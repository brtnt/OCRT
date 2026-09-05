#!/usr/bin/env python3
"""GPU 무결성 게이트.

같은 grid를 CPU(numpy)와 GPU(CuPy)로 각각 실행해 ρ_TOA/ρ_R/ρ_R+A/Rrs를
7유효자리(상대오차 1e-7)로 대조한다. 배치화의 부동소수 순서차는 수렴차수를
±1 어긋나게 할 수 있어 비트 일치가 아니라 7유효자리를 무결성 기준으로 둔다.

CuPy가 설치된 로컬 CUDA 환경에서 실행한다. 샌드박스에는 GPU가 없어 이 게이트를
통과시키는 것이 로컬 검증의 필수 단계다.

사용법:
    python gpu_gate.py --grid small_grid.csv --data data
    (권장: 4~16행 정도의 소규모 grid)
"""
import argparse
import os
import sys
import time

import numpy as np

# produce_grid의 내부 유틸을 재사용한다.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)


def _run(grid_path, data_dir, use_gpu, bands, chunk, n_mu_water, nt_atm,
         m_max, max_it_water):
    """한 backend로 grid를 풀어 {(case_id, band): (rTOA, rR, rRA, Rrs)}를 반환한다."""
    os.environ['OCRT_PY_GPU'] = '1' if use_gpu else '0'
    # backend는 import 시점에 OCRT_PY_GPU를 읽으므로, 하위 모듈을 매번 새로
    # 로드하기 위해 관련 모듈을 캐시에서 제거한다.
    for m in list(sys.modules):
        if m.startswith('ocrt_py') or m == 'produce_grid':
            del sys.modules[m]
    import produce_grid as PG
    from ocrt_py import batch_driver as BD
    from ocrt_py.constituent import OCRTConstituentModel

    rows = PG._read_grid(grid_path)
    bands = bands or PG.BANDS_NM
    mie_cache = PG._mie_cache(rows, data_dir)
    cm = OCRTConstituentModel(data_dir, 'micro', 'red_clay')
    from ocrt_py.absorption import Absorption
    ab = Absorption(os.path.join(data_dir, 'afgl_atm'),
                    os.path.join(data_dir, 'xsec'))

    out = {}
    for band in bands:
        for chunk_rows in PG._chunks(rows, chunk):
            c2 = PG._band_cases(chunk_rows, band, mie_cache, 'r2')
            c3 = PG._band_cases(chunk_rows, band, mie_cache, 'r3')
            c1 = PG._band_cases(chunk_rows, band, mie_cache, 'r1')
            rR = BD.solve_r2_grid(c2, pressure_hpa=1013.25,
                                  n_mu_gl=n_mu_water, absorption=ab)
            rRA = BD.solve_r3_grid(c3, mie_cache, n_mu_gl=n_mu_water,
                                   nt=nt_atm, m_max=m_max,
                                   max_iterations=100, absorption=ab)
            rTOA, rrs, _tr, _iop = BD.solve_r1_grid(c1, mie_cache, cm,
                                         n_mu_water=n_mu_water, nt_atm=nt_atm,
                                         fourier_m_max=m_max, max_it_atm=100,
                                         max_it_water=max_it_water,
                                         absorption=ab)
            for k, r in enumerate(chunk_rows):
                out[(r['case_id'], band)] = (float(rTOA[k]), float(rR[k]),
                                             float(rRA[k]), float(rrs[k]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--grid', required=True, help='소규모 grid CSV (4~16행 권장)')
    ap.add_argument('--data', default='data')
    ap.add_argument('--bands', default='443,555,865',
                    help='대조할 밴드 nm (기본 443,555,865)')
    ap.add_argument('--chunk', type=int, default=8)
    ap.add_argument('--n-mu-water', type=int, default=24)
    ap.add_argument('--nt-atm', type=int, default=400)
    ap.add_argument('--m-max', type=int, default=16)
    ap.add_argument('--max-it-water', type=int, default=500)
    ap.add_argument('--tol', type=float, default=1e-7,
                    help='통과 기준 상대오차 (기본 1e-7 = 7유효자리)')
    ap.add_argument('--gpu-only', action='store_true',
                    help='CPU 대조를 생략하고 GPU 시간만 측정한다(큰 B 처리량 측정용). '
                         'CPU 기준 실행이 큰 B에서 지나치게 느리므로, 처리량만 볼 때 사용')
    args = ap.parse_args()
    bands = [int(x) for x in args.bands.split(',')]

    kw = dict(bands=bands, chunk=args.chunk, n_mu_water=args.n_mu_water,
              nt_atm=args.nt_atm, m_max=args.m_max,
              max_it_water=args.max_it_water)

    # CuPy 사용 가능 여부 확인 (GPU를 먼저 실행해 GPU 동작을 바로 볼 수 있게 한다)
    try:
        import cupy  # noqa: F401
    except Exception as e:
        print(f'[gpu_gate] CuPy import 실패: {e}')
        print('  CuPy를 설치해야 GPU 게이트를 실행할 수 있다.')
        print('  예: pip install cupy-cuda12x')
        sys.exit(2)

    print('[gpu_gate] GPU(CuPy) 실행 중...')
    t0 = time.time()
    gpu = _run(args.grid, args.data, True, **kw)
    t_gpu = time.time() - t0
    n_cb = max(len(gpu), 1)
    print(f'  GPU 완료: {len(gpu)} 케이스, {t_gpu:.1f}s '
          f'(평균 {t_gpu / n_cb * 1000:.0f} ms/케이스, chunk(B)={args.chunk})')

    if args.gpu_only:
        print('[gpu_gate] --gpu-only 모드로 CPU 대조를 생략하였다. '
              'chunk가 커질수록 케이스당 시간이 감소하면 GPU 활용률이 향상된 것이다.')
        return

    print('[gpu_gate] CPU(numpy) 기준 대조 실행 중...')
    t0 = time.time()
    cpu = _run(args.grid, args.data, False, **kw)
    t_cpu = time.time() - t0
    print(f'  CPU 완료: {len(cpu)} 케이스, {t_cpu:.1f}s')

    # 대조
    labels = ['rho_TOA', 'rho_R', 'rho_RpA', 'Rrs']
    worst = [0.0, 0.0, 0.0, 0.0]
    worst_key = [None, None, None, None]
    n = 0
    for key in cpu:
        if key not in gpu:
            continue
        n += 1
        for j in range(4):
            a = cpu[key][j]
            b = gpu[key][j]
            denom = abs(a) if abs(a) > 0 else 1.0
            rel = abs(a - b) / denom
            if rel > worst[j]:
                worst[j] = rel
                worst_key[j] = key

    print(f'\n[gpu_gate] 대조 결과 ({n} 케이스, 기준 rel < {args.tol:.0e})')
    all_pass = True
    for j in range(4):
        status = 'PASS' if worst[j] < args.tol else 'FAIL'
        if worst[j] >= args.tol:
            all_pass = False
        wk = worst_key[j]
        print(f'  {labels[j]:8s}: 최대 rel = {worst[j]:.2e}  [{status}]'
              + (f'  (worst: case {wk[0]} band {wk[1]})' if wk else ''))
    print(f'\n속도: CPU {t_cpu:.1f}s / GPU {t_gpu:.1f}s'
          + (f' (가속 {t_cpu / t_gpu:.1f}x)' if t_gpu > 0 else ''))
    print('결과:', 'PASS (7유효자리 통과)' if all_pass else 'FAIL (기준 초과)')
    sys.exit(0 if all_pass else 1)


if __name__ == '__main__':
    main()
