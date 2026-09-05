#!/usr/bin/env python3
"""GOCI-III 편광 논문 grid 재생성.

설계 사양(2026-07-18 확정):
- 1000케이스, 각 행이 완결 조건.
- 해수 IOP는 대양->연안 covary: 로그공간 상관 다변량 정규 + 범위 truncation.
  TSM 0.1-20 g/m3, Chl 0.1-5 mg/m3, aCDOM440 0.01-0.1 /m.
  로그공간 Pearson 상관 TSM-Chl 0.85, Chl-aCDOM 0.75, TSM-aCDOM 0.77.
- 에어로졸 6종(C50/T50/M80C/M95C/M98C/O99) 균등 랜덤.
- AOD865 0.05-0.3 균등, 풍속 1-10 m/s 균등.
- 기하 sza 0-65, vza 0-55, raa 0-180 중 sunglint-free
  (Cox-Munk 직접 glint rho_g < 5e-4; 거울점 회피).
- 해수는 covary, 대기/기하는 독립.

주의: 원본 grid의 설계 스크립트와 난수 시드가 이전 세션에서 소실되어, 원본과
비트 동일한 재현은 불가능하다. 이 스크립트는 위 사양(범위·상관·sunglint 기준)을
만족하는 새 샘플을 생성한다. 시드를 고정해 이 스크립트 자체는 재현 가능하다.
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ocrt_py.surface import direct_sunglint_rho

AEROSOLS = ['C50', 'T50', 'M80C', 'M95C', 'M98C', 'O99']

# 해수 IOP 범위 (min, max)
TSM_RANGE = (0.1, 20.0)     # g/m3
CHL_RANGE = (0.1, 5.0)      # mg/m3
ACDOM_RANGE = (0.01, 0.1)   # /m

# 로그공간 Pearson 상관 (순서: TSM, Chl, aCDOM)
CORR = np.array([
    [1.00, 0.85, 0.77],
    [0.85, 1.00, 0.75],
    [0.77, 0.75, 1.00],
])

# 기하 범위
SZA_MAX = 65.0
VZA_MAX = 55.0
RAA_MAX = 180.0

# sunglint-free 기준
GLINT_THRESH = 5.0e-4
N_WATER_GLINT = 1.34     # 판정용 대표 굴절률
SIGMA_TYPE = 1
Q_CONV = 1


def _logspace_params(lo, hi, k_sigma=2.0):
    """범위 [lo, hi]를 로그공간에서 ±k_sigma로 덮는 (mu, sigma) (자연로그)."""
    lmu = 0.5 * (np.log(lo) + np.log(hi))
    lsig = (np.log(hi) - np.log(lo)) / (2.0 * k_sigma)
    return lmu, lsig


def sample_iop(n, rng):
    """covary 해수 IOP (TSM, Chl, aCDOM) 샘플. 범위 밖은 재샘플(truncation)."""
    mu = np.zeros(3)
    sig = np.zeros(3)
    for j, (lo, hi) in enumerate([TSM_RANGE, CHL_RANGE, ACDOM_RANGE]):
        mu[j], sig[j] = _logspace_params(lo, hi)
    # 로그공간 공분산 = D @ CORR @ D  (D=diag(sig))
    D = np.diag(sig)
    cov = D @ CORR @ D
    L = np.linalg.cholesky(cov)

    ranges = [TSM_RANGE, CHL_RANGE, ACDOM_RANGE]
    out = np.zeros((n, 3))
    filled = 0
    while filled < n:
        # 여유있게 뽑고 truncation
        m = (n - filled) * 3
        z = rng.standard_normal((m, 3))
        logx = mu[None, :] + z @ L.T
        x = np.exp(logx)
        ok = np.ones(m, dtype=bool)
        for j, (lo, hi) in enumerate(ranges):
            ok &= (x[:, j] >= lo) & (x[:, j] <= hi)
        good = x[ok]
        take = min(len(good), n - filled)
        out[filled:filled + take] = good[:take]
        filled += take
    return out[:, 0], out[:, 1], out[:, 2]   # TSM, Chl, aCDOM


def sample_geometry(n, wind, rng):
    """sunglint-free 기하 (sza, vza, raa) 샘플. glint rho_g < 5e-4 재샘플."""
    sza = np.zeros(n); vza = np.zeros(n); raa = np.zeros(n)
    for i in range(n):
        while True:
            s = rng.uniform(0.0, SZA_MAX)
            v = rng.uniform(0.0, VZA_MAX)
            r = rng.uniform(0.0, RAA_MAX)
            mu0 = np.cos(np.radians(s))
            muv = np.cos(np.radians(v))
            phi_v = np.radians(r - 180.0)   # driver 관례
            # 순수 표면 glint 반사도 (tau=0)
            rg = direct_sunglint_rho(muv, phi_v, mu0, 0.0, wind[i],
                                     SIGMA_TYPE, N_WATER_GLINT, 0.0, Q_CONV)
            if abs(rg[0]) < GLINT_THRESH:
                sza[i] = s; vza[i] = v; raa[i] = r
                break
    return sza, vza, raa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=1000, help='케이스 수 (기본 1000)')
    ap.add_argument('--seed', type=int, default=20260718, help='난수 시드')
    ap.add_argument('--out', default='full_grid_design_v1.csv')
    ap.add_argument('--iop-out', default='iop_grid_design_v1.csv',
                    help='해수 IOP만 담은 부속 CSV')
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n = args.n

    # 대기/기하 (독립)
    aer = rng.choice(AEROSOLS, size=n)
    aod865 = rng.uniform(0.05, 0.30, size=n)
    wind = rng.uniform(1.0, 10.0, size=n)

    # 해수 IOP (covary)
    tsm, chl, acdom = sample_iop(n, rng)

    # 기하 (sunglint-free)
    sza, vza, raa = sample_geometry(n, wind, rng)

    # full grid
    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case_id', 'sza', 'vza', 'raa', 'wind',
                    'aerosol', 'aod865', 'chl', 'tsm', 'acdom440'])
        for i in range(n):
            w.writerow([f'G{i + 1:04d}',
                        f'{sza[i]:.4f}', f'{vza[i]:.4f}', f'{raa[i]:.4f}',
                        f'{wind[i]:.4f}', aer[i], f'{aod865[i]:.5f}',
                        f'{chl[i]:.5f}', f'{tsm[i]:.5f}', f'{acdom[i]:.6f}'])

    # iop 부속 (검증·기록용)
    with open(args.iop_out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['case_id', 'TSM_g_m3', 'Chl_mg_m3', 'aCDOM440_m_inv'])
        for i in range(n):
            w.writerow([f'G{i + 1:04d}', f'{tsm[i]:.5f}',
                        f'{chl[i]:.5f}', f'{acdom[i]:.6f}'])

    # 통계 요약
    print(f'[make_grid] {n} 케이스 생성 (seed={args.seed}) -> {args.out}')
    logs = np.log(np.vstack([tsm, chl, acdom]))
    C = np.corrcoef(logs)
    print('  해수 IOP 로그공간 상관 (목표 vs 실현):')
    print(f'    TSM-Chl:   {CORR[0,1]:.2f} vs {C[0,1]:.3f}')
    print(f'    Chl-aCDOM: {CORR[1,2]:.2f} vs {C[1,2]:.3f}')
    print(f'    TSM-aCDOM: {CORR[0,2]:.2f} vs {C[0,2]:.3f}')
    print('  범위 확인:')
    print(f'    TSM   [{tsm.min():.3f}, {tsm.max():.3f}]  (목표 {TSM_RANGE})')
    print(f'    Chl   [{chl.min():.3f}, {chl.max():.3f}]  (목표 {CHL_RANGE})')
    print(f'    aCDOM [{acdom.min():.4f}, {acdom.max():.4f}]  (목표 {ACDOM_RANGE})')
    print(f'    sza   [{sza.min():.2f}, {sza.max():.2f}]  vza [{vza.min():.2f}, {vza.max():.2f}]  raa [{raa.min():.2f}, {raa.max():.2f}]')
    print(f'    aod865[{aod865.min():.3f}, {aod865.max():.3f}]  wind [{wind.min():.2f}, {wind.max():.2f}]')
    from collections import Counter
    print(f'  에어로졸 분포: {dict(Counter(aer))}')


if __name__ == '__main__':
    main()
