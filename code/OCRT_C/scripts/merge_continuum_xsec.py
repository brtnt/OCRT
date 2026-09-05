#!/usr/bin/env python3
"""
merge_continuum_xsec.py — layered LBL xsec 테이블에 실측 UV-Vis continuum 합산

배경:
  HITRAN line-by-line(LBL)로 생성한 layered xsec_o3.dat / xsec_no2.dat 에는
  가시광 전자전이 continuum(O3 Chappuis/Huggins, NO2 visible)이 없다(4~5자릿수 과소).
  continuum 은 line 흡수와 독립적으로 더해지는 성분이므로, 실측 continuum
  cross section 을 층별 온도로 내삽하여 LBL 테이블에 합산한다:
      sigma_total(layer, λ) = sigma_LBL(layer, λ) + sigma_cont(T_layer, λ)

필요한 실측 데이터(사용자가 다운로드하여 준비):
  O3 : Serdyuchenko / Gorshelev et al. 2014 (AMT 7, 609 & 625).
       213–1100 nm, 온도 11개(193–293 K, 10 K 간격). IUP Bremen 또는
       HITRAN UV-Vis cross-section 섹션에서 배포.
       보통 단일 파일: 1열 wavelength[nm] + 11열 sigma[cm^2](193K..293K).
  NO2: Vandaele et al. 1998 (JQSRT 59, 171). 238–1000 nm, 220 K / 294 K 두 파일.
       HITRAN UV-Vis 또는 BIRA-IASB 배포. 보통 1열 wavenumber[cm^-1] + 1열 sigma.

사용 예:
  # O3 (단일 파일, 온도열 11개)
  python3 merge_continuum_xsec.py \
      --lbl inputs/xsec/xsec_o3.dat \
      --cont serdyuchenko_o3.dat --temps 193,203,213,223,233,243,253,263,273,283,293 \
      --out inputs/xsec/xsec_o3.dat.new

  # NO2 (온도별 파일 2개: 파일경로:온도K)
  python3 merge_continuum_xsec.py \
      --lbl inputs/xsec/xsec_no2.dat \
      --cont-t no2_220K.dat:220 --cont-t no2_294K.dat:294 \
      --out inputs/xsec/xsec_no2.dat.new

규칙:
  * 내삽은 전부 선형(파장 boxcar 평균 + 온도 선형). nearest 금지.
  * 온도 격자 밖은 경계값 고정(clamp)하고 경고 출력(외삽 금지).
  * continuum 데이터의 파장범위 밖은 0 으로 두고 해당 구간을 경고 출력.
  * 합산 결과 음수 발생 시 오류로 중단(비물리).
"""
import sys, argparse
import numpy as np


# ---------------------------------------------------------------- I/O helpers
def read_layered_lbl(path):
    """layered xsec 파일: 헤더(#) + 1행 파장격자 + n_layers행 sigma. 층 T는 헤더에서."""
    header, rows, layer_T = [], [], []
    for line in open(path):
        if line.startswith('#'):
            header.append(line.rstrip('\r\n'))
            if 'layer' in line and 'T=' in line and 'z=' in line:
                try:
                    layer_T.append(float(line.split('T=')[1].split('K')[0]))
                except (IndexError, ValueError):
                    pass
            continue
        if line.strip():
            rows.append(line.split())
    wl = np.array(rows[0], dtype=float)
    sig = np.array(rows[1:], dtype=float)
    if len(layer_T) != sig.shape[0]:
        sys.exit(f'error: 헤더의 층 T 개수({len(layer_T)}) != 데이터 층 수({sig.shape[0]}). '
                 f'"# layer k: z=.. T=..K" 헤더가 있는 layered 파일이 필요하다.')
    return header, wl, sig, np.array(layer_T)


def read_cont_table(path):
    """continuum 파일: 주석(#,;,%,문자행) 무시, 수치행만. 1열 wl 또는 wavenumber + sigma열(들)."""
    rows = []
    for line in open(path, errors='replace'):
        t = line.replace(',', ' ').split()
        if not t:
            continue
        try:
            rows.append([float(x) for x in t])
        except ValueError:
            continue                       # 헤더/텍스트 행
    if not rows:
        sys.exit(f'error: {path} 에서 수치 데이터를 찾지 못했다.')
    ncol = min(len(r) for r in rows)
    arr = np.array([r[:ncol] for r in rows], dtype=float)
    x = arr[:, 0]
    if np.nanmedian(x) > 2000.0:           # wavenumber cm^-1 로 판단 → nm 변환
        arr[:, 0] = 1.0e7 / x
    arr = arr[np.argsort(arr[:, 0])]
    return arr[:, 0], arr[:, 1:]


# ---------------------------------------------------------------- resampling
def boxcar_to_grid(wl_src, sig_src, wl_dst, half_width=0.5):
    """소스 고분해능 sigma 를 목적 격자에 ±half_width nm boxcar 평균으로 재표본.
    창 안에 소스점이 없으면(성긴 데이터) 선형 내삽으로 대체. nearest 미사용."""
    out = np.zeros_like(wl_dst)
    idx_lo = np.searchsorted(wl_src, wl_dst - half_width, side='left')
    idx_hi = np.searchsorted(wl_src, wl_dst + half_width, side='right')
    lin = np.interp(wl_dst, wl_src, sig_src)          # 성긴 구간 fallback (선형)
    for i, (lo, hi) in enumerate(zip(idx_lo, idx_hi)):
        out[i] = sig_src[lo:hi].mean() if hi - lo >= 2 else lin[i]
    inside = (wl_dst >= wl_src[0] - half_width) & (wl_dst <= wl_src[-1] + half_width)
    out[~inside] = 0.0
    return out, inside


def interp_T_masked(sig_by_T, avail, temps, T_layer):
    """파장별 가용 온도 부분집합 내 선형 내삽(부분집합 경계 밖은 clamp).
    sig_by_T, avail: (n_T, n_wl). 반환: (sigma(n_wl), clamp 발생 여부).
    데이터가 없는 파장(가용 온도 0개)은 0. nearest 미사용 — 가용 격자 내 선형."""
    temps = np.asarray(temps, dtype=float)
    order = np.argsort(temps)
    temps, sig_by_T, avail = temps[order], sig_by_T[order], avail[order]
    n_T, n_wl = sig_by_T.shape
    out = np.zeros(n_wl)
    clamped_any = False
    # 가용 마스크 패턴별로 묶어 벡터화 (범위가 온도별 연속 구간이라 패턴 수 적음)
    keys = {}
    for j in range(n_wl):
        keys.setdefault(tuple(avail[:, j]), []).append(j)
    for mask, cols in keys.items():
        idx = [k for k in range(n_T) if mask[k]]
        cols = np.array(cols)
        if not idx:
            continue                                  # 데이터 전무 → 0
        t_sub = temps[idx]
        s_sub = sig_by_T[np.ix_(idx, cols)]
        if T_layer <= t_sub[0]:
            out[cols] = s_sub[0]; clamped_any = True
        elif T_layer >= t_sub[-1]:
            out[cols] = s_sub[-1]; clamped_any = True
        else:
            j2 = int(np.searchsorted(t_sub, T_layer)); j1 = j2 - 1
            f = (T_layer - t_sub[j1]) / (t_sub[j2] - t_sub[j1])
            out[cols] = (1.0 - f) * s_sub[j1] + f * s_sub[j2]
    return out, clamped_any


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lbl', required=True, help='입력 layered LBL xsec 파일')
    ap.add_argument('--out', required=True, help='출력 파일 (합산본)')
    ap.add_argument('--cont', help='continuum 단일 파일(1열 wl/wn + 온도별 sigma 열들)')
    ap.add_argument('--temps', help='--cont 의 sigma 열 온도 목록 [K], 예: 193,203,...,293')
    ap.add_argument('--cont-t', action='append', default=[],
                    help='온도별 개별 파일, "경로:온도K" 형식. 반복 지정 가능')
    args = ap.parse_args()

    header, wl, sig_lbl, layer_T = read_layered_lbl(args.lbl)
    print(f'LBL: {args.lbl}  wl {wl[0]:.0f}-{wl[-1]:.0f}nm x {sig_lbl.shape[0]}층, '
          f'층 T {layer_T.min():.1f}-{layer_T.max():.1f}K')

    # continuum 로드 → (temps, sigma_on_dst_grid)
    if args.cont:
        if not args.temps:
            sys.exit('error: --cont 사용 시 --temps 필수 (sigma 열 온도 목록)')
        temps = [float(t) for t in args.temps.split(',')]
        wl_c, sig_c = read_cont_table(args.cont)
        if sig_c.shape[1] < len(temps):
            sys.exit(f'error: continuum sigma 열 {sig_c.shape[1]}개 < 지정 온도 {len(temps)}개')
        sig_c = sig_c[:, :len(temps)]
        srcs = [(temps[k], wl_c, sig_c[:, k]) for k in range(len(temps))]
    elif args.cont_t:
        srcs = []
        for spec in args.cont_t:
            path, T = spec.rsplit(':', 1)
            wl_c, sig_c = read_cont_table(path)
            srcs.append((float(T), wl_c, sig_c[:, 0]))
    else:
        sys.exit('error: --cont 또는 --cont-t 중 하나는 필요하다.')

    temps = [s[0] for s in srcs]
    grid_sig, grid_avail = [], []
    for T, wl_c, s_c in srcs:
        # 실측 잡음 음수 클리핑: sigma 는 물리적으로 비음(측정잡음 규모 확인 후 적용;
        # O3 min -1.7e-24 = Chappuis 정점의 3.4e-4배, NO2 min ~ 정점의 1% 대역날개)
        n_neg = int((s_c < 0).sum())
        if n_neg:
            print(f'  T={T:.0f}K: 원천 음수 {n_neg}점(측정잡음) → 0 클리핑 '
                  f'(최솟값 {s_c.min():.2e})')
            s_c = np.maximum(s_c, 0.0)
        g, inside = boxcar_to_grid(wl_c, s_c, wl)
        grid_sig.append(g); grid_avail.append(inside)
        print(f'  continuum T={T:.0f}K: 원천 {wl_c[0]:.1f}-{wl_c[-1]:.1f}nm, '
              f'{len(wl_c)}점 → 1nm boxcar 재표본')
    grid_sig = np.array(grid_sig); grid_avail = np.array(grid_avail)

    none_cov = ~grid_avail.any(axis=0)
    if none_cov.any():
        miss = wl[none_cov]
        print(f'  WARNING: 모든 온도에서 미포함 파장 {miss.min():.0f}-{miss.max():.0f}nm '
              f'({none_cov.sum()}점) → continuum=0 유지')
    part_cov = grid_avail.any(axis=0) & ~grid_avail.all(axis=0)
    if part_cov.any():
        miss = wl[part_cov]
        print(f'  NOTE: 일부 온도만 포함하는 파장 {miss.min():.0f}-{miss.max():.0f}nm '
              f'({part_cov.sum()}점) → 가용 온도 부분집합 내 선형 내삽/clamp')

    clamped = 0
    sig_out = sig_lbl.copy()
    for k, Tk in enumerate(layer_T):
        s_T, was_clamped = interp_T_masked(grid_sig, grid_avail, temps, Tk)
        clamped += int(was_clamped)
        sig_out[k] += s_T
    if clamped:
        print(f'  NOTE: 층 {clamped}개에서 온도 clamp 발생(가용격자 밖, 외삽 금지)')

    if (sig_out < 0).any():
        sys.exit('error: 합산 결과에 음수 sigma 발생 — 입력 데이터 확인 필요 (비물리)')

    with open(args.out, 'w') as f:
        for h in header:
            f.write(h + '\n')
        f.write(f'# CONTINUUM MERGED by merge_continuum_xsec.py\n')
        f.write(f'#   sources: ' + '; '.join(
            (f'{args.cont} (T={args.temps}K)' if args.cont else
             ', '.join(args.cont_t))) + '\n')
        f.write(f'#   sigma_total = sigma_LBL + sigma_continuum(T_layer, lambda); '
                f'linear-in-T, 1nm boxcar-in-lambda; no nearest\n')
        f.write('  '.join(f'{x:.4f}' for x in wl) + '\n')
        for k in range(sig_out.shape[0]):
            f.write('  '.join(f'{x:.6e}' for x in sig_out[k]) + '\n')
    print(f'출력: {args.out}')

    # sanity: 대표 파장 sigma 출력 (지상층)
    for lam in (440.0, 555.0, 602.0):
        j = int(np.argmin(np.abs(wl - lam)))
        print(f'  sanity sigma(ground, {wl[j]:.0f}nm): '
              f'LBL={sig_lbl[0, j]:.3e} → merged={sig_out[0, j]:.3e} cm^2')


if __name__ == '__main__':
    main()
