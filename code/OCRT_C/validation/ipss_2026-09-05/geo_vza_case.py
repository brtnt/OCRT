#!/usr/bin/env python3
"""실제 기하로 표층 VZA vs TOA VZA 계산 (2026-09-05).

표적 T: 36.0 N, 128.25 E (지표)
센서 S: 0.0 N, 128.25 E (같은 경도, 적도) — 고도 (a) 정지궤도 35786 km, (b) 700 km

정의
  gamma          지심각 (OT 와 OS 사이 각) = 같은 경도이므로 위도차
  theta_surface  표적의 국소 천정각 = 위성 L1B senz
  theta_TOA      같은 직선이 r = Re + H 를 지나는 점의 국소 천정각
                 r sin(theta) = const (직선 불변량)  =>  sin(theta_TOA) = Re/(Re+H) sin(theta_surface)
"""
import numpy as np

RE = 6371.0          # OCRT 모델 지구 반경
H = 100.0            # OCRT 모델 TOA 고도
D = np.degrees
R = np.radians


def theta_surface(gamma_deg, h_sat):
    """표적에서 위성을 볼 때의 국소 천정각. 90도 초과면 지평선 아래(관측 불가)."""
    g = R(gamma_deg)
    Rs = RE + h_sat
    d = np.sqrt(RE ** 2 + Rs ** 2 - 2 * RE * Rs * np.cos(g))
    return D(np.arccos(np.clip((Rs * np.cos(g) - RE) / d, -1.0, 1.0))), d


def theta_toa(th_s_deg, H_km=H):
    s = RE / (RE + H_km) * np.sin(R(th_s_deg))
    return D(np.arcsin(np.clip(s, -1.0, 1.0)))


def horizon_gamma(h_sat):
    return D(np.arccos(RE / (RE + h_sat)))


def nadir_angle(gamma_deg, h_sat):
    """위성에서 본 천저각(look angle)."""
    g = R(gamma_deg)
    Rs = RE + h_sat
    d = np.sqrt(RE ** 2 + Rs ** 2 - 2 * RE * Rs * np.cos(g))
    return D(np.arcsin(np.clip(RE * np.sin(g) / d, -1.0, 1.0)))


GAMMA = 36.0
print("=" * 78)
print(f"표적 36.0N 128.25E / 센서 0.0N 128.25E  =>  지심각 gamma = {GAMMA:.2f} deg")
print(f"구면 지구 Re = {RE:.0f} km, 모델 TOA H = {H:.0f} km")
print("=" * 78)

for h_sat, name in ((35786.0, "정지궤도 (GEO)"), (700.0, "저궤도 700 km")):
    hg = horizon_gamma(h_sat)
    ths, d = theta_surface(GAMMA, h_sat)
    print(f"\n[{name}]  h_sat = {h_sat:.0f} km   (Rs = {RE + h_sat:.0f} km)")
    print(f"  지평선 한계 지심각          = {hg:8.3f} deg")
    if GAMMA >= hg:
        print(f"  gamma = {GAMMA:.2f} deg > {hg:.3f} deg  ==>  표적이 지평선 아래. 관측 불가.")
        print(f"  (형식적으로 계산하면 theta_surface = {ths:.3f} deg > 90 deg)")
        print(f"  이 고도에서 36.0N 을 보려면 위성이 위도 {90.0 - 0:.0f}..? — "
              f"지심각 {hg:.2f} deg 안에 있어야 하므로 위성 위도가 "
              f"{GAMMA - hg:.2f} deg 이상이어야 한다.")
        continue
    tht = theta_toa(ths)
    print(f"  슬랜트 거리 |TS|            = {d:10.1f} km")
    print(f"  위성 천저각 (look angle)    = {nadir_angle(GAMMA, h_sat):8.3f} deg")
    print(f"  theta_surface (L1B senz)    = {ths:8.4f} deg")
    print(f"  theta_TOA (h = 100 km)      = {tht:8.4f} deg")
    print(f"  차이 Delta theta            = {ths - tht:8.4f} deg")

# ---------------------------------------------------------------- 핵심 대조
print("\n" + "=" * 78)
print("핵심: Delta theta 는 theta_surface 만의 함수다 (h_sat 무관)")
print("=" * 78)
ths_geo, _ = theta_surface(GAMMA, 35786.0)
print(f"{'h_sat [km]':>12} {'gamma [deg]':>12} {'th_surface':>11} {'th_TOA':>9} {'Delta':>8}")
for h_sat in (700.0, 2000.0, 5000.0, 35786.0):
    # 같은 theta_surface 를 만드는 gamma 를 역산
    lo, hi = 1e-6, horizon_gamma(h_sat) - 1e-6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        t, _d = theta_surface(mid, h_sat)
        if t < ths_geo:
            lo = mid
        else:
            hi = mid
    g = 0.5 * (lo + hi)
    t, _d = theta_surface(g, h_sat)
    print(f"{h_sat:12.0f} {g:12.3f} {t:11.4f} {theta_toa(t):9.4f} {t - theta_toa(t):8.4f}")
print("=> theta_surface 가 같으면 고도가 700 km 든 GEO 든 Delta theta 가 완전히 같다.")

# ------------------------------------------------- 사용자 직관이 맞는 부분
print("\n" + "=" * 78)
print(f"사용자 직관이 맞는 부분: gamma = {GAMMA:.1f} deg 고정하고 센서를 올리면")
print("theta_surface 가 gamma 로 수렴하고, 그에 따라 Delta theta 도 줄어든다")
print("=" * 78)
print(f"{'h_sat [km]':>12} {'th_surface':>11} {'th_TOA':>9} {'Delta':>8}  비고")
rows = []
for h_sat in (1600.0, 2000.0, 3000.0, 5000.0, 10000.0, 20000.0, 35786.0, 1.0e7):
    hg = horizon_gamma(h_sat)
    if GAMMA >= hg:
        print(f"{h_sat:12.0f} {'—':>11} {'—':>9} {'—':>8}  지평선 아래 (한계 {hg:.2f} deg)")
        continue
    t, _d = theta_surface(GAMMA, h_sat)
    tt = theta_toa(t)
    tag = "무한 극한" if h_sat > 1e6 else ""
    print(f"{h_sat:12.0f} {t:11.4f} {tt:9.4f} {t - tt:8.4f}  {tag}")
    rows.append((h_sat, t, tt))
print(f"{'':12} {'':11} {'':9} {'':8}  h_sat -> inf 에서 theta_surface -> gamma = {GAMMA:.1f} deg")
print(f"36.0 N 이 보이기 시작하는 최소 고도 = "
      f"{RE / np.cos(R(GAMMA)) - RE:.0f} km")

# ------------------------------------------------- Delta theta(theta_surface)
print("\n" + "=" * 78)
print("Delta theta(theta_surface), H = 100 km — 고도는 어디에도 들어가지 않는다")
print("=" * 78)
print(f"{'th_surface':>11} {'th_TOA':>9} {'Delta':>8}")
for t in (0, 10, 20, 30, 36, 40, 50, 55, 60, 70, 80):
    print(f"{t:11.1f} {theta_toa(t):9.4f} {t - theta_toa(t):8.4f}")

# WGS84 타원체 보정 크기
f = 1.0 / 298.257223563
phi_c = D(np.arctan((1 - f) ** 2 * np.tan(R(36.0))))
print(f"\n[타원체 참고] WGS84 측지위도 36.000 N -> 지심위도 {phi_c:.4f} N "
      f"(gamma 가 {36.0 - phi_c:.3f} deg 작아진다)")
t_sph, _ = theta_surface(36.0, 35786.0)
t_ell, _ = theta_surface(phi_c, 35786.0)
print(f"  GEO theta_surface: 구면 {t_sph:.4f} deg / 지심위도 사용 {t_ell:.4f} deg "
      f"(차 {t_sph - t_ell:.4f} deg)")
print(f"  그에 따른 Delta theta: {t_sph - theta_toa(t_sph):.4f} vs "
      f"{t_ell - theta_toa(t_ell):.4f} deg — 결론에 영향 없음")
