#!/usr/bin/env python3
"""rt_ipss.c 와 무관한 독립 검산 — 순수 Python 3D 수치적분.

질문: 지상(화소) 앵커에서 VZA 55, SZA 0 에 남는 -0.24 % 는 물리인가 오류인가.
방법: 구면쉘 단일산란 적분을 좌표 벡터로 직접 수치적분(레이-셸 교차 루틴 없이
      촘촘한 Simpson), 평면평행 폐형식과 비교.  위상함수·E0 는 비에서 상쇄되므로 생략.
"""
import numpy as np

RE, H = 6371.0, 100.0


def beta_of_h(h, tau_t, Hs):
    """지수 대기, 연직 tau = tau_t 로 정규화. Hs<=0 이면 균질."""
    if Hs > 0:
        b0 = tau_t / (Hs * (1.0 - np.exp(-H / Hs)))
        return b0 * np.exp(-h / Hs)
    return np.full_like(h, tau_t / H)


def tau_along(p0, d, tau_t, Hs, n=4001):
    """점 p0 에서 단위벡터 d 방향으로 대기 밖(|r|>Re+H)까지의 광학두께.
    지표(|r|<Re)에 닿으면 inf (차폐)."""
    # 구간 끝: |p0 + t d| = Re + H 의 양의 근
    b = p0 @ d
    c = p0 @ p0 - (RE + H) ** 2
    disc = b * b - c
    t_exit = -b + np.sqrt(disc)
    # 지표 충돌 검사: |p0 + t d| = Re 의 근이 (0, t_exit) 안에 있으면 차폐
    c2 = p0 @ p0 - RE ** 2
    disc2 = b * b - c2
    if disc2 > 0:
        t1 = -b - np.sqrt(disc2)
        if 1e-9 < t1 < t_exit:
            return np.inf
    t = np.linspace(0.0, t_exit, n)
    r = p0[None, :] + t[:, None] * d[None, :]
    h = np.linalg.norm(r, axis=1) - RE
    f = beta_of_h(h, tau_t, Hs)
    return np.trapezoid(f, t)


def I1_sphere(vza, sza, raa, tau_t, Hs, n_los=2001):
    """화소 (0,0,Re) 앵커. 태양은 화소에서 국소 천정각 sza 인 평행빔."""
    tv, ts, fp = np.radians(vza), np.radians(sza), np.radians(180.0 - raa)
    u = np.array([np.sin(tv) * np.cos(fp), np.sin(tv) * np.sin(fp), np.cos(tv)])
    to_sun = np.array([np.sin(ts), 0.0, np.cos(ts)])   # 화소에서 태양을 향하는 단위벡터
    P = np.array([0.0, 0.0, RE])
    b = P @ u; c = P @ P - (RE + H) ** 2
    s_max = -b + np.sqrt(b * b - c)
    s = np.linspace(0.0, s_max, n_los)
    r = P[None, :] + s[:, None] * u[None, :]
    h = np.linalg.norm(r, axis=1) - RE
    beta = beta_of_h(h, tau_t, Hs)
    # 시선 감쇠: 산란점 s 에서 TOA(s_max) 까지
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (beta[1:] + beta[:-1]) * np.diff(s))])
    tau_view = cum[-1] - cum
    # 태양 감쇠: 각 점에서 to_sun 방향으로 대기 밖까지 (독립 적분)
    tau_sun = np.array([tau_along(r[i], to_sun, tau_t, Hs) for i in range(n_los)])
    integrand = beta * np.exp(-tau_view - tau_sun)
    integrand[~np.isfinite(integrand)] = 0.0
    return np.trapezoid(integrand, s)


def I1_pp(vza, sza, tau_t):
    mv, ms = np.cos(np.radians(vza)), np.cos(np.radians(sza))
    C = 1.0 / mv + 1.0 / ms
    return (1.0 / (mv * C)) * (1.0 - np.exp(-C * tau_t))


if __name__ == "__main__":
    tau = 0.0935
    print("독립 Python 적분  —  dev = pp/sph − 1  [%]   (rt_ipss.c 미사용)")
    print(f"{'VZA':>5} {'SZA':>5} {'Hs[km]':>7} {'I1_sph':>12} {'I1_pp':>12} {'dev[%]':>9}")
    for vza, sza, Hs in ((55, 0, 8.0), (55, 0, 2.0), (55, 0, 0.5), (55, 0, 0.1),
                         (40, 0, 8.0), (70, 0, 8.0), (55, 80, 8.0), (41.7788, 0, 8.0)):
        a = I1_sphere(vza, sza, 90.0, tau, Hs)
        b = I1_pp(vza, sza, tau)
        print(f"{vza:5.1f} {sza:5.0f} {Hs:7.1f} {a:12.6e} {b:12.6e} {100*(b/a-1):+9.3f}")
    # 소각 해석 추정치: 질량이 고도 Hs 에 있을 때 경로인자 비 sec(theta(Hs))/sec(vza)
    th = np.degrees(np.arcsin(RE / (RE + 8.0) * np.sin(np.radians(55))))
    print(f"\n해석 추정 (질량이 8 km 한 층에 있다고 볼 때): sec(θ(8km))/sec55 − 1 = "
          f"{np.cos(np.radians(55))/np.cos(np.radians(th))-1:+.4%}   (θ(8km)={th:.3f}°)")
