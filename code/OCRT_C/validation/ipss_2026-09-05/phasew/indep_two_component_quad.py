#!/usr/bin/env python3
"""독립 검산 v3 — 2성분(Rayleigh + 에어로졸) 위상함수 가중 κ.  rt_ipss.c 미사용.

κ = I1_sph / I1_pp,  I1 = ∫ [β_R(h) P_R(Θ) + ω_A β_A(h) P_A(Θ)] exp(−τ_view − τ_sun) ds
  · Rayleigh: 지수 8 km, τ_R = 0.0935, depol 0.0279
  · 에어로졸: 지수 2 km, τ_A = 0.1, ω_A = 0.95, Henyey-Greenstein g = 0.7
  · 소멸(감쇠)은 β_R + β_A (흡수 포함), 원천은 산란분만.
  · Θ 는 평행빔 + 직선 시선이라 경로상 일정: rt_ipss.c 와 같은 벡터 규약
    V0 = (sin θs, 0, −cos θs), u = (sin θv cos φp, sin θv sin φp, cos θv), φp = 180° − raa.
같은 프로파일에서 위상함수를 뺀(phase-free) κ 도 같이 내서 두 정의의 차이를 보인다.
scipy quad, 화소(0,0,Re) 앵커, 지표 차폐 포함."""
import numpy as np
from scipy.integrate import quad

RE, H = 6371.0, 100.0
TAU_R, HS_R, DEPOL = 0.0935, 8.0, 0.0279
TAU_A, HS_A, OM_A, G_HG = 0.1, 2.0, 0.95, 0.7


def bexp(h, tau, Hs):
    b0 = tau / (Hs * (1.0 - np.exp(-H / Hs)))
    return b0 * np.exp(-h / Hs)


def beta_R(h): return bexp(h, TAU_R, HS_R)
def beta_A(h): return bexp(h, TAU_A, HS_A)
def beta_ext(h): return beta_R(h) + beta_A(h)


def P_R(ct):
    b2 = (1.0 - DEPOL) / (2.0 + DEPOL)
    return 1.0 + b2 * 0.5 * (3.0 * ct * ct - 1.0)


def P_A(ct):
    g = G_HG
    return (1.0 - g * g) / (1.0 + g * g - 2.0 * g * ct) ** 1.5


def t_exit(p0, d):
    b = p0 @ d; c = p0 @ p0 - (RE + H) ** 2
    return -b + np.sqrt(b * b - c)


def hits_ground(p0, d, tmax):
    b = p0 @ d; c = p0 @ p0 - RE ** 2
    disc = b * b - c
    if disc <= 0: return False
    t1 = -b - np.sqrt(disc)
    return 1e-9 < t1 < tmax


def tau_along(p0, d):
    tm = t_exit(p0, d)
    if hits_ground(p0, d, tm): return np.inf
    f = lambda t: beta_ext(np.linalg.norm(p0 + t * d) - RE)
    pts = [x for x in (1.0, 2.0, 4.0, 8.0, 16.0, 40.0) if 0 < x < tm]
    v, _ = quad(f, 0.0, tm, points=pts, limit=500, epsabs=0, epsrel=1e-10)
    return v


def I1_sphere(vza, sza, raa, phase):
    tv, ts, fp = np.radians(vza), np.radians(sza), np.radians(180.0 - raa)
    u = np.array([np.sin(tv) * np.cos(fp), np.sin(tv) * np.sin(fp), np.cos(tv)])
    to_sun = np.array([-np.sin(ts), 0.0, np.cos(ts)])          # rt_ipss.c 규약
    ct = -(to_sun @ u)                                          # V0 · u
    PR, PA = (P_R(ct), P_A(ct)) if phase else (1.0, 1.0)
    P = np.array([0.0, 0.0, RE])
    s_max = t_exit(P, u)

    def tau_view(s):
        f = lambda x: beta_ext(np.linalg.norm(P + x * u) - RE)
        v, _ = quad(f, s, s_max, limit=500, epsabs=0, epsrel=1e-10)
        return v

    def integrand(s):
        r = P + s * u
        h = np.linalg.norm(r) - RE
        tsun = tau_along(r, to_sun)
        if not np.isfinite(tsun): return 0.0
        src = beta_R(h) * PR + OM_A * beta_A(h) * PA
        return src * np.exp(-tau_view(s) - tsun)

    mv = np.cos(tv)
    pts = [x / mv for x in (1.0, 2.0, 4.0, 8.0, 16.0, 40.0) if 0 < x / mv < s_max]
    v, _ = quad(integrand, 0.0, s_max, points=pts, limit=500, epsabs=0, epsrel=1e-9)
    return v, ct


def I1_pp(vza, sza, ct, phase):
    mv, ms = np.cos(np.radians(vza)), np.cos(np.radians(sza))
    C = 1.0 / mv + 1.0 / ms
    PR, PA = (P_R(ct), P_A(ct)) if phase else (1.0, 1.0)
    tau_above = lambda h: quad(beta_ext, h, H, limit=500, epsabs=0, epsrel=1e-10)[0]
    f = lambda h: (beta_R(h) * PR + OM_A * beta_A(h) * PA) * np.exp(-C * tau_above(h)) / mv
    pts = [1.0, 2.0, 4.0, 8.0, 16.0, 40.0]
    v, _ = quad(f, 0.0, H, points=pts, limit=500, epsabs=0, epsrel=1e-9)
    return v


if __name__ == "__main__":
    cases = [(55, 0, 90), (55, 40, 90), (55, 80, 90), (70, 80, 0), (70, 80, 180), (40, 60, 45), (0, 40, 0)]
    print("vza,sza,raa,cos_theta,P_R,P_A,kappa_phase,kappa_free")
    for vza, sza, raa in cases:
        a, ct = I1_sphere(vza, sza, raa, True)
        b = I1_pp(vza, sza, ct, True)
        a0, _ = I1_sphere(vza, sza, raa, False)
        b0 = I1_pp(vza, sza, ct, False)
        print(f"{vza},{sza},{raa},{ct:.9f},{P_R(ct):.6f},{P_A(ct):.6f},{a/b:.9f},{a0/b0:.9f}")
