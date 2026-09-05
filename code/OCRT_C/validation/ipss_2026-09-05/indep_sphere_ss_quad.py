#!/usr/bin/env python3
"""독립 검산 v2 — scipy 적응 구적(quad). 작은 척도높이에서도 수렴하도록 함.
rt_ipss.c 와 코드를 전혀 공유하지 않는다."""
import numpy as np
from scipy.integrate import quad

RE, H = 6371.0, 100.0


def beta_of_h(h, tau_t, Hs):
    if Hs > 0:
        b0 = tau_t / (Hs * (1.0 - np.exp(-H / Hs)))
        return b0 * np.exp(-h / Hs)
    return tau_t / H


def t_exit(p0, d):
    b = p0 @ d; c = p0 @ p0 - (RE + H) ** 2
    return -b + np.sqrt(b * b - c)


def hits_ground(p0, d, tmax):
    b = p0 @ d; c = p0 @ p0 - RE ** 2
    disc = b * b - c
    if disc <= 0: return False
    t1 = -b - np.sqrt(disc)
    return 1e-9 < t1 < tmax


def tau_along(p0, d, tau_t, Hs):
    tm = t_exit(p0, d)
    if hits_ground(p0, d, tm): return np.inf
    f = lambda t: beta_of_h(np.linalg.norm(p0 + t * d) - RE, tau_t, Hs)
    pts = [x for x in (Hs * 0.5, Hs, 3 * Hs, 10 * Hs) if 0 < x < tm] if Hs > 0 else None
    v, _ = quad(f, 0.0, tm, points=pts, limit=400, epsabs=0, epsrel=1e-10)
    return v


def I1_sphere(vza, sza, raa, tau_t, Hs):
    tv, ts, fp = np.radians(vza), np.radians(sza), np.radians(180.0 - raa)
    u = np.array([np.sin(tv) * np.cos(fp), np.sin(tv) * np.sin(fp), np.cos(tv)])
    to_sun = np.array([np.sin(ts), 0.0, np.cos(ts)])
    P = np.array([0.0, 0.0, RE])
    s_max = t_exit(P, u)

    def tau_view(s):   # 산란점 s 에서 TOA 까지 시선 감쇠
        f = lambda x: beta_of_h(np.linalg.norm(P + x * u) - RE, tau_t, Hs)
        v, _ = quad(f, s, s_max, limit=400, epsabs=0, epsrel=1e-10)
        return v

    def integrand(s):
        r = P + s * u
        h = np.linalg.norm(r) - RE
        ts_ = tau_along(r, to_sun, tau_t, Hs)
        if not np.isfinite(ts_): return 0.0
        return beta_of_h(h, tau_t, Hs) * np.exp(-tau_view(s) - ts_)

    mv = np.cos(tv)
    pts = [x / mv for x in (Hs * 0.5, Hs, 3 * Hs, 10 * Hs, 30 * Hs) if 0 < x / mv < s_max] if Hs > 0 else None
    v, _ = quad(integrand, 0.0, s_max, points=pts, limit=400, epsabs=0, epsrel=1e-9)
    return v


def I1_pp(vza, sza, tau_t):
    mv, ms = np.cos(np.radians(vza)), np.cos(np.radians(sza))
    C = 1.0 / mv + 1.0 / ms
    return (1.0 / (mv * C)) * (1.0 - np.exp(-C * tau_t))


if __name__ == "__main__":
    tau = 0.0935
    print("독립 Python 적응구적  —  dev = pp/sph − 1 [%]  (rt_ipss.c 미사용)")
    print(f"{'VZA':>7} {'SZA':>4} {'Hs[km]':>7} {'dev[%]':>9}   rt_ipss.c 값")
    ref = {(55, 0, 8.0): "+0.240", (55, 0, 2.0): "+0.060", (55, 0, 0.5): "+0.015",
           (55, 0, 0.1): "+0.004", (40, 0, 8.0): "+0.085", (70, 0, 8.0): "+0.830",
           (55, 80, 8.0): "-0.591", (41.7788, 0, 8.0): "+0.096", (55, 0, 0.0): "+1.317"}
    for vza, sza, Hs in ((55, 0, 8.0), (55, 0, 2.0), (55, 0, 0.5), (55, 0, 0.1),
                         (40, 0, 8.0), (70, 0, 8.0), (55, 80, 8.0), (41.7788, 0, 8.0),
                         (55, 0, 0.0)):
        a = I1_sphere(vza, sza, 90.0, tau, Hs)
        b = I1_pp(vza, sza, tau)
        print(f"{vza:7.4f} {sza:4.0f} {Hs:7.1f} {100*(b/a-1):+9.3f}   {ref.get((vza,sza,Hs),'')}")
