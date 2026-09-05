"""OCRT v1.19 Python port — quadrature + angular kernels.

Faithful transliteration of:
  - src/rt_quadrature.c  : rt_quadrature_gauss_legendre_pos
  - src/rt_kernel.c      : rt_legendre_compute(_pol), rt_kernel_phase_fourier(_pol),
                           rt_kernel_phase_fourier_aerosol_full
Storage convention mirrors the C code: direction index j in [-n_mu, +n_mu]
is mapped to numpy index j + n_mu.  plm/rrl/rtl are (l_max+1, 2n+1) arrays;
phase_fourier_m / gr / gt / arr / art / att are (n_mu+1, 2n+1) arrays.
"""
import numpy as np
from math import sqrt, cos, pi, fabs, pow as fpow


def gauss_legendre_pos(n_mu):
    """rt_quadrature_gauss_legendre_pos: positive half of GL(2*n_mu) on [-1,1].
    Newton iteration with Tricomi initial guess, exactly as the C code."""
    N = 2 * n_mu
    tol = 1.0e-15
    mu = np.empty(n_mu)
    w = np.empty(n_mu)

    def legendre_eval(t):
        p0, p1 = 1.0, t
        for n in range(1, N):
            p2 = ((2.0 * n + 1.0) * t * p1 - n * p0) / (n + 1)
            p0, p1 = p1, p2
        return p1, p0  # P_N, P_{N-1}

    for k in range(1, n_mu + 1):
        t = cos(pi * (k - 0.25) / (N + 0.5))
        converged = False
        for _ in range(100):
            pn, pn_m1 = legendre_eval(t)
            dp = N * (t * pn - pn_m1) / (t * t - 1.0)
            dt = -pn / dp
            t += dt
            if fabs(dt) < tol * (fabs(t) + 1.0):
                converged = True
                break
        if not converged:
            raise RuntimeError("GL Newton did not converge")
        pn, pn_m1 = legendre_eval(t)
        dp = N * (t * pn - pn_m1) / (t * t - 1.0)
        w_full = 2.0 / ((1.0 - t * t) * dp * dp)
        idx = n_mu - k
        mu[idx] = t
        w[idx] = w_full
    return mu, w


class LegendreWorkspace:
    """rt_legendre_workspace_t equivalent.  Index helper: dir j -> j+n_mu."""

    def __init__(self, n_mu, l_max):
        self.n_mu = n_mu
        self.l_max = l_max
        dirs = 2 * n_mu + 1
        self.plm = np.zeros((l_max + 1, dirs))
        self.rrl = np.zeros((l_max + 1, dirs))
        self.rtl = np.zeros((l_max + 1, dirs))
        self.pfm = np.zeros((n_mu + 1, dirs))   # phase_fourier_m
        self.gr = np.zeros((n_mu + 1, dirs))
        self.gt = np.zeros((n_mu + 1, dirs))
        self.arr = np.zeros((n_mu + 1, dirs))
        self.art = np.zeros((n_mu + 1, dirs))
        self.att = np.zeros((n_mu + 1, dirs))


def legendre_compute(ws, rm, m, xpl_out=None):
    """rt_legendre_compute.  rm: signed mu array of length 2n+1 (index j+n_mu).
    Fills ws.plm; returns xpl (= plm[2] row copy) like atm->xpl."""
    n_mu, l_max = ws.n_mu, ws.l_max
    if l_max < m or l_max < 2:
        raise ValueError("l_max too small")
    plm = ws.plm
    c = rm  # full signed vector

    plm[:min(m, l_max + 1), :] = 0.0

    if m == 0:
        plm[0, :] = 1.0
        plm[1, :] = c
        plm[2, :] = 0.5 * (3.0 * c * c - 1.0)
    elif m == 1:
        sqrt3 = sqrt(3.0)
        x = 1.0 - c * c
        plm[1, :] = np.sqrt(0.5 * x)
        plm[2, :] = c * plm[1, :] * sqrt3
    else:
        a = 1.0
        for i in range(1, m + 1):
            a *= sqrt((i + m) / i) * 0.5
        xx = 1.0 - c * c
        plm[m - 1, :] = 0.0
        plm[m, :] = a * np.power(xx, 0.5 * m)

    l_start = 2 if m < 2 else m
    for l in range(l_start, l_max):
        a_rec = (2.0 * l + 1.0) / sqrt((l + m + 1) * (l - m + 1))
        b_rec = sqrt((l + m) * (l - m)) / (2.0 * l + 1.0)
        plm[l + 1, :] = a_rec * (c * plm[l, :] - b_rec * plm[l - 1, :])

    xpl = plm[2, :].copy()
    if xpl_out is not None:
        xpl_out[:] = xpl
    return xpl


def legendre_compute_pol(ws, rm, m):
    """rt_legendre_compute_pol: plm + rrl/rtl.  Returns (xpl, xrl, xtl)."""
    xpl = legendre_compute(ws, rm, m)
    n_mu, l_max = ws.n_mu, ws.l_max
    rrl, rtl = ws.rrl, ws.rtl
    c = rm

    lim = min(max(m, 2), l_max + 1)
    rrl[:lim, :] = 0.0
    rtl[:lim, :] = 0.0

    if m == 0:
        xx = 1.0 - c * c
        rrl[2, :] = 3.0 * xx / (2.0 * sqrt(6.0))
        rtl[2, :] = 0.0
    elif m == 1:
        x = 1.0 - c * c
        sx = np.sqrt(x)
        rrl[2, :] = -c * sx * 0.5
        rtl[2, :] = -sx * 0.5
    else:
        a = 1.0
        for i in range(1, m + 1):
            a *= sqrt((i + m) / i) * 0.5
        b = a * sqrt(m / (m + 1)) * sqrt((m - 1) / (m + 2))
        xx = 1.0 - c * c
        pxx = np.power(xx, 0.5 * m - 1.0)
        rrl[m, :] = b * (1.0 + c * c) * pxx
        rtl[m, :] = b * 2.0 * c * pxx

    l_start = 2 if m < 2 else m
    for l in range(l_start, l_max):
        d = (l + 1) * (2 * l + 1) / sqrt((l + 3) * (l - 1) * (l + m + 1) * (l - m + 1))
        e = sqrt((l + 2) * (l - 2) * (l + m) * (l - m)) / (l * (2 * l + 1))
        f = 2.0 * m / (l * (l + 1))
        rrl[l + 1, :] = d * (c * rrl[l, :] - f * rtl[l, :] - e * rrl[l - 1, :])
        rtl[l + 1, :] = d * (c * rtl[l, :] - f * rrl[l, :] - e * rtl[l - 1, :])

    xrl = rrl[2, :].copy()
    xtl = rtl[2, :].copy()
    return xpl, xrl, xtl


def kernel_phase_fourier(ws, m, betal):
    """rt_kernel_phase_fourier: pfm[j][k] = sum_l plm[l][j] plm[l][k] betal[l]."""
    n_mu, l_max = ws.n_mu, ws.l_max
    if betal is None:
        ws.pfm[:, :] = 0.0
        return
    b = np.asarray(betal[: l_max + 1])
    P = ws.plm[m: l_max + 1, :]          # (L, dirs)
    Pw = P * b[m:, None]                  # weight by betal
    # rows j = 0..n_mu correspond to dirs index j+n_mu
    ws.pfm[:, :] = P[:, n_mu:].T @ Pw    # (n_mu+1, dirs)


def kernel_phase_fourier_pol(ws, m, gammal):
    n_mu, l_max = ws.n_mu, ws.l_max
    if gammal is None:
        ws.gr[:, :] = 0.0
        ws.gt[:, :] = 0.0
        return
    g = np.asarray(gammal[: l_max + 1])
    P = ws.plm[m: l_max + 1, :]
    R = ws.rrl[m: l_max + 1, :]
    T = ws.rtl[m: l_max + 1, :]
    Pj = P[:, n_mu:]                      # (L, n_mu+1)
    ws.gr[:, :] = Pj.T @ (R * g[m:, None])
    ws.gt[:, :] = Pj.T @ (T * g[m:, None])


def kernel_phase_fourier_aerosol_full(ws, m, alphal, zetal):
    n_mu, l_max = ws.n_mu, ws.l_max
    if alphal is None or zetal is None:
        ws.arr[:, :] = 0.0
        ws.art[:, :] = 0.0
        ws.att[:, :] = 0.0
        return
    a = np.asarray(alphal[: l_max + 1])[m:]
    z = np.asarray(zetal[: l_max + 1])[m:]
    R = ws.rrl[m: l_max + 1, :]
    T = ws.rtl[m: l_max + 1, :]
    Rj, Tj = R[:, n_mu:], T[:, n_mu:]
    # att = tj*tk*a + rj*rk*z ; arr = tj*tk*z + rj*rk*a ; art = tj*rk*a + rj*tk*z
    ws.att[:, :] = Tj.T @ (T * a[:, None]) + Rj.T @ (R * z[:, None])
    ws.arr[:, :] = Tj.T @ (T * z[:, None]) + Rj.T @ (R * a[:, None])
    ws.art[:, :] = Tj.T @ (R * a[:, None]) + Rj.T @ (T * z[:, None])
