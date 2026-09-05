# -*- coding: utf-8 -*-
"""Aerosol mie file reader and Greek-coefficient decomposition (Phase C).

read_mie: OPAC-layout .mie parser (mie_io.c read_mie_file).
  Line 1: n_ang.  Spectral table (n_wl rows x 6 cols after wavelength).
  Phase blocks P11, P12, P33 (each n_ang rows of angle + n_phase_wl values).
"""
import os
import warnings
from dataclasses import dataclass

import numpy as np

from .spectral_contract import require_wavelength, require_query_in_table
from .backend import xp


def _floats(line):
    out = []
    for tok in line.split():
        try:
            out.append(float(tok))
        except ValueError:
            return out if out else []  # first non-float ends the numeric run only if we have some
    return out


def _all_floats(line):
    """Parse a line as floats; return [] if ANY token is non-numeric
    (matches C parse_floats which reads a fixed column count)."""
    toks = line.split()
    if not toks:
        return []
    out = []
    for tok in toks:
        try:
            out.append(float(tok))
        except ValueError:
            return []
    return out


class MieData:
    __slots__ = ("n_ang", "n_wl", "n_phase_wl", "wavelengths",
                 "phase_wavelengths", "angles", "spectral", "P11", "P12", "P33",
                 "source_path", "source_mtime_ns", "source_size", "cache_key")


# MIE ANGLE-GRID / TRUNCATION POLICY
#
# Exact FR631 is the canonical raw-phase scientific representation and broad
# truncation is OFF by default. A tested 0--0.005 degree local cap was
# negligible, but OSOAA-style broad truncation was not and can change
# directional Rrs/rrs by several percent. Non-FR631 is not synonymous with
# coarse: recognize only exact FR631 and the known legacy 361 x 0.5-degree
# grid; treat every other custom grid as unvalidated and never auto-enable
# truncation.
#
# DOC-REF:
#   docs/OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21.md
#   validation/artifacts/README_KO.md
_MIE_GRID_WARNED_KEYS = set()


def _fr631_expected_descending():
    return np.concatenate([
        np.arange(180.0, 19.999999, -0.5),
        np.arange(19.9, 4.999999, -0.1),
        np.arange(4.95, 0.999999, -0.05),
        np.arange(0.98, 0.199999, -0.02),
        np.arange(0.195, -0.000001, -0.005),
    ])


def classify_mie_angle_grid(angles, atol=1.0e-8):
    """Return ``FR631``, ``legacy-361-uniform`` or
    ``noncanonical-unvalidated``. Classification is descriptive only and never
    changes truncation or any other physics option."""
    a = np.asarray(angles, float)
    fr = _fr631_expected_descending()
    if len(a) == 631 and (
        np.allclose(a, fr, rtol=0.0, atol=atol)
        or np.allclose(a, fr[::-1], rtol=0.0, atol=atol)
    ):
        return "FR631"
    legacy = np.arange(180.0, -0.000001, -0.5)
    if len(a) == 361 and (
        np.allclose(a, legacy, rtol=0.0, atol=atol)
        or np.allclose(a, legacy[::-1], rtol=0.0, atol=atol)
    ):
        return "legacy-361-uniform"
    return "noncanonical-unvalidated"


def _warn_mie_angle_grid_once(source_path, cache_key, angles):
    cls = classify_mie_angle_grid(angles)
    if cls == "FR631" or cache_key in _MIE_GRID_WARNED_KEYS:
        return
    _MIE_GRID_WARNED_KEYS.add(cache_key)
    a = np.asarray(angles, float)
    positive = a[a > 0.0]
    min_positive = float(positive.min()) if positive.size else float("nan")
    gaps = np.abs(np.diff(a))
    forward_pairs = np.minimum(a[:-1], a[1:]) <= 20.0
    max_forward_gap = float(gaps[forward_pairs].max()) if np.any(forward_pairs) else float("nan")
    if cls == "legacy-361-uniform":
        message = (
            f"legacy 361-angle 0.5-degree Mie grid detected in '{source_path}'. "
            "The forward peak may be under-resolved. Replace it with canonical "
            "FR631 when possible. Broad truncation is an advanced fallback that "
            "changes effective phase, b, omega, tau and source function; it is "
            "never enabled automatically."
        )
    else:
        message = (
            f"non-canonical Mie angle grid detected in '{source_path}' "
            f"(n={len(a)}, min-positive={min_positive:.9g} deg, "
            f"max-gap<=20deg={max_forward_gap:.9g} deg). Forward adequacy is "
            "unvalidated. Do not assume every non-FR631 grid is coarse; inspect "
            "the grid and transport convergence or regenerate FR631. Truncation "
            "is not enabled automatically."
        )
    warnings.warn(message, RuntimeWarning, stacklevel=2)


def read_mie(path):
    """Return MieData. wavelengths in um. spectral columns:
    0 Nor_Ext_Co, 1 Nor_Sca_Co, 2 SSA, 3 Asymm, 4 Extinct_Co, 5 Scatter_Co."""
    source_path = os.path.realpath(os.fspath(path))
    st = os.stat(source_path)
    with open(source_path, encoding="utf-8", errors="ignore") as fh:
        raw = fh.readlines()

    # Line 1: n_ang (first line whose first token is numeric)
    i = 0
    n_ang = 0
    while i < len(raw):
        s = raw[i].strip()
        if s:
            tok0 = s.split()[0]
            try:
                n_ang = int(float(tok0))
                i += 1
                break
            except ValueError:
                pass
        i += 1
    if n_ang <= 0:
        raise ValueError("mie: n_ang not found")

    # Spectral table: rows with exactly 7 floats, until "Phase Function".
    wl = []
    spec = []
    while i < len(raw):
        if "Phase Function" in raw[i]:
            break
        vals = _all_floats(raw[i])
        if len(vals) == 7:
            wl.append(vals[0])
            spec.append(vals[1:7])
        i += 1
    n_wl = len(wl)
    if n_wl < 1:
        raise ValueError("mie: no spectral rows")

    # Phase wavelength grid from the TETA header (line after "Phase Function (P11)").
    # The header line is "TETA <wl1> <wl2> ...": non-numeric first token.
    phase_wl = None
    j = i
    while j < len(raw):
        if "TETA" in raw[j]:
            toks = raw[j].split()
            cand = []
            for t in toks:
                try:
                    cand.append(float(t))
                except ValueError:
                    continue
            if cand:
                phase_wl = cand
            break
        j += 1
    if phase_wl is None:
        phase_wl = list(wl)
    n_phase_wl = len(phase_wl)

    # Phase blocks: skip to first row with (n_phase_wl+1) floats, read n_ang rows.
    blocks = []
    angles = None
    while i < len(raw) and len(blocks) < 3:
        while i < len(raw):
            vals = _all_floats(raw[i])
            if len(vals) == n_phase_wl + 1:
                break
            i += 1
        if i >= len(raw):
            break
        blk = np.zeros((n_ang, n_phase_wl))
        angs = np.zeros(n_ang)
        ok = 0
        for k in range(n_ang):
            if i + k >= len(raw):
                break
            vals = _all_floats(raw[i + k])
            if len(vals) < n_phase_wl + 1:
                break
            angs[k] = vals[0]
            blk[k, :] = vals[1:n_phase_wl + 1]
            ok += 1
        if ok != n_ang:
            break
        if angles is None:
            angles = angs
        blocks.append(blk)
        i += n_ang
    if len(blocks) < 3:
        raise ValueError("mie: expected 3 phase blocks, found %d" % len(blocks))

    m = MieData()
    m.n_ang = n_ang
    m.n_wl = n_wl
    m.n_phase_wl = n_phase_wl
    m.wavelengths = np.array(wl)
    m.phase_wavelengths = np.array(phase_wl)
    m.angles = angles
    m.spectral = np.array(spec)            # (n_wl, 6)
    m.P11, m.P12, m.P33 = blocks           # each (n_ang, n_phase_wl)
    m.source_path = source_path
    m.source_mtime_ns = int(st.st_mtime_ns)
    m.source_size = int(st.st_size)
    m.cache_key = (source_path, m.source_mtime_ns, m.source_size)
    _warn_mie_angle_grid_once(source_path, m.cache_key, m.angles)
    return m


# ---- PCHIP (Fritsch-Carlson 1980, matches shared/numerics.c) ----
class Pchip:
    def __init__(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        n = len(x)
        if n < 2:
            raise ValueError("PCHIP requires at least two points")
        h = np.diff(x)
        if np.any(h <= 0.0):
            raise ValueError("PCHIP x grid must be strictly increasing")
        d = np.diff(y) / h
        m = np.zeros(n)
        if n == 2:
            m[:] = d[0]
            self.x = x; self.y = y; self.m = m; self.n = n
            return
        for k in range(1, n - 1):
            if d[k - 1] * d[k] <= 0.0:
                m[k] = 0.0
            else:
                w1 = 2.0 * h[k] + h[k - 1]
                w2 = h[k] + 2.0 * h[k - 1]
                m[k] = (w1 + w2) / (w1 / d[k - 1] + w2 / d[k])
        m[0] = ((2.0 * h[0] + h[1]) * d[0] - h[0] * d[1]) / (h[0] + h[1])
        if m[0] * d[0] <= 0.0:
            m[0] = 0.0
        elif (d[0] * d[1] <= 0.0) and (abs(m[0]) > 3.0 * abs(d[0])):
            m[0] = 3.0 * d[0]
        nn = n - 1
        m[nn] = ((2.0 * h[nn - 1] + h[nn - 2]) * d[nn - 1] - h[nn - 1] * d[nn - 2]) \
            / (h[nn - 1] + h[nn - 2])
        if m[nn] * d[nn - 1] <= 0.0:
            m[nn] = 0.0
        elif (d[nn - 1] * d[nn - 2] <= 0.0) and (abs(m[nn]) > 3.0 * abs(d[nn - 1])):
            m[nn] = 3.0 * d[nn - 1]
        self.x = x; self.y = y; self.m = m; self.n = n

    def eval_vec(self, xq):
        xq = np.asarray(xq, float)
        x = self.x; n = self.n
        lo = np.searchsorted(x, xq, side='right') - 1
        lo = np.clip(lo, 0, n - 2)
        hi = lo + 1
        h = x[hi] - x[lo]
        t = (xq - x[lo]) / h
        t2 = t * t; t3 = t2 * t
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2
        out = h00 * self.y[lo] + h10 * h * self.m[lo] + \
            h01 * self.y[hi] + h11 * h * self.m[hi]
        # shared/numerics.c::pchip_interp clamps outside the tabulated range.
        out = np.where(xq <= x[0], self.y[0], out)
        out = np.where(xq >= x[-1], self.y[-1], out)
        return out


def _ascending(theta_deg, P):
    """OPAC files are 180->0; return theta ascending + P reordered."""
    t = np.asarray(theta_deg, float)
    p = np.asarray(P, float)
    if t[0] > t[-1]:
        t = t[::-1].copy(); p = p[::-1].copy()
    order = np.argsort(t, kind='stable')
    return t[order], p[order]


def compute_legendre_moments(P, theta_deg, n_max, fine=4096):
    """chi_l = (2l+1)/2 * int P(mu) P_l(mu) dmu (rt_aerosol_compute_legendre_moments).
    PCHIP in theta, dense grid theta=180*(fine-1-i)/(fine-1) (mu ascending),
    adjacent-trapezoid in mu, Bonnet recurrence."""
    ts, Ps = _ascending(theta_deg, P)
    sp = Pchip(ts, Ps)
    i = np.arange(fine)
    theta_i = 180.0 * (fine - 1 - i) / (fine - 1)
    mu = np.cos(theta_i * np.pi / 180.0)
    Pf = sp.eval_vec(theta_i)
    chi = np.zeros(n_max)
    Plm1 = np.zeros(fine); Pl = np.ones(fine)
    dmu = mu[1:] - mu[:-1]
    for l in range(n_max):
        integ = Pf * Pl
        s = np.sum(0.5 * dmu * (integ[1:] + integ[:-1]))
        chi[l] = (2.0 * l + 1.0) * 0.5 * s
        if l + 1 < n_max:
            nxt = ((2.0 * l + 1.0) * mu * Pl - l * Plm1) / (l + 1.0)
            Plm1 = Pl; Pl = nxt
    return chi


def _p2_matrix(mu, L_max):
    """legendre_P2_recurrence for all mu: returns (len(mu), L_max+1)."""
    fine = len(mu)
    pol = np.zeros((fine, L_max + 1))
    if L_max >= 2:
        pol[:, 2] = 3.0 * (1.0 - mu * mu) / (2.0 * np.sqrt(6.0))
        for k in range(2, L_max):
            d = (2.0 * k + 1.0) / np.sqrt((k + 3.0) * (k - 1.0))
            e = np.sqrt((k + 2.0) * (k - 2.0)) / (2.0 * k + 1.0)
            pol[:, k + 1] = d * (mu * pol[:, k] - e * pol[:, k - 1])
    return pol


def compute_p2_moments(P12, theta_deg, L_max, fine=4096):
    """gamma_l = (2l+1)/2 * int P12(mu) P2_l(mu) dmu (compute_p2_moments).
    Central-difference trapezoid in mu."""
    ts, Ps = _ascending(theta_deg, P12)
    sp = Pchip(ts, Ps)
    dtheta = (ts[-1] - ts[0]) / (fine - 1)
    i = np.arange(fine)
    th = ts[0] + i * dtheta
    Pth = sp.eval_vec(th)
    mu_i = np.cos(th * np.pi / 180.0)      # descending as th ascends
    mu_grid = mu_i[::-1].copy()            # ascending
    integ = Pth[::-1].copy()
    P2 = _p2_matrix(mu_grid, L_max)        # (fine, L_max+1)
    dmu = np.empty(fine)
    dmu[0] = (mu_grid[1] - mu_grid[0]) * 0.5
    dmu[-1] = (mu_grid[-1] - mu_grid[-2]) * 0.5
    dmu[1:-1] = (mu_grid[2:] - mu_grid[:-2]) * 0.5
    gam = (integ[:, None] * P2 * dmu[:, None]).sum(axis=0)
    gam *= (2.0 * np.arange(L_max + 1) + 1.0) * 0.5
    return gam


def compute_vector_legendre(P11, P12, P33, theta_deg, L_max, fine=4096):
    """Return betal,gammal,alphal,zetal (rt_aerosol_compute_vector_legendre).
    P11/P12/P33 are the phase values at one wavelength; theta_deg the angle grid."""
    betal = compute_legendre_moments(P11, theta_deg, L_max + 1, fine)
    gammal = compute_p2_moments(P12, theta_deg, L_max, fine)
    deltal = compute_legendre_moments(P33, theta_deg, L_max + 1, fine)
    alphal = np.zeros(L_max + 1)
    zetal = np.zeros(L_max + 1)
    for ii in range(2, L_max + 1):
        co1 = 4.0 * (2.0 * ii + 1.0) / (ii * (ii - 1.0) * (ii + 1.0) * (ii + 2.0))
        co2 = ii * (ii - 1.0) / ((ii + 1.0) * (ii + 2.0))
        co3 = co2 * deltal[ii]
        co2_b = co2 * betal[ii]
        nn = ii // 2
        mm = (ii - 1) // 2
        som1 = som2 = som3 = som4 = 0.0
        for j in range(1, nn + 1):
            c2 = (ii - 1.0) * (ii - 1.0) - 3.0 * (2 * j - 1) * (ii - j)
            idx = ii - 2 * j
            som1 += c2 * betal[idx]
            som2 += c2 * deltal[idx]
        for j in range(0, mm + 1):
            c2 = (ii - 1.0) * (ii - 1.0) - 3.0 * j * (2 * ii - 2 * j - 1)
            idx = ii - 2 * j - 1
            som3 += c2 * betal[idx]
            som4 += c2 * deltal[idx]
        zetal[ii] = co3 - co1 * (som2 - som3)
        alphal[ii] = co2_b - co1 * (som1 - som4)
    b0 = betal[0]
    if b0 != 0.0:
        alphal = alphal / b0
        betal = betal / b0
        gammal = gammal / b0
        zetal = zetal / b0
    return betal, gammal, alphal, zetal


def loglin_truncate(P11, P12, P33, theta_deg, mu1=0.8, mu2=0.94, thresh=0.1):
    """rt_aerosol_loglin_truncate: real-angle log10-linear forward-peak
    truncation.  Returns (P11_out, P12_out, P33_out, A) in the ORIGINAL angle
    ordering.  A is the truncated forward fraction (0 if no truncation)."""
    from math import acos, log10, pi as PI
    t = np.asarray(theta_deg, float)
    P11a = np.asarray(P11, float)
    P12a = np.asarray(P12, float)
    P33a = np.asarray(P33, float)
    N = len(t)
    ascending = bool(t[0] <= t[-1])
    if ascending:
        order = np.arange(N)
    else:
        order = np.arange(N)[::-1]
    ts = t[order]
    P11s = P11a[order]

    theta1 = acos(mu1) * 180.0 / PI
    theta2 = acos(mu2) * 180.0 / PI
    K_hi = 0
    while K_hi < N and ts[K_hi] < theta1:
        K_hi += 1
    KK_hi = 0
    while KK_hi < N and ts[KK_hi] < theta2:
        KK_hi += 1

    def _no_trunc():
        return P11a.copy(), P12a.copy(), P33a.copy(), 0.0

    if K_hi >= N or K_hi == 0 or KK_hi >= N or KK_hi == 0 or K_hi <= KK_hi + 1:
        return _no_trunc()
    K_lo = K_hi - 1
    KK_lo = KK_hi - 1
    tK = (theta1 - ts[K_lo]) / (ts[K_hi] - ts[K_lo])
    tKK = (theta2 - ts[KK_lo]) / (ts[KK_hi] - ts[KK_lo])
    P11_at_t1 = (1.0 - tK) * P11s[K_lo] + tK * P11s[K_hi]
    P11_at_t2 = (1.0 - tKK) * P11s[KK_lo] + tKK * P11s[KK_hi]
    KK = KK_hi
    if P11_at_t1 <= 0.0 or P11_at_t2 <= 0.0:
        return _no_trunc()

    t1r = theta1 * PI / 180.0
    t2r = theta2 * PI / 180.0
    AA = (log10(P11_at_t2) - log10(P11_at_t1)) / (t2r - t1r)
    X1 = log10(P11_at_t2)
    P11_trunc_s = P11s.copy()
    for i in range(KK):
        tir = ts[i] * PI / 180.0
        P11_trunc_s[i] = 10.0 ** (X1 + AA * (tir - t2r))

    chi0 = compute_legendre_moments(P11_trunc_s, ts, 1, 3601)[0]
    A_trunc = 2.0 * (1.0 - chi0)
    if A_trunc < thresh:
        return _no_trunc()

    P11_out = np.empty(N)
    P12_out = np.empty(N)
    P33_out = np.empty(N)
    for i in range(N):
        oi = order[i]
        ratio = (P11_trunc_s[i] / P11s[i]) if P11s[i] > 0.0 else 1.0
        P11_out[oi] = P11_trunc_s[i]
        P12_out[oi] = P12a[oi] * ratio
        P33_out[oi] = P33a[oi] * ratio
    return P11_out, P12_out, P33_out, A_trunc


def _aer_p11_at_costheta(th_asc, P_asc, cth):
    """Linear interp of P11 at theta=acos(cth) on ascending theta grid (deg).
    Matches aer_phase_p11_at_costheta (endpoint clamp).  Backend-agnostic
    (xp) so it runs on the GPU when cth is a device array."""
    th_asc = xp.asarray(th_asc); P_asc = xp.asarray(P_asc)
    cth = xp.clip(xp.asarray(cth), -1.0, 1.0)
    theta = xp.arccos(cth) * 180.0 / xp.pi
    idx = xp.searchsorted(th_asc, theta, side='right') - 1
    idx = xp.clip(idx, 0, len(th_asc) - 2)
    lo = th_asc[idx]
    hi = th_asc[idx + 1]
    frac = (theta - lo) / (hi - lo)
    out = P_asc[idx] + frac * (P_asc[idx + 1] - P_asc[idx])
    out = xp.where(theta <= th_asc[0], P_asc[0], out)
    out = xp.where(theta >= th_asc[-1], P_asc[-1], out)
    return out


def aerosol_value_pfm(atm, m, th_asc, P_asc, nphi=720):
    """aerosol_value_phase_fourier: pfm[j][k] = (1/nphi) sum_q P11(cth) cos(m phi),
    cth = mu_j mu_k + sj sk cos(phi).  th_asc ascending (0..180 deg), P_asc the
    NORMALIZED P11 (divided by beta0).  Returns (n_mu+1, dirs), row 0 = solar."""
    n_mu = atm.n_mu
    dirs = 2 * n_mu + 1
    mu_j = atm.rm[n_mu:2 * n_mu + 1]          # rm[0..n_mu] in C offset = solar+positive
    mu_k = atm.rm[0:2 * n_mu + 1]             # rm[-n_mu..n_mu]
    sj = xp.sqrt(xp.maximum(0.0, 1.0 - mu_j * mu_j))
    sk = xp.sqrt(xp.maximum(0.0, 1.0 - mu_k * mu_k))
    phi = 2.0 * xp.pi * (xp.arange(nphi) + 0.5) / nphi
    cosmphi = xp.cos(m * phi)
    cth = (mu_j[:, None, None] * mu_k[None, :, None]
           + sj[:, None, None] * sk[None, :, None] * xp.cos(phi)[None, None, :])
    P = _aer_p11_at_costheta(th_asc, P_asc, cth)     # (n_mu+1, dirs, nphi)
    pfm = xp.einsum('jkp,p->jk', P, cosmphi) / nphi
    return pfm


def _aer_p11_at_costheta_batch(th_asc, P_asc, cth):
    """Batched _aer_p11_at_costheta.  th_asc (N,) shared ascending grid; P_asc
    (B, N) per-case P11; cth (B, ...) arbitrary trailing shape.  Returns cth's
    shape.  Per-case gather along the angle axis via take_along_axis; bit-
    identical to the per-case linear interp with endpoint clamp."""
    th_asc = xp.asarray(th_asc); P_asc = xp.asarray(P_asc)
    cth = xp.clip(xp.asarray(cth), -1.0, 1.0)
    theta = xp.arccos(cth) * 180.0 / xp.pi                 # (B, ...)
    B = P_asc.shape[0]
    N = th_asc.shape[0]
    idx = xp.searchsorted(th_asc, theta, side='right') - 1
    idx = xp.clip(idx, 0, N - 2)                            # (B, ...)
    lo = th_asc[idx]; hi = th_asc[idx + 1]                 # 1-D grid indexed by nD idx
    frac = (theta - lo) / (hi - lo)
    idx_flat = idx.reshape(B, -1)                          # (B, T)
    P0 = xp.take_along_axis(P_asc, idx_flat, axis=1).reshape(idx.shape)
    P1 = xp.take_along_axis(P_asc, idx_flat + 1, axis=1).reshape(idx.shape)
    out = P0 + frac * (P1 - P0)
    bshape = (B,) + (1,) * (theta.ndim - 1)
    p0 = P_asc[:, 0].reshape(bshape)
    pN = P_asc[:, -1].reshape(bshape)
    out = xp.where(theta <= th_asc[0], p0, out)
    out = xp.where(theta >= th_asc[-1], pN, out)
    return out


def aerosol_value_pfm_batch(bA, m, th_asc, P_asc, nphi=720):
    """Batched aerosol_value_phase_fourier over the case axis.  bA is a BatchAtm
    (ring = bA.rm (B, dirs)); th_asc (N,) shared; P_asc (B, N) per-case NORMALIZED
    P11.  Returns pfm (B, n_mu+1, dirs).  Bit-identical to per-case
    aerosol_value_pfm."""
    n_mu = bA.n_mu
    dirs = 2 * n_mu + 1
    rm = xp.asarray(bA.rm)                                  # (B, dirs)
    mu_j = rm[:, n_mu:2 * n_mu + 1]                         # (B, n_mu+1)
    mu_k = rm[:, 0:2 * n_mu + 1]                            # (B, dirs)
    sj = xp.sqrt(xp.maximum(0.0, 1.0 - mu_j * mu_j))
    sk = xp.sqrt(xp.maximum(0.0, 1.0 - mu_k * mu_k))
    phi = 2.0 * xp.pi * (xp.arange(nphi) + 0.5) / nphi
    cosmphi = xp.cos(m * phi)
    cth = (mu_j[:, :, None, None] * mu_k[:, None, :, None]
           + sj[:, :, None, None] * sk[:, None, :, None]
           * xp.cos(phi)[None, None, None, :])              # (B, n_mu+1, dirs, nphi)
    P = _aer_p11_at_costheta_batch(th_asc, P_asc, cth)      # (B, n_mu+1, dirs, nphi)
    pfm = xp.einsum('bjkp,p->bjk', P, cosmphi) / nphi       # (B, n_mu+1, dirs)
    return pfm


def aod_at_wavelength(mie, wl_nm, user_aod_ref, aod_ref_nm=555.0):
    """AODREF (v1.2): AOD(wl) = AOD(ref) * NorExt(wl)/NorExt(ref).
    NorExt = .mie spectral column 0, linear interp over wavelengths (um).
    At wl == ref the ratio is exactly 1 (C rt_aerosol_runtime.c:38-50)."""
    wl_nm = require_wavelength(wl_nm, context='aerosol AOD target')
    aod_ref_nm = require_wavelength(aod_ref_nm, context='aerosol AOD reference')
    wl_um = wl_nm / 1000.0
    ref_um = aod_ref_nm / 1000.0
    require_query_in_table(wl_um, mie.wavelengths, context='aerosol bulk grid', unit='um')
    require_query_in_table(ref_um, mie.wavelengths, context='aerosol bulk grid', unit='um')
    ext_t = float(np.interp(wl_um, mie.wavelengths, mie.spectral[:, 0]))
    ext_r = float(np.interp(ref_um, mie.wavelengths, mie.spectral[:, 0]))
    if not (ext_t > 0.0 and ext_r > 0.0):
        raise ValueError('non-positive NorExt in AODREF scaling')
    return user_aod_ref * (ext_t / ext_r)


def phase_matrix_at_wavelength_linear(mie, wl_nm):
    """Interpolate P11/P12/P33 in wavelength with the FR631 consumption
    contract: piecewise-LINEAR over the delivered 1-nm phase wavelength grid.
    Integer-nm queries resolve to exact column lookup (grid nodes), matching
    the C reader (mie_io.c phase_linear_eval / rt_aerosol_runtime prepare).
    DOC-REF: OCRT_330_1100_FR631_INTEGRATION_GUIDE_KO.md §3.2-3.3.
    """
    wl_nm = require_wavelength(wl_nm, context='aerosol phase')
    wl_um = float(wl_nm) / 1000.0
    x = np.asarray(mie.phase_wavelengths, float)
    require_query_in_table(wl_um, x, context='aerosol phase grid', unit='um')
    j = int(np.searchsorted(x, wl_um))
    if j <= 0:
        w0, j0, j1 = 1.0, 0, 0
    elif j >= x.size:
        w0, j0, j1 = 1.0, x.size - 1, x.size - 1
    else:
        j0, j1 = j - 1, j
        d = x[j1] - x[j0]
        w0 = (x[j1] - wl_um) / d if d > 0 else 1.0
    if abs(x[min(j1, x.size - 1)] - wl_um) < 1e-12:
        w0, j0 = 0.0, j1
    out = []
    for block in (mie.P11, mie.P12, mie.P33):
        B = np.asarray(block, float)
        out.append(w0 * B[:, j0] + (1.0 - w0) * B[:, j1])
    return tuple(out)


def phase_matrix_at_wavelength_pchip(mie, wl_nm):
    """Legacy-named compatibility wrapper. The FR631 contract fixes phase
    wavelength interpolation to LINEAR; this delegates to
    phase_matrix_at_wavelength_linear (same pattern as C mie_io.c
    mie_phase_nodes_at_wavelength_pchip -> _linear). 2026-08-23 cleanup.
    DOC-REF: OCRT_MIGRATION_STATUS §6(a).
    """
    return phase_matrix_at_wavelength_linear(mie, wl_nm)


def _readonly(a):
    out = np.asarray(a, float).copy()
    out.setflags(write=False)
    return out


@dataclass(frozen=True)
class AerosolRuntime:
    """Immutable aerosol object shared by every geometry in a full-grid run.

    It is the Python counterpart of C ``rt_aerosol_input_t`` plus runtime
    diagnostics.  Geometry-independent work (Mie spectral interpolation,
    forward-peak truncation, optical-depth/SSA transform and vector Legendre
    construction) is completed exactly once before the VZA/RAA loop.
    """
    source_key: tuple
    wavelength_nm: float
    aod_ref: float
    aod_ref_nm: float
    aod_target: float
    extinction_ratio: float
    ssa_raw: float
    tau_a_eff: float
    ssa_a_eff: float
    truncation_A: float
    L_max: int
    theta_asc: np.ndarray
    P11_norm_asc: np.ndarray
    P12_norm_asc: np.ndarray
    P33_norm_asc: np.ndarray
    betal: np.ndarray
    gammal: np.ndarray
    alphal: np.ndarray
    zetal: np.ndarray
    interpolation: str = "pchip-wavelength"
    truncation: str = "loglin"


def prepare_aerosol_runtime(mie, wavelength_nm, user_aod_ref,
                             aod_ref_nm=865.0, L_max=80,
                             loglin_mu1=0.8, loglin_mu2=0.94,
                             loglin_threshold=0.1):
    """Prepare one immutable aerosol runtime object.

    Parameters mirror C ``rt_aerosol_runtime_prepare`` for the production
    log-linear path.  ``user_aod_ref`` may be zero for the AOD=0 control test;
    negative AOD is rejected.  AOD/SSA interpolation is linear and phase
    P11/P12/P33 wavelength interpolation is PCHIP, exactly as in the C update.
    """
    wavelength_nm = require_wavelength(
        wavelength_nm, context='aerosol runtime target')
    user_aod_ref = float(user_aod_ref)
    aod_ref_nm = require_wavelength(
        aod_ref_nm, context='aerosol runtime reference')
    L_max = int(L_max)
    if mie is None:
        raise ValueError("MieData is required to prepare an aerosol runtime")
    if user_aod_ref < 0.0:
        raise ValueError("AOD must be non-negative")
    if aod_ref_nm <= 0.0:
        raise ValueError("AOD reference wavelength must be positive")
    if L_max < 2:
        raise ValueError("L_max must be >= 2")

    aod_target = aod_at_wavelength(mie, wavelength_nm, user_aod_ref,
                                   aod_ref_nm)
    extinction_ratio = (aod_target / user_aod_ref) if user_aod_ref > 0.0 else \
        (aod_at_wavelength(mie, wavelength_nm, 1.0, aod_ref_nm))
    wl_um = wavelength_nm / 1000.0
    ssa_raw = float(np.interp(wl_um, mie.wavelengths, mie.spectral[:, 2]))
    P11, P12, P33 = phase_matrix_at_wavelength_pchip(mie, wavelength_nm)
    P11t, P12t, P33t, trunc_A = loglin_truncate(
        P11, P12, P33, mie.angles,
        mu1=float(loglin_mu1), mu2=float(loglin_mu2),
        thresh=float(loglin_threshold))

    tau_a_eff = aod_target * (1.0 - 0.5 * ssa_raw * trunc_A)
    den = 1.0 - 0.5 * ssa_raw * trunc_A
    ssa_a_eff = ((1.0 - 0.5 * trunc_A) * ssa_raw / den
                 if trunc_A > 0.0 else ssa_raw)

    betal, gammal, alphal, zetal = compute_vector_legendre(
        P11t, P12t, P33t, mie.angles, L_max, 4096)
    beta0 = compute_legendre_moments(P11t, mie.angles, 1, 4096)[0]
    if beta0 == 0.0:
        raise ValueError("zero beta0 after aerosol phase truncation")
    theta_asc, P11_asc = _ascending(mie.angles, P11t)
    _, P12_asc = _ascending(mie.angles, P12t)
    _, P33_asc = _ascending(mie.angles, P33t)

    source_key = getattr(mie, 'cache_key', ('memory', id(mie)))
    return AerosolRuntime(
        source_key=tuple(source_key),
        wavelength_nm=wavelength_nm,
        aod_ref=user_aod_ref,
        aod_ref_nm=aod_ref_nm,
        aod_target=float(aod_target),
        extinction_ratio=float(extinction_ratio),
        ssa_raw=float(ssa_raw),
        tau_a_eff=float(tau_a_eff),
        ssa_a_eff=float(ssa_a_eff),
        truncation_A=float(trunc_A),
        L_max=L_max,
        theta_asc=_readonly(theta_asc),
        P11_norm_asc=_readonly(P11_asc / beta0),
        P12_norm_asc=_readonly(P12_asc / beta0),
        P33_norm_asc=_readonly(P33_asc / beta0),
        betal=_readonly(betal),
        gammal=_readonly(gammal),
        alphal=_readonly(alphal),
        zetal=_readonly(zetal),
    )


# ---- gauss-quadrature vector Legendre (rt_aerosol_compute_vector_legendre_gauss) ----
def _clamped_cubic_spline_build(x, y):
    """clamped_cubic_spline_build: end slopes from first/last finite differences
    (Numerical-Recipes tridiagonal solve). Returns y2 (second derivatives)."""
    n = len(x)
    y2 = np.zeros(n)
    u = np.zeros(n - 1)
    dx0 = x[1] - x[0]
    dxn = x[n - 1] - x[n - 2]
    yp1 = (y[1] - y[0]) / dx0
    ypn = (y[n - 1] - y[n - 2]) / dxn
    y2[0] = -0.5
    u[0] = (3.0 / dx0) * ((y[1] - y[0]) / dx0 - yp1)
    for i in range(1, n - 1):
        dxm = x[i] - x[i - 1]
        dxp = x[i + 1] - x[i]
        sig = dxm / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        dd = ((y[i + 1] - y[i]) / dxp - (y[i] - y[i - 1]) / dxm)
        u[i] = (6.0 * dd / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p
    qn = 0.5
    un = (3.0 / dxn) * (ypn - (y[n - 1] - y[n - 2]) / dxn)
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)
    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
    return y2


def _clamped_cubic_spline_eval(x, y, y2, xq):
    """Vectorized clamped_cubic_spline_eval with edge clamping."""
    xq = np.asarray(xq, float)
    n = len(x)
    out = np.empty_like(xq)
    below = xq <= x[0]
    above = xq >= x[-1]
    mid = ~(below | above)
    out[below] = y[0]
    out[above] = y[-1]
    if np.any(mid):
        xm = xq[mid]
        khi = np.searchsorted(x, xm, side='right')
        khi = np.clip(khi, 1, n - 1)
        klo = khi - 1
        h = x[khi] - x[klo]
        a = (x[khi] - xm) / h
        b = (xm - x[klo]) / h
        out[mid] = (a * y[klo] + b * y[khi] +
                    ((a**3 - a) * y2[klo] + (b**3 - b) * y2[khi]) * (h * h) / 6.0)
    return out


def _build_mu_table_unique(P, theta_deg):
    """mu = cos(theta) sorted ascending, duplicate mu collapsed (keep latest)."""
    mu = np.cos(np.asarray(theta_deg, float) * np.pi / 180.0)
    mu = np.clip(mu, -1.0, 1.0)
    order = np.argsort(mu, kind='stable')
    xs = mu[order]
    ys = np.asarray(P, float)[order]
    keep_x = [xs[0]]
    keep_y = [ys[0]]
    for i in range(1, len(xs)):
        if abs(xs[i] - keep_x[-1]) < 1e-14:
            keep_y[-1] = ys[i]
        else:
            keep_x.append(xs[i])
            keep_y.append(ys[i])
    return np.array(keep_x), np.array(keep_y)


def _legendre_P_all(l_max, xq):
    """P_l(xq) for l=0..l_max (Bonnet recurrence)."""
    P = np.zeros(l_max + 1)
    P[0] = 1.0
    if l_max >= 1:
        P[1] = xq
    for l in range(1, l_max):
        P[l + 1] = ((2.0 * l + 1.0) * xq * P[l] - l * P[l - 1]) / (l + 1.0)
    return P


def _spin2_m0(l_max, xq):
    """rt_phase_spin2_m0: spin-2 associated basis P^l_{0,2} at m=0."""
    out = np.zeros(l_max + 1)
    if l_max < 2:
        return out
    out[2] = 3.0 * (1.0 - xq * xq) / (2.0 * np.sqrt(6.0))
    for ell in range(2, l_max):
        scale = (2.0 * ell + 1.0) / np.sqrt((ell - 1.0) * (ell + 3.0))
        back = np.sqrt((ell - 2.0) * (ell + 2.0)) / (2.0 * ell + 1.0)
        out[ell + 1] = scale * (xq * out[ell] - back * out[ell - 1])
    return out


def _complete_spherical(l_max, beta11, beta22, delta33):
    """rt_phase_complete_spherical -> (alpha, zeta)."""
    alpha = np.zeros(l_max + 1)
    zeta = np.zeros(l_max + 1)
    for ell in range(2, l_max + 1):
        e = float(ell)
        lower = e * (e - 1.0) / ((e + 1.0) * (e + 2.0))
        coup = 4.0 * (2.0 * e + 1.0) / (e * (e - 1.0) * (e + 1.0) * (e + 2.0))
        base = (e - 1.0) * (e - 1.0)
        sb = sd = ob = od = 0.0
        for order in range(ell - 1, -1, -1):
            gap = ell - order
            w = base - 1.5 * (gap - 1) * (ell + order)
            if (gap & 1) == 0:
                sb += w * beta11[order]
                sd += w * delta33[order]
            else:
                ob += w * beta22[order]
                od += w * delta33[order]
        zeta[ell] = lower * delta33[ell] - coup * (sd - ob)
        alpha[ell] = lower * beta11[ell] - coup * (sb - od)
    return alpha, zeta


def compute_vector_legendre_gauss(P11, P12, P33, theta_deg, L_max, mie_n_mu=60):
    """rt_aerosol_compute_vector_legendre_gauss (moment_mode=1, C default).
    Clamped cubic spline in mu + Gauss-Legendre quadrature over +/-mu.
    Returns betal,gammal,alphal,zetal normalized by beta0."""
    from .kernel import gauss_legendre_pos
    x11, y11 = _build_mu_table_unique(P11, theta_deg)
    x12, y12 = _build_mu_table_unique(P12, theta_deg)
    x33, y33 = _build_mu_table_unique(P33, theta_deg)
    y2_11 = _clamped_cubic_spline_build(x11, y11)
    y2_12 = _clamped_cubic_spline_build(x12, y12)
    y2_33 = _clamped_cubic_spline_build(x33, y33)
    mu, w = gauss_legendre_pos(mie_n_mu)     # positive nodes of GL(2*mie_n_mu)

    betal = np.zeros(L_max + 1)
    beta22 = np.zeros(L_max + 1)
    deltal = np.zeros(L_max + 1)
    gammal = np.zeros(L_max + 1)
    for ih in range(mie_n_mu):
        for sgn in (-1, 1):
            xq = -mu[ih] if sgn < 0 else mu[ih]
            wt = w[ih]
            p11 = float(_clamped_cubic_spline_eval(x11, y11, y2_11, xq))
            p12 = float(_clamped_cubic_spline_eval(x12, y12, y2_12, xq))
            p33 = float(_clamped_cubic_spline_eval(x33, y33, y2_33, xq))
            PL = _legendre_P_all(L_max, xq)
            P2 = _spin2_m0(L_max, xq)
            c = (2.0 * np.arange(L_max + 1) + 1.0) * 0.5 * wt
            betal += c * p11 * PL
            beta22 += c * p11 * PL          # spherical external: P22 = P11
            deltal += c * p33 * PL
            gammal += c * p12 * P2
    alphal, zetal = _complete_spherical(L_max, betal, beta22, deltal)
    b0 = betal[0]
    if b0 == 0.0 or not np.isfinite(b0):
        raise ValueError('beta0 non-normalizable in gauss vector Legendre')
    return betal / b0, gammal / b0, alphal / b0, zetal / b0


def water_component_phase_moments_gauss(mie, wl_nm, L_max=200, mie_n_mu=60):
    """water_component_phase_moments (moment_mode=1): eval_aerosol_phase at the
    mie angle grid for wl_nm, then gauss vector Legendre.  Also returns the
    PCHIP-integral bb/b ratio (build_aerosol_interpolators)."""
    from .constituent import bb_b_ratio as _bbr
    wl_um = wl_nm / 1000.0
    order = np.argsort(mie.angles, kind='stable')
    ang = mie.angles[order]
    def _interp_block(B):
        Br = B[order, :]
        return np.array([np.interp(wl_um, mie.phase_wavelengths, Br[k, :])
                         for k in range(len(ang))])
    p11 = _interp_block(mie.P11)
    p12 = _interp_block(mie.P12)
    p33 = _interp_block(mie.P33)
    be, ga, al, ze = compute_vector_legendre_gauss(p11, p12, p33, ang, L_max, mie_n_mu)
    return be, ga, al, ze, _bbr(mie, wl_nm)
