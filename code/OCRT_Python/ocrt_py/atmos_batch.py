"""Batched (case-axis b) Legendre + phase-kernel helpers for the atmosphere
SOS, GPU-capable via backend.xp.

The single-case functions in kernel.py operate on a shared ring rm (length
2*n_mu+1) and build plm/rrl/rtl of shape (L+1, dirs) and phase kernels
pfm/gr/gt/arr/art/att of shape (n_mu+1, dirs).  For our production grid the
ring differs per case (solar/view slots are case-specific), so every quantity
carries a leading batch axis b:

    rm    (B, dirs)
    plm   (B, L+1, dirs)   rrl/rtl (B, L+1, dirs)
    pfm   (B, n_mu+1, dirs)  gr/gt/arr/art/att (B, n_mu+1, dirs)

Greek coefficients (betal/gammal/alphal/zetal) are (B, L+1) since aerosol
model / IOP differ per case.  l-recurrences run sequentially over l but the
b and dirs axes are vectorized; there is no per-case Python loop.

Verified elsewhere against the single-case kernel.py path (bit-identical on
numpy).
"""
import numpy as np
from math import sqrt, pi, exp

from .backend import xp, to_np

M_PI = pi


class BatchLegWorkspace:
    """Batched analogue of kernel.LegendreWorkspace.  Arrays are backend.xp."""
    def __init__(self, B, n_mu, l_max):
        self.B = B
        self.n_mu = n_mu
        self.l_max = l_max
        dirs = 2 * n_mu + 1
        self.dirs = dirs
        self.plm = xp.zeros((B, l_max + 1, dirs))
        self.rrl = xp.zeros((B, l_max + 1, dirs))
        self.rtl = xp.zeros((B, l_max + 1, dirs))
        self.pfm = xp.zeros((B, n_mu + 1, dirs))
        self.gr = xp.zeros((B, n_mu + 1, dirs))
        self.gt = xp.zeros((B, n_mu + 1, dirs))
        self.arr = xp.zeros((B, n_mu + 1, dirs))
        self.art = xp.zeros((B, n_mu + 1, dirs))
        self.att = xp.zeros((B, n_mu + 1, dirs))


def legendre_compute_batch(ws, rm, m):
    """rt_legendre_compute over a batch.  rm (B, dirs) signed mu.  Fills
    ws.plm; returns xpl = plm[:, 2, :] (B, dirs)."""
    n_mu, l_max = ws.n_mu, ws.l_max
    if l_max < m or l_max < 2:
        raise ValueError("l_max too small")
    plm = ws.plm
    c = rm                                      # (B, dirs)
    plm[:] = 0.0
    if m == 0:
        plm[:, 0, :] = 1.0
        plm[:, 1, :] = c
        plm[:, 2, :] = 0.5 * (3.0 * c * c - 1.0)
    elif m == 1:
        sqrt3 = sqrt(3.0)
        x = 1.0 - c * c
        plm[:, 1, :] = xp.sqrt(0.5 * x)
        plm[:, 2, :] = c * plm[:, 1, :] * sqrt3
    else:
        a = 1.0
        for i in range(1, m + 1):
            a *= sqrt((i + m) / i) * 0.5
        xx = 1.0 - c * c
        plm[:, m - 1, :] = 0.0
        plm[:, m, :] = a * xp.power(xx, 0.5 * m)
    l_start = 2 if m < 2 else m
    for l in range(l_start, l_max):
        a_rec = (2.0 * l + 1.0) / sqrt((l + m + 1) * (l - m + 1))
        b_rec = sqrt((l + m) * (l - m)) / (2.0 * l + 1.0)
        plm[:, l + 1, :] = a_rec * (c * plm[:, l, :] - b_rec * plm[:, l - 1, :])
    return plm[:, 2, :].copy()


def legendre_compute_pol_batch(ws, rm, m):
    """rt_legendre_compute_pol over a batch.  Returns (xpl, xrl, xtl), each
    (B, dirs).  Fills ws.plm/rrl/rtl."""
    xpl = legendre_compute_batch(ws, rm, m)
    n_mu, l_max = ws.n_mu, ws.l_max
    rrl, rtl = ws.rrl, ws.rtl
    c = rm
    lim = min(max(m, 2), l_max + 1)
    rrl[:, :lim, :] = 0.0
    rtl[:, :lim, :] = 0.0
    if m == 0:
        xx = 1.0 - c * c
        rrl[:, 2, :] = 3.0 * xx / (2.0 * sqrt(6.0))
        rtl[:, 2, :] = 0.0
    elif m == 1:
        x = 1.0 - c * c
        sx = xp.sqrt(x)
        rrl[:, 2, :] = -c * sx * 0.5
        rtl[:, 2, :] = -sx * 0.5
    else:
        a = 1.0
        for i in range(1, m + 1):
            a *= sqrt((i + m) / i) * 0.5
        b = a * sqrt(m / (m + 1)) * sqrt((m - 1) / (m + 2))
        xx = 1.0 - c * c
        pxx = xp.power(xx, 0.5 * m - 1.0)
        rrl[:, m, :] = b * (1.0 + c * c) * pxx
        rtl[:, m, :] = b * 2.0 * c * pxx
    l_start = 2 if m < 2 else m
    for l in range(l_start, l_max):
        d = (l + 1) * (2 * l + 1) / sqrt((l + 3) * (l - 1) * (l + m + 1) * (l - m + 1))
        e = sqrt((l + 2) * (l - 2) * (l + m) * (l - m)) / (l * (2 * l + 1))
        f = 2.0 * m / (l * (l + 1))
        rrl[:, l + 1, :] = d * (c * rrl[:, l, :] - f * rtl[:, l, :] - e * rrl[:, l - 1, :])
        rtl[:, l + 1, :] = d * (c * rtl[:, l, :] - f * rrl[:, l, :] - e * rtl[:, l - 1, :])
    xrl = rrl[:, 2, :].copy()
    xtl = rtl[:, 2, :].copy()
    return xpl, xrl, xtl


def kernel_phase_fourier_batch(ws, m, betal):
    """pfm[b,j,k] = sum_l plm[b,l,j] plm[b,l,k] betal[b,l].  betal (B, L+1)."""
    n_mu, l_max = ws.n_mu, ws.l_max
    if betal is None:
        ws.pfm[:] = 0.0
        return
    b = betal[:, : l_max + 1]                       # (B, L+1)
    P = ws.plm[:, m: l_max + 1, :]                   # (B, L, dirs)
    Pw = P * b[:, m:, None]                          # weight
    Pj = P[:, :, n_mu:]                              # (B, L, n_mu+1)
    # pfm[b] = Pj[b].T @ Pw[b]  -> (B, n_mu+1, dirs)
    ws.pfm[:] = xp.einsum('blj,bld->bjd', Pj, Pw)


def kernel_phase_fourier_pol_batch(ws, m, gammal):
    n_mu, l_max = ws.n_mu, ws.l_max
    if gammal is None:
        ws.gr[:] = 0.0
        ws.gt[:] = 0.0
        return
    g = gammal[:, : l_max + 1]
    P = ws.plm[:, m: l_max + 1, :]
    R = ws.rrl[:, m: l_max + 1, :]
    T = ws.rtl[:, m: l_max + 1, :]
    Pj = P[:, :, n_mu:]                              # (B, L, n_mu+1)
    ws.gr[:] = xp.einsum('blj,bld->bjd', Pj, R * g[:, m:, None])
    ws.gt[:] = xp.einsum('blj,bld->bjd', Pj, T * g[:, m:, None])


def kernel_phase_fourier_aerosol_full_batch(ws, m, alphal, zetal):
    n_mu, l_max = ws.n_mu, ws.l_max
    if alphal is None or zetal is None:
        ws.arr[:] = 0.0
        ws.art[:] = 0.0
        ws.att[:] = 0.0
        return
    a = alphal[:, : l_max + 1][:, m:]
    z = zetal[:, : l_max + 1][:, m:]
    R = ws.rrl[:, m: l_max + 1, :]
    T = ws.rtl[:, m: l_max + 1, :]
    Rj, Tj = R[:, :, n_mu:], T[:, :, n_mu:]
    ws.att[:] = (xp.einsum('blj,bld->bjd', Tj, T * a[:, :, None]) +
                 xp.einsum('blj,bld->bjd', Rj, R * z[:, :, None]))
    ws.arr[:] = (xp.einsum('blj,bld->bjd', Tj, T * z[:, :, None]) +
                 xp.einsum('blj,bld->bjd', Rj, R * a[:, :, None]))
    ws.art[:] = (xp.einsum('blj,bld->bjd', Tj, R * a[:, :, None]) +
                 xp.einsum('blj,bld->bjd', Rj, T * z[:, :, None]))


# ---- batched medium + primary source + integrate ---------------------------
class BatchAtm:
    """Batched in-atmosphere medium.  Scalars beta0/beta2/gamma2/alpha2 are
    Rayleigh (case-independent).  Per-case arrays carry axis b."""
    __slots__ = ('B', 'nt', 'n_mu', 'beta0', 'beta2', 'gamma2', 'alpha2',
                 'h', 'ch', 'xdel', 'ydel', 'rm', 'gb', 'mu_pos', 'mu_sun',
                 'betal', 'gammal', 'alphal', 'zetal', 'L_max', 'beam_q',
                 'bq_active')

    def to_xp(self):
        """Move every per-case array onto the active backend (GPU when
        OCRT_PY_GPU=1).  Needed for media stacked with np.stack so they don't
        mix numpy arrays with CuPy workspace arrays."""
        for name in ('h', 'ch', 'xdel', 'ydel', 'rm', 'gb', 'mu_pos',
                     'betal', 'gammal', 'alphal', 'zetal', 'beam_q'):
            v = getattr(self, name, None)
            if v is not None and not np.isscalar(v):
                setattr(self, name, xp.asarray(v))
        return self


def build_inwater_atm_batch(nt_case, tau_total, omega_w, omega_particle,
                            mu_sun_water, rings, weights, betal, gammal,
                            alphal, zetal, ntw, Lmix, beam_q):
    """Batched in-water medium build over the case axis.  Replaces the per-case
    build_inwater_atm + _pad_water_medium + xp.stack loop with a single set of
    vectorized ops (no per-case GPU launch).  All per-case inputs are (B,)
    except rings/weights (B, n_mu_w), Greek (B, Lmix+1) and beam_q (B,).  Bottom
    identity padding is baked in via the k <= nt_case mask (k > nt_case -> h held
    constant so dtau=0, ch/xdel/ydel=0 so no source).  Bit-identical to the
    per-case path.  Returns a ready BatchAtm."""
    from math import sqrt
    nt_case = xp.asarray(nt_case, dtype=float)
    tau_total = xp.asarray(tau_total, dtype=float)
    omega_w = xp.asarray(omega_w, dtype=float)
    omega_particle = xp.asarray(omega_particle, dtype=float)
    mu_sun_water = xp.asarray(mu_sun_water, dtype=float)
    rings = xp.asarray(rings, dtype=float)
    weights = xp.asarray(weights, dtype=float)
    B, n_mu_w = rings.shape
    delta_w = 0.039
    ron_w = 2.0 * (1.0 - delta_w) / (2.0 + delta_w)
    k = xp.arange(ntw + 1, dtype=float)[None, :]           # (1, ntw+1)
    ntc = nt_case[:, None]                                 # (B, 1)
    valid = k <= ntc
    h = xp.minimum(k, ntc) * (tau_total / nt_case)[:, None]
    ch = xp.where(valid, 0.5 * xp.exp(-h / mu_sun_water[:, None]), 0.0)
    xdel = xp.where(valid, omega_particle[:, None], 0.0)
    ydel = xp.where(valid, omega_w[:, None], 0.0)
    zeros_col = xp.zeros((B, 1))
    rm = xp.concatenate([-xp.flip(rings, axis=1), -mu_sun_water[:, None], rings], axis=1)
    gb = xp.concatenate([xp.flip(weights, axis=1), zeros_col, weights], axis=1)
    bAw = BatchAtm()
    bAw.B = B; bAw.nt = ntw; bAw.n_mu = n_mu_w; bAw.L_max = Lmix
    bAw.beta0 = 1.0; bAw.beta2 = 0.5 * ron_w
    bAw.gamma2 = -ron_w * sqrt(1.5); bAw.alpha2 = 3.0 * ron_w
    bAw.h = h; bAw.ch = ch; bAw.xdel = xdel; bAw.ydel = ydel
    bAw.rm = rm; bAw.gb = gb; bAw.mu_pos = rings; bAw.mu_sun = mu_sun_water
    bAw.betal = xp.asarray(betal); bAw.gammal = xp.asarray(gammal)
    bAw.alphal = xp.asarray(alphal); bAw.zetal = xp.asarray(zetal)
    bAw.beam_q = xp.asarray(beam_q, dtype=float)
    bAw.bq_active = bool(xp.any(bAw.beam_q != 0.0))
    return bAw


def build_atm_aerosol_batch(nt, tau_R, tau_a, ssa_a, mu_sun, rings, weights,
                            betal, gammal, alphal, zetal, L_max, aer_h_km=2.0):
    """Batched build_atm_aerosol over the case axis.  tau_R is shared within a
    band; tau_a/ssa_a/mu_sun are (B,); rings/weights are (B, n_mu_atm) (the GL
    ring with the view node inserted); Greek is (B, L_max+1).  Returns
    (bAa, z_km_level) with z_km_level (B, nt+1).  The 64-step altitude bisection
    and the US62 / exponential aerosol column profiles run vectorized over the
    (B, nt) layer grid, matching the per-case build_atm_aerosol bit for bit.
    beam_q=0 (unpolarized solar beam).  The view-node ring index (vj) is left to
    the caller since it is a host scalar."""
    from math import sqrt
    from .atmos import DEPOL, _tau_a_above_exp
    from . import molprofile as mp
    tau_a = xp.asarray(tau_a, dtype=float)
    ssa_a = xp.asarray(ssa_a, dtype=float)
    mu_sun = xp.asarray(mu_sun, dtype=float)
    rings = xp.asarray(rings, dtype=float)
    weights = xp.asarray(weights, dtype=float)
    B, n_mu = rings.shape
    tau_total = tau_R + tau_a                              # (B,)
    ron = 2.0 * (1.0 - DEPOL) / (2.0 + DEPOL)
    z_top = 100.0
    k_arr = xp.arange(nt + 1, dtype=float)[None, :]        # (1, nt+1)
    h = k_arr * tau_total[:, None] / nt                    # (B, nt+1); (k*tau)/nt
    ch = 0.5 * xp.exp(-h / mu_sun[:, None])
    xdel = xp.zeros((B, nt + 1)); ydel = xp.zeros((B, nt + 1))
    if tau_R == 0.0:
        xdel[:, 0] = ssa_a; ydel[:, 0] = 0.0
    else:
        xdel[:, 0] = 0.0; ydel[:, 0] = 1.0
    # altitude bisection: find z_j with tau_ray_above(z)+tau_a_above(z)=h[j]
    target = h[:, 1:nt + 1]                                # (B, nt)
    z_lo = xp.zeros((B, nt)); z_hi = xp.full((B, nt), z_top)
    for _ in range(64):
        z_mid = 0.5 * (z_lo + z_hi)
        tau_above = (tau_R * mp.us62_grid_fraction_above(z_mid)
                     + _tau_a_above_exp(z_mid, tau_a[:, None], aer_h_km, z_top))
        above = tau_above > target
        z_lo = xp.where(above, z_mid, z_lo)
        z_hi = xp.where(above, z_hi, z_mid)
    z_j = 0.5 * (z_lo + z_hi)
    z_j[:, nt - 1] = 0.0
    z_km_level = xp.zeros((B, nt + 1))
    z_km_level[:, 0] = z_top
    z_km_level[:, 1:] = z_j
    ray_full = tau_R * mp.us62_grid_fraction_above(z_km_level)      # (B, nt+1)
    ea_full = _tau_a_above_exp(z_km_level, tau_a[:, None], aer_h_km, z_top)
    dt_ray = ray_full[:, 1:] - ray_full[:, :-1]                     # (B, nt)
    dt_aer = ea_full[:, 1:] - ea_full[:, :-1]
    dt = dt_ray + dt_aer
    pos = dt > 0.0
    ydel_else = 1.0 if tau_R > 0.0 else 0.0
    xdel[:, 1:] = xp.where(pos, dt_aer * ssa_a[:, None] / dt, 0.0)
    ydel[:, 1:] = xp.where(pos, dt_ray / dt, ydel_else)
    rm = xp.concatenate([-xp.flip(rings, axis=1), -mu_sun[:, None], rings], axis=1)
    gb = xp.concatenate([xp.flip(weights, axis=1), xp.zeros((B, 1)), weights], axis=1)
    bAa = BatchAtm()
    bAa.B = B; bAa.nt = nt; bAa.n_mu = n_mu; bAa.L_max = L_max; bAa.beam_q = 0.0
    bAa.beta0 = 1.0; bAa.beta2 = 0.5 * ron
    bAa.gamma2 = -ron * sqrt(1.5); bAa.alpha2 = 3.0 * ron
    bAa.h = h; bAa.ch = ch; bAa.xdel = xdel; bAa.ydel = ydel
    bAa.rm = rm; bAa.gb = gb; bAa.mu_pos = rings; bAa.mu_sun = mu_sun
    bAa.betal = xp.asarray(betal); bAa.gammal = xp.asarray(gammal)
    bAa.alphal = xp.asarray(alphal); bAa.zetal = xp.asarray(zetal)
    bAa.bq_active = False
    return bAa, z_km_level


def build_atm_rayleigh_batch(nt, tau_R, mu_sun, rings, weights):
    """Batched build_atm_rayleigh (no aerosol) over the case axis.  tau_R shared
    within a band; mu_sun (B,); rings/weights (B, n_mu) [GL + view node].
    Returns (bA, z_km_level (B, nt+1)).  z_km_level uses the US62 direct inverse
    (us62_altitude_from_grid_fraction), matching the per-case build_atm_rayleigh
    bit for bit (no aerosol bisection)."""
    from math import sqrt
    from .atmos import DEPOL
    from . import molprofile as mp
    mu_sun = xp.asarray(mu_sun, dtype=float)
    rings = xp.asarray(rings, dtype=float)
    weights = xp.asarray(weights, dtype=float)
    B, n_mu = rings.shape
    ron = 2.0 * (1.0 - DEPOL) / (2.0 + DEPOL)
    k_arr = xp.arange(nt + 1, dtype=float)[None, :]
    h = xp.broadcast_to(k_arr * tau_R / nt, (B, nt + 1)) * xp.ones((B, 1))  # (B,nt+1)
    ch = 0.5 * xp.exp(-h / mu_sun[:, None])
    xdel = xp.zeros((B, nt + 1)); ydel = xp.ones((B, nt + 1))
    z_top = 100.0
    z_km = xp.zeros((B, nt + 1))
    z_km[:, 0] = z_top
    z_km[:, 1:] = mp.us62_altitude_from_grid_fraction(h[:, 1:] / tau_R)
    rm = xp.concatenate([-xp.flip(rings, axis=1), -mu_sun[:, None], rings], axis=1)
    gb = xp.concatenate([xp.flip(weights, axis=1), xp.zeros((B, 1)), weights], axis=1)
    bA = BatchAtm()
    bA.B = B; bA.nt = nt; bA.n_mu = n_mu; bA.L_max = 2; bA.beam_q = 0.0
    bA.beta0 = 1.0; bA.beta2 = 0.5 * ron
    bA.gamma2 = -ron * sqrt(1.5); bA.alpha2 = 3.0 * ron
    bA.h = h; bA.ch = ch; bA.xdel = xdel; bA.ydel = ydel
    bA.rm = rm; bA.gb = gb; bA.mu_pos = rings; bA.mu_sun = mu_sun
    bA.bq_active = False
    return bA, z_km


def apply_gas_absorption_batch(bAa, absorption, wl_nm, z_km_level):
    """Batched rt_atm_apply_gas_absorption.  Modifies bAa.h/ch/xdel/ydel in
    place over the case axis and returns tau_abs_total (B,).  z_km_level is
    (B, nt+1).  Per-layer absorption is computed by flattening the (B, nt) layer
    grid to a 1-D vector for the array-safe column/xsec helpers, then reshaping.
    Bit-identical to the per-case apply_gas_absorption."""
    from .absorption import (afgl_cumulative_column_vec,
                             xsec_interp_layer_z_vec, N_GAS)
    nt = bAa.nt
    at = absorption.atm
    z_km_level = xp.asarray(z_km_level, dtype=float)
    B = z_km_level.shape[0]
    z_hi = z_km_level[:, :nt]                              # (B, nt)
    z_lo = z_km_level[:, 1:nt + 1]
    zlo_s = xp.minimum(z_lo, z_hi); zhi_s = xp.maximum(z_lo, z_hi)
    z_mid = 0.5 * (zlo_s + zhi_s)
    zlo_f = zlo_s.reshape(-1); zhi_f = zhi_s.reshape(-1); zmid_f = z_mid.reshape(-1)
    tau_abs_layer = xp.zeros((B, nt))
    for g in range(N_GAS):
        xs = absorption.xsec[g]
        if not xs.available:
            continue
        Ndef = at.column_default[g]
        N_layer = (afgl_cumulative_column_vec(at, g, zlo_f)
                   - afgl_cumulative_column_vec(at, g, zhi_f)).reshape(B, nt)
        if absorption.overridden[g] and Ndef > 0.0:
            N = absorption.column_eff[g] * (N_layer / Ndef)
        else:
            N = N_layer
        sigma = xsec_interp_layer_z_vec(xs, at.z_km, zmid_f, wl_nm).reshape(B, nt)
        tau_abs_layer = tau_abs_layer + sigma * N
    h_old = bAa.h.copy() if hasattr(bAa.h, 'copy') else xp.array(bAa.h)
    cum_abs = xp.concatenate([xp.zeros((B, 1)), xp.cumsum(tau_abs_layer, axis=1)], axis=1)
    bAa.h = h_old + cum_abs
    bAa.ch = 0.5 * xp.exp(-bAa.h / bAa.mu_sun[:, None])
    dt_old = h_old[:, 1:] - h_old[:, :-1]
    dt_new = dt_old + tau_abs_layer
    ratio = xp.where((dt_new > 0.0) & (dt_old > 0.0), dt_old / dt_new, 1.0)
    bAa.xdel[:, 1:] = bAa.xdel[:, 1:] * ratio
    bAa.ydel[:, 1:] = bAa.ydel[:, 1:] * ratio
    return tau_abs_layer.sum(axis=1)


def primary_source_batch(atm, m, xpl, pfm, ws=None, xrl=None):
    """Batched rt_solver_primary_source.  xpl (B, dirs); pfm (B, n_mu+1, dirs).
    beam_q (B,) may be nonzero (in-water Fresnel-polarized solar beam): then
    the aerosol bq term needs ws.plm/rrl and atm.gammal, and the Rayleigh bq
    term needs xrl.  Returns src (B, nt+1, dirs)."""
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    dirs = 2 * n_mu + 1
    beta0_m = atm.beta0 if m == 0 else 0.0
    xpl0 = xpl[:, n_mu][:, None]                        # (B,1)
    sa_ray = beta0_m + atm.beta2 * xpl * xpl0           # (B, dirs)
    sa_aer = pfm[:, 0, :]                               # (B, dirs) solar row
    ch = atm.ch[:, :, None]
    xd = atm.xdel[:, :, None]
    yd = atm.ydel[:, :, None]
    base = ch * (xd * sa_aer[:, None, :] + yd * sa_ray[:, None, :])

    beam_q = getattr(atm, 'beam_q', 0.0)
    bq_active = getattr(atm, 'bq_active', None)
    if bq_active is None:
        bq_active = (np.any(np.asarray(to_np(beam_q)) != 0.0)
                     if not np.isscalar(beam_q) else (beam_q != 0.0))
    if bq_active:
        bqv = beam_q if not np.isscalar(beam_q) else xp.full((B,), beam_q)
        # aerosol bq: sum_l plm[l][j] rrl[l][n_mu] gammal[l], l=m..L_max
        gl = atm.gammal
        Lg = gl.shape[1] - 1
        P = ws.plm[:, m:Lg + 1, :]                      # (B, Lp, dirs)
        Rn = ws.rrl[:, m:Lg + 1, n_mu]                  # (B, Lp)
        bq_aer = xp.einsum('bld,bl->bd', P, Rn * gl[:, m:Lg + 1])
        # Rayleigh bq: gamma2*xpl*xrl0 (m<=2)
        if m <= 2:
            xrl0 = xrl[:, n_mu][:, None]
            bq_ray = atm.gamma2 * xpl * xrl0
        else:
            bq_ray = xp.zeros((B, dirs))
        bqf = bqv[:, None, None]
        base = base + ch * bqf * (xd * bq_aer[:, None, :] + yd * bq_ray[:, None, :])
    base[:, :, n_mu] = xp.nan
    return base


def primary_source_pol_batch(atm, m, xrl, xtl, xpl, gr, gt, ws=None):
    """Batched rt_solver_primary_source_pol.  beam_q (B,) may be nonzero: the
    aerosol bq_q/bq_u terms need ws.rrl/rtl and atm.alphal/zetal.  Returns
    (src_Q, src_U), each (B, nt+1, dirs)."""
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    dirs = 2 * n_mu + 1
    xpl0 = xpl[:, n_mu][:, None]                        # (B,1)
    rayleigh_active = (m <= 2)
    if rayleigh_active:
        sb_ray = atm.gamma2 * xrl * xpl0               # (B, dirs)
        sc_ray = atm.gamma2 * xtl * xpl0
    else:
        sb_ray = xp.zeros((B, dirs))
        sc_ray = xp.zeros((B, dirs))
    sb_aer = gr[:, 0, :]
    sc_aer = gt[:, 0, :]
    ch = atm.ch[:, :, None]
    xd = atm.xdel[:, :, None]
    yd = atm.ydel[:, :, None]
    src_q = +ch * (xd * sb_aer[:, None, :] + yd * sb_ray[:, None, :])
    src_u = -ch * (xd * sc_aer[:, None, :] + yd * sc_ray[:, None, :])

    beam_q = getattr(atm, 'beam_q', 0.0)
    bq_active = getattr(atm, 'bq_active', None)
    if bq_active is None:
        bq_active = (np.any(np.asarray(to_np(beam_q)) != 0.0)
                     if not np.isscalar(beam_q) else (beam_q != 0.0))
    if bq_active:
        bqv = beam_q if not np.isscalar(beam_q) else xp.full((B,), beam_q)
        al = atm.alphal; ze = atm.zetal
        La = al.shape[1] - 1
        R = ws.rrl[:, m:La + 1, :]                      # (B, Lp, dirs) rj
        T = ws.rtl[:, m:La + 1, :]                      # tj
        r0 = ws.rrl[:, m:La + 1, n_mu]                  # (B, Lp)
        t0 = ws.rtl[:, m:La + 1, n_mu]
        a = al[:, m:La + 1]; z = ze[:, m:La + 1]
        bq_q = (xp.einsum('bl,bld->bd', t0 * z, T) +
                xp.einsum('bl,bld->bd', r0 * a, R))
        bq_u = (xp.einsum('bl,bld->bd', r0 * a, T) +
                xp.einsum('bl,bld->bd', t0 * z, R))
        if rayleigh_active:
            xrl0 = xrl[:, n_mu][:, None]
            bq_q_ray = atm.alpha2 * xrl * xrl0
            bq_u_ray = atm.alpha2 * xtl * xrl0
        else:
            bq_q_ray = xp.zeros((B, dirs))
            bq_u_ray = xp.zeros((B, dirs))
        bqf = bqv[:, None, None]
        src_q = src_q + ch * bqf * (xd * bq_q[:, None, :] + yd * bq_q_ray[:, None, :])
        src_u = src_u - ch * bqf * (xd * bq_u[:, None, :] + yd * bq_u_ray[:, None, :])
    src_q[:, :, n_mu] = xp.nan
    src_u[:, :, n_mu] = xp.nan
    return src_q, src_u


def _bcs_scan_fwd(c, contrib, x0):
    """Solve the forward linear recurrence x[k] = c[k-1]*x[k-1] + contrib[k-1]
    for k=1..nt with x[0]=x0.  c: (B, nt, n_mu) shared transmission; contrib and
    x0 may carry a trailing Stokes axis S (contrib (B, nt, n_mu[, S]), x0
    (B, n_mu[, S])).  Returns x (B, nt+1, n_mu[, S]).

    Division-free affine-map Hillis-Steele prefix scan (see below); a stays in
    [0,1] and b stays bounded, so stable for arbitrary optical depth.  log2(nt)
    full-array steps replace nt sequential launches.  When contrib carries the
    Stokes axis, I/Q/U ride the same scan (c is shared), cutting the scan-call
    count by 3x with bit-identical results (no cross-Stokes reassociation)."""
    B, nt, n_mu = c.shape
    has_s = (contrib.ndim == 4)
    a = c.copy(); b = contrib.copy()
    d = 1
    while d < nt:
        a_sh = xp.concatenate([xp.ones((B, d, n_mu)), a[:, :-d, :]], axis=1)
        if has_s:
            S = b.shape[-1]
            b_sh = xp.concatenate([xp.zeros((B, d, n_mu, S)), b[:, :-d]], axis=1)
            b = a[..., None] * b_sh + b
        else:
            b_sh = xp.concatenate([xp.zeros((B, d, n_mu)), b[:, :-d, :]], axis=1)
            b = a * b_sh + b
        a = a_sh * a
        d *= 2
    if has_s:
        x = xp.empty((B, nt + 1, n_mu, b.shape[-1]))
        x[:, 0] = x0
        x[:, 1:] = a[..., None] * x0[:, None] + b
    else:
        x = xp.empty((B, nt + 1, n_mu))
        x[:, 0, :] = x0
        x[:, 1:, :] = a * x0[:, None, :] + b
    return x


def integrate_bcs_batch(atm, src, surface_bc=None, top_down_bc=None):
    """Batched rt_solver_integrate_bcs (LINEAR).  src (B, nt+1, dirs) or, with a
    trailing Stokes axis, (B, nt+1, dirs, S); surface_bc/top_down_bc (B, n_mu),
    (B, n_mu, S) or None.  Returns rad shaped like src (solar slot NaN).  Per-
    layer sweeps are vectorized linear-recurrence scans; OCRT_PY_SEQ_BCS=1 forces
    the bit-identical sequential loop (validation)."""
    import os
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    dirs = 2 * n_mu + 1
    squeeze = (src.ndim == 3)
    if squeeze:
        src = src[..., None]                            # (B, nt+1, dirs, 1)
    S = src.shape[-1]
    h = atm.h; mu = atm.mu_pos
    dtau = h[:, 1:] - h[:, :-1]
    dtau_safe = xp.where(dtau > 1e-300, dtau, 1.0)
    rad = xp.full((B, nt + 1, dirs, S), xp.nan)
    c = xp.exp(-dtau[:, :, None] / mu[:, None, :])      # (B, nt, n_mu)
    seq = os.environ.get('OCRT_PY_SEQ_BCS', '0') == '1'
    hlo = h[:, :-1, None, None]; hhi = h[:, 1:, None, None]
    mub = mu[:, None, :, None]; dts = dtau_safe[:, :, None, None]
    cS = c[:, :, :, None]

    def _bc(x):
        if x is None:
            return xp.zeros((B, n_mu, S))
        x = xp.asarray(x, dtype=float)
        return x[..., None] if x.ndim == 2 else x

    # upward (j>0): sweep k=nt-1..0
    J = src[:, :, n_mu + 1:]                            # (B, nt+1, n_mu, S)
    a = (J[:, 1:] - J[:, :-1]) / dts
    b = J[:, :-1] - a * hlo
    contrib = ((1.0 - cS) * (b + a * mub) + a * (hlo - hhi * cS)) * 0.5
    I0 = _bc(surface_bc)
    if seq:
        I = I0.copy()
        rad[:, nt, n_mu + 1:] = I
        for k in range(nt - 1, -1, -1):
            I = cS[:, k] * I + contrib[:, k]
            rad[:, k, n_mu + 1:] = I
    else:
        xup = _bcs_scan_fwd(xp.flip(c, axis=1), xp.flip(contrib, axis=1), I0)
        rad[:, :, n_mu + 1:] = xp.flip(xup, axis=1)

    # downward (j<0): abs-indexed columns j=-1..-n_mu; storage n_mu - j_abs
    Jd = src[:, :, :n_mu][:, :, ::-1]                   # (B, nt+1, n_mu, S) abs
    a = (Jd[:, 1:] - Jd[:, :-1]) / dts
    b = Jd[:, 1:] - a * hhi
    contrib = ((1.0 - cS) * (b + a * (-mub)) + a * (hhi - hlo * cS)) * 0.5
    I0d = _bc(top_down_bc)
    if seq:
        I = I0d.copy()
        rad[:, 0, :n_mu] = I[:, ::-1]
        for k in range(1, nt + 1):
            I = cS[:, k - 1] * I + contrib[:, k - 1]
            rad[:, k, :n_mu] = I[:, ::-1]
    else:
        xdn = _bcs_scan_fwd(c, contrib, I0d)
        rad[:, :, :n_mu] = xdn[:, :, ::-1]
    if squeeze:
        rad = rad[..., 0]
    return rad


def integrate_bcs_iqu(atm, src_i, src_q, src_u,
                      bc_i=None, bc_q=None, bc_u=None, top=False):
    """Run integrate_bcs_batch for I, Q, U at once (Stokes stacked on the last
    axis) so the shared-transmission scan runs a third as often.  bc_* are
    surface BCs unless top=True, then top-down BCs.  Returns (rad_i,rad_q,rad_u).
    Bit-identical to three separate calls (Stokes are independent here)."""
    src = xp.stack([src_i, src_q, src_u], axis=-1)      # (B, nt+1, dirs, 3)
    if bc_i is None and bc_q is None and bc_u is None:
        bc = None
    else:
        n_mu = atm.n_mu
        z = xp.zeros((atm.B, n_mu))
        def _b(x): return z if x is None else xp.asarray(x, dtype=float)
        bc = xp.stack([_b(bc_i), _b(bc_q), _b(bc_u)], axis=-1)   # (B, n_mu, 3)
    rad = (integrate_bcs_batch(atm, src, None, bc) if top
           else integrate_bcs_batch(atm, src, bc, None))
    return rad[..., 0], rad[..., 1], rad[..., 2]


# ---- batched source build (sos_build_source_pol) ---------------------------
def sos_build_source_pol_batch(atm, m, kt, I_prev, Q_prev, U_prev,
                               xpl, xrl, xtl):
    """Batched sos_build_source_pol.  kt = (pfm, gr, gt, arr, art, att), each
    (B, n_mu+1, dirs).  I/Q/U_prev (B, nt+1, dirs).  xpl/xrl/xtl (B, dirs).
    Returns J_I, J_Q, J_U each (B, nt+1, dirs)."""
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    dirs = 2 * n_mu + 1
    beta0_m = atm.beta0 if m == 0 else 0.0
    beta2, gamma2, alpha2 = atm.beta2, atm.gamma2, atm.alpha2
    ray = (m <= 2)
    P, GR, GT, AR, AT, ATT = kt

    jp = np.arange(1, n_mu + 1)
    pk = n_mu + jp                          # +k / +jp column indices
    mk = n_mu - jp                          # -k / -jp column indices
    xpv = xpl[:, pk]                         # (B, n_mu)  xpl[+j]
    xrv = xrl[:, pk]
    xtv = xtl[:, pk]
    ypv = xpl[:, mk]                         # xpl[-j]
    yrv = xrl[:, mk]
    ytv = xtl[:, mk]
    w = atm.gb[:, pk]                        # (B, n_mu)

    def Tp(X, cols):
        # X (B, n_mu+1, dirs); rows jp, given cols -> (B, jp, k) then T -> (B,k,jp)
        return xp.transpose(X[:, jp][:, :, cols], (0, 2, 1))

    def Dr(X, cols):
        # rows kk(=jp range), given cols -> (B, k, jp) directly
        return X[:, jp][:, :, cols]

    a = [None] * 18
    a[0] = Tp(P, pk);  a[1] = Tp(P, mk)
    a[2] = Tp(GT, pk); a[3] = Tp(GT, mk)
    a[4] = Dr(GT, pk); a[5] = Dr(GT, mk)
    a[6] = Tp(GR, pk); a[7] = Tp(GR, mk)
    a[8] = Dr(GR, pk); a[9] = Dr(GR, mk)
    a[10] = Tp(AR, pk); a[11] = Tp(AR, mk)
    a[12] = Tp(AT, pk); a[13] = Tp(AT, mk)
    a[14] = Dr(AT, pk); a[15] = Dr(AT, mk)
    a[16] = Tp(ATT, pk); a[17] = Tp(ATT, mk)

    r = [None] * 18
    if ray:
        # outer products (B, k, jp): row index = first factor, col = second
        def out_kj(u_row, v_col):
            return u_row[:, :, None] * v_col[:, None, :]
        r[0] = beta0_m + beta2 * out_kj(xpv, xpv)
        r[1] = beta0_m + beta2 * out_kj(ypv, xpv)
        r[2] = gamma2 * out_kj(xtv, xpv)
        r[3] = gamma2 * out_kj(ytv, xpv)
        r[4] = gamma2 * out_kj(xpv, xtv)
        r[5] = gamma2 * out_kj(xpv, ytv)
        r[6] = gamma2 * out_kj(xrv, xpv)
        r[7] = gamma2 * out_kj(yrv, xpv)
        r[8] = gamma2 * out_kj(xpv, xrv)
        r[9] = gamma2 * out_kj(xpv, yrv)
        r[10] = alpha2 * out_kj(xrv, xrv)
        r[11] = alpha2 * out_kj(yrv, xrv)
        r[12] = alpha2 * out_kj(xrv, xtv)
        r[13] = alpha2 * out_kj(yrv, xtv)
        r[14] = alpha2 * out_kj(xtv, xrv)
        r[15] = alpha2 * out_kj(xtv, yrv)
        r[16] = alpha2 * out_kj(xtv, xtv)
        r[17] = alpha2 * out_kj(ytv, xtv)

    x = atm.xdel                            # (B, nt+1)
    y = atm.ydel

    def fld(F, cols):
        return xp.transpose(F[:, :, cols], (0, 2, 1))   # (B, jp, nk)
    Ip = fld(I_prev, pk); Im = fld(I_prev, mk)
    Qp = fld(Q_prev, pk); Qm = fld(Q_prev, mk)
    Up = fld(U_prev, pk); Um = fld(U_prev, mk)
    wIp = w[:, :, None] * Ip; wIm = w[:, :, None] * Im
    wQp = w[:, :, None] * Qp; wQm = w[:, :, None] * Qm
    wUp = w[:, :, None] * Up; wUm = w[:, :, None] * Um

    def mixdot(e, F):
        out = (a[e] @ F) * x[:, None, :]
        if ray:
            out = out + (r[e] @ F) * y[:, None, :]
        return out

    aI2 = mixdot(0, wIp) + mixdot(8, wQp) - mixdot(4, wUp) \
        + mixdot(1, wIm) + mixdot(9, wQm) - mixdot(5, wUm)
    aI1 = mixdot(1, wIp) + mixdot(9, wQp) + mixdot(5, wUp) \
        + mixdot(0, wIm) + mixdot(8, wQm) + mixdot(4, wUm)
    aQ2 = mixdot(6, wIp) + mixdot(10, wQp) - mixdot(12, wUp) \
        + mixdot(7, wIm) + mixdot(11, wQm) + mixdot(13, wUm)
    aQ1 = mixdot(7, wIp) + mixdot(11, wQp) - mixdot(13, wUp) \
        + mixdot(6, wIm) + mixdot(10, wQm) + mixdot(12, wUm)
    aU2 = -(mixdot(2, wIp) + mixdot(14, wQp) - mixdot(16, wUp)) \
        - (-mixdot(3, wIm) + mixdot(15, wQm) - mixdot(17, wUm))
    aU1 = -(mixdot(3, wIp) - mixdot(15, wQp) - mixdot(17, wUp)) \
        - (-mixdot(2, wIm) - mixdot(14, wQm) - mixdot(16, wUm))

    J_I = xp.full((B, nt + 1, dirs), xp.nan)
    J_Q = xp.full((B, nt + 1, dirs), xp.nan)
    J_U = xp.full((B, nt + 1, dirs), xp.nan)
    # aX* are (B, k, nk); place with k over the layer axis, direction over ±jp
    J_I[:, :, pk] = xp.transpose(aI2, (0, 2, 1))
    J_I[:, :, mk] = xp.transpose(aI1, (0, 2, 1))
    J_Q[:, :, pk] = xp.transpose(aQ2, (0, 2, 1))
    J_Q[:, :, mk] = xp.transpose(aQ1, (0, 2, 1))
    J_U[:, :, pk] = xp.transpose(aU2, (0, 2, 1))
    J_U[:, :, mk] = xp.transpose(aU1, (0, 2, 1))
    return J_I, J_Q, J_U


# ---- batched order iteration (black + Cox-Munk surface) --------------------
def _field_max_abs_batch(F, n_mu):
    """nanmax over directions excluding the solar slot, per case -> (B,)."""
    A = xp.abs(F)
    A2 = xp.concatenate([A[:, :, :n_mu], A[:, :, n_mu + 1:]], axis=2)
    return xp.nanmax(A2, axis=(1, 2))


def sos_atm_batch(atm, m, kt, prim_I, prim_Q, prim_U, xpl, xrl, xtl,
                  max_iterations, tolerance, R_m=None, mu_factor=None,
                  surf_seed=None, ext_init=None):
    """Batched atm vector SOS.  R_m None -> black surface; else Cox-Munk
    diffuse-reflection m-mode kernel R_m (B, n_mu, n_mu, 9).  surf_seed /
    ext_init (B, nt+1, dirs) or None.  Returns tot_I/Q/U (B, nt+1, dirs),
    n_orders (B,), conv (B,), resid (B,)."""
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    use_surface = R_m is not None
    if surf_seed is not None:
        sI, sQ, sU = surf_seed
        tot_I = prim_I + sI; tot_Q = prim_Q + sQ; tot_U = prim_U + sU
    else:
        tot_I = prim_I.copy(); tot_Q = prim_Q.copy(); tot_U = prim_U.copy()
    if ext_init is not None:
        eI, eQ, eU = ext_init
        tot_I = tot_I + eI; tot_Q = tot_Q + eQ; tot_U = tot_U + eU

    n_ord = xp.ones((B,), dtype=int)
    conv = xp.zeros((B,), dtype=bool)
    resid = xp.zeros((B,))
    if max_iterations == 1:
        if surf_seed is not None:
            tot_I = tot_I - sI; tot_Q = tot_Q - sQ; tot_U = tot_U - sU
        return tot_I, tot_Q, tot_U, n_ord, xp.ones((B,), dtype=bool), resid

    I_prev, Q_prev, U_prev = tot_I.copy(), tot_Q.copy(), tot_U.copy()
    if use_surface:
        mu_pos = atm.mu_pos
        w_pos = atm.gb[:, n_mu + 1:]
        mu_w = mu_factor * mu_pos * w_pos                 # (B, n_mu)
        RII = R_m[:, :, :, 0]; RIQ = R_m[:, :, :, 1]; RIU = R_m[:, :, :, 2]
        RQI = R_m[:, :, :, 3]; RQQ = R_m[:, :, :, 4]; RQU = R_m[:, :, :, 5]
        RUI = R_m[:, :, :, 6]; RUQ = R_m[:, :, :, 7]; RUU = R_m[:, :, :, 8]

    active = xp.ones((B,), dtype=bool)
    for it in range(2, max_iterations + 1):
        # Periodic convergence check (every 8 orders).  Converged cases are
        # masked to zero contribution (am below), so extra orders before the
        # next check add exactly 0 -> bit-identical result.  Removes the
        # per-order device->host sync so the CPU keeps queuing kernels.
        if (it % 8 == 0) and (not bool(active.any())):
            break
        J_I, J_Q, J_U = sos_build_source_pol_batch(
            atm, m, kt, I_prev, Q_prev, U_prev, xpl, xrl, xtl)
        if use_surface:
            Idn = I_prev[:, nt, :n_mu][:, ::-1]           # (B, n_mu)
            Qdn = Q_prev[:, nt, :n_mu][:, ::-1]
            Udn = U_prev[:, nt, :n_mu][:, ::-1]
            wIdn = Idn * mu_w; wQdn = Qdn * mu_w; wUdn = Udn * mu_w
            surf_I = (xp.einsum('bjk,bk->bj', RII, wIdn) +
                      xp.einsum('bjk,bk->bj', RIQ, wQdn) +
                      xp.einsum('bjk,bk->bj', RIU, wUdn))
            surf_Q = (xp.einsum('bjk,bk->bj', RQI, wIdn) +
                      xp.einsum('bjk,bk->bj', RQQ, wQdn) +
                      xp.einsum('bjk,bk->bj', RQU, wUdn))
            surf_U = (xp.einsum('bjk,bk->bj', RUI, wIdn) +
                      xp.einsum('bjk,bk->bj', RUQ, wQdn) +
                      xp.einsum('bjk,bk->bj', RUU, wUdn))
            I_curr, Q_curr, U_curr = integrate_bcs_iqu(
                atm, J_I, J_Q, J_U, surf_I, surf_Q, surf_U)
        else:
            I_curr, Q_curr, U_curr = integrate_bcs_iqu(atm, J_I, J_Q, J_U)
        am = active[:, None, None]
        tot_I = tot_I + xp.where(am, xp.where(xp.isnan(I_curr), 0.0, I_curr), 0.0)
        tot_I[:, :, n_mu] = xp.nan
        tot_Q = tot_Q + xp.where(am, xp.where(xp.isnan(Q_curr), 0.0, Q_curr), 0.0)
        tot_Q[:, :, n_mu] = xp.nan
        tot_U = tot_U + xp.where(am, xp.where(xp.isnan(U_curr), 0.0, U_curr), 0.0)
        tot_U[:, :, n_mu] = xp.nan
        max_in = xp.maximum(xp.maximum(_field_max_abs_batch(I_curr, n_mu),
                                       _field_max_abs_batch(Q_curr, n_mu)),
                            _field_max_abs_batch(U_curr, n_mu))
        max_t = xp.maximum(xp.maximum(_field_max_abs_batch(tot_I, n_mu),
                                      _field_max_abs_batch(tot_Q, n_mu)),
                           _field_max_abs_batch(tot_U, n_mu))
        r = xp.where(max_t > 0.0, max_in / xp.where(max_t > 0.0, max_t, 1.0), 0.0)
        n_ord = xp.where(active, it, n_ord)
        resid = xp.where(active, r, resid)
        newly = active & (r < tolerance)
        conv = conv | newly
        active = active & ~newly
        I_prev, Q_prev, U_prev = I_curr, Q_curr, U_curr

    if surf_seed is not None:
        tot_I = tot_I - sI; tot_I[:, :, n_mu] = xp.nan
        tot_Q = tot_Q - sQ; tot_Q[:, :, n_mu] = xp.nan
        tot_U = tot_U - sU; tot_U[:, :, n_mu] = xp.nan
    return tot_I, tot_Q, tot_U, n_ord, conv, resid


def reconstruct_phi_batch(per_m, m_max, dphi):
    """per_m (B, m_max+1); dphi (B,).  cos reconstruction -> (B,)."""
    base = dphi + M_PI
    I = per_m[:, 0].copy()
    for m in range(1, m_max + 1):
        I = I + 2.0 * per_m[:, m] * xp.cos(m * base)
    return I


def reconstruct_phi_sin_batch(per_m, m_max, dphi):
    base = dphi + M_PI
    U = xp.zeros(per_m.shape[0])
    for m in range(1, m_max + 1):
        U = U + 2.0 * per_m[:, m] * xp.sin(m * base)
    return U


# ---- batched atmosphere orchestrator (R2 coxmunk / R3 aerosol+coxmunk) ------
def _build_batch_medium(atms, n_mu, nt, L_max, greek=False):
    """Stack per-case single-case atm objects into a BatchAtm."""
    B = len(atms)
    bA = BatchAtm()
    bA.B = B; bA.nt = nt; bA.n_mu = n_mu; bA.L_max = L_max; bA.beam_q = 0.0
    a0 = atms[0]
    bA.beta0 = a0.beta0; bA.beta2 = a0.beta2
    bA.gamma2 = a0.gamma2; bA.alpha2 = a0.alpha2
    bA.h = xp.stack([xp.asarray(a.h) for a in atms])
    bA.ch = xp.stack([xp.asarray(a.ch) for a in atms])
    bA.xdel = xp.stack([xp.asarray(a.xdel) for a in atms])
    bA.ydel = xp.stack([xp.asarray(a.ydel) for a in atms])
    bA.rm = xp.stack([xp.asarray(a.rm) for a in atms])
    bA.gb = xp.stack([xp.asarray(a.gb) for a in atms])
    bA.mu_pos = xp.stack([xp.asarray(a.rm[n_mu + 1:]) for a in atms])
    if greek:
        bA.betal = xp.stack([xp.asarray(a.betal) for a in atms])
        bA.gammal = xp.stack([xp.asarray(a.gammal) for a in atms])
        bA.alphal = xp.stack([xp.asarray(a.alphal) for a in atms])
        bA.zetal = xp.stack([xp.asarray(a.zetal) for a in atms])
    return bA


def _rough_fresnel_surf_kernels(bA, mu_suns, m, wind_list, n_mu, nt,
                          sigma_type, n_water, q_convention, n_phi_quad):
    """Batched air-side R_m, R_solar and surf_seed for mode m.  Uses the batched
    Cox-Munk kernel (surface_batch.fourier_kernel_coxmunk_batch) so the whole
    case axis is one GPU pass with no per-case loop.  bA is a BatchAtm
    (ring = bA.mu_pos, heights = bA.h)."""
    from . import surface_batch as SB
    B = bA.B
    dirs = 2 * n_mu + 1
    m_factor = pi if m == 0 else 0.5 * pi
    rings = xp.asarray(bA.mu_pos)                             # (B, n_mu)
    winds = xp.asarray(wind_list, dtype=float)
    mu_suns_arr = xp.asarray(mu_suns, dtype=float)
    # R_m (B, n_mu, n_mu, 9) and R_solar (B, n_mu, 1, 9) in two batched passes
    R_m = SB.fourier_kernel_coxmunk_batch(rings, rings, m, n_phi_quad, winds,
                                          sigma_type, n_water, q_convention,
                                          osoaa_sign_fix=True)
    Rs = SB.fourier_kernel_coxmunk_batch(rings, mu_suns_arr[:, None], m,
                                         n_phi_quad, winds, sigma_type, n_water,
                                         q_convention, osoaa_sign_fix=True)
    # surf_seed: direct solar Fresnel reflection, up-attenuated per level
    h = xp.asarray(bA.h)                                      # (B, nt+1)
    tau_total = h[:, nt][:, None]                             # (B,1)
    trans_sun = xp.exp(-tau_total / mu_suns_arr[:, None])     # (B,1)
    dtau = tau_total - h                                      # (B, nt+1)
    trans = trans_sun[:, :, None] * xp.exp(-dtau[:, :, None] / rings[:, None, :])
    factor = m_factor * mu_suns_arr[:, None, None] * trans    # (B, nt+1, n_mu)
    sI = xp.zeros((B, nt + 1, dirs))
    sQ = xp.zeros((B, nt + 1, dirs))
    sU = xp.zeros((B, nt + 1, dirs))
    sI[:, :, n_mu + 1:] = Rs[:, :, 0, 0][:, None, :] * factor
    sQ[:, :, n_mu + 1:] = Rs[:, :, 0, 3][:, None, :] * factor
    sU[:, :, n_mu + 1:] = Rs[:, :, 0, 6][:, None, :] * factor
    sI[:, :, n_mu] = xp.nan; sQ[:, :, n_mu] = xp.nan; sU[:, :, n_mu] = xp.nan
    return (R_m, (sI, sQ, sU))


def solve_atm_batch(atms, vjs, mu_suns, dphis, wind_list, n_mu, nt,
                    m_max, max_iterations, tolerance, mode='black_fresnel_ocean',
                    L_max=2, sigma_type=1, n_water=1.34, q_convention=1,
                    n_phi_quad=1024, bottom_source=None,
                    use_value_kernel=False, aer_p11=None):
    """Batched atmosphere SOS.  atms: per-case single-case medium objects
    (build_atm_rayleigh for mode 'black_fresnel_ocean', build_atm_aerosol for
    'aerosol_black_fresnel_ocean').  Returns rho_I/Q/U (B,), plus per-m boa I (B, m+1, n_mu)."""
    if mode == 'coxmunk':
        mode = 'black_fresnel_ocean'
    elif mode == 'aerosol_coxmunk':
        mode = 'aerosol_black_fresnel_ocean'
    """Batched atmosphere SOS.  atms is either a list of per-case single-case
    medium objects (build_atm_rayleigh / build_atm_aerosol), which are stacked
    here, or a pre-built BatchAtm (used directly, no per-case build/stack).
    Returns rho_I/Q/U (B,), plus per-m boa I (B, m+1, n_mu)."""
    greek = (mode == 'aerosol_black_fresnel_ocean')
    if isinstance(atms, BatchAtm):
        bA = atms; B = bA.B; atms_list = None
    else:
        B = len(atms); atms_list = atms
        bA = _build_batch_medium(atms, n_mu, nt, max(m_max, L_max), greek=greek)
    dirs = 2 * n_mu + 1
    l_max_ws = max(m_max, L_max)
    zero_kt = tuple(xp.zeros((B, n_mu + 1, dirs)) for _ in range(6))
    mu_factor = (2.0 * pi)   # placeholder; per-m below

    I_pm = xp.zeros((B, m_max + 1))
    Q_pm = xp.zeros((B, m_max + 1))
    U_pm = xp.zeros((B, m_max + 1))
    boa_I = xp.zeros((B, m_max + 1, n_mu))
    boa_Q = xp.zeros((B, m_max + 1, n_mu))
    boa_U = xp.zeros((B, m_max + 1, n_mu))
    n_orders = xp.ones((B,), dtype=int)
    vj_arr = np.asarray(vjs)
    # bottom_source (B, m+1, n_mu) upward per m for atm pass-2
    bsrc = bottom_source is not None
    if bsrc:
        bI, bQ, bU = bottom_source
        bI = xp.asarray(bI); bQ = xp.asarray(bQ); bU = xp.asarray(bU)  # GPU-safe
    mu_pos_b = bA.mu_pos                                # (B, n_mu)
    h_b = bA.h                                          # (B, nt+1)
    tau_tot_b = h_b[:, nt][:, None]                     # (B,1)

    # value-kernel P11 stacked once (shared angle grid, per-case P11 values)
    if greek and use_value_kernel and aer_p11 is not None:
        _aer_th_common = xp.asarray(aer_p11[0][0])
        _aer_P_stack = xp.stack([xp.asarray(aer_p11[b][1]) for b in range(B)])

    for m in range(0, m_max + 1):
        wsb = BatchLegWorkspace(B, n_mu, l_max_ws)
        xpl, xrl, xtl = legendre_compute_pol_batch(wsb, bA.rm, m)
        if greek:
            kernel_phase_fourier_batch(wsb, m, bA.betal)
            kernel_phase_fourier_pol_batch(wsb, m, bA.gammal)
            kernel_phase_fourier_aerosol_full_batch(wsb, m, bA.alphal, bA.zetal)
            pfm = wsb.pfm
            if use_value_kernel and aer_p11 is not None:
                # C aerosol auto-dispatch: I-channel uses the Gibbs-free value
                # kernel (P11 azimuth-direct); Q/U keep the Greek gr/gt/arr/att.
                # Batched over the case axis (shared angle grid, per-case P11).
                from .aerosol import aerosol_value_pfm_batch as _vpfm_batch
                pfm = _vpfm_batch(bA, m, _aer_th_common, _aer_P_stack, 720)
            kt = (pfm, wsb.gr, wsb.gt, wsb.arr, wsb.art, wsb.att)
        else:
            kt = zero_kt
        if bsrc:
            # C3 pass-2: solar primary zeroed; injected water-leaving enters as
            # an up-attenuated initial field at every level.
            prim_i = xp.zeros((B, nt + 1, dirs))
            prim_q = xp.zeros((B, nt + 1, dirs))
            prim_u = xp.zeros((B, nt + 1, dirs))
            prim_i[:, :, n_mu] = xp.nan
            prim_q[:, :, n_mu] = xp.nan
            prim_u[:, :, n_mu] = xp.nan
            # atten[b,i_lev,+k] = exp(-(tau_tot - h)/mu_pos)
            atten = xp.exp(-(tau_tot_b[:, :, None] - h_b[:, :, None]) /
                           mu_pos_b[:, None, :])          # (B, nt+1, n_mu)
            eI = xp.zeros((B, nt + 1, dirs)); eQ = xp.zeros((B, nt + 1, dirs)); eU = xp.zeros((B, nt + 1, dirs))
            eI[:, :, n_mu + 1:] = bI[:, m, :][:, None, :] * atten
            eQ[:, :, n_mu + 1:] = bQ[:, m, :][:, None, :] * atten
            eU[:, :, n_mu + 1:] = bU[:, m, :][:, None, :] * atten
            R_m, _ = _rough_fresnel_surf_kernels(bA, mu_suns, m, wind_list, n_mu, nt,
                                           sigma_type, n_water, q_convention, n_phi_quad)
            mu_fac = (2.0 * pi) if m == 0 else pi
            ti, tq, tu, no, cv, rs = sos_atm_batch(
                bA, m, kt, prim_i, prim_q, prim_u, xpl, xrl, xtl,
                max_iterations, tolerance, R_m=R_m, mu_factor=mu_fac,
                surf_seed=None, ext_init=(eI, eQ, eU))
        else:
            src_i = primary_source_batch(bA, m, xpl, kt[0], ws=wsb, xrl=xrl)
            src_q, src_u = primary_source_pol_batch(bA, m, xrl, xtl, xpl,
                                                    kt[1], kt[2], ws=wsb)
            prim_i, prim_q, prim_u = integrate_bcs_iqu(bA, src_i, src_q, src_u)
            R_m, surf_seed = _rough_fresnel_surf_kernels(
                bA, mu_suns, m, wind_list, n_mu, nt,
                sigma_type, n_water, q_convention, n_phi_quad)
            mu_fac = (2.0 * pi) if m == 0 else pi
            ti, tq, tu, no, cv, rs = sos_atm_batch(
                bA, m, kt, prim_i, prim_q, prim_u, xpl, xrl, xtl,
                max_iterations, tolerance, R_m=R_m, mu_factor=mu_fac,
                surf_seed=surf_seed)
        idx = n_mu + vj_arr
        rows = np.arange(B)
        I_pm[:, m] = ti[rows, 0, idx]
        Q_pm[:, m] = tq[rows, 0, idx]
        U_pm[:, m] = tu[rows, 0, idx]
        boa_I[:, m, :] = ti[:, nt, n_mu - 1::-1]
        boa_Q[:, m, :] = tq[:, nt, n_mu - 1::-1]
        boa_U[:, m, :] = tu[:, nt, n_mu - 1::-1]
        n_orders = xp.maximum(n_orders, no)

    I_TOA = reconstruct_phi_batch(I_pm, m_max, xp.asarray(dphis))
    Q_TOA = reconstruct_phi_batch(Q_pm, m_max, xp.asarray(dphis))
    U_TOA = reconstruct_phi_sin_batch(U_pm, m_max, xp.asarray(dphis))
    mus = xp.asarray(mu_suns)
    rho_I = I_TOA / mus
    rho_Q = Q_TOA / mus
    rho_U = U_TOA / mus
    # T_diff_dn_hemi = (2pi sum mu w boa_I[m=0]) / (pi mu_sun) ; view node w=0
    w_pos_b = bA.gb[:, n_mu + 1:]                       # (B, n_mu)
    F_dn = xp.sum(mu_pos_b * w_pos_b * boa_I[:, 0, :], axis=1) * 2.0 * M_PI
    T_diff_dn_hemi = F_dn / (M_PI * mus)
    return dict(rho_I=rho_I, rho_Q=rho_Q, rho_U=rho_U,
                I_TOA=I_TOA, Q_TOA=Q_TOA, U_TOA=U_TOA,
                boa_I=boa_I, boa_Q=boa_Q, boa_U=boa_U,
                mu_quad=mu_pos_b, w_quad=w_pos_b,
                T_diff_dn_hemi=T_diff_dn_hemi, n_orders=n_orders)


# ---- batched in-water order iteration (internal reflection, sos_pol_intrefl_rough)
def sos_water_intrefl_batch(atm, m, kt, prim_I, prim_Q, prim_U, xpl, xrl, xtl,
                            Rww_K, max_iterations, tolerance):
    """Batched sos_pol_intrefl_rough.  Rww_K (B, n_mu, n_mu, 9) m-mode water
    Cox-Munk internal-reflection kernel.  Top BC at k=0: upwelling field
    reflected downward.  Returns tot_I/Q/U (B, nt+1, dirs), n_ord/conv/resid."""
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    msign = 1.0  # specular water reflection preserves azimuth; no (-1)^m
    tot_I = prim_I.copy(); tot_Q = prim_Q.copy(); tot_U = prim_U.copy()
    n_ord = xp.ones((B,), dtype=int)
    conv = xp.zeros((B,), dtype=bool)
    resid = xp.zeros((B,))
    if max_iterations == 1:
        return tot_I, tot_Q, tot_U, n_ord, xp.ones((B,), dtype=bool), resid
    I_prev, Q_prev, U_prev = prim_I.copy(), prim_Q.copy(), prim_U.copy()
    mu_pos = atm.mu_pos
    w_pos = atm.gb[:, n_mu + 1:]
    az = (2.0 * M_PI) if m == 0 else M_PI
    wk = az * mu_pos * w_pos                          # (B, n_mu)
    K0 = Rww_K[:, :, :, 0]; K1 = Rww_K[:, :, :, 1]; K2 = Rww_K[:, :, :, 2]
    K3 = Rww_K[:, :, :, 3]; K4 = Rww_K[:, :, :, 4]; K5 = Rww_K[:, :, :, 5]
    K6 = Rww_K[:, :, :, 6]; K7 = Rww_K[:, :, :, 7]; K8 = Rww_K[:, :, :, 8]
    active = xp.ones((B,), dtype=bool)
    for it in range(2, max_iterations + 1):
        # Periodic convergence check (every 8 orders) — see atm loop; masked
        # converged cases add exactly 0, so result is bit-identical and the
        # per-order device->host sync is removed.
        if (it % 8 == 0) and (not bool(active.any())):
            break
        J_I, J_Q, J_U = sos_build_source_pol_batch(
            atm, m, kt, I_prev, Q_prev, U_prev, xpl, xrl, xtl)
        I_up = I_prev[:, 0, n_mu + 1:]                # (B, n_mu) upwelling at k=0
        Q_up = Q_prev[:, 0, n_mu + 1:]
        U_up = U_prev[:, 0, n_mu + 1:]
        wIu = I_up * wk; wQu = Q_up * wk; wUu = U_up * wk
        tb_I = msign * (xp.einsum('bjk,bk->bj', K0, wIu) +
                        xp.einsum('bjk,bk->bj', K1, wQu) +
                        xp.einsum('bjk,bk->bj', K2, wUu))
        tb_Q = msign * (xp.einsum('bjk,bk->bj', K3, wIu) +
                        xp.einsum('bjk,bk->bj', K4, wQu) +
                        xp.einsum('bjk,bk->bj', K5, wUu))
        tb_U = msign * (xp.einsum('bjk,bk->bj', K6, wIu) +
                        xp.einsum('bjk,bk->bj', K7, wQu) +
                        xp.einsum('bjk,bk->bj', K8, wUu))
        I_curr, Q_curr, U_curr = integrate_bcs_iqu(
            atm, J_I, J_Q, J_U, tb_I, tb_Q, tb_U, top=True)
        am = active[:, None, None]
        tot_I = tot_I + xp.where(am, xp.where(xp.isnan(I_curr), 0.0, I_curr), 0.0)
        tot_I[:, :, n_mu] = xp.nan
        tot_Q = tot_Q + xp.where(am, xp.where(xp.isnan(Q_curr), 0.0, Q_curr), 0.0)
        tot_Q[:, :, n_mu] = xp.nan
        tot_U = tot_U + xp.where(am, xp.where(xp.isnan(U_curr), 0.0, U_curr), 0.0)
        tot_U[:, :, n_mu] = xp.nan
        max_in = xp.maximum(xp.maximum(_field_max_abs_batch(I_curr, n_mu),
                                       _field_max_abs_batch(Q_curr, n_mu)),
                            _field_max_abs_batch(U_curr, n_mu))
        max_t = xp.maximum(xp.maximum(_field_max_abs_batch(tot_I, n_mu),
                                      _field_max_abs_batch(tot_Q, n_mu)),
                           _field_max_abs_batch(tot_U, n_mu))
        r = xp.where(max_t > 0.0, max_in / xp.where(max_t > 0.0, max_t, 1.0), 0.0)
        n_ord = xp.where(active, it, n_ord)
        resid = xp.where(active, r, resid)
        newly = active & (r < tolerance)
        conv = conv | newly
        active = active & ~newly
        I_prev, Q_prev, U_prev = I_curr, Q_curr, U_curr
    return tot_I, tot_Q, tot_U, n_ord, conv, resid


# ---- batched air->water forward coupling (couple_atm_to_water, wind>0) ------
def couple_atm_to_water_batch(boa_I, boa_Q, boa_U, mu_atm_list, m_max,
                              n_water, q_convention, mu_water_pos_list,
                              w_atm_list, wind_list, sigma_type,
                              n_phi_kernel=128):
    """Batched rt_air_water_couple_atm_to_water (Cox-Munk Fourier, wind>0).
    boa_{I,Q,U} (B, m+1, n_a); mu_atm_list/w_atm_list/mu_water_pos_list are
    per-case lists; wind_list (B,).  Returns cI/cQ/cU (B, m+1, n_w)."""
    from . import surface_batch as SB
    B = len(mu_atm_list)
    n_w = len(mu_water_pos_list[0])
    cI = xp.zeros((B, m_max + 1, n_w))
    cQ = xp.zeros((B, m_max + 1, n_w))
    cU = xp.zeros((B, m_max + 1, n_w))
    mu_a = xp.stack([xp.asarray(x, dtype=float) for x in mu_atm_list])   # (B, n_a)
    w_a = xp.stack([xp.asarray(x, dtype=float) for x in w_atm_list])
    mu_w_stack = xp.stack([xp.asarray(x, dtype=float) for x in mu_water_pos_list])  # (B, n_w)
    winds_arr = xp.asarray(wind_list, dtype=float)
    for m in range(m_max + 1):
        az = (2.0 * M_PI) if m == 0 else M_PI
        # batched air->water T_aw kernel over the whole group (no per-case loop)
        Tm = SB.fourier_kernel_T_aw_batch(mu_w_stack, mu_a, m, n_phi_kernel,
                                          winds_arr, sigma_type, n_water,
                                          q_convention)                 # (B,n_w,n_a,9)
        wj = az * mu_a * w_a                            # (B, n_a)
        Ia = boa_I[:, m, :]; Qa = boa_Q[:, m, :]; Ua = boa_U[:, m, :]
        wIa = wj * Ia; wQa = wj * Qa; wUa = wj * Ua
        cI[:, m, :] = (xp.einsum('bij,bj->bi', Tm[..., 0], wIa) +
                       xp.einsum('bij,bj->bi', Tm[..., 1], wQa) +
                       xp.einsum('bij,bj->bi', Tm[..., 2], wUa))
        cQ[:, m, :] = (xp.einsum('bij,bj->bi', Tm[..., 3], wIa) +
                       xp.einsum('bij,bj->bi', Tm[..., 4], wQa) +
                       xp.einsum('bij,bj->bi', Tm[..., 5], wUa))
        cU[:, m, :] = (xp.einsum('bij,bj->bi', Tm[..., 6], wIa) +
                       xp.einsum('bij,bj->bi', Tm[..., 7], wQa) +
                       xp.einsum('bij,bj->bi', Tm[..., 8], wUa))
    return cI, cQ, cU


# ---- batched water->atm reverse coupling (couple_water_to_atm) --------------
def _aw_interp_batch(mu_grid, vals, targets):
    """Batched aw_interp_on_unsorted for ascending, no-dup rings (N<=64 after
    the C clamp).  mu_grid (B,N) ascending; vals (B,M2,N); targets (B,M).
    Linear interp with linear extrapolation at both ends, matching
    _linear_interp_ascending exactly.  Returns (B,M2,M).  The interp bracket
    depends only on (B,targets), so vals is gathered along its ring axis."""
    B, N = mu_grid.shape
    if N > 64:                                     # C: defensive clamp to first 64
        mu_grid = mu_grid[:, :64]
        vals = vals[:, :, :64]
        N = 64
    M = targets.shape[1]
    M2 = vals.shape[1]
    c = (mu_grid[:, :, None] < targets[:, None, :]).sum(axis=1)   # (B,M)
    hi = xp.clip(c, 1, N - 1)
    lo = hi - 1
    mu0 = xp.take_along_axis(mu_grid, lo, axis=1)                 # (B,M)
    mu1 = xp.take_along_axis(mu_grid, hi, axis=1)
    t = (targets - mu0) / (mu1 - mu0)                            # (B,M)
    lo_bc = xp.ascontiguousarray(xp.broadcast_to(lo[:, None, :], (B, M2, M)))
    hi_bc = xp.ascontiguousarray(xp.broadcast_to(hi[:, None, :], (B, M2, M)))
    v0 = xp.take_along_axis(vals, lo_bc, axis=2)                 # (B,M2,M)
    v1 = xp.take_along_axis(vals, hi_bc, axis=2)
    return (1.0 - t)[:, None, :] * v0 + t[:, None, :] * v1


def _diffuse_interp_batch(ext_mu, ext_vals, targets):
    """Clamped linear interp (t in [0,1]) of ext_vals (B,N) sampled on ext_mu
    (B,N ascending) onto targets (B,M).  Matches _add_diffuse_top_primary:
    ring nodes coinciding with an ext_mu node fall out as t=0/1 (exact), and
    off-grid nodes use clamped interpolation (no extrapolation).  Returns (B,M)."""
    B, N = ext_mu.shape
    cnt = (ext_mu[:, :, None] <= targets[:, None, :]).sum(axis=1)   # (B,M) <= (non-strict)
    lo = xp.clip(cnt - 1, 0, N - 2)
    hi = lo + 1
    mu0 = xp.take_along_axis(ext_mu, lo, axis=1)
    mu1 = xp.take_along_axis(ext_mu, hi, axis=1)
    denom = xp.where(mu1 != mu0, mu1 - mu0, 1.0)
    t = xp.clip((targets - mu0) / denom, 0.0, 1.0)
    v0 = xp.take_along_axis(ext_vals, lo, axis=1)
    v1 = xp.take_along_axis(ext_vals, hi, axis=1)
    return (1.0 - t) * v0 + t * v1


def couple_water_to_atm_batch(I_w, Q_w, U_w, mu_water_pos_list, m_max,
                              n_water, q_convention, mu_atm_pos_list):
    """Batched rt_air_water_couple_water_to_atm (flat reverse Snell + direct
    T_wa Mueller + aw_interp), fully on the active backend (no host round-trip).
    I/Q/U_w (B, m+1, n_w); per-case ring lists (uniform sizes within a group).
    Returns wlI/wlQ/wlU (B, m+1, n_a)."""
    from . import surface_batch as SB
    mu_wr = xp.stack([xp.asarray(x, dtype=float) for x in mu_water_pos_list])  # (B,n_w)
    mu_a = xp.stack([xp.asarray(x, dtype=float) for x in mu_atm_pos_list])     # (B,n_a)
    n2 = n_water * n_water
    sin2_air = 1.0 - mu_a * mu_a
    mu_w = xp.sqrt(xp.maximum(0.0, 1.0 - sin2_air / n2))          # (B,n_a) water-side
    MT = SB._fresnel_T_mueller_vec(mu_w, n_water, 1.0, q_convention)  # 9 x (B,n_a)
    Iint = _aw_interp_batch(mu_wr, I_w, mu_w)                     # (B,m+1,n_a)
    Qint = _aw_interp_batch(mu_wr, Q_w, mu_w)
    Uint = _aw_interp_batch(mu_wr, U_w, mu_w)
    T0 = MT[0][:, None, :]; T1 = MT[1][:, None, :]; T2 = MT[2][:, None, :]
    T3 = MT[3][:, None, :]; T4 = MT[4][:, None, :]; T5 = MT[5][:, None, :]
    T6 = MT[6][:, None, :]; T7 = MT[7][:, None, :]; T8 = MT[8][:, None, :]
    wlI = T0 * Iint + T1 * Qint + T2 * Uint
    wlQ = T3 * Iint + T4 * Qint + T5 * Uint
    wlU = T6 * Iint + T7 * Qint + T8 * Uint
    return wlI, wlQ, wlU


# ---- batched diffuse-top primary injection (_add_diffuse_top_primary) -------
def add_diffuse_top_primary_batch(atm, m, ext_I, ext_Q, ext_U, ext_mu_list,
                                  kt, xpl, xrl, xtl, src_i, src_q, src_u):
    """Batched ocrt_add_diffuse_top_primary (DT-HIOM TABLE kernel).  Adds the
    first-order scattering source of the injected downwelling diffuse field to
    src_{i,q,u} (B, nt+1, dirs) in place.  ext_{I,Q,U} (B, n_ext); ext_mu_list
    per-case ext sample mu (== ring positive nodes here, exact match)."""
    B, nt, n_mu = atm.B, atm.nt, atm.n_mu
    dirs = 2 * n_mu + 1
    ray_on = (m <= 2)
    beta0_m = atm.beta0 if m == 0 else 0.0
    beta2, gamma2, alpha2 = atm.beta2, atm.gamma2, atm.alpha2
    P, GR, GT, AR, AT, ATT = kt
    gb = atm.gb                                        # (B, dirs)
    jp = np.arange(1, n_mu + 1)
    pk = n_mu + jp; mk = n_mu - jp

    # AI[b,c] via per-case ext_mu match (exact) or clamped linear interp,
    # replicating _add_diffuse_top_primary (ext_mu may be a subgrid of the ring).
    # AI[b,c] = 2*gb[b,c] * ext_interp(ext_mu -> ring node c), clamped linear
    # interp.  Ring nodes coinciding with an ext_mu node (all GL nodes here, as
    # ext_mu is the GL subgrid) fall out exactly; inserted special nodes use the
    # clamped interpolation.  Fully on the backend (no per-mode host round-trip).
    ring_pos = atm.rm[:, n_mu + 1:]                    # (B, n_mu) target ring nodes
    ext_mu_b = xp.stack([xp.asarray(e, dtype=float) for e in ext_mu_list])  # (B, n_ext)
    g2 = 2.0 * gb[:, n_mu + 1:]                        # (B, n_mu)
    AI = g2 * _diffuse_interp_batch(ext_mu_b, ext_I, ring_pos)
    AQ = (g2 * _diffuse_interp_batch(ext_mu_b, ext_Q, ring_pos)
          if ext_Q is not None else xp.zeros((B, n_mu)))
    AU = (g2 * _diffuse_interp_batch(ext_mu_b, ext_U, ring_pos)
          if ext_U is not None else xp.zeros((B, n_mu)))

    # per incident column c=1..n_mu: replicate the single-case row construction,
    # vectorized over batch b and layers.
    for cc_i, c in enumerate(range(1, n_mu + 1)):
        col_p = n_mu + c
        col_m = n_mu - c
        # kernel rows (B, dirs) for this incident column c
        rI = xp.zeros((B, dirs)); rQ = xp.zeros((B, dirs)); rU = xp.zeros((B, dirs))
        rIq = xp.zeros((B, dirs)); rQq = xp.zeros((B, dirs)); rUq = xp.zeros((B, dirs))
        rIu = xp.zeros((B, dirs)); rQu = xp.zeros((B, dirs)); rUu = xp.zeros((B, dirs))
        for j in range(1, n_mu + 1):
            # +j
            ip = n_mu + j
            rI[:, ip] = P[:, c, n_mu - j]
            rQ[:, ip] = GR[:, c, n_mu - j]
            rU[:, ip] = -GT[:, c, n_mu - j]
            rIq[:, ip] = GR[:, j, n_mu - c]
            rQq[:, ip] = AR[:, c, n_mu - j]
            rUq[:, ip] = AT[:, j, n_mu - c]
            rIu[:, ip] = -GT[:, j, n_mu - c]
            rQu[:, ip] = AT[:, c, n_mu - j]
            rUu[:, ip] = -ATT[:, c, n_mu - j]
            # -j (a=j)
            im = n_mu - j
            rI[:, im] = P[:, c, n_mu + j]
            rQ[:, im] = GR[:, c, n_mu + j]
            rU[:, im] = -GT[:, c, n_mu + j]
            rIq[:, im] = GR[:, j, n_mu + c]
            rQq[:, im] = AR[:, c, n_mu + j]
            rUq[:, im] = -AT[:, j, n_mu + c]
            rIu[:, im] = GT[:, j, n_mu + c]
            rQu[:, im] = AT[:, c, n_mu + j]
            rUu[:, im] = -ATT[:, c, n_mu + j]
        # Rayleigh rows (B, dirs)
        xplc = xpl[:, col_m][:, None]                  # xpl[-c]
        xrlc = xrl[:, col_m][:, None]
        xtlc = xtl[:, col_m][:, None]
        raI = beta0_m + beta2 * xpl * xplc
        raQ = (gamma2 * xrl * xplc) if ray_on else xp.zeros((B, dirs))
        raU = (gamma2 * xtl * xplc) if ray_on else xp.zeros((B, dirs))
        rqI = (gamma2 * xpl * xrlc) if ray_on else xp.zeros((B, dirs))
        rqQ = (alpha2 * xrl * xrlc) if ray_on else xp.zeros((B, dirs))
        rqU = (alpha2 * xtl * xrlc) if ray_on else xp.zeros((B, dirs))
        ruI = (gamma2 * xpl * xtlc) if ray_on else xp.zeros((B, dirs))
        ruQ = (alpha2 * xrl * xtlc) if ray_on else xp.zeros((B, dirs))
        ruU = (alpha2 * xtl * xtlc) if ray_on else xp.zeros((B, dirs))
        # poison solar slot
        for arr_ in (rI, rQ, rU, rIq, rQq, rUq, rIu, rQu, rUu,
                     raI, raQ, raU, rqI, rqQ, rqU, ruI, ruQ, ruU):
            arr_[:, n_mu] = 0.0
        AIc = AI[:, cc_i][:, None]                     # (B,1)
        AQc = AQ[:, cc_i][:, None]
        AUc = AU[:, cc_i][:, None]
        mu_c = atm.mu_pos[:, cc_i][:, None]            # (B,1)
        # per-layer ch_c (B, nt+1)
        ch_c = 0.5 * xp.exp(-atm.h / mu_c)             # (B, nt+1)
        xd = atm.xdel; yd = atm.ydel                   # (B, nt+1)
        # contribution to src (B, nt+1, dirs)
        contribI = (AIc[:, :, None] * (xd[:, :, None] * rI[:, None, :] + yd[:, :, None] * raI[:, None, :])
                    + AQc[:, :, None] * (xd[:, :, None] * rIq[:, None, :] + yd[:, :, None] * rqI[:, None, :])
                    + AUc[:, :, None] * (xd[:, :, None] * rIu[:, None, :] + yd[:, :, None] * ruI[:, None, :]))
        contribQ = (AIc[:, :, None] * (xd[:, :, None] * rQ[:, None, :] + yd[:, :, None] * raQ[:, None, :])
                    + AQc[:, :, None] * (xd[:, :, None] * rQq[:, None, :] + yd[:, :, None] * rqQ[:, None, :])
                    + AUc[:, :, None] * (xd[:, :, None] * rQu[:, None, :] + yd[:, :, None] * ruQ[:, None, :]))
        contribU = (AIc[:, :, None] * (xd[:, :, None] * rU[:, None, :] + yd[:, :, None] * raU[:, None, :])
                    + AQc[:, :, None] * (xd[:, :, None] * rUq[:, None, :] + yd[:, :, None] * rqU[:, None, :])
                    + AUc[:, :, None] * (xd[:, :, None] * rUu[:, None, :] + yd[:, :, None] * ruU[:, None, :]))
        src_i += ch_c[:, :, None] * contribI
        src_q += ch_c[:, :, None] * contribQ
        src_u += -ch_c[:, :, None] * contribU
    src_i[:, :, n_mu] = xp.nan
    src_q[:, :, n_mu] = xp.nan
    src_u[:, :, n_mu] = xp.nan
    return src_i, src_q, src_u


# ---- batched water solve (constituent, ext_top diffuse-top) for R1 ----------
def solve_water_batch_r1(atms_w, mu_sun_water_list, mu_view_water_list,
                         Rww_b, m_max, max_it_water, tol_water,
                         F_sun_water_list, ext_top=None,
                         raa_list=None, n_water=1.34, q_convention=1,
                         ext_list=None, ext_w_list=None):
    """Batched constituent water SOS for R1.  F_sun_water_list (B,) scales the
    exported fields by f_scale = F_sun_water/pi (matching sos_pure_fixed_bulk).
    If raa_list is given, also reconstructs the already-solved I/Q/U fields
    at 0- and transmits the same node fields through T_wa to 0+.  This adds
    only vectorized Fourier reconstruction and six scalar normalizations in
    the caller; no SOS order, kernel build, quadrature, or coupling pass is
    added."""
    bA = atms_w
    B, nt, n_mu = bA.B, bA.nt, bA.n_mu
    dirs = 2 * n_mu + 1
    Iup = xp.zeros((B, m_max + 1, n_mu))
    Qup = xp.zeros((B, m_max + 1, n_mu))
    Uup = xp.zeros((B, m_max + 1, n_mu))
    Iview = xp.zeros((B, m_max + 1))
    Qview = xp.zeros((B, m_max + 1))
    Uview = xp.zeros((B, m_max + 1))
    n_orders = xp.ones((B,), dtype=int)
    # view node index per case (exact match to mu_view_water), vectorized
    mu_view = xp.asarray(mu_view_water_list)
    pos = bA.rm[:, n_mu + 1:]                            # (B, n_mu) positive nodes
    match = xp.abs(pos - mu_view[:, None]) < 1e-13       # (B, n_mu)
    has = match.any(axis=1)
    jex = xp.where(has, xp.argmax(match, axis=1) + 1, -1)   # (B,) 1-based ring idx or -1
    if ext_top is not None:
        dtI, dtQ, dtU, ext_mu_list = ext_top
    # Fourier-mode early exit (commit #14) per case
    m_amp_max = xp.zeros((B,))
    streak = xp.zeros((B,), dtype=int)
    exited = xp.zeros((B,), dtype=bool)
    rows = xp.arange(B)
    ti_m0 = None                                          # m=0 field for Kd(0-)
    for m in range(m_max + 1):
        wsb = BatchLegWorkspace(B, n_mu, max(bA.L_max, 2))
        xpl, xrl, xtl = legendre_compute_pol_batch(wsb, bA.rm, m)
        kernel_phase_fourier_batch(wsb, m, bA.betal)
        kernel_phase_fourier_pol_batch(wsb, m, bA.gammal)
        kernel_phase_fourier_aerosol_full_batch(wsb, m, bA.alphal, bA.zetal)
        kt = (wsb.pfm, wsb.gr, wsb.gt, wsb.arr, wsb.art, wsb.att)
        src_i = primary_source_batch(bA, m, xpl, kt[0], ws=wsb, xrl=xrl)
        src_q, src_u = primary_source_pol_batch(bA, m, xrl, xtl, xpl, kt[1], kt[2], ws=wsb)
        if ext_top is not None and m <= m_max:
            add_diffuse_top_primary_batch(
                bA, m, dtI[:, m, :], dtQ[:, m, :], dtU[:, m, :], ext_mu_list,
                kt, xpl, xrl, xtl, src_i, src_q, src_u)
        prim_i, prim_q, prim_u = integrate_bcs_iqu(bA, src_i, src_q, src_u)
        ti, tq, tu, no, cv, rs = sos_water_intrefl_batch(
            bA, m, kt, prim_i, prim_q, prim_u, xpl, xrl, xtl,
            Rww_b[m], max_it_water, tol_water)
        if m == 0:
            ti_m0 = ti
        # 7e view extraction (only for non-exited cases; exited stay 0)
        iv = ti[rows, 0, n_mu + jex]
        qv = tq[rows, 0, n_mu + jex]
        uv = tu[rows, 0, n_mu + jex]
        keep = ~exited
        Iview[:, m] = xp.where(keep, iv, 0.0)
        Qview[:, m] = xp.where(keep, qv, 0.0)
        Uview[:, m] = xp.where(keep, uv, 0.0)
        # early-exit decision (commit #14): view amp, streak
        amp = xp.maximum(xp.maximum(xp.abs(iv), xp.abs(qv)), xp.abs(uv))
        m_amp_max = xp.where(keep, xp.maximum(m_amp_max, amp), m_amp_max)
        cond = keep & (m >= 3) & (m_amp_max > 0.0)
        small = cond & (amp < tol_water * m_amp_max)
        streak = xp.where(small, streak + 1, xp.where(cond & ~small, 0, streak))
        newly_exit = cond & (streak >= 2)
        # 7e-bis node store: only cases not exited AND not exiting this mode
        store = (~exited) & (~newly_exit)
        sm = store[:, None]
        Iup[:, m, :] = xp.where(sm, ti[:, 0, n_mu + 1:], 0.0)
        Qup[:, m, :] = xp.where(sm, tq[:, 0, n_mu + 1:], 0.0)
        Uup[:, m, :] = xp.where(sm, tu[:, 0, n_mu + 1:], 0.0)
        n_orders = xp.where(keep, xp.maximum(n_orders, no), n_orders)
        exited = exited | newly_exit
    fsc = (xp.asarray(F_sun_water_list) / M_PI)[:, None, None]   # (B,1,1)
    fscv = (xp.asarray(F_sun_water_list) / M_PI)[:, None]        # (B,1)
    out = dict(I_up_per_m=Iup * fsc, Q_up_per_m=Qup * fsc, U_up_per_m=Uup * fsc,
               I_view=Iview * fscv, Q_view=Qview * fscv, U_view=Uview * fscv,
               jex=jex, n_orders=n_orders)
    # ---- Kd(0-), Ed(0-) : m=0 하향조도의 표층 로그미분 ---------------------------
    # Ed(z) = 직달빔 + 확산장(ti) 반구조도 + 하늘광 미산란 투과분.
    # 하늘광(확산-top dtI)의 미산란 투과분은 SOS 장 ti 에 들어있지 않다
    # (add_diffuse_top 은 하늘광의 산란 소스만 넣는다). 이 항을 빠뜨리면 근표층
    # Ed 가 비물리적으로 증가해 Kd 가 음수가 된다. 미산란 투과분을 명시적으로
    # 더하면 Ed 프로파일이 단조감소하고 Kd 가 물리적(양수, ~a+bb)이 된다.
    Kd0 = xp.zeros((B,)); Ed0m = xp.zeros((B,))
    if ti_m0 is not None and ext_list is not None:
        mu_pos = bA.mu_pos                                 # (B, n_mu) 양의 링 절점
        w_ring = bA.gb[:, n_mu + 1:]                       # (B, n_mu) 양의 링 가중치
        musw = bA.mu_sun                                   # (B,) 수중 태양 코사인
        Fsw = xp.asarray(F_sun_water_list)                # (B,) 수면 직하 빔
        fscl = Fsw / M_PI                                  # 확산장 스케일(내보내기와 동일)

        def _ed_diff(lvl):                                 # ti 확산장 하향 반구조도
            D = ti_m0[:, lvl, :n_mu][:, ::-1]             # (B, n_mu) abs 하향
            return 2.0 * M_PI * fscl * xp.sum(D * mu_pos * w_ring, axis=1)

        # 하늘광 미산란 투과분: dtI(mu_c)*exp(-h/mu_c) 를 GL 절점·가중치로 반구적분
        if ext_top is not None and ext_w_list is not None:
            emu = xp.asarray(ext_top[3][0])               # (n_ext,) = mu_wg (GL 절점)
            ewt = xp.asarray(ext_w_list[0])               # (n_ext,) = w_wg (GL 가중치)
            dtI0 = ext_top[0][:, 0, :]                     # (B, n_ext) m=0 하늘광 라디언스

            def _ed_sky(lvl):
                att = xp.exp(-bA.h[:, lvl:lvl + 1] / emu[None, :])         # (B, n_ext)
                return 2.0 * M_PI * fscl * xp.sum(
                    dtI0 * att * emu[None, :] * ewt[None, :], axis=1)      # (B,)
        else:
            def _ed_sky(lvl):
                return xp.zeros((B,))

        h1 = bA.h[:, 1]                                    # level 1 광학깊이 (h[0]=0)
        Ed0 = Fsw * musw + _ed_diff(0) + _ed_sky(0)        # 직달+확산+하늘광 미산란
        Ed1 = Fsw * musw * xp.exp(-h1 / musw) + _ed_diff(1) + _ed_sky(1)
        extv = xp.asarray(ext_list)                       # (B,) a+b
        z1 = xp.where(extv > 0.0, h1 / xp.where(extv > 0.0, extv, 1.0), 0.0)
        good = (Ed0 > 0.0) & (Ed1 > 0.0) & (z1 > 0.0)
        Kd0 = xp.where(good, -xp.log(xp.where(good, Ed1, 1.0)
                                     / xp.where(good, Ed0, 1.0))
                             / xp.where(z1 > 0.0, z1, 1.0), 0.0)
        Ed0m = Ed0
    out['Kd0minus'] = Kd0
    out['Ed_0minus'] = Ed0m
    # ---- Lu_0plus (I_air): node interp at mu_w -> T_wa -> reconstruct -------
    if raa_list is not None:
        from .water import raa_to_water_view_phi
        from . import surface_batch as SB
        mu_ring = bA.rm[:, n_mu + 1:]                       # (B, n_w) xp
        mu_w = xp.asarray([mu_view_water_list[b] for b in range(B)])   # (B,)
        M = SB._fresnel_T_mueller_vec(mu_w, n_water, 1.0, q_convention)  # 9 x (B,)
        tgt = mu_w[:, None]                                 # (B,1)
        Iwm = _aw_interp_batch(mu_ring, Iup, tgt)[:, :, 0]  # (B, m+1) raw (pre f_scale)
        Qwm = _aw_interp_batch(mu_ring, Qup, tgt)[:, :, 0]
        Uwm = _aw_interp_batch(mu_ring, Uup, tgt)[:, :, 0]
        Ia = M[0][:, None] * Iwm + M[1][:, None] * Qwm + M[2][:, None] * Uwm  # (B,m+1)
        Qa = M[3][:, None] * Iwm + M[4][:, None] * Qwm + M[5][:, None] * Uwm
        Ua = M[6][:, None] * Iwm + M[7][:, None] * Qwm + M[8][:, None] * Uwm
        phi_v = xp.asarray([raa_to_water_view_phi(raa_list[b]) for b in range(B)])
        pa = phi_v + M_PI                                   # (B,)
        m_arr = xp.arange(m_max + 1)                        # (m+1,)
        fac = xp.where(m_arr == 0, 1.0, 2.0)               # (m+1,)
        cos_mp = xp.cos(m_arr[None, :] * pa[:, None])       # (B, m+1)
        sin_mp = xp.sin(m_arr[None, :] * pa[:, None])
        fs = xp.asarray(F_sun_water_list) / M_PI            # (B,)
        # Just below the surface: use the exact view-node per-mode fields that
        # were already retained for convergence and output extraction.
        out['Lu_0minus'] = xp.sum(fac[None, :] * Iview * cos_mp, axis=1) * fs
        out['Qu_0minus'] = xp.sum(fac[None, :] * Qview * cos_mp, axis=1) * fs
        out['Uu_0minus'] = xp.sum(-fac[None, :] * Uview * sin_mp, axis=1) * fs
        # Just above the surface: same solved node field after T_wa.
        out['Lu_0plus'] = xp.sum(fac[None, :] * Ia * cos_mp, axis=1) * fs
        out['Qu_0plus'] = xp.sum(fac[None, :] * Qa * cos_mp, axis=1) * fs
        out['Uu_0plus'] = xp.sum(-fac[None, :] * Ua * sin_mp, axis=1) * fs
    return out
