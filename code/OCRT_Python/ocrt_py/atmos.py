"""Atmospheric Rayleigh SOS (Phase A: black surface; Phase B: Cox-Munk ocean, rho_R). Gas absorption optional.

Port of the atmosphere branch of rt_solve_case_pol_impl (rt_solver.c) for the
Rayleigh-only, black-surface case.  Reuses the shared SOS engine already ported
for the in-water side (kernel.py legendre/phase-fourier, sos.py primary sources,
integrate_bcs, sos_build_source_pol, sos_atm_black, reconstruct_phi).

Only the atmosphere BUILD differs from water:
  - depol = 0.0279 (6SV fixed), ron = 2(1-d)/(2+d)
  - Greek coeffs beta0=1, beta2=ron/2, gamma2=-ron*sqrt(1.5), alpha2=3*ron
  - layer profile xdel=0 (no aerosol), ydel=1 (Rayleigh active)
  - solar/direct beam is unpolarized: beam_q = 0
  - ring = GL(n_mu) + zero-weight view node (sampled at TOA only)

C reference (production defaults, Rayleigh-only): n_mu=24, n_layers=40,
fourier_m_max=2, sos.max_iterations=20, tolerance=1e-7, LINEAR integration.
TOA output: rho = I_TOA / mu_solar (F=pi internal), dphi = raa*pi/180,
reconstruct base = dphi + pi.
"""
import numpy as np
from math import cos, pi, exp, sqrt

from . import kernel as kn
from . import sos as S
from . import molprofile as mp
from .rayleigh import tau_rayleigh_bodhaine
from .backend import xp
from .spectral_contract import require_wavelength

DEPOL = 0.0279


def build_atm_rayleigh(nt, tau_R, mu_sun, mu_view, n_mu_gl):
    """Build the Rayleigh atmosphere medium (rt_atm_build_rayleigh) with a
    zero-weight view node appended to the GL ring.  Returns (atm, vj) where
    vj is the POSITIVE ring index of the view node."""
    gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_gl)
    u = S.UAngles(list(gl_mu), list(gl_w))
    u.add(mu_view)
    mu_ring = xp.asarray(u.mu)
    w_ring = xp.asarray(u.w)
    n_mu = len(mu_ring)

    atm = S.Atm(nt, n_mu)
    ron = 2.0 * (1.0 - DEPOL) / (2.0 + DEPOL)
    atm.mu_sun = mu_sun
    atm.beta0 = 1.0
    atm.beta2 = 0.5 * ron
    atm.gamma2 = -ron * sqrt(1.5)
    atm.alpha2 = 3.0 * ron
    atm.beam_q = 0.0                 # unpolarized solar beam

    k_arr = xp.arange(nt + 1, dtype=float)
    atm.h = k_arr * tau_R / nt
    atm.ch = 0.5 * xp.exp(-atm.h / mu_sun)
    atm.xdel = xp.zeros(nt + 1)      # no aerosol
    atm.ydel = xp.ones(nt + 1)       # Rayleigh active

    atm.rm[n_mu] = -mu_sun           # solar slot j=0
    atm.gb[n_mu] = 0.0
    atm.rm[n_mu + 1:] = mu_ring
    atm.rm[:n_mu] = -xp.flip(mu_ring)
    atm.gb[n_mu + 1:] = w_ring
    atm.gb[:n_mu] = xp.flip(w_ring)

    # z_km_level from the US Standard Atmosphere 1962 cumulative molecular
    # column (rt_atm_build_rayleigh, OCRT v1.2):
    #   h[k]/tau_R = F_mol_above(z_k),  F normalized to 1 at sea level and
    #   0 at the US62 model top.  k=0 retains the 100 km AFGL integration
    #   cap so absorbing-gas columns above the US62 top stay in layer 0.
    z_cap = 100.0
    z_km_level = xp.zeros(nt + 1)
    z_km_level[1:] = mp.us62_altitude_from_grid_fraction(atm.h[1:] / tau_R)
    z_km_level[0] = z_cap

    vj = int(xp.argmin(xp.abs(mu_ring - mu_view))) + 1   # positive index
    return atm, vj, z_km_level


def apply_gas_absorption(atm, absorption, wl_nm, z_km_level):
    """rt_atm_apply_gas_absorption: add per-layer gas absorption OD to the
    Rayleigh atmosphere.  Modifies atm.h, atm.ch, atm.xdel, atm.ydel in place
    and returns tau_abs_total.  z_km_level from build_atm_rayleigh.  Array-safe
    (layer axis vectorized; gas axis is a short loop over N_GAS)."""
    from .absorption import (afgl_cumulative_column_vec,
                             xsec_interp_layer_z_vec, N_GAS)
    nt = atm.n_layers
    ab = absorption
    at = ab.atm
    z_km_level = xp.asarray(z_km_level, dtype=float)

    # Step 1: per-layer absorption OD over the SOS layer bands (layer axis
    # vectorized; z_lo<=z_hi enforced by min/max)
    z_hi = z_km_level[:nt]
    z_lo = z_km_level[1:nt + 1]
    zlo_s = xp.minimum(z_lo, z_hi)
    zhi_s = xp.maximum(z_lo, z_hi)
    z_mid = 0.5 * (zlo_s + zhi_s)
    tau_abs_layer = xp.zeros(nt)
    for g in range(N_GAS):
        xs = ab.xsec[g]
        if not xs.available:
            continue
        Ndef = at.column_default[g]
        N_layer = (afgl_cumulative_column_vec(at, g, zlo_s)
                   - afgl_cumulative_column_vec(at, g, zhi_s))
        if ab.overridden[g] and Ndef > 0.0:
            N = ab.column_eff[g] * (N_layer / Ndef)
        else:
            N = N_layer
        sigma = xsec_interp_layer_z_vec(xs, at.z_km, z_mid, wl_nm)
        tau_abs_layer = tau_abs_layer + sigma * N

    # Step 2: rebuild h, ch cumulating gas absorption
    h_old = atm.h.copy()
    cum_abs = xp.concatenate([xp.zeros(1), xp.cumsum(tau_abs_layer)])
    atm.h = h_old + cum_abs
    atm.ch = 0.5 * xp.exp(-atm.h / atm.mu_sun)

    # Step 3: re-normalize xdel, ydel (scattering ODs unchanged, denom grows)
    dt_old = h_old[1:] - h_old[:-1]
    dt_new = dt_old + tau_abs_layer
    ratio = xp.where((dt_new > 0.0) & (dt_old > 0.0), dt_old / dt_new, 1.0)
    atm.xdel[1:] = atm.xdel[1:] * ratio
    atm.ydel[1:] = atm.ydel[1:] * ratio

    atm.tau_total = atm.h[nt]
    return float(tau_abs_layer.sum())


def solve_atm_rayleigh_black(sza_deg, vza_deg, raa_deg, wl_nm,
                             pressure_hpa=1013.25, n_mu_gl=24, nt=40,
                             m_max=2, max_iterations=20, tolerance=1.0e-7,
                             absorption=None):
    """Return (rho_I, tau_R, n_orders_per_m).  Rayleigh-only, black surface.
    If absorption (an Absorption instance) is given, gas absorption is applied
    to match the C default (gas-on) behavior."""
    wl_nm = require_wavelength(wl_nm, context='Rayleigh atmosphere')
    mu_sun = cos(sza_deg * pi / 180.0)
    mu_view = cos(vza_deg * pi / 180.0)
    tau_R = tau_rayleigh_bodhaine(wl_nm, pressure_hpa, 45.0, 0.0, 360.0)

    atm, vj, z_km_level = build_atm_rayleigh(nt, tau_R, mu_sun, mu_view, n_mu_gl)
    if absorption is not None:
        apply_gas_absorption(atm, absorption, wl_nm, z_km_level)
    n_mu = atm.n_mu
    dirs = 2 * n_mu + 1
    l_max = max(m_max, 2)
    ws = kn.LegendreWorkspace(n_mu, l_max)

    zero = np.zeros((n_mu + 1, dirs))
    kt = S.KernelTables(zero, zero, zero, zero, zero, zero)

    I_per_m = np.zeros(m_max + 1)
    Q_per_m = np.zeros(m_max + 1)
    U_per_m = np.zeros(m_max + 1)
    n_orders = []

    for m in range(0, m_max + 1):
        xpl, xrl, xtl = kn.legendre_compute_pol(ws, atm.rm, m)
        atm.xpl = xpl
        atm.xrl = xrl
        atm.xtl = xtl

        src_i = S.primary_source(atm, m, ws, zero)
        src_q, src_u = S.primary_source_pol(atm, m, ws, zero, zero)
        prim_i = S.integrate_bcs(atm, src_i)
        prim_q = S.integrate_bcs(atm, src_q)
        prim_u = S.integrate_bcs(atm, src_u)

        tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_black(
            atm, m, kt, prim_i, prim_q, prim_u, max_iterations, tolerance)
        n_orders.append(n_ord)

        I_per_m[m] = tot_i[0, n_mu + vj]
        Q_per_m[m] = tot_q[0, n_mu + vj]
        U_per_m[m] = tot_u[0, n_mu + vj]

    dphi = raa_deg * pi / 180.0
    I_TOA = S.reconstruct_phi(I_per_m, m_max, dphi)
    rho_I = I_TOA / mu_sun if mu_sun > 0.0 else 0.0
    return rho_I, tau_R, n_orders


def solve_atm_black_fresnel_ocean(sza_deg, vza_deg, raa_deg, wl_nm, wind_ms,
                      pressure_hpa=1013.25, n_mu_gl=24, nt=40, m_max=2,
                      max_iterations=20, tolerance=1.0e-7, n_water=1.34,
                      sigma_type=1, q_convention=1, n_phi_quad=1024,
                      absorption=None):
    """ρ_R: Rayleigh atmosphere over a rough Fresnel interface with a black ocean below (no water-leaving),
    direct sunglint decoupled.  Matches C --surface black_fresnel_ocean --aod 0
    --decouple-sunglint.  I component returned."""
    wl_nm = require_wavelength(wl_nm, context='black Fresnel atmosphere')
    from . import surface as SF
    mu_sun = cos(sza_deg * pi / 180.0)
    mu_view = cos(vza_deg * pi / 180.0)
    tau_R = tau_rayleigh_bodhaine(wl_nm, pressure_hpa, 45.0, 0.0, 360.0)

    atm, vj, z_km_level = build_atm_rayleigh(nt, tau_R, mu_sun, mu_view, n_mu_gl)
    if absorption is not None:
        apply_gas_absorption(atm, absorption, wl_nm, z_km_level)
    n_mu = atm.n_mu
    dirs = 2 * n_mu + 1
    l_max = max(m_max, 2)
    ws = kn.LegendreWorkspace(n_mu, l_max)
    mu_pos = atm.rm[n_mu + 1:].copy()          # (n_mu,), includes view node

    zero = np.zeros((n_mu + 1, dirs))
    kt = S.KernelTables(zero, zero, zero, zero, zero, zero)

    I_per_m = np.zeros(m_max + 1)
    Q_per_m = np.zeros(m_max + 1)
    U_per_m = np.zeros(m_max + 1)
    n_orders = []

    for m in range(0, m_max + 1):
        xpl, xrl, xtl = kn.legendre_compute_pol(ws, atm.rm, m)
        atm.xpl = xpl
        atm.xrl = xrl
        atm.xtl = xtl

        # air-side rough-Fresnel m-mode surface kernel R_m (n_mu,n_mu,9)
        R_m = SF.fourier_kernel(SF.R_rough_fresnel_trig, mu_pos, mu_pos, m,
                                n_phi_quad, wind_ms, sigma_type, n_water,
                                q_convention, osoaa_sign_fix=True)
        mu_factor = (2.0 * pi) if m == 0 else pi   # rough-Fresnel

        src_i = S.primary_source(atm, m, ws, zero)
        src_q, src_u = S.primary_source_pol(atm, m, ws, zero, zero)
        prim_i = S.integrate_bcs(atm, src_i)
        prim_q = S.integrate_bcs(atm, src_q)
        prim_u = S.integrate_bcs(atm, src_u)

        # r_first: first-order direct-sun surface reflection (surf_seed).
        # R_solar[k_o] = kernel from solar direction to upward mu_k.
        R_solar = SF.fourier_kernel(SF.R_rough_fresnel_trig, mu_pos,
                                    np.array([mu_sun]), m, n_phi_quad, wind_ms,
                                    sigma_type, n_water, q_convention,
                                    osoaa_sign_fix=True)          # (n_mu,1,9)
        tau_total = atm.h[nt]
        trans_sun = exp(-tau_total / mu_sun)
        m_factor = pi if m == 0 else 0.5 * pi
        sI = np.zeros((nt + 1, dirs))
        sQ = np.zeros((nt + 1, dirs))
        sU = np.zeros((nt + 1, dirs))
        for i_lev in range(nt + 1):
            dtau_up = tau_total - atm.h[i_lev]
            trans = trans_sun * np.exp(-dtau_up / mu_pos)         # (n_mu,)
            factor = m_factor * mu_sun * trans                    # (n_mu,)
            sI[i_lev, n_mu + 1:] = R_solar[:, 0, 0] * factor
            sQ[i_lev, n_mu + 1:] = R_solar[:, 0, 3] * factor
            sU[i_lev, n_mu + 1:] = R_solar[:, 0, 6] * factor
        sI[:, n_mu] = np.nan
        sQ[:, n_mu] = np.nan
        sU[:, n_mu] = np.nan

        tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_surface(
            atm, m, kt, prim_i, prim_q, prim_u, R_m, mu_factor,
            max_iterations, tolerance, surf_seed=(sI, sQ, sU))
        n_orders.append(n_ord)

        I_per_m[m] = tot_i[0, n_mu + vj]
        Q_per_m[m] = tot_q[0, n_mu + vj]
        U_per_m[m] = tot_u[0, n_mu + vj]

    dphi = raa_deg * pi / 180.0
    I_TOA = S.reconstruct_phi(I_per_m, m_max, dphi)
    rho_I = I_TOA / mu_sun if mu_sun > 0.0 else 0.0
    return rho_I, tau_R, n_orders


# 6SV ODA550.f an23 (V=23 km) aerosol number density (rt_atm.c).
_Z_AER = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                   20,21,22,23,24,25,30,35,40,45,50,70,100,99999.0])
_AN23 = np.array([2.828e+03,1.244e+03,5.371e+02,2.256e+02,1.192e+02,
                  8.987e+01,6.337e+01,5.890e+01,6.069e+01,5.818e+01,
                  5.675e+01,5.317e+01,5.585e+01,5.156e+01,5.048e+01,
                  4.744e+01,4.511e+01,4.458e+01,4.314e+01,3.634e+01,
                  2.667e+01,1.933e+01,1.455e+01,1.113e+01,8.826e+00,
                  7.429e+00,2.238e+00,5.890e-01,1.550e-01,4.082e-02,
                  1.078e-02,5.550e-05,1.969e-08,0.0])
_AN23_TOTAL = float(np.sum(0.5 * (_AN23[:32] + _AN23[1:33]) * (_Z_AER[1:33] - _Z_AER[:32])))


def _tau_a_above(z, tau_a):
    """TAU_A_ABOVE: tau_a * (an23 mass above z / total), trapezoid."""
    if z >= _Z_AER[33]:
        return 0.0
    mass_above = 0.0
    k = 0
    while k < 32 and _Z_AER[k + 1] <= z:
        k += 1
    if k < 32 and _Z_AER[k] <= z < _Z_AER[k + 1]:
        n_at_z = _AN23[k] + (_AN23[k + 1] - _AN23[k]) * \
            (z - _Z_AER[k]) / (_Z_AER[k + 1] - _Z_AER[k])
        mass_above += 0.5 * (n_at_z + _AN23[k + 1]) * (_Z_AER[k + 1] - z)
        for kk in range(k + 1, 32):
            mass_above += 0.5 * (_AN23[kk] + _AN23[kk + 1]) * (_Z_AER[kk + 1] - _Z_AER[kk])
    return tau_a * mass_above / _AN23_TOTAL


def _tau_a_above_exp(z, tau_a, ha, z_top=100.0):
    """TAU_A_ABOVE_EXP: single-exponential exp(-z/ha) aerosol mass above z.
    Matches C when aer_h_km>0 (default 2.0, OSOAA -AP.HA).  Array-safe."""
    zz = xp.asarray(z, dtype=float)
    scalar = (zz.ndim == 0)
    total = ha * (1.0 - exp(-z_top / ha))
    mass_above = ha * (xp.exp(-zz / ha) - exp(-z_top / ha))
    res = tau_a * mass_above / total
    res = xp.where(zz >= z_top, 0.0, res)
    return float(res) if scalar else res


def build_atm_aerosol(nt, tau_R, tau_a, ssa_a, mu_sun, mu_view, n_mu_gl,
                      betal, gammal, alphal, zetal, L_max, aer_h_km=2.0):
    """rt_atm_build_aerosol_rayleigh: Rayleigh (tau_total=tau_R+tau_a) then
    overwrite xdel/ydel with the aerosol vertical profile.  aer_h_km>0 (default
    2.0, OSOAA -AP.HA) uses a single exponential exp(-z/ha); aer_h_km=0 uses the
    6SV an23 table.  Stores aerosol Greek coefficients on atm."""
    tau_total = tau_R + tau_a
    atm, vj, z_km_level = build_atm_rayleigh(nt, tau_total, mu_sun, mu_view, n_mu_gl)
    z_top = 100.0
    ha = aer_h_km

    def tau_a_above(z):
        if ha > 0.0:
            return _tau_a_above_exp(z, tau_a, ha, z_top)
        return _tau_a_above(z, tau_a)

    def tau_ray_above(z):
        # OCRT v1.2: Rayleigh follows the US62 cumulative molecular column.
        return tau_R * mp.us62_grid_fraction_above(z)

    if tau_R == 0.0 and tau_a > 0.0:
        atm.xdel[0] = ssa_a
        atm.ydel[0] = 0.0
    else:
        atm.xdel[0] = 0.0
        atm.ydel[0] = 1.0
    z_km_level[0] = z_top

    # 모든 층 독립 bisection (상한을 z_top 으로 통일; 64회면 근이 순차판과
    # ~1e-17 수준으로 일치, 최종 rho 무영향).  j==nt 은 z=0 명시.
    target = atm.h[1:nt + 1]
    z_lo = xp.zeros(nt)
    z_hi = xp.full(nt, z_top)
    for _ in range(64):
        z_mid = 0.5 * (z_lo + z_hi)
        tau_above = tau_ray_above(z_mid) + tau_a_above(z_mid)
        above = tau_above > target
        z_lo = xp.where(above, z_mid, z_lo)
        z_hi = xp.where(above, z_hi, z_mid)
    z_j = 0.5 * (z_lo + z_hi)
    z_j[nt - 1] = 0.0
    z_km_level[1:] = z_j
    ray_full = tau_ray_above(z_km_level)
    ea_full = tau_a_above(z_km_level)
    dt_ray = ray_full[1:] - ray_full[:-1]
    dt_aer = ea_full[1:] - ea_full[:-1]
    dt = dt_ray + dt_aer
    pos = dt > 0.0
    ydel_else = 1.0 if tau_R > 0.0 else 0.0
    atm.xdel[1:] = xp.where(pos, dt_aer * ssa_a / dt, 0.0)
    atm.ydel[1:] = xp.where(pos, dt_ray / dt, ydel_else)

    atm.ssa = (tau_R + tau_a * ssa_a) / tau_total
    atm.betal = xp.asarray(betal[:L_max + 1], dtype=float)
    atm.gammal = xp.asarray(gammal[:L_max + 1], dtype=float)
    atm.alphal = xp.asarray(alphal[:L_max + 1], dtype=float)
    atm.zetal = xp.asarray(zetal[:L_max + 1], dtype=float)
    atm.L_max = L_max
    return atm, vj, z_km_level




def build_atm_rayleigh_views(nt, tau_R, mu_sun, mu_views, n_mu_gl):
    """Multi-view analogue of :func:`build_atm_rayleigh`.

    Every requested view cosine is inserted as a zero-weight ordinate in one
    shared angular ring.  The atmosphere SOS field is therefore solved once
    and sampled exactly for all requested VZAs.
    """
    gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_gl)
    u = S.UAngles(list(gl_mu), list(gl_w))
    mu_views = [float(x) for x in mu_views]
    for mv in mu_views:
        u.add(mv)
    mu_ring = xp.asarray(u.mu)
    w_ring = xp.asarray(u.w)
    n_mu = len(mu_ring)

    atm = S.Atm(nt, n_mu)
    ron = 2.0 * (1.0 - DEPOL) / (2.0 + DEPOL)
    atm.mu_sun = mu_sun
    atm.beta0 = 1.0
    atm.beta2 = 0.5 * ron
    atm.gamma2 = -ron * sqrt(1.5)
    atm.alpha2 = 3.0 * ron
    atm.beam_q = 0.0

    k_arr = xp.arange(nt + 1, dtype=float)
    atm.h = k_arr * tau_R / nt
    atm.ch = 0.5 * xp.exp(-atm.h / mu_sun)
    atm.xdel = xp.zeros(nt + 1)
    atm.ydel = xp.ones(nt + 1)
    atm.rm[n_mu] = -mu_sun
    atm.gb[n_mu] = 0.0
    atm.rm[n_mu + 1:] = mu_ring
    atm.rm[:n_mu] = -xp.flip(mu_ring)
    atm.gb[n_mu + 1:] = w_ring
    atm.gb[:n_mu] = xp.flip(w_ring)

    z_cap = 100.0
    z_km_level = xp.zeros(nt + 1)
    z_km_level[1:] = mp.us62_altitude_from_grid_fraction(atm.h[1:] / tau_R)
    z_km_level[0] = z_cap
    vjs = []
    for mv in mu_views:
        vjs.append(int(xp.argmin(xp.abs(mu_ring - mv))) + 1)
    return atm, vjs, z_km_level


def build_atm_aerosol_views(nt, tau_R, tau_a, ssa_a, mu_sun, mu_views,
                             n_mu_gl, betal, gammal, alphal, zetal, L_max,
                             aer_h_km=2.0):
    """Multi-view aerosol/Rayleigh atmosphere on one shared angular ring."""
    tau_total = tau_R + tau_a
    atm, vjs, z_km_level = build_atm_rayleigh_views(
        nt, tau_total, mu_sun, mu_views, n_mu_gl)
    z_top = 100.0
    ha = aer_h_km

    def tau_a_above(z):
        if ha > 0.0:
            return _tau_a_above_exp(z, tau_a, ha, z_top)
        return _tau_a_above(z, tau_a)

    def tau_ray_above(z):
        return tau_R * mp.us62_grid_fraction_above(z)

    if tau_R == 0.0 and tau_a > 0.0:
        atm.xdel[0] = ssa_a
        atm.ydel[0] = 0.0
    else:
        atm.xdel[0] = 0.0
        atm.ydel[0] = 1.0
    z_km_level[0] = z_top

    target = atm.h[1:nt + 1]
    z_lo = xp.zeros(nt)
    z_hi = xp.full(nt, z_top)
    for _ in range(64):
        z_mid = 0.5 * (z_lo + z_hi)
        tau_above = tau_ray_above(z_mid) + tau_a_above(z_mid)
        above = tau_above > target
        z_lo = xp.where(above, z_mid, z_lo)
        z_hi = xp.where(above, z_hi, z_mid)
    z_j = 0.5 * (z_lo + z_hi)
    z_j[nt - 1] = 0.0
    z_km_level[1:] = z_j
    ray_full = tau_ray_above(z_km_level)
    ea_full = tau_a_above(z_km_level)
    dt_ray = ray_full[1:] - ray_full[:-1]
    dt_aer = ea_full[1:] - ea_full[:-1]
    dt = dt_ray + dt_aer
    pos = dt > 0.0
    ydel_else = 1.0 if tau_R > 0.0 else 0.0
    atm.xdel[1:] = xp.where(pos, dt_aer * ssa_a / dt, 0.0)
    atm.ydel[1:] = xp.where(pos, dt_ray / dt, ydel_else)
    atm.ssa = (tau_R + tau_a * ssa_a) / tau_total
    atm.betal = xp.asarray(betal[:L_max + 1], dtype=float)
    atm.gammal = xp.asarray(gammal[:L_max + 1], dtype=float)
    atm.alphal = xp.asarray(alphal[:L_max + 1], dtype=float)
    atm.zetal = xp.asarray(zetal[:L_max + 1], dtype=float)
    atm.L_max = L_max
    return atm, vjs, z_km_level

def solve_atm_aerosol_black_fresnel_ocean(sza_deg, vza_deg, raa_deg, wl_nm, wind_ms,
                              tau_a, ssa_a, betal, gammal, alphal, zetal,
                              L_max=80, pressure_hpa=1013.25, n_mu_gl=24,
                              nt=400, m_max=16, max_iterations=100,
                              tolerance=1.0e-7, n_water=1.34, sigma_type=1,
                              q_convention=1, n_phi_quad=1024, absorption=None):
    """rho_R+A: Rayleigh+aerosol atmosphere over a rough Fresnel interface with a black ocean below,
    direct sunglint decoupled.  Matches C --surface black_fresnel_ocean --aod>0 --mie ...
    --decouple-sunglint.  I component returned."""
    wl_nm = require_wavelength(wl_nm, context='aerosol atmosphere')
    from . import surface as SF
    mu_sun = cos(sza_deg * pi / 180.0)
    mu_view = cos(vza_deg * pi / 180.0)
    tau_R = tau_rayleigh_bodhaine(wl_nm, pressure_hpa, 45.0, 0.0, 360.0)

    atm, vj, z_km_level = build_atm_aerosol(nt, tau_R, tau_a, ssa_a, mu_sun,
                                            mu_view, n_mu_gl, betal, gammal,
                                            alphal, zetal, L_max)
    if absorption is not None:
        apply_gas_absorption(atm, absorption, wl_nm, z_km_level)
    n_mu = atm.n_mu
    dirs = 2 * n_mu + 1
    l_max = max(m_max, L_max)
    ws = kn.LegendreWorkspace(n_mu, l_max)
    mu_pos = atm.rm[n_mu + 1:].copy()

    I_per_m = np.zeros(m_max + 1)
    n_orders = []

    for m in range(0, m_max + 1):
        xpl, xrl, xtl = kn.legendre_compute_pol(ws, atm.rm, m)
        atm.xpl = xpl
        atm.xrl = xrl
        atm.xtl = xtl
        # aerosol phase Fourier kernels (betal/gammal/alphal/zetal) -> KernelTables
        kn.kernel_phase_fourier(ws, m, atm.betal)
        kn.kernel_phase_fourier_pol(ws, m, atm.gammal)
        kn.kernel_phase_fourier_aerosol_full(ws, m, atm.alphal, atm.zetal)
        kt = S.KernelTables(ws.pfm, ws.gr, ws.gt, ws.arr, ws.art, ws.att)

        R_m = SF.fourier_kernel(SF.R_rough_fresnel_trig, mu_pos, mu_pos, m,
                                n_phi_quad, wind_ms, sigma_type, n_water,
                                q_convention, osoaa_sign_fix=True)
        mu_factor = (2.0 * pi) if m == 0 else pi

        zero = np.zeros((n_mu + 1, dirs))
        src_i = S.primary_source(atm, m, ws, kt.pfm)
        src_q, src_u = S.primary_source_pol(atm, m, ws, kt.gr, kt.gt)
        prim_i = S.integrate_bcs(atm, src_i)
        prim_q = S.integrate_bcs(atm, src_q)
        prim_u = S.integrate_bcs(atm, src_u)

        # r_first (direct-sun surface reflection, surf_seed)
        R_solar = SF.fourier_kernel(SF.R_rough_fresnel_trig, mu_pos,
                                    np.array([mu_sun]), m, n_phi_quad, wind_ms,
                                    sigma_type, n_water, q_convention,
                                    osoaa_sign_fix=True)
        tau_total = atm.h[nt]
        trans_sun = exp(-tau_total / mu_sun)
        m_factor = pi if m == 0 else 0.5 * pi
        sI = np.zeros((nt + 1, dirs)); sQ = np.zeros((nt + 1, dirs)); sU = np.zeros((nt + 1, dirs))
        for i_lev in range(nt + 1):
            dtau_up = tau_total - atm.h[i_lev]
            trans = trans_sun * np.exp(-dtau_up / mu_pos)
            factor = m_factor * mu_sun * trans
            sI[i_lev, n_mu + 1:] = R_solar[:, 0, 0] * factor
            sQ[i_lev, n_mu + 1:] = R_solar[:, 0, 3] * factor
            sU[i_lev, n_mu + 1:] = R_solar[:, 0, 6] * factor
        sI[:, n_mu] = np.nan; sQ[:, n_mu] = np.nan; sU[:, n_mu] = np.nan

        tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_surface(
            atm, m, kt, prim_i, prim_q, prim_u, R_m, mu_factor,
            max_iterations, tolerance, surf_seed=(sI, sQ, sU))
        n_orders.append(n_ord)
        I_per_m[m] = tot_i[0, n_mu + vj]

    dphi = raa_deg * pi / 180.0
    I_TOA = S.reconstruct_phi(I_per_m, m_max, dphi)
    rho_I = I_TOA / mu_sun if mu_sun > 0.0 else 0.0
    return rho_I, tau_R, n_orders


def solve_aerosol_black_fresnel_ocean_value(sza_deg, vza_deg, raa_deg, wl_nm, wind_ms, mie,
                                user_aod, L_max=80, pressure_hpa=1013.25,
                                n_mu_gl=24, nt=400, m_max=16, max_iterations=100,
                                tolerance=1.0e-7, n_water=1.34, sigma_type=1,
                                q_convention=1, n_phi_quad=1024, nphi_value=720,
                                surface='black_fresnel_ocean', absorption=None, aer_h_km=2.0,
                                return_boa=False, use_value_kernel=True,
                                bottom_source=None, full_result=False,
                                aerosol_runtime=None):
    """rho_R+A with the C production aerosol protocol: loglin forward-peak
    truncation + value-kernel scalar I phase (pfm) + Greek-coefficient Q/U
    (gr/gt/arr/att).  Matches C aerosol-on auto-dispatch.  I returned."""
    wl_nm = require_wavelength(wl_nm, context='aerosol value atmosphere')
    from . import aerosol as AER
    from . import surface as SF
    mu_sun = cos(sza_deg * pi / 180.0)
    mu_view = cos(vza_deg * pi / 180.0)
    # One immutable aerosol object may be prepared before a full-grid loop and
    # passed to every cell.  Backward-compatible callers that provide raw Mie
    # data still use the same centralized preparation path, with user_aod
    # interpreted as AOD at the current calculation wavelength.
    if aerosol_runtime is None:
        if mie is None:
            if user_aod > 0.0:
                raise ValueError(
                    "Aerosol runtime is required when AOD > 0 in ocean full-grid mode.")
            raise ValueError("MieData or AerosolRuntime is required")
        aerosol_runtime = AER.prepare_aerosol_runtime(
            mie, wl_nm, user_aod, aod_ref_nm=wl_nm, L_max=L_max)
    else:
        if abs(float(aerosol_runtime.wavelength_nm) - float(wl_nm)) > 1.0e-10:
            raise ValueError("AerosolRuntime wavelength does not match solver wavelength")
        if int(aerosol_runtime.L_max) != int(L_max):
            raise ValueError("AerosolRuntime L_max does not match solver L_max")
        if abs(float(aerosol_runtime.aod_target) - float(user_aod)) > \
                1.0e-12 * max(1.0, abs(float(user_aod))):
            raise ValueError("AerosolRuntime AOD does not match solver AOD")

    A = aerosol_runtime.truncation_A
    tau_a_eff = aerosol_runtime.tau_a_eff
    ssa_a_eff = aerosol_runtime.ssa_a_eff
    ts = aerosol_runtime.theta_asc
    P_norm_asc = aerosol_runtime.P11_norm_asc
    be = aerosol_runtime.betal
    ga = aerosol_runtime.gammal
    al = aerosol_runtime.alphal
    ze = aerosol_runtime.zetal

    tau_R = tau_rayleigh_bodhaine(wl_nm, pressure_hpa, 45.0, 0.0, 360.0)
    atm, vj, z_km_level = build_atm_aerosol(nt, tau_R, tau_a_eff, ssa_a_eff, mu_sun,
                                            mu_view, n_mu_gl, be, ga, al, ze, L_max, aer_h_km=aer_h_km)
    tau_scat_total = atm.h[nt]           # tau_R + tau_a_eff (before gas OD)
    if absorption is not None:
        apply_gas_absorption(atm, absorption, wl_nm, z_km_level)
    n_mu = atm.n_mu
    dirs = 2 * n_mu + 1
    l_max = max(m_max, L_max)
    ws = kn.LegendreWorkspace(n_mu, l_max)
    mu_pos = atm.rm[n_mu + 1:].copy()

    if surface == 'coxmunk':
        import warnings
        warnings.warn("surface='coxmunk' is deprecated; use "
                      "surface='black_fresnel_ocean'. The mode denotes a rough "
                      "Fresnel interface over a black ocean; slope variance is "
                      "selected separately.", DeprecationWarning, stacklevel=2)
        surface = 'black_fresnel_ocean'
    if surface not in ('black_fresnel_ocean', 'black'):
        raise ValueError("surface must be 'black_fresnel_ocean' or 'black'")
    use_surface = (surface == 'black_fresnel_ocean')
    I_per_m = np.zeros(m_max + 1)
    Q_per_m = np.zeros(m_max + 1)
    U_per_m = np.zeros(m_max + 1)
    # BOA downward skylight per m per positive node (for ocean coupling, Phase D)
    boa_I = np.zeros((m_max + 1, n_mu))
    boa_Q = np.zeros((m_max + 1, n_mu))
    boa_U = np.zeros((m_max + 1, n_mu))

    for m in range(0, m_max + 1):
        xpl, xrl, xtl = kn.legendre_compute_pol(ws, atm.rm, m)
        atm.xpl = xpl; atm.xrl = xrl; atm.xtl = xtl
        # value-kernel scalar I phase (overwrites pfm); Greek Q/U kernels
        if use_value_kernel:
            pfm_val = AER.aerosol_value_pfm(atm, m, ts, P_norm_asc, nphi_value)
            ws.pfm[:, :] = pfm_val
        else:
            kn.kernel_phase_fourier(ws, m, atm.betal)
        kn.kernel_phase_fourier_pol(ws, m, atm.gammal)
        kn.kernel_phase_fourier_aerosol_full(ws, m, atm.alphal, atm.zetal)
        kt = S.KernelTables(ws.pfm, ws.gr, ws.gt, ws.arr, ws.art, ws.att)

        if bottom_source is None:
            src_i = S.primary_source(atm, m, ws, kt.pfm)
            src_q, src_u = S.primary_source_pol(atm, m, ws, kt.gr, kt.gt)
            prim_i = S.integrate_bcs(atm, src_i)
            prim_q = S.integrate_bcs(atm, src_q)
            prim_u = S.integrate_bcs(atm, src_u)

        if use_surface:
            R_m = SF.fourier_kernel(SF.R_rough_fresnel_trig, mu_pos, mu_pos, m,
                                    n_phi_quad, wind_ms, sigma_type, n_water,
                                    q_convention, osoaa_sign_fix=True)
            mu_factor = (2.0 * pi) if m == 0 else pi
            if bottom_source is not None:
                # C3 pass-2 (bottom_source_only): solar primary zeroed; the
                # injected water-leaving 0+ field enters as an up-attenuated
                # initial field at every level (add_external_bottom_source).
                bI, bQ, bU = bottom_source
                tau_total = atm.h[nt]
                z0 = np.zeros((nt + 1, dirs))
                prim_i = z0.copy(); prim_q = z0.copy(); prim_u = z0.copy()
                prim_i[:, n_mu] = np.nan
                prim_q[:, n_mu] = np.nan
                prim_u[:, n_mu] = np.nan
                eI = np.zeros((nt + 1, dirs))
                eQ = np.zeros((nt + 1, dirs))
                eU = np.zeros((nt + 1, dirs))
                for i_lev in range(nt + 1):
                    atten = np.exp(-(tau_total - atm.h[i_lev]) / mu_pos)
                    eI[i_lev, n_mu + 1:] = bI[m] * atten
                    eQ[i_lev, n_mu + 1:] = bQ[m] * atten
                    eU[i_lev, n_mu + 1:] = bU[m] * atten
                tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_surface(
                    atm, m, kt, prim_i, prim_q, prim_u, R_m, mu_factor,
                    max_iterations, tolerance, surf_seed=None,
                    ext_init=(eI, eQ, eU))
            else:
                R_solar = SF.fourier_kernel(SF.R_rough_fresnel_trig, mu_pos,
                                            np.array([mu_sun]), m, n_phi_quad, wind_ms,
                                            sigma_type, n_water, q_convention,
                                            osoaa_sign_fix=True)
                tau_total = atm.h[nt]
                trans_sun = exp(-tau_total / mu_sun)
                m_factor = pi if m == 0 else 0.5 * pi
                sI = np.zeros((nt + 1, dirs)); sQ = np.zeros((nt + 1, dirs)); sU = np.zeros((nt + 1, dirs))
                for i_lev in range(nt + 1):
                    trans = trans_sun * np.exp(-(tau_total - atm.h[i_lev]) / mu_pos)
                    factor = m_factor * mu_sun * trans
                    sI[i_lev, n_mu + 1:] = R_solar[:, 0, 0] * factor
                    sQ[i_lev, n_mu + 1:] = R_solar[:, 0, 3] * factor
                    sU[i_lev, n_mu + 1:] = R_solar[:, 0, 6] * factor
                sI[:, n_mu] = np.nan; sQ[:, n_mu] = np.nan; sU[:, n_mu] = np.nan
                tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_surface(
                    atm, m, kt, prim_i, prim_q, prim_u, R_m, mu_factor,
                    max_iterations, tolerance, surf_seed=(sI, sQ, sU))
        else:
            tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_black(
                atm, m, kt, prim_i, prim_q, prim_u, max_iterations, tolerance)
        I_per_m[m] = tot_i[0, n_mu + vj]
        Q_per_m[m] = tot_q[0, n_mu + vj]
        U_per_m[m] = tot_u[0, n_mu + vj]
        # BOA (k=nt) downward: idx_dn = n_mu - jp, jp=1..n_mu (C rt_solver.c:815)
        boa_I[m, :] = tot_i[nt, n_mu - 1::-1]
        boa_Q[m, :] = tot_q[nt, n_mu - 1::-1]
        boa_U[m, :] = tot_u[nt, n_mu - 1::-1]

    dphi = raa_deg * pi / 180.0
    I_TOA = S.reconstruct_phi(I_per_m, m_max, dphi)
    Q_TOA = S.reconstruct_phi(Q_per_m, m_max, dphi)
    U_TOA = S.reconstruct_phi_sin(U_per_m, m_max, dphi)
    rho_I = I_TOA / mu_sun if mu_sun > 0.0 else 0.0
    rho_Q = Q_TOA / mu_sun if mu_sun > 0.0 else 0.0
    rho_U = U_TOA / mu_sun if mu_sun > 0.0 else 0.0
    # T_diff_dn_hemi = (2pi sum mu w I_dn^{m=0}) / (pi mu_sun)  (view node w=0)
    w_pos = atm.gb[n_mu + 1:]
    F_dn_diff = float((mu_pos * w_pos * boa_I[0]).sum()) * 2.0 * pi
    T_diff_dn_hemi = F_dn_diff / (pi * mu_sun) if mu_sun > 0.0 else 0.0
    tau_abs_total = getattr(atm, 'tau_abs_total', 0.0)
    boa = {'I_per_m': boa_I, 'Q_per_m': boa_Q, 'U_per_m': boa_U,
           'mu_quad_pos': mu_pos.copy(), 'w_quad_pos': w_pos.copy(),
           'tau_atm_total': tau_scat_total + tau_abs_total,
           'mu_sun_air': mu_sun}
    if full_result:
        return {'rho_I': rho_I, 'rho_Q': rho_Q, 'rho_U': rho_U,
                'I_TOA': I_TOA, 'Q_TOA': Q_TOA, 'U_TOA': U_TOA,
                'T_diff_dn_hemi': T_diff_dn_hemi,
                'tau_R': tau_R, 'tau_a_eff': tau_a_eff, 'A': A,
                'tau_scat_total': tau_scat_total,
                'tau_abs_total': tau_abs_total, 'boa': boa}
    if return_boa:
        return rho_I, rho_Q, rho_U, tau_R, A, boa
    return rho_I, rho_Q, rho_U, tau_R, A


def solve_aerosol_black_fresnel_ocean_value_lut(
        sza_deg, vza_values, raa_values, wl_nm, wind_ms, mie, user_aod,
        L_max=80, pressure_hpa=1013.25, n_mu_gl=24, nt=400, m_max=16,
        max_iterations=100, tolerance=1.0e-7, n_water=1.34, sigma_type=1,
        q_convention=1, n_phi_quad=1024, nphi_value=720,
        surface='black_fresnel_ocean', absorption=None, aer_h_km=2.0,
        use_value_kernel=True, bottom_source=None, aerosol_runtime=None,
        surface_kernel_cache=None):
    """Native all-VZA/all-RAA atmosphere LUT for one case and wavelength.

    The SOS field for each Fourier order is solved once on a shared ring that
    contains every requested view as a zero-weight node.  RAA dependence is
    reconstructed from the stored Fourier samples only.
    """
    wl_nm = require_wavelength(wl_nm, context='aerosol atmosphere LUT')
    from . import aerosol as AER
    from . import surface as SF
    vza_values = [float(v) for v in vza_values]
    raa_values = [float(r) for r in raa_values]
    if not vza_values or not raa_values:
        raise ValueError('vza_values and raa_values must be non-empty')
    mu_sun = cos(sza_deg * pi / 180.0)
    mu_views = [cos(v * pi / 180.0) for v in vza_values]

    if aerosol_runtime is None:
        if mie is None:
            if user_aod > 0.0:
                raise ValueError('Aerosol runtime is required when AOD > 0')
            raise ValueError('MieData or AerosolRuntime is required')
        aerosol_runtime = AER.prepare_aerosol_runtime(
            mie, wl_nm, user_aod, aod_ref_nm=wl_nm, L_max=L_max)
    else:
        if abs(float(aerosol_runtime.wavelength_nm) - float(wl_nm)) > 1.0e-10:
            raise ValueError('AerosolRuntime wavelength does not match solver wavelength')
        if int(aerosol_runtime.L_max) != int(L_max):
            raise ValueError('AerosolRuntime L_max does not match solver L_max')
        if abs(float(aerosol_runtime.aod_target) - float(user_aod)) > \
                1.0e-12 * max(1.0, abs(float(user_aod))):
            raise ValueError('AerosolRuntime AOD does not match solver AOD')

    A = aerosol_runtime.truncation_A
    tau_a_eff = aerosol_runtime.tau_a_eff
    ssa_a_eff = aerosol_runtime.ssa_a_eff
    ts = aerosol_runtime.theta_asc
    P_norm_asc = aerosol_runtime.P11_norm_asc
    be = aerosol_runtime.betal
    ga = aerosol_runtime.gammal
    al = aerosol_runtime.alphal
    ze = aerosol_runtime.zetal

    tau_R = tau_rayleigh_bodhaine(wl_nm, pressure_hpa, 45.0, 0.0, 360.0)
    atm, vjs, z_km_level = build_atm_aerosol_views(
        nt, tau_R, tau_a_eff, ssa_a_eff, mu_sun, mu_views, n_mu_gl,
        be, ga, al, ze, L_max, aer_h_km=aer_h_km)
    tau_scat_total = atm.h[nt]
    if absorption is not None:
        apply_gas_absorption(atm, absorption, wl_nm, z_km_level)
    n_mu = atm.n_mu
    dirs = 2 * n_mu + 1
    l_max = max(m_max, L_max)
    ws = kn.LegendreWorkspace(n_mu, l_max)
    mu_pos = atm.rm[n_mu + 1:].copy()

    if surface == 'coxmunk':
        surface = 'black_fresnel_ocean'
    if surface not in ('black_fresnel_ocean', 'black'):
        raise ValueError("surface must be 'black_fresnel_ocean' or 'black'")
    use_surface = (surface == 'black_fresnel_ocean')

    def cached_surface_kernel(mu_o, mu_i, mode):
        """Exact per-case cache shared by atmosphere pass 1 and pass 2.

        The caller supplies a fresh dictionary for one case/wavelength.  A
        full byte key is nevertheless retained so accidental cross-case reuse
        cannot return a false hit.  Stored arrays are private copies; hits also
        return a copy because the SOS boundary routine is not part of the cache
        ownership contract.
        """
        mu_o_arr = np.ascontiguousarray(mu_o, dtype=np.float64)
        mu_i_arr = np.ascontiguousarray(mu_i, dtype=np.float64)
        key = (int(mode), int(n_phi_quad), float(wind_ms), int(sigma_type),
               float(n_water), int(q_convention),
               mu_o_arr.shape, mu_o_arr.tobytes(),
               mu_i_arr.shape, mu_i_arr.tobytes())
        if surface_kernel_cache is not None:
            cached = surface_kernel_cache.get(key)
            if cached is not None:
                return cached.copy()
        value = SF.fourier_kernel(
            SF.R_rough_fresnel_trig, mu_o_arr, mu_i_arr, mode,
            n_phi_quad, wind_ms, sigma_type, n_water, q_convention,
            osoaa_sign_fix=True)
        if surface_kernel_cache is not None:
            surface_kernel_cache[key] = value.copy()
        return value

    nv = len(vza_values)
    nr = len(raa_values)
    I_pm = np.zeros((m_max + 1, nv))
    Q_pm = np.zeros((m_max + 1, nv))
    U_pm = np.zeros((m_max + 1, nv))
    boa_I = np.zeros((m_max + 1, n_mu))
    boa_Q = np.zeros((m_max + 1, n_mu))
    boa_U = np.zeros((m_max + 1, n_mu))
    n_orders = 1

    for m in range(0, m_max + 1):
        xpl, xrl, xtl = kn.legendre_compute_pol(ws, atm.rm, m)
        atm.xpl = xpl; atm.xrl = xrl; atm.xtl = xtl
        if use_value_kernel:
            pfm_val = AER.aerosol_value_pfm(atm, m, ts, P_norm_asc, nphi_value)
            ws.pfm[:, :] = pfm_val
        else:
            kn.kernel_phase_fourier(ws, m, atm.betal)
        kn.kernel_phase_fourier_pol(ws, m, atm.gammal)
        kn.kernel_phase_fourier_aerosol_full(ws, m, atm.alphal, atm.zetal)
        kt = S.KernelTables(ws.pfm, ws.gr, ws.gt, ws.arr, ws.art, ws.att)

        if bottom_source is None:
            src_i = S.primary_source(atm, m, ws, kt.pfm)
            src_q, src_u = S.primary_source_pol(atm, m, ws, kt.gr, kt.gt)
            prim_i = S.integrate_bcs(atm, src_i)
            prim_q = S.integrate_bcs(atm, src_q)
            prim_u = S.integrate_bcs(atm, src_u)

        if use_surface:
            R_m = cached_surface_kernel(mu_pos, mu_pos, m)
            mu_factor = (2.0 * pi) if m == 0 else pi
            if bottom_source is not None:
                bI, bQ, bU = bottom_source
                tau_total = atm.h[nt]
                z0 = np.zeros((nt + 1, dirs))
                prim_i = z0.copy(); prim_q = z0.copy(); prim_u = z0.copy()
                prim_i[:, n_mu] = np.nan; prim_q[:, n_mu] = np.nan; prim_u[:, n_mu] = np.nan
                eI = np.zeros((nt + 1, dirs)); eQ = np.zeros((nt + 1, dirs)); eU = np.zeros((nt + 1, dirs))
                for i_lev in range(nt + 1):
                    atten = np.exp(-(tau_total - atm.h[i_lev]) / mu_pos)
                    eI[i_lev, n_mu + 1:] = bI[m] * atten
                    eQ[i_lev, n_mu + 1:] = bQ[m] * atten
                    eU[i_lev, n_mu + 1:] = bU[m] * atten
                tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_surface(
                    atm, m, kt, prim_i, prim_q, prim_u, R_m, mu_factor,
                    max_iterations, tolerance, surf_seed=None,
                    ext_init=(eI, eQ, eU))
            else:
                R_solar = cached_surface_kernel(
                    mu_pos, np.array([mu_sun]), m)
                tau_total = atm.h[nt]
                trans_sun = exp(-tau_total / mu_sun)
                m_factor = pi if m == 0 else 0.5 * pi
                sI = np.zeros((nt + 1, dirs)); sQ = np.zeros((nt + 1, dirs)); sU = np.zeros((nt + 1, dirs))
                for i_lev in range(nt + 1):
                    trans = trans_sun * np.exp(-(tau_total - atm.h[i_lev]) / mu_pos)
                    factor = m_factor * mu_sun * trans
                    sI[i_lev, n_mu + 1:] = R_solar[:, 0, 0] * factor
                    sQ[i_lev, n_mu + 1:] = R_solar[:, 0, 3] * factor
                    sU[i_lev, n_mu + 1:] = R_solar[:, 0, 6] * factor
                sI[:, n_mu] = np.nan; sQ[:, n_mu] = np.nan; sU[:, n_mu] = np.nan
                tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_surface(
                    atm, m, kt, prim_i, prim_q, prim_u, R_m, mu_factor,
                    max_iterations, tolerance, surf_seed=(sI, sQ, sU))
        else:
            tot_i, tot_q, tot_u, n_ord, conv, resid = S.sos_atm_black(
                atm, m, kt, prim_i, prim_q, prim_u, max_iterations, tolerance)
        n_orders = max(n_orders, int(n_ord))
        for iv, vj in enumerate(vjs):
            I_pm[m, iv] = tot_i[0, n_mu + vj]
            Q_pm[m, iv] = tot_q[0, n_mu + vj]
            U_pm[m, iv] = tot_u[0, n_mu + vj]
        boa_I[m, :] = tot_i[nt, n_mu - 1::-1]
        boa_Q[m, :] = tot_q[nt, n_mu - 1::-1]
        boa_U[m, :] = tot_u[nt, n_mu - 1::-1]

    I_grid = np.zeros((nv, nr)); Q_grid = np.zeros((nv, nr)); U_grid = np.zeros((nv, nr))
    for iv in range(nv):
        for ir, raa in enumerate(raa_values):
            dphi = raa * pi / 180.0
            I_grid[iv, ir] = S.reconstruct_phi(I_pm[:, iv], m_max, dphi)
            Q_grid[iv, ir] = S.reconstruct_phi(Q_pm[:, iv], m_max, dphi)
            U_grid[iv, ir] = S.reconstruct_phi_sin(U_pm[:, iv], m_max, dphi)
    w_pos = atm.gb[n_mu + 1:]
    F_dn_diff = float((mu_pos * w_pos * boa_I[0]).sum()) * 2.0 * pi
    T_diff_dn_hemi = F_dn_diff / (pi * mu_sun) if mu_sun > 0.0 else 0.0
    tau_abs_total = getattr(atm, 'tau_abs_total', 0.0)
    boa = {'I_per_m': boa_I, 'Q_per_m': boa_Q, 'U_per_m': boa_U,
           'mu_quad_pos': mu_pos.copy(), 'w_quad_pos': w_pos.copy(),
           'tau_atm_total': tau_scat_total + tau_abs_total,
           'mu_sun_air': mu_sun}
    return {
        'I_TOA': I_grid, 'Q_TOA': Q_grid, 'U_TOA': U_grid,
        'rho_I': I_grid / mu_sun, 'rho_Q': Q_grid / mu_sun,
        'rho_U': U_grid / mu_sun, 'I_per_m': I_pm, 'Q_per_m': Q_pm,
        'U_per_m': U_pm, 'T_diff_dn_hemi': T_diff_dn_hemi,
        'tau_R': tau_R, 'tau_a_eff': tau_a_eff, 'A': A,
        'tau_scat_total': tau_scat_total, 'tau_abs_total': tau_abs_total,
        'boa': boa, 'n_orders': n_orders, 'vza_values': vza_values,
        'raa_values': raa_values,
    }


# Deprecated high-level API aliases.  These preserve old scripts while all new
# code uses the physical boundary-condition name rather than a slope-model name.
def solve_atm_coxmunk(*args, **kwargs):
    import warnings
    warnings.warn("solve_atm_coxmunk is deprecated; use "
                  "solve_atm_black_fresnel_ocean", DeprecationWarning, stacklevel=2)
    return solve_atm_black_fresnel_ocean(*args, **kwargs)

def solve_atm_aerosol_coxmunk(*args, **kwargs):
    import warnings
    warnings.warn("solve_atm_aerosol_coxmunk is deprecated; use "
                  "solve_atm_aerosol_black_fresnel_ocean", DeprecationWarning, stacklevel=2)
    return solve_atm_aerosol_black_fresnel_ocean(*args, **kwargs)

def solve_aerosol_coxmunk_value(*args, **kwargs):
    import warnings
    warnings.warn("solve_aerosol_coxmunk_value is deprecated; use "
                  "solve_aerosol_black_fresnel_ocean_value",
                  DeprecationWarning, stacklevel=2)
    return solve_aerosol_black_fresnel_ocean_value(*args, **kwargs)
