"""rt_water_rt_sos_pure transliteration — fixed-bulk lut-deltam value-kernel
path (the Tier-0 production route).  Not a full port of every diagnostic
branch; env-gated diagnostics and unused modes are omitted, production
algebra is verbatim."""

# CANONICAL FR631 PARTICLE POLICY
# Raw FR631 with broad truncation OFF is the scientific reference. A local
# 0--0.005 degree cap was negligible in validation, whereas broad truncation
# can materially change b/omega/tau and directional Rrs/rrs. Keep detailed
# evidence in validation CSV/figure campaigns and synchronize this comment
# with docs/OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21.md.

import os
import numpy as np
from math import pi, sqrt, exp, cos, sin, log, acos

from . import kernel as kn
from . import surface as sf
from . import phase as ph
from . import sos
from . import value_phase as vp
from .rww import rww_coxmunk_allm, R_ww_coxmunk_direct

M_PI = pi

_VALUE_POL_NADIR_WARNED = False

def _warn_value_kernel_exact_nadir(view_vza_values):
    """Warn without changing geometry or solver data.

    The direct polarized value kernel has a known exact-nadir output-extraction
    singularity.  The validated operational workaround is VZA=0.001 deg
    (safe range 0.001--0.01 deg); angles below roughly 1e-6 deg can round back
    to the singular branch in float64.
    """
    global _VALUE_POL_NADIR_WARNED
    if _VALUE_POL_NADIR_WARNED:
        return
    if any(abs(float(v)) <= 1.0e-15 for v in view_vza_values):
        import warnings
        warnings.warn(
            'OCRT_WATER_VALUE_KERNEL_POL at exact VZA=0 deg has a known '
            'polarized output-extraction singularity; use VZA=0.001 deg '
            '(validated safe range 0.001--0.01 deg; avoid <1e-6 deg). '
            'The solver input is not modified automatically.',
            RuntimeWarning, stacklevel=3)
        _VALUE_POL_NADIR_WARNED = True


def snell_down(mu_air, n_water):
    sin_a_sq = max(0.0, 1.0 - mu_air * mu_air)
    sin_w_sq = sin_a_sq / (n_water * n_water)
    if sin_w_sq >= 1.0:
        return 0.0
    return sqrt(1.0 - sin_w_sq)


def raa_to_water_view_phi(raa_deg):
    return (180.0 - raa_deg) * M_PI / 180.0


def _linear_interp_ascending(mu_grid, vals, mu_t):
    n = len(mu_grid)
    if n < 2:
        return vals[0] if n == 1 else 0.0
    hi = 1
    while hi < n and mu_grid[hi] < mu_t:
        hi += 1
    lo = hi - 1
    if hi >= n:
        hi = n - 1
        lo = n - 2
    if lo < 0:
        lo, hi = 0, 1
    t = (mu_t - mu_grid[lo]) / (mu_grid[hi] - mu_grid[lo])
    return (1.0 - t) * vals[lo] + t * vals[hi]


def aw_interp_on_unsorted(mu, vals, mu_t):
    # C verbatim: defensive clamp to the first 64 entries BEFORE sorting.
    # For rings > 64 nodes (LUT-mode n_mu_water=96) this truncates the
    # interpolation grid to the 64 lowest-index ring entries, so high-mu
    # targets are linearly EXTRAPOLATED from the two highest kept nodes.
    # Faithful to rt_air_water_coupling.c:aw_interp_on_unsorted (n>64 -> 64).
    if len(mu) > 64:
        mu = mu[:64]
        vals = vals[:64]
    order = sorted(range(len(mu)), key=lambda i: mu[i])
    ms, tv = [], []
    for i in order:
        if ms and abs(mu[i] - ms[-1]) < 1e-12:
            continue
        ms.append(mu[i])
        tv.append(vals[i])
    return _linear_interp_ascending(ms, tv, mu_t)


class WaterResult:
    pass


def _add_diffuse_top_primary(atm, ws, m, ext_I, ext_Q, ext_U, ext_mu,
                             src_i, src_q, src_u):
    """ocrt_add_diffuse_top_primary (DT-HIOM TABLE kernel default).

    Adds the first-order scattering source of the externally injected
    downwelling diffuse field (per mode m, sampled on ext_mu) to the SOS
    source arrays.  Each positive ring column c acts as a downward beam
    exp(-h/mu_c); the scattering kernel rows are read from the solver's own
    per-m tables (value pfm + Greek gr/gt/arr/art/att), which is the C
    default (OCRT_EXTTOP_KERNEL unset).  The water-Rayleigh part rides
    xdel/ydel per layer (inert for fixed-bulk where ydel=0).
    """
    nt = atm.n_layers
    n_mu = atm.n_mu
    ray_on = 1 if m <= 2 else 0
    beta0_m = atm.beta0 if m == 0 else 0.0
    beta2 = atm.beta2
    gamma2 = atm.gamma2
    alpha2 = atm.alpha2
    ext_n = len(ext_mu)
    pfm, gr, gt = ws.pfm, ws.gr, ws.gt
    arr, art, att = ws.arr, ws.art, ws.att

    for c in range(1, n_mu + 1):
        mu_c = atm.rm[n_mu + c]
        # map ring column -> incident sample (exact match, else linear interp)
        AI = AQ = AU = 0.0
        hit = False
        for e in range(ext_n):
            if abs(ext_mu[e] - mu_c) < 1e-12:
                AI = 2.0 * atm.gb[n_mu + c] * ext_I[e]
                AQ = 2.0 * atm.gb[n_mu + c] * ext_Q[e] if ext_Q is not None else 0.0
                AU = 2.0 * atm.gb[n_mu + c] * ext_U[e] if ext_U is not None else 0.0
                hit = True
                break
        if not hit and AI == 0.0 and AQ == 0.0 and AU == 0.0:
            # DTP-FIX: inserted (view/slot) columns take clamped linear interp
            if ext_n >= 2:
                lo = 0
                while lo + 1 < ext_n - 1 and ext_mu[lo + 1] <= mu_c:
                    lo += 1
                m0_, m1_ = ext_mu[lo], ext_mu[lo + 1]
                t_ = (mu_c - m0_) / (m1_ - m0_) if m1_ != m0_ else 0.0
                if t_ < 0.0:
                    t_ = 0.0
                elif t_ > 1.0:
                    t_ = 1.0
                AI = 2.0 * atm.gb[n_mu + c] * ((1.0 - t_) * ext_I[lo] + t_ * ext_I[lo + 1])
                if ext_Q is not None:
                    AQ = 2.0 * atm.gb[n_mu + c] * ((1.0 - t_) * ext_Q[lo] + t_ * ext_Q[lo + 1])
                if ext_U is not None:
                    AU = 2.0 * atm.gb[n_mu + c] * ((1.0 - t_) * ext_U[lo] + t_ * ext_U[lo + 1])
            if AI == 0.0 and AQ == 0.0 and AU == 0.0:
                continue
        cc = -c                            # downward incident column
        # kernel rows via the solver's per-m tables (TABLE default)
        js = np.arange(-n_mu, n_mu + 1)
        rowI = np.zeros(2 * n_mu + 1); rowQ = np.zeros(2 * n_mu + 1); rowU = np.zeros(2 * n_mu + 1)
        rowIq = np.zeros(2 * n_mu + 1); rowQq = np.zeros(2 * n_mu + 1); rowUq = np.zeros(2 * n_mu + 1)
        rowIu = np.zeros(2 * n_mu + 1); rowQu = np.zeros(2 * n_mu + 1); rowUu = np.zeros(2 * n_mu + 1)
        for j in js:
            if j > 0:
                sI = pfm[c, n_mu - j]
                sQ = gr[c, n_mu - j]
                sU = -gt[c, n_mu - j]
                qI = gr[j, n_mu - c]
                qQ = arr[c, n_mu - j]
                qU = +art[j, n_mu - c]
                uI = -gt[j, n_mu - c]
                uQ = +art[c, n_mu - j]
                uU = -att[c, n_mu - j]
            elif j < 0:
                a = -j
                sI = pfm[c, n_mu + a]
                sQ = gr[c, n_mu + a]
                sU = -gt[c, n_mu + a]
                qI = gr[a, n_mu + c]
                qQ = arr[c, n_mu + a]
                qU = -art[a, n_mu + c]
                uI = +gt[a, n_mu + c]
                uQ = +art[c, n_mu + a]
                uU = -att[c, n_mu + a]
            else:
                sI = sQ = sU = qI = qQ = qU = uI = uQ = uU = 0.0
            idx = j + n_mu
            rowI[idx] = sI; rowQ[idx] = sQ; rowU[idx] = sU
            rowIq[idx] = qI; rowQq[idx] = qQ; rowUq[idx] = qU
            rowIu[idx] = uI; rowQu[idx] = uQ; rowUu[idx] = uU
        xplc = atm.xpl[cc + n_mu]
        xrlc = atm.xrl[cc + n_mu]
        xtlc = atm.xtl[cc + n_mu]
        xpl_j = atm.xpl
        xrl_j = atm.xrl
        xtl_j = atm.xtl
        raI = beta0_m + beta2 * xpl_j * xplc
        raQ = (gamma2 * xrl_j * xplc) if ray_on else np.zeros_like(xpl_j)
        raU = (gamma2 * xtl_j * xplc) if ray_on else np.zeros_like(xpl_j)
        rqI = (gamma2 * xpl_j * xrlc) if ray_on else np.zeros_like(xpl_j)
        rqQ = (alpha2 * xrl_j * xrlc) if ray_on else np.zeros_like(xpl_j)
        rqU = (alpha2 * xtl_j * xrlc) if ray_on else np.zeros_like(xpl_j)
        ruI = (gamma2 * xpl_j * xtlc) if ray_on else np.zeros_like(xpl_j)
        ruQ = (alpha2 * xrl_j * xtlc) if ray_on else np.zeros_like(xpl_j)
        ruU = (alpha2 * xtl_j * xtlc) if ray_on else np.zeros_like(xpl_j)
        mask = np.ones(2 * n_mu + 1, dtype=bool)
        mask[n_mu] = False               # solar slot stays poisoned
        for k in range(nt + 1):
            ch_c = 0.5 * exp(-atm.h[k] / mu_c)
            xd = atm.xdel[k]
            yd = atm.ydel[k]
            src_i[k, mask] += ch_c * (AI * (xd * rowI[mask] + yd * raI[mask])
                                      + AQ * (xd * rowIq[mask] + yd * rqI[mask])
                                      + AU * (xd * rowIu[mask] + yd * ruI[mask]))
            src_q[k, mask] += +ch_c * (AI * (xd * rowQ[mask] + yd * raQ[mask])
                                       + AQ * (xd * rowQq[mask] + yd * rqQ[mask])
                                       + AU * (xd * rowQu[mask] + yd * ruQ[mask]))
            src_u[k, mask] += -ch_c * (AI * (xd * rowU[mask] + yd * raU[mask])
                                       + AQ * (xd * rowUq[mask] + yd * rqU[mask])
                                       + AU * (xd * rowUu[mask] + yd * ruU[mask]))


def sos_pure_fixed_bulk(sza_deg_air, vza_deg_air, raa_deg, lambda_nm,
                        n_water, F_sun, a_tot, b_tot, bb_tot,
                        phase_lut_path, phase_case_id, phase_wl_nm,
                        wind_speed, sigma_type, q_convention,
                        n_mu_water=48, m_max_water=30, n_layers_water=300,
                        tau_max_target=20.0, max_tau_max_target=0.0,
                        max_z_max_m=0.0, depth_bottom_tol=1.0e-8,
                        layer_dtau_target=0.05, max_iterations=10000,
                        tolerance=1.0e-7, view_as_node=True,
                        fixed_bulk_lmax=30, phase_nphi=720,
                        ext_top_I=None, ext_top_Q=None, ext_top_U=None,
                        ext_top_mu=None, ext_top_m_max=-1, constituent=None,
                        view_vza_values=None, raa_values=None, return_lut=False):
    """Fixed-bulk (direct IOP + phase LUT, phase_model=3, value kernel).

    If ``constituent`` is a dict with keys betal/gammal/alphal/zetal/b_w/
    b_particle, the OCRT/CCRR constituent path is used instead: pure-water
    Rayleigh slot (omega_w=b_w/ext, depol 0.039) + particle Greek-coefficient
    slot (omega_particle=b_particle/ext), Greek kernel_phase_fourier (no value
    kernel, no forward truncation, b_rt=b_tot)."""
    res = WaterResult()
    use_constituent = constituent is not None

    mu_sun_air = cos(sza_deg_air * M_PI / 180.0)
    if view_vza_values is None:
        view_vza_values = [float(vza_deg_air)]
    else:
        view_vza_values = [float(v) for v in view_vza_values]
    if raa_values is None:
        raa_values = [float(raa_deg)]
    else:
        raa_values = [float(r) for r in raa_values]
    if not view_vza_values or not raa_values:
        raise ValueError('view_vza_values and raa_values must be non-empty')
    mu_view_air_list = [cos(v * M_PI / 180.0) for v in view_vza_values]
    mu_sun_water = snell_down(mu_sun_air, n_water)
    mu_view_water_list = [snell_down(mu, n_water) for mu in mu_view_air_list]
    mu_view_air = mu_view_air_list[0]
    mu_view_water = mu_view_water_list[0]
    n_view = len(view_vza_values)

    # --- phase setup ---
    if use_constituent:
        # OCRT/CCRR constituent: Greek-coefficient particle phase + pure-water
        # Rayleigh slot.  No forward truncation (water_mie_truncation_mode=0).
        c_betal = np.asarray(constituent['betal'], float)
        c_gammal = np.asarray(constituent['gammal'], float)
        c_alphal = np.asarray(constituent['alphal'], float)
        c_zetal = np.asarray(constituent['zetal'], float)
        b_w = float(constituent['b_w'])
        b_particle = float(constituent['b_particle'])
        b_rt = b_tot
        ext = a_tot + b_rt
        omega = b_rt / ext if ext > 0.0 else 0.0
        omega_w = (b_w / ext) if ext > 0.0 else 0.0
        omega_particle = (b_particle / ext) if ext > 0.0 else 0.0
        L = len(c_betal) - 1
        particle_Lmax = L
        particle_betal = c_betal
        tab = None
    else:
        tab = ph.phase_lut_load(phase_lut_path, phase_case_id, phase_wl_nm)
        A_cap = ph.phase_table_osoaa_cap(tab, 6.0, 3.0)
        delta_f = A_cap
        b_rt = b_tot * (1.0 - delta_f)
        ext = a_tot + b_rt
        omega = b_rt / ext if ext > 0.0 else 0.0
        b_particle = b_rt
        b_w = 0.0
        omega_w = 0.0
        omega_particle = (b_particle / ext) if ext > 0.0 else 0.0

        # particle moments for value kernel: betal=[1,0..0], L
        L = fixed_bulk_lmax if 0 < fixed_bulk_lmax <= 200 else 30
        particle_Lmax = L
        particle_betal = np.zeros(L + 1)
        particle_betal[0] = 1.0

    # --- depth policy ---
    tau_base = tau_max_target if tau_max_target > 0.0 else 20.0
    g_particle = 0.0
    if particle_Lmax >= 1:
        g_particle = particle_betal[1] / 3.0
        g_particle = min(0.999, max(-0.999, g_particle))
    phase_b_sum = b_particle + b_w
    g_eff = (b_particle * g_particle) / phase_b_sum if phase_b_sum > 0.0 else 0.0
    g_eff = min(0.999, max(-0.999, g_eff))
    z_req = tau_base / ext if ext > 0.0 else 0.0
    tol_depth = min(1.0e-2, max(1.0e-30, depth_bottom_tol))
    efolds = -log(tol_depth)
    transport_ext = a_tot + b_rt * max(0.0, 1.0 - g_eff)
    if transport_ext > 0.0:
        z_tr = efolds / transport_ext
        if z_tr > z_req:
            z_req = z_tr
    tau_max = ext * z_req if ext > 0.0 else tau_base
    if max_tau_max_target > 0.0 and tau_max > max_tau_max_target:
        tau_max = max_tau_max_target
        z_req = tau_max / ext if ext > 0.0 else z_req
    if max_z_max_m > 0.0 and z_req > max_z_max_m:
        z_req = max_z_max_m
        tau_max = ext * z_req if ext > 0.0 else tau_max
    z_max = tau_max / ext if ext > 0.0 else z_req

    nt = n_layers_water
    if layer_dtau_target > 0.0 and tau_max > 0.0:
        n_need = int(np.ceil(tau_max / layer_dtau_target - 1e-12))
        # C: ceil(tau/dtau)
        n_need = int(np.ceil(tau_max / layer_dtau_target))
        if n_need > nt:
            nt = n_need
    nt = max(1, nt)

    # --- ring (rt_uangles): GL(n_mu_water) + view node + 3 zero slots ---
    gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_water)
    u = sos.UAngles(gl_mu, gl_w)
    if view_as_node:
        for mu_vw in mu_view_water_list:
            u.add(mu_vw)
    u.add(1.0)
    u.add(mu_sun_water)
    u.add(mu_sun_air)
    n_mu = len(u.mu)

    if use_constituent:
        atm = sos.build_inwater_atm(nt, n_mu, tau_max, omega_w, omega_particle,
                                    mu_sun_water, n_mu, particle_Lmax,
                                    particle_betal, c_gammal, c_alphal, c_zetal,
                                    u.mu, u.w)
    else:
        atm = sos.build_inwater_atm(nt, n_mu, tau_max, omega_w, omega_particle,
                                    mu_sun_water, n_mu, particle_Lmax,
                                    particle_betal, None, None, None,
                                    u.mu, u.w)
    # ring already ascending sorted (UAngles); rm/gb set by build via u lists
    dirs = 2 * n_mu + 1

    # workspace
    ws_l_max = max(atm.L_max, 2)
    ws = kn.LegendreWorkspace(n_mu, ws_l_max)
    m_loop_max = min(m_max_water, ws_l_max)

    # beam_q (flat T_aw Fresnel at air-side solar zenith)
    M_T_aw = sf.T_aw(mu_sun_air, n_water, q_convention)
    atm.beam_q = (M_T_aw[3] / M_T_aw[0]) if M_T_aw[0] != 0.0 else 0.0

    # Value-kernel gates are strictly opt-in.  Gate OFF retains the
    # authoritative coefficient/scalar path bit-for-bit.
    rm_map = {j: atm.rm[j + n_mu] for j in range(-n_mu, n_mu + 1)}
    value_spline = os.environ.get('OCRT_VALUE_PHASE_SPLINE', '0') != '0'
    value_pol = os.environ.get('OCRT_WATER_VALUE_KERNEL_POL', '0') != '0'
    if value_pol:
        _warn_value_kernel_exact_nadir(view_vza_values)
    norm_nmu = int(os.environ.get('OCRT_WATER_MIE_MOMENT_N_MU', '400'))
    value_nphi = int(os.environ.get('OCRT_WATER_VALUE_NPHI', str(phase_nphi)))
    pol_allm = None
    if use_constituent:
        fb_allm = None
        raw = constituent.get('value_phase') if isinstance(constituent, dict) else None
        if value_pol and raw is not None:
            interp = vp.ValuePhaseInterp(raw['theta_deg'], raw['p11'],
                                         raw.get('p12'), raw.get('p33'),
                                         raw.get('norm_n_mu', norm_nmu))
            pol_allm = vp.vector_fourier_allm(atm.rm, n_mu, m_loop_max + 1,
                                              interp, value_nphi)
    else:
        if value_spline:
            interp = vp.ValuePhaseInterp(tab.theta, tab.p11, norm_n_mu=norm_nmu)
            fb_allm = vp.scalar_fourier_allm(atm.rm, n_mu, m_loop_max + 1,
                                             interp, value_nphi)
        else:
            fb_allm = ph.direct_phase_fourier_allm(rm_map, n_mu, m_loop_max + 1, tab,
                                                   phase_nphi)

    # Rww rough-Fresnel microfacet all-m Fourier kernel (wind>0)
    n_phi_ww = max(176, 2 * m_max_water + 16, 128)
    mu_pos = atm.rm[n_mu + 1:].copy()
    Rww_allm = rww_coxmunk_allm(mu_pos, mu_pos, m_loop_max + 1, n_phi_ww,
                                wind_speed, sigma_type, n_water, q_convention)

    Mmodes = m_max_water + 1
    I_m_view_all = np.zeros((n_view, Mmodes))
    Q_m_view_all = np.zeros((n_view, Mmodes))
    U_m_view_all = np.zeros((n_view, Mmodes))
    Iss_m_view_all = np.zeros((n_view, Mmodes))
    I_m_node = np.zeros((Mmodes, n_mu))
    Q_m_node = np.zeros((Mmodes, n_mu))
    U_m_node = np.zeros((Mmodes, n_mu))
    I_m_pos = np.zeros(n_mu)
    I_m_neg = np.zeros(n_mu)
    Q_m_pos = np.zeros(n_mu)
    I_m0_lvl1 = np.zeros(dirs)
    bottom_Ed = bottom_Eu = 0.0

    max_orders_seen = 0
    all_conv = 1
    worst_resid = 0.0
    m_amp_max = np.zeros(n_view)
    m_small_streak = np.zeros(n_view, dtype=int)
    view_done = np.zeros(n_view, dtype=bool)
    view_m_hi = np.full(n_view, -1, dtype=int)
    m_last = -1

    for m in range(0, m_loop_max + 1):
        xpl, xrl, xtl = kn.legendre_compute_pol(ws, atm.rm, m)
        atm.xpl = xpl
        atm.xrl = xrl
        atm.xtl = xtl
        kn.kernel_phase_fourier(ws, m, atm.betal)
        # fixed-bulk overwrites pfm with the value kernel; constituent keeps the
        # Greek-coefficient pfm (rt_kernel_phase_fourier path, use_value_kernel=0)
        if not use_constituent:
            ws.pfm[:, :] = fb_allm[m]
        kn.kernel_phase_fourier_pol(ws, m, atm.gammal)
        kn.kernel_phase_fourier_aerosol_full(ws, m, atm.alphal, atm.zetal)
        if pol_allm is not None:
            ws.pfm[:, :] = pol_allm[0][m]
            ws.gr[:, :]  = pol_allm[1][m]
            ws.gt[:, :]  = pol_allm[2][m]
            ws.arr[:, :] = pol_allm[3][m]
            ws.art[:, :] = pol_allm[4][m]
            ws.att[:, :] = pol_allm[5][m]
        kt = sos.KernelTables(ws.pfm, ws.gr, ws.gt, ws.arr, ws.art, ws.att)

        src_i = sos.primary_source(atm, m, ws, ws.pfm)
        src_q, src_u = sos.primary_source_pol(atm, m, ws, ws.gr, ws.gt)
        if (ext_top_I is not None and ext_top_mu is not None and
                m <= ext_top_m_max):
            _add_diffuse_top_primary(
                atm, ws, m, ext_top_I[m],
                ext_top_Q[m] if ext_top_Q is not None else None,
                ext_top_U[m] if ext_top_U is not None else None,
                ext_top_mu, src_i, src_q, src_u)
        prim_i = sos.integrate_bcs(atm, src_i)
        prim_q = sos.integrate_bcs(atm, src_q)
        prim_u = sos.integrate_bcs(atm, src_u)

        tot_i, tot_q, tot_u, n_ord, conv, resid = sos.sos_pol_intrefl_rough(
            atm, m, kt, prim_i, prim_q, prim_u, Rww_allm[m],
            max_iterations, tolerance)

        # Exact target extraction for every requested VZA.  Each view has
        # its own Fourier early-exit state; the shared SOS field continues
        # until every requested view has converged.
        for iv, mu_t in enumerate(mu_view_water_list):
            if mu_t >= atm.rm[n_mu + n_mu]:
                jl, jh = n_mu - 1, n_mu
            elif mu_t <= atm.rm[n_mu + 1]:
                jl, jh = 1, 2
            else:
                jl, jh = 1, 2
                for jp in range(1, n_mu):
                    if atm.rm[n_mu + jp] <= mu_t < atm.rm[n_mu + jp + 1]:
                        jl, jh = jp, jp + 1
                        break
            ml_, mh_ = atm.rm[n_mu + jl], atm.rm[n_mu + jh]
            whi = (mu_t - ml_) / (mh_ - ml_) if mh_ != ml_ else 0.0
            wlo = 1.0 - whi
            Iss_m_view_all[iv, m] = (wlo * prim_i[0, n_mu + jl] +
                                      whi * prim_i[0, n_mu + jh])
            j_exact = -1
            for jp in range(1, n_mu + 1):
                if abs(atm.rm[n_mu + jp] - mu_t) < 1e-13:
                    j_exact = jp
                    break
            if j_exact > 0:
                I_m_view_all[iv, m] = tot_i[0, n_mu + j_exact]
                Q_m_view_all[iv, m] = tot_q[0, n_mu + j_exact]
                U_m_view_all[iv, m] = tot_u[0, n_mu + j_exact]
            else:
                I_m_view_all[iv, m] = wlo * tot_i[0, n_mu + jl] + whi * tot_i[0, n_mu + jh]
                Q_m_view_all[iv, m] = wlo * tot_q[0, n_mu + jl] + whi * tot_q[0, n_mu + jh]
                U_m_view_all[iv, m] = wlo * tot_u[0, n_mu + jl] + whi * tot_u[0, n_mu + jh]
            if not view_done[iv]:
                amp = max(abs(I_m_view_all[iv, m]), abs(Q_m_view_all[iv, m]),
                          abs(U_m_view_all[iv, m]))
                if amp > m_amp_max[iv]:
                    m_amp_max[iv] = amp
                if m >= 3 and m_amp_max[iv] > 0.0:
                    if amp < tolerance * m_amp_max[iv]:
                        m_small_streak[iv] += 1
                    else:
                        m_small_streak[iv] = 0
                    if m_small_streak[iv] >= 2:
                        view_done[iv] = True
                        view_m_hi[iv] = m
        m_last = m
        # Preserve the single-view ordering: the trigger mode contributes to
        # the target view but not to the exported nodal coupling field.
        if bool(view_done.all()):
            break

        # 7e-bis: per-node field
        I_m_node[m, :] = tot_i[0, n_mu + 1:]
        Q_m_node[m, :] = tot_q[0, n_mu + 1:]
        U_m_node[m, :] = tot_u[0, n_mu + 1:]

        if m == 0:
            I_m_pos[:] = tot_i[0, n_mu + 1:]
            I_m_neg[:] = tot_i[0, :n_mu][::-1]
            Q_m_pos[:] = tot_q[0, n_mu + 1:]
            I_m0_lvl1[:] = tot_i[1, :]
            mu_p = atm.rm[n_mu + 1:]
            w_p = atm.gb[n_mu + 1:]
            bottom_Eu = float((tot_i[nt, n_mu + 1:] * mu_p * w_p).sum()) * 2.0 * M_PI
            bottom_Ed = float((tot_i[nt, :n_mu][::-1] * mu_p * w_p).sum()) * 2.0 * M_PI

    # 8. all-view/all-RAA reconstruction.  The target-view series uses the
    # trigger mode inclusively; water->air nodal coupling uses modes below the
    # trigger, matching the authoritative single-view ordering.
    for iv in range(n_view):
        if view_m_hi[iv] < 0:
            view_m_hi[iv] = m_last
    n_raa = len(raa_values)
    I_diffuse_grid = np.zeros((n_view, n_raa))
    Q_diffuse_grid = np.zeros((n_view, n_raa))
    U_diffuse_grid = np.zeros((n_view, n_raa))
    for iv in range(n_view):
        hi = int(view_m_hi[iv])
        for ir, raa_v in enumerate(raa_values):
            dphi_rad = raa_to_water_view_phi(raa_v)
            I_diffuse_grid[iv, ir] = sos.reconstruct_phi(
                I_m_view_all[iv, :hi + 1], hi, dphi_rad)
            Q_diffuse_grid[iv, ir] = sos.reconstruct_phi(
                Q_m_view_all[iv, :hi + 1], hi, dphi_rad)
            U_diffuse_grid[iv, ir] = -sos.reconstruct_phi_sin(
                U_m_view_all[iv, :hi + 1], hi, dphi_rad)
    I_diffuse_view = I_diffuse_grid[0, 0]
    Q_diffuse_view = Q_diffuse_grid[0, 0]
    U_diffuse_view = U_diffuse_grid[0, 0]

    # 9. hemispheric integrals (m=0)
    mu_p = atm.rm[n_mu + 1:]
    w_p = atm.gb[n_mu + 1:]
    Eu_diffuse = float((I_m_pos * mu_p * w_p).sum()) * 2.0 * M_PI
    Ed_diffuse = float((I_m_neg * mu_p * w_p).sum()) * 2.0 * M_PI
    # Ed_internal_reflect zeroed because Rww_K feedback active
    Ed_internal_reflect = 0.0

    Ed_direct = F_sun * mu_sun_water
    f_scale = F_sun / M_PI
    Ed_diffuse_SI = (Ed_diffuse + Ed_internal_reflect) * f_scale
    Eu_diffuse_SI = Eu_diffuse * f_scale
    Ed_total = Ed_direct + Ed_diffuse_SI
    Eu_total = Eu_diffuse_SI
    I_view = I_diffuse_view * f_scale
    Q_view = Q_diffuse_view * f_scale
    U_view = U_diffuse_view * f_scale

    # 12. Kd, Ku
    Ed_diff_lvl1 = 0.0
    Eu_diff_lvl1 = 0.0
    for jp in range(1, n_mu + 1):
        mu_j = atm.rm[n_mu + jp]
        w_j = atm.gb[n_mu + jp]
        Ed_diff_lvl1 += I_m0_lvl1[n_mu - jp] * mu_j * w_j
        Eu_diff_lvl1 += I_m0_lvl1[n_mu + jp] * mu_j * w_j
    Ed_diff_lvl1 *= 2.0 * M_PI * f_scale
    Eu_diff_lvl1 *= 2.0 * M_PI * f_scale
    tau_lvl1 = tau_max / nt
    Ed_direct_lvl1 = F_sun * mu_sun_water * exp(-tau_lvl1 / mu_sun_water)
    Ed_total_lvl1 = Ed_direct_lvl1 + Ed_diff_lvl1
    Eu_total_lvl1 = Eu_diff_lvl1
    z_lvl1 = tau_lvl1 / ext
    Kd = -log(Ed_total_lvl1 / Ed_total) / max(z_lvl1, 1e-12) \
        if (Ed_total_lvl1 > 0.0 and Ed_total > 0.0) else 0.0
    Ku = -log(Eu_total_lvl1 / Eu_total) / max(z_lvl1, 1e-12) \
        if (Eu_total_lvl1 > 0.0 and Eu_total > 0.0) else 0.0

    # 13. above-water via water->air coupling (wind>0 path, production:
    # OCRT_TWA_MATRIX unset -> analytic flat Snell/Fresnel per-m coupling,
    # mb_re_escape no-op).
    Mc = m_loop_max + 1
    I_air_grid = np.zeros((n_view, n_raa))
    Q_air_grid = np.zeros((n_view, n_raa))
    U_air_grid = np.zeros((n_view, n_raa))
    mu_list = list(mu_pos)
    for iv, mu_a in enumerate(mu_view_air_list):
        n2 = n_water * n_water
        sin2_air = 1.0 - mu_a * mu_a
        mu_w = sqrt(max(0.0, 1.0 - sin2_air / n2))
        M_T_wa = sf.T_wa(mu_w, n_water, q_convention)
        cpl_hi = max(-1, int(view_m_hi[iv]) - 1)
        Ia = np.zeros(Mc); Qa = np.zeros(Mc); Ua = np.zeros(Mc)
        for m in range(min(Mc - 1, cpl_hi) + 1):
            I_w = aw_interp_on_unsorted(mu_list, list(I_m_node[m]), mu_w)
            Q_w = aw_interp_on_unsorted(mu_list, list(Q_m_node[m]), mu_w)
            U_w = aw_interp_on_unsorted(mu_list, list(U_m_node[m]), mu_w)
            Ia[m] = M_T_wa[0] * I_w + M_T_wa[1] * Q_w + M_T_wa[2] * U_w
            Qa[m] = M_T_wa[3] * I_w + M_T_wa[4] * Q_w + M_T_wa[5] * U_w
            Ua[m] = M_T_wa[6] * I_w + M_T_wa[7] * Q_w + M_T_wa[8] * U_w
        for ir, raa_v in enumerate(raa_values):
            phi_v = raa_to_water_view_phi(raa_v)
            LI = LQ = LU = 0.0
            for m in range(min(Mc - 1, cpl_hi) + 1):
                f = 1.0 if m == 0 else 2.0
                pa = phi_v + M_PI
                c = cos(m * pa); sn = sin(m * pa)
                LI += f * Ia[m] * c
                LQ += f * Qa[m] * c
                LU += -f * Ua[m] * sn
            I_air_grid[iv, ir] = LI * f_scale
            Q_air_grid[iv, ir] = LQ * f_scale
            U_air_grid[iv, ir] = LU * f_scale
    I_air = I_air_grid[0, 0]
    Q_air = Q_air_grid[0, 0]
    U_air = U_air_grid[0, 0]

    Ed_0plus_air = F_sun * mu_sun_air

    res.max_orders_used = max_orders_seen
    res.all_converged = all_conv
    res.I_0minus_view = I_view
    res.Q_0minus_view = Q_view
    res.U_0minus_view = U_view
    res.I_0plus_view = I_air
    res.Q_0plus_view = Q_air
    res.U_0plus_view = U_air
    res.Ed_0minus_water = Ed_total
    res.Eu_0minus_water = Eu_total
    res.Ed_0plus_air = Ed_0plus_air
    res.r_rs_0minus = I_view / Ed_total if Ed_total > 0.0 else 0.0
    res.R_rs_0plus = I_air / Ed_0plus_air if Ed_0plus_air > 0.0 else 0.0
    res.Kd_0minus = Kd
    res.Ku_0minus = Ku
    res.a_total_used = a_tot
    res.b_total_used = b_tot
    res.bb_total_used = bb_tot
    res.omega_used = omega
    res.tau_max_used = tau_max
    res.z_max_used = z_max
    res.mu_sun_air = mu_sun_air
    res.mu_sun_water = mu_sun_water
    res.delta_f = 0.0 if use_constituent else delta_f
    res.b_rt = b_rt
    res.ext = ext
    # B.5 export: z=0- upwelling Fourier field (physical = internal x f_scale)
    res.I_up_per_m = I_m_node * f_scale
    res.Q_up_per_m = Q_m_node * f_scale
    res.U_up_per_m = U_m_node * f_scale
    res.mu_water_pos = np.array(mu_pos)
    res.w_water_pos = atm.gb[n_mu + 1:].copy()
    res.n_mu_water_filled = n_mu
    res.m_max_filled = m_max_water
    if return_lut or n_view > 1 or len(raa_values) > 1:
        res.lut = {
            'vza_values': list(view_vza_values), 'raa_values': list(raa_values),
            'I_0minus': I_diffuse_grid * f_scale,
            'Q_0minus': Q_diffuse_grid * f_scale,
            'U_0minus': U_diffuse_grid * f_scale,
            'I_0plus': I_air_grid, 'Q_0plus': Q_air_grid, 'U_0plus': U_air_grid,
            'I_view_per_m': I_m_view_all.copy(),
            'Q_view_per_m': Q_m_view_all.copy(),
            'U_view_per_m': U_m_view_all.copy(),
            'view_m_hi': view_m_hi.copy(),
        }
    return res
