"""Batch orchestrator for GOCI-III production: groups grid cases by common
angular ring (n_mu) and runs R2 (rho_R) / R3 (rho_R+A) / R1 (rho_TOA) as
GPU-capable batches via atmos_batch.

Grouping: a case's ring is GL(n_mu_water) + solar/view specials; the special
nodes are dedup-dropped if they coincide with a GL node (UAngles), so n_mu can
differ per case.  Cases sharing the same (n_mu_atm) — and for R1 also
(n_mu_water_ring, nt_water) — are batched together; the common case (no
coincidence) forms one large group.
"""
import numpy as np
from math import cos, pi, exp
import math

from .backend import xp, to_np
from . import aerosol as AER
from . import kernel as kn
from . import sos as sos1
from .atmos import build_atm_aerosol, build_atm_rayleigh, tau_rayleigh_bodhaine
from . import atmos_batch as AB

M_PI = pi


_AER_PHASE_CACHE = {}


# v1.11 (2026-08-29): per-case CDOM / detritus spectral-slope controls.
# A case dict may carry 'cdom_slope', 'det_a440' and 'det_slope'; absent keys
# fall back to the C defaults (rt_types.h: S_cdom 0.014 nm^-1 Bricaud 1981,
# detritus_slope 0.0109 nm^-1 Bricaud-Stramski, detritus a440 0), so existing
# callers keep bit-identical results.
def _case_slopes(c):
    from .constituent import CDOM_SLOPE_DEFAULT, ORGANIC_DETRITUS_SLOPE_DEFAULT
    return (float(c.get('cdom_slope', CDOM_SLOPE_DEFAULT)),
            float(c.get('det_a440', 0.0)),
            float(c.get('det_slope', ORGANIC_DETRITUS_SLOPE_DEFAULT)))


def _aer_greek(mie, wl_nm, aod, L_max, aer_name=None):
    """Per-case wavelength interp + loglin truncation + Greek coefficients.
    Returns (betal, gammal, alphal, zetal, tau_a_eff, ssa_a_eff, ts, P_norm_asc)
    where the last two feed the value kernel (ascending angle, P11/beta0).

    The phase part (be/ga/al/ze/ssa/A/ts/P_norm_asc) depends only on the
    physical Mie file identity, wavelength and numerical options -- not AOD.
    The cache key therefore uses the absolute file/stat fingerprint populated
    by ``read_mie`` rather than the user-facing aerosol name."""
    mie_key = getattr(mie, 'cache_key', ('memory', id(mie)))
    key = (tuple(mie_key), round(float(wl_nm), 6), int(L_max),
           'pchip-wavelength', 'loglin', 0.8, 0.94, 0.1)
    cached = _AER_PHASE_CACHE.get(key)
    if cached is not None:
        be, ga, al, ze, ssa, A, ts, P_norm_asc = cached
    else:
        wl_um = wl_nm / 1000.0
        P11, P12, P33 = AER.phase_matrix_at_wavelength_pchip(mie, wl_nm)
        ssa = float(np.interp(wl_um, mie.wavelengths, mie.spectral[:, 2]))
        P11t, P12t, P33t, A = AER.loglin_truncate(P11, P12, P33, mie.angles)
        be, ga, al, ze = AER.compute_vector_legendre(P11t, P12t, P33t, mie.angles, L_max, 4096)
        beta0 = AER.compute_legendre_moments(P11t, mie.angles, 1, 4096)[0]
        ts, P11_ts = AER._ascending(mie.angles, P11t)
        P_norm_asc = P11_ts / beta0
        _AER_PHASE_CACHE[key] = (be, ga, al, ze, ssa, A, ts, P_norm_asc)
    tau_a_eff = aod * (1.0 - 0.5 * ssa * A)
    ssa_a_eff = (1.0 - 0.5 * A) * ssa / (1.0 - 0.5 * ssa * A) if A > 0.0 else ssa
    return be, ga, al, ze, tau_a_eff, ssa_a_eff, ts, P_norm_asc


def solve_r3_grid(cases, mie_cache, L_max=80, pressure_hpa=1013.25,
                  n_mu_gl=24, nt=400, m_max=16, max_iterations=100,
                  tolerance=1.0e-7, n_water=1.34, sigma_type=1,
                  q_convention=1, n_phi_quad=1024, aer_h_km=2.0,
                  absorption=None):
    """rho_R+A for a list of cases.  Each case dict: sza,vza,raa,wl,wind,aer,aod.
    mie_cache: {aer_name: Mie}.  absorption: Absorption or None (gas OD applied
    to each medium).  Returns rho_I array aligned with cases."""
    B = len(cases)
    # per-case host scheduling only: aerosol Greek + atm ring + scalars
    prepped = []
    for c in cases:
        mie = mie_cache[c['aer']]
        be, ga, al, ze, tau_a, ssa_e, ts, P_asc = _aer_greek(mie, c['wl'], c['aod'], L_max, c.get('aer'))
        tauR = tau_rayleigh_bodhaine(c['wl'], pressure_hpa, 45.0, 0.0, 360.0)
        mva = cos(c['vza'] * pi / 180.0)
        gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_gl)
        u = sos1.UAngles(list(gl_mu), list(gl_w)); u.add(mva)
        vj = int(np.argmin(np.abs(np.asarray(u.mu) - mva))) + 1
        prepped.append(dict(c=c, be=be, ga=ga, al=al, ze=ze, tau_a=tau_a,
                            ssa_e=ssa_e, ts=ts, P_asc=P_asc, tauR=tauR,
                            uw=u, vj=vj, n_mu=len(u.mu)))
    groups = {}
    for i, p in enumerate(prepped):
        groups.setdefault(p['n_mu'], []).append(i)
    rho = np.zeros(B)
    for n_mu, idxs in groups.items():
        rings = np.stack([np.asarray(prepped[i]['uw'].mu) for i in idxs])
        wts = np.stack([np.asarray(prepped[i]['uw'].w) for i in idxs])
        tau_a_arr = [prepped[i]['tau_a'] for i in idxs]
        ssa_arr = [prepped[i]['ssa_e'] for i in idxs]
        msa_arr = [cos(cases[i]['sza'] * pi / 180.0) for i in idxs]
        tauR0 = prepped[idxs[0]]['tauR']          # shared within band
        wl0 = cases[idxs[0]]['wl']
        aBE = np.stack([np.asarray(prepped[i]['be'])[:L_max + 1] for i in idxs])
        aGA = np.stack([np.asarray(prepped[i]['ga'])[:L_max + 1] for i in idxs])
        aAL = np.stack([np.asarray(prepped[i]['al'])[:L_max + 1] for i in idxs])
        aZE = np.stack([np.asarray(prepped[i]['ze'])[:L_max + 1] for i in idxs])
        bAa, z_km = AB.build_atm_aerosol_batch(nt, tauR0, tau_a_arr, ssa_arr,
                                               msa_arr, rings, wts, aBE, aGA,
                                               aAL, aZE, L_max, aer_h_km=aer_h_km)
        if absorption is not None:
            AB.apply_gas_absorption_batch(bAa, absorption, wl0, z_km)
        vjs = [prepped[i]['vj'] for i in idxs]
        aer_p11 = [(prepped[i]['ts'], prepped[i]['P_asc']) for i in idxs]
        dphis = [cases[i]['raa'] * pi / 180.0 for i in idxs]
        winds = [cases[i]['wind'] for i in idxs]
        res = AB.solve_atm_batch(bAa, vjs, msa_arr, dphis, winds, n_mu, nt,
                                 m_max, max_iterations, tolerance,
                                 mode='aerosol_black_fresnel_ocean', L_max=L_max,
                                 sigma_type=sigma_type, n_water=n_water,
                                 q_convention=q_convention, n_phi_quad=n_phi_quad,
                                 use_value_kernel=True, aer_p11=aer_p11)
        rI = to_np(res['rho_I'])
        for k, i in enumerate(idxs):
            rho[i] = rI[k]
    return rho


def solve_r2_grid(cases, L_max=2, pressure_hpa=1013.25, n_mu_gl=24, nt=40,
                  m_max=2, max_iterations=20, tolerance=1.0e-7, n_water=1.34,
                  sigma_type=1, q_convention=1, n_phi_quad=1024, absorption=None):
    """rho_R (Rayleigh + Cox-Munk) for a list of cases.  Returns rho_I array."""
    B = len(cases)
    prepped = []
    for c in cases:
        tauR = tau_rayleigh_bodhaine(c['wl'], pressure_hpa, 45.0, 0.0, 360.0)
        mva = cos(c['vza'] * pi / 180.0)
        gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_gl)
        u = sos1.UAngles(list(gl_mu), list(gl_w)); u.add(mva)
        vj = int(np.argmin(np.abs(np.asarray(u.mu) - mva))) + 1
        prepped.append(dict(c=c, tauR=tauR, uw=u, vj=vj, n_mu=len(u.mu)))
    groups = {}
    for i, p in enumerate(prepped):
        groups.setdefault(p['n_mu'], []).append(i)
    rho = np.zeros(B)
    for n_mu, idxs in groups.items():
        rings = np.stack([np.asarray(prepped[i]['uw'].mu) for i in idxs])
        wts = np.stack([np.asarray(prepped[i]['uw'].w) for i in idxs])
        msa_arr = [cos(cases[i]['sza'] * pi / 180.0) for i in idxs]
        tauR0 = prepped[idxs[0]]['tauR']          # shared within band
        wl0 = cases[idxs[0]]['wl']
        bA, z_km = AB.build_atm_rayleigh_batch(nt, tauR0, msa_arr, rings, wts)
        if absorption is not None:
            AB.apply_gas_absorption_batch(bA, absorption, wl0, z_km)
        vjs = [prepped[i]['vj'] for i in idxs]
        dphis = [cases[i]['raa'] * pi / 180.0 for i in idxs]
        winds = [cases[i]['wind'] for i in idxs]
        res = AB.solve_atm_batch(bA, vjs, msa_arr, dphis, winds, n_mu, nt,
                                 m_max, max_iterations, tolerance,
                                 mode='black_fresnel_ocean', L_max=L_max,
                                 sigma_type=sigma_type, n_water=n_water,
                                 q_convention=q_convention, n_phi_quad=n_phi_quad)
        rI = to_np(res['rho_I'])
        for k, i in enumerate(idxs):
            rho[i] = rI[k]
    return rho


# ---- R1 (rho_TOA): 3-pass coupled batch ------------------------------------
from .water import snell_down
from .rww import rww_coxmunk_allm, T_aw_coxmunk_direct
from . import surface as sf
from .constituent import OCRTConstituentModel


def _water_depth_policy(a_tot, b_tot, betal, b_particle, b_w,
                        tau_max_target=15.0, max_tau_max_target=200.0,
                        max_z_max_m=200.0, depth_bottom_tol=1.0e-8,
                        n_layers_water=300, layer_dtau_target=0.05):
    """Replicate sos_pure_fixed_bulk depth policy exactly (constituent path:
    b_rt=b_tot, ext=a+b_tot)."""
    b_rt = b_tot
    ext = a_tot + b_rt
    tau_base = tau_max_target if tau_max_target > 0.0 else 20.0
    particle_Lmax = len(betal) - 1
    g_particle = 0.0
    if particle_Lmax >= 1:
        g_particle = betal[1] / 3.0
        g_particle = min(0.999, max(-0.999, g_particle))
    phase_b_sum = b_particle + b_w
    g_eff = (b_particle * g_particle) / phase_b_sum if phase_b_sum > 0.0 else 0.0
    g_eff = min(0.999, max(-0.999, g_eff))
    z_req = tau_base / ext if ext > 0.0 else 0.0
    tol_depth = min(1.0e-2, max(1.0e-30, depth_bottom_tol))
    efolds = -math.log(tol_depth)
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
    nt = n_layers_water
    n_need = nt
    if layer_dtau_target > 0.0 and tau_max > 0.0:
        n_need = int(np.ceil(tau_max / layer_dtau_target))
    ow = (b_w / ext) if ext > 0.0 else 0.0
    op = (b_particle / ext) if ext > 0.0 else 0.0
    return ext, ow, op, tau_max, n_need


TRANS_KEYS = ['T_dir_dn', 'T_diff_dn_hemi', 'T_total_dn_hemi',
              'T_dir_up_view', 'T_diff_up_view', 'T_total_up_view',
              'TOA_water_signal_I', 'T_up_rt_valid',
              'Rrs_I', 'Rrs_Q', 'Rrs_U', 'rrs_I', 'rrs_Q', 'rrs_U']

IOP_KEYS = ['a_total', 'b_total', 'bb_total', 'a_w', 'b_w', 'bb_w',
            'a_chl', 'a_phyto_detritus', 'b_phyto_detritus', 'bb_phyto_detritus',
            'a_dom', 'a_min', 'b_min', 'bb_min', 'Kd0minus']


def solve_r1_grid(cases, mie_cache, constituent_model, L_max=80, Lmix=200,
                  pressure_hpa=1013.25, n_mu_water=12, nt_atm=400, fourier_m_max=16,
                  max_it_atm=100, tol_atm=1e-7, n_layers_water_min=300,
                  max_it_water=500, tol_water=1e-7, n_water=1.34, sigma_type=1,
                  q_convention=1, n_phi_quad=1024, n_phi_taw=128, F_sun_TOA=1.0,
                  aer_h_km=2.0, absorption=None):
    """rho_TOA (coupled) for a list of cases.  Each case dict adds chl,tsm,ad
    (constituent IOP).  constituent_model: OCRTConstituentModel.  absorption:
    Absorption or None (gas OD applied to pass-1/pass-2 medium; tau_abs_col
    folded into tau_atm).  Groups by (n_mu_atm, n_mu_water_ring); common
    nt_water = max(nt_case) per group with bottom identity padding.
    Returns rho_I aligned with cases."""
    from .atmos import apply_gas_absorption
    B = len(cases)
    # ---- batched IOP + cached water moments (band level: all cases share wl) -
    # Within a produce_grid band call every case has the same wl.  The IOP is a
    # single vectorized pass (bit-identical to per-case evaluate).  The water
    # component moments depend on wl only, so they are computed ONCE per
    # component (via the per-case decomposition) and reused for every case; the
    # per-case weighted sum keeps the same accumulation order, so the result is
    # bit-identical to the per-case path while dropping B-1 redundant
    # decompositions per component.
    _iop_np = None; _mom_cache = None
    if len(cases) > 0:
        _wl0 = cases[0]['wl']
        _cs0 = _case_slopes(cases[0])
        # the batched IOP path shares one wavelength AND one slope set; mixed
        # per-case slopes fall through to the per-case evaluate() below.
        if (all(c['wl'] == _wl0 for c in cases) and
                all(_case_slopes(c) == _cs0 for c in cases)):
            _chl = np.array([c['chl'] for c in cases], dtype=float)
            _tsm = np.array([c['tsm'] for c in cases], dtype=float)
            _acd = np.array([c['ad'] for c in cases], dtype=float)
            _iopb = constituent_model.evaluate_batch(
                _wl0, _chl, _tsm, _acd, cdom_slope=_cs0[0],
                det_a440=_cs0[1], det_slope=_cs0[2])
            _iop_np = {k: to_np(v) for k, v in _iopb.items()}
            from .phase_batch import water_component_phase_moments_gauss_vec
            # Stage-2 constituent Chl scattering is disabled.  Preserve a
            # three-slot layout for deterministic mixing, but never prepare
            # the phytoplankton phase moments.
            _mom_cache = [None]
            for mie_c in (constituent_model.det_mie, constituent_model.min_mie):
                bb_, gg, aa, zz, _ = water_component_phase_moments_gauss_vec(
                    mie_c, _wl0, Lmix, 400)
                _mom_cache.append((bb_, gg, aa, zz))
    # ---- per-case prep: atm medium (pass-1), water medium, coupling scalars --
    prep = []
    for i, c in enumerate(cases):
        mie = mie_cache[c['aer']]
        be, ga, al, ze, tau_a, ssa_e, _ts, _pasc = _aer_greek(mie, c['wl'], c['aod'], L_max, c.get('aer'))
        tauR = tau_rayleigh_bodhaine(c['wl'], pressure_hpa, 45.0, 0.0, 360.0)
        mu_sun_air = cos(c['sza'] * pi / 180.0); mu_v_air = cos(c['vza'] * pi / 180.0)
        # atm ring (GL + view node); the medium is built batched per group in
        # _r1_group, so here we only assemble host scheduling data.
        gla_mu, gla_w = kn.gauss_legendre_pos(n_mu_water)
        atm_uw = sos1.UAngles(list(gla_mu), list(gla_w))
        atm_uw.add(mu_v_air)
        vj = int(np.argmin(np.abs(np.asarray(atm_uw.mu) - mu_v_air))) + 1
        n_mu_atm = len(atm_uw.mu)
        tau_abs_col = 0.0
        if absorption is not None:
            tau_abs_col = float(absorption.tau_total(c['wl']))
        # constituent IOP + Greek (batched IOP + cached moments if available)
        if _iop_np is not None:
            r = {k: _iop_np[k][i] for k in _iop_np}
            bp = r['b_phyto'] + r['b_det'] + r['b_min']
            be2 = np.zeros(Lmix + 1); ga2 = np.zeros(Lmix + 1)
            al2 = np.zeros(Lmix + 1); ze2 = np.zeros(Lmix + 1)
            if bp > 0.0:
                for moments, bc in zip(
                        _mom_cache, (r['b_phyto'], r['b_det'], r['b_min'])):
                    if moments is None or not (bc > 0.0):
                        continue
                    bb_, gg, aa, zz = moments
                    wgt = bc / bp
                    be2 += wgt * bb_; ga2 += wgt * gg
                    al2 += wgt * aa; ze2 += wgt * zz
        else:
            _cs = _case_slopes(c)
            r = constituent_model.evaluate(c['wl'], c['chl'], c['tsm'], c['ad'],
                                           cdom_slope=_cs[0], det_a440=_cs[1],
                                           det_slope=_cs[2])
            bp = r['b_phyto'] + r['b_det'] + r['b_min']
            be2 = np.zeros(Lmix + 1); ga2 = np.zeros(Lmix + 1)
            al2 = np.zeros(Lmix + 1); ze2 = np.zeros(Lmix + 1)
            from .phase_batch import water_component_phase_moments_gauss_vec
            if bp > 0.0:
                for mie_c, bc in [(constituent_model.det_mie, r['b_det']),
                                  (constituent_model.min_mie, r['b_min'])]:
                    if not (bc > 0.0):
                        continue
                    bb_, gg, aa, zz, _ = water_component_phase_moments_gauss_vec(
                        mie_c, c['wl'], Lmix, 400)
                    wgt = bc / bp
                    be2 += wgt * bb_; ga2 += wgt * gg
                    al2 += wgt * aa; ze2 += wgt * zz
        mu_sun_water = snell_down(mu_sun_air, n_water)
        mu_view_water = snell_down(mu_v_air, n_water)
        ext, ow, op, tau_max, n_need = _water_depth_policy(r['a'], r['b'], be2, bp, r['b_w'])
        nt_case = max(n_layers_water_min, n_need)
        gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_water)
        uw = sos1.UAngles(gl_mu, gl_w)
        uw.add(mu_view_water); uw.add(1.0); uw.add(mu_sun_water); uw.add(mu_sun_air)
        n_mu_w = len(uw.mu)
        prep.append(dict(c=c, atm_uw=atm_uw, vj=vj, n_mu_atm=n_mu_atm,
                         aer_be=be, aer_ga=ga, aer_al=al, aer_ze=ze,
                         ssa_e=ssa_e, tauR=tauR, tau_a=tau_a, tau_abs_col=tau_abs_col,
                         mu_sun_air=mu_sun_air, mu_v_air=mu_v_air,
                         mu_sun_water=mu_sun_water, mu_view_water=mu_view_water,
                         r=r, bp=bp, be2=be2, ga2=ga2, al2=al2, ze2=ze2,
                         ow=ow, op=op, tau_max=tau_max, n_need=n_need, nt_case=nt_case,
                         n_mu_w=n_mu_w, uw=uw))
    # ---- group by (n_mu_atm, n_mu_w); pad water to common ntw = max nt_case --
    groups = {}
    for i, p in enumerate(prep):
        groups.setdefault((p['n_mu_atm'], p['n_mu_w']), []).append(i)
    rho = np.zeros(B)
    rrs = np.zeros(B)
    trans = {k: np.zeros(B) for k in TRANS_KEYS}
    iop = {k: np.zeros(B) for k in IOP_KEYS}
    for (n_mu_atm, n_mu_w), idxs in groups.items():
        ntw = max(prep[i]['nt_case'] for i in idxs)
        _r1_group(prep, idxs, n_mu_atm, n_mu_w, ntw, Lmix, nt_atm, fourier_m_max,
                  max_it_atm, tol_atm, max_it_water, tol_water, n_water, sigma_type,
                  q_convention, n_phi_quad, n_phi_taw, F_sun_TOA, L_max, n_mu_water,
                  absorption, aer_h_km, rho, rrs, trans, iop)
    return rho, rrs, trans, iop


def _pad_water_medium(atmw, ntw):
    """Pad an in-water medium at the bottom with identity layers up to ntw
    (h held constant -> dtau=0 -> transmission=1; xdel/ydel=0 -> no source)."""
    nt0 = atmw.n_layers
    if nt0 >= ntw:
        return atmw
    h = xp.zeros(ntw + 1); h[:nt0 + 1] = atmw.h; h[nt0 + 1:] = atmw.h[nt0]
    ch = xp.zeros(ntw + 1); ch[:nt0 + 1] = atmw.ch
    xd = xp.zeros(ntw + 1); xd[:nt0 + 1] = atmw.xdel
    yd = xp.zeros(ntw + 1); yd[:nt0 + 1] = atmw.ydel
    atmw.h = h; atmw.ch = ch; atmw.xdel = xd; atmw.ydel = yd
    atmw.n_layers = ntw
    return atmw


def _r1_group(prep, idxs, n_mu_atm, n_mu_w, ntw, Lmix, nt_atm, fmm, mia, ta,
              miw, tw, n_water, sigma_type, q_convention, n_phi_quad, n_phi_taw,
              F_sun_TOA, L_max, n_mu_water, absorption, aer_h_km, rho_out, rrs_out,
              trans_out, iop_out):
    """Run one common-ring group through the 3-pass coupled pipeline.
    n_mu_water: plain GL count used for the coupling grid.
    rrs_out is retained as the backward-compatible Rrs0plus-I array.
    All Rrs(0+) and rrs(0-) I/Q/U components are also written to trans_out.
    They are algebraic views of fields already solved by the R1 pipeline."""
    G = len(idxs)
    # ---- batched in-atmosphere medium build (single GPU pass, no per-case) --
    # atm rings differ per case (view node) but share size within a group.
    atm_rings = np.stack([np.asarray(prep[i]['atm_uw'].mu) for i in idxs])   # (G, n_mu_atm)
    atm_wts = np.stack([np.asarray(prep[i]['atm_uw'].w) for i in idxs])
    tau_a_arr = [prep[i]['tau_a'] for i in idxs]
    ssa_arr = [prep[i]['ssa_e'] for i in idxs]
    msa_arr = [prep[i]['mu_sun_air'] for i in idxs]
    tauR0 = prep[idxs[0]]['tauR']          # shared within band (wl, pressure fixed)
    wl0 = prep[idxs[0]]['c']['wl']
    aBE = np.stack([np.asarray(prep[i]['aer_be'])[:L_max + 1] for i in idxs])
    aGA = np.stack([np.asarray(prep[i]['aer_ga'])[:L_max + 1] for i in idxs])
    aAL = np.stack([np.asarray(prep[i]['aer_al'])[:L_max + 1] for i in idxs])
    aZE = np.stack([np.asarray(prep[i]['aer_ze'])[:L_max + 1] for i in idxs])
    bAa, z_km_b = AB.build_atm_aerosol_batch(nt_atm, tauR0, tau_a_arr, ssa_arr,
                                             msa_arr, atm_rings, atm_wts,
                                             aBE, aGA, aAL, aZE, L_max, aer_h_km=aer_h_km)
    if absorption is not None:
        AB.apply_gas_absorption_batch(bAa, absorption, wl0, z_km_b)
    vjs = [prep[i]['vj'] for i in idxs]
    mu_suns = [prep[i]['mu_sun_air'] for i in idxs]
    dphis = [prep[i]['c']['raa'] * pi / 180.0 for i in idxs]
    winds = [prep[i]['c']['wind'] for i in idxs]
    n_mu = n_mu_atm

    # ---- pass-1 atm (aerosol + Cox-Munk), full boa export ------------------
    res1 = AB.solve_atm_batch(bAa, vjs, mu_suns, dphis, winds, n_mu, nt_atm,
                              fmm, mia, ta, mode='aerosol_black_fresnel_ocean', L_max=L_max,
                              n_phi_quad=n_phi_quad)
    atm_path_I = res1['I_TOA']                        # (G,) xp
    boaI = res1['boa_I']; boaQ = res1['boa_Q']; boaU = res1['boa_U']
    mu_arr = res1['mu_quad']; w_arr = res1['w_quad']  # (G, n_a) xp (atm ring)
    mu_atm_list = [mu_arr[k] for k in range(G)]        # xp per-case (kernel builders)
    w_atm_list = [w_arr[k] for k in range(G)]

    # ---- forward coupling (Cox-Munk Fourier), pi-normalized snapshot -------
    mu_wg, w_wg = kn.gauss_legendre_pos(n_mu_water)
    cIb, cQb, cUb = AB.couple_atm_to_water_batch(
        boaI, boaQ, boaU, mu_atm_list, fmm, n_water, q_convention,
        [np.asarray(mu_wg)] * G, w_atm_list, winds, sigma_type, n_phi_taw)

    # ---- per-case setup scalars (host, no GPU read) + vectorized Ed/rescale -
    mus_arr = np.array([prep[i]['mu_sun_air'] for i in idxs])
    musw_arr = np.array([prep[i]['mu_sun_water'] for i in idxs])
    tau_atm_arr = np.array([prep[i]['tauR'] + prep[i]['tau_a'] + prep[i]['tau_abs_col']
                            for i in idxs])
    Fboa_arr = np.exp(-tau_atm_arr / mus_arr)
    Taw_arr = np.array([min(1.0, max(0.0, T_aw_coxmunk_direct(
        prep[i]['mu_sun_air'], n_water, prep[i]['c']['wind'], sigma_type)))
        for i in idxs])
    Fw_arr = Fboa_arr * mus_arr * Taw_arr / musw_arr
    dtsc_arr = np.where(Fw_arr > 0.0, F_sun_TOA / np.where(Fw_arr > 0.0, Fw_arr, 1.0), 0.0)
    # Ed_diff_BOA (sky diffuse downwelling, hemispheric integral of boa m=0), xp
    ed_diff = 2.0 * F_sun_TOA * xp.sum(w_arr * mu_arr * boaI[:, 0, :], axis=1)   # (G,)
    Ed_0plus = xp.asarray(Fboa_arr * mus_arr) + ed_diff                          # (G,)
    dtsc_xp = xp.asarray(dtsc_arr)[:, None, None]
    dtI = cIb * dtsc_xp; dtQ = cQb * dtsc_xp; dtU = cUb * dtsc_xp                 # xp
    # DTPSIGN: atmosphere BOA modes use phi+pi, water diffuse-top uses phi.
    # Flip odd m at the hand-off only; keep the coupling operator itself raw.
    from . import coupling as _cpl_basis
    _cpl_basis.apply_diffuse_top_basis_sign_inplace(dtI, dtQ, dtU, mode_axis=1)
    F_sun_water = Fw_arr                                # (G,) host, passed as list

    # ---- collect per-case water inputs (host scheduling only, no GPU) -------
    rings_list = [np.asarray(prep[i]['uw'].mu) for i in idxs]
    wts_list = [np.asarray(prep[i]['uw'].w) for i in idxs]
    ntc_list = [prep[i]['nt_case'] for i in idxs]
    tau_list = [prep[i]['tau_max'] for i in idxs]
    ow_list = [prep[i]['ow'] for i in idxs]
    op_list = [prep[i]['op'] for i in idxs]
    msw_list = [prep[i]['mu_sun_water'] for i in idxs]
    beamq_list = []
    for i in idxs:
        M = sf.T_aw(prep[i]['mu_sun_air'], n_water, q_convention)
        beamq_list.append((M[3] / M[0]) if M[0] != 0.0 else 0.0)
    BE = np.stack([np.asarray(prep[i]['be2'])[:Lmix + 1] for i in idxs])
    GA = np.stack([np.asarray(prep[i]['ga2'])[:Lmix + 1] for i in idxs])
    AL = np.stack([np.asarray(prep[i]['al2'])[:Lmix + 1] for i in idxs])
    ZE = np.stack([np.asarray(prep[i]['ze2'])[:Lmix + 1] for i in idxs])
    mu_stack = np.stack(rings_list)                                    # (G, n_mu_w)
    wt_stack = np.stack(wts_list)
    # ---- batched surface reflection kernel over the group ------------------
    # Rings differ per case (view/sun nodes) but share size within a group, so
    # rww is a single batched pass, chunked so the intermediate
    # (chunk, n, n, n_phi, 9) Mueller buffer stays bounded on the GPU.
    from .rww_batch import rww_coxmunk_allm_batch
    n_phi_rww = max(176, 2 * fmm + 16, 128)
    ws_stack = np.array([prep[i]['c']['wind'] for i in idxs], dtype=float)
    _RWW_CHUNK = 16
    _rc = []
    for s in range(0, G, _RWW_CHUNK):
        e = min(s + _RWW_CHUNK, G)
        _rc.append(rww_coxmunk_allm_batch(mu_stack[s:e], mu_stack[s:e], fmm + 1,
                   n_phi_rww, ws_stack[s:e], sigma_type, n_water, q_convention))
    Rww_all = _rc[0] if len(_rc) == 1 else xp.concatenate(_rc, axis=0)    # (G,M,n,n,9)
    # ---- batched in-water medium build (single GPU pass, no per-case loop) --
    bAw = AB.build_inwater_atm_batch(ntc_list, tau_list, ow_list, op_list,
                                     msw_list, mu_stack, wt_stack, BE, GA, AL, ZE,
                                     ntw, Lmix, beamq_list)
    Rww_b = [Rww_all[:, m] for m in range(fmm + 1)]   # (G, n, n, 9) per mode m
    mu_sun_water_list = [prep[i]['mu_sun_water'] for i in idxs]
    mu_view_water_list = [prep[i]['mu_view_water'] for i in idxs]
    raa_list = [prep[i]['c']['raa'] for i in idxs]
    resW = AB.solve_water_batch_r1(bAw, mu_sun_water_list, mu_view_water_list,
                                   Rww_b, fmm, miw, tw, list(F_sun_water),
                                   ext_top=(dtI, dtQ, dtU, [np.asarray(mu_wg)] * G),
                                   raa_list=raa_list, n_water=n_water,
                                   q_convention=q_convention,
                                   ext_list=[prep[i]['r']['a'] + prep[i]['r']['b']
                                             for i in idxs],
                                   ext_w_list=[np.asarray(w_wg)] * G)
    Iup = resW['I_up_per_m']; Qup = resW['Q_up_per_m']; Uup = resW['U_up_per_m']   # xp
    Lu_0plus = resW['Lu_0plus']    # (G,) xp — water-leaving radiance above surface
    Qu_0plus = resW['Qu_0plus']; Uu_0plus = resW['Uu_0plus']
    Lu_0minus = resW['Lu_0minus']; Qu_0minus = resW['Qu_0minus']; Uu_0minus = resW['Uu_0minus']

    # ---- reverse coupling + pass-2 atm (bottom source) ---------------------
    wlI, wlQ, wlU = AB.couple_water_to_atm_batch(Iup, Qup, Uup, rings_list,
                                                 fmm, n_water, q_convention, mu_atm_list)
    res2 = AB.solve_atm_batch(bAa, vjs, mu_suns, dphis, winds, n_mu, nt_atm,
                              fmm, mia, ta, mode='aerosol_black_fresnel_ocean', L_max=L_max,
                              n_phi_quad=n_phi_quad,
                              bottom_source=(wlI, wlQ, wlU))       # xp
    TOA_wl = res2['I_TOA']                              # (G,) xp

    # ---- assembly (vectorized; single host readback of rho/rrs) ------------
    mus_xp = xp.asarray(mus_arr)
    rho_arr = atm_path_I / mus_xp + TOA_wl * pi / (F_sun_TOA * mus_xp)           # (G,) xp
    den0p = xp.where(Ed_0plus != 0.0, Ed_0plus, 1.0)
    den0m = xp.where(resW['Ed_0minus'] != 0.0, resW['Ed_0minus'], 1.0)
    RrsI_arr = xp.where(Ed_0plus > 0.0, Lu_0plus / den0p, 0.0)
    RrsQ_arr = xp.where(Ed_0plus != 0.0, Qu_0plus / den0p, 0.0)
    RrsU_arr = xp.where(Ed_0plus != 0.0, Uu_0plus / den0p, 0.0)
    rrsI_arr = xp.where(resW['Ed_0minus'] > 0.0, Lu_0minus / den0m, 0.0)
    rrsQ_arr = xp.where(resW['Ed_0minus'] != 0.0, Qu_0minus / den0m, 0.0)
    rrsU_arr = xp.where(resW['Ed_0minus'] != 0.0, Uu_0minus / den0m, 0.0)
    rho_np = to_np(rho_arr); rrs_np = to_np(RrsI_arr)
    RrsQ_np = to_np(RrsQ_arr); RrsU_np = to_np(RrsU_arr)
    rrsI_np = to_np(rrsI_arr); rrsQ_np = to_np(rrsQ_arr); rrsU_np = to_np(rrsU_arr)
    # ---- RT 유도 대기투과율 (C v1.2 rt-derived-transmittance fix와 동일 정의) --
    # 하향: T_dir_dn=exp(-tau/mu_sun), T_diff_dn_hemi=대기 SOS BOA 확산 반구적분,
    #       T_total_dn_hemi=Ed(0+)/(F_sun·mu_sun)=직달+확산.
    # 상향: T_dir_up_view=exp(-tau/mu_view), T_total_up_view=(TOA 수광신호)/Lu(0+)
    #       =RT 유효 상향투과율(상반성 근사 대신 실제 bottom-source pass 비율),
    #       T_diff_up_view=T_total_up_view-T_dir_up_view(각도재분배로 음수 가능).
    # TOA_water_signal_I=I_TOA_total-I_atm_path=TOA 수광 radiance(=TOA_wl).
    # T_up_rt_valid=1(수광 source 존재). Lu(0+)=0이면 상향은 NaN, 플래그 0.
    # rho_wn=t_rho_w/(T_total_dn_hemi·T_total_up_view)=표층 water-leaving 반사도.
    Ed0 = to_np(Ed_0plus); Lu0 = to_np(Lu_0plus); TOAwl = to_np(TOA_wl)
    Kd_np = to_np(resW['Kd0minus'])
    Tdiff_dn = to_np(res1['T_diff_dn_hemi'])
    muv_arr = np.array([prep[i]['mu_v_air'] for i in idxs])
    T_dir_dn_np   = np.asarray(Fboa_arr, dtype=float)          # exp(-tau/mu_sun)
    T_totdn_np    = Ed0 / (F_sun_TOA * mus_arr)                # direct + diffuse
    T_dir_up_np   = np.exp(-tau_atm_arr / muv_arr)             # exp(-tau/mu_view)
    valid = np.abs(Lu0) > 0.0
    T_totup_np    = np.where(valid, TOAwl / np.where(valid, Lu0, 1.0), np.nan)
    T_diff_up_np  = T_totup_np - T_dir_up_np
    for k, i in enumerate(idxs):
        rho_out[i] = rho_np[k]
        rrs_out[i] = rrs_np[k]
        trans_out['T_dir_dn'][i]           = T_dir_dn_np[k]
        trans_out['T_diff_dn_hemi'][i]     = Tdiff_dn[k]
        trans_out['T_total_dn_hemi'][i]    = T_totdn_np[k]
        trans_out['T_dir_up_view'][i]      = T_dir_up_np[k]
        trans_out['T_diff_up_view'][i]     = T_diff_up_np[k]
        trans_out['T_total_up_view'][i]    = T_totup_np[k]
        trans_out['TOA_water_signal_I'][i] = TOAwl[k]
        trans_out['T_up_rt_valid'][i]      = 1.0 if valid[k] else 0.0
        trans_out['Rrs_I'][i] = rrs_np[k]
        trans_out['Rrs_Q'][i] = RrsQ_np[k]
        trans_out['Rrs_U'][i] = RrsU_np[k]
        trans_out['rrs_I'][i] = rrsI_np[k]
        trans_out['rrs_Q'][i] = rrsQ_np[k]
        trans_out['rrs_U'][i] = rrsU_np[k]
        # 구성성분 IOP (C rt_water_rt: a_pig=a_phyto+a_det; a_chl=a_phyto 별도)
        ri = prep[i]['r']
        iop_out['a_total'][i]  = ri['a'];   iop_out['b_total'][i]  = ri['b']
        iop_out['bb_total'][i] = ri['bb']
        iop_out['a_w'][i]  = ri['a_w'];  iop_out['b_w'][i]  = ri['b_w']
        iop_out['bb_w'][i] = ri['bb_w']
        iop_out['a_chl'][i]            = ri['a_phyto']
        iop_out['a_phyto_detritus'][i]  = ri['a_phyto'] + ri['a_det']
        iop_out['b_phyto_detritus'][i]  = ri['b_phyto'] + ri['b_det']
        iop_out['bb_phyto_detritus'][i] = ri['bb_phyto'] + ri['bb_det']
        iop_out['a_dom'][i] = ri['a_cdom']
        iop_out['a_min'][i]  = ri['a_min'];  iop_out['b_min'][i]  = ri['b_min']
        iop_out['bb_min'][i] = ri['bb_min']
        iop_out['Kd0minus'][i] = Kd_np[k]
