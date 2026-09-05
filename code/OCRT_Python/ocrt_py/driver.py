"""rt_solve_case_ocean + rt_io_print_single transliteration for the
no-atmosphere (pressure 0, aod 0) fixed-bulk route."""
from math import pi, cos, exp
import os

from . import surface as sf
from .rww import T_aw_coxmunk_direct
from .absorption import Absorption
from .water import sos_pure_fixed_bulk, snell_down
from .spectral_contract import require_wavelength

M_PI = pi


def solve_case_ocean_coupled(sza_deg, vza_deg, raa_deg, wavelength_nm,
                             n_water, wind_speed, a_tot, b_tot, bb_tot,
                             phase_lut_path, mie, user_aod,
                             pressure_hpa=1013.25, absorption=None,
                             F_sun_TOA=1.0,
                             n_mu_water=48, fourier_m_max=16, nt_atm=400,
                             L_max=80, max_it_atm=100, tol_atm=1.0e-7,
                             m_max_water=30, n_layers_water=300,
                             max_it_water=160, tol_water=1.0e-7,
                             sigma_type=1, q_convention=1, n_phi_quad=1024,
                             n_phi_taw=128, constituent=None,
                             aerosol_runtime=None):
    """rt_solve_case_ocean, pristine default path (diffuse-top + C3d), for
    the fixed-bulk IOP branch with aerosol atmosphere and wind>0.

    Flow: atm pass-1 (COXMUNK, shared-N ring, boa export) -> forward Cox-Munk
    Fourier coupling -> pi-normalized diffuse-top snapshot (x F_TOA/F_sun_w)
    -> one water solve (direct beam + injected diffuse top) -> flat reverse
    coupling of the 0- per-m upwelling -> atm pass-2 (bottom_source_only)
    -> assembly.  Skylight first-cut terms are zero in diffuse-top mode.
    Sunglint is always decoupled in rho_I and added back in rho_I_glint.
    """
    wavelength_nm = require_wavelength(wavelength_nm, context='coupled ocean case')
    if wind_speed <= 0.0:
        raise NotImplementedError('flat-surface (wind=0) pass-2 not ported')
    if user_aod > 0.0 and aerosol_runtime is None and mie is None:
        raise ValueError(
            'Aerosol runtime is required when AOD > 0 in ocean full-grid mode.')
    from .atmos import solve_aerosol_black_fresnel_ocean_value
    from . import coupling as cpl
    from .kernel import gauss_legendre_pos
    import numpy as np

    mu_sun_air = cos(sza_deg * M_PI / 180.0)
    mu_v_air = cos(vza_deg * M_PI / 180.0)

    # ---- atm pass-1: path radiance + BOA skylight export ------------------
    r1 = solve_aerosol_black_fresnel_ocean_value(
        sza_deg, vza_deg, raa_deg, wavelength_nm, wind_speed, mie, user_aod,
        L_max=L_max, pressure_hpa=pressure_hpa, n_mu_gl=n_mu_water,
        nt=nt_atm, m_max=fourier_m_max, max_iterations=max_it_atm,
        tolerance=tol_atm, n_water=n_water, sigma_type=sigma_type,
        q_convention=q_convention, n_phi_quad=n_phi_quad,
        surface='black_fresnel_ocean', absorption=absorption, aer_h_km=2.0,
        use_value_kernel=True, full_result=True,
        aerosol_runtime=aerosol_runtime)
    atm_path_I = r1['I_TOA']
    atm_path_Q = r1['Q_TOA']
    atm_path_U = r1['U_TOA']
    boa = r1['boa']
    mu_atm = boa['mu_quad_pos']
    w_atm = boa['w_quad_pos']

    # ---- Stage 2a: direct-beam attenuation --------------------------------
    tau_abs_col = absorption.tau_total(wavelength_nm) if absorption is not None else 0.0
    tau_atm = r1['tau_R'] + r1['tau_a_eff'] + tau_abs_col
    F_sun_BOA = F_sun_TOA
    if mu_sun_air > 0.0 and tau_atm > 0.0:
        F_sun_BOA = F_sun_TOA * exp(-tau_atm / mu_sun_air)
    T_dir_beam = exp(-tau_atm / mu_sun_air) if (mu_sun_air > 0.0 and tau_atm > 0.0) else 1.0

    # ---- coupling grid + forward coupling (wind>0: kernel path) -----------
    mu_wg, w_wg = gauss_legendre_pos(n_mu_water)   # plain [0,1] GL (wind>0)
    cI, cQ, cU = cpl.couple_atm_to_water(
        boa['I_per_m'], boa['Q_per_m'], boa['U_per_m'], mu_atm,
        fourier_m_max, n_water, q_convention, np.asarray(mu_wg), w_atm,
        wind_speed, sigma_type, n_phi_kernel=n_phi_taw)
    dtI, dtQ, dtU = cI.copy(), cQ.copy(), cU.copy()   # pi-normalized snapshot
    cpl.apply_diffuse_top_basis_sign_inplace(dtI, dtQ, dtU, mode_axis=0)

    # ---- Ed integrals on the unified air ring -----------------------------
    Ed_diff_water_from_atm = 0.0
    Ed_diff_BOA_air = 0.0
    for j in range(len(mu_atm)):
        Td = T_aw_coxmunk_direct(mu_atm[j], n_water, wind_speed, sigma_type)
        muL = w_atm[j] * mu_atm[j] * boa['I_per_m'][0, j]
        Ed_diff_water_from_atm += muL * Td
        Ed_diff_BOA_air += muL
    Ed_diff_water_from_atm *= 2.0 * M_PI * (F_sun_TOA / M_PI)
    Ed_diff_BOA_air *= 2.0 * M_PI * (F_sun_TOA / M_PI)

    # ---- F_sun_water + diffuse-top rescale --------------------------------
    mu_sun_water = snell_down(mu_sun_air, n_water)
    F_sun_water = F_sun_BOA
    if mu_sun_air > 0.0 and mu_sun_water > 0.0:
        T_aw = T_aw_coxmunk_direct(mu_sun_air, n_water, wind_speed, sigma_type)
        T_aw = min(1.0, max(0.0, T_aw))
        F_sun_water = F_sun_BOA * mu_sun_air * T_aw / mu_sun_water
        if F_sun_water > 0.0:
            dtsc = F_sun_TOA / F_sun_water
            dtI *= dtsc
            dtQ *= dtsc
            dtU *= dtsc

    # ---- water solve (direct beam + injected diffuse top) -----------------
    w = sos_pure_fixed_bulk(sza_deg, vza_deg, raa_deg, wavelength_nm,
                            n_water, F_sun_water, a_tot, b_tot, bb_tot,
                            phase_lut_path, None, 0.0,
                            wind_speed, sigma_type, q_convention,
                            n_mu_water=n_mu_water, m_max_water=m_max_water,
                            n_layers_water=n_layers_water,
                            max_iterations=max_it_water, tolerance=tol_water,
                            ext_top_I=dtI, ext_top_Q=dtQ, ext_top_U=dtU,
                            ext_top_mu=np.asarray(mu_wg),
                            ext_top_m_max=fourier_m_max, constituent=constituent)

    # ---- C3d: reverse coupling + atm pass-2 (bottom source only) ----------
    mmax = min(fourier_m_max, w.m_max_filled)
    wlI, wlQ, wlU = cpl.couple_water_to_atm(
        w.I_up_per_m, w.Q_up_per_m, w.U_up_per_m, w.mu_water_pos, mmax,
        n_water, q_convention, mu_atm)
    r2 = solve_aerosol_black_fresnel_ocean_value(
        sza_deg, vza_deg, raa_deg, wavelength_nm, wind_speed, mie, user_aod,
        L_max=L_max, pressure_hpa=pressure_hpa, n_mu_gl=n_mu_water,
        nt=nt_atm, m_max=fourier_m_max, max_iterations=max_it_atm,
        tolerance=tol_atm, n_water=n_water, sigma_type=sigma_type,
        q_convention=q_convention, n_phi_quad=n_phi_quad,
        surface='black_fresnel_ocean', absorption=absorption, aer_h_km=2.0,
        use_value_kernel=True, bottom_source=(wlI, wlQ, wlU),
        full_result=True, aerosol_runtime=aerosol_runtime)
    TOA_wl_direct_I = r2['I_TOA']
    TOA_wl_direct_Q = r2['Q_TOA']
    TOA_wl_direct_U = r2['U_TOA']

    # ---- assembly (diffuse-top: sky first-cut terms are zero) -------------
    Ed_0plus_air = F_sun_BOA * mu_sun_air if mu_sun_air > 0.0 else w.Ed_0plus_air
    Lu_total_below = w.I_0minus_view
    Qu_total_below = w.Q_0minus_view
    Uu_total_below = w.U_0minus_view
    Ed_total_below = w.Ed_0minus_water + Ed_diff_water_from_atm
    Ed_above_air = Ed_0plus_air + Ed_diff_BOA_air
    Lu_total_above = w.I_0plus_view
    Qu_total_above = w.Q_0plus_view
    Uu_total_above = w.U_0plus_view

    mu_v_for_T = mu_v_air if mu_v_air > 0.0 else 1.0
    T_atm_dir_up_view = exp(-tau_atm / mu_v_for_T) if tau_atm > 0.0 else 1.0
    T_atm_diff_up_view = r1['T_diff_dn_hemi']
    T_atm_up_view = T_atm_dir_up_view + T_atm_diff_up_view

    I_TOA_total = atm_path_I + TOA_wl_direct_I
    Q_TOA_total = atm_path_Q + TOA_wl_direct_Q
    U_TOA_total = atm_path_U + TOA_wl_direct_U

    rho_factor = (M_PI / (F_sun_TOA * mu_sun_air)) if (F_sun_TOA > 0.0 and mu_sun_air > 0.0) else 0.0
    atm_refl_norm = (1.0 / mu_sun_air) if mu_sun_air > 0.0 else 0.0
    rho_I = atm_path_I * atm_refl_norm + TOA_wl_direct_I * rho_factor
    rho_Q = atm_path_Q * atm_refl_norm + TOA_wl_direct_Q * rho_factor
    rho_U = atm_path_U * atm_refl_norm + TOA_wl_direct_U * rho_factor

    # ---- dual-output sunglint --------------------------------------------
    rho_glint = [0.0, 0.0, 0.0]
    if wind_speed > 0.0:
        phi_v_g = (raa_deg - 180.0) * M_PI / 180.0
        mu_v_g = mu_v_air if vza_deg > 0.0 else 1.0
        rg = sf.direct_sunglint_rho(mu_v_g, phi_v_g, mu_sun_air, 0.0,
                                    wind_speed, sigma_type, n_water,
                                    tau_atm, q_convention)
        rho_glint = [rg[0], rg[1], rg[2]]

    return dict(
        rho_I=rho_I, rho_Q=rho_Q, rho_U=rho_U,
        rho_glint_direct_I=rho_glint[0], rho_glint_direct_Q=rho_glint[1],
        rho_glint_direct_U=-rho_glint[2],
        rho_I_glint=rho_I + rho_glint[0],
        rho_Q_glint=rho_Q + rho_glint[1],
        rho_U_glint=rho_U - rho_glint[2],
        I_TOA=I_TOA_total, Q_TOA=Q_TOA_total, U_TOA=U_TOA_total,
        rho_atm_path_I=atm_path_I * atm_refl_norm,
        rho_atm_path_Q=atm_path_Q * atm_refl_norm,
        rho_atm_path_U=atm_path_U * atm_refl_norm,
        rho_water_direct_I=TOA_wl_direct_I * rho_factor,
        rho_water_total_I=TOA_wl_direct_I * rho_factor,
        T_dir_dn=T_dir_beam, T_diff_dn_hemi=r1['T_diff_dn_hemi'],
        T_dir_up_view=T_atm_dir_up_view, T_diff_up_hemi=T_atm_diff_up_view,
        T_total_up_view=T_atm_up_view,
        Lu_0plus=Lu_total_above, Qu_0plus=Qu_total_above,
        Uu_0plus=Uu_total_above,
        Lu_0minus=Lu_total_below, Qu_0minus=Qu_total_below,
        Uu_0minus=Uu_total_below,
        Ed_0plus=Ed_above_air, Ed_0minus=Ed_total_below,
        Eu_0minus=w.Eu_0minus_water,
        Rrs0plus=(Lu_total_above / Ed_above_air) if Ed_above_air > 0 else 0.0,
        Rrs0plus_Q=(Qu_total_above / Ed_above_air) if Ed_above_air != 0 else 0.0,
        Rrs0plus_U=(Uu_total_above / Ed_above_air) if Ed_above_air != 0 else 0.0,
        rrs0minus=(Lu_total_below / Ed_total_below) if Ed_total_below > 0 else 0.0,
        rrs0minus_Q=(Qu_total_below / Ed_total_below) if Ed_total_below != 0 else 0.0,
        rrs0minus_U=(Uu_total_below / Ed_total_below) if Ed_total_below != 0 else 0.0,
        Kd=w.Kd_0minus, Ku=w.Ku_0minus,
        orders=w.max_orders_used, conv=w.all_converged,
        atm_path_I=atm_path_I, TOA_wl_direct_I=TOA_wl_direct_I,
        Ed_diff_water_from_atm=Ed_diff_water_from_atm,
        Ed_diff_BOA_air=Ed_diff_BOA_air,
        F_sun_BOA=F_sun_BOA, F_sun_water=F_sun_water, tau_atm=tau_atm,
    )


def solve_case_ocean_coupled_lut(
        sza_deg, vza_values, raa_values, wavelength_nm, n_water, wind_speed,
        a_tot, b_tot, bb_tot, phase_lut_path, mie, user_aod,
        pressure_hpa=1013.25, absorption=None, F_sun_TOA=1.0,
        n_mu_water=48, fourier_m_max=16, nt_atm=400, L_max=80,
        max_it_atm=100, tol_atm=1.0e-7, m_max_water=30,
        n_layers_water=300, max_it_water=160, tol_water=1.0e-7,
        sigma_type=1, q_convention=1, n_phi_quad=1024, n_phi_taw=128,
        constituent=None, aerosol_runtime=None):
    """Native coupled all-VZA/all-RAA LUT for one physical case.

    Atmosphere pass-1, in-water SOS, reverse coupling, and atmosphere pass-2
    are each evaluated once.  Requested view directions are zero-weight nodes
    in the shared angular rings; RAA dependence is Fourier reconstruction only.
    """
    wavelength_nm = require_wavelength(wavelength_nm, context='coupled ocean LUT')
    if wind_speed <= 0.0:
        raise NotImplementedError('flat-surface (wind=0) LUT not ported')
    if user_aod > 0.0 and aerosol_runtime is None and mie is None:
        raise ValueError('Aerosol runtime is required when AOD > 0')
    from .atmos import solve_aerosol_black_fresnel_ocean_value_lut
    from . import coupling as cpl
    from .kernel import gauss_legendre_pos
    import numpy as np

    vza_values = [float(v) for v in vza_values]
    raa_values = [float(r) for r in raa_values]
    nv, nr = len(vza_values), len(raa_values)
    mu_sun_air = cos(sza_deg * M_PI / 180.0)
    mu_v_air = np.cos(np.deg2rad(vza_values))

    # One exact rough-Fresnel Fourier cache is shared by atmosphere pass 1
    # and pass 2.  The expensive azimuth integration depends on the angular
    # ring and surface options, not on the injected water-leaving boundary.
    # The opt-out is retained only for byte-for-byte regression testing.
    surface_kernel_cache = None if os.environ.get(
        'OCRT_PY_AIR_SURFACE_CACHE_OFF') else {}

    # pass-1: one atmosphere field for every requested view and RAA
    r1 = solve_aerosol_black_fresnel_ocean_value_lut(
        sza_deg, vza_values, raa_values, wavelength_nm, wind_speed, mie,
        user_aod, L_max=L_max, pressure_hpa=pressure_hpa,
        n_mu_gl=n_mu_water, nt=nt_atm, m_max=fourier_m_max,
        max_iterations=max_it_atm, tolerance=tol_atm, n_water=n_water,
        sigma_type=sigma_type, q_convention=q_convention,
        n_phi_quad=n_phi_quad, surface='black_fresnel_ocean',
        absorption=absorption, aer_h_km=2.0, use_value_kernel=True,
        aerosol_runtime=aerosol_runtime,
        surface_kernel_cache=surface_kernel_cache)
    boa = r1['boa']
    mu_atm = boa['mu_quad_pos']; w_atm = boa['w_quad_pos']

    tau_abs_col = absorption.tau_total(wavelength_nm) if absorption is not None else 0.0
    tau_atm = r1['tau_R'] + r1['tau_a_eff'] + tau_abs_col
    F_sun_BOA = F_sun_TOA
    if mu_sun_air > 0.0 and tau_atm > 0.0:
        F_sun_BOA = F_sun_TOA * exp(-tau_atm / mu_sun_air)
    T_dir_beam = exp(-tau_atm / mu_sun_air) if (mu_sun_air > 0.0 and tau_atm > 0.0) else 1.0

    mu_wg, w_wg = gauss_legendre_pos(n_mu_water)
    cI, cQ, cU = cpl.couple_atm_to_water(
        boa['I_per_m'], boa['Q_per_m'], boa['U_per_m'], mu_atm,
        fourier_m_max, n_water, q_convention, np.asarray(mu_wg), w_atm,
        wind_speed, sigma_type, n_phi_kernel=n_phi_taw)
    dtI, dtQ, dtU = cI.copy(), cQ.copy(), cU.copy()
    cpl.apply_diffuse_top_basis_sign_inplace(dtI, dtQ, dtU, mode_axis=0)

    Ed_diff_water_from_atm = 0.0
    Ed_diff_BOA_air = 0.0
    for j in range(len(mu_atm)):
        Td = T_aw_coxmunk_direct(mu_atm[j], n_water, wind_speed, sigma_type)
        muL = w_atm[j] * mu_atm[j] * boa['I_per_m'][0, j]
        Ed_diff_water_from_atm += muL * Td
        Ed_diff_BOA_air += muL
    Ed_diff_water_from_atm *= 2.0 * M_PI * (F_sun_TOA / M_PI)
    Ed_diff_BOA_air *= 2.0 * M_PI * (F_sun_TOA / M_PI)

    mu_sun_water = snell_down(mu_sun_air, n_water)
    F_sun_water = F_sun_BOA
    if mu_sun_air > 0.0 and mu_sun_water > 0.0:
        T_aw = T_aw_coxmunk_direct(mu_sun_air, n_water, wind_speed, sigma_type)
        T_aw = min(1.0, max(0.0, T_aw))
        F_sun_water = F_sun_BOA * mu_sun_air * T_aw / mu_sun_water
        if F_sun_water > 0.0:
            dtsc = F_sun_TOA / F_sun_water
            dtI *= dtsc; dtQ *= dtsc; dtU *= dtsc

    # one in-water SOS field with every requested VZA inserted as a node
    w = sos_pure_fixed_bulk(
        sza_deg, vza_values[0], raa_values[0], wavelength_nm,
        n_water, F_sun_water, a_tot, b_tot, bb_tot,
        phase_lut_path, None, 0.0, wind_speed, sigma_type, q_convention,
        n_mu_water=n_mu_water, m_max_water=m_max_water,
        n_layers_water=n_layers_water, max_iterations=max_it_water,
        tolerance=tol_water, ext_top_I=dtI, ext_top_Q=dtQ, ext_top_U=dtU,
        ext_top_mu=np.asarray(mu_wg), ext_top_m_max=fourier_m_max,
        constituent=constituent, view_vza_values=vza_values,
        raa_values=raa_values, return_lut=True)

    mmax = min(fourier_m_max, w.m_max_filled)
    wlI, wlQ, wlU = cpl.couple_water_to_atm(
        w.I_up_per_m, w.Q_up_per_m, w.U_up_per_m, w.mu_water_pos,
        mmax, n_water, q_convention, mu_atm)
    # one pass-2 atmosphere field for all requested views/RAAs
    r2 = solve_aerosol_black_fresnel_ocean_value_lut(
        sza_deg, vza_values, raa_values, wavelength_nm, wind_speed, mie,
        user_aod, L_max=L_max, pressure_hpa=pressure_hpa,
        n_mu_gl=n_mu_water, nt=nt_atm, m_max=fourier_m_max,
        max_iterations=max_it_atm, tolerance=tol_atm, n_water=n_water,
        sigma_type=sigma_type, q_convention=q_convention,
        n_phi_quad=n_phi_quad, surface='black_fresnel_ocean',
        absorption=absorption, aer_h_km=2.0, use_value_kernel=True,
        bottom_source=(wlI, wlQ, wlU), aerosol_runtime=aerosol_runtime,
        surface_kernel_cache=surface_kernel_cache)

    Ed_0plus_air = F_sun_BOA * mu_sun_air if mu_sun_air > 0.0 else w.Ed_0plus_air
    Ed_total_below = w.Ed_0minus_water + Ed_diff_water_from_atm
    Ed_above_air = Ed_0plus_air + Ed_diff_BOA_air
    LmI = w.lut['I_0minus']; LmQ = w.lut['Q_0minus']; LmU = w.lut['U_0minus']
    LpI = w.lut['I_0plus'];  LpQ = w.lut['Q_0plus'];  LpU = w.lut['U_0plus']
    atmI, atmQ, atmU = r1['I_TOA'], r1['Q_TOA'], r1['U_TOA']
    wlTI, wlTQ, wlTU = r2['I_TOA'], r2['Q_TOA'], r2['U_TOA']

    rho_factor = (M_PI / (F_sun_TOA * mu_sun_air)) if (F_sun_TOA > 0.0 and mu_sun_air > 0.0) else 0.0
    atm_norm = (1.0 / mu_sun_air) if mu_sun_air > 0.0 else 0.0
    rho_I = atmI * atm_norm + wlTI * rho_factor
    rho_Q = atmQ * atm_norm + wlTQ * rho_factor
    rho_U = atmU * atm_norm + wlTU * rho_factor
    T_dir_up = np.exp(-tau_atm / np.maximum(mu_v_air, 1.0e-300))[:, None]
    T_diff_up = float(r1['T_diff_dn_hemi'])

    glI = np.zeros((nv, nr)); glQ = np.zeros((nv, nr)); glU = np.zeros((nv, nr))
    if wind_speed > 0.0:
        for iv, vza in enumerate(vza_values):
            for ir, raa in enumerate(raa_values):
                phi_v_g = (raa - 180.0) * M_PI / 180.0
                mu_v_g = mu_v_air[iv] if vza > 0.0 else 1.0
                rg = sf.direct_sunglint_rho(mu_v_g, phi_v_g, mu_sun_air, 0.0,
                                            wind_speed, sigma_type, n_water,
                                            tau_atm, q_convention)
                glI[iv, ir], glQ[iv, ir], glU[iv, ir] = rg[0], rg[1], -rg[2]

    return dict(
        rho_I=rho_I, rho_Q=rho_Q, rho_U=rho_U,
        rho_I_glint=rho_I+glI, rho_Q_glint=rho_Q+glQ, rho_U_glint=rho_U+glU,
        Lu_0plus=LpI, Qu_0plus=LpQ, Uu_0plus=LpU,
        Lu_0minus=LmI, Qu_0minus=LmQ, Uu_0minus=LmU,
        Ed_0plus=Ed_above_air, Ed_0minus=Ed_total_below,
        Rrs0plus=LpI/Ed_above_air if Ed_above_air else np.zeros_like(LpI),
        Rrs0plus_Q=LpQ/Ed_above_air if Ed_above_air else np.zeros_like(LpQ),
        Rrs0plus_U=LpU/Ed_above_air if Ed_above_air else np.zeros_like(LpU),
        rrs0minus=LmI/Ed_total_below if Ed_total_below else np.zeros_like(LmI),
        rrs0minus_Q=LmQ/Ed_total_below if Ed_total_below else np.zeros_like(LmQ),
        rrs0minus_U=LmU/Ed_total_below if Ed_total_below else np.zeros_like(LmU),
        orders=w.max_orders_used, conv=w.all_converged,
        T_dir_dn=T_dir_beam, T_diff_dn_hemi=r1['T_diff_dn_hemi'],
        T_dir_up_view=T_dir_up, T_diff_up_hemi=T_diff_up,
        Kd=w.Kd_0minus, Ku=w.Ku_0minus,
        diagnostics={'atmos_pass1_solves': 1, 'water_solves': 1,
                     'atmos_pass2_solves': 1, 'cell_solver_calls': 0,
                     'vza_count': nv, 'raa_count': nr,
                     'air_surface_cache_enabled':
                         surface_kernel_cache is not None,
                     'air_surface_cache_entries':
                         len(surface_kernel_cache)
                         if surface_kernel_cache is not None else 0})


def solve_case_ocean_tier0(sza_deg, vza_deg, raa_deg, wavelength_nm,
                           n_water, wind_speed, a_tot, b_tot, bb_tot,
                           phase_lut_path, inputs_dir,
                           pressure_hpa=0.0, aod=0.0,
                           sigma_type=1, q_convention=1,
                           F_sun_TOA=1.0, decouple_sunglint=True):
    if pressure_hpa > 0.0 or aod > 0.0:
        raise NotImplementedError('atm SOS branch not ported (pressure/aod>0)')

    mu_sun_air = cos(sza_deg * M_PI / 180.0)
    tau_atm = 0.0
    ab = Absorption(inputs_dir + '/afgl_atm', inputs_dir + '/xsec')
    tau_abs_col = ab.tau_total(wavelength_nm)
    tau_atm += tau_abs_col

    F_sun_BOA = F_sun_TOA
    if mu_sun_air > 0.0 and tau_atm > 0.0:
        F_sun_BOA = F_sun_TOA * exp(-tau_atm / mu_sun_air)

    # F_sun_water (option A, Cox-Munk direct T_aw for wind>0)
    mu_sun_water = snell_down(mu_sun_air, n_water)
    if wind_speed > 0.0:
        T_aw = T_aw_coxmunk_direct(mu_sun_air, n_water, wind_speed, sigma_type)
    else:
        M_R = sf.R_aa(mu_sun_air, n_water, q_convention)
        T_aw = 1.0 - M_R[0]
    T_aw = min(1.0, max(0.0, T_aw))
    F_sun_water = F_sun_BOA * mu_sun_air * T_aw / mu_sun_water

    w = sos_pure_fixed_bulk(sza_deg, vza_deg, raa_deg, wavelength_nm,
                            n_water, F_sun_water, a_tot, b_tot, bb_tot,
                            phase_lut_path, None, 0.0,
                            wind_speed, sigma_type, q_convention)

    # Ed_0plus_air 재정정
    Ed_above_air = F_sun_BOA * mu_sun_air

    # C assembly adds the (here zero) atm-diffuse superposition terms; the
    # additions also normalize any -0.0 from the sin reconstruction.
    Lu_total_below = w.I_0minus_view + 0.0
    Qu_total_below = w.Q_0minus_view + 0.0
    Uu_total_below = w.U_0minus_view + 0.0
    Ed_total_below = w.Ed_0minus_water + 0.0
    Lu_total_above = w.I_0plus_view + 0.0
    Qu_total_above = w.Q_0plus_view + 0.0
    Uu_total_above = w.U_0plus_view + 0.0

    mu_v_air = cos(vza_deg * M_PI / 180.0)
    # no atm branch -> T_atm_up = 1.0, atm_path = 0
    T_dir_up = 1.0
    T_up = 1.0
    atm_path_I = atm_path_Q = atm_path_U = 0.0

    rho_factor = (M_PI / (F_sun_TOA * mu_sun_air)) if (F_sun_TOA > 0.0 and mu_sun_air > 0.0) else 0.0
    atm_refl_norm = (1.0 / mu_sun_air) if mu_sun_air > 0.0 else 0.0
    rho_I = atm_path_I * atm_refl_norm + (Lu_total_above * T_up) * rho_factor
    rho_Q = atm_path_Q * atm_refl_norm + (Qu_total_above * T_dir_up) * rho_factor
    rho_U = atm_path_U * atm_refl_norm + (Uu_total_above * T_dir_up) * rho_factor

    # sunglint (always populated)
    rho_glint = [0.0, 0.0, 0.0]
    if wind_speed > 0.0:
        phi_v_g = (raa_deg - 180.0) * M_PI / 180.0
        mu_v_g = cos(vza_deg * M_PI / 180.0) if vza_deg > 0.0 else 1.0
        rg = sf.direct_sunglint_rho(mu_v_g, phi_v_g, mu_sun_air, 0.0,
                                    wind_speed, sigma_type, n_water,
                                    tau_atm, q_convention)
        rho_glint = [rg[0], rg[1], rg[2]]
    rho_I_glint = rho_I + rho_glint[0]
    rho_Q_glint = rho_Q + rho_glint[1]
    rho_U_glint = rho_U - rho_glint[2]

    out = dict(
        rho_I=rho_I, rho_Q=rho_Q, rho_U=rho_U,
        rho_I_glint=rho_I_glint, rho_Q_glint=rho_Q_glint,
        rho_U_glint=rho_U_glint,
        Lu_0plus=Lu_total_above, Qu_0plus=Qu_total_above,
        Uu_0plus=Uu_total_above,
        Lu_0minus=Lu_total_below, Qu_0minus=Qu_total_below,
        Uu_0minus=Uu_total_below,
        Ed_0plus=Ed_above_air, Ed_0minus=Ed_total_below,
        Rrs0plus=(Lu_total_above / Ed_above_air) if Ed_above_air > 0 else 0.0,
        Rrs0plus_Q=(Qu_total_above / Ed_above_air) if Ed_above_air != 0 else 0.0,
        Rrs0plus_U=(Uu_total_above / Ed_above_air) if Ed_above_air != 0 else 0.0,
        rrs0minus=(Lu_total_below / Ed_total_below) if Ed_total_below > 0 else 0.0,
        rrs0minus_Q=(Qu_total_below / Ed_total_below) if Ed_total_below != 0 else 0.0,
        rrs0minus_U=(Uu_total_below / Ed_total_below) if Ed_total_below != 0 else 0.0,
        Kd=w.Kd_0minus, Ku=w.Ku_0minus,
        a_total=w.a_total_used, b_total=w.b_total_used,
        bb_total=w.bb_total_used, omega_total=w.omega_used,
        orders=w.max_orders_used, conv=w.all_converged,
        tau_abs_col=tau_abs_col,
    )
    return out


def print_single_ocean(cs, r):
    """rt_io_print_single OCEAN branch."""
    wavelength_nm = require_wavelength(wavelength_nm, context='tier-0 ocean case')
    line = (
        "%.10e %.10e %.10e  "
        "rho_I_glint=%.10e rho_Q_glint=%.10e rho_U_glint=%.10e  "
        "Lu0plus=%.6e Qu0plus=%.6e Uu0plus=%.6e  "
        "Lu0minus=%.6e Qu0minus=%.6e Uu0minus=%.6e Ed0plus=%.6e Ed0minus=%.6e  "
        "Rrs0plus=%.6e Rrs0plus_Q=%.6e Rrs0plus_U=%.6e "
        "rrs0minus=%.6e rrs0minus_Q=%.6e rrs0minus_U=%.6e Kd=%.6f Ku=%.6f  "
        "a_w=%.8e b_w=%.8e bb_w=%.8e a_dom=%.8e "
        "a_pig=%.8e b_pig=%.8e bb_pig=%.8e "
        "a_min=%.8e b_min=%.8e bb_min=%.8e "
        "a_total=%.8e b_total=%.8e bb_total=%.8e omega_total=%.6f  "
        "sza=%.3f vza=%.3f raa=%.3f wl=%.3f orders=%d conv=%d"
    ) % (
        r['rho_I'], r['rho_Q'], r['rho_U'],
        r['rho_I_glint'], r['rho_Q_glint'], r['rho_U_glint'],
        r['Lu_0plus'], r['Qu_0plus'], r['Uu_0plus'],
        r['Lu_0minus'], r['Qu_0minus'], r['Uu_0minus'],
        r['Ed_0plus'], r['Ed_0minus'],
        r['Rrs0plus'], r['Rrs0plus_Q'], r['Rrs0plus_U'],
        r['rrs0minus'], r['rrs0minus_Q'], r['rrs0minus_U'],
        r['Kd'], r['Ku'],
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        r['a_total'], r['b_total'], r['bb_total'], r['omega_total'],
        cs['sza'], cs['vza'], cs['raa'], cs['wl'],
        r['orders'], r['conv'],
    )
    print(line)
