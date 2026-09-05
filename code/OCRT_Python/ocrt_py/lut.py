"""LUT full-grid driver (Python analogue of C --output-full-grid), batched
and GPU-capable.  Emits the same CSV row format as ocean_rrs_full_grid."""
import numpy as np
import sys
from math import pi, cos, sin, exp, log, floor

from . import surface as sf
from .absorption import Absorption
from .rww import T_aw_coxmunk_direct
from .water import snell_down, aw_interp_on_unsorted
from . import lutbatch as lb

M_PI = pi

CSV_HEADER = (
    "case_kind,sza_deg,wavelength_nm,wind_speed_ms,n_mu,mu_index,mu_view,"
    "vza_deg,raa_deg,"
    "Ed0plus_air,Lu0plus_I,Lu0plus_Q,Lu0plus_U,Rrs_I,Rrs_Q,Rrs_U,"
    "Ed0minus_water,Lu0minus_I,Lu0minus_Q,Lu0minus_U,rrs_I,rrs_Q,rrs_U,"
    "TOA_rho_I,TOA_rho_Q,TOA_rho_U,a_total,b_total,bb_total,omega_water,"
    "water_orders,water_converged")


def run_lut_grid(case_specs, sza_deg, n_water, inputs_dir,
                 vza_step_deg=2.5, vza_max_deg=85.0, raa_step_deg=2.5,
                 n_mu_water=96, m_max_water=30, max_iterations=10000,
                 tolerance=1.0e-7, sigma_type=1, q_convention=1,
                 F_sun_TOA=1.0, n_mu_header=24, out_csv=None,
                 phase_nphi=720, progress=False):
    """case_specs: list of dicts with keys a, b, bb, wind, wl, phase_lut
    (+optional phase_case_id, phase_wl).  All cases share sza & n_water
    (batch precondition: identical angular ring).  Atmosphere-free route
    (pressure 0, aod 0) — the only route currently ported."""
    n_vza = int(floor(vza_max_deg / vza_step_deg + 1e-9)) + 1
    vza_list = [iv * vza_step_deg for iv in range(n_vza)]
    n_raa = max(1, int(floor(360.0 / raa_step_deg + 0.5)))
    raa_list = [ir * raa_step_deg for ir in range(n_raa)]

    ab = Absorption(inputs_dir + '/afgl_atm', inputs_dir + '/xsec')
    mu_sun_air = cos(sza_deg * M_PI / 180.0)
    mu_sun_water = snell_down(mu_sun_air, n_water)

    cases = []
    F_list = []
    meta = []
    n_spec = len(case_specs)
    for i_spec, spec in enumerate(case_specs):
        if progress:
            sys.stderr.write('\r  [prep] %d/%d cases   ' % (i_spec + 1, n_spec))
            sys.stderr.flush()
        c = lb.prepare_case(spec['a'], spec['b'], spec['bb'], spec['wind'],
                            spec['wl'], spec['phase_lut'],
                            spec.get('phase_case_id'),
                            spec.get('phase_wl', 0.0))
        tau_abs = ab.tau_total(spec['wl'])
        F_BOA = F_sun_TOA * exp(-tau_abs / mu_sun_air) if tau_abs > 0.0 else F_sun_TOA
        if spec['wind'] > 0.0:
            T_aw = T_aw_coxmunk_direct(mu_sun_air, n_water, spec['wind'],
                                       sigma_type)
        else:
            T_aw = 1.0 - sf.R_aa(mu_sun_air, n_water, q_convention)[0]
        T_aw = min(1.0, max(0.0, T_aw))
        c.F_sun = F_BOA * mu_sun_air * T_aw / mu_sun_water
        cases.append(c)
        F_list.append(c.F_sun)
        meta.append(dict(F_BOA=F_BOA, tau_abs=tau_abs))

    if progress:
        sys.stderr.write('\r  [prep] %d/%d cases  done\n' % (n_spec, n_spec))
        sys.stderr.flush()

    res = lb.solve_batch(cases, sza_deg, vza_list, n_water,
                         n_mu_water=n_mu_water, m_max_water=m_max_water,
                         max_iterations=max_iterations, tolerance=tolerance,
                         sigma_type=sigma_type, q_convention=q_convention,
                         phase_nphi=phase_nphi, F_sun_list=F_list,
                         first_view_idx=0, progress=progress)

    rows = []
    mu_pos = res.mu_pos
    w_pos = res.w_pos
    n_mu = res.n_mu
    Mst = m_max_water + 1
    n_case = len(cases)
    for b, c in enumerate(cases):
        if progress:
            sys.stderr.write('\r  [assemble] %d/%d cases   ' % (b + 1, n_case))
            sys.stderr.flush()
        F_sun = c.F_sun
        f_scale = F_sun / M_PI
        Ed_direct = F_sun * mu_sun_water
        Ed_diffuse = float((res.I_neg0[b] * mu_pos * w_pos).sum()) * 2.0 * M_PI
        Ed_total = Ed_direct + Ed_diffuse * f_scale
        Ed_above = meta[b]['F_BOA'] * mu_sun_air
        # Kd/Ku (view-independent)
        Ed_l1 = Eu_l1 = 0.0
        for jp in range(1, n_mu + 1):
            Ed_l1 += res.I_lvl1[b, n_mu - jp] * mu_pos[jp - 1] * w_pos[jp - 1]
            Eu_l1 += res.I_lvl1[b, n_mu + jp] * mu_pos[jp - 1] * w_pos[jp - 1]
        Ed_l1 *= 2.0 * M_PI * f_scale
        Eu_l1 *= 2.0 * M_PI * f_scale
        tau_l1 = c.tau_max / c.nt
        Ed_dir_l1 = F_sun * mu_sun_water * exp(-tau_l1 / mu_sun_water)
        Ed_tot_l1 = Ed_dir_l1 + Ed_l1
        Eu0 = float((res.I_node[b, 0] * mu_pos * w_pos).sum()) * 2.0 * M_PI * f_scale
        z_l1 = tau_l1 / c.ext
        Kd = -log(Ed_tot_l1 / Ed_total) / max(z_l1, 1e-12) \
            if (Ed_tot_l1 > 0.0 and Ed_total > 0.0) else 0.0
        Ku = -log(Eu_l1 / Eu0) / max(z_l1, 1e-12) \
            if (Eu_l1 > 0.0 and Eu0 > 0.0) else 0.0

        rho_factor = M_PI / (F_sun_TOA * mu_sun_air)
        m_avail = res.m_exit[b] if res.exited[b] else (Mst - 1)
        for iv, vza in enumerate(vza_list):
            jn = res.node_of_view[iv]
            mu_va = cos(vza * M_PI / 180.0)
            mu_wv = snell_down(mu_va, n_water)
            M_T_wa = sf.T_wa(mu_wv, n_water, q_convention)
            # --- C grid-cache replay parity (rt_water_rt.c commit #14/#16):
            # every row RE-RUNS the m loop over the cached modes with ITS OWN
            # view amplitude, so the early-exit point is view-dependent.
            # Ordering per mode: (7e) view extraction -> exit check/break ->
            # (7e-bis) node store.  Hence for a row that breaks at mode mb:
            #   view reconstruction uses modes 0..mb (inclusive),
            #   the water->air coupling (I_m_node) uses modes 0..mb-1.
            # A row that never breaks uses all cached modes for both.
            amp_max = 0.0
            streak = 0
            mb = None                       # triggering mode for this row
            for m in range(m_avail + 1):
                amp = abs(res.I_node[b, m, jn])   # Q/U view values are 0 here
                if amp > amp_max:
                    amp_max = amp
                if m >= 3 and amp_max > 0.0:
                    if amp < tolerance * amp_max:
                        streak += 1
                    else:
                        streak = 0
                    if streak >= 2:
                        mb = m
                        break
            m_view_hi = mb if mb is not None else m_avail    # inclusive
            m_cpl_hi = (mb - 1) if mb is not None else m_avail
            Im_view = np.zeros(Mst)
            Im_view[:m_view_hi + 1] = res.I_node[b, :m_view_hi + 1, jn]
            # per-m above-water via analytic flat coupling at the node
            mu_list = list(mu_pos)
            Ia_m = np.zeros(Mst)
            Qa_m = np.zeros(Mst)
            Ua_m = np.zeros(Mst)
            for m in range(m_cpl_hi + 1):
                I_w = aw_interp_on_unsorted(mu_list, list(res.I_node[b, m]),
                                            mu_wv)
                # Q_w = U_w = 0 in this route; full Mueller rows on I_w
                Ia_m[m] = M_T_wa[0] * I_w
                Qa_m[m] = M_T_wa[3] * I_w
                Ua_m[m] = M_T_wa[6] * I_w
            for ir, raa in enumerate(raa_list):
                dphi = (180.0 - raa) * M_PI / 180.0
                base = dphi + M_PI
                Lu_below = Im_view[0]
                Lu_above = Ia_m[0]
                Qu_above = Qa_m[0]
                Uu_above = 0.0
                for m in range(1, Mst):
                    cf = 2.0 * cos(m * base)
                    sf_ = 2.0 * sin(m * base)
                    Lu_below += cf * Im_view[m]
                    Lu_above += cf * Ia_m[m]
                    Qu_above += cf * Qa_m[m]
                    Uu_above += -sf_ * Ua_m[m]
                Lu_below = Lu_below * f_scale + 0.0
                Lu_above = Lu_above * f_scale + 0.0
                Qu_above = Qu_above * f_scale + 0.0
                Uu_above = Uu_above * f_scale + 0.0
                Rrs = Lu_above / Ed_above if Ed_above > 0.0 else 0.0
                Rq = Qu_above / Ed_above if Ed_above > 0.0 else 0.0
                Ru = Uu_above / Ed_above if Ed_above > 0.0 else 0.0
                rrs = Lu_below / Ed_total if Ed_total > 0.0 else 0.0
                rho_I = Lu_above * rho_factor
                rho_Q = Qu_above * rho_factor
                rho_U = Uu_above * rho_factor
                rows.append((
                    "ocean_rrs_grid,%.10g,%.10g,%.10g,%d,%d,%.12g,%.10g,%.10g,"
                    "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,"
                    "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,"
                    "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%d,%d"
                ) % (
                    sza_deg, c.wl_nm, c.wind, n_mu_header, iv, mu_va, vza, raa,
                    Ed_above, Lu_above, Qu_above, Uu_above,
                    Rrs, Rq, Ru,
                    Ed_total, Lu_below, 0.0, 0.0,
                    rrs, 0.0, 0.0,
                    rho_I, rho_Q, rho_U,
                    c.a_tot, c.b_tot, c.bb_tot, c.omega,
                    res.n_orders[b], res.all_conv[b]))

    if progress:
        sys.stderr.write('\r  [assemble] %d/%d cases  done\n' % (n_case, n_case))
        sys.stderr.flush()

    text = CSV_HEADER + "\n" + "\n".join(rows) + "\n"
    if out_csv:
        with open(out_csv, 'w', encoding='utf-8', newline='') as f:
            f.write(text)
    return text


AEROSOL_FULLGRID_FIELDS = [
    'case_kind', 'sza_deg', 'wavelength_nm', 'wind_speed_ms',
    'mu_index', 'mu_view', 'vza_deg', 'raa_deg',
    'Ed0plus_air', 'Lu0plus_I', 'Lu0plus_Q', 'Lu0plus_U',
    'Rrs_I', 'Rrs_Q', 'Rrs_U',
    'Ed0minus_water', 'Lu0minus_I', 'Lu0minus_Q', 'Lu0minus_U',
    'rrs_I', 'rrs_Q', 'rrs_U',
    'TOA_rho_I', 'TOA_rho_Q', 'TOA_rho_U',
    'a_total', 'b_total', 'bb_total', 'omega_water',
    'water_orders', 'water_converged',
    'aod_ref', 'aod_ref_nm', 'aod_target', 'tau_a_eff', 'ssa_a_eff',
]


def run_ocean_aerosol_full_grid(
        *, sza_deg, wavelength_nm, n_water, wind_speed,
        a_tot, b_tot, bb_tot, phase_lut_path, mie,
        aod_ref, aod_ref_nm, vza_values, raa_values,
        pressure_hpa=1013.25, absorption=None, F_sun_TOA=1.0,
        n_mu_water=48, fourier_m_max=16, nt_atm=400, L_max=80,
        max_it_atm=100, tol_atm=1.0e-7,
        m_max_water=30, n_layers_water=300,
        max_it_water=160, tol_water=1.0e-7,
        sigma_type=1, q_convention=1, n_phi_quad=1024,
        n_phi_taw=128, constituent=None, out_csv=None,
        progress=False):
    """Atmosphere-ocean full-grid path with one explicit aerosol runtime.

    This is the Python counterpart of the C full-grid aerosol-object fix.  The
    Mie file is assumed to have been read once by the caller; the immutable
    ``AerosolRuntime`` is prepared once here and passed explicitly to every
    VZA/RAA cell and to both atmospheric passes inside each coupled solve.
    """
    import csv
    from .aerosol import prepare_aerosol_runtime
    from .driver import solve_case_ocean_coupled, solve_case_ocean_coupled_lut

    if mie is None:
        if float(aod_ref) > 0.0:
            raise ValueError(
                'Aerosol runtime is required when AOD > 0 in ocean full-grid mode.')
        raise ValueError('MieData is required for deterministic AOD=0 control runs')

    runtime = prepare_aerosol_runtime(
        mie, wavelength_nm, aod_ref, aod_ref_nm=aod_ref_nm, L_max=L_max)
    vza_values = [float(v) for v in vza_values]
    raa_values = [float(r) for r in raa_values]
    rows = []
    n_cells = len(vza_values) * len(raa_values)
    omega = b_tot / (a_tot + b_tot) if (a_tot + b_tot) > 0.0 else 0.0
    native = not bool(__import__('os').environ.get('OCRT_PY_NATIVE_COUPLED_LUT_OFF'))
    if native:
        result = solve_case_ocean_coupled_lut(
            sza_deg, vza_values, raa_values, wavelength_nm,
            n_water=n_water, wind_speed=wind_speed,
            a_tot=a_tot, b_tot=b_tot, bb_tot=bb_tot,
            phase_lut_path=phase_lut_path,
            mie=None, user_aod=runtime.aod_target,
            pressure_hpa=pressure_hpa, absorption=absorption,
            F_sun_TOA=F_sun_TOA, n_mu_water=n_mu_water,
            fourier_m_max=fourier_m_max, nt_atm=nt_atm,
            L_max=L_max, max_it_atm=max_it_atm, tol_atm=tol_atm,
            m_max_water=m_max_water, n_layers_water=n_layers_water,
            max_it_water=max_it_water, tol_water=tol_water,
            sigma_type=sigma_type, q_convention=q_convention,
            n_phi_quad=n_phi_quad, n_phi_taw=n_phi_taw,
            constituent=constituent, aerosol_runtime=runtime)
        for iv, vza in enumerate(vza_values):
            for ir, raa in enumerate(raa_values):
                rows.append(dict(
                    case_kind='ocean_rrs_grid', sza_deg=float(sza_deg),
                    wavelength_nm=float(wavelength_nm), wind_speed_ms=float(wind_speed),
                    mu_index=iv, mu_view=cos(vza * M_PI / 180.0),
                    vza_deg=vza, raa_deg=raa,
                    Ed0plus_air=float(result['Ed_0plus']),
                    Lu0plus_I=float(result['Lu_0plus'][iv, ir]),
                    Lu0plus_Q=float(result['Qu_0plus'][iv, ir]),
                    Lu0plus_U=float(result['Uu_0plus'][iv, ir]),
                    Rrs_I=float(result['Rrs0plus'][iv, ir]),
                    Rrs_Q=float(result['Rrs0plus_Q'][iv, ir]),
                    Rrs_U=float(result['Rrs0plus_U'][iv, ir]),
                    Ed0minus_water=float(result['Ed_0minus']),
                    Lu0minus_I=float(result['Lu_0minus'][iv, ir]),
                    Lu0minus_Q=float(result['Qu_0minus'][iv, ir]),
                    Lu0minus_U=float(result['Uu_0minus'][iv, ir]),
                    rrs_I=float(result['rrs0minus'][iv, ir]),
                    rrs_Q=float(result['rrs0minus_Q'][iv, ir]),
                    rrs_U=float(result['rrs0minus_U'][iv, ir]),
                    TOA_rho_I=float(result['rho_I'][iv, ir]),
                    TOA_rho_Q=float(result['rho_Q'][iv, ir]),
                    TOA_rho_U=float(result['rho_U'][iv, ir]),
                    a_total=float(a_tot), b_total=float(b_tot), bb_total=float(bb_tot),
                    omega_water=omega, water_orders=result['orders'],
                    water_converged=result['conv'],
                    aod_ref=runtime.aod_ref, aod_ref_nm=runtime.aod_ref_nm,
                    aod_target=runtime.aod_target, tau_a_eff=runtime.tau_a_eff,
                    ssa_a_eff=runtime.ssa_a_eff))
        diagnostics = dict(result['diagnostics'])
        diagnostics.update(mie_read_count=1, aerosol_prepare_count=1,
                           runtime_object_id=id(runtime), native_lut=True)
    else:
        cell_no = 0
        for iv, vza in enumerate(vza_values):
            for raa in raa_values:
                cell_no += 1
                if progress:
                    sys.stderr.write('\r  [py full-grid legacy] %d/%d' % (cell_no, n_cells))
                    sys.stderr.flush()
                result = solve_case_ocean_coupled(
                    sza_deg, vza, raa, wavelength_nm,
                    n_water=n_water, wind_speed=wind_speed,
                    a_tot=a_tot, b_tot=b_tot, bb_tot=bb_tot,
                    phase_lut_path=phase_lut_path,
                    mie=None, user_aod=runtime.aod_target,
                    pressure_hpa=pressure_hpa, absorption=absorption,
                    F_sun_TOA=F_sun_TOA, n_mu_water=n_mu_water,
                    fourier_m_max=fourier_m_max, nt_atm=nt_atm,
                    L_max=L_max, max_it_atm=max_it_atm, tol_atm=tol_atm,
                    m_max_water=m_max_water, n_layers_water=n_layers_water,
                    max_it_water=max_it_water, tol_water=tol_water,
                    sigma_type=sigma_type, q_convention=q_convention,
                    n_phi_quad=n_phi_quad, n_phi_taw=n_phi_taw,
                    constituent=constituent, aerosol_runtime=runtime)
                rows.append(dict(
                    case_kind='ocean_rrs_grid', sza_deg=float(sza_deg),
                    wavelength_nm=float(wavelength_nm), wind_speed_ms=float(wind_speed),
                    mu_index=iv, mu_view=cos(vza * M_PI / 180.0),
                    vza_deg=vza, raa_deg=raa,
                    Ed0plus_air=result['Ed_0plus'], Lu0plus_I=result['Lu_0plus'],
                    Lu0plus_Q=result['Qu_0plus'], Lu0plus_U=result['Uu_0plus'],
                    Rrs_I=result['Rrs0plus'], Rrs_Q=result['Rrs0plus_Q'],
                    Rrs_U=result['Rrs0plus_U'], Ed0minus_water=result['Ed_0minus'],
                    Lu0minus_I=result['Lu_0minus'], Lu0minus_Q=result['Qu_0minus'],
                    Lu0minus_U=result['Uu_0minus'], rrs_I=result['rrs0minus'],
                    rrs_Q=result['rrs0minus_Q'], rrs_U=result['rrs0minus_U'],
                    TOA_rho_I=result['rho_I'], TOA_rho_Q=result['rho_Q'],
                    TOA_rho_U=result['rho_U'], a_total=float(a_tot),
                    b_total=float(b_tot), bb_total=float(bb_tot),
                    omega_water=omega, water_orders=result['orders'],
                    water_converged=result['conv'], aod_ref=runtime.aod_ref,
                    aod_ref_nm=runtime.aod_ref_nm, aod_target=runtime.aod_target,
                    tau_a_eff=runtime.tau_a_eff, ssa_a_eff=runtime.ssa_a_eff))
        diagnostics = {'mie_read_count': 1, 'aerosol_prepare_count': 1,
                       'cell_solver_count': n_cells, 'runtime_object_id': id(runtime),
                       'native_lut': False}
    if progress:
        sys.stderr.write('\n')
        sys.stderr.flush()
    if out_csv:
        with open(out_csv, 'w', newline='', encoding='utf-8') as fp:
            writer = csv.DictWriter(fp, fieldnames=AEROSOL_FULLGRID_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
    return {
        'rows': rows,
        'aerosol_runtime': runtime,
        'diagnostics': diagnostics,
    }
