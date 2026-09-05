import os

import numpy as np

from ocrt_py import surface
from ocrt_py.aerosol import read_mie
from ocrt_py.constituent import OCRTConstituentModel
from ocrt_py.driver import solve_case_ocean_coupled_lut
from ocrt_solve import _constituent_greek


def _kwargs():
    root = os.path.dirname(os.path.dirname(__file__))
    data = os.path.join(root, 'data')
    mie = read_mie(os.path.join(data, 'C50.mie'))
    cm = OCRTConstituentModel(data, 'micro', 'red_clay')
    rr, constituent = _constituent_greek(
        cm, 490.0, 0.0, 3.5, 0.04, cdom_slope=0.014)
    return dict(
        sza_deg=30.0,
        vza_values=[0.0, 30.0],
        raa_values=[0.0, 90.0],
        wavelength_nm=490.0,
        n_water=1.34,
        wind_speed=3.0,
        a_tot=rr['a'],
        b_tot=rr['b'],
        bb_tot=rr['bb'],
        phase_lut_path=None,
        mie=mie,
        user_aod=0.2,
        pressure_hpa=1013.25,
        absorption=None,
        n_mu_water=4,
        fourier_m_max=2,
        nt_atm=8,
        L_max=20,
        max_it_atm=4,
        tol_atm=1.0e-5,
        m_max_water=2,
        n_layers_water=8,
        max_it_water=4,
        tol_water=1.0e-5,
        sigma_type=1,
        q_convention=1,
        n_phi_quad=32,
        n_phi_taw=16,
        constituent=constituent,
    )


def _public_arrays(result):
    keys = (
        'rho_I', 'rho_Q', 'rho_U',
        'rho_I_glint', 'rho_Q_glint', 'rho_U_glint',
        'Lu_0plus', 'Qu_0plus', 'Uu_0plus',
        'Lu_0minus', 'Qu_0minus', 'Uu_0minus',
        'Rrs0plus', 'Rrs0plus_Q', 'Rrs0plus_U',
        'rrs0minus', 'rrs0minus_Q', 'rrs0minus_U',
    )
    return tuple(np.asarray(result[key]) for key in keys)


def test_pass1_pass2_surface_fourier_cache_is_exact(monkeypatch):
    calls = {'count': 0}
    original = surface.fourier_kernel

    def counted(*args, **kwargs):
        calls['count'] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(surface, 'fourier_kernel', counted)
    monkeypatch.delenv('OCRT_PY_AIR_SURFACE_CACHE_OFF', raising=False)
    cached = solve_case_ocean_coupled_lut(**_kwargs())
    cached_calls = calls['count']

    calls['count'] = 0
    monkeypatch.setenv('OCRT_PY_AIR_SURFACE_CACHE_OFF', '1')
    uncached = solve_case_ocean_coupled_lut(**_kwargs())
    uncached_calls = calls['count']

    for lhs, rhs in zip(_public_arrays(cached), _public_arrays(uncached)):
        assert np.array_equal(lhs, rhs)

    m_count = _kwargs()['fourier_m_max'] + 1
    assert cached_calls == 2 * m_count
    assert uncached_calls == 3 * m_count
    assert cached['diagnostics']['air_surface_cache_enabled'] is True
    assert cached['diagnostics']['air_surface_cache_entries'] == 2 * m_count
    assert uncached['diagnostics']['air_surface_cache_enabled'] is False
    assert uncached['diagnostics']['air_surface_cache_entries'] == 0
