import os

from ocrt_py.aerosol import read_mie
from ocrt_py.constituent import OCRTConstituentModel
from ocrt_py.lut import run_ocean_aerosol_full_grid
from ocrt_solve import _constituent_greek


def _kwargs():
    data = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    mie = read_mie(os.path.join(data, 'C50.mie'))
    cm = OCRTConstituentModel(data, 'micro', 'red_clay')
    rr, constituent = _constituent_greek(
        cm, 490.0, 0.0, 3.5, 0.04, cdom_slope=0.014)
    return dict(
        sza_deg=30.0, wavelength_nm=490.0, n_water=1.34,
        wind_speed=3.0, a_tot=rr['a'], b_tot=rr['b'], bb_tot=rr['bb'],
        phase_lut_path=None, mie=mie, aod_ref=0.2, aod_ref_nm=865.0,
        vza_values=[0.0, 30.0], raa_values=[0.0, 90.0],
        pressure_hpa=1013.25, absorption=None, n_mu_water=4,
        fourier_m_max=2, nt_atm=8, L_max=20, max_it_atm=4,
        tol_atm=1.0e-5, m_max_water=2, n_layers_water=8,
        max_it_water=4, tol_water=1.0e-5, sigma_type=1,
        q_convention=1, n_phi_quad=32, n_phi_taw=16,
        constituent=constituent, progress=False)


def test_native_coupled_lut_matches_authoritative_cell_replay(monkeypatch):
    kw = _kwargs()
    monkeypatch.delenv('OCRT_PY_NATIVE_COUPLED_LUT_OFF', raising=False)
    native = run_ocean_aerosol_full_grid(**kw)
    monkeypatch.setenv('OCRT_PY_NATIVE_COUPLED_LUT_OFF', '1')
    legacy = run_ocean_aerosol_full_grid(**kw)
    assert native['rows'] == legacy['rows']
    d = native['diagnostics']
    assert d['atmos_pass1_solves'] == 1
    assert d['water_solves'] == 1
    assert d['atmos_pass2_solves'] == 1
    assert d['cell_solver_calls'] == 0
    assert legacy['diagnostics']['cell_solver_count'] == 4
