import os
import numpy as np

from ocrt_py.aerosol import read_mie
from ocrt_py.constituent import OCRTConstituentModel
from ocrt_py.lut import run_ocean_aerosol_full_grid
from ocrt_solve import _constituent_greek


def _phase_at(mie, wavelength_nm):
    order = np.argsort(mie.angles, kind='stable')
    theta = np.asarray(mie.angles)[order]
    wl_um = wavelength_nm / 1000.0
    def interp(block):
        b = np.asarray(block)[order, :]
        return np.array([np.interp(wl_um, mie.phase_wavelengths, b[i])
                         for i in range(len(theta))])
    return dict(theta_deg=theta, p11=interp(mie.P11),
                p12=interp(mie.P12), p33=interp(mie.P33), norm_n_mu=400)


def _kwargs():
    data = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    aer = read_mie(os.path.join(data, 'C50.mie'))
    cm = OCRTConstituentModel(data, 'micro', 'red_clay')
    rr, constituent = _constituent_greek(
        cm, 490.0, 0.0, 3.5, 0.04, L=6, nmg=20, cdom_slope=0.014)
    constituent = dict(constituent)
    constituent['value_phase'] = _phase_at(cm.min_mie, 490.0)
    return dict(
        sza_deg=30.0, wavelength_nm=490.0, n_water=1.34,
        wind_speed=3.0, a_tot=rr['a'], b_tot=rr['b'], bb_tot=rr['bb'],
        phase_lut_path=None, mie=aer, aod_ref=0.2, aod_ref_nm=865.0,
        vza_values=[30.0], raa_values=[0.0, 90.0],
        pressure_hpa=1013.25, absorption=None, n_mu_water=2,
        fourier_m_max=1, nt_atm=2, L_max=6, max_it_atm=2,
        tol_atm=1.0e-4, m_max_water=1, n_layers_water=2,
        max_it_water=2, tol_water=1.0e-4, sigma_type=1,
        q_convention=1, n_phi_quad=32, n_phi_taw=16,
        constituent=constituent, progress=False)


def test_polar_value_kernel_preserves_native_coupled_lut(monkeypatch):
    kw = _kwargs()
    monkeypatch.setenv('OCRT_WATER_VALUE_KERNEL_POL', '1')
    monkeypatch.setenv('OCRT_WATER_VALUE_NPHI', '64')
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
