import numpy as np
from math import cos, pi

from ocrt_py import value_phase as vp
from ocrt_py.aerosol import compute_vector_legendre_gauss
from ocrt_py.kernel import (
    gauss_legendre_pos, LegendreWorkspace, legendre_compute_pol,
    kernel_phase_fourier, kernel_phase_fourier_pol,
    kernel_phase_fourier_aerosol_full,
)
from ocrt_py.sos import build_inwater_atm


def _rayleigh_like():
    theta = np.linspace(0.0, 180.0, 361)
    mu = np.cos(np.deg2rad(theta))
    p11 = 0.75 * (1.0 + mu * mu)
    p12 = -0.75 * (1.0 - mu * mu)
    p33 = 1.5 * mu
    return theta, p11, p12, p33


def test_direct_polar_value_kernel_matches_coefficient_path():
    theta, p11, p12, p33 = _rayleigh_like()
    l_max, n_mu, m_count = 10, 6, 5
    beta, gamma, alpha, zeta = compute_vector_legendre_gauss(
        p11, p12, p33, theta, l_max, 400)
    gl_mu, gl_w = gauss_legendre_pos(n_mu)
    atm = build_inwater_atm(
        1, n_mu, 0.1, 0.0, 1.0, cos(30.0*pi/180.0),
        n_mu, l_max, beta, gamma, alpha, zeta, gl_mu, gl_w)
    phase = vp.ValuePhaseInterp(theta, p11, p12, p33, 400)
    direct = vp.vector_fourier_allm(
        np.asarray(atm.rm), n_mu, m_count, phase, 4096)
    ws = LegendreWorkspace(n_mu, l_max)
    max_abs = 0.0
    for m in range(m_count):
        legendre_compute_pol(ws, np.asarray(atm.rm), m)
        kernel_phase_fourier(ws, m, beta)
        kernel_phase_fourier_pol(ws, m, gamma)
        kernel_phase_fourier_aerosol_full(ws, m, alpha, zeta)
        coeff = (ws.pfm, ws.gr, ws.gt, ws.arr, ws.art, ws.att)
        for c in range(6):
            max_abs = max(max_abs, float(np.max(np.abs(coeff[c] - direct[c][m]))))
    assert max_abs <= 8.0e-12


def test_scalar_value_kernel_uses_same_spline_and_gauss_normalization():
    theta, p11, p12, p33 = _rayleigh_like()
    l_max, n_mu, m_count = 10, 6, 5
    beta, gamma, alpha, zeta = compute_vector_legendre_gauss(
        p11, p12, p33, theta, l_max, 400)
    gl_mu, gl_w = gauss_legendre_pos(n_mu)
    atm = build_inwater_atm(
        1, n_mu, 0.1, 0.0, 1.0, cos(30.0*pi/180.0),
        n_mu, l_max, beta, gamma, alpha, zeta, gl_mu, gl_w)
    phase = vp.ValuePhaseInterp(theta, p11, norm_n_mu=400)
    direct = vp.scalar_fourier_allm(
        np.asarray(atm.rm), n_mu, m_count, phase, 4096)
    ws = LegendreWorkspace(n_mu, l_max)
    max_abs = 0.0
    for m in range(m_count):
        legendre_compute_pol(ws, np.asarray(atm.rm), m)
        kernel_phase_fourier(ws, m, beta)
        max_abs = max(max_abs, float(np.max(np.abs(ws.pfm - direct[m]))))
    assert max_abs <= 8.0e-12
