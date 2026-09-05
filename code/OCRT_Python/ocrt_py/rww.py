"""Vectorized rough-Fresnel microfacet kernels (bit-level algebra identical to surface.py
scalars; used for speed in the m-mode Fourier tables)."""
import numpy as np
from math import pi, sqrt, exp
from .surface import slope_variance

M_PI = pi
_COS_ENTRY = np.array([1, 1, 0, 1, 1, 0, 0, 0, 1], dtype=bool)


def _sancer_lambda_vec(mu, sigma_len):
    lam = np.zeros_like(mu)
    mask = mu < 1.0 - 1e-10
    sin_t = np.sqrt(np.maximum(1e-20, 1.0 - mu * mu))
    nu = np.where(mask, (mu / np.where(sin_t > 0, sin_t, 1.0)) / sigma_len, 100.0)
    ok = mask & (nu <= 20.0)
    from scipy.special import erfc as _erfc
    lam_v = 0.5 * (np.exp(-nu * nu) / (nu * np.sqrt(M_PI)) - _erfc(nu))
    lam = np.where(ok, lam_v, 0.0)
    return lam


def _fresnel_R_mueller_vec(cos_omega, n1, n2, q_convention):
    """Vector Fresnel reflection Mueller entries (M00,M01,M22) with TIR
    identity.  Returns (M00, M01, M22, is_tir)."""
    mu_i = cos_omega
    sin_i_sq = np.maximum(0.0, 1.0 - mu_i * mu_i)
    sin_t_sq = (n1 / n2) ** 2 * sin_i_sq
    tir = sin_t_sq >= 1.0
    mu_t = np.sqrt(np.maximum(0.0, 1.0 - np.minimum(sin_t_sq, 1.0)))
    denom_s = n1 * mu_i + n2 * mu_t
    denom_p = n2 * mu_i + n1 * mu_t
    rs = np.where(tir, 0.0, (n1 * mu_i - n2 * mu_t) / np.where(denom_s != 0, denom_s, 1.0))
    rp = np.where(tir, 0.0, (n2 * mu_i - n1 * mu_t) / np.where(denom_p != 0, denom_p, 1.0))
    Rs = rs * rs
    Rp = rp * rp
    Qk = 0.5 * (Rp - Rs) if q_convention == 1 else 0.5 * (Rs - Rp)
    M00 = np.where(tir, 1.0, 0.5 * (Rs + Rp))
    M01 = np.where(tir, 0.0, Qk)
    M22 = np.where(tir, 1.0, rs * rp)
    return M00, M01, M22


def rww_coxmunk_allm(mu_o, mu_i, m_count, n_phi, ws, sigma_type, n_water,
                     q_convention):
    """All-m Fourier coefficients of surface_R_ww_coxmunk_trig over
    midpoint phi grid.  Returns K[m, jo, ji, 9]."""
    n_o = len(mu_o)
    n_i = len(mu_i)
    dphi = 2.0 * M_PI / n_phi
    phi = (np.arange(n_phi) + 0.5) * dphi
    cphi = np.cos(phi)
    sphi = np.sin(phi)

    sigma_sq = slope_variance(ws, sigma_type)
    sigma_len = sqrt(sigma_sq)

    MO, MI, PH = np.meshgrid(mu_o, mu_i, phi, indexing='ij')
    CP = np.cos(PH)
    SP = np.sin(PH)
    mu_in_actual = -MI
    s1 = np.sqrt(np.maximum(0.0, 1.0 - MO * MO))
    s2 = np.sqrt(np.maximum(0.0, 1.0 - mu_in_actual * mu_in_actual))
    cT = np.clip(MO * mu_in_actual + s1 * s2 * CP, -1.0, 1.0)
    cos_omega = np.sqrt(np.maximum(0.0, 0.5 * (1.0 - cT)))
    cos_beta = np.where(cos_omega > 1e-6, (MO + MI) / (2.0 * np.where(cos_omega > 0, cos_omega, 1.0)), 0.0)
    valid = (cos_omega > 1e-6) & (cos_beta > 1e-6)
    sT = np.sqrt(np.maximum(0.0, 1.0 - cT * cT))
    denom = np.maximum(1e-12, sT)
    s2_safe = np.maximum(1e-12, s2)
    s1_safe = np.maximum(1e-12, s1)
    ci1 = (MO - mu_in_actual * cT) / (denom * s2_safe)
    si1 = s1 * SP / denom
    n1n = np.sqrt(ci1 * ci1 + si1 * si1)
    small1 = n1n < 1e-12
    ci1 = np.where(small1, 1.0, ci1 / np.where(small1, 1.0, n1n))
    si1 = np.where(small1, 0.0, si1 / np.where(small1, 1.0, n1n))
    ci2 = (mu_in_actual - MO * cT) / (denom * s1_safe)
    si2 = s2 * SP / denom
    n2n = np.sqrt(ci2 * ci2 + si2 * si2)
    small2 = n2n < 1e-12
    ci2 = np.where(small2, 1.0, ci2 / np.where(small2, 1.0, n2n))
    si2 = np.where(small2, 0.0, si2 / np.where(small2, 1.0, n2n))

    M00, M01, M22 = _fresnel_R_mueller_vec(cos_omega, n_water, 1.0, q_convention)

    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / np.maximum(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * sigma_sq)) * np.exp(-tan_beta_sq / sigma_sq)
    lam_i = _sancer_lambda_vec(MI, sigma_len)
    lam_v = _sancer_lambda_vec(MO, sigma_len)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    brdf = (P_slope * S_bi) / (4.0 * MI * MO * cos_beta_sq * cos_beta_sq)
    brdf = np.where(valid, brdf, 0.0)

    # R = L2 · MF · L1 · brdf; MF = [[M00,M01,0],[M01,M00,0],[0,0,M22]]
    c2a = ci1 * ci1 - si1 * si1
    s2a = 2.0 * ci1 * si1
    c2b = ci2 * ci2 - si2 * si2
    s2b = 2.0 * ci2 * si2
    # MF@L1:
    # row0: [M00, M01*c2a, -M01*s2a]
    # row1: [M01, M00*c2a, -M00*s2a]
    # row2: [0,   M22*s2a,  M22*c2a]
    # L2@(MF@L1):
    R = np.empty(MO.shape + (9,))
    R[..., 0] = M00
    R[..., 1] = M01 * c2a
    R[..., 2] = -M01 * s2a
    R[..., 3] = c2b * M01
    R[..., 4] = c2b * M00 * c2a - s2b * M22 * s2a
    R[..., 5] = -c2b * M00 * s2a - s2b * M22 * c2a
    R[..., 6] = s2b * M01
    R[..., 7] = s2b * M00 * c2a + c2b * M22 * s2a
    R[..., 8] = -s2b * M00 * s2a + c2b * M22 * c2a
    R *= brdf[..., None]
    R = np.where(valid[..., None], R, 0.0)

    ms = np.arange(m_count)
    cm = np.cos(ms[:, None] * phi[None, :])      # (M, P)
    sm = np.sin(ms[:, None] * phi[None, :])
    W = np.where(_COS_ENTRY[None, :, None], cm[:, None, :], sm[:, None, :])  # (M,9,P)
    out = np.einsum('oipe,mep->moie', R, W)
    norm = np.where(ms == 0, 1.0 / (2.0 * M_PI), 1.0 / M_PI)
    out *= (dphi * norm)[:, None, None, None]
    return out


def T_aw_coxmunk_direct(mu_sun_air, n_water, wind_speed, sigma_type):
    """surface_T_aw_coxmunk_direct: slope-integrated air->water irradiance
    transmittance (trapezoid 101x101)."""
    if mu_sun_air <= 1e-9:
        return 0.0
    sig2 = slope_variance(wind_speed, sigma_type)
    sin_s = sqrt(max(0.0, 1.0 - mu_sun_air * mu_sun_air))
    if sig2 <= 1e-10:
        cti0 = mu_sun_air
        stt0 = sin_s / n_water
        ctt0 = sqrt(max(0.0, 1.0 - stt0 * stt0))
        rs0 = (cti0 - n_water * ctt0) / (cti0 + n_water * ctt0)
        rp0 = (n_water * cti0 - ctt0) / (n_water * cti0 + ctt0)
        return 1.0 - 0.5 * (rs0 * rs0 + rp0 * rp0)
    sig = sqrt(sig2)
    L = 6.0 * sig
    N = 101
    dz = 2.0 * L / (N - 1)
    z = -L + np.arange(N) * dz
    wt = np.ones(N)
    wt[0] = wt[-1] = 0.5
    ZX, ZY = np.meshgrid(z, z, indexing='ij')
    WX, WY = np.meshgrid(wt, wt, indexing='ij')
    proj = mu_sun_air + ZX * sin_s
    r2 = ZX * ZX + ZY * ZY
    cti = np.where(proj > 0.0, proj / np.sqrt(1.0 + r2), 0.0)
    sti2 = np.maximum(0.0, 1.0 - cti * cti)
    stt = np.sqrt(sti2) / n_water
    ctt = np.sqrt(np.maximum(0.0, 1.0 - stt * stt))
    rs = (cti - n_water * ctt) / (cti + n_water * ctt)
    rp = (n_water * cti - ctt) / (n_water * cti + ctt)
    TF = 1.0 - 0.5 * (rs * rs + rp * rp)
    p = (1.0 / (M_PI * sig2)) * np.exp(-r2 / sig2)
    integrand = np.where(proj > 0.0, WX * WY * p * proj * TF, 0.0)
    acc = integrand.sum() * dz * dz
    T = acc / mu_sun_air
    return min(1.0, max(0.0, T))


def R_ww_coxmunk_direct(mu_up_water, n_water, wind_speed, sigma_type):
    """surface_R_ww_coxmunk_direct (first-encounter integral)."""
    if mu_up_water <= 1e-9:
        return 1.0
    s2v = slope_variance(wind_speed, sigma_type)
    sin_u = sqrt(max(0.0, 1.0 - mu_up_water * mu_up_water))

    def rwa_flat(cti):
        sti2 = np.maximum(0.0, 1.0 - cti * cti)
        sair = n_water * np.sqrt(sti2)
        tir = sair >= 1.0
        cair = np.sqrt(np.maximum(0.0, 1.0 - np.minimum(sair, 1.0) ** 2))
        rs = (n_water * cti - cair) / (n_water * cti + cair)
        rp = (cti - n_water * cair) / (cti + n_water * cair)
        return np.where(tir, 1.0, 0.5 * (rs * rs + rp * rp))

    if s2v <= 1e-10:
        return float(rwa_flat(np.array([mu_up_water]))[0])
    sig = sqrt(s2v)
    L = 6.0 * sig
    N = 101
    dz = 2.0 * L / (N - 1)
    z = -L + np.arange(N) * dz
    wt = np.ones(N)
    wt[0] = wt[-1] = 0.5
    ZX, ZY = np.meshgrid(z, z, indexing='ij')
    WX, WY = np.meshgrid(wt, wt, indexing='ij')
    proj = mu_up_water + ZX * sin_u
    r2 = ZX * ZX + ZY * ZY
    cti = np.where(proj > 0.0, proj / np.sqrt(1.0 + r2), 0.0)
    Rwa = rwa_flat(cti)
    p = (1.0 / (M_PI * s2v)) * np.exp(-r2 / s2v)
    integ = np.where(proj > 0.0, WX * WY * p * proj * Rwa, 0.0)
    acc = integ.sum() * dz * dz
    return min(1.0, max(0.0, acc / mu_up_water))
