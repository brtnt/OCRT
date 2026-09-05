"""Batched (case-axis) water-water rough-Fresnel microfacet reflection Fourier kernel.

rww_coxmunk_allm_batch stacks the per-case ring (mu_o, mu_i) and wind over a
leading B axis so all cases' surface kernels compute in one backend pass,
eliminating the per-case host loop.  Bit-identical to rww.rww_coxmunk_allm
applied case by case.
"""
import numpy as np
from math import pi as M_PI
from .backend import xp, erfc

_COS_ENTRY = np.array([1, 1, 0, 1, 1, 0, 0, 0, 1], dtype=bool)


def _slope_variance_batch(ws, sigma_type):
    ws = xp.asarray(ws, dtype=float)
    if sigma_type == 1:
        return 0.003 + 0.00512 * xp.maximum(0.01, ws)
    sig = 0.0731 * xp.sqrt(xp.maximum(0.0, ws))
    return xp.maximum(sig * sig, 1e-10)


def _sancer_lambda_batch(mu, sigma_len):
    """Sancer shadowing lambda (broadcasting).  mu any shape, sigma_len
    broadcastable.  Uses backend erfc (scipy on numpy, cupyx on GPU)."""
    mask = mu < 1.0 - 1e-10
    sin_t = xp.sqrt(xp.maximum(1e-20, 1.0 - mu * mu))
    nu = xp.where(mask, (mu / xp.where(sin_t > 0, sin_t, 1.0)) / sigma_len, 100.0)
    ok = mask & (nu <= 20.0)
    lam_v = 0.5 * (xp.exp(-nu * nu) / (nu * xp.sqrt(M_PI)) - erfc(nu))
    return xp.where(ok, lam_v, 0.0)


def _fresnel_R_mueller_batch(cos_omega, n1, n2, q_convention):
    mu_i = cos_omega
    sin_i_sq = xp.maximum(0.0, 1.0 - mu_i * mu_i)
    sin_t_sq = (n1 / n2) ** 2 * sin_i_sq
    tir = sin_t_sq >= 1.0
    mu_t = xp.sqrt(xp.maximum(0.0, 1.0 - xp.minimum(sin_t_sq, 1.0)))
    denom_s = n1 * mu_i + n2 * mu_t
    denom_p = n2 * mu_i + n1 * mu_t
    rs = xp.where(tir, 0.0, (n1 * mu_i - n2 * mu_t) / xp.where(denom_s != 0, denom_s, 1.0))
    rp = xp.where(tir, 0.0, (n2 * mu_i - n1 * mu_t) / xp.where(denom_p != 0, denom_p, 1.0))
    Rs = rs * rs
    Rp = rp * rp
    Qk = 0.5 * (Rp - Rs) if q_convention == 1 else 0.5 * (Rs - Rp)
    M00 = xp.where(tir, 1.0, 0.5 * (Rs + Rp))
    M01 = xp.where(tir, 0.0, Qk)
    M22 = xp.where(tir, 1.0, rs * rp)
    return M00, M01, M22


def rww_coxmunk_allm_batch(mu_o, mu_i, m_count, n_phi, ws, sigma_type,
                           n_water, q_convention):
    """Batched all-m water-water rough-Fresnel microfacet Fourier kernel over the case axis.

    mu_o (B, n_o), mu_i (B, n_i): per-case rings.  ws (B,): per-case wind.
    Returns K (B, M, n_o, n_i, 9).  Bit-identical to rww_coxmunk_allm per case.

    Memory note: the intermediate R is (B, n_o, n_i, n_phi, 9); large B*n_phi
    can be heavy, so the caller controls B via the chunk size."""
    mu_o = xp.asarray(mu_o, dtype=float)
    mu_i = xp.asarray(mu_i, dtype=float)
    B, n_o = mu_o.shape
    n_i = mu_i.shape[1]
    dphi = 2.0 * M_PI / n_phi
    phi = (xp.arange(n_phi, dtype=float) + 0.5) * dphi        # (P,)
    sigma_sq = _slope_variance_batch(ws, sigma_type)          # (B,)
    sigma_len = xp.sqrt(sigma_sq)
    ss = sigma_sq[:, None, None, None]
    sl = sigma_len[:, None, None, None]

    MO = mu_o[:, :, None, None]                               # (B, n_o, 1, 1)
    MI = mu_i[:, None, :, None]                               # (B, 1, n_i, 1)
    CP = xp.cos(phi)[None, None, None, :]
    SP = xp.sin(phi)[None, None, None, :]
    mu_in_actual = -MI
    s1 = xp.sqrt(xp.maximum(0.0, 1.0 - MO * MO))
    s2 = xp.sqrt(xp.maximum(0.0, 1.0 - mu_in_actual * mu_in_actual))
    cT = xp.clip(MO * mu_in_actual + s1 * s2 * CP, -1.0, 1.0)
    cos_omega = xp.sqrt(xp.maximum(0.0, 0.5 * (1.0 - cT)))
    cos_beta = xp.where(cos_omega > 1e-6,
                        (MO + MI) / (2.0 * xp.where(cos_omega > 0, cos_omega, 1.0)), 0.0)
    valid = (cos_omega > 1e-6) & (cos_beta > 1e-6)
    sT = xp.sqrt(xp.maximum(0.0, 1.0 - cT * cT))
    denom = xp.maximum(1e-12, sT)
    s2_safe = xp.maximum(1e-12, s2)
    s1_safe = xp.maximum(1e-12, s1)
    ci1 = (MO - mu_in_actual * cT) / (denom * s2_safe)
    si1 = s1 * SP / denom
    n1n = xp.sqrt(ci1 * ci1 + si1 * si1)
    small1 = n1n < 1e-12
    ci1 = xp.where(small1, 1.0, ci1 / xp.where(small1, 1.0, n1n))
    si1 = xp.where(small1, 0.0, si1 / xp.where(small1, 1.0, n1n))
    ci2 = (mu_in_actual - MO * cT) / (denom * s1_safe)
    si2 = s2 * SP / denom
    n2n = xp.sqrt(ci2 * ci2 + si2 * si2)
    small2 = n2n < 1e-12
    ci2 = xp.where(small2, 1.0, ci2 / xp.where(small2, 1.0, n2n))
    si2 = xp.where(small2, 0.0, si2 / xp.where(small2, 1.0, n2n))

    M00, M01, M22 = _fresnel_R_mueller_batch(cos_omega, n_water, 1.0, q_convention)

    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / xp.maximum(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * ss)) * xp.exp(-tan_beta_sq / ss)
    lam_i = _sancer_lambda_batch(MI, sl)
    lam_v = _sancer_lambda_batch(MO, sl)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    brdf = (P_slope * S_bi) / (4.0 * MI * MO * cos_beta_sq * cos_beta_sq)
    brdf = xp.where(valid, brdf, 0.0)

    c2a = ci1 * ci1 - si1 * si1
    s2a = 2.0 * ci1 * si1
    c2b = ci2 * ci2 - si2 * si2
    s2b = 2.0 * ci2 * si2
    R = xp.stack([
        M00 + 0.0 * c2a,                                     # broadcast to full shape
        M01 * c2a,
        -M01 * s2a,
        c2b * M01,
        c2b * M00 * c2a - s2b * M22 * s2a,
        -c2b * M00 * s2a - s2b * M22 * c2a,
        s2b * M01,
        s2b * M00 * c2a + c2b * M22 * s2a,
        -s2b * M00 * s2a + c2b * M22 * c2a,
    ], axis=-1)                                              # (B, n_o, n_i, P, 9)
    R = R * brdf[..., None]
    R = xp.where(valid[..., None], R, 0.0)

    ms = xp.arange(m_count)
    cm = xp.cos(ms[:, None] * phi[None, :])                 # (M, P)
    sm = xp.sin(ms[:, None] * phi[None, :])
    ce = xp.asarray(_COS_ENTRY)
    W = xp.where(ce[None, :, None], cm[:, None, :], sm[:, None, :])   # (M, 9, P)
    out = xp.einsum('boipe,mep->bmoie', R, W)               # (B, M, n_o, n_i, 9)
    norm = xp.where(ms == 0, 1.0 / (2.0 * M_PI), 1.0 / M_PI)
    out = out * (dphi * norm)[None, :, None, None, None]
    return out
