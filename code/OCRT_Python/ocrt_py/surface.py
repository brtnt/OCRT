"""OCRT shared/surface.c + rt_air_water.c transliteration."""
import numpy as np
from math import sqrt, exp, cos, sin, erfc, pi, fabs

M_PI = pi


SLOPE_MODEL_NAKAJIMA_TANAKA = 0
SLOPE_MODEL_OCRT_FLOOR = 1


def slope_variance(wind_speed_ms, sigma_type):
    """Return the isotropic rough-surface slope variance.

    sigma_type=1 is the OCRT floor law, not an unmodified Cox-Munk model:
    sigma^2 = 0.003 + 0.00512 * max(0.01, wind speed).  sigma_type=0 is
    the Nakajima-Tanaka relation.  The surface boundary condition is named
    separately as black_fresnel_ocean.
    """
    if sigma_type == SLOPE_MODEL_OCRT_FLOOR:
        return 0.003 + 0.00512 * max(0.01, wind_speed_ms)
    sigma = 0.0731 * sqrt(max(0.0, wind_speed_ms))
    sig2 = sigma * sigma
    return max(sig2, 1e-10)


class FresnelCore:
    __slots__ = ("mu_i", "mu_t", "rs", "rp", "is_TIR")

    def __init__(self, mu_i, n1, n2):
        self.mu_i = mu_i
        sin_i_sq = max(0.0, 1.0 - mu_i * mu_i)
        sin_t_sq = (n1 / n2) * (n1 / n2) * sin_i_sq
        if sin_t_sq >= 1.0:
            self.is_TIR = True
            self.mu_t = 0.0
            self.rs = 0.0
            self.rp = 0.0
            return
        self.is_TIR = False
        self.mu_t = sqrt(1.0 - sin_t_sq)
        self.rs = (n1 * mu_i - n2 * self.mu_t) / (n1 * mu_i + n2 * self.mu_t)
        self.rp = (n2 * mu_i - n1 * self.mu_t) / (n2 * mu_i + n1 * self.mu_t)


def fresnel_R_mueller(fc, q_convention):
    MR = np.zeros(9)
    if fc.is_TIR:
        MR[0] = MR[4] = MR[8] = 1.0
        return MR
    Rs = fc.rs * fc.rs
    Rp = fc.rp * fc.rp
    Qk = 0.5 * (Rp - Rs) if q_convention == 1 else 0.5 * (Rs - Rp)
    MR[0] = 0.5 * (Rs + Rp)
    MR[1] = Qk
    MR[3] = Qk
    MR[4] = 0.5 * (Rs + Rp)
    MR[8] = fc.rs * fc.rp
    return MR


def fresnel_T_mueller_radiance(fc, n1, n2, q_convention):
    MT = np.zeros(9)
    if fc.is_TIR:
        return MT
    ts = 1.0 + fc.rs
    tp = (n1 / n2) * (1.0 + fc.rp)
    flux_factor = (n2 * fc.mu_t) / (n1 * fc.mu_i)
    n_ratio_sq = (n2 * n2) / (n1 * n1)
    rad = n_ratio_sq * flux_factor
    Ts = rad * ts * ts
    Tp = rad * tp * tp
    Qk = 0.5 * (Tp - Ts) if q_convention == 1 else 0.5 * (Ts - Tp)
    MT[0] = 0.5 * (Ts + Tp)
    MT[1] = Qk
    MT[3] = Qk
    MT[4] = 0.5 * (Ts + Tp)
    MT[8] = rad * ts * tp
    return MT


def flat_fresnel_R_general(mu_i, n1, n2, q_convention):
    return fresnel_R_mueller(FresnelCore(mu_i, n1, n2), q_convention)


def flat_fresnel_T(mu_i, n1, n2, q_convention):
    fc = FresnelCore(mu_i, n1, n2)
    return fresnel_T_mueller_radiance(fc, n1, n2, q_convention)


# rt_air_water wrappers
def T_aw(mu_a, n_water, q_convention):
    return flat_fresnel_T(mu_a, 1.0, n_water, q_convention)


def T_wa(mu_w, n_water, q_convention):
    return flat_fresnel_T(mu_w, n_water, 1.0, q_convention)


def R_aa(mu_a, n_water, q_convention):
    return flat_fresnel_R_general(mu_a, 1.0, n_water, q_convention)


def R_ww(mu_w, n_water, q_convention):
    return flat_fresnel_R_general(mu_w, n_water, 1.0, q_convention)


def mu_refracted_down(mu_a, n_water):
    sin_a_sq = max(0.0, 1.0 - mu_a * mu_a)
    sin_w_sq = sin_a_sq / (n_water * n_water)
    if sin_w_sq >= 1.0:
        return 0.0
    return sqrt(1.0 - sin_w_sq)


def mu_critical(n_water):
    inv = 1.0 / (n_water * n_water)
    if inv >= 1.0:
        return 0.0
    return sqrt(1.0 - inv)


def build_rotation_L(ci, si):
    c2 = ci * ci - si * si
    s2 = 2.0 * ci * si
    L = np.zeros(9)
    L[0] = 1.0
    L[4] = c2
    L[5] = -s2
    L[7] = s2
    L[8] = c2
    return L


def mat3_mul(A, B):
    return (A.reshape(3, 3) @ B.reshape(3, 3)).reshape(9)


def _sancer_lambda(mu, sigma_len):
    if mu < 1.0 - 1e-10:
        sin_t = sqrt(max(1e-20, 1.0 - mu * mu))
        nu = (mu / sin_t) / sigma_len
        if nu <= 20.0:
            return 0.5 * (exp(-nu * nu) / (nu * sqrt(M_PI)) - erfc(nu))
    return 0.0


def R_coxmunk_trig(mu_out, mu_in_dn, cos_phi, sin_phi, ws, sigma_type,
                   n_water, q_convention):
    """air-side rough Fresnel reflection Mueller over a black ocean (3x3)."""
    R = np.zeros(9)
    sigma_sq = slope_variance(ws, sigma_type)
    mu_in_actual = -mu_in_dn
    s1 = sqrt(max(0.0, 1.0 - mu_out * mu_out))
    s2 = sqrt(max(0.0, 1.0 - mu_in_actual * mu_in_actual))
    cT = mu_out * mu_in_actual + s1 * s2 * cos_phi
    cT = min(1.0, max(-1.0, cT))
    cos_omega = sqrt(max(0.0, 0.5 * (1.0 - cT)))
    if cos_omega <= 1e-6:
        return R
    cos_beta = (mu_out + mu_in_dn) / (2.0 * cos_omega)
    if cos_beta <= 1e-6:
        return R
    sT = sqrt(max(0.0, 1.0 - cT * cT))
    denom = max(1e-12, sT)
    s2_safe = max(1e-12, s2)
    s1_safe = max(1e-12, s1)
    ci1 = (mu_out - mu_in_actual * cT) / (denom * s2_safe)
    si1 = s1 * sin_phi / denom
    n1 = sqrt(ci1 * ci1 + si1 * si1)
    if n1 < 1e-12:
        ci1, si1 = 1.0, 0.0
    else:
        ci1 /= n1
        si1 /= n1
    ci2 = (mu_in_actual - mu_out * cT) / (denom * s1_safe)
    si2 = s2 * sin_phi / denom
    n2 = sqrt(ci2 * ci2 + si2 * si2)
    if n2 < 1e-12:
        ci2, si2 = 1.0, 0.0
    else:
        ci2 /= n2
        si2 /= n2
    MF = fresnel_R_mueller(FresnelCore(cos_omega, 1.0, n_water), q_convention)
    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / max(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * sigma_sq)) * exp(-tan_beta_sq / sigma_sq)
    sigma_len = sqrt(sigma_sq)
    lam_i = _sancer_lambda(mu_in_dn, sigma_len)
    lam_v = _sancer_lambda(mu_out, sigma_len)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    brdf = (P_slope * S_bi) / (4.0 * mu_in_dn * mu_out * cos_beta_sq * cos_beta_sq)
    L1 = build_rotation_L(ci1, si1)
    L2 = build_rotation_L(ci2, si2)
    return mat3_mul(L2, mat3_mul(MF, L1)) * brdf


def R_ww_coxmunk_trig(mu_out, mu_in_dn, cos_phi, sin_phi, ws, sigma_type,
                      n_water, q_convention):
    """water-side rough Fresnel internal-reflection Mueller."""
    R = np.zeros(9)
    sigma_sq = slope_variance(ws, sigma_type)
    mu_in_actual = -mu_in_dn
    s1 = sqrt(max(0.0, 1.0 - mu_out * mu_out))
    s2 = sqrt(max(0.0, 1.0 - mu_in_actual * mu_in_actual))
    cT = mu_out * mu_in_actual + s1 * s2 * cos_phi
    cT = min(1.0, max(-1.0, cT))
    cos_omega = sqrt(max(0.0, 0.5 * (1.0 - cT)))
    if cos_omega <= 1e-6:
        return R
    cos_beta = (mu_out + mu_in_dn) / (2.0 * cos_omega)
    if cos_beta <= 1e-6:
        return R
    sT = sqrt(max(0.0, 1.0 - cT * cT))
    denom = max(1e-12, sT)
    s2_safe = max(1e-12, s2)
    s1_safe = max(1e-12, s1)
    ci1 = (mu_out - mu_in_actual * cT) / (denom * s2_safe)
    si1 = s1 * sin_phi / denom
    n1 = sqrt(ci1 * ci1 + si1 * si1)
    if n1 < 1e-12:
        ci1, si1 = 1.0, 0.0
    else:
        ci1 /= n1
        si1 /= n1
    ci2 = (mu_in_actual - mu_out * cT) / (denom * s1_safe)
    si2 = s2 * sin_phi / denom
    n2 = sqrt(ci2 * ci2 + si2 * si2)
    if n2 < 1e-12:
        ci2, si2 = 1.0, 0.0
    else:
        ci2 /= n2
        si2 /= n2
    MF = flat_fresnel_R_general(cos_omega, n_water, 1.0, q_convention)
    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / max(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * sigma_sq)) * exp(-tan_beta_sq / sigma_sq)
    sigma_len = sqrt(sigma_sq)
    lam_i = _sancer_lambda(mu_in_dn, sigma_len)
    lam_v = _sancer_lambda(mu_out, sigma_len)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    brdf = (P_slope * S_bi) / (4.0 * mu_in_dn * mu_out * cos_beta_sq * cos_beta_sq)
    L1 = build_rotation_L(ci1, si1)
    L2 = build_rotation_L(ci2, si2)
    return mat3_mul(L2, mat3_mul(MF, L1)) * brdf


def T_coxmunk_trig(mu_out_air, mu_in_water, cos_phi, sin_phi, ws, sigma_type,
                   n_water, q_convention):
    """water->air rough Fresnel BTDF Mueller (3x3)."""
    T = np.zeros(9)
    mu1 = mu_in_water
    mu2 = mu_out_air
    if mu1 <= 1e-9 or mu2 <= 1e-9:
        return T
    sigma_sq = slope_variance(ws, sigma_type)
    n = n_water
    s1 = sqrt(max(0.0, 1.0 - mu1 * mu1))
    s2 = sqrt(max(0.0, 1.0 - mu2 * mu2))
    cosPsi = mu1 * mu2 + s1 * s2 * cos_phi
    cosPsi = min(1.0, max(-1.0, cosPsi))
    Nsq = n * n + 1.0 - 2.0 * n * cosPsi
    if Nsq <= 1e-12:
        return T
    Nmag = sqrt(Nsq)
    cos_beta = (n * mu1 - mu2) / Nmag
    cos_omega_w = (n - cosPsi) / Nmag
    if cos_beta <= 1e-6 or cos_omega_w <= 1e-6:
        return T
    fc = FresnelCore(cos_omega_w, n, 1.0)
    if fc.is_TIR:
        return T
    mu_af = fc.mu_t
    ts = 1.0 + fc.rs
    tp = n * (1.0 + fc.rp)
    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / max(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * sigma_sq)) * exp(-tan_beta_sq / sigma_sq)
    sigma_len = sqrt(sigma_sq)
    lam_i = _sancer_lambda(mu1, sigma_len)
    lam_v = _sancer_lambda(mu2, sigma_len)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    common = (mu_af * mu_af) / (n * mu1 * mu2 * Nsq) * \
        (P_slope / (cos_beta_sq * cos_beta_sq)) * S_bi
    ft_s = common * ts * ts
    ft_p = common * tp * tp
    ft_sp = common * ts * tp
    MT = np.zeros(9)
    Qk = 0.5 * (ft_p - ft_s) if q_convention == 1 else 0.5 * (ft_s - ft_p)
    MT[0] = 0.5 * (ft_s + ft_p)
    MT[1] = Qk
    MT[3] = Qk
    MT[4] = 0.5 * (ft_s + ft_p)
    MT[8] = ft_sp
    sT = sqrt(max(0.0, 1.0 - cosPsi * cosPsi))
    denom = max(1e-12, sT)
    s1_safe = max(1e-12, s1)
    s2_safe = max(1e-12, s2)
    ci1 = (mu2 - mu1 * cosPsi) / (denom * s1_safe)
    si1 = s2 * sin_phi / denom
    r1 = sqrt(ci1 * ci1 + si1 * si1)
    if r1 < 1e-12:
        ci1, si1 = 1.0, 0.0
    else:
        ci1 /= r1
        si1 /= r1
    ci2 = (mu1 - mu2 * cosPsi) / (denom * s2_safe)
    si2 = s1 * sin_phi / denom
    r2 = sqrt(ci2 * ci2 + si2 * si2)
    if r2 < 1e-12:
        ci2, si2 = 1.0, 0.0
    else:
        ci2 /= r2
        si2 /= r2
    L1 = build_rotation_L(ci1, si1)
    L2 = build_rotation_L(ci2, si2)
    return mat3_mul(L2, mat3_mul(MT, L1))


_COS_ENTRY = np.array([1, 1, 0, 1, 1, 0, 0, 0, 1], dtype=bool)


def fourier_kernel(trig_fn, mu_o, mu_i, m, n_phi_quad, ws, sigma_type,
                   n_water, q_convention, osoaa_sign_fix=False):
    """Generic phi-quadrature Fourier m-mode kernel (matches the C loops).
    Returns array (n_o, n_i, 9).  osoaa_sign_fix flips entries 2/5 for m>0
    (only used by the air-side surface_coxmunk_fourier_kernel)."""
    n_o, n_i = len(mu_o), len(mu_i)
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    dphi = 2.0 * M_PI / n_phi_quad
    out = np.zeros((n_o, n_i, 9))
    phis = (np.arange(n_phi_quad) + 0.5) * dphi
    cphis = np.cos(phis)
    sphis = np.sin(phis)
    cm = np.cos(m * phis)
    sm = np.sin(m * phis)
    w = np.where(_COS_ENTRY[None, :], cm[:, None], sm[:, None])  # (nphi, 9)
    for jo in range(n_o):
        mo = mu_o[jo]
        for ji in range(n_i):
            mi = mu_i[ji]
            acc = np.zeros(9)
            for p in range(n_phi_quad):
                Rl = trig_fn(mo, mi, cphis[p], sphis[p], ws, sigma_type,
                             n_water, q_convention)
                acc += Rl * w[p]
            out[jo, ji, :] = acc * (dphi * norm)
            if osoaa_sign_fix and m > 0:
                out[jo, ji, 2] = -out[jo, ji, 2]
                out[jo, ji, 5] = -out[jo, ji, 5]
    return out


def R_ww_coxmunk_fourier_kernel(mu_o, mu_i, m, n_phi_quad, ws, sigma_type,
                                n_water, q_convention):
    return fourier_kernel(R_ww_coxmunk_trig, mu_o, mu_i, m, n_phi_quad, ws,
                          sigma_type, n_water, q_convention)


def T_wa_coxmunk_fourier_kernel(mu_o, mu_i, m, n_phi_quad, ws, sigma_type,
                                n_water, q_convention):
    return fourier_kernel(T_coxmunk_trig, mu_o, mu_i, m, n_phi_quad, ws,
                          sigma_type, n_water, q_convention)


def direct_sunglint_rho(mu_v, phi_v, mu_0, phi_0, ws, sigma_type, n_water,
                        tau_total, q_convention):
    """surface_direct_sunglint_rho (ws>0 rough-Fresnel microfacet branch; ws=0 handled too)."""
    dphi = phi_v - phi_0
    cos_dphi = cos(dphi)
    sin_dphi = sin(dphi)
    T_dn = exp(-tau_total / mu_0)
    T_up = exp(-tau_total / mu_v)
    trans = T_dn * T_up
    if ws < 1e-6:
        sigma_fov_deg = 0.5
        sigma_fov_rad = sigma_fov_deg * M_PI / 180.0
        sin_theta_0 = sqrt(max(0.0, 1.0 - mu_0 * mu_0))
        sigma_mu = sigma_fov_rad * max(sin_theta_0, 1e-3)
        sigma_phi = sigma_fov_rad
        dphi_eff = dphi
        while dphi_eff > M_PI:
            dphi_eff -= 2.0 * M_PI
        while dphi_eff < -M_PI:
            dphi_eff += 2.0 * M_PI
        dmu = mu_v - mu_0
        expo = -0.5 * ((dmu * dmu) / (sigma_mu * sigma_mu) +
                       (dphi_eff * dphi_eff) / (sigma_phi * sigma_phi))
        if expo < -50.0:
            return np.zeros(3)
        gw = exp(expo)
        MF = flat_fresnel_R_general(mu_0, 1.0, n_water, q_convention)
        sin2_theta_0 = 2.0 * mu_0 * sin_theta_0
        denom = max(sin2_theta_0 * mu_0 * mu_v, 1e-12)
        rho_total = 4.0 * M_PI * trans / denom
        rho_local = rho_total * gw
        return np.array([MF[0], MF[3], MF[6]]) * rho_local
    R = R_coxmunk_trig(mu_v, mu_0, cos_dphi, sin_dphi, ws, sigma_type,
                       n_water, q_convention)
    return np.array([M_PI * R[0] * trans, M_PI * R[3] * trans,
                     M_PI * R[6] * trans])


# =============================================================================
# Phase D (2026-07-19): air->water rough Fresnel transmission (forward coupling).
# Transliterations of shared/surface.c: surface_aw_microfacet_geometry_,
# surface_T_aw_coxmunk_trig, surface_T_aw_coxmunk_btdf_scalar,
# surface_T_aw_coxmunk_fourier_kernel.
# =============================================================================

def _aw_microfacet_geometry(mu_i, mu_o, cphi, sphi, n_water, slope_var):
    """Half-vector + Beckmann density for air->water transmission.
    Returns (si, so, ih, oh, hz, D, denom_sq) or None (C returns -1)."""
    if slope_var <= 1e-10 or mu_i <= 1e-9 or mu_o <= 1e-9:
        return None
    si = sqrt(max(0.0, 1.0 - mu_i * mu_i))
    so = sqrt(max(0.0, 1.0 - mu_o * mu_o))
    ix, iy, iz = si, 0.0, mu_i
    ox, oy, oz = so * cphi, so * sphi, -mu_o
    hx = -(ix + n_water * ox)
    hy = -(iy + n_water * oy)
    hz = -(iz + n_water * oz)
    hnorm = sqrt(hx * hx + hy * hy + hz * hz)
    if hnorm < 1e-12:
        return None
    hx /= hnorm; hy /= hnorm; hz /= hnorm
    if hz < 0.0:
        hx, hy, hz = -hx, -hy, -hz
    ih = ix * hx + iy * hy + iz * hz
    oh = ox * hx + oy * hy + oz * hz
    if ih <= 0.0 or hz <= 1e-9:
        return None
    tan2 = (1.0 - hz * hz) / (hz * hz)
    D = exp(-tan2 / slope_var) / (M_PI * slope_var * hz * hz * hz * hz)
    denom = ih + n_water * oh
    denom *= denom
    return si, so, ih, oh, hz, D, denom


def T_aw_coxmunk_trig(mu_i_air, mu_o_water, cphi, sphi, ws, sigma_type,
                      n_water, q_convention):
    """surface_T_aw_coxmunk_trig: air->water polarized microfacet transmission
    Mueller (3x3, flat 9-vector).  Includes the n^2 radiance factor."""
    T = np.zeros(9)
    s2v = slope_variance(ws, sigma_type)
    mu_i, mu_o = mu_i_air, mu_o_water
    n = n_water
    g = _aw_microfacet_geometry(mu_i, mu_o, cphi, sphi, n, s2v)
    if g is None:
        return T
    si, so, ih, oh, hz, D, denom_h = g
    # air(1)->water(n) Fresnel at local incidence
    cti = ih
    sti = sqrt(max(0.0, 1.0 - cti * cti))
    stt = sti / n
    if stt >= 1.0:
        Ts = 0.0
        Tp = 0.0
    else:
        ctt = sqrt(max(0.0, 1.0 - stt * stt))
        rs = (cti - n * ctt) / (cti + n * ctt)
        rp = (n * cti - ctt) / (n * cti + ctt)
        Ts = 1.0 - rs * rs
        Tp = 1.0 - rp * rp
    if denom_h < 1e-12:
        return T
    common = (abs(ih) * abs(oh)) / (mu_i * mu_o) * (n * n * D) / denom_h
    ft_s = common * Ts
    ft_p = common * Tp
    ft_sp = common * sqrt(max(0.0, Ts * Tp))

    MT = np.zeros(9)
    Q_kernel = 0.5 * (ft_p - ft_s) if q_convention == 1 else 0.5 * (ft_s - ft_p)
    MT[0 * 3 + 0] = 0.5 * (ft_s + ft_p)
    MT[0 * 3 + 1] = Q_kernel
    MT[1 * 3 + 0] = Q_kernel
    MT[1 * 3 + 1] = 0.5 * (ft_s + ft_p)
    MT[2 * 3 + 2] = ft_sp

    # Mishchenko meridian rotations, roles: incoming = air, outgoing = water
    cosPsi = mu_i * mu_o + si * so * cphi
    if cosPsi > 1.0:
        cosPsi = 1.0
    elif cosPsi < -1.0:
        cosPsi = -1.0
    sT = sqrt(max(0.0, 1.0 - cosPsi * cosPsi))
    dnm = max(1e-12, sT)
    s1_safe = max(1e-12, si)
    s2_safe = max(1e-12, so)
    ci1 = (mu_o - mu_i * cosPsi) / (dnm * s1_safe)
    si1 = so * sphi / dnm
    r1 = sqrt(ci1 * ci1 + si1 * si1)
    if r1 < 1e-12:
        ci1, si1 = 1.0, 0.0
    else:
        ci1 /= r1; si1 /= r1
    ci2 = (mu_i - mu_o * cosPsi) / (dnm * s2_safe)
    si2 = si * sphi / dnm
    r2 = sqrt(ci2 * ci2 + si2 * si2)
    if r2 < 1e-12:
        ci2, si2 = 1.0, 0.0
    else:
        ci2 /= r2; si2 /= r2
    L1 = build_rotation_L(ci1, si1)
    L2 = build_rotation_L(ci2, si2)
    MTL1 = mat3_mul(MT, L1)
    OUT = mat3_mul(L2, MTL1)
    return OUT


def T_aw_coxmunk_btdf_scalar(mu_i, mu_o, dphi, n_water, wind_speed, sigma_type):
    """surface_T_aw_coxmunk_btdf_scalar: scalar air->water microfacet BTDF."""
    s2 = slope_variance(wind_speed, sigma_type)
    g = _aw_microfacet_geometry(mu_i, mu_o, cos(dphi), sin(dphi), n_water, s2)
    if g is None:
        return 0.0
    si, so, ih, oh, hz, D, denom = g
    n = n_water
    cti = ih
    sti = sqrt(max(0.0, 1.0 - cti * cti))
    stt = sti / n
    if stt >= 1.0:
        F = 1.0
    else:
        ctt = sqrt(max(0.0, 1.0 - stt * stt))
        rs = (cti - n * ctt) / (cti + n * ctt)
        rp = (n * cti - ctt) / (n * cti + ctt)
        F = 0.5 * (rs * rs + rp * rp)
    if denom < 1e-12:
        return 0.0
    ft = (abs(ih) * abs(oh)) / (mu_i * mu_o) * (n * n * (1.0 - F) * D) / denom
    return ft if ft > 0.0 else 0.0


def T_aw_coxmunk_fourier_kernel(mu_o_water, mu_i_air, m, n_phi_quad, ws,
                                sigma_type, n_water, q_convention):
    """surface_T_aw_coxmunk_fourier_kernel: Fourier-mode air->water rough-Fresnel microfacet
    transmission Mueller table.  Returns (n_o, n_i, 9); loop order and the
    trig-argument order (air incident first) follow the C source verbatim."""
    n_o = len(mu_o_water)
    n_i = len(mu_i_air)
    if slope_variance(ws, sigma_type) <= 1e-10:
        return None                      # flat (C returns -2)
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    cos_entry = np.array([1, 1, 0, 1, 1, 0, 0, 0, 1])
    dphi = 2.0 * M_PI / n_phi_quad
    T_m = np.zeros((n_o, n_i, 9))
    for j_o in range(n_o):
        mw = mu_o_water[j_o]
        for j_i in range(n_i):
            ma = mu_i_air[j_i]
            acc = np.zeros(9)
            for p in range(n_phi_quad):
                phi = (p + 0.5) * dphi
                cm = cos(m * phi)
                sm = sin(m * phi)
                Tl = T_aw_coxmunk_trig(ma, mw, cos(phi), sin(phi), ws,
                                       sigma_type, n_water, q_convention)
                w = np.where(cos_entry == 1, cm, sm)
                acc += Tl * w
            T_m[j_o, j_i, :] = acc * (norm * dphi)
    return T_m


# Canonical low-level names.  The historical ``*_coxmunk_*`` symbols remain
# source-compatible, but the operator is a rough Fresnel microfacet boundary
# whose slope-variance law is selected independently.
R_rough_fresnel_trig = R_coxmunk_trig
R_ww_rough_fresnel_trig = R_ww_coxmunk_trig
T_wa_rough_fresnel_trig = T_coxmunk_trig
T_aw_rough_fresnel_trig = T_aw_coxmunk_trig
T_aw_rough_fresnel_btdf_scalar = T_aw_coxmunk_btdf_scalar
R_ww_rough_fresnel_fourier_kernel = R_ww_coxmunk_fourier_kernel
T_wa_rough_fresnel_fourier_kernel = T_wa_coxmunk_fourier_kernel
T_aw_rough_fresnel_fourier_kernel = T_aw_coxmunk_fourier_kernel
