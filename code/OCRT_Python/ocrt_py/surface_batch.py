"""Vectorized rough-Fresnel microfacet surface kernels (GPU-friendly, no Python phi loop).

surface.py computes the Fourier m-mode surface kernel with a triple Python
loop (n_o x n_i x n_phi), calling scalar R_coxmunk_trig per (mu_o, mu_i, phi).
For n_phi=1024 and many angles this is the dominant CPU bottleneck and never
reaches the GPU.

Here the whole (n_o, n_i, n_phi) grid is broadcast into array ops so CuPy runs
it on the GPU in one shot.  Conditional early-returns (cos_omega/cos_beta
cutoffs, rotation normalization, Fresnel TIR) become xp.where masks; the 3x3
Mueller products become einsum; the phi integral is a weighted sum over the
phi axis.  Bit-for-bit against surface.fourier_kernel(R_coxmunk_trig, ...).
"""
from math import pi

from .backend import xp

M_PI = pi

# entry i uses cos(m*phi) when True, sin(m*phi) when False (== surface._COS_ENTRY)
_COS_ENTRY = [True, True, False, True, True, False, False, False, True]


def _erfc(x):
    """erfc on the active backend."""
    try:
        if xp.__name__ == 'cupy':
            from cupyx.scipy.special import erfc as _e
        else:
            from scipy.special import erfc as _e
        return _e(x)
    except Exception:
        # scipy 없으면 numpy 폴백 (느리지만 정확)
        import numpy as _np
        from scipy.special import erfc as _e
        return xp.asarray(_e(xp.asnumpy(x) if xp.__name__ == 'cupy' else x))


def _slope_variance(ws, sigma_type):
    """slope_variance (scalar wind).  sigma_type 1은 OCRT floor slope law(그리드 자료생산 기본)이다."""
    if sigma_type == 1:
        return 0.003 + 0.00512 * max(0.01, ws)
    sig = 0.0731 * (max(0.0, ws) ** 0.5)
    return max(sig * sig, 1e-10)


def _fresnel_rs_rp(mu_i, n1, n2):
    """FresnelCore rs/rp/is_TIR, array mu_i.  Returns (rs, rp, tir_mask)."""
    sin_i_sq = xp.maximum(0.0, 1.0 - mu_i * mu_i)
    sin_t_sq = (n1 / n2) * (n1 / n2) * sin_i_sq
    tir = sin_t_sq >= 1.0
    mu_t = xp.sqrt(xp.maximum(0.0, 1.0 - sin_t_sq))       # TIR 위치는 0 취급
    denom_s = n1 * mu_i + n2 * mu_t
    denom_p = n2 * mu_i + n1 * mu_t
    ds = xp.where(xp.abs(denom_s) < 1e-300, 1.0, denom_s)
    dp = xp.where(xp.abs(denom_p) < 1e-300, 1.0, denom_p)
    rs = (n1 * mu_i - n2 * mu_t) / ds
    rp = (n2 * mu_i - n1 * mu_t) / dp
    rs = xp.where(tir, 0.0, rs)
    rp = xp.where(tir, 0.0, rp)
    return rs, rp, tir


def _fresnel_mueller(mu_i, n1, n2, q_conv):
    """fresnel_R_mueller as 9 arrays (each broadcast shape).  Returns list of 9."""
    rs, rp, tir = _fresnel_rs_rp(mu_i, n1, n2)
    Rs = rs * rs
    Rp = rp * rp
    Qk = (0.5 * (Rp - Rs)) if q_conv == 1 else (0.5 * (Rs - Rp))
    half = 0.5 * (Rs + Rp)
    rr = rs * rp
    z = xp.zeros_like(mu_i)
    one = xp.ones_like(mu_i)
    # non-TIR values
    M = [half, Qk, z, Qk, half, z, z, z, rr]
    # TIR: identity diag(1,1,1)
    M[0] = xp.where(tir, one, M[0])
    M[1] = xp.where(tir, z, M[1])
    M[3] = xp.where(tir, z, M[3])
    M[4] = xp.where(tir, one, M[4])
    M[8] = xp.where(tir, one, M[8])
    return M


def _fresnel_T_mueller_vec(mu_i, n1, n2, q_conv):
    """fresnel_T_mueller_radiance as 9 arrays (radiance-convention flat-Fresnel
    transmission Mueller).  TIR positions -> 0.  Array mu_i."""
    rs, rp, tir = _fresnel_rs_rp(mu_i, n1, n2)
    sin_i_sq = xp.maximum(0.0, 1.0 - mu_i * mu_i)
    sin_t_sq = (n1 / n2) * (n1 / n2) * sin_i_sq
    mu_t = xp.sqrt(xp.maximum(0.0, 1.0 - sin_t_sq))
    ts = 1.0 + rs
    tp = (n1 / n2) * (1.0 + rp)
    mu_i_safe = xp.where(mu_i > 0.0, mu_i, 1.0)
    flux = (n2 * mu_t) / (n1 * mu_i_safe)
    rad = (n2 * n2) / (n1 * n1) * flux
    Ts = rad * ts * ts
    Tp = rad * tp * tp
    Qk = (0.5 * (Tp - Ts)) if q_conv == 1 else (0.5 * (Ts - Tp))
    half = 0.5 * (Ts + Tp)
    z = xp.zeros_like(mu_i)
    M = [half, Qk, z, Qk, half, z, z, z, rad * ts * tp]
    for e in range(9):
        M[e] = xp.where(tir, 0.0, M[e])
    return M


def _rotation_L(ci, si):
    """build_rotation_L as 9 arrays.  L[0]=1, L[4]=c2, L[5]=-s2, L[7]=s2, L[8]=c2."""
    c2 = ci * ci - si * si
    s2 = 2.0 * ci * si
    z = xp.zeros_like(ci)
    one = xp.ones_like(ci)
    return [one, z, z, z, c2, -s2, z, s2, c2]


def _mat3_mul(A, B):
    """3x3 x 3x3 for 9-lists A, B (each element broadcast array).  Returns 9-list."""
    C = [None] * 9
    for r in range(3):
        for c in range(3):
            acc = A[r * 3 + 0] * B[0 * 3 + c] + A[r * 3 + 1] * B[1 * 3 + c] \
                + A[r * 3 + 2] * B[2 * 3 + c]
            C[r * 3 + c] = acc
    return C


def _sancer_lambda(mu, sigma_len):
    """_sancer_lambda, array mu.  mu>=1-1e-10 또는 nu>20 이면 0."""
    valid = mu < (1.0 - 1e-10)
    sin_t = xp.sqrt(xp.maximum(1e-20, 1.0 - mu * mu))
    nu = (mu / sin_t) / sigma_len
    nu_ok = valid & (nu <= 20.0)
    nu_safe = xp.where(nu_ok, nu, 1.0)                    # 0 나눗셈 회피
    val = 0.5 * (xp.exp(-nu_safe * nu_safe) / (nu_safe * (M_PI ** 0.5))
                 - _erfc(nu_safe))
    return xp.where(nu_ok, val, 0.0)


def fourier_kernel_coxmunk_vec(mu_o, mu_i, m, n_phi_quad, ws, sigma_type,
                               n_water, q_convention, osoaa_sign_fix=False):
    """Vectorized air-side rough-Fresnel microfacet reflection Fourier m-mode kernel.

    Bit-for-bit replacement for
        surface.fourier_kernel(surface.R_coxmunk_trig, mu_o, mu_i, m,
                               n_phi_quad, ws, sigma_type, n_water,
                               q_convention, osoaa_sign_fix)
    Returns (n_o, n_i, 9).  mu_o, mu_i are 1-D arrays (any backend).
    """
    mu_o = xp.asarray(mu_o, dtype=float)
    mu_i = xp.asarray(mu_i, dtype=float)
    n_o = mu_o.shape[0]
    n_i = mu_i.shape[0]
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    dphi = 2.0 * M_PI / n_phi_quad
    sigma_sq = _slope_variance(ws, sigma_type)
    sigma_len = sigma_sq ** 0.5

    # phi quadrature nodes / weights per entry
    p = xp.arange(n_phi_quad, dtype=float)
    phis = (p + 0.5) * dphi
    cphi = xp.cos(phis)                                   # (P,)
    sphi = xp.sin(phis)
    cm = xp.cos(m * phis)
    sm = xp.sin(m * phis)

    # broadcast grid: o (n_o,1,1), i (1,n_i,1), p (1,1,P)
    mo = mu_o[:, None, None]
    mi_dn = mu_i[None, :, None]
    cph = cphi[None, None, :]
    sph = sphi[None, None, :]

    mi_act = -mi_dn
    s1 = xp.sqrt(xp.maximum(0.0, 1.0 - mo * mo))
    s2 = xp.sqrt(xp.maximum(0.0, 1.0 - mi_act * mi_act))
    cT = mo * mi_act + s1 * s2 * cph
    cT = xp.minimum(1.0, xp.maximum(-1.0, cT))
    cos_omega = xp.sqrt(xp.maximum(0.0, 0.5 * (1.0 - cT)))
    m_omega = cos_omega > 1e-6                            # else R=0

    co_safe = xp.where(m_omega, cos_omega, 1.0)
    cos_beta = (mo + mi_dn) / (2.0 * co_safe)
    m_beta = cos_beta > 1e-6                              # else R=0
    valid = m_omega & m_beta

    sT = xp.sqrt(xp.maximum(0.0, 1.0 - cT * cT))
    denom = xp.maximum(1e-12, sT)
    s2_safe = xp.maximum(1e-12, s2)
    s1_safe = xp.maximum(1e-12, s1)

    ci1 = (mo - mi_act * cT) / (denom * s2_safe)
    si1 = s1 * sph / denom
    n1 = xp.sqrt(ci1 * ci1 + si1 * si1)
    small1 = n1 < 1e-12
    n1s = xp.where(small1, 1.0, n1)
    ci1 = xp.where(small1, 1.0, ci1 / n1s)
    si1 = xp.where(small1, 0.0, si1 / n1s)

    ci2 = (mi_act - mo * cT) / (denom * s1_safe)
    si2 = s2 * sph / denom
    n2 = xp.sqrt(ci2 * ci2 + si2 * si2)
    small2 = n2 < 1e-12
    n2s = xp.where(small2, 1.0, n2)
    ci2 = xp.where(small2, 1.0, ci2 / n2s)
    si2 = xp.where(small2, 0.0, si2 / n2s)

    MF = _fresnel_mueller(cos_omega, 1.0, n_water, q_convention)
    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / xp.maximum(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * sigma_sq)) * xp.exp(-tan_beta_sq / sigma_sq)
    lam_i = _sancer_lambda(mi_dn * xp.ones_like(cT), sigma_len)
    lam_v = _sancer_lambda(mo * xp.ones_like(cT), sigma_len)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    brdf = (P_slope * S_bi) / (4.0 * mi_dn * mo * cos_beta_sq * cos_beta_sq)

    L1 = _rotation_L(ci1, si1)
    L2 = _rotation_L(ci2, si2)
    R = _mat3_mul(L2, _mat3_mul(MF, L1))                  # 9-list, each (n_o,n_i,P)

    # apply brdf and validity mask
    zero = xp.zeros_like(cT)
    for e in range(9):
        R[e] = xp.where(valid, R[e] * brdf, zero)

    # phi integral: out[o,i,e] = sum_p R[e][o,i,p] * w_p[e] * dphi * norm
    # w_p[e] = cm[p] if COS entry else sm[p]
    out = xp.zeros((n_o, n_i, 9))
    for e in range(9):
        w = cm if _COS_ENTRY[e] else sm                  # (P,)
        acc = xp.sum(R[e] * w[None, None, :], axis=2)    # (n_o,n_i)
        out = _set_slice(out, e, acc * (dphi * norm))

    if osoaa_sign_fix and m > 0:
        out = _set_slice(out, 2, -out[:, :, 2])
        out = _set_slice(out, 5, -out[:, :, 5])
    return out


def fourier_kernel_coxmunk_batch(mu_o, mu_i, m, n_phi_quad, winds, sigma_type,
                                 n_water, q_convention, osoaa_sign_fix=False):
    """Batched air-side rough-Fresnel microfacet reflection Fourier m-mode kernel over the case
    axis.  mu_o (B, n_o), mu_i (B, n_i): per-case angle sets.  winds (B,):
    per-case wind.  Returns (B, n_o, n_i, 9).  Bit-identical to the per-case
    fourier_kernel_coxmunk_vec; the leading B axis rides through the shared
    element-wise Fresnel / rotation / Sancer helpers."""
    from .rww_batch import _slope_variance_batch
    mu_o = xp.asarray(mu_o, dtype=float)
    mu_i = xp.asarray(mu_i, dtype=float)
    B, n_o = mu_o.shape
    n_i = mu_i.shape[1]
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    dphi = 2.0 * M_PI / n_phi_quad
    sigma_sq = _slope_variance_batch(winds, sigma_type)     # (B,)
    sigma_len = xp.sqrt(sigma_sq)
    ss = sigma_sq[:, None, None, None]
    sl = sigma_len[:, None, None, None]
    p = xp.arange(n_phi_quad, dtype=float)
    phis = (p + 0.5) * dphi
    cphi = xp.cos(phis); sphi = xp.sin(phis)
    cm = xp.cos(m * phis); sm = xp.sin(m * phis)
    mo = mu_o[:, :, None, None]                             # (B,n_o,1,1)
    mi_dn = mu_i[:, None, :, None]                          # (B,1,n_i,1)
    cph = cphi[None, None, None, :]
    sph = sphi[None, None, None, :]
    mi_act = -mi_dn
    s1 = xp.sqrt(xp.maximum(0.0, 1.0 - mo * mo))
    s2 = xp.sqrt(xp.maximum(0.0, 1.0 - mi_act * mi_act))
    cT = mo * mi_act + s1 * s2 * cph
    cT = xp.minimum(1.0, xp.maximum(-1.0, cT))
    cos_omega = xp.sqrt(xp.maximum(0.0, 0.5 * (1.0 - cT)))
    m_omega = cos_omega > 1e-6
    co_safe = xp.where(m_omega, cos_omega, 1.0)
    cos_beta = (mo + mi_dn) / (2.0 * co_safe)
    m_beta = cos_beta > 1e-6
    valid = m_omega & m_beta
    sT = xp.sqrt(xp.maximum(0.0, 1.0 - cT * cT))
    denom = xp.maximum(1e-12, sT)
    s2_safe = xp.maximum(1e-12, s2)
    s1_safe = xp.maximum(1e-12, s1)
    ci1 = (mo - mi_act * cT) / (denom * s2_safe)
    si1 = s1 * sph / denom
    n1 = xp.sqrt(ci1 * ci1 + si1 * si1)
    small1 = n1 < 1e-12; n1s = xp.where(small1, 1.0, n1)
    ci1 = xp.where(small1, 1.0, ci1 / n1s); si1 = xp.where(small1, 0.0, si1 / n1s)
    ci2 = (mi_act - mo * cT) / (denom * s1_safe)
    si2 = s2 * sph / denom
    n2 = xp.sqrt(ci2 * ci2 + si2 * si2)
    small2 = n2 < 1e-12; n2s = xp.where(small2, 1.0, n2)
    ci2 = xp.where(small2, 1.0, ci2 / n2s); si2 = xp.where(small2, 0.0, si2 / n2s)
    MF = _fresnel_mueller(cos_omega, 1.0, n_water, q_convention)
    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / xp.maximum(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * ss)) * xp.exp(-tan_beta_sq / ss)
    lam_i = _sancer_lambda(mi_dn * xp.ones_like(cT), sl)
    lam_v = _sancer_lambda(mo * xp.ones_like(cT), sl)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    brdf = (P_slope * S_bi) / (4.0 * mi_dn * mo * cos_beta_sq * cos_beta_sq)
    L1 = _rotation_L(ci1, si1); L2 = _rotation_L(ci2, si2)
    R = _mat3_mul(L2, _mat3_mul(MF, L1))
    zero = xp.zeros_like(cT)
    for e in range(9):
        R[e] = xp.where(valid, R[e] * brdf, zero)
    outs = []
    for e in range(9):
        w = cm if _COS_ENTRY[e] else sm
        acc = xp.sum(R[e] * w[None, None, None, :], axis=3)   # (B,n_o,n_i)
        outs.append(acc * (dphi * norm))
    out = xp.stack(outs, axis=-1)                            # (B,n_o,n_i,9)
    if osoaa_sign_fix and m > 0:
        out[..., 2] = -out[..., 2]
        out[..., 5] = -out[..., 5]
    return out


def _set_slice(out, e, val):
    """out[:, :, e] = val  (works for numpy and cupy)."""
    out[:, :, e] = val
    return out


def _aw_microfacet_geometry_vec(mu_i, mu_o, cph, sph, n_water, slope_var):
    """_aw_microfacet_geometry, broadcast arrays.  mu_i (air), mu_o (water),
    cph/sph broadcast.  Returns (si, so, ih, oh, hz, D, denom, valid_mask).
    valid False 이면 그 지점 T=0."""
    si = xp.sqrt(xp.maximum(0.0, 1.0 - mu_i * mu_i))
    so = xp.sqrt(xp.maximum(0.0, 1.0 - mu_o * mu_o))
    ix, iy, iz = si, xp.zeros_like(si), mu_i * xp.ones_like(si)
    ox = so * cph
    oy = so * sph
    oz = -mu_o * xp.ones_like(so)
    hx = -(ix + n_water * ox)
    hy = -(iy + n_water * oy)
    hz = -(iz + n_water * oz)
    hnorm = xp.sqrt(hx * hx + hy * hy + hz * hz)
    ok = hnorm >= 1e-12
    hn = xp.where(ok, hnorm, 1.0)
    hx = hx / hn; hy = hy / hn; hz = hz / hn
    flip = hz < 0.0
    hx = xp.where(flip, -hx, hx)
    hy = xp.where(flip, -hy, hy)
    hz = xp.where(flip, -hz, hz)
    ih = ix * hx + iy * hy + iz * hz
    oh = ox * hx + oy * hy + oz * hz
    ok = ok & (ih > 0.0) & (hz > 1e-9)
    hz_s = xp.where(ok, hz, 1.0)
    tan2 = (1.0 - hz_s * hz_s) / (hz_s * hz_s)
    D = xp.exp(-tan2 / slope_var) / (M_PI * slope_var * hz_s ** 4)
    denom = ih + n_water * oh
    denom = denom * denom
    # base validity requires slope_var, mu_i, mu_o > cutoffs (scalar-ish via arrays)
    ok = ok & (mu_i > 1e-9) & (mu_o > 1e-9) & (denom >= 1e-12)
    return si, so, ih, oh, hz, D, denom, ok


def fourier_kernel_T_aw_vec(mu_o_water, mu_i_air, m, n_phi_quad, ws,
                            sigma_type, n_water, q_convention):
    """Vectorized air->water rough-Fresnel microfacet transmission Fourier m-mode kernel.

    Bit-for-bit replacement for surface.T_aw_coxmunk_fourier_kernel.
    Returns (n_o, n_i, 9); n_o over mu_o_water, n_i over mu_i_air.
    """
    mu_o = xp.asarray(mu_o_water, dtype=float)
    mu_i = xp.asarray(mu_i_air, dtype=float)
    n_o = mu_o.shape[0]
    n_i = mu_i.shape[0]
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    dphi = 2.0 * M_PI / n_phi_quad
    s2v = _slope_variance(ws, sigma_type)
    n = n_water

    p = xp.arange(n_phi_quad, dtype=float)
    phis = (p + 0.5) * dphi
    cph = xp.cos(phis)[None, None, :]
    sph = xp.sin(phis)[None, None, :]
    cm = xp.cos(m * phis)
    sm = xp.sin(m * phis)

    mi = mu_i[None, :, None]                    # air, (1,n_i,1)
    mo = mu_o[:, None, None]                    # water, (n_o,1,1)

    si, so, ih, oh, hz, D, denom_h, valid = \
        _aw_microfacet_geometry_vec(mi, mo, cph, sph, n, s2v)

    # air(1)->water(n) Fresnel transmission at local incidence
    cti = ih
    sti = xp.sqrt(xp.maximum(0.0, 1.0 - cti * cti))
    stt = sti / n
    tir = stt >= 1.0
    ctt = xp.sqrt(xp.maximum(0.0, 1.0 - stt * stt))
    ds = cti + n * ctt
    dp = n * cti + ctt
    ds = xp.where(xp.abs(ds) < 1e-300, 1.0, ds)
    dp = xp.where(xp.abs(dp) < 1e-300, 1.0, dp)
    rs = (cti - n * ctt) / ds
    rp = (n * cti - ctt) / dp
    Ts = xp.where(tir, 0.0, 1.0 - rs * rs)
    Tp = xp.where(tir, 0.0, 1.0 - rp * rp)

    dh = xp.where(valid, denom_h, 1.0)
    mi_s = xp.where(valid, mi, 1.0)
    mo_s = xp.where(valid, mo, 1.0)
    common = (xp.abs(ih) * xp.abs(oh)) / (mi_s * mo_s) * (n * n * D) / dh
    ft_s = common * Ts
    ft_p = common * Tp
    ft_sp = common * xp.sqrt(xp.maximum(0.0, Ts * Tp))

    Qk = (0.5 * (ft_p - ft_s)) if q_convention == 1 else (0.5 * (ft_s - ft_p))
    half = 0.5 * (ft_s + ft_p)
    z = xp.zeros_like(half)
    MT = [half, Qk, z, Qk, half, z, z, z, ft_sp]

    # Mishchenko meridian rotations (incoming=air mi, outgoing=water mo)
    cosPsi = xp.minimum(1.0, xp.maximum(-1.0, mi * mo + si * so * cph))
    sT = xp.sqrt(xp.maximum(0.0, 1.0 - cosPsi * cosPsi))
    dnm = xp.maximum(1e-12, sT)
    s1_safe = xp.maximum(1e-12, si)
    s2_safe = xp.maximum(1e-12, so)
    ci1 = (mo - mi * cosPsi) / (dnm * s1_safe)
    si1 = so * sph / dnm
    r1 = xp.sqrt(ci1 * ci1 + si1 * si1)
    sm1 = r1 < 1e-12
    r1s = xp.where(sm1, 1.0, r1)
    ci1 = xp.where(sm1, 1.0, ci1 / r1s)
    si1 = xp.where(sm1, 0.0, si1 / r1s)
    ci2 = (mi - mo * cosPsi) / (dnm * s2_safe)
    si2 = si * sph / dnm
    r2 = xp.sqrt(ci2 * ci2 + si2 * si2)
    sm2 = r2 < 1e-12
    r2s = xp.where(sm2, 1.0, r2)
    ci2 = xp.where(sm2, 1.0, ci2 / r2s)
    si2 = xp.where(sm2, 0.0, si2 / r2s)
    L1 = _rotation_L(ci1, si1)
    L2 = _rotation_L(ci2, si2)
    R = _mat3_mul(L2, _mat3_mul(MT, L1))

    zero = xp.zeros_like(cosPsi)
    for e in range(9):
        R[e] = xp.where(valid, R[e], zero)

    out = xp.zeros((n_o, n_i, 9))
    for e in range(9):
        w = cm if _COS_ENTRY[e] else sm
        out[:, :, e] = xp.sum(R[e] * w[None, None, :], axis=2) * (dphi * norm)
    return out


def fourier_kernel_T_aw_batch(mu_o_water, mu_i_air, m, n_phi_quad, winds,
                              sigma_type, n_water, q_convention):
    """Batched air->water rough-Fresnel microfacet transmission Fourier m-mode kernel over the
    case axis.  mu_o_water (B, n_o), mu_i_air (B, n_i): per-case rings.  winds
    (B,): per-case wind.  Returns (B, n_o, n_i, 9).  Bit-identical to the
    per-case fourier_kernel_T_aw_vec."""
    from .rww_batch import _slope_variance_batch
    mu_o = xp.asarray(mu_o_water, dtype=float)
    mu_i = xp.asarray(mu_i_air, dtype=float)
    B, n_o = mu_o.shape
    n_i = mu_i.shape[1]
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    dphi = 2.0 * M_PI / n_phi_quad
    s2v = _slope_variance_batch(winds, sigma_type)[:, None, None, None]   # (B,1,1,1)
    n = n_water
    p = xp.arange(n_phi_quad, dtype=float)
    phis = (p + 0.5) * dphi
    cph = xp.cos(phis)[None, None, None, :]
    sph = xp.sin(phis)[None, None, None, :]
    cm = xp.cos(m * phis); sm = xp.sin(m * phis)
    mi = mu_i[:, None, :, None]                  # air (B,1,n_i,1)
    mo = mu_o[:, :, None, None]                  # water (B,n_o,1,1)
    si, so, ih, oh, hz, D, denom_h, valid = \
        _aw_microfacet_geometry_vec(mi, mo, cph, sph, n, s2v)
    cti = ih
    sti = xp.sqrt(xp.maximum(0.0, 1.0 - cti * cti))
    stt = sti / n
    tir = stt >= 1.0
    ctt = xp.sqrt(xp.maximum(0.0, 1.0 - stt * stt))
    ds = cti + n * ctt; dp = n * cti + ctt
    ds = xp.where(xp.abs(ds) < 1e-300, 1.0, ds)
    dp = xp.where(xp.abs(dp) < 1e-300, 1.0, dp)
    rs = (cti - n * ctt) / ds; rp = (n * cti - ctt) / dp
    Ts = xp.where(tir, 0.0, 1.0 - rs * rs); Tp = xp.where(tir, 0.0, 1.0 - rp * rp)
    dh = xp.where(valid, denom_h, 1.0)
    mi_s = xp.where(valid, mi, 1.0); mo_s = xp.where(valid, mo, 1.0)
    common = (xp.abs(ih) * xp.abs(oh)) / (mi_s * mo_s) * (n * n * D) / dh
    ft_s = common * Ts; ft_p = common * Tp
    ft_sp = common * xp.sqrt(xp.maximum(0.0, Ts * Tp))
    Qk = (0.5 * (ft_p - ft_s)) if q_convention == 1 else (0.5 * (ft_s - ft_p))
    half = 0.5 * (ft_s + ft_p)
    z = xp.zeros_like(half)
    MT = [half, Qk, z, Qk, half, z, z, z, ft_sp]
    cosPsi = xp.minimum(1.0, xp.maximum(-1.0, mi * mo + si * so * cph))
    sT = xp.sqrt(xp.maximum(0.0, 1.0 - cosPsi * cosPsi))
    dnm = xp.maximum(1e-12, sT)
    s1_safe = xp.maximum(1e-12, si); s2_safe = xp.maximum(1e-12, so)
    ci1 = (mo - mi * cosPsi) / (dnm * s1_safe); si1 = so * sph / dnm
    r1 = xp.sqrt(ci1 * ci1 + si1 * si1); sm1 = r1 < 1e-12; r1s = xp.where(sm1, 1.0, r1)
    ci1 = xp.where(sm1, 1.0, ci1 / r1s); si1 = xp.where(sm1, 0.0, si1 / r1s)
    ci2 = (mi - mo * cosPsi) / (dnm * s2_safe); si2 = si * sph / dnm
    r2 = xp.sqrt(ci2 * ci2 + si2 * si2); sm2 = r2 < 1e-12; r2s = xp.where(sm2, 1.0, r2)
    ci2 = xp.where(sm2, 1.0, ci2 / r2s); si2 = xp.where(sm2, 0.0, si2 / r2s)
    L1 = _rotation_L(ci1, si1); L2 = _rotation_L(ci2, si2)
    R = _mat3_mul(L2, _mat3_mul(MT, L1))
    zero = xp.zeros_like(cosPsi)
    for e in range(9):
        R[e] = xp.where(valid, R[e], zero)
    outs = []
    for e in range(9):
        w = cm if _COS_ENTRY[e] else sm
        outs.append(xp.sum(R[e] * w[None, None, None, :], axis=3) * (dphi * norm))
    return xp.stack(outs, axis=-1)               # (B, n_o, n_i, 9)


def fourier_kernel_T_wa_vec(mu_o_air, mu_i_water, m, n_phi_quad, ws,
                            sigma_type, n_water, q_convention):
    """Vectorized water->air rough-Fresnel microfacet BTDF Fourier m-mode kernel.

    Bit-for-bit replacement for
        surface.fourier_kernel(surface.T_coxmunk_trig, mu_o, mu_i, m, ...)
    (== T_wa_coxmunk_fourier_kernel).  Returns (n_o, n_i, 9); n_o over
    mu_o_air, n_i over mu_i_water.
    """
    mu_o = xp.asarray(mu_o_air, dtype=float)       # air
    mu_i = xp.asarray(mu_i_water, dtype=float)     # water
    n_o = mu_o.shape[0]
    n_i = mu_i.shape[0]
    norm = (1.0 / (2.0 * M_PI)) if m == 0 else (1.0 / M_PI)
    dphi = 2.0 * M_PI / n_phi_quad
    sigma_sq = _slope_variance(ws, sigma_type)
    sigma_len = sigma_sq ** 0.5
    n = n_water

    p = xp.arange(n_phi_quad, dtype=float)
    phis = (p + 0.5) * dphi
    cph = xp.cos(phis)[None, None, :]
    sph = xp.sin(phis)[None, None, :]
    cm = xp.cos(m * phis)
    sm = xp.sin(m * phis)

    mu2 = mu_o[:, None, None]                       # air out
    mu1 = mu_i[None, :, None]                        # water in
    s1 = xp.sqrt(xp.maximum(0.0, 1.0 - mu1 * mu1))
    s2 = xp.sqrt(xp.maximum(0.0, 1.0 - mu2 * mu2))
    cosPsi = xp.minimum(1.0, xp.maximum(-1.0, mu1 * mu2 + s1 * s2 * cph))
    Nsq = n * n + 1.0 - 2.0 * n * cosPsi
    ok = (mu1 > 1e-9) & (mu2 > 1e-9) & (Nsq > 1e-12)
    Nsq_s = xp.where(ok, Nsq, 1.0)
    Nmag = xp.sqrt(Nsq_s)
    cos_beta = (n * mu1 - mu2) / Nmag
    cos_omega_w = (n - cosPsi) / Nmag
    ok = ok & (cos_beta > 1e-6) & (cos_omega_w > 1e-6)

    # water(n)->air(1) Fresnel at cos_omega_w
    rs, rp, tir = _fresnel_rs_rp(cos_omega_w, n, 1.0)
    ok = ok & (~tir)
    ts = 1.0 + rs
    tp = n * (1.0 + rp)
    cos_beta_sq = cos_beta * cos_beta
    tan_beta_sq = (1.0 - cos_beta_sq) / xp.maximum(cos_beta_sq, 1e-8)
    P_slope = (1.0 / (M_PI * sigma_sq)) * xp.exp(-tan_beta_sq / sigma_sq)
    lam_i = _sancer_lambda(mu1 * xp.ones_like(cosPsi), sigma_len)
    lam_v = _sancer_lambda(mu2 * xp.ones_like(cosPsi), sigma_len)
    S_bi = 1.0 / (1.0 + lam_i + lam_v)
    # mu_af = fc.mu_t (transmitted cosine) from FresnelCore(cos_omega_w, n, 1)
    sin_i_sq = xp.maximum(0.0, 1.0 - cos_omega_w * cos_omega_w)
    sin_t_sq = (n / 1.0) * (n / 1.0) * sin_i_sq
    mu_af = xp.sqrt(xp.maximum(0.0, 1.0 - sin_t_sq))
    mu1_s = xp.where(ok, mu1, 1.0)
    mu2_s = xp.where(ok, mu2, 1.0)
    cbs = xp.where(ok, cos_beta_sq, 1.0)
    common = (mu_af * mu_af) / (n * mu1_s * mu2_s * Nsq_s) \
        * (P_slope / (cbs * cbs)) * S_bi
    ft_s = common * ts * ts
    ft_p = common * tp * tp
    ft_sp = common * ts * tp

    Qk = (0.5 * (ft_p - ft_s)) if q_convention == 1 else (0.5 * (ft_s - ft_p))
    half = 0.5 * (ft_s + ft_p)
    z = xp.zeros_like(half)
    MT = [half, Qk, z, Qk, half, z, z, z, ft_sp]

    sT = xp.sqrt(xp.maximum(0.0, 1.0 - cosPsi * cosPsi))
    denom = xp.maximum(1e-12, sT)
    s1_safe = xp.maximum(1e-12, s1)
    s2_safe = xp.maximum(1e-12, s2)
    ci1 = (mu2 - mu1 * cosPsi) / (denom * s1_safe)
    si1 = s2 * sph / denom
    r1 = xp.sqrt(ci1 * ci1 + si1 * si1)
    sm1 = r1 < 1e-12
    r1s = xp.where(sm1, 1.0, r1)
    ci1 = xp.where(sm1, 1.0, ci1 / r1s)
    si1 = xp.where(sm1, 0.0, si1 / r1s)
    ci2 = (mu1 - mu2 * cosPsi) / (denom * s2_safe)
    si2 = s1 * sph / denom
    r2 = xp.sqrt(ci2 * ci2 + si2 * si2)
    sm2 = r2 < 1e-12
    r2s = xp.where(sm2, 1.0, r2)
    ci2 = xp.where(sm2, 1.0, ci2 / r2s)
    si2 = xp.where(sm2, 0.0, si2 / r2s)
    L1 = _rotation_L(ci1, si1)
    L2 = _rotation_L(ci2, si2)
    R = _mat3_mul(L2, _mat3_mul(MT, L1))

    zero = xp.zeros_like(cosPsi)
    for e in range(9):
        R[e] = xp.where(ok, R[e], zero)

    out = xp.zeros((n_o, n_i, 9))
    for e in range(9):
        w = cm if _COS_ENTRY[e] else sm
        out[:, :, e] = xp.sum(R[e] * w[None, None, :], axis=2) * (dphi * norm)
    return out
