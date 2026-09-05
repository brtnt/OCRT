"""Vectorized vector-Legendre phase moments (GPU-friendly).

aerosol.compute_vector_legendre_gauss loops over 2*mie_n_mu (=800) Gauss
nodes, and _complete_spherical loops over ell x order (~2e4).  Both are pure
CPU (numpy) and cost ~1.5s per water case (nmg=400), which starves the GPU.

Here the node loop is turned into an array axis: the Legendre / spin-2
recurrences keep the sequential l axis but evaluate all nodes at once, the
moment accumulation is an einsum over nodes, and _complete_spherical becomes a
triangular matrix product.  Bit-for-bit against the scalar version on CPU.
"""
from .backend import xp


def _legendre_P_all_vec(l_max, xq):
    """P_l(xq) for l=0..l_max, xq array (n,).  Returns (n, l_max+1).
    Bonnet recurrence; l sequential, node axis vectorized."""
    n = xq.shape[0]
    P = xp.zeros((n, l_max + 1))
    P[:, 0] = 1.0
    if l_max >= 1:
        P[:, 1] = xq
    for l in range(1, l_max):
        P[:, l + 1] = ((2.0 * l + 1.0) * xq * P[:, l] - l * P[:, l - 1]) / (l + 1.0)
    return P


def _spin2_m0_vec(l_max, xq):
    """spin-2 basis P^l_{0,2} at m=0, xq array (n,).  Returns (n, l_max+1)."""
    n = xq.shape[0]
    out = xp.zeros((n, l_max + 1))
    if l_max < 2:
        return out
    out[:, 2] = 3.0 * (1.0 - xq * xq) / (2.0 * (6.0 ** 0.5))
    for ell in range(2, l_max):
        scale = (2.0 * ell + 1.0) / ((ell - 1.0) * (ell + 3.0)) ** 0.5
        back = ((ell - 2.0) * (ell + 2.0)) ** 0.5 / (2.0 * ell + 1.0)
        out[:, ell + 1] = scale * (xq * out[:, ell] - back * out[:, ell - 1])
    return out


def _clamped_cubic_spline_build_batch(x, Y):
    """Batched clamped_cubic_spline_build.  x (m,) shared grid; Y (B, m) per-case
    values -> Y2 (B, m) second derivatives.  The tridiagonal multipliers depend
    on x only (shared); only the RHS (u) and the back-substitution carry the
    case axis.  Bit-identical to per-case _clamped_cubic_spline_build."""
    n = int(x.shape[0])
    B = Y.shape[0]
    dx0 = x[1] - x[0]; dxn = x[n - 1] - x[n - 2]
    yp1 = (Y[:, 1] - Y[:, 0]) / dx0            # (B,)
    ypn = (Y[:, n - 1] - Y[:, n - 2]) / dxn
    M = xp.zeros(n)                            # x-only forward multiplier
    M[0] = -0.5
    U = xp.zeros((B, n))
    U[:, 0] = (3.0 / dx0) * ((Y[:, 1] - Y[:, 0]) / dx0 - yp1)
    for i in range(1, n - 1):
        dxm = x[i] - x[i - 1]; dxp = x[i + 1] - x[i]
        sig = dxm / (x[i + 1] - x[i - 1])
        p = sig * M[i - 1] + 2.0               # x-only
        M[i] = (sig - 1.0) / p
        dd = ((Y[:, i + 1] - Y[:, i]) / dxp - (Y[:, i] - Y[:, i - 1]) / dxm)
        U[:, i] = (6.0 * dd / (x[i + 1] - x[i - 1]) - sig * U[:, i - 1]) / p
    qn = 0.5
    un = (3.0 / dxn) * (ypn - (Y[:, n - 1] - Y[:, n - 2]) / dxn)
    Y2 = xp.zeros((B, n))
    Y2[:, n - 1] = (un - qn * U[:, n - 2]) / (qn * M[n - 2] + 1.0)
    for k in range(n - 2, -1, -1):
        Y2[:, k] = M[k] * Y2[:, k + 1] + U[:, k]
    return Y2


def _spline_eval_batch(x, Y, Y2, xq):
    """Batched _clamped_cubic_spline_eval.  x (m,), xq (n,) shared; Y, Y2 (B, m)
    per case -> (B, n).  Bracket indices depend on x,xq only (shared); the
    interpolation gathers Y/Y2 per case."""
    n = int(x.shape[0])
    B = Y.shape[0]
    below = xq <= x[0]
    above = xq >= x[-1]
    mid = ~(below | above)
    khi = xp.searchsorted(x, xq, side='right')
    khi = xp.clip(khi, 1, n - 1)
    klo = khi - 1
    h = x[khi] - x[klo]
    h = xp.where(h == 0, 1.0, h)
    a = (x[khi] - xq) / h                       # (n,) shared
    b = (xq - x[klo]) / h
    ca = (a ** 3 - a) * (h * h) / 6.0           # (n,)
    cb = (b ** 3 - b) * (h * h) / 6.0
    Yklo = Y[:, klo]; Ykhi = Y[:, khi]          # (B, n) gathers
    Y2klo = Y2[:, klo]; Y2khi = Y2[:, khi]
    val = (a[None, :] * Yklo + b[None, :] * Ykhi
           + ca[None, :] * Y2klo + cb[None, :] * Y2khi)     # (B, n)
    out = xp.where(below[None, :], Y[:, 0:1], 0.0)
    out = xp.where(above[None, :], Y[:, -1:], out)
    return xp.where(mid[None, :], val, out)


def water_component_moments_batch(mie, wls, L_max=200, mie_n_mu=400):
    """Batched water_component_phase_moments_gauss over cases.

    mie: one fixed component model (phyto/det/min).  wls: (B,) wavelengths (nm).
    Returns betal,gammal,alphal,zetal each (B, L+1) on the backend.  The angle
    grid / mu dedup / spline bracket are wavelength-independent (shared); only
    the phase values carry the case axis.  Bit-identical to the per-case path
    applied case by case."""
    import numpy as _np
    from .kernel import gauss_legendre_pos
    order = _np.argsort(mie.angles, kind='stable')      # angle order (shared)
    ang = mie.angles[order]
    mu_ang = _np.clip(_np.cos(ang * _np.pi / 180.0), -1.0, 1.0)
    order2 = _np.argsort(mu_ang, kind='stable')          # mu order within angle-order
    mus = mu_ang[order2]
    keep = []                                            # last of each dup group
    for i in range(len(mus)):
        if i + 1 < len(mus) and abs(mus[i + 1] - mus[i]) < 1e-14:
            continue
        keep.append(i)
    keep = _np.asarray(keep)
    x_shared = mus[keep]                                 # deduped mu grid (shared)

    pw = mie.phase_wavelengths                           # (W,)
    wl_um = _np.asarray(wls, dtype=float) / 1000.0       # (B,)
    idx = _np.clip(_np.searchsorted(pw, wl_um) - 1, 0, len(pw) - 2)
    lo = idx; hi = idx + 1
    frac = _np.clip((wl_um - pw[lo]) / (pw[hi] - pw[lo]), 0.0, 1.0)   # (B,) np.interp clamp

    def _phase_at_nodes(Pmat):
        # (n_ang_full, W) -> angle order -> interp at B wls -> mu-order+dedup -> (B, n_keep)
        Pr = Pmat[order, :]                              # (n_ang, W)
        Yb = (1.0 - frac)[:, None] * Pr[:, lo].T + frac[:, None] * Pr[:, hi].T   # (B, n_ang)
        Yb = Yb[:, order2][:, keep]                      # (B, n_keep) mu-order + dedup
        return xp.asarray(Yb)

    Y11 = _phase_at_nodes(mie.P11)
    Y12 = _phase_at_nodes(mie.P12)
    Y33 = _phase_at_nodes(mie.P33)
    xs = xp.asarray(x_shared)
    Y2_11 = _clamped_cubic_spline_build_batch(xs, Y11)
    Y2_12 = _clamped_cubic_spline_build_batch(xs, Y12)
    Y2_33 = _clamped_cubic_spline_build_batch(xs, Y33)
    mu, w = gauss_legendre_pos(mie_n_mu)
    mu = xp.asarray(mu); w = xp.asarray(w)
    xq = xp.concatenate([-mu, mu])
    p11 = _spline_eval_batch(xs, Y11, Y2_11, xq)         # (B, n_nodes)
    p12 = _spline_eval_batch(xs, Y12, Y2_12, xq)
    p33 = _spline_eval_batch(xs, Y33, Y2_33, xq)
    return compute_vector_legendre_gauss_batch(p11, p12, p33, mu, w, L_max)


def water_greek_band_batch(components, wl_nm, weights, L_max=200, mie_n_mu=400):
    """Band-level constituent Greek (shared wl).  Component moments are computed
    once at wl (they depend on wl only, not on the case), then the per-case Greek
    is weights (B, n_comp) @ moments (n_comp, L+1).  Returns be2,ga2,al2,ze2 each
    (B, L+1).  Identical result to water_greek_batch with wls=[wl]*B but computes
    each component's decomposition once instead of B times."""
    mom_be = []; mom_ga = []; mom_al = []; mom_ze = []
    for mie in components:
        be, ga, al, ze = water_component_moments_batch(mie, [wl_nm], L_max, mie_n_mu)
        mom_be.append(be[0]); mom_ga.append(ga[0])
        mom_al.append(al[0]); mom_ze.append(ze[0])
    Mbe = xp.stack(mom_be); Mga = xp.stack(mom_ga)      # (n_comp, L+1)
    Mal = xp.stack(mom_al); Mze = xp.stack(mom_ze)
    W = xp.asarray(weights)                             # (B, n_comp)
    return W @ Mbe, W @ Mga, W @ Mal, W @ Mze


def water_greek_batch(components, wls, weights, L_max=200, mie_n_mu=400):
    """Batched constituent (water) Greek: weighted sum over components.

    components: list of fixed component mie models (e.g. [phyto, det, min]).
    wls: (B,) wavelengths.  weights: (B, n_comp) = b_component / b_particle per
    case.  Returns be2,ga2,al2,ze2 each (B, L+1) on the backend."""
    B = len(wls)
    be2 = xp.zeros((B, L_max + 1)); ga2 = xp.zeros((B, L_max + 1))
    al2 = xp.zeros((B, L_max + 1)); ze2 = xp.zeros((B, L_max + 1))
    for j, mie in enumerate(components):
        be, ga, al, ze = water_component_moments_batch(mie, wls, L_max, mie_n_mu)
        wgt = xp.asarray(weights)[:, j][:, None]         # (B,1)
        be2 = be2 + wgt * be; ga2 = ga2 + wgt * ga
        al2 = al2 + wgt * al; ze2 = ze2 + wgt * ze
    return be2, ga2, al2, ze2


def _spline_eval_vec(x, y, y2, xq):
    """_clamped_cubic_spline_eval, xp arrays.  x,y,y2 (m,), xq (n,)."""
    n = len(x)
    out = xp.empty_like(xq)
    below = xq <= x[0]
    above = xq >= x[-1]
    mid = ~(below | above)
    out = xp.where(below, y[0], out)
    out = xp.where(above, y[-1], out)
    # searchsorted for all nodes; masked assembly (branch-free for GPU)
    khi = xp.searchsorted(x, xq, side='right')
    khi = xp.clip(khi, 1, n - 1)
    klo = khi - 1
    h = x[khi] - x[klo]
    h = xp.where(h == 0, 1.0, h)
    a = (x[khi] - xq) / h
    b = (xq - x[klo]) / h
    val = (a * y[klo] + b * y[khi] +
           ((a ** 3 - a) * y2[klo] + (b ** 3 - b) * y2[khi]) * (h * h) / 6.0)
    return xp.where(mid, val, out)


def _complete_spherical_vec(l_max, beta11, beta22, delta33):
    """_complete_spherical as triangular matrix products.  Returns alpha, zeta.

    zeta[ell] = lower*delta33[ell] - coup*(sd - ob)
    alpha[ell]= lower*beta11[ell] - coup*(sb - od)
    with sb/sd summed over order<ell where (ell-order) even, ob/od where odd.
    """
    L = l_max
    ellv = xp.arange(L + 1, dtype=float)[:, None]     # (L+1,1)
    ordv = xp.arange(L + 1, dtype=float)[None, :]     # (1,L+1)
    gap = ellv - ordv
    tri = (ordv < ellv) & (ellv >= 2.0)               # order<ell, ell>=2
    base = (ellv - 1.0) ** 2
    W = base - 1.5 * (gap - 1.0) * (ellv + ordv)
    even = tri & (((gap.astype(xp.int64)) % 2) == 0)
    odd = tri & (((gap.astype(xp.int64)) % 2) == 1)
    We = xp.where(even, W, 0.0)                        # (L+1,L+1)
    Wo = xp.where(odd, W, 0.0)
    sb = We @ beta11
    sd = We @ delta33
    ob = Wo @ beta22
    od = Wo @ delta33
    e = xp.arange(L + 1, dtype=float)
    denom_l = (e + 1.0) * (e + 2.0)
    denom_c = e * (e - 1.0) * (e + 1.0) * (e + 2.0)
    safe_l = xp.where(denom_l == 0, 1.0, denom_l)
    safe_c = xp.where(denom_c == 0, 1.0, denom_c)
    lower = e * (e - 1.0) / safe_l
    coup = 4.0 * (2.0 * e + 1.0) / safe_c
    zeta = lower * delta33 - coup * (sd - ob)
    alpha = lower * beta11 - coup * (sb - od)
    mask = e >= 2.0
    zeta = xp.where(mask, zeta, 0.0)
    alpha = xp.where(mask, alpha, 0.0)
    return alpha, zeta


def _complete_spherical_batch(l_max, beta11, beta22, delta33):
    """Batched _complete_spherical.  Inputs (B, L+1); returns alpha, zeta (B, L+1).
    Same math as _complete_spherical_vec but the triangular products use
    `x @ W.T` so a leading case axis broadcasts (x @ W.T == W @ x for 1-D x)."""
    L = l_max
    ellv = xp.arange(L + 1, dtype=float)[:, None]
    ordv = xp.arange(L + 1, dtype=float)[None, :]
    gap = ellv - ordv
    tri = (ordv < ellv) & (ellv >= 2.0)
    base = (ellv - 1.0) ** 2
    W = base - 1.5 * (gap - 1.0) * (ellv + ordv)
    even = tri & (((gap.astype(xp.int64)) % 2) == 0)
    odd = tri & (((gap.astype(xp.int64)) % 2) == 1)
    We = xp.where(even, W, 0.0)
    Wo = xp.where(odd, W, 0.0)
    sb = beta11 @ We.T
    sd = delta33 @ We.T
    ob = beta22 @ Wo.T
    od = delta33 @ Wo.T
    e = xp.arange(L + 1, dtype=float)
    denom_l = (e + 1.0) * (e + 2.0)
    denom_c = e * (e - 1.0) * (e + 1.0) * (e + 2.0)
    safe_l = xp.where(denom_l == 0, 1.0, denom_l)
    safe_c = xp.where(denom_c == 0, 1.0, denom_c)
    lower = e * (e - 1.0) / safe_l
    coup = 4.0 * (2.0 * e + 1.0) / safe_c
    zeta = lower * delta33 - coup * (sd - ob)
    alpha = lower * beta11 - coup * (sb - od)
    mask = e >= 2.0
    zeta = xp.where(mask, zeta, 0.0)
    alpha = xp.where(mask, alpha, 0.0)
    return alpha, zeta


def compute_vector_legendre_gauss_batch(p11_B, p12_B, p33_B, mu, w, l_max):
    """Batched Legendre projection over the case axis.

    p{11,12,33}_B: (B, 2*nmg) phase values at the GL nodes xq=[-mu,mu], stacked
    over cases.  mu,w: GL positive nodes/weights (case-independent, shared).
    Returns betal,gammal,alphal,zetal each (B, L+1), normalized by per-case
    beta0.  The Legendre tables (PL,P2) and coeffs (c) are built once; only the
    (B,n)@(n,L+1) projections carry the case axis, so nothing returns to host."""
    xq = xp.concatenate([-mu, mu])                    # (n,) shared
    wt = xp.concatenate([w, w])                       # (n,)
    PL = _legendre_P_all_vec(l_max, xq)               # (n, L+1) shared
    P2 = _spin2_m0_vec(l_max, xq)                     # (n, L+1) shared
    c = (2.0 * xp.arange(l_max + 1) + 1.0) * 0.5      # (L+1,)
    wp11 = wt[None, :] * p11_B                        # (B, n)
    wp12 = wt[None, :] * p12_B
    wp33 = wt[None, :] * p33_B
    betal = c[None, :] * (wp11 @ PL)                  # (B, L+1)
    deltal = c[None, :] * (wp33 @ PL)
    gammal = c[None, :] * (wp12 @ P2)
    alphal, zetal = _complete_spherical_batch(l_max, betal, betal, deltal)
    b0 = betal[:, 0:1]                                # (B, 1)
    return betal / b0, gammal / b0, alphal / b0, zetal / b0


def compute_vector_legendre_gauss_vec(x11, y11, y2_11, x12, y12, y2_12,
                                      x33, y33, y2_33, mu, w, l_max):
    """Vectorized compute_vector_legendre_gauss body.

    Spline tables (x,y,y2) are prebuilt (CPU, cheap); mu,w are GL positive
    nodes.  All arrays already on the active backend.  Returns
    betal, gammal, alphal, zetal (normalized by beta0).
    """
    xq = xp.concatenate([-mu, mu])                    # (2*nmg,)
    wt = xp.concatenate([w, w])
    p11 = _spline_eval_vec(x11, y11, y2_11, xq)
    p12 = _spline_eval_vec(x12, y12, y2_12, xq)
    p33 = _spline_eval_vec(x33, y33, y2_33, xq)
    PL = _legendre_P_all_vec(l_max, xq)               # (n, L+1)
    P2 = _spin2_m0_vec(l_max, xq)
    c = (2.0 * xp.arange(l_max + 1) + 1.0) * 0.5      # (L+1,)
    betal = c * ((wt * p11) @ PL)                     # (L+1,)
    deltal = c * ((wt * p33) @ PL)
    gammal = c * ((wt * p12) @ P2)
    beta22 = betal                                    # spherical: P22 = P11
    alphal, zetal = _complete_spherical_vec(l_max, betal, beta22, deltal)
    b0 = betal[0]
    return betal / b0, gammal / b0, alphal / b0, zetal / b0


def water_component_phase_moments_gauss_vec(mie, wl_nm, L_max=200, mie_n_mu=400):
    """GPU-vectorized water_component_phase_moments_gauss.  Spline tables are
    built on CPU (small angle grid); the node-sum / recurrences run on the
    active backend.  Returns numpy arrays (small (L+1,) results) so the caller's
    CPU accumulation is unchanged."""
    import numpy as np
    from .aerosol import _build_mu_table_unique, _clamped_cubic_spline_build
    from .constituent import bb_b_ratio as _bbr
    from .kernel import gauss_legendre_pos
    from .backend import to_np
    wl_um = wl_nm / 1000.0
    order = np.argsort(mie.angles, kind='stable')
    ang = mie.angles[order]

    def _blk(B):
        Br = B[order, :]
        return np.array([np.interp(wl_um, mie.phase_wavelengths, Br[k, :])
                         for k in range(len(ang))])
    p11 = _blk(mie.P11); p12 = _blk(mie.P12); p33 = _blk(mie.P33)
    x11, y11 = _build_mu_table_unique(p11, ang); y2_11 = _clamped_cubic_spline_build(x11, y11)
    x12, y12 = _build_mu_table_unique(p12, ang); y2_12 = _clamped_cubic_spline_build(x12, y12)
    x33, y33 = _build_mu_table_unique(p33, ang); y2_33 = _clamped_cubic_spline_build(x33, y33)
    mu, w = gauss_legendre_pos(mie_n_mu)
    be, ga, al, ze = compute_vector_legendre_gauss_vec(
        xp.asarray(x11), xp.asarray(y11), xp.asarray(y2_11),
        xp.asarray(x12), xp.asarray(y12), xp.asarray(y2_12),
        xp.asarray(x33), xp.asarray(y33), xp.asarray(y2_33),
        xp.asarray(mu), xp.asarray(w), L_max)
    return to_np(be), to_np(ga), to_np(al), to_np(ze), _bbr(mie, wl_nm)
