"""rt_solver.c SOS engine transliteration (vectorized with numpy).

Field layout: X[k, jidx] with k=0..nt (levels), jidx=j+n_mu for signed
direction j in [-n_mu..n_mu]; j=0 = solar slot (NaN-poisoned, excluded).
"""
import numpy as np
from math import pi, exp, cos, sin

from .backend import xp

F_SOLAR_PI = pi


class Atm:
    """In-water medium built by build_inwater_atm (rt_water_rt.c)."""

    def __init__(self, nt, n_mu):
        self.n_layers = nt
        self.n_mu = n_mu
        dirs = 2 * n_mu + 1
        self.h = xp.zeros(nt + 1)
        self.ch = xp.zeros(nt + 1)
        self.xdel = xp.zeros(nt + 1)
        self.ydel = xp.zeros(nt + 1)
        self.rm = xp.zeros(dirs)   # index j+n_mu
        self.gb = xp.zeros(dirs)
        self.xpl = xp.zeros(dirs)
        self.xrl = xp.zeros(dirs)
        self.xtl = xp.zeros(dirs)
        self.beta0 = 1.0
        self.beta2 = 0.0
        self.gamma2 = 0.0
        self.alpha2 = 0.0
        self.beam_q = 0.0
        self.L_max = 0
        self.betal = None
        self.gammal = None
        self.alphal = None
        self.zetal = None
        self.mu_sun = 0.0

    def R(self, j):
        return self.rm[j + self.n_mu]

    def G(self, j):
        return self.gb[j + self.n_mu]


def build_inwater_atm(nt, n_mu_alloc, tau_total, omega_w, omega_particle,
                      mu_sun_water, n_mu_quadrature, particle_Lmax,
                      betal, gammal, alphal, zetal, gl_nodes, gl_weights):
    """build_inwater_atm (rt_water_rt.c:2115).  gl_nodes/gl_weights =
    GL positive nodes for n_mu_quadrature (precomputed by caller)."""
    from math import sqrt
    atm = Atm(nt, n_mu_alloc)
    delta_w = 0.039
    ron_w = 2.0 * (1.0 - delta_w) / (2.0 + delta_w)
    atm.mu_sun = mu_sun_water
    atm.beta0 = 1.0
    atm.beta2 = 0.5 * ron_w
    atm.gamma2 = -ron_w * sqrt(1.5)
    atm.alpha2 = 3.0 * ron_w
    k_arr = xp.arange(nt + 1, dtype=float)
    atm.h = k_arr * tau_total / nt
    atm.ch = 0.5 * xp.exp(-atm.h / mu_sun_water)
    atm.xdel = xp.full(nt + 1, omega_particle)
    atm.ydel = xp.full(nt + 1, omega_w)
    n_mu = n_mu_alloc
    nq = n_mu_quadrature
    atm.rm[n_mu] = -mu_sun_water   # solar slot j=0
    atm.gb[n_mu] = 0.0
    gln = xp.asarray(gl_nodes[:nq])
    glw = xp.asarray(gl_weights[:nq])
    atm.rm[n_mu + 1:n_mu + 1 + nq] = gln
    atm.rm[n_mu - nq:n_mu] = -xp.flip(gln)
    atm.gb[n_mu + 1:n_mu + 1 + nq] = glw
    atm.gb[n_mu - nq:n_mu] = xp.flip(glw)
    if omega_particle > 0.0 and particle_Lmax >= 0 and betal is not None:
        L = particle_Lmax
        atm.L_max = L
        atm.betal = xp.asarray(betal[:L + 1], dtype=float)
        atm.gammal = xp.zeros(L + 1) if gammal is None else xp.asarray(gammal[:L + 1])
        atm.alphal = xp.zeros(L + 1) if alphal is None else xp.asarray(alphal[:L + 1])
        atm.zetal = xp.zeros(L + 1) if zetal is None else xp.asarray(zetal[:L + 1])
    return atm


def set_ring(atm, mu_ring, w_ring):
    """Overwrite positive/negative nodes with the unified angle table (sorted
    GL + zero-weight specials).  mu_ring ascending, len = new n_mu."""
    n_new = len(mu_ring)
    atm.n_mu = n_new
    # rm/gb arrays were allocated for n_mu_alloc >= n_new; rebuild them fresh
    dirs = 2 * n_new + 1
    mu_ring = xp.asarray(mu_ring)
    w_ring = xp.asarray(w_ring)
    rm = xp.zeros(dirs)
    gb = xp.zeros(dirs)
    rm[n_new] = -atm.mu_sun
    rm[n_new + 1:] = mu_ring
    rm[:n_new] = -xp.flip(mu_ring)
    gb[n_new + 1:] = w_ring
    gb[:n_new] = xp.flip(w_ring)
    atm.rm = rm
    atm.gb = gb
    atm.xpl = xp.zeros(dirs)
    atm.xrl = xp.zeros(dirs)
    atm.xtl = xp.zeros(dirs)


class UAngles:
    """rt_angles_unified: GL(n) core + sorted-insert/dedup zero-weight specials."""

    def __init__(self, gl_nodes, gl_weights):
        self.mu = list(gl_nodes)
        self.w = list(gl_weights)

    def add(self, mu_new):
        for i, m in enumerate(self.mu):
            if abs(m - mu_new) < 1e-12:
                return i
        pos = len(self.mu)
        for i, m in enumerate(self.mu):
            if m > mu_new:
                pos = i
                break
        self.mu.insert(pos, mu_new)
        self.w.insert(pos, 0.0)
        return pos


def primary_source(atm, m, ws, pfm):
    """rt_solver_primary_source: returns src_I[k, jidx].  pfm = value/moment
    kernel phase_fourier_m[j0][k] with j0=0..n_mu (row 0 = solar), second axis
    signed offset (2n_mu+1).  ws = LegendreWorkspace (plm/rrl needed only when
    beam_q!=0 with gammal)."""
    nt, n_mu = atm.n_layers, atm.n_mu
    dirs = 2 * n_mu + 1
    beta0_m = atm.beta0 if m == 0 else 0.0
    beta2 = atm.beta2
    xpl0 = atm.xpl[n_mu]
    beam_q = atm.beam_q
    rayleigh_active = (m <= 2)
    gamma2 = atm.gamma2
    xrl0 = atm.xrl[n_mu]

    bq_i_aer = np.zeros(dirs)
    if beam_q != 0.0 and atm.gammal is not None:
        l_max = ws.l_max
        for j in range(-n_mu, n_mu + 1):
            s = 0.0
            for l in range(m, l_max + 1):
                s += ws.plm[l][j + n_mu] * ws.rrl[l][n_mu] * atm.gammal[l] \
                    if l <= atm.L_max else 0.0
            bq_i_aer[j + n_mu] = s

    sa_ray = beta0_m + beta2 * atm.xpl * xpl0            # (dirs,)
    sa_aer = pfm[0, :]                                    # solar row, (dirs,)
    base = np.empty((nt + 1, dirs))
    for k in range(nt + 1):
        val = atm.ch[k] * (atm.xdel[k] * sa_aer + atm.ydel[k] * sa_ray)
        if beam_q != 0.0:
            bq_ray = (gamma2 * atm.xpl * xrl0) if rayleigh_active else 0.0
            val = val + atm.ch[k] * beam_q * (atm.xdel[k] * bq_i_aer +
                                              atm.ydel[k] * bq_ray)
        base[k] = val
    base[:, n_mu] = np.nan
    return base


def primary_source_pol(atm, m, ws, gr_pol, gt_pol):
    """rt_solver_primary_source_pol -> (src_Q, src_U).  gr_pol/gt_pol are
    kernel tables [j0=0..n_mu][signed offset]."""
    nt, n_mu = atm.n_layers, atm.n_mu
    dirs = 2 * n_mu + 1
    gamma2 = atm.gamma2
    xpl0 = atm.xpl[n_mu]
    rayleigh_active = (m <= 2)
    beam_q = atm.beam_q
    alpha2 = atm.alpha2
    xrl0 = atm.xrl[n_mu]

    bq_q = np.zeros(dirs)
    bq_u = np.zeros(dirs)
    if beam_q != 0.0 and atm.alphal is not None and atm.zetal is not None:
        l_max = min(ws.l_max, atm.L_max)
        for j in range(-n_mu, n_mu + 1):
            sq = su = 0.0
            for l in range(m, l_max + 1):
                rj = ws.rrl[l][j + n_mu]
                tj = ws.rtl[l][j + n_mu]
                r0 = ws.rrl[l][n_mu]
                t0 = ws.rtl[l][n_mu]
                a = atm.alphal[l]
                z = atm.zetal[l]
                sq += t0 * tj * z + r0 * rj * a
                su += tj * r0 * a + rj * t0 * z
            bq_q[j + n_mu] = sq
            bq_u[j + n_mu] = su

    sb_ray = (gamma2 * atm.xrl * xpl0) if rayleigh_active else np.zeros(dirs)
    sc_ray = (gamma2 * atm.xtl * xpl0) if rayleigh_active else np.zeros(dirs)
    sb_aer = gr_pol[0, :]
    sc_aer = gt_pol[0, :]
    src_q = np.empty((nt + 1, dirs))
    src_u = np.empty((nt + 1, dirs))
    for k in range(nt + 1):
        vq = +atm.ch[k] * (atm.xdel[k] * sb_aer + atm.ydel[k] * sb_ray)
        vu = -atm.ch[k] * (atm.xdel[k] * sc_aer + atm.ydel[k] * sc_ray)
        if beam_q != 0.0:
            bq_q_ray = (alpha2 * atm.xrl * xrl0) if rayleigh_active else 0.0
            bq_u_ray = (alpha2 * atm.xtl * xrl0) if rayleigh_active else 0.0
            vq = vq + atm.ch[k] * beam_q * (atm.xdel[k] * bq_q + atm.ydel[k] * bq_q_ray)
            vu = vu - atm.ch[k] * beam_q * (atm.xdel[k] * bq_u + atm.ydel[k] * bq_u_ray)
        src_q[k] = vq
        src_u[k] = vu
    src_q[:, n_mu] = np.nan
    src_u[:, n_mu] = np.nan
    return src_q, src_u


def integrate_bcs(atm, src, surface_bc=None, top_down_bc=None):
    """rt_solver_integrate_bcs, LINEAR method, vectorized over directions.
    src[k, jidx]; surface_bc[j-1] for j>0 (bottom), top_down_bc[|j|-1] for j<0.
    Returns rad[k, jidx] (solar slot NaN)."""
    nt, n_mu = atm.n_layers, atm.n_mu
    dirs = 2 * n_mu + 1
    rad = np.full((nt + 1, dirs), np.nan)
    h = atm.h
    mu_pos = atm.rm[n_mu + 1:]           # (n_mu,) ascending
    # --- upward (j>0), sweep k = nt-1 -> 0 ---
    J = src[:, n_mu + 1:]                # (nt+1, n_mu)
    bc = surface_bc if surface_bc is not None else np.zeros(n_mu)
    I = bc.copy().astype(float)
    rad[nt, n_mu + 1:] = I
    dtau = np.diff(h)                    # (nt,) h[k+1]-h[k]
    c_up = np.exp(-dtau[:, None] / mu_pos[None, :])       # (nt, n_mu)
    a_up = (J[1:, :] - J[:-1, :]) / dtau[:, None]         # per layer k (J_next-J_now)/dtau
    b_up = J[:-1, :] - a_up * h[:-1, None]
    contrib_up = ((1.0 - c_up) * (b_up + a_up * mu_pos[None, :]) +
                  a_up * (h[:-1, None] - h[1:, None] * c_up)) * 0.5
    for k in range(nt - 1, -1, -1):
        I = c_up[k] * I + contrib_up[k]
        rad[k, n_mu + 1:] = I
    # --- downward (j<0), sweep k = 1 -> nt ---
    # signed mu = -mu_pos for node |j|; storage index n_mu - j_abs
    Jd = src[:, :n_mu][:, ::-1]          # columns ordered j=-1..-n_mu -> use abs index
    # Jd[:, j_abs-1] corresponds to direction -j_abs
    bc_dn = top_down_bc if top_down_bc is not None else np.zeros(n_mu)
    I = bc_dn.copy().astype(float)
    rad[0, :n_mu] = I[::-1]
    c_dn = np.exp(-dtau[:, None] / mu_pos[None, :])
    # layer k (1..nt): J_now = src[k], J_next = src[k-1]; a=(J_now-J_next)/dtau
    a_dn = (Jd[1:, :] - Jd[:-1, :]) / dtau[:, None]
    b_dn = Jd[1:, :] - a_dn * h[1:, None]
    mu_signed = -mu_pos
    contrib_dn = ((1.0 - c_dn) * (b_dn + a_dn * mu_signed[None, :]) +
                  a_dn * (h[1:, None] - h[:-1, None] * c_dn)) * 0.5
    for k in range(1, nt + 1):
        I = c_dn[k - 1] * I + contrib_dn[k - 1]
        rad[k, :n_mu] = I[::-1]
    return rad


class KernelTables:
    """Per-m kernel tables consumed by sos_build_source_pol.
    Arrays indexed [j0 in 0..n_mu, signed offset]."""

    def __init__(self, pfm, gr, gt, arr, art, att):
        self.pfm = pfm
        self.gr = gr
        self.gt = gt
        self.arr = arr
        self.art = art
        self.att = att


def sos_build_source_pol(atm, m, kt, I_prev, Q_prev, U_prev):
    """sos_build_source_pol (fastk expressions).  Returns (J_I, J_Q, J_U)."""
    nt, n_mu = atm.n_layers, atm.n_mu
    dirs = 2 * n_mu + 1
    beta0_m = atm.beta0 if m == 0 else 0.0
    beta2 = atm.beta2
    gamma2 = atm.gamma2
    alpha2 = atm.alpha2
    ray = (m <= 2)

    idx_p = np.arange(1, n_mu + 1)       # jp / k positive indices
    xp = atm.xpl[n_mu + idx_p]           # xpl[+j]
    xr = atm.xrl[n_mu + idx_p]
    xt = atm.xtl[n_mu + idx_p]
    yp = atm.xpl[n_mu - idx_p]           # xpl[-j]
    yr = atm.xrl[n_mu - idx_p]
    yt = atm.xtl[n_mu - idx_p]
    w = atm.gb[n_mu + idx_p]             # weights z

    # aerosol tables a0..a17: shape (k, jp) with k rows, jp cols
    P = kt.pfm
    GR = kt.gr
    GT = kt.gt
    AR = kt.arr
    AT = kt.art
    ATT = kt.att
    jp = idx_p
    kk = idx_p
    # helper index arrays
    off_p = n_mu + jp                   # +jp offset
    off_m = n_mu - jp
    a = np.empty((18, n_mu, n_mu))       # (entry, k, jp)
    a[0] = P[jp][:, off_p].T if False else P[np.ix_(jp, off_p)].T  # pfm[+jp][+k] -> (k,jp)
    # Build carefully: table[X][Y] means first index X in 0..n_mu, second signed Y.
    a[0] = P[np.ix_(jp, n_mu + kk)].T    # pfm[+jp][+k]
    a[1] = P[np.ix_(jp, n_mu - kk)].T    # pfm[+jp][-k]
    a[2] = GT[np.ix_(jp, n_mu + kk)].T   # gt[+jp][+k]
    a[3] = GT[np.ix_(jp, n_mu - kk)].T
    a[4] = GT[np.ix_(kk, n_mu + jp)]     # gt[+k][+jp] -> (k,jp)
    a[5] = GT[np.ix_(kk, n_mu - jp)]
    a[6] = GR[np.ix_(jp, n_mu + kk)].T
    a[7] = GR[np.ix_(jp, n_mu - kk)].T
    a[8] = GR[np.ix_(kk, n_mu + jp)]
    a[9] = GR[np.ix_(kk, n_mu - jp)]
    a[10] = AR[np.ix_(jp, n_mu + kk)].T
    a[11] = AR[np.ix_(jp, n_mu - kk)].T
    a[12] = AT[np.ix_(jp, n_mu + kk)].T
    a[13] = AT[np.ix_(jp, n_mu - kk)].T
    a[14] = AT[np.ix_(kk, n_mu + jp)]
    a[15] = AT[np.ix_(kk, n_mu - jp)]
    a[16] = ATT[np.ix_(jp, n_mu + kk)].T
    a[17] = ATT[np.ix_(jp, n_mu - kk)].T

    r = np.zeros((18, n_mu, n_mu))
    if ray:
        # outer products (k, jp)
        r[0] = beta0_m + beta2 * np.outer(xp, xp).T * 0 + beta2 * xp[None, :] * xp[:, None]
        r[0] = beta0_m + beta2 * xp[None, :] * xp[:, None]      # xpj*xpk -> (k,jp): xp[k]*xp[jp]
        r[1] = beta0_m + beta2 * xp[None, :] * yp[:, None]      # xpj*ypk
        r[2] = gamma2 * xp[None, :] * xt[:, None]               # xpj*xtk
        r[3] = gamma2 * xp[None, :] * yt[:, None]
        r[4] = gamma2 * xp[:, None] * xt[None, :]               # xpk*xtj
        r[5] = gamma2 * xp[:, None] * yt[None, :]
        r[6] = gamma2 * xp[None, :] * xr[:, None]
        r[7] = gamma2 * xp[None, :] * yr[:, None]
        r[8] = gamma2 * xp[:, None] * xr[None, :]
        r[9] = gamma2 * xp[:, None] * yr[None, :]
        r[10] = alpha2 * xr[None, :] * xr[:, None]
        r[11] = alpha2 * xr[None, :] * yr[:, None]
        r[12] = alpha2 * xt[None, :] * xr[:, None]
        r[13] = alpha2 * xt[None, :] * yr[:, None]
        r[14] = alpha2 * xt[:, None] * xr[None, :]
        r[15] = alpha2 * xt[:, None] * yr[None, :]
        r[16] = alpha2 * xt[None, :] * xt[:, None]
        r[17] = alpha2 * xt[None, :] * yt[:, None]

    x = atm.xdel                          # (nt+1,)
    y = atm.ydel
    # eff[entry, kk_layer, k, jp] would be huge; but x,y are layer-scalars here.
    # Follow the fastk algebra: eff = a*x + r*y with x,y per layer kk.
    # Fields: transposed planes (jp, kk_layer)
    Ip = I_prev[:, n_mu + jp].T          # (jp, nk)
    Im = I_prev[:, n_mu - jp].T
    Qp = Q_prev[:, n_mu + jp].T
    Qm = Q_prev[:, n_mu - jp].T
    Up = U_prev[:, n_mu + jp].T
    Um = U_prev[:, n_mu - jp].T

    nk = nt + 1
    # For each layer kk: kernel K(entry,k,jp) mixed by (x[kk], y[kk]).
    # aI2[k,kk] = sum_jp w[jp] * ( ... ). Use einsum with x,y broadcast.
    # Compose per-entry contraction terms:
    def C(e):
        # returns (k, jp) aerosol table and rayleigh table for entry e
        return a[e], r[e]

    # weighted fields (jp, nk)
    wIp = w[:, None] * Ip
    wIm = w[:, None] * Im
    wQp = w[:, None] * Qp
    wQm = w[:, None] * Qm
    wUp = w[:, None] * Um * 0 + w[:, None] * Up
    wUm = w[:, None] * Um

    def mixdot(e, Fp, Fm=None):
        """sum_jp K_e(k,jp)*(x[kk]*a + y[kk]*r) * F(jp,kk) ->
        returns (k, nk) = x*(a@F) + y*(r@F)."""
        ae, re = C(e)
        out = (ae @ Fp) * x[None, :]
        if ray:
            out += (re @ Fp) * y[None, :]
        return out

    # accumulate per fastk expressions
    aI2 = mixdot(0, wIp) + mixdot(8, wQp) - mixdot(4, wUp) \
        + mixdot(1, wIm) + mixdot(9, wQm) - mixdot(5, wUm)
    aI1 = mixdot(1, wIp) + mixdot(9, wQp) + mixdot(5, wUp) \
        + mixdot(0, wIm) + mixdot(8, wQm) + mixdot(4, wUm)
    aQ2 = mixdot(6, wIp) + mixdot(10, wQp) - mixdot(12, wUp) \
        + mixdot(7, wIm) + mixdot(11, wQm) + mixdot(13, wUm)
    aQ1 = mixdot(7, wIp) + mixdot(11, wQp) - mixdot(13, wUp) \
        + mixdot(6, wIm) + mixdot(10, wQm) + mixdot(12, wUm)
    aU2 = -(mixdot(2, wIp) + mixdot(14, wQp) - mixdot(16, wUp)) \
        - (-mixdot(3, wIm) + mixdot(15, wQm) - mixdot(17, wUm))
    aU1 = -(mixdot(3, wIp) - mixdot(15, wQp) - mixdot(17, wUp)) \
        - (-mixdot(2, wIm) - mixdot(14, wQm) - mixdot(16, wUm))

    J_I = np.full((nk, dirs), np.nan)
    J_Q = np.full((nk, dirs), np.nan)
    J_U = np.full((nk, dirs), np.nan)
    J_I[:, n_mu + idx_p] = aI2.T
    J_I[:, n_mu - idx_p] = aI1.T
    J_Q[:, n_mu + idx_p] = aQ2.T
    J_Q[:, n_mu - idx_p] = aQ1.T
    J_U[:, n_mu + idx_p] = aU2.T
    J_U[:, n_mu - idx_p] = aU1.T
    return J_I, J_Q, J_U


def field_max_abs(F, n_mu):
    A = np.abs(F)
    A = np.concatenate([A[:, :n_mu], A[:, n_mu + 1:]], axis=1)
    return np.nanmax(A)


def sos_pol_intrefl_rough(atm, m, kt, prim_I, prim_Q, prim_U, Rww_K,
                          max_iterations, tolerance):
    """rt_solver_sos_pol_intrefl_rough.  Rww_K[jp-1, kp-1, 9] m-mode kernel."""
    nt, n_mu = atm.n_layers, atm.n_mu
    msign = 1.0  # specular water reflection preserves azimuth; no (-1)^m
    tot_I = prim_I.copy()
    tot_Q = prim_Q.copy()
    tot_U = prim_U.copy()
    if max_iterations == 1:
        return tot_I, tot_Q, tot_U, 1, True, 0.0
    I_prev, Q_prev, U_prev = prim_I.copy(), prim_Q.copy(), prim_U.copy()
    n = 1
    residual = 0.0
    converged = False
    mu_pos = atm.rm[n_mu + 1:]
    w_pos = atm.gb[n_mu + 1:]
    az = (2.0 * F_SOLAR_PI) if m == 0 else F_SOLAR_PI
    wk = az * mu_pos * w_pos                     # (kp,)
    K = Rww_K                                    # (jp, kp, 9)
    for it in range(2, max_iterations + 1):
        J_I, J_Q, J_U = sos_build_source_pol(atm, m, kt, I_prev, Q_prev, U_prev)
        # top BC from previous-order upwelling at k=0
        I_up = I_prev[0, n_mu + 1:]
        Q_up = Q_prev[0, n_mu + 1:]
        U_up = U_prev[0, n_mu + 1:]
        SI = K[:, :, 0] * I_up[None, :] + K[:, :, 1] * Q_up[None, :] + K[:, :, 2] * U_up[None, :]
        SQ = K[:, :, 3] * I_up[None, :] + K[:, :, 4] * Q_up[None, :] + K[:, :, 5] * U_up[None, :]
        SU = K[:, :, 6] * I_up[None, :] + K[:, :, 7] * Q_up[None, :] + K[:, :, 8] * U_up[None, :]
        tb_I = msign * (SI * wk[None, :]).sum(axis=1)
        tb_Q = msign * (SQ * wk[None, :]).sum(axis=1)
        tb_U = msign * (SU * wk[None, :]).sum(axis=1)
        I_curr = integrate_bcs(atm, J_I, None, tb_I)
        Q_curr = integrate_bcs(atm, J_Q, None, tb_Q)
        U_curr = integrate_bcs(atm, J_U, None, tb_U)
        tot_I += np.nan_to_num(I_curr, nan=0.0) * 0 + np.where(np.isnan(I_curr), 0.0, I_curr)
        tot_I[:, n_mu] = np.nan
        tot_Q += np.where(np.isnan(Q_curr), 0.0, Q_curr)
        tot_Q[:, n_mu] = np.nan
        tot_U += np.where(np.isnan(U_curr), 0.0, U_curr)
        tot_U[:, n_mu] = np.nan
        n = it
        max_in = max(field_max_abs(I_curr, n_mu), field_max_abs(Q_curr, n_mu),
                     field_max_abs(U_curr, n_mu))
        max_t = max(field_max_abs(tot_I, n_mu), field_max_abs(tot_Q, n_mu),
                    field_max_abs(tot_U, n_mu))
        residual = (max_in / max_t) if max_t > 0.0 else 0.0
        if residual < tolerance:
            converged = True
            break
        I_prev, Q_prev, U_prev = I_curr, Q_curr, U_curr
    return tot_I, tot_Q, tot_U, n, converged, residual


def sos_atm_black(atm, m, kt, prim_I, prim_Q, prim_U,
                  max_iterations, tolerance):
    """atm-only vector SOS (rt_solver_sos_pol) with a BLACK surface:
    no surface reflection BC, no top-of-atmosphere diffuse down.
    g_tail = 1 (OCRT_M1_MS_TAIL_MODE off = production baseline)."""
    nt, n_mu = atm.n_layers, atm.n_mu
    tot_I = prim_I.copy()
    tot_Q = prim_Q.copy()
    tot_U = prim_U.copy()
    if max_iterations == 1:
        return tot_I, tot_Q, tot_U, 1, True, 0.0
    I_prev, Q_prev, U_prev = prim_I.copy(), prim_Q.copy(), prim_U.copy()
    n = 1
    residual = 0.0
    converged = False
    for it in range(2, max_iterations + 1):
        J_I, J_Q, J_U = sos_build_source_pol(atm, m, kt, I_prev, Q_prev, U_prev)
        I_curr = integrate_bcs(atm, J_I, None, None)   # black: no surface BC
        Q_curr = integrate_bcs(atm, J_Q, None, None)
        U_curr = integrate_bcs(atm, J_U, None, None)
        tot_I += np.where(np.isnan(I_curr), 0.0, I_curr); tot_I[:, n_mu] = np.nan
        tot_Q += np.where(np.isnan(Q_curr), 0.0, Q_curr); tot_Q[:, n_mu] = np.nan
        tot_U += np.where(np.isnan(U_curr), 0.0, U_curr); tot_U[:, n_mu] = np.nan
        n = it
        max_in = max(field_max_abs(I_curr, n_mu), field_max_abs(Q_curr, n_mu),
                     field_max_abs(U_curr, n_mu))
        max_t = max(field_max_abs(tot_I, n_mu), field_max_abs(tot_Q, n_mu),
                    field_max_abs(tot_U, n_mu))
        residual = (max_in / max_t) if max_t > 0.0 else 0.0
        if residual < tolerance:
            converged = True
            break
        I_prev, Q_prev, U_prev = I_curr, Q_curr, U_curr
    return tot_I, tot_Q, tot_U, n, converged, residual


def sos_atm_surface(atm, m, kt, prim_I, prim_Q, prim_U, R_m, mu_factor,
                    max_iterations, tolerance, surf_seed=None, ext_init=None):
    """atm SOS with a reflecting surface (rt_solver_sos_pol_with_surface,
    Cox-Munk diffuse-reflection path).  Surface BC at k=nt: downward field
    reflected upward via the m-mode Mueller kernel R_m (shape (n_mu,n_mu,9)).
    surf_seed = (sI,sQ,sU) is the first-order direct-sun surface reflection
    (r_first): added to the seed and total, then subtracted at the end to get
    the sunglint-decoupled (wo_SG) result.  g_tail = 1.  mu_factor = 2pi (m=0)
    or pi (m>=1)."""
    nt, n_mu = atm.n_layers, atm.n_mu
    if surf_seed is not None:
        sI, sQ, sU = surf_seed
        tot_I = prim_I + sI
        tot_Q = prim_Q + sQ
        tot_U = prim_U + sU
    else:
        tot_I = prim_I.copy()
        tot_Q = prim_Q.copy()
        tot_U = prim_U.copy()
    if ext_init is not None:
        # add_external_bottom_source: seed AND total, never subtracted.
        eI, eQ, eU = ext_init
        tot_I = tot_I + eI
        tot_Q = tot_Q + eQ
        tot_U = tot_U + eU
    if max_iterations == 1:
        out_I, out_Q, out_U = tot_I, tot_Q, tot_U
        if surf_seed is not None:
            out_I = out_I - sI; out_Q = out_Q - sQ; out_U = out_U - sU
        return out_I, out_Q, out_U, 1, True, 0.0
    I_prev, Q_prev, U_prev = tot_I.copy(), tot_Q.copy(), tot_U.copy()
    n = 1
    residual = 0.0
    converged = False
    mu_pos = atm.rm[n_mu + 1:]
    w_pos = atm.gb[n_mu + 1:]
    mu_w = mu_factor * mu_pos * w_pos          # (n_mu,); view node w=0 -> 0
    RII = R_m[:, :, 0]; RIQ = R_m[:, :, 1]; RIU = R_m[:, :, 2]
    RQI = R_m[:, :, 3]; RQQ = R_m[:, :, 4]; RQU = R_m[:, :, 5]
    RUI = R_m[:, :, 6]; RUQ = R_m[:, :, 7]; RUU = R_m[:, :, 8]
    for it in range(2, max_iterations + 1):
        J_I, J_Q, J_U = sos_build_source_pol(atm, m, kt, I_prev, Q_prev, U_prev)
        Idn = I_prev[nt, :n_mu][::-1]
        Qdn = Q_prev[nt, :n_mu][::-1]
        Udn = U_prev[nt, :n_mu][::-1]
        SI = RII * Idn[None, :] + RIQ * Qdn[None, :] + RIU * Udn[None, :]
        SQ = RQI * Idn[None, :] + RQQ * Qdn[None, :] + RQU * Udn[None, :]
        SU = RUI * Idn[None, :] + RUQ * Qdn[None, :] + RUU * Udn[None, :]
        surf_I = (SI * mu_w[None, :]).sum(axis=1)
        surf_Q = (SQ * mu_w[None, :]).sum(axis=1)
        surf_U = (SU * mu_w[None, :]).sum(axis=1)
        I_curr = integrate_bcs(atm, J_I, surf_I, None)
        Q_curr = integrate_bcs(atm, J_Q, surf_Q, None)
        U_curr = integrate_bcs(atm, J_U, surf_U, None)
        tot_I += np.where(np.isnan(I_curr), 0.0, I_curr); tot_I[:, n_mu] = np.nan
        tot_Q += np.where(np.isnan(Q_curr), 0.0, Q_curr); tot_Q[:, n_mu] = np.nan
        tot_U += np.where(np.isnan(U_curr), 0.0, U_curr); tot_U[:, n_mu] = np.nan
        n = it
        max_in = max(field_max_abs(I_curr, n_mu), field_max_abs(Q_curr, n_mu),
                     field_max_abs(U_curr, n_mu))
        max_t = max(field_max_abs(tot_I, n_mu), field_max_abs(tot_Q, n_mu),
                    field_max_abs(tot_U, n_mu))
        residual = (max_in / max_t) if max_t > 0.0 else 0.0
        if residual < tolerance:
            converged = True
            break
        I_prev, Q_prev, U_prev = I_curr, Q_curr, U_curr

    # wo_SG: subtract the first-order direct surface reflection (sunglint decoupled)
    if surf_seed is not None:
        tot_I = tot_I - sI; tot_I[:, n_mu] = np.nan
        tot_Q = tot_Q - sQ; tot_Q[:, n_mu] = np.nan
        tot_U = tot_U - sU; tot_U[:, n_mu] = np.nan
    return tot_I, tot_Q, tot_U, n, converged, residual


def reconstruct_phi(per_m, m_max, dphi_rad):
    base = dphi_rad + F_SOLAR_PI
    I = per_m[0]
    for m in range(1, m_max + 1):
        I += 2.0 * per_m[m] * cos(m * base)
    return I


def reconstruct_phi_sin(per_m, m_max, dphi_rad):
    base = dphi_rad + F_SOLAR_PI
    U = 0.0
    for m in range(1, m_max + 1):
        U += 2.0 * per_m[m] * sin(m * base)
    return U
