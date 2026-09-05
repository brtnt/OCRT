"""Batched in-water SOS for LUT (full-grid) mode, GPU-capable via backend.xp.

Scope (asserted, exact for this route):
  fixed-bulk value-kernel path (phase LUT + OSOAA cap), so the polarized
  coupling kernels (gr/gt/arr/att) and the water-Rayleigh slot (ydel) are
  identically zero.  Consequences used here, each verified against the C
  single-case port (byte-identical Tier-0):
    - J_Q = J_U = 0 for every order -> upward Q/U = 0 exactly;
    - downward Q/U per order = top-BC * exp(-h/mu) (zero-source integrate);
    - their contribution to the convergence residual attains its max at k=0,
      so tracking the BC sums reproduces the C residual exactly;
    - the SOS I-dynamics never see Q/U (all couplings zero).

Batch axis b: cases sharing (sza, n_water, view-vza list) -> identical
angular ring; wavelength / IOPs / wind / depth grid vary per case.  Cases
with fewer layers are padded at the bottom with identity layers
(c=1, contrib=0), which is exact.
"""
import numpy as np
import sys
from math import pi, sqrt, exp, log, cos, sin

from .backend import xp, to_np
from . import kernel as kn
from . import surface as sf
from . import phase as ph
from . import sos as sos1
from .rww import rww_coxmunk_allm

M_PI = pi


class BatchCase:
    """Per-case scalars/tables prepared on CPU."""
    __slots__ = ('a_tot', 'b_tot', 'bb_tot', 'wind', 'wl_nm', 'tab',
                 'delta_f', 'b_rt', 'ext', 'omega', 'omega_particle',
                 'tau_max', 'z_max', 'nt', 'dtau', 'F_sun', 'wind_idx')


def prepare_case(a_tot, b_tot, bb_tot, wind, wl_nm, phase_lut_path,
                 phase_case_id=None, phase_wl_nm=0.0,
                 tau_max_target=15.0, max_tau_max_target=200.0,
                 max_z_max_m=200.0, depth_bottom_tol=1.0e-8,
                 layer_dtau_target=0.05, n_layers_water=300,
                 fixed_bulk_lmax=30):
    """Cap + depth policy, identical to water.sos_pure_fixed_bulk."""
    c = BatchCase()
    c.a_tot, c.b_tot, c.bb_tot, c.wind, c.wl_nm = a_tot, b_tot, bb_tot, wind, wl_nm
    tab = ph.phase_lut_load(phase_lut_path, phase_case_id, phase_wl_nm)
    c.delta_f = ph.phase_table_osoaa_cap(tab, 6.0, 3.0)
    c.tab = tab
    c.b_rt = b_tot * (1.0 - c.delta_f)
    c.ext = a_tot + c.b_rt
    c.omega = c.b_rt / c.ext if c.ext > 0.0 else 0.0
    c.omega_particle = c.omega
    L = fixed_bulk_lmax if 0 < fixed_bulk_lmax <= 200 else 30
    # value kernel: betal=[1,0,...]; g_particle = betal[1]/3 = 0
    g_eff = 0.0
    tau_base = tau_max_target if tau_max_target > 0.0 else 20.0
    z_req = tau_base / c.ext if c.ext > 0.0 else 0.0
    tol_depth = min(1.0e-2, max(1.0e-30, depth_bottom_tol))
    efolds = -log(tol_depth)
    transport_ext = a_tot + c.b_rt * max(0.0, 1.0 - g_eff)
    if transport_ext > 0.0:
        z_tr = efolds / transport_ext
        if z_tr > z_req:
            z_req = z_tr
    tau_max = c.ext * z_req if c.ext > 0.0 else tau_base
    if max_tau_max_target > 0.0 and tau_max > max_tau_max_target:
        tau_max = max_tau_max_target
        z_req = tau_max / c.ext if c.ext > 0.0 else z_req
    if max_z_max_m > 0.0 and z_req > max_z_max_m:
        z_req = max_z_max_m
        tau_max = c.ext * z_req if c.ext > 0.0 else tau_max
    c.tau_max = tau_max
    c.z_max = tau_max / c.ext if c.ext > 0.0 else z_req
    nt = n_layers_water
    if layer_dtau_target > 0.0 and tau_max > 0.0:
        n_need = int(np.ceil(tau_max / layer_dtau_target))
        if n_need > nt:
            nt = n_need
    c.nt = max(1, nt)
    c.dtau = tau_max / c.nt
    return c


class BatchResult:
    pass


def solve_batch(cases, sza_deg, view_vza_list, n_water,
                n_mu_water=96, m_max_water=30,
                max_iterations=10000, tolerance=1.0e-7,
                sigma_type=1, q_convention=1, phase_nphi=720,
                F_sun_list=None, first_view_idx=0, progress=False):
    """Batched fixed-bulk value-kernel water SOS.

    cases: list[BatchCase]; F_sun_list: per-case F_sun_water.
    Returns BatchResult with per-case per-m node fields and diagnostics.
    """
    B = len(cases)
    mu_sun_air = cos(sza_deg * M_PI / 180.0)
    from .water import snell_down
    mu_sun_water = snell_down(mu_sun_air, n_water)

    # --- shared ring: GL(n_mu_water) + view-list snell nodes + 3 zero slots
    gl_mu, gl_w = kn.gauss_legendre_pos(n_mu_water)
    u = sos1.UAngles(gl_mu, gl_w)
    view_nodes = []
    for v in view_vza_list:
        mu_wv = snell_down(cos(v * M_PI / 180.0), n_water)
        view_nodes.append(mu_wv)
        u.add(mu_wv)
    u.add(1.0)
    u.add(mu_sun_water)
    u.add(mu_sun_air)
    n_mu = len(u.mu)
    dirs = 2 * n_mu + 1
    mu_pos = np.array(u.mu)
    w_pos = np.array(u.w)
    rm = np.zeros(dirs)
    gb = np.zeros(dirs)
    rm[n_mu] = -mu_sun_water
    rm[n_mu + 1:] = mu_pos
    rm[:n_mu] = -mu_pos[::-1]
    gb[n_mu + 1:] = w_pos
    gb[:n_mu] = w_pos[::-1]
    node_of_view = []
    for mu_wv in view_nodes:
        j = int(np.argmin(np.abs(mu_pos - mu_wv)))
        assert abs(mu_pos[j] - mu_wv) < 1e-13
        node_of_view.append(j)

    # --- per-case value kernel planes (M, n_mu+1, dirs), CPU prep
    Mloop = m_max_water + 1
    rm_map = {j: rm[j + n_mu] for j in range(-n_mu, n_mu + 1)}
    pfm_all = np.empty((B, Mloop, n_mu + 1, dirs))
    for b, c in enumerate(cases):
        pfm_all[b] = ph.direct_phase_fourier_allm(rm_map, n_mu, Mloop,
                                                  c.tab, phase_nphi)

    # --- Rww kernels per unique wind (rows 0/3/6 only: Q_up=U_up=0)
    winds = sorted({c.wind for c in cases})
    K036 = {}
    for wv in winds:
        K = rww_coxmunk_allm(mu_pos, mu_pos, Mloop, max(176, 2 * m_max_water + 16),
                             wv, sigma_type, n_water, q_convention)
        K036[wv] = np.ascontiguousarray(K[..., [0, 3, 6]])   # (M, n, n, 3)
    for c in cases:
        c.wind_idx = winds.index(c.wind)

    # --- per-case vertical grids (padded to ntmax)
    ntmax = max(c.nt for c in cases)
    h = np.zeros((B, ntmax + 1))
    ch = np.zeros((B, ntmax + 1))
    real_layer = np.zeros((B, ntmax), dtype=bool)
    dtau_b = np.array([c.dtau for c in cases])
    xd = np.array([c.omega_particle for c in cases])
    for b, c in enumerate(cases):
        ks = np.arange(c.nt + 1)
        h[b, :c.nt + 1] = ks * c.tau_max / c.nt
        h[b, c.nt + 1:] = h[b, c.nt]
        ch[b, :c.nt + 1] = 0.5 * np.exp(-h[b, :c.nt + 1] / mu_sun_water)
        ch[b, c.nt + 1:] = 0.0            # padded-source rows unused
        real_layer[b, :c.nt] = True

    # xp transfers
    Xh = xp.asarray(h)
    Xch = xp.asarray(ch)
    Xreal = xp.asarray(real_layer)
    level_real = np.zeros((B, ntmax + 1), dtype=bool)
    for b, c in enumerate(cases):
        level_real[b, :c.nt + 1] = True
    Xlevel = xp.asarray(level_real)
    Xdtau = xp.asarray(dtau_b)[:, None, None]
    Xxd = xp.asarray(xd)[:, None, None]
    Xmu = xp.asarray(mu_pos)
    Xw = xp.asarray(w_pos)
    c_lay = xp.exp(-(Xh[:, 1:] - Xh[:, :-1])[:, :, None] /
                   Xmu[None, None, :])
    c_lay = xp.where(Xreal[:, :, None], c_lay, 1.0)

    def integrate_up(J, bc):
        """Upward sweep; J (B, nt+1, n) positive-direction source rows.
        Padded J rows must be 0."""
        a = (J[:, 1:, :] - J[:, :-1, :]) / Xdtau
        bcoef = J[:, :-1, :] - a * Xh[:, :-1, None]
        contrib = ((1.0 - c_lay) * (bcoef + a * Xmu[None, None, :]) +
                   a * (Xh[:, :-1, None] - Xh[:, 1:, None] * c_lay)) * 0.5
        contrib = xp.where(Xreal[:, :, None], contrib, 0.0)
        rad = xp.empty_like(J)
        I = bc
        rad[:, ntmax, :] = I
        for k in range(ntmax - 1, -1, -1):
            I = c_lay[:, k, :] * I + contrib[:, k, :]
            rad[:, k, :] = I
        return rad

    def integrate_dn(J, bc):
        """Downward sweep; J (B, nt+1, n) rows for direction -j (abs index)."""
        a = (J[:, 1:, :] - J[:, :-1, :]) / Xdtau
        bcoef = J[:, 1:, :] - a * Xh[:, 1:, None]
        contrib = ((1.0 - c_lay) * (bcoef + a * (-Xmu)[None, None, :]) +
                   a * (Xh[:, 1:, None] - Xh[:, :-1, None] * c_lay)) * 0.5
        contrib = xp.where(Xreal[:, :, None], contrib, 0.0)
        rad = xp.empty_like(J)
        I = bc
        rad[:, 0, :] = I
        for k in range(1, ntmax + 1):
            I = c_lay[:, k - 1, :] * I + contrib[:, k - 1, :]
            rad[:, k, :] = I
        return rad

    # legendre workspace (shared ring); assert zero polarized couplings
    ws = kn.LegendreWorkspace(n_mu, max(30, 2))

    Mst = Mloop
    I_node = np.zeros((B, Mst, n_mu))     # k=0 upward per m
    I_neg0 = np.zeros((B, n_mu))          # m=0 downward at k=0
    I_lvl1 = np.zeros((B, dirs))          # m=0 all dirs at k=1
    n_orders = np.zeros((B,), dtype=int)
    all_conv = np.ones((B,), dtype=int)
    m_exit = np.full((B,), Mst, dtype=int)
    amp_max = np.zeros((B,))
    streak = np.zeros((B,), dtype=int)
    exited = np.zeros((B,), dtype=bool)

    az0 = 2.0 * M_PI
    az1 = M_PI

    for m in range(0, Mst):
        if exited.all():
            break
        if progress:
            sys.stderr.write('\r  [SOS] mode %2d/%d  (%d cases)   ' % (m, Mst - 1, B))
            sys.stderr.flush()
        msign = 1.0  # specular water reflection preserves azimuth; no (-1)^m

        pfm_m = xp.asarray(pfm_all[:, m])            # (B, n_mu+1, dirs)
        A0 = pfm_m[:, 1:, n_mu + 1:]                  # P[jp, +k] -> (B, jp, k)
        A1 = pfm_m[:, 1:, :n_mu][:, :, ::-1]          # P[jp, -k] (k=1..n)
        A0T = xp.ascontiguousarray(xp.transpose(A0, (0, 2, 1)))   # (B, k, jp)
        A1T = xp.ascontiguousarray(xp.transpose(A1, (0, 2, 1)))
        sol_p = pfm_m[:, 0, n_mu + 1:]                # solar row +j
        sol_m = pfm_m[:, 0, :n_mu][:, ::-1]           # solar row -j

        # primary source (ydel=0, gammal=0): src = ch * xd * pfm[0]
        src_p = Xch[:, :, None] * Xxd * sol_p[:, None, :]
        src_m = Xch[:, :, None] * Xxd * sol_m[:, None, :]
        prim_p = integrate_up(src_p, xp.zeros((B, n_mu)))
        prim_m = integrate_dn(src_m, xp.zeros((B, n_mu)))

        Kw = xp.asarray(np.stack([K036[winds[c.wind_idx]][m] for c in cases]))
        wk = (az0 if m == 0 else az1) * Xmu * Xw      # (n,)

        tot_p = prim_p.copy()
        tot_m = prim_m.copy()
        SQ = xp.zeros((B, n_mu))    # sum of Q top-BCs over orders
        SU = xp.zeros((B, n_mu))
        active = xp.asarray(~exited)
        conv_mask = xp.zeros((B,), dtype=bool)
        n_ord = xp.ones((B,), dtype=int)
        resid_last = xp.zeros((B,))

        Ip, Im = prim_p, prim_m
        for it in range(2, max_iterations + 1):
            if not bool(active.any()):
                break
            # source build: J+ = xd*(A0T@wIp + A1T@wIm); J- = xd*(A1T@wIp + A0T@wIm)
            wIp = xp.transpose(Ip, (0, 2, 1)) * Xw[None, :, None]   # (B, jp, nk)
            wIm = xp.transpose(Im, (0, 2, 1)) * Xw[None, :, None]
            Jp = xp.transpose(A0T @ wIp + A1T @ wIm, (0, 2, 1)) * Xxd
            Jm = xp.transpose(A1T @ wIp + A0T @ wIm, (0, 2, 1)) * Xxd
            Jp = xp.where(Xlevel[:, :, None], Jp, 0.0)
            Jm = xp.where(Xlevel[:, :, None], Jm, 0.0)
            # top BC from previous-order upwelling at k=0 (I only; Q_up=U_up=0)
            I_up = Ip[:, 0, :]
            tb_I = msign * xp.einsum('bjk,bk->bj', Kw[..., 0], I_up * wk[None, :])
            tb_Q = msign * xp.einsum('bjk,bk->bj', Kw[..., 1], I_up * wk[None, :])
            tb_U = msign * xp.einsum('bjk,bk->bj', Kw[..., 2], I_up * wk[None, :])
            Icur_p = integrate_up(Jp, xp.zeros((B, n_mu)))
            Icur_m = integrate_dn(Jm, tb_I)
            am = active[:, None, None]
            tot_p = tot_p + xp.where(am, Icur_p, 0.0)
            tot_m = tot_m + xp.where(am, Icur_m, 0.0)
            SQ = SQ + xp.where(active[:, None], tb_Q, 0.0)
            SU = SU + xp.where(active[:, None], tb_U, 0.0)
            # residual (per case): includes Q/U via BC maxima (exact, decay<=1)
            max_in = xp.maximum(
                xp.abs(Icur_p).max(axis=(1, 2)),
                xp.abs(Icur_m).max(axis=(1, 2)))
            max_in = xp.maximum(max_in, xp.abs(tb_Q).max(axis=1))
            max_in = xp.maximum(max_in, xp.abs(tb_U).max(axis=1))
            max_t = xp.maximum(
                xp.abs(tot_p).max(axis=(1, 2)),
                xp.abs(tot_m).max(axis=(1, 2)))
            max_t = xp.maximum(max_t, xp.abs(SQ).max(axis=1))
            max_t = xp.maximum(max_t, xp.abs(SU).max(axis=1))
            resid = xp.where(max_t > 0.0, max_in / xp.where(max_t > 0, max_t, 1.0), 0.0)
            n_ord = xp.where(active, it, n_ord)
            resid_last = xp.where(active, resid, resid_last)
            newly = active & (resid < tolerance)
            conv_mask = conv_mask | newly
            active = active & ~newly
            Ip, Im = Icur_p, Icur_m

        n_ord_np = to_np(n_ord)
        conv_np = to_np(conv_mask)
        for b in range(B):
            if exited[b]:
                continue
            if n_ord_np[b] > n_orders[b]:
                n_orders[b] = int(n_ord_np[b])
            if not conv_np[b]:
                all_conv[b] = 0

        tp0 = to_np(tot_p[:, 0, :])
        I_node_m = tp0
        for b in range(B):
            if not exited[b]:
                I_node[b, m, :] = I_node_m[b]
        if m == 0:
            tm0 = to_np(tot_m[:, 0, :])
            lvl1_p = to_np(tot_p[:, 1, :])
            lvl1_m = to_np(tot_m[:, 1, :])
            for b in range(B):
                I_neg0[b] = tm0[b]
                I_lvl1[b, n_mu + 1:] = lvl1_p[b]
                I_lvl1[b, :n_mu] = lvl1_m[b][::-1]

        # early exit per case (amp at first_view node; Q/U views are 0)
        jv = node_of_view[first_view_idx]
        for b in range(B):
            if exited[b]:
                continue
            amp = abs(I_node[b, m, jv])
            if amp > amp_max[b]:
                amp_max[b] = amp
            if m >= 3 and amp_max[b] > 0.0:
                if amp < tolerance * amp_max[b]:
                    streak[b] += 1
                else:
                    streak[b] = 0
                if streak[b] >= 2:
                    exited[b] = True
                    m_exit[b] = m

    if progress:
        sys.stderr.write('\r  [SOS] done (%d modes, %d cases)        \n'
                         % (m + 1, B))
        sys.stderr.flush()

    res = BatchResult()
    res.mu_pos = mu_pos
    res.w_pos = w_pos
    res.n_mu = n_mu
    res.node_of_view = node_of_view
    res.mu_sun_air = mu_sun_air
    res.mu_sun_water = mu_sun_water
    res.I_node = I_node
    res.I_neg0 = I_neg0
    res.I_lvl1 = I_lvl1
    res.n_orders = n_orders
    res.all_conv = all_conv
    # C parity (rt_water_rt.c commit #14 ordering): the early-exit check runs
    # AFTER the view extraction (7e) but BEFORE the node store (7e-bis), so the
    # mode that TRIGGERS the exit is present in the view reconstruction but
    # absent from I_m_node (the water->air coupling input).  I_node here is the
    # inclusive store; consumers that mirror the coupling must zero mode
    # m_exit[b] for exited cases (see lut.run_lut_grid).
    res.m_exit = m_exit
    res.exited = exited
    res.ntmax = ntmax
    return res
