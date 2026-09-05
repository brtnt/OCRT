"""rt_air_water_coupling.c transliteration (OCRT v1.2, pristine paths).

FORWARD  couple_atm_to_water : atm BOA downwelling Fourier field -> refracted
         in-water downwelling field on the water GL grid, per mode m.
         wind>0 with atm ring weights -> full polarized rough-Fresnel microfacet air->water
         Fourier transmission operator; otherwise flat Snell + T_aw Mueller
         hybrid with the FIX-B m=0 rough-Fresnel microfacet BTDF override for wind>0.
REVERSE  couple_water_to_atm : in-water z=0- upwelling field -> transmitted
         0+ field on the atm mu grid (flat reverse Snell + T_wa Mueller).
         The OCRT_MB_CLOSURE multi-bounce re-escape is opt-in in C (default
         off) and is not ported.

Interp rule: linear only, 2-point edge extrapolation, sort+dedup at the
consumer boundary (aw_interp_on_unsorted, shared with water.py).
"""
import numpy as np
from math import pi, sqrt, cos, sin, acos, exp

from . import surface as sf
from .water import aw_interp_on_unsorted

M_PI = pi


def apply_diffuse_top_basis_sign_inplace(*fields, mode_axis=0):
    """Convert atmosphere Fourier basis ``phi+pi`` to water ``phi`` basis.

    The air-water transmission operator preserves physical azimuth, but the
    stored atmosphere BOA coefficients are reconstructed in the legacy
    ``phi+pi`` basis whereas the in-water diffuse-top source uses
    ``cos/sin(m*phi)``.  Therefore only odd Fourier modes require a sign flip.
    Apply this at the hand-off snapshot, not inside the source injector, so
    source/operator closure tests remain identities.  Works in-place for
    NumPy or CuPy arrays.
    """
    for field in fields:
        if field is None:
            continue
        sl = [slice(None)] * field.ndim
        sl[mode_axis] = slice(1, None, 2)
        field[tuple(sl)] *= -1.0
    return fields


def couple_atm_to_water(boa_I, boa_Q, boa_U, mu_atm, m_max, n_water,
                        q_convention, mu_water_pos, w_atm_pos,
                        wind_speed, sigma_type, n_phi_kernel=128):
    """rt_air_water_couple_atm_to_water.

    boa_{I,Q,U} : (m_max+1, n_a) BOA downwelling Fourier field on mu_atm
    mu_atm      : (n_a,) atm positive mu ring (ascending; may hold a
                  zero-weight view node)
    w_atm_pos   : (n_a,) ring quadrature weights, or None -> legacy hybrid
    mu_water_pos: (n_w,) water GL grid (coupling grid)
    Returns (cI, cQ, cU) each (m_max+1, n_w) in the boa normalization.
    """
    n_a = len(mu_atm)
    n_w = len(mu_water_pos)
    cI = np.zeros((m_max + 1, n_w))
    cQ = np.zeros((m_max + 1, n_w))
    cU = np.zeros((m_max + 1, n_w))

    if w_atm_pos is not None and wind_speed > 0.0:
        # Full polarized rough-Fresnel microfacet air->water Fourier transmission operator.
        #   S_w^m(i) = C_m sum_j mu_a_j w_j T_aw^m[i][j] . S_air^m(j),
        #   C_0 = 2*pi, C_{m>0} = pi, no (-1)^m (transmission).
        for m in range(m_max + 1):
            Tm = sf.T_aw_coxmunk_fourier_kernel(mu_water_pos, mu_atm, m,
                                                n_phi_kernel, wind_speed,
                                                sigma_type, n_water,
                                                q_convention)
            if Tm is None:
                break                    # kernel failure -> legacy hybrid
            az = (2.0 * M_PI) if m == 0 else M_PI
            Ia = boa_I[m]; Qa = boa_Q[m]; Ua = boa_U[m]
            for i_w in range(n_w):
                Iw = Qw = Uw = 0.0
                for j_a in range(n_a):
                    wj = az * mu_atm[j_a] * w_atm_pos[j_a]
                    if wj == 0.0:
                        continue
                    T = Tm[i_w, j_a]
                    Iw += wj * (T[0] * Ia[j_a] + T[1] * Qa[j_a] + T[2] * Ua[j_a])
                    Qw += wj * (T[3] * Ia[j_a] + T[4] * Qa[j_a] + T[5] * Ua[j_a])
                    Uw += wj * (T[6] * Ia[j_a] + T[7] * Qa[j_a] + T[8] * Ua[j_a])
                cI[m, i_w] = Iw
                cQ[m, i_w] = Qw
                cU[m, i_w] = Uw
        else:
            return cI, cQ, cU
        # fell out of the for via break -> continue into the legacy hybrid
        cI[:, :] = 0.0; cQ[:, :] = 0.0; cU[:, :] = 0.0

    # ---- legacy hybrid: flat Snell + linear interp + T_aw Mueller ----------
    mu_list = list(mu_atm)
    for j_w in range(n_w):
        mu_w = mu_water_pos[j_w]
        if mu_w <= 0.0 or mu_w > 1.0:
            continue
        sin2_water = 1.0 - mu_w * mu_w
        sin2_air = n_water * n_water * sin2_water
        if sin2_air >= 1.0:
            continue                     # TIR cone: exact zeros
        mu_air = sqrt(1.0 - sin2_air)
        M_T_aw = sf.T_aw(mu_air, n_water, q_convention)
        for m in range(m_max + 1):
            I_air = aw_interp_on_unsorted(mu_list, list(boa_I[m]), mu_air)
            Q_air = aw_interp_on_unsorted(mu_list, list(boa_Q[m]), mu_air)
            U_air = aw_interp_on_unsorted(mu_list, list(boa_U[m]), mu_air)
            cI[m, j_w] = M_T_aw[0] * I_air + M_T_aw[1] * Q_air + M_T_aw[2] * U_air
            cQ[m, j_w] = M_T_aw[3] * I_air + M_T_aw[4] * Q_air + M_T_aw[5] * U_air
            cU[m, j_w] = M_T_aw[6] * I_air + M_T_aw[7] * Q_air + M_T_aw[8] * U_air

    # FIX-B: wind>0 -> override the m=0 intensity column with the air->water
    # rough-Fresnel microfacet microfacet BTDF integral of the m=0 skylight.
    if wind_speed > 0.0:
        Nfa = 192
        Nph = 96
        from .kernel import gauss_legendre_pos
        famu, fawt = gauss_legendre_pos(Nfa)
        dph = 2.0 * M_PI / Nph
        I_atm0 = list(boa_I[0])
        for j_w in range(n_w):
            mu_w = mu_water_pos[j_w]
            if mu_w <= 1e-9 or mu_w > 1.0:
                continue
            Lw = 0.0
            for a in range(Nfa):
                mua = famu[a]
                Lair0 = aw_interp_on_unsorted(mu_list, I_atm0, mua)
                K0 = 0.0
                for p in range(Nph):
                    K0 += sf.T_aw_coxmunk_btdf_scalar(mua, mu_w, dph * p,
                                                      n_water, wind_speed,
                                                      sigma_type)
                K0 *= dph
                Lw += fawt[a] * mua * K0 * Lair0
            cI[0, j_w] = Lw
    return cI, cQ, cU


def couple_water_to_atm(I_water_per_m, Q_water_per_m, U_water_per_m,
                        mu_water_pos, m_max, n_water, q_convention,
                        mu_atm_pos):
    """rt_air_water_couple_water_to_atm (flat reverse Snell + T_wa Mueller;
    OCRT_MB_CLOSURE off).  Inputs (m_max+1, n_w); returns (m_max+1, n_a)."""
    n_w = len(mu_water_pos)
    n_a = len(mu_atm_pos)
    n2 = n_water * n_water
    wlI = np.zeros((m_max + 1, n_a))
    wlQ = np.zeros((m_max + 1, n_a))
    wlU = np.zeros((m_max + 1, n_a))
    mu_list = list(mu_water_pos)
    for j_a in range(n_a):
        mu_a = mu_atm_pos[j_a]
        sin2_air = 1.0 - mu_a * mu_a
        mu_w = sqrt(max(0.0, 1.0 - sin2_air / n2))
        M_T_wa = sf.T_wa(mu_w, n_water, q_convention)
        for m in range(m_max + 1):
            I_w = aw_interp_on_unsorted(mu_list, list(I_water_per_m[m]), mu_w)
            Q_w = aw_interp_on_unsorted(mu_list, list(Q_water_per_m[m]), mu_w)
            U_w = aw_interp_on_unsorted(mu_list, list(U_water_per_m[m]), mu_w)
            wlI[m, j_a] = M_T_wa[0] * I_w + M_T_wa[1] * Q_w + M_T_wa[2] * U_w
            wlQ[m, j_a] = M_T_wa[3] * I_w + M_T_wa[4] * Q_w + M_T_wa[5] * U_w
            wlU[m, j_a] = M_T_wa[6] * I_w + M_T_wa[7] * Q_w + M_T_wa[8] * U_w
    return wlI, wlQ, wlU
