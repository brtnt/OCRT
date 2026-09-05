/* rt_air_water_coupling.h — atmospheric BOA → in-water boundary mapping
 *
 * Phase B.4 Stage 2b step2: Given the atmospheric SOS BOA downward Stokes
 * field at air-side quadrature nodes, produce the in-water inbound Stokes
 * field at water-side quadrature nodes by combining:
 *   1) Snell refraction: μ_air = sqrt(1 - n_w² × (1 - μ_water²))
 *      Always real for downward refraction (sin(θ_water) ≤ 1/n_w).
 *   2) Linear interpolation in μ_air across atm quadrature grid (TODO 1 will
 *      revisit higher-order interp, e.g. PCHIP, once OSOAA benchmark is set).
 *   3) T_aw radiance Mueller (water-side incidence — same μ as the in-water
 *      direction we're filling, by reversibility of the refraction relation).
 *
 * Output is per Fourier m mode, per positive water μ node.
 * The m=0 mode reconstruction gives the hemispheric diffuse Ed_diff_water:
 *   Ed_diff_water = 2π Σ w_j_water μ_j_water · I^{m=0}_water_inbound(μ_j_water)
 */
#ifndef OCRT_V2_RT_AIR_WATER_COUPLING_H
#define OCRT_V2_RT_AIR_WATER_COUPLING_H

#include "rt_types.h"

/* Coupling output struct.
 * Pre-allocated by caller (helper available below).
 * Layout: field[m * n_mu_water + (j-1)] for m ∈ [0, m_max], j ∈ [1, n_mu_water]
 */
typedef struct {
    int     m_max;
    int     n_mu_water;
    double *mu_water_pos;         /* [n_mu_water]: positive water-side μ grid */
    double *I_inwater_per_m;      /* [(m_max+1) × n_mu_water] in-water inbound I^m */
    double *Q_inwater_per_m;      /* same for Q */
    double *U_inwater_per_m;      /* same for U */
} rt_aw_coupled_field_t;

int  rt_aw_coupled_field_alloc(rt_aw_coupled_field_t *f, int m_max, int n_mu_water);
void rt_aw_coupled_field_free (rt_aw_coupled_field_t *f);

/* Main coupling function.
 *
 * Inputs:
 *   boa_export   — atmospheric BOA downward Stokes field (from
 *                  rt_solve_case_pol_for_ocean)
 *   n_water      — water refractive index (real, ≈ 1.34)
 *   q_convention — 0 legacy or 1 Mishchenko (passed to T_aw_radiance)
 *   mu_water_pos — [n_mu_water] positive water-side μ grid (caller-provided;
 *                  typically Gauss-Legendre pos nodes used by in-water SOS)
 *   n_mu_water   — number of water-side positive μ nodes
 *
 * Output:
 *   coupled      — pre-allocated with matching m_max and n_mu_water;
 *                  fields filled with in-water inbound Stokes per m per μ_w.
 *
 * Returns 0 on success, negative on failure.
 */
int rt_air_water_couple_atm_to_water(const rt_atm_boa_export_t *boa_export,
                                       double n_water,
                                       int q_convention,
                                       const double *mu_water_pos,
                                       int n_mu_water,
                                       const double *w_atm_pos,
                                       double wind_speed,
                                       int sigma_type,
                                       rt_aw_coupled_field_t *coupled);

/* Reverse coupling: in-water z=0⁻ upwelling Fourier field → above-water z=0⁺
 * field on the atmosphere μ-node grid (atm SOS bottom source). FLAT interface
 * (per-mode specular + T_wa); rough (Cox-Munk water→air BTDF) = Milestone 2b.
 * Output buffers {I,Q,U}_air_per_m are CALLER-allocated, each size
 * (m_max+1)*n_mu_atm, indexed [m*n_mu_atm + j_a]. See the .c for full notes. */
int rt_air_water_couple_water_to_atm(const double *I_water_per_m,
                                       const double *Q_water_per_m,
                                       const double *U_water_per_m,
                                       const double *mu_water_pos,
                                       const double *w_water_pos,
                                       int n_mu_water, int m_max,
                                       double n_water, int q_convention,
                                       double wind_speed, int sigma_type,
                                       const double *mu_atm_pos, int n_mu_atm,
                                       const double *w_atm_pos,
                                       double *I_air_per_m,
                                       double *Q_air_per_m,
                                       double *U_air_per_m);

void mb_re_escape_add_m0(const double *I_water_per_m,
                         const double *mu_water_pos, const double *w_water,
                         int n_mu_water, double n_water, double wind_speed,
                         const double *mu_atm_pos, int n_mu_atm,
                         const double *w_atm_pos, double *I_air_per_m);

#endif /* OCRT_V2_RT_AIR_WATER_COUPLING_H */
