/* ============================================================================
 * rt_air_water.h — Air-water interface vector Mueller transmission API
 *
 * Phase B.2 implementation (2026-05-22).
 *
 * 본 모듈은 atmosphere-ocean coupled vector RT의 *interface coupling 층*을
 * 제공한다. shared/surface.h 의 generalized Fresnel API
 * (surface_flat_fresnel_R_matrix_general, surface_flat_fresnel_T_matrix)를
 * thin wrapper로 호출하여, OCRT의 in-water solver와 atmospheric solver
 * 사이의 boundary coupling을 단순화한다.
 *
 * Symmetry with existing surface code:
 *   기존 surface_flat_fresnel_matrix_raw 는 air-side reflection at flat water
 *   surface (n1=1, n2=n_water) 만 처리하였다. 본 모듈은 같은 Fresnel core를
 *   재사용하여 4가지 방향 (R/T, air/water side)를 모두 처리한다.
 *
 * Reference:
 *   - Mishchenko, Travis & Lacis 2002, Cambridge UP
 *   - Hovenier, van der Mee & Domke 2004, Kluwer
 *   - Mobley 1994, "Light and Water", Ch. 4.5
 * ============================================================================ */

#ifndef OCRT_RT_AIR_WATER_H
#define OCRT_RT_AIR_WATER_H

#ifdef __cplusplus
extern "C" {
#endif

/* ===========================================================================
 * Snell refraction helpers
 * =========================================================================== */

/* Refraction: air → water
 *   mu_a    : abs cos of incidence angle in air (>= 0)
 *   n_water : water refractive index
 * Returns mu_w = sqrt(1 - sin²(theta_w)), where sin(theta_w) = sin(theta_a)/n_w
 * No TIR in this direction (n_air < n_water). Result in (0, 1].
 */
double rt_air_water_mu_refracted_down(double mu_a, double n_water);

/* Refraction: water → air
 *   mu_w    : abs cos of incidence angle in water (>= 0)
 *   n_water : water refractive index
 * Returns mu_a if no TIR, else -1.0.
 * TIR occurs when sin(theta_w) > 1/n_water (i.e. mu_w < mu_critical).
 */
double rt_air_water_mu_refracted_up(double mu_w, double n_water);

/* Critical-angle cosine in water medium.
 *   Returns mu_w_critical = sqrt(1 - 1/n_water²)
 *   For n_water = 1.34: mu_w_critical ≈ 0.6614, theta_critical ≈ 48.61°
 */
double rt_air_water_mu_critical(double n_water);

/* TIR check: returns 1 if total internal reflection at water-air interface
 *   for given water-side incidence cosine, else 0.
 */
int rt_air_water_is_TIR(double mu_w, double n_water);

/* ===========================================================================
 * Air ↔ Water Mueller matrices (3x3, meridian basis)
 *
 * These are thin wrappers around shared/surface.h generalized API. The
 * underlying Fresnel core is identical to that used by the existing
 * surface_flat_fresnel_matrix_raw (air-side reflection), guaranteeing
 * numerical consistency between reflection and transmission paths.
 * =========================================================================== */

/* Air → Water transmission Mueller matrix (radiance form)
 *   mu_a         : abs cos of incidence angle in air (>= 0)
 *   n_water      : water refractive index
 *   q_convention : 0 legacy, 1 Mishchenko (matches OCRT default)
 *   M_T[9]       : output, row-major 3x3
 * Includes:
 *   - Snell refraction (internal)
 *   - Fresnel transmission Mueller (s/p combined)
 *   - Radiance n² law:  T_radiance = (n_w/1)² · T_flux
 * No TIR possible (air → water always transmits).
 */
void rt_air_water_T_aw(double mu_a, double n_water,
                       int q_convention, double *M_T);

/* Water → Air transmission Mueller matrix (radiance form)
 *   mu_w         : abs cos of incidence angle in water (>= 0)
 *   n_water      : water refractive index
 *   q_convention : 0 legacy, 1 Mishchenko
 *   M_T[9]       : output, row-major 3x3
 * TIR branch (mu_w < mu_critical): returns zero matrix.
 * Otherwise applies Snell + Fresnel transmission with n² law for upward
 * transmission:  T_radiance = (1/n_w)² · T_flux  (n_ratio < 1).
 */
void rt_air_water_T_wa(double mu_w, double n_water,
                       int q_convention, double *M_T);

/* Air-side reflection Mueller (existing surface_flat_fresnel_matrix_raw alias)
 *   Provided here for API completeness; identical to surface_flat_fresnel_matrix_raw.
 */
void rt_air_water_R_aa(double mu_a, double n_water,
                       int q_convention, double *M_R);

/* Water-side internal reflection Mueller matrix
 *   mu_w         : abs cos of incidence angle in water (>= 0)
 *   n_water      : water refractive index
 *   q_convention : 0 legacy, 1 Mishchenko
 *   M_R[9]       : output, row-major 3x3
 * Includes TIR branch: returns identity matrix when mu_w < mu_critical.
 * Non-TIR branch returns standard Fresnel reflection with (n1, n2) =
 * (n_water, 1).
 *
 * Limitation (B.2): TIR phase shift between s and p amplitudes
 * (which would populate off-diagonal Mueller elements for polarized light)
 * is not modelled in this implementation; identity is used as a first-order
 * approximation. For pure water Rayleigh-like scattering this effect is
 * small. Revisit if validation against OSOAA shows discrepancy at high
 * water-side incidence angles.
 */
void rt_air_water_R_ww(double mu_w, double n_water,
                       int q_convention, double *M_R);

/* ===========================================================================
 * Diagnostics
 * =========================================================================== */

/* Energy conservation check at flat interface.
 *   For (mu_i, n1, n2), verify:
 *     |1 - (R_X_flux + T_X_flux)| summed over X in {s, p}
 *   Returns residual; expect < 1e-12 for non-TIR.
 *   Returns 0.0 for TIR (R=1 by definition, no transmission).
 *
 *   Wrapper around surface_flat_fresnel_check_energy.
 */
double rt_air_water_check_energy(double mu_i, double n1, double n2);

#ifdef __cplusplus
}
#endif

#endif /* OCRT_RT_AIR_WATER_H */
