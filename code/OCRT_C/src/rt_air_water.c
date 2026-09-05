/* ============================================================================
 * rt_air_water.c — Air-water interface Mueller transmission implementation
 *
 * Phase B.2 (2026-05-22).
 *
 * Thin wrappers around shared/surface.h generalized Fresnel API. The Fresnel
 * core (Snell + s/p amplitudes) is implemented file-static in surface.c and
 * shared between reflection and transmission paths — this guarantees
 * numerical consistency.
 * ============================================================================ */

#include "rt_air_water.h"
#include "shared/surface.h"
#include "shared/mat3.h"
#include <math.h>

/* ---------------------------------------------------------------------------
 * Snell helpers
 * --------------------------------------------------------------------------- */

double rt_air_water_mu_critical(double n_water) {
    /* sin(theta_c) = 1/n_water  =>  mu_c = sqrt(1 - 1/n_water²) */
    double inv_n_sq = 1.0 / (n_water * n_water);
    if (inv_n_sq >= 1.0) return 0.0;
    return sqrt(1.0 - inv_n_sq);
}

int rt_air_water_is_TIR(double mu_w, double n_water) {
    return (mu_w < rt_air_water_mu_critical(n_water)) ? 1 : 0;
}

double rt_air_water_mu_refracted_down(double mu_a, double n_water) {
    /* Snell: sin(theta_w) = sin(theta_a)/n_water, n_water > 1 → no TIR. */
    double sin_a_sq = fmax(0.0, 1.0 - mu_a * mu_a);
    double sin_w_sq = sin_a_sq / (n_water * n_water);
    if (sin_w_sq >= 1.0) return 0.0;  /* numerical edge (n_water > 1 shouldn't hit) */
    return sqrt(1.0 - sin_w_sq);
}

double rt_air_water_mu_refracted_up(double mu_w, double n_water) {
    /* Snell: sin(theta_a) = n_water · sin(theta_w). TIR if > 1. */
    double sin_w_sq = fmax(0.0, 1.0 - mu_w * mu_w);
    double sin_a_sq = n_water * n_water * sin_w_sq;
    if (sin_a_sq >= 1.0) return -1.0;  /* TIR signal */
    return sqrt(1.0 - sin_a_sq);
}

/* ---------------------------------------------------------------------------
 * Mueller matrix wrappers
 * --------------------------------------------------------------------------- */

void rt_air_water_T_aw(double mu_a, double n_water,
                       int q_convention, double *M_T) {
    /* (n1, n2) = (1, n_water): downward refraction across air-water surface */
    surface_flat_fresnel_T_matrix(mu_a, 1.0, n_water, q_convention, M_T);
}

void rt_air_water_T_wa(double mu_w, double n_water,
                       int q_convention, double *M_T) {
    /* (n1, n2) = (n_water, 1): upward refraction; TIR returns zero matrix */
    surface_flat_fresnel_T_matrix(mu_w, n_water, 1.0, q_convention, M_T);
}

void rt_air_water_R_aa(double mu_a, double n_water,
                       int q_convention, double *M_R) {
    /* Identical to surface_flat_fresnel_matrix_raw (legacy alias) */
    surface_flat_fresnel_R_matrix_general(mu_a, 1.0, n_water, q_convention, M_R);
}

void rt_air_water_R_ww(double mu_w, double n_water,
                       int q_convention, double *M_R) {
    /* (n1, n2) = (n_water, 1): water-side reflection; TIR → identity */
    surface_flat_fresnel_R_matrix_general(mu_w, n_water, 1.0, q_convention, M_R);
}

double rt_air_water_check_energy(double mu_i, double n1, double n2) {
    return surface_flat_fresnel_check_energy(mu_i, n1, n2);
}
