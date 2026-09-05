#ifndef OCRT_V2_RT_AEROSOL_RUNTIME_H
#define OCRT_V2_RT_AEROSOL_RUNTIME_H

#include "rt_solver.h"
#include "shared/mie_io.h"

/* Runtime adapter from .mie tables to rt_aerosol_input_t.
 *
 * This module centralizes the aerosol path used by both single-case and
 * batch execution:
 *   1) wavelength interpolation of SSA/P11/P12/P33,
 *   2) optional forward-peak truncation,
 *   3) optional NT optical-depth transform,
 *   4) vector Legendre moment construction.
 */
typedef struct {
    int    aer_L_max;        /* Legendre max order, e.g. 80 */
    double theta_cut_deg;    /* delta-fit cutoff; 0 disables */
    int    delta_m_N;        /* delta-M N; 0 disables */
    int    apply_nt_tau;     /* apply Wiscombe NT τ/SSA transform if f_fwd > 0 */
    int    use_loglin_trunc; /* real-angle log10-linear forward-peak truncation
                              * + matching τ/SSA renormalization (Potter 1970-style
                              * truncation, mass-conserving variant) */
    double loglin_mu1;       /* anchor cos(θ_1); default 0.8 → 36.87° */
    double loglin_mu2;       /* anchor cos(θ_2); default 0.94 → 19.95° */
    double loglin_threshold; /* skip truncation if A < this; default 0.1 */
    int    use_value_kernel; /* 1 = build SOS phase_fourier_m from angle-space
                              * (loglin-capped) phase VALUES (Gibbs-free), not from
                              * Legendre moments. Default policy set in main.c. */
} rt_aerosol_runtime_options_t;

typedef struct {
    double f_fwd;
    double ssa_raw;
    double tau_a_eff;
    double ssa_a_eff;
    double aod_ref;
    double aod_ref_nm;
    double aod_target;
    double extinction_ratio;
} rt_aerosol_runtime_diag_t;

int rt_aerosol_runtime_prepare(const mie_data_t *mie,
                               double wavelength_nm,
                               double user_aod_ref,
                               double aod_ref_nm,
                               const rt_aerosol_runtime_options_t *ropts,
                               rt_aerosol_input_t *aer,
                               rt_aerosol_runtime_diag_t *diag);

void rt_aerosol_runtime_free(rt_aerosol_input_t *aer);

#endif /* OCRT_V2_RT_AEROSOL_RUNTIME_H */
