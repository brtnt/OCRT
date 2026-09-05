#ifndef OCRT_V2_RT_IO_H
#define OCRT_V2_RT_IO_H

#include <stdio.h>
#include "rt_types.h"

/* Optional aerosol batch configuration.
 * When mie_path is non-NULL and the CSV contains an aod_555 or aod_865 column,
 * rt_io_run_batch computes wavelength-interpolated aerosol Legendre moments
 * per row and calls rt_solve_case_pol_aerosol().
 */
typedef struct {
    const char *mie_path;
    int         aer_L_max;
    double      theta_cut_deg;
    int         delta_m_N;
    int         apply_nt_tau;
    int         use_loglin_trunc;
} rt_batch_aerosol_config_t;

/* Batch CSV driver (Step 10 + Step 7 vector extension).
 *
 * Input CSV (header + rows):
 *   case_id, sza_deg, vza_deg, raa_deg, wavelength_nm[, tau_R_ref, rho_I_ref]
 *
 * The reference columns are optional.  Header column count determines
 * which columns are present (5 = no reference, 6 = + tau_R_ref only,
 * 7 = + tau_R_ref + rho_I_ref).
 *
 * Output CSV — scalar mode (opts->vector_mode = 0):
 *   case_id, sza_deg, vza_deg, raa_deg, wavelength_nm,
 *   rho_I_v2, rho_I_ref, delta_abs, delta_rel_pct, delta_rel_abs_pct,
 *   tau_R_v2, tau_R_ref, n_orders_m0, n_orders_m1, n_orders_m2,
 *   walltime_ms
 *
 * Output CSV — vector mode (opts->vector_mode = 1):
 *   ... (same 16 columns above) ..., rho_Q_v2, rho_U_v2
 *
 *   The Q/U columns are APPENDED at the end so the column-index
 *   positions of the legacy 16 columns remain unchanged.  Downstream
 *   tools that parse rho_I_v2 by column index keep working.  Reference
 *   Q/U columns are not yet wired (Step 8 — 1200-case validation will
 *   add `rho_Q_ref`, `rho_U_ref` columns in a separate path).
 *
 * delta_rel_pct  = 100 * (rho_v2 - rho_ref) / rho_ref     (signed)
 * delta_rel_abs_pct = |delta_rel_pct|                     (PASS gauge)
 *
 * Per-row execution is parallelized via OpenMP; output rows are
 * collected into an in-memory array first and written in input order
 * after all cases finish.
 *
 * Solver dispatch: opts->vector_mode == 0 calls rt_solve_case (scalar);
 * == 1 calls rt_solve_case_pol (vector I/Q/U).  Both share the same
 * tau_R, sos_orders_per_m, walltime_ms diagnostics.
 *
 * Returns 0 on completion (independent of PASS rate -- caller decides),
 * -1 on file I/O error or malformed CSV header, -2 if any individual
 * solver call returns -2 (NaN propagation).
 *
 * opts pass NULL to use defaults (n_mu=24, n_layers=40, m_max=2,
 * Bodhaine, LINEAR, vector_mode=0).
 */
/* cs_template (NULL → black surface, no aerosol, Rayleigh-only — original
 * default): provides case-level surface BC, wind speed, water refractive
 * index, Cox-Munk sigma type, q convention, sunglint decoupling flag.
 * These are applied uniformly to every row of the batch CSV; per-row
 * geometry/wavelength/tau_R_ref still come from the CSV columns. */
int rt_io_run_batch(const char *in_csv, const char *out_csv,
                    const rt_options_t *opts,
                    const rt_case_t *cs_template,
                    const rt_batch_aerosol_config_t *aer_cfg);

/* Print a single result line to stdout in stable format. */
void rt_io_print_single(const rt_case_t *cs, const rt_result_t *res);

#endif /* OCRT_V2_RT_IO_H */
