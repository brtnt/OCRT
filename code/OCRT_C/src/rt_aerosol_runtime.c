#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rt_aerosol_runtime.h"
#include "rt_aerosol.h"
#include "shared/numerics.h"

static void free_phase(double *a, double *b, double *c) {
    free(a); free(b); free(c);
}

int rt_aerosol_runtime_prepare(const mie_data_t *mie,
                               double wavelength_nm,
                               double user_aod_ref,
                               double aod_ref_nm,
                               const rt_aerosol_runtime_options_t *ropts,
                               rt_aerosol_input_t *aer,
                               rt_aerosol_runtime_diag_t *diag) {
    if (!mie || !ropts || !aer) return -1;
    if (mie->n_wl <= 0 || mie->n_ang <= 0 || user_aod_ref <= 0.0 ||
        aod_ref_nm <= 0.0 ||
        !mie_data_covers_wavelength(mie, wavelength_nm / 1000.0) ||
        !mie_data_covers_bulk_wavelength(mie, aod_ref_nm / 1000.0)) return -1;
    if (ropts->aer_L_max < 2) return -1;

    memset(aer, 0, sizeof(*aer));
    if (diag) memset(diag, 0, sizeof(*diag));

    const double wl_um = wavelength_nm / 1000.0;
    int wl_row_n = (mie->n_phase_wl > mie->n_wl) ? mie->n_phase_wl : mie->n_wl;
    double *wl_row = (double*)malloc((size_t)wl_row_n * sizeof(double));
    double *P11_w  = (double*)malloc((size_t)mie->n_ang * sizeof(double));
    double *P12_w  = (double*)malloc((size_t)mie->n_ang * sizeof(double));
    double *P33_w  = (double*)malloc((size_t)mie->n_ang * sizeof(double));
    if (!wl_row || !P11_w || !P12_w || !P33_w) {
        free(wl_row); free_phase(P11_w, P12_w, P33_w);
        return -1;
    }

    /* User AOD is specified at a fixed reference wavelength. Convert it to
     * the current calculation wavelength with the .mie spectral extinction
     * ratio. Column 0 is Nor_Ext_Co; absolute Extinct_Co would give the same
     * ratio. This keeps one physical aerosol column across a spectral run. */
    for (int iw = 0; iw < mie->n_wl; ++iw) wl_row[iw] = SPECTRAL_AT(mie, iw, 0);
    const double ext_target = linterp(mie->wavelengths, wl_row, mie->n_wl, wl_um);
    const double ext_ref = linterp(mie->wavelengths, wl_row, mie->n_wl, aod_ref_nm / 1000.0);
    if (!(ext_target > 0.0 && ext_ref > 0.0)) {
        free(wl_row); free_phase(P11_w, P12_w, P33_w);
        return -1;
    }
    const double extinction_ratio = ext_target / ext_ref;
    const double user_aod = user_aod_ref * extinction_ratio;

    for (int iw = 0; iw < mie->n_wl; ++iw) wl_row[iw] = SPECTRAL_AT(mie, iw, 2);
    const double ssa_a = linterp(mie->wavelengths, wl_row, mie->n_wl, wl_um);

    /* Phase-matrix wavelength interpolation follows OSOAA conv_local:
     * shape-preserving cubic PCHIP for P11/P12/P33.  AOD extinction and SSA
     * retain the pre-existing linear spectral interpolation. */
    for (int ia = 0; ia < mie->n_ang; ++ia) {
        for (int iw = 0; iw < mie->n_phase_wl; ++iw) wl_row[iw] = BLOCK_AT(mie->P11, mie->n_phase_wl, ia, iw);
        P11_w[ia] = pchip_interp(mie->phase_wavelengths, wl_row, mie->n_phase_wl, wl_um);
        for (int iw = 0; iw < mie->n_phase_wl; ++iw) wl_row[iw] = BLOCK_AT(mie->P12, mie->n_phase_wl, ia, iw);
        P12_w[ia] = pchip_interp(mie->phase_wavelengths, wl_row, mie->n_phase_wl, wl_um);
        for (int iw = 0; iw < mie->n_phase_wl; ++iw) wl_row[iw] = BLOCK_AT(mie->P33, mie->n_phase_wl, ia, iw);
        P33_w[ia] = pchip_interp(mie->phase_wavelengths, wl_row, mie->n_phase_wl, wl_um);
    }
    free(wl_row); wl_row = NULL;

    double f_fwd = 0.0;
    if (ropts->use_loglin_trunc) {
        /* Real-angle log10-linear forward-peak truncation.
         * Mass-conserving variant: layer τ_sca shrinks by (1 − A/2), and
         * the τ_ext / SSA renormalization recovers the actual optical
         * properties of the truncated phase function (see formulas below). */
        double *P11_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        double *P12_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        double *P33_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        if (!P11_t || !P12_t || !P33_t) {
            free_phase(P11_w, P12_w, P33_w);
            free_phase(P11_t, P12_t, P33_t);
            return -1;
        }
        double mu1 = (ropts->loglin_mu1 > 0.0) ? ropts->loglin_mu1 : 0.8;
        double mu2 = (ropts->loglin_mu2 > 0.0) ? ropts->loglin_mu2 : 0.94;
        double thresh = (ropts->loglin_threshold > 0.0) ? ropts->loglin_threshold : 0.1;
        if (rt_aerosol_loglin_truncate(P11_w, P12_w, P33_w,
                                       mie->angles, mie->n_ang,
                                       mu1, mu2, thresh,
                                       P11_t, P12_t, P33_t, &f_fwd) != 0) {
            free_phase(P11_w, P12_w, P33_w);
            free_phase(P11_t, P12_t, P33_t);
            return -1;
        }
        free_phase(P11_w, P12_w, P33_w);
        P11_w = P11_t; P12_w = P12_t; P33_w = P33_t;
        /* NOTE: f_fwd here is the truncation coefficient A (range 0–2),
         * not the Wiscombe forward fraction. τ/SSA renormalization below
         * uses the matching mass-conserving formulas. */
    } else if (ropts->delta_m_N > 0) {
        double *P11_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        double *P12_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        double *P33_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        if (!P11_t || !P12_t || !P33_t) {
            free_phase(P11_w, P12_w, P33_w);
            free_phase(P11_t, P12_t, P33_t);
            return -1;
        }
        if (rt_aerosol_delta_m_truncate(P11_w, P12_w, P33_w,
                                        mie->angles, mie->n_ang,
                                        ropts->delta_m_N,
                                        P11_t, P12_t, P33_t, &f_fwd) != 0) {
            free_phase(P11_w, P12_w, P33_w);
            free_phase(P11_t, P12_t, P33_t);
            return -1;
        }
        free_phase(P11_w, P12_w, P33_w);
        P11_w = P11_t; P12_w = P12_t; P33_w = P33_t;
    } else if (ropts->theta_cut_deg > 0.0) {
        double *P11_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        double *P12_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        double *P33_t = (double*)malloc((size_t)mie->n_ang * sizeof(double));
        if (!P11_t || !P12_t || !P33_t) {
            free_phase(P11_w, P12_w, P33_w);
            free_phase(P11_t, P12_t, P33_t);
            return -1;
        }
        if (rt_aerosol_delta_fit_truncate(P11_w, P12_w, P33_w,
                                          mie->angles, mie->n_ang,
                                          ropts->theta_cut_deg,
                                          P11_t, P12_t, P33_t, &f_fwd) != 0) {
            free_phase(P11_w, P12_w, P33_w);
            free_phase(P11_t, P12_t, P33_t);
            return -1;
        }
        free_phase(P11_w, P12_w, P33_w);
        P11_w = P11_t; P12_w = P12_t; P33_w = P33_t;
    }

    double tau_a_eff = user_aod;
    double ssa_a_eff = ssa_a;
    if (ropts->use_loglin_trunc && f_fwd > 0.0) {
        /* Mass-conserving τ/SSA renormalization for real-angle log10-linear
         * truncation. f_fwd here is the truncation coefficient A.
         *
         *   τ_sca' = τ_sca · (1 − A/2)            (scattering depth shrinks)
         *   τ_ext' = τ_ext · (1 − ω · A/2)        (extinction follows)
         *   ω'     = (1 − A/2)·ω / (1 − ω·A/2)    (SSA adjustment)
         *
         * Derivation: starting from τ_sca' = ω·τ_ext · (1 − A/2) and
         * τ_ext' = τ_sca' / ω', one gets the closed-form ω' above. */
        double A = f_fwd;
        tau_a_eff = user_aod * (1.0 - 0.5 * ssa_a * A);
        ssa_a_eff = (1.0 - 0.5 * A) * ssa_a / (1.0 - 0.5 * ssa_a * A);
    } else if (ropts->apply_nt_tau && f_fwd > 0.0) {
        /* Wiscombe/NT-style similarity transform for delta-M / delta-fit
         * truncation (f_fwd = Wiscombe forward fraction f in [0,1)).
         *
         * *** SUSPECTED BUG (v1.09 annotation pass, 2026-07-04;
         *     BUG_SUSPECTS.md S-001) ***
         * The standard Wiscombe (1977) similarity relations are
         *   tau' = (1 - omega * f) * tau
         *   omega' = (1 - f) * omega / (1 - omega * f)
         * The omega' line below IS the standard form, but the tau' line is
         *   tau' = (1 - f/2) * tau        <- factor 0.5, and NO omega
         * which is neither Wiscombe nor self-consistent with the omega'
         * line: the implied scattering depth omega'*tau' equals the correct
         * (1-f)*omega*tau only at omega = 1/2.  Compare the loglin branch
         * above, which implements the correct pattern (with f = A/2).
         * EXPOSURE: only the optional --nt-tau flag (default off); the
         * consistency-harness goldens never enable it, and production
         * aerosol runs are frozen by the value-kernel open defect.
         * FIXED in v1.09 commit #4 (2026-07-05).  Validation record
         * (BUG_SUSPECTS S-001 items 1-9): absorption-invariance analytic
         * proof; M80C internal A/B vs untruncated reference (error halved
         * to a third); OSOAA truncation-invariance ratio (1.36-1.45 -> 
         * 1.13-1.19 against the 0.997 OSOAA baseline); matched absolute
         * point +20.0% -> +0.3% (vza30).  The whole delta-fit/delta-M +
         * nt-tau family is DEBUG-gated at the CLI (non-validated research
         * path); the remaining ~+-10% truncation-shape residual is NOT a
         * bookkeeping error and is out of scope here.
         * Pre-fix formula, for archaeology:  tau' = tau * (1 - f/2). */
        tau_a_eff = user_aod * (1.0 - f_fwd * ssa_a);   /* Wiscombe tau' */
        ssa_a_eff = (1.0 - f_fwd) * ssa_a / (1.0 - f_fwd * ssa_a);
    }

    double *betal  = (double*)calloc((size_t)(ropts->aer_L_max + 1), sizeof(double));
    double *gammal = (double*)calloc((size_t)(ropts->aer_L_max + 1), sizeof(double));
    double *alphal = (double*)calloc((size_t)(ropts->aer_L_max + 1), sizeof(double));
    double *zetal  = (double*)calloc((size_t)(ropts->aer_L_max + 1), sizeof(double));
    if (!betal || !gammal || !alphal || !zetal) {
        free_phase(P11_w, P12_w, P33_w);
        free(betal); free(gammal); free(alphal); free(zetal);
        return -1;
    }

    if (rt_aerosol_compute_vector_legendre(P11_w, P12_w, P33_w,
                                           mie->angles, mie->n_ang,
                                           ropts->aer_L_max, 4096,
                                           betal, gammal, alphal, zetal) != 0) {
        free_phase(P11_w, P12_w, P33_w);
        free(betal); free(gammal); free(alphal); free(zetal);
        return -1;
    }

    /* Normalize the angle-space phase for the value kernel by
     * beta0 = (1/2) int P11 dmu.  The moment kernel divides ALL Legendre
     * moments (alpha,beta,gamma,zeta) by beta0 = betal[0] (rt_aerosol.c:684-690),
     * so its effective phase integrates to 1.  The value kernel reconstructs the
     * raw truncated phase directly, which integrates to beta0 (=0.758 for
     * loglin-truncated M80C@412), making the value-kernel source term uniformly
     * too small by exactly that factor.  P11/P12/P33_w are absolute matrix
     * elements (not ratios), so all three scale by 1/beta0 together -- this is the
     * same delta-M renormalization the moment path applies.  Use the IDENTICAL
     * dense-PCHIP fine_grid_points (4096) as the moment call above so beta0 is
     * bit-consistent with betal[0]. */
    {
        double chi0[1];
        if (rt_aerosol_compute_legendre_moments(P11_w, mie->angles, mie->n_ang,
                                                1, 4096, chi0) == 0
            && chi0[0] != 0.0) {
            double inv_b0 = 1.0 / chi0[0];
            for (int i = 0; i < mie->n_ang; ++i) {
                P11_w[i] *= inv_b0;
                P12_w[i] *= inv_b0;
                P33_w[i] *= inv_b0;
            }
        }
    }

    /* Keep the (capped) angle-space phase for the value kernel; do NOT free.
     * theta_phase is copied (mie->angles is owned by mie); P11/P12/P33_phase
     * take ownership of the (possibly truncated) P11_w/P12_w/P33_w buffers,
     * which are released in rt_aerosol_runtime_free.
     *
     * CRITICAL: .mie angle grids are DESCENDING (180 -> 0 deg, OPAC convention),
     * but aer_phase_p11_at_costheta() in rt_solver.c does an ascending binary
     * search.  Store theta + P11/P12/P33 in ASCENDING order (reverse when the
     * input is descending) so the lookup is correct.  Without this the lookup
     * short-circuits to P[0] = P11(180 deg) for every scattering angle, which
     * replaces the aerosol phase with an isotropic constant. */
    {
        int na = mie->n_ang;
        int descending = (na >= 2 && mie->angles[0] > mie->angles[na - 1]);
        double *th = (double*)malloc((size_t)na * sizeof(double));
        if (th) {
            for (int i = 0; i < na; ++i)
                th[i] = descending ? mie->angles[na - 1 - i] : mie->angles[i];
        }
        if (descending) {
            for (int i = 0; i < na / 2; ++i) {
                int r = na - 1 - i; double t;
                t = P11_w[i]; P11_w[i] = P11_w[r]; P11_w[r] = t;
                t = P12_w[i]; P12_w[i] = P12_w[r]; P12_w[r] = t;
                t = P33_w[i]; P33_w[i] = P33_w[r]; P33_w[r] = t;
            }
        }
        aer->n_ang_phase = na;
        aer->theta_phase = th;
        aer->P11_phase = P11_w;
        aer->P12_phase = P12_w;
        aer->P33_phase = P33_w;
    }
    aer->use_value_kernel = ropts->use_value_kernel;

    aer->tau_a = tau_a_eff;
    aer->ssa_a = ssa_a_eff;
    aer->L_max = ropts->aer_L_max;
    aer->betal_aer = betal;
    aer->gammal_aer = gammal;
    aer->alphal_aer = alphal;
    aer->zetal_aer = zetal;

    if (diag) {
        diag->f_fwd = f_fwd;
        diag->ssa_raw = ssa_a;
        diag->tau_a_eff = tau_a_eff;
        diag->ssa_a_eff = ssa_a_eff;
        diag->aod_ref = user_aod_ref;
        diag->aod_ref_nm = aod_ref_nm;
        diag->aod_target = user_aod;
        diag->extinction_ratio = extinction_ratio;
    }
    return 0;
}

void rt_aerosol_runtime_free(rt_aerosol_input_t *aer) {
    if (!aer) return;
    free((void*)aer->betal_aer);
    free((void*)aer->gammal_aer);
    free((void*)aer->alphal_aer);
    free((void*)aer->zetal_aer);
    free((void*)aer->theta_phase);
    free((void*)aer->P11_phase);
    free((void*)aer->P12_phase);
    free((void*)aer->P33_phase);
    memset(aer, 0, sizeof(*aer));
}
