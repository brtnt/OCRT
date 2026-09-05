/* =============================================================================
 * rt_aerosol.c — Phase-matrix Greek-moment production + forward-peak truncation
 * =============================================================================
 *
 * ROLE (producer side of the kernels in rt_kernel.c)
 *   Converts tabulated Mueller elements P11(theta), P12(theta), P33(theta)
 *   (from tabulated .mie or external phase data) into the Greek expansion moments
 *   consumed by the SOS kernels, and applies forward-peak truncation:
 *
 *     beta_l  : a1 = P11 moments in ordinary Legendre P_l
 *     gamma_l : b1 = P12 moments in the generalized basis P2_l (see below)
 *     alpha_l : a2 = P22 moments   } completed from the ordinary and
 *     zeta_l  : a3 = P33 moments   } spin-2 moment identities
 *     delta_l : a4-side helper = P33 moments in ordinary P_l (internal)
 *
 *   Moment conventions: chi_l ("betal") = (2l+1) * g_l with g_l the
 *   normalized moment (g_0 = 1 after the final beta_0 division); this is
 *   exactly the normalization rt_kernel_phase_fourier() expects.
 *
 * BASIS CONSISTENCY (VERIFIED, closes the note left in rt_kernel.c header)
 *   The P2_l recurrence used for gamma_l below is ALGEBRAICALLY IDENTICAL
 *   to the m = 0 rsl recurrence in rt_kernel.c:
 *     seed:  P2_2(mu) = 3(1-mu^2)/(2 sqrt(6))  ==  rsl(2,mu) at m=0
 *     coeffs: d = (2k+1)/sqrt((k+3)(k-1))
 *             == [(l+1)(2l+1)/sqrt((l+3)(l-1)(l+1)^2)]  (kernel d at m=0)
 *             e = sqrt((k+2)(k-2))/(2k+1)
 *             == [sqrt((l+2)(l-2)) * l / (l(2l+1))]     (kernel e at m=0)
 *     and the kernel's rtl coupling term vanishes at m=0 (f = 2m/... = 0,
 *     rtl(l,mu)|_{m=0} = 0), so the decoupled recurrence here is exact.
 *   Hence gamma_l produced here contracts consistently with the spin-2
 *   angular basis used by the Fourier kernel.  Alpha/zeta are completed by
 *   the published reciprocal spherical-particle moment identities.
 *
 * SPHERICAL-PARTICLE ASSUMPTION
 *   The diagonal-moment completion uses beta_l(P11) where P22 moments are
 *   needed, i.e. it assumes P22 == P11.  For Mie SPHERES this is exact
 *   (a2 = a1 and a4 = a3 identically); for nonspherical particles it is an
 *   approximation.  All current OCRT aerosol/hydrosol phase inputs are Mie.
 *
 * INPUT TABLE CONSTRAINT
 *   theta tables MUST span the full [0, 180] deg range.
 *   rt_aerosol_compute_legendre_moments() integrates on a fixed dense
 *   0..180 grid and pchip_eval EXTRAPOLATES with the edge cubic outside the
 *   table — for a forward-peaked P11 with a table starting above 0 deg that
 *   extrapolation is uncontrolled.  compute_p2_moments() by contrast
 *   integrates only over the tabulated range.  With full-range tables
 *   (verified for all production .mie inputs: 0.000..180.000) the two
 *   integrators are consistent; narrower tables would make beta/delta and
 *   gamma mutually inconsistent — do not feed partial-range tables.
 *   [BUG_SUSPECTS.md S-002]
 *
 * IMPLEMENTATION
 *   OCRT-owned interpolation, quadrature, moment projection and truncation.
 *   Dependencies: shared/numerics.h, rt_quadrature.h, rt_phase_moments.h.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rt_aerosol.h"
#include "shared/numerics.h"
#include "rt_quadrature.h"
#include "rt_phase_moments.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ======================================================================== */
/* Ordinary-Legendre moments of a scalar phase table  chi_l, l = 0..n_max-1.
 *
 *   chi_l = (2l+1) * (1/2) * Integral_{-1}^{+1} P(mu) P_l(mu) dmu
 *
 * Pipeline: (1) sort table ascending in theta; (2) PCHIP interpolant in
 * theta (shape-preserving — no overshoot on the forward peak, unlike a
 * cubic spline); (3) evaluate on a dense UNIFORM-IN-THETA grid over the
 * full 0..180 deg (see file-header constraint about table coverage);
 * (4) trapezoid in mu (non-uniform mu spacing handled by the pairwise
 * rule); (5) Bonnet recurrence (l+1)P_{l+1} = (2l+1) mu P_l - l P_{l-1}
 * carried across the l loop (one array pass per order).
 * chi_0 equals the normalization integral: chi_0 = 1 for a properly
 * normalized phase function ((1/2) int P dmu = 1). */

/* v1.09 commit #5 (BUG_SUSPECTS S-002): input-range guard.
 *
 * WHY: the intensity-moment integrator (compute_legendre_moments) fills a
 * fixed dense 0..180 grid and EXTRAPOLATES outside the table (edge cubic,
 * unbounded on a forward peak), while the polarization-moment integrator
 * (compute_p2_moments) integrates only the tabulated range.  A table not
 * covering [0,180] therefore yields intensity and polarization coefficient
 * sets taken from effectively DIFFERENT phase curves (measured on M80C:
 * a 5-deg clip distorts beta_l by up to -10.7% and gamma_l by +-2%, in
 * different directions).  Rather than harmonizing the integrators (a
 * numerics change touching validated goldens), we REJECT partial tables at
 * every public entry: production files all span [0.000, 180.000] exactly,
 * so this guard is inert for spec-conforming inputs.
 * Tolerance 0.5 deg admits slightly coarse endpoint grids while catching
 * any real clipping. */
static int aerosol_table_full_range_ok(const double *theta_deg, int N,
                                       const char *who) {
    if (!theta_deg || N < 2) return 0;
    double tmin = theta_deg[0], tmax = theta_deg[0];
    for (int i = 1; i < N; i++) {
        if (theta_deg[i] < tmin) tmin = theta_deg[i];
        if (theta_deg[i] > tmax) tmax = theta_deg[i];
    }
    if (tmin > 0.5 || tmax < 179.5) {
        fprintf(stderr, "error: %s: phase table covers [%.3f, %.3f] deg; "
                "full [0,180] required (S-002: intensity/polarization moment "
                "integrators disagree on partial tables; see "
                "BUG_SUSPECTS.md S-002)\n", who, tmin, tmax);
        return 0;
    }
    return 1;
}

int rt_aerosol_compute_legendre_moments(const double* P_values, const double* theta_deg, int N,
                                      int n_max, int fine_grid_points, double* chi_out) {
    if (!aerosol_table_full_range_ok(theta_deg, N, "compute_legendre_moments")) return -1;
    /* Ensure theta ascending */
    int ascending = 1;
    for (int i = 1; i < N; ++i) if (theta_deg[i] < theta_deg[i - 1]) { ascending = 0; break; }
    double* theta_sorted = (double*)malloc(N * sizeof(double));
    double* P_sorted     = (double*)malloc(N * sizeof(double));
    if (!theta_sorted || !P_sorted) { free(theta_sorted); free(P_sorted); return -1; }
    if (ascending) {
        memcpy(theta_sorted, theta_deg, N * sizeof(double));
        memcpy(P_sorted,     P_values,  N * sizeof(double));
    } else {
        /* Simple reverse (common case for OPAC 180->0 ordering) */
        for (int i = 0; i < N; ++i) {
            theta_sorted[i] = theta_deg[N - 1 - i];
            P_sorted[i]     = P_values[N - 1 - i];
        }
        /* Still may not be strictly ascending; insertion sort if needed */
        for (int i = 1; i < N; ++i) {
            if (theta_sorted[i] < theta_sorted[i - 1]) {
                double t = theta_sorted[i], v = P_sorted[i]; int j = i - 1;
                while (j >= 0 && theta_sorted[j] > t) {
                    theta_sorted[j + 1] = theta_sorted[j];
                    P_sorted[j + 1] = P_sorted[j]; j--;
                }
                theta_sorted[j + 1] = t; P_sorted[j + 1] = v;
            }
        }
    }

    pchip_t interp;
    if (pchip_build(&interp, theta_sorted, P_sorted, N) != 0) {
        free(theta_sorted); free(P_sorted); return -1;
    }

    /* Dense theta grid, evaluate P, sort by mu ascending */
    double* mu_asc = (double*)malloc(fine_grid_points * sizeof(double));
    double* P_fine = (double*)malloc(fine_grid_points * sizeof(double));
    if (!mu_asc || !P_fine) {
        pchip_free(&interp); free(theta_sorted); free(P_sorted);
        free(mu_asc); free(P_fine); return -1;
    }
    for (int i = 0; i < fine_grid_points; ++i) {
        /* theta descending in index so mu ascending */
        double theta_i = 180.0 * (fine_grid_points - 1 - i) / (double)(fine_grid_points - 1);
        mu_asc[i] = cos(theta_i * M_PI / 180.0);
        P_fine[i] = pchip_eval(&interp, theta_i);
    }

    /* chi_l = (2l+1)*(1/2)*int P P_l dmu, trapezoid */
    double* Pl_vals = (double*)malloc(fine_grid_points * sizeof(double));
    if (!Pl_vals) { pchip_free(&interp); free(theta_sorted); free(P_sorted);
                    free(mu_asc); free(P_fine); return -1; }

    /* Recurrence-based evaluation: carry P_{l-1}, P_l arrays across l loop.
     * Init so that at start of iteration l, `Pl` holds P_l.
     *   Init: Plm1 = P_{-1} (dummy, zero), Pl = P_0 = 1.
     * Advance after each iteration produces P_{l+1} correctly. */
    double* Plm1 = (double*)malloc(fine_grid_points * sizeof(double));
    double* Pl   = (double*)malloc(fine_grid_points * sizeof(double));
    if (!Plm1 || !Pl) { free(Plm1); free(Pl); free(Pl_vals); pchip_free(&interp);
                        free(theta_sorted); free(P_sorted); free(mu_asc); free(P_fine); return -1; }
    for (int i = 0; i < fine_grid_points; ++i) { Plm1[i] = 0.0; Pl[i] = 1.0; }

    for (int l = 0; l < n_max; ++l) {
        double* integrand = Pl_vals;
        for (int i = 0; i < fine_grid_points; ++i) integrand[i] = P_fine[i] * Pl[i];
        double sum = 0.0;
        for (int i = 1; i < fine_grid_points; ++i) {
            sum += 0.5 * (mu_asc[i] - mu_asc[i - 1]) * (integrand[i] + integrand[i - 1]);
        }
        chi_out[l] = (2.0 * l + 1.0) * 0.5 * sum;

        /* Advance recurrence: P_{l+1} = [(2l+1) x P_l - l P_{l-1}] / (l+1) */
        if (l + 1 < n_max) {
            for (int i = 0; i < fine_grid_points; ++i) {
                double next = ((2.0 * l + 1.0) * mu_asc[i] * Pl[i] - (double)l * Plm1[i])
                              / (double)(l + 1);
                Plm1[i] = Pl[i];
                Pl[i]   = next;
            }
        }
    }

    free(Plm1); free(Pl); free(Pl_vals);
    pchip_free(&interp);
    free(theta_sorted); free(P_sorted); free(mu_asc); free(P_fine);
    return 0;
}

/* ======================================================================== */
/* Roof-cap ("delta-fit" family, Hu & Stamnes 2000 style) truncation.
 *
 * Construction, for a cut angle theta_c:
 *   P_cap        = P11(theta_c)            (linear interp AT theta_c;
 *                                           v1.02 grid-independence fix)
 *   P11_capped   = P_cap for theta < theta_c, P11 elsewhere (flat roof)
 *   f            = (1/2) Int_0^{theta_c} (P11 - P_cap) sin(theta) dtheta
 *                = fraction of phase-function mass removed by the roof
 *                  (the boundary term at theta_c vanishes by construction)
 *   P11_out      = P11_capped / (1 - f)    (renormalize to (1/2)IntP dmu = 1)
 *   P12/P33      : capped at the SAME theta_c with their own interpolated
 *                  cap values, same 1/(1-f) renormalization.
 * Caller compensates optics via the standard similarity relations
 * (tau' = (1 - w*f) tau, w' = (1-f) w / (1 - w*f)) in rt_aerosol_runtime.
 * Degenerate cases (cut outside table, f <= 0 or >= 0.999) pass inputs
 * through unchanged with f = 0. */
int rt_aerosol_delta_fit_truncate(const double* P11_in, const double* P12_in, const double* P33_in,
                                const double* theta_deg, int N,
                                double theta_cut_deg,
                                double* P11_out, double* P12_out, double* P33_out,
                                double* f_out) {
    /* Sort ascending if needed */
    int ascending = 1;
    for (int i = 1; i < N; ++i) if (theta_deg[i] < theta_deg[i - 1]) { ascending = 0; break; }

    double* theta_s = (double*)malloc(N * sizeof(double));
    double* P11_s   = (double*)malloc(N * sizeof(double));
    double* P12_s   = (double*)malloc(N * sizeof(double));
    double* P33_s   = (double*)malloc(N * sizeof(double));
    int* idx = (int*)malloc(N * sizeof(int));
    if (!theta_s || !P11_s || !P12_s || !P33_s || !idx) {
        free(theta_s); free(P11_s); free(P12_s); free(P33_s); free(idx);
        return -1;
    }
    if (ascending) {
        for (int i = 0; i < N; ++i) {
            theta_s[i] = theta_deg[i]; idx[i] = i;
            P11_s[i] = P11_in[i]; P12_s[i] = P12_in[i]; P33_s[i] = P33_in[i];
        }
    } else {
        for (int i = 0; i < N; ++i) {
            int j = N - 1 - i;
            theta_s[i] = theta_deg[j]; idx[i] = j;
            P11_s[i] = P11_in[j]; P12_s[i] = P12_in[j]; P33_s[i] = P33_in[j];
        }
    }

    /* Find bracket: theta_s[idx_lo] <= theta_cut <= theta_s[idx_hi] (idx_hi = idx_lo+1).
     * v1.02 fix: interpolate P_cap AT theta_cut (linear), NOT use tabulated row.
     * Previous version (nearest tabulated row) introduced grid dependence. */
    int idx_hi = 0;
    while (idx_hi < N && theta_s[idx_hi] < theta_cut_deg) idx_hi++;
    int idx_cut = idx_hi;  /* keep variable name for downstream */
    if (idx_cut >= N - 1 || idx_cut < 2) {
        /* No meaningful truncation -- copy inputs to outputs */
        for (int i = 0; i < N; ++i) {
            P11_out[i] = P11_in[i]; P12_out[i] = P12_in[i]; P33_out[i] = P33_in[i];
        }
        *f_out = 0.0;
        free(theta_s); free(P11_s); free(P12_s); free(P33_s); free(idx);
        return 0;
    }

    /* Linear interpolation of P at exact theta_cut */
    int idx_lo = idx_hi - 1;
    double t = (theta_cut_deg - theta_s[idx_lo]) /
               (theta_s[idx_hi] - theta_s[idx_lo]);
    double P11_cap = (1.0 - t) * P11_s[idx_lo] + t * P11_s[idx_hi];
    double P12_cap = (1.0 - t) * P12_s[idx_lo] + t * P12_s[idx_hi];
    double P33_cap = (1.0 - t) * P33_s[idx_lo] + t * P33_s[idx_hi];

    /* Capped arrays (roof variant): flat cap from 0 to idx_cut */
    double* P11_capped = (double*)malloc(N * sizeof(double));
    double* P12_capped = (double*)malloc(N * sizeof(double));
    double* P33_capped = (double*)malloc(N * sizeof(double));
    if (!P11_capped || !P12_capped || !P33_capped) {
        free(P11_capped); free(P12_capped); free(P33_capped);
        free(theta_s); free(P11_s); free(P12_s); free(P33_s); free(idx);
        return -1;
    }
    for (int i = 0; i < N; ++i) {
        if (i < idx_cut) {
            P11_capped[i] = P11_cap;
            P12_capped[i] = P12_cap;
            P33_capped[i] = P33_cap;
        } else {
            P11_capped[i] = P11_s[i];
            P12_capped[i] = P12_s[i];
            P33_capped[i] = P33_s[i];
        }
    }

    /* Forward peak mass: f = (1/2) int_0^theta_cut (P11_orig - P11_cap) sin(theta) dtheta
     * v1.02 fix: include theta_cut as explicit integration endpoint via interpolated
     * P11 at theta_cut (which equals P11_cap by construction → boundary integrand = 0).
     * This makes the integral independent of where theta_cut falls relative to the
     * tabulated grid. */
    double f = 0.0;
    for (int i = 1; i < idx_cut; ++i) {
        double th_a = theta_s[i - 1] * M_PI / 180.0;
        double th_b = theta_s[i]     * M_PI / 180.0;
        double g_a = (P11_s[i - 1] - P11_cap) * sin(th_a);
        double g_b = (P11_s[i]     - P11_cap) * sin(th_b);
        f += 0.5 * (th_b - th_a) * (g_a + g_b);
    }
    /* Last interval: theta_s[idx_lo] -> theta_cut (interpolated endpoint).
     * At theta_cut, P11_orig(theta_cut) = P11_cap by interpolation -> integrand = 0. */
    {
        double th_a = theta_s[idx_cut - 1] * M_PI / 180.0;
        double th_b = theta_cut_deg * M_PI / 180.0;
        double g_a = (P11_s[idx_cut - 1] - P11_cap) * sin(th_a);
        double g_b = 0.0;  /* P11_orig(theta_cut) - P11_cap = 0 */
        f += 0.5 * (th_b - th_a) * (g_a + g_b);
    }
    f *= 0.5;

    if (f <= 0.0 || f >= 0.999) {
        for (int i = 0; i < N; ++i) {
            P11_out[i] = P11_in[i]; P12_out[i] = P12_in[i]; P33_out[i] = P33_in[i];
        }
        *f_out = 0.0;
        free(P11_capped); free(P12_capped); free(P33_capped);
        free(theta_s); free(P11_s); free(P12_s); free(P33_s); free(idx);
        return 0;
    }

    /* Normalize by (1-f) and restore original ordering */
    double inv_denom = 1.0 / (1.0 - f);
    for (int i = 0; i < N; ++i) {
        int orig_i = idx[i];
        P11_out[orig_i] = P11_capped[i] * inv_denom;
        P12_out[orig_i] = P12_capped[i] * inv_denom;
        P33_out[orig_i] = P33_capped[i] * inv_denom;
    }
    *f_out = f;

    free(P11_capped); free(P12_capped); free(P33_capped);
    free(theta_s); free(P11_s); free(P12_s); free(P33_s); free(idx);
    return 0;
}

/* ======================================================================== */
/* Real-angle log10-linear forward-peak truncation.
 *
 * Algorithm:
 *   1. Pick two anchor angles θ_1, θ_2 at cos θ = mu1, mu2 (defaults
 *      0.8 → 36.87°, 0.94 → 19.95°).
 *   2. Replace P11 in the forward region θ < θ_2 with a linear segment
 *      in log10(P11) vs θ, extrapolated through the two anchor values.
 *   3. Compute the truncation coefficient A = 2·(1 − β₀_truncated),
 *      where β₀ is the zeroth Legendre moment (normalization) of the
 *      truncated P11.  A is the fraction of forward-peak mass removed.
 *   4. If A is below a small threshold (default 0.1) the truncation is
 *      cancelled — appropriate for phase functions without a strong
 *      forward peak (e.g., highly absorbing aerosols).
 *   5. Scale P12, P22, P33 by the pointwise ratio P11_truncated / P11_orig.
 *
 * This forward-peak handling differs from the two truncation schemes
 * already implemented above:
 *
 *   - delta-fit (Hu & Stamnes 2000): replaces the forward peak with a
 *     constant roof cap and renormalizes the entire phase matrix by
 *     1/(1−f).
 *   - delta-M (Wiscombe 1977): operates in Legendre space, separates
 *     the forward peak into a δ-function plus a smooth residual, and
 *     truncates Legendre coefficients above order N.
 *
 * The present operator works in real (θ) space, preserves the smooth
 * shape of P11 outside the forward cone, and produces a truncation
 * coefficient A directly tied to the fraction of phase-function mass
 * removed.  The caller renormalizes layer optical thickness and single-
 * scattering albedo accordingly (see rt_aerosol_runtime.c).
 *
 * References:
 *   Potter, J. F. (1970). The delta function approximation in radiative
 *     transfer theory. J. Atmos. Sci., 27, 943–949.
 *   Wiscombe, W. J. (1977). The delta-M method: rapid yet accurate
 *     radiative flux calculations for strongly asymmetric phase
 *     functions. J. Atmos. Sci., 34, 1408–1422.
 *   Hu, Y.-X., and Stamnes, K. (2000). An accurate parameterization of
 *     the radiative properties of water clouds. J. Atmos. Sci., 57,
 *     1568–1591.
 *
 * Numerically cross-checked against an independent vector RT
 * implementation (Chami et al. 2015, Opt. Express 23, 27829) on the
 * FresnelOcean aerosol benchmark used by this project (case 721:
 * reference I = 0.11699, this implementation = 0.11782, Δ = +0.71%).
 *
 * f_out returns A (truncation coefficient, typically 0.5–1.5 for
 * forward-peaked aerosols); semantically distinct from the f_out of
 * delta-fit and delta-M (Wiscombe forward fraction, range 0–1).
 */
int rt_aerosol_loglin_truncate(const double* P11_in, const double* P12_in, const double* P33_in,
                               const double* theta_deg, int N,
                               double mu1, double mu2,
                               double trunc_threshold,
                               double* P11_out, double* P12_out, double* P33_out,
                               double* f_out) {
    /* Sort ascending in theta if needed */
    int ascending = 1;
    for (int i = 1; i < N; ++i) if (theta_deg[i] < theta_deg[i - 1]) { ascending = 0; break; }

    double* theta_s = (double*)malloc(N * sizeof(double));
    double* P11_s   = (double*)malloc(N * sizeof(double));
    int* idx        = (int*)malloc(N * sizeof(int));
    if (!theta_s || !P11_s || !idx) {
        free(theta_s); free(P11_s); free(idx);
        return -1;
    }
    if (ascending) {
        for (int i = 0; i < N; ++i) { theta_s[i] = theta_deg[i]; P11_s[i] = P11_in[i]; idx[i] = i; }
    } else {
        for (int i = 0; i < N; ++i) {
            int j = N - 1 - i;
            theta_s[i] = theta_deg[j]; P11_s[i] = P11_in[j]; idx[i] = j;
        }
    }

    /* Thresholds: theta1 = acos(MU1), theta2 = acos(MU2). MU1=0.8 → 36.87°, MU2=0.94 → 19.95°.
     * Forward region: theta < theta2 (smaller theta = more forward). */
    double theta1_deg = acos(mu1) * 180.0 / M_PI;   /* ~36.87° */
    double theta2_deg = acos(mu2) * 180.0 / M_PI;   /* ~19.95° */

    /* v1.02 fix: linear interpolation of P11 AT theta1_deg, theta2_deg (not at
     * nearest tabulated row). Previous version used P11_s[K], P11_s[KK] (closest
     * <= row) which introduced grid dependence. */

    /* Find bracket indices for theta1_deg, theta2_deg */
    int K_hi = 0;
    while (K_hi < N && theta_s[K_hi] < theta1_deg) K_hi++;
    int KK_hi = 0;
    while (KK_hi < N && theta_s[KK_hi] < theta2_deg) KK_hi++;

    /* Validity: theta1 > theta2, both within tabulated range, well-separated */
    if (K_hi >= N || K_hi == 0 || KK_hi >= N || KK_hi == 0 || K_hi <= KK_hi + 1) {
        /* No valid anchors → no truncation */
        for (int i = 0; i < N; ++i) { P11_out[i] = P11_in[i]; P12_out[i] = P12_in[i]; P33_out[i] = P33_in[i]; }
        *f_out = 0.0;
        free(theta_s); free(P11_s); free(idx);
        return 0;
    }

    /* Linear-interpolated P11 values AT theta1 and theta2 (exact anchor angles) */
    int K_lo  = K_hi  - 1;
    int KK_lo = KK_hi - 1;
    double t_K  = (theta1_deg - theta_s[K_lo])  / (theta_s[K_hi]  - theta_s[K_lo]);
    double t_KK = (theta2_deg - theta_s[KK_lo]) / (theta_s[KK_hi] - theta_s[KK_lo]);
    double P11_at_theta1 = (1.0 - t_K)  * P11_s[K_lo]  + t_K  * P11_s[K_hi];
    double P11_at_theta2 = (1.0 - t_KK) * P11_s[KK_lo] + t_KK * P11_s[KK_hi];

    /* Transition index KK = first tabulated row with theta_s[KK] >= theta2 */
    int KK = KK_hi;

    if (P11_at_theta1 <= 0.0 || P11_at_theta2 <= 0.0) {
        for (int i = 0; i < N; ++i) { P11_out[i] = P11_in[i]; P12_out[i] = P12_in[i]; P33_out[i] = P33_in[i]; }
        *f_out = 0.0;
        free(theta_s); free(P11_s); free(idx);
        return 0;
    }

    /* Slope AA in log10(P11) vs theta (radians) based on INTERPOLATED anchor values */
    double theta1_rad = theta1_deg * M_PI / 180.0;
    double theta2_rad = theta2_deg * M_PI / 180.0;
    double AA = (log10(P11_at_theta2) - log10(P11_at_theta1)) / (theta2_rad - theta1_rad);

    /* Build truncated P11 on the sorted grid:
     *   for i < KK (theta_s[i] < theta2): P11_new = 10^(log10(P11_at_theta2) + AA*(theta_i - theta2_rad))
     *   for i >= KK: unchanged.
     * v1.02 fix: anchor X1 = log10(P11_at_theta2) is now INTERPOLATED at exact
     * theta2, not at tabulated row P11_s[KK]. Grid-independent. */
    double* P11_trunc_s = (double*)malloc(N * sizeof(double));
    if (!P11_trunc_s) {
        free(theta_s); free(P11_s); free(idx);
        return -1;
    }
    double X1 = log10(P11_at_theta2);
    for (int i = 0; i < N; ++i) {
        if (i < KK) {
            double theta_i_rad = theta_s[i] * M_PI / 180.0;
            double logp = X1 + AA * (theta_i_rad - theta2_rad);
            P11_trunc_s[i] = pow(10.0, logp);
        } else {
            P11_trunc_s[i] = P11_s[i];
        }
    }

    /* Compute BETA11(0) of truncated P11.
     * BETA11(0) = (1/2) ∫_{-1}^{+1} P11(μ) dμ = (1/2) ∫_0^π P11(θ) sin(θ) dθ
     * → use chi_0 from rt_aerosol_compute_legendre_moments with n_max=1. */
    double chi0_trunc = 0.0;
    int rc = rt_aerosol_compute_legendre_moments(P11_trunc_s, theta_s, N,
                                                 1, 3601, &chi0_trunc);
    if (rc != 0) {
        free(P11_trunc_s); free(theta_s); free(P11_s); free(idx);
        return -1;
    }

    double A_trunc = 2.0 * (1.0 - chi0_trunc);

    /* Cancel if too weak */
    if (A_trunc < trunc_threshold) {
        for (int i = 0; i < N; ++i) { P11_out[i] = P11_in[i]; P12_out[i] = P12_in[i]; P33_out[i] = P33_in[i]; }
        *f_out = 0.0;
        free(P11_trunc_s); free(theta_s); free(P11_s); free(idx);
        return 0;
    }

    /* Apply truncation: restore original ordering, scale P12/P33 by ratio. */
    for (int i = 0; i < N; ++i) {
        int orig_i = idx[i];
        double ratio = (P11_s[i] > 0.0) ? (P11_trunc_s[i] / P11_s[i]) : 1.0;
        P11_out[orig_i] = P11_trunc_s[i];
        P12_out[orig_i] = P12_in[orig_i] * ratio;
        P33_out[orig_i] = P33_in[orig_i] * ratio;
    }
    *f_out = A_trunc;

    free(P11_trunc_s); free(theta_s); free(P11_s); free(idx);
    return 0;
}

/* ======================================================================== */
/* delta-M truncation (Wiscombe 1977) in Legendre space, order N_trunc.
 *
 *   g_l  = chi_l / (2l+1)                (normalized moments, g_0 = 1)
 *   f    = g_{N_trunc}                   (forward delta-function fraction)
 *   g'_l = (g_l - f) / (1 - f),  l < N_trunc     [Wiscombe Eq. (5)]
 *   i.e. chi'_l = (chi_l - f (2l+1)) / (1 - f)
 * and the truncated P11 is re-synthesized as sum_l chi'_l P_l(mu)
 * (clipped at 0 against Gibbs undershoot).
 *
 * P12/P33 handling — DOCUMENTED HEURISTIC, not standard theory:
 * the forward delta mass assigned to P12/P33 is scaled by the LAST-MOMENT
 * ratio r = chi_X[N_trunc] / chi_11[N_trunc], i.e. each element is assumed
 * to carry the same forward-peak shape as P11 weighted by its high-order
 * moment content.  For P12 (which vanishes at exact forward scattering for
 * spheres) this over-assigns delta mass in principle; kept because it
 * matches the validated v1/Python reference implementation bit-for-bit.
 * NOTE: a "+48% high-omega 412 nm" defect once attributed to this
 * scheme's water-RT application was RETRACTED (no measurement record;
 * S-004, 2026-07-05).  The P12/P33 heuristic above remains a documented
 * approximation on its own merits; do not modify before validation
 * matrix item #5 completes. */
int rt_aerosol_delta_m_truncate(const double* P11_in, const double* P12_in, const double* P33_in,
                              const double* theta_deg, int N, int N_trunc,
                              double* P11_out, double* P12_out, double* P33_out,
                              double* f_out) {
    double* chi = (double*)malloc((N_trunc + 1) * sizeof(double));
    if (!chi) return -1;
    /* Use fine_grid_points=3601 for integration accuracy, matches Python */
    if (rt_aerosol_compute_legendre_moments(P11_in, theta_deg, N, N_trunc + 1, 3601, chi) != 0) {
        free(chi); return -1;
    }
    double f = chi[N_trunc] / (2.0 * N_trunc + 1.0);
    if (f <= 0.0 || f >= 0.999) {
        for (int i = 0; i < N; ++i) {
            P11_out[i] = P11_in[i]; P12_out[i] = P12_in[i]; P33_out[i] = P33_in[i];
        }
        *f_out = 0.0; free(chi); return 0;
    }

    /* Truncated moments for P11 */
    double* chi_t = (double*)malloc(N_trunc * sizeof(double));
    if (!chi_t) { free(chi); return -1; }
    for (int l = 0; l < N_trunc; ++l) {
        chi_t[l] = (chi[l] - f * (2.0 * l + 1.0)) / (1.0 - f);
    }

    /* FIX: Also truncate P12 and P33 in Legendre space using same f.
     * Standard vector delta-M (Hu-Stamnes 2000): subtract delta-mass forward peak
     * from all Mueller elements consistently. */
    double* chi12 = (double*)malloc((N_trunc + 1) * sizeof(double));
    double* chi33 = (double*)malloc((N_trunc + 1) * sizeof(double));
    double* chi12_t = (double*)malloc(N_trunc * sizeof(double));
    double* chi33_t = (double*)malloc(N_trunc * sizeof(double));
    if (!chi12 || !chi33 || !chi12_t || !chi33_t) {
        free(chi); free(chi_t); free(chi12); free(chi33); free(chi12_t); free(chi33_t);
        return -1;
    }
    if (rt_aerosol_compute_legendre_moments(P12_in, theta_deg, N, N_trunc + 1, 3601, chi12) != 0 ||
        rt_aerosol_compute_legendre_moments(P33_in, theta_deg, N, N_trunc + 1, 3601, chi33) != 0) {
        free(chi); free(chi_t); free(chi12); free(chi33); free(chi12_t); free(chi33_t);
        return -1;
    }
    /* Forward-peak ratio of P12, P33 inferred from last-moment ratio relative
     * to P11. Guard against chi[N_trunc] ≈ 0 (matches Python vrt_solver.py
     * delta_m_truncate behavior). In practice chi[N_trunc] > 0 is guaranteed
     * by the (f > 0) check above, but the explicit guard makes intent clear. */
    double r12 = 0.0, r33 = 0.0;
    if (fabs(chi[N_trunc]) > 1e-300) {
        r12 = chi12[N_trunc] / chi[N_trunc];
        r33 = chi33[N_trunc] / chi[N_trunc];
    }
    for (int l = 0; l < N_trunc; ++l) {
        chi12_t[l] = (chi12[l] - f * (2.0 * l + 1.0) * r12) / (1.0 - f);
        chi33_t[l] = (chi33[l] - f * (2.0 * l + 1.0) * r33) / (1.0 - f);
    }

    /* Reconstruct P11', P12', P33' on input theta grid via Sigma chi_l Pl(mu) */
    double* Plbuf = (double*)malloc(N_trunc * sizeof(double));
    if (!Plbuf) {
        free(chi); free(chi_t); free(chi12); free(chi33); free(chi12_t); free(chi33_t);
        return -1;
    }
    for (int i = 0; i < N; ++i) {
        double mu = cos(theta_deg[i] * M_PI / 180.0);
        legendre_P_all(N_trunc, mu, Plbuf);
        double sum_p11 = 0.0, sum_p12 = 0.0, sum_p33 = 0.0;
        for (int l = 0; l < N_trunc; ++l) {
            sum_p11 += chi_t[l]   * Plbuf[l];
            sum_p12 += chi12_t[l] * Plbuf[l];
            sum_p33 += chi33_t[l] * Plbuf[l];
        }
        if (sum_p11 < 0.0) sum_p11 = 0.0;
        P11_out[i] = sum_p11;
        P12_out[i] = sum_p12;
        P33_out[i] = sum_p33;
    }
    *f_out = f;

    free(Plbuf); free(chi_t); free(chi); free(chi12); free(chi33);
    free(chi12_t); free(chi33_t);
    return 0;
}

/* ========================================================================
 * Vector phase-matrix moment expansion
 * ======================================================================== */


/* gamma_l = (2l+1)/2 * Integral P12(mu) P2_l(mu) dmu on a PCHIP-densified
 * grid — same numerical pattern as compute_legendre_moments but with the
 * generalized basis P2_l (== rt_kernel rsl at m=0; see file header for the
 * verified identity).  NOTE the dense grid here spans only the TABULATED
 * theta range (vs fixed 0..180 in compute_legendre_moments) — consistent
 * only for full-range tables (file-header constraint).  The pointwise
 * weights (half-spans) are the trapezoid rule written per node, handling
 * the non-uniform mu spacing of a uniform-theta grid. */
static int compute_p2_moments(const double* P12_values, const double* theta_deg, int N,
                                int L_max, int fine_grid_points, double* gammal_out) {
    /* Sort ascending if needed */
    int ascending = 1;
    for (int i = 1; i < N; ++i) if (theta_deg[i] < theta_deg[i - 1]) { ascending = 0; break; }
    double* theta_sorted = (double*)malloc(N * sizeof(double));
    double* P_sorted     = (double*)malloc(N * sizeof(double));
    if (!theta_sorted || !P_sorted) { free(theta_sorted); free(P_sorted); return -1; }
    if (ascending) {
        memcpy(theta_sorted, theta_deg, N * sizeof(double));
        memcpy(P_sorted,     P12_values, N * sizeof(double));
    } else {
        for (int i = 0; i < N; ++i) {
            theta_sorted[i] = theta_deg[N - 1 - i];
            P_sorted[i]     = P12_values[N - 1 - i];
        }
    }

    /* PCHIP fit in θ */
    pchip_t spline;
    if (pchip_build(&spline, theta_sorted, P_sorted, N) != 0) {
        free(theta_sorted); free(P_sorted); return -1;
    }

    /* Dense uniform θ grid + trapezoid in μ */
    double* pol = (double*)malloc((L_max + 1) * sizeof(double));
    double* integrand = (double*)malloc(fine_grid_points * sizeof(double));
    double* mu_grid   = (double*)malloc(fine_grid_points * sizeof(double));
    if (!pol || !integrand || !mu_grid) {
        free(pol); free(integrand); free(mu_grid);
        pchip_free(&spline);
        free(theta_sorted); free(P_sorted);
        return -1;
    }

    /* Initialize γ_l = 0 */
    for (int l = 0; l <= L_max; ++l) gammal_out[l] = 0.0;

    /* For each Legendre order l, integrate ∫_{-1}^{+1} P12(μ) P²_l(μ) dμ
     * Trapezoid in μ (descending as θ ascends), so we go from θ=π to θ=0
     * in PCHIP eval but reverse-fill so μ ascends. */
    double dtheta = (theta_sorted[N - 1] - theta_sorted[0]) / (fine_grid_points - 1);

    /* Precompute μ_grid and P12 values (one PCHIP eval per fine-grid point) */
    for (int i = 0; i < fine_grid_points; ++i) {
        double th = theta_sorted[0] + i * dtheta;
        double P_th = pchip_eval(&spline, th);
        double mu = cos(th * M_PI / 180.0);
        /* Reverse-fill: smallest μ corresponds to largest θ */
        int j = fine_grid_points - 1 - i;
        mu_grid[j]   = mu;
        integrand[j] = P_th;     /* P12 at this μ — multiplied by P²_l later */
    }

    /* For each l, compute trapezoid of P12 × P²_l in μ */
    for (int i = 0; i < fine_grid_points; ++i) {
        rt_phase_spin2_m0(L_max, mu_grid[i], pol);
        double dmu = (i == 0)
                     ? (mu_grid[1] - mu_grid[0]) * 0.5
                     : (i == fine_grid_points - 1)
                       ? (mu_grid[fine_grid_points - 1] - mu_grid[fine_grid_points - 2]) * 0.5
                       : (mu_grid[i + 1] - mu_grid[i - 1]) * 0.5;
        for (int l = 0; l <= L_max; ++l) {
            gammal_out[l] += integrand[i] * pol[l] * dmu;
        }
    }

    /* Normalization: γ_l = (2l+1)/2 × ∫ */
    for (int l = 0; l <= L_max; ++l) {
        gammal_out[l] *= (2.0 * l + 1.0) * 0.5;
    }

    free(pol); free(integrand); free(mu_grid);
    pchip_free(&spline);
    free(theta_sorted); free(P_sorted);
    return 0;
}

/* ======================================================================== */
int rt_aerosol_compute_vector_legendre(
    const double* P11, const double* P12, const double* P33,
    const double* theta_deg, int N,
    int L_max, int fine_grid_points,
    double* betal_out, double* gammal_out,
    double* alphal_out, double* zetal_out)
{
    if (!aerosol_table_full_range_ok(theta_deg, N, "compute_vector_legendre")) return -1;
    if (!P11 || !P12 || !P33 || !theta_deg || N < 3 || L_max < 2) return -1;
    if (!betal_out || !gammal_out || !alphal_out || !zetal_out)   return -1;

    /* β_l of P11 */
    if (rt_aerosol_compute_legendre_moments(P11, theta_deg, N, L_max + 1,
                                            fine_grid_points, betal_out) != 0) return -1;

    /* γ_l of P12 (Mishchenko P²_l) */
    if (compute_p2_moments(P12, theta_deg, N, L_max, fine_grid_points, gammal_out) != 0) return -1;

    /* Ordinary Legendre moments of P33.  Together with P11=P22 these
     * determine the two diagonal spin-2 coefficient families. */
    double* p33_moments = (double*)calloc(L_max + 1, sizeof(double));
    if (!p33_moments) return -1;
    if (rt_aerosol_compute_legendre_moments(P33, theta_deg, N, L_max + 1,
                                            fine_grid_points, p33_moments) != 0) {
        free(p33_moments); return -1;
    }

    if (rt_phase_complete_spherical(L_max, betal_out, betal_out, p33_moments,
                                    alphal_out, zetal_out) != 0) {
        free(p33_moments);
        return -1;
    }

    if (rt_phase_normalize_vector(L_max, betal_out, gammal_out,
                                  alphal_out, zetal_out) != 0) {
        free(p33_moments);
        return -1;
    }

    free(p33_moments);
    return 0;
}


/* ========================================================================
 * Gauss-quadrature vector moment expansion for external phase tables
 * ======================================================================== */
typedef struct { double x; double y; } rt_xy_pair_t;
static int rt_xy_pair_cmp(const void *a, const void *b) {
    const rt_xy_pair_t *pa = (const rt_xy_pair_t*)a;
    const rt_xy_pair_t *pb = (const rt_xy_pair_t*)b;
    return (pa->x > pb->x) - (pa->x < pb->x);
}

static int build_mu_table_unique(const double *P, const double *theta_deg, int N,
                                 double **x_out, double **y_out, int *M_out) {
    if (!P || !theta_deg || N < 3 || !x_out || !y_out || !M_out) return -1;
    rt_xy_pair_t *tmp = (rt_xy_pair_t*)calloc((size_t)N, sizeof(rt_xy_pair_t));
    if (!tmp) return -1;
    for (int i = 0; i < N; ++i) {
        double mu = cos(theta_deg[i] * M_PI / 180.0);
        if (mu < -1.0) mu = -1.0;
        if (mu >  1.0) mu =  1.0;
        tmp[i].x = mu;
        tmp[i].y = P[i];
    }
    qsort(tmp, (size_t)N, sizeof(rt_xy_pair_t), rt_xy_pair_cmp);
    double *x = (double*)malloc((size_t)N * sizeof(double));
    double *y = (double*)malloc((size_t)N * sizeof(double));
    if (!x || !y) { free(tmp); free(x); free(y); return -1; }
    int M = 0;
    for (int i = 0; i < N; ++i) {
        if (M > 0 && fabs(tmp[i].x - x[M-1]) < 1e-14) {
            /* If duplicate μ appears, keep the latest value.  Valid phase
             * tables normally have unique angles; this branch is defensive. */
            y[M-1] = tmp[i].y;
        } else {
            x[M] = tmp[i].x;
            y[M] = tmp[i].y;
            ++M;
        }
    }
    free(tmp);
    if (M < 3) { free(x); free(y); return -1; }
    *x_out = x; *y_out = y; *M_out = M;
    return 0;
}

static int clamped_cubic_spline_build(const double *x, const double *y, int n, double *y2) {
    if (!x || !y || !y2 || n < 3) return -1;
    double *u = (double*)calloc((size_t)(n - 1), sizeof(double));
    if (!u) return -1;
    double dx0 = x[1] - x[0];
    double dxn = x[n-1] - x[n-2];
    if (!(dx0 > 0.0) || !(dxn > 0.0)) { free(u); return -1; }
    double yp1 = (y[1] - y[0]) / dx0;
    double ypn = (y[n-1] - y[n-2]) / dxn;
    y2[0] = -0.5;
    u[0] = (3.0 / dx0) * ((y[1] - y[0]) / dx0 - yp1);
    for (int i = 1; i < n - 1; ++i) {
        double dxm = x[i] - x[i-1];
        double dxp = x[i+1] - x[i];
        if (!(dxm > 0.0) || !(dxp > 0.0)) { free(u); return -1; }
        double sig = dxm / (x[i+1] - x[i-1]);
        double p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        double dd = ((y[i+1] - y[i]) / dxp - (y[i] - y[i-1]) / dxm);
        u[i] = (6.0 * dd / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }
    double qn = 0.5;
    double un = (3.0 / dxn) * (ypn - (y[n-1] - y[n-2]) / dxn);
    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
    for (int k = n - 2; k >= 0; --k) {
        y2[k] = y2[k] * y2[k+1] + u[k];
    }
    free(u);
    return 0;
}

static double clamped_cubic_spline_eval(const double *x, const double *y, const double *y2, int n, double xq) {
    if (xq <= x[0]) return y[0];
    if (xq >= x[n-1]) return y[n-1];
    int klo = 0, khi = n - 1;
    while (khi - klo > 1) {
        int k = (khi + klo) >> 1;
        if (x[k] > xq) khi = k; else klo = k;
    }
    double h = x[khi] - x[klo];
    if (!(h > 0.0)) return y[klo];
    double a = (x[khi] - xq) / h;
    double b = (xq - x[klo]) / h;
    return a * y[klo] + b * y[khi] +
           ((a*a*a - a) * y2[klo] + (b*b*b - b) * y2[khi]) * (h*h) / 6.0;
}

int rt_aerosol_compute_vector_legendre_gauss(
    const double* P11, const double* P12, const double* P33,
    const double* theta_deg, int N,
    int L_max, int mie_n_mu,
    double* betal_out, double* gammal_out,
    double* alphal_out, double* zetal_out)
{
    if (!aerosol_table_full_range_ok(theta_deg, N, "compute_vector_legendre_gauss")) return -1;
    if (!P11 || !P12 || !P33 || !theta_deg || N < 3 || L_max < 2 || mie_n_mu < 2) return -1;
    if (!betal_out || !gammal_out || !alphal_out || !zetal_out) return -1;

    double *x11=NULL, *y11=NULL, *x12=NULL, *y12=NULL, *x33=NULL, *y33=NULL;
    int n11=0, n12=0, n33=0;
    if (build_mu_table_unique(P11, theta_deg, N, &x11, &y11, &n11) != 0) return -1;
    if (build_mu_table_unique(P12, theta_deg, N, &x12, &y12, &n12) != 0) { free(x11); free(y11); return -1; }
    if (build_mu_table_unique(P33, theta_deg, N, &x33, &y33, &n33) != 0) { free(x11); free(y11); free(x12); free(y12); return -1; }

    double *y2_11=(double*)calloc((size_t)n11,sizeof(double));
    double *y2_12=(double*)calloc((size_t)n12,sizeof(double));
    double *y2_33=(double*)calloc((size_t)n33,sizeof(double));
    double *mu=(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *w =(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *PL=(double*)malloc((size_t)(L_max+1)*sizeof(double));
    double *P2=(double*)malloc((size_t)(L_max+1)*sizeof(double));
    double *beta22=(double*)calloc((size_t)(L_max+1),sizeof(double));
    double *deltal=(double*)calloc((size_t)(L_max+1),sizeof(double));
    if (!y2_11 || !y2_12 || !y2_33 || !mu || !w || !PL || !P2 || !beta22 || !deltal) {
        free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
        free(y2_11); free(y2_12); free(y2_33); free(mu); free(w); free(PL); free(P2); free(beta22); free(deltal);
        return -1;
    }
    if (clamped_cubic_spline_build(x11,y11,n11,y2_11)!=0 ||
        clamped_cubic_spline_build(x12,y12,n12,y2_12)!=0 ||
        clamped_cubic_spline_build(x33,y33,n33,y2_33)!=0 ||
        rt_quadrature_gauss_legendre_pos(mie_n_mu, mu, w) != 0) {
        free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
        free(y2_11); free(y2_12); free(y2_33); free(mu); free(w); free(PL); free(P2); free(beta22); free(deltal);
        return -1;
    }

    for (int l=0; l<=L_max; ++l) {
        betal_out[l]=0.0; gammal_out[l]=0.0; alphal_out[l]=0.0; zetal_out[l]=0.0;
    }

    for (int ih=0; ih<mie_n_mu; ++ih) {
        for (int sgn=-1; sgn<=1; sgn+=2) {
            double xq = (sgn < 0) ? -mu[ih] : mu[ih];
            double wt = w[ih];
            double p11 = clamped_cubic_spline_eval(x11,y11,y2_11,n11,xq);
            double p12 = clamped_cubic_spline_eval(x12,y12,y2_12,n12,xq);
            double p33 = clamped_cubic_spline_eval(x33,y33,y2_33,n33,xq);
            /* v1.09 S-006 fix: legendre_P_all fills out[0..n_max-1] ("count"
             * convention; see the reconstruction caller at ~line 619 which
             * allocates exactly N_trunc).  This caller allocates L_max+1 and
             * consumes PL[L_max], so it must request L_max+1 terms.  The old
             * call left PL[L_max] UNWRITTEN (malloc garbage): first calls in a
             * fresh heap read 0.0 (beta_Lmax silently lost), later calls read
             * stale heap doubles (beta_Lmax = deterministic garbage) — the
             * case-to-case drift of S-006.  The spin-2 basis routine
             * already fills pol[L_max] (max-degree convention), so only the
             * beta moments were affected. */
            legendre_P_all(L_max + 1, xq, PL);
            rt_phase_spin2_m0(L_max, xq, P2);
            for (int l=0; l<=L_max; ++l) {
                double c = (2.0*l + 1.0) * 0.5 * wt;
                betal_out[l] += c * p11 * PL[l];
                beta22[l]    += c * p11 * PL[l];  /* spherical external phase: P22 = P11 */
                deltal[l]    += c * p33 * PL[l];
                gammal_out[l]+= c * p12 * P2[l];
            }
        }
    }

    if (rt_phase_complete_spherical(L_max, betal_out, beta22, deltal,
                                    alphal_out, zetal_out) != 0 ||
        rt_phase_normalize_vector(L_max, betal_out, gammal_out,
                                  alphal_out, zetal_out) != 0) {
        free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
        free(y2_11); free(y2_12); free(y2_33); free(mu); free(w);
        free(PL); free(P2); free(beta22); free(deltal);
        return -1;
    }

    free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
    free(y2_11); free(y2_12); free(y2_33); free(mu); free(w); free(PL); free(P2); free(beta22); free(deltal);
    return 0;
}


/* Gauss-grid forward-peak smoothing followed by vector moments.
 *
 * Differences from rt_aerosol_loglin_truncate(): this routine deliberately
 * operates directly on the selected positive Gauss nodes.  The two
 * truncation anchors are the effective Gauss nodes just below MU1_TRONCA and
 * MU2_TRONCA, not interpolated exact-angle anchors.  P12/P33 are scaled by
 * P11_truncated/P11_raw pointwise before moment decomposition. */
int rt_aerosol_compute_vector_legendre_gauss_truncated(
    const double* P11, const double* P12, const double* P33,
    const double* theta_deg, int N,
    int L_max, int mie_n_mu,
    double mu1_tronca, double mu2_tronca, double trunc_threshold,
    double* betal_out, double* gammal_out,
    double* alphal_out, double* zetal_out,
    double* A_out)
{
    if (!aerosol_table_full_range_ok(theta_deg, N, "compute_vector_legendre_gauss_truncated")) return -1;
    if (A_out) *A_out = 0.0;
    if (!P11 || !P12 || !P33 || !theta_deg || N < 3 || L_max < 2 || mie_n_mu < 2) return -1;
    if (!betal_out || !gammal_out || !alphal_out || !zetal_out) return -1;

    double *x11=NULL, *y11=NULL, *x12=NULL, *y12=NULL, *x33=NULL, *y33=NULL;
    int n11=0, n12=0, n33=0;
    if (build_mu_table_unique(P11, theta_deg, N, &x11, &y11, &n11) != 0) return -1;
    if (build_mu_table_unique(P12, theta_deg, N, &x12, &y12, &n12) != 0) { free(x11); free(y11); return -1; }
    if (build_mu_table_unique(P33, theta_deg, N, &x33, &y33, &n33) != 0) { free(x11); free(y11); free(x12); free(y12); return -1; }

    double *y2_11=(double*)calloc((size_t)n11,sizeof(double));
    double *y2_12=(double*)calloc((size_t)n12,sizeof(double));
    double *y2_33=(double*)calloc((size_t)n33,sizeof(double));
    double *mu=(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *w =(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *p11p=(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *p12p=(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *p33p=(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *p11t=(double*)malloc((size_t)mie_n_mu*sizeof(double));
    double *P11n=(double*)malloc((size_t)(2*mie_n_mu)*sizeof(double));
    double *P12n=(double*)malloc((size_t)(2*mie_n_mu)*sizeof(double));
    double *P33n=(double*)malloc((size_t)(2*mie_n_mu)*sizeof(double));
    double *thn =(double*)malloc((size_t)(2*mie_n_mu)*sizeof(double));
    if (!y2_11 || !y2_12 || !y2_33 || !mu || !w || !p11p || !p12p || !p33p || !p11t || !P11n || !P12n || !P33n || !thn) {
        free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
        free(y2_11); free(y2_12); free(y2_33); free(mu); free(w);
        free(p11p); free(p12p); free(p33p); free(p11t); free(P11n); free(P12n); free(P33n); free(thn);
        return -1;
    }
    if (clamped_cubic_spline_build(x11,y11,n11,y2_11)!=0 ||
        clamped_cubic_spline_build(x12,y12,n12,y2_12)!=0 ||
        clamped_cubic_spline_build(x33,y33,n33,y2_33)!=0 ||
        rt_quadrature_gauss_legendre_pos(mie_n_mu, mu, w) != 0) {
        free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
        free(y2_11); free(y2_12); free(y2_33); free(mu); free(w);
        free(p11p); free(p12p); free(p33p); free(p11t); free(P11n); free(P12n); free(P33n); free(thn);
        return -1;
    }

    for (int i=0; i<mie_n_mu; ++i) {
        p11p[i] = clamped_cubic_spline_eval(x11,y11,y2_11,n11,mu[i]);
        p12p[i] = clamped_cubic_spline_eval(x12,y12,y2_12,n12,mu[i]);
        p33p[i] = clamped_cubic_spline_eval(x33,y33,y2_33,n33,mu[i]);
        p11t[i] = p11p[i];
    }

    int K=-1, KK=-1;
    for (int i=0; i<mie_n_mu; ++i) { if (mu[i] > mu1_tronca) { K=i-1; break; } }
    for (int i=0; i<mie_n_mu; ++i) { if (mu[i] > mu2_tronca) { KK=i-1; break; } }
    int trunc_applied = 0;
    double A = 0.0;
    if (K >= 0 && KK >= 0 && KK < mie_n_mu-1 && K < KK && p11p[K] > 0.0 && p11p[KK] > 0.0) {
        const double thK  = acos(mu[K]);
        const double thKK = acos(mu[KK]);
        double AA = (log10(p11p[KK]) - log10(p11p[K])) / (thKK - thK);
        double X1 = log10(p11p[KK]);
        double X2 = thKK;
        for (int i=KK+1; i<mie_n_mu; ++i) {
            double logp = X1 + AA * (acos(mu[i]) - X2);
            p11t[i] = pow(10.0, logp);
        }
        double beta0 = 0.0;
        for (int i=0; i<mie_n_mu; ++i) {
            double ppos = p11t[i];
            double pneg = clamped_cubic_spline_eval(x11,y11,y2_11,n11,-mu[i]);
            beta0 += 0.5 * w[i] * (ppos + pneg);
        }
        A = 2.0 * (1.0 - beta0);
        if (A >= trunc_threshold) trunc_applied = 1;
        else A = 0.0;
    }

    int n2 = 0;
    /* Build full -1..+1 node table in ascending theta_deg (0 -> 180) so the
     * generic Gauss moment routine can be reused.  Positive mu first
     * means forward angles first. */
    for (int i=mie_n_mu-1; i>=0; --i) {
        double ratio = (trunc_applied && p11p[i] > 0.0) ? (p11t[i] / p11p[i]) : 1.0;
        thn[n2] = acos(mu[i]) * 180.0 / M_PI;
        P11n[n2] = trunc_applied ? p11t[i] : p11p[i];
        P12n[n2] = p12p[i] * ratio;
        P33n[n2] = p33p[i] * ratio;
        n2++;
    }
    for (int i=0; i<mie_n_mu; ++i) {
        double xq = -mu[i];
        double p11 = clamped_cubic_spline_eval(x11,y11,y2_11,n11,xq);
        double p12 = clamped_cubic_spline_eval(x12,y12,y2_12,n12,xq);
        double p33 = clamped_cubic_spline_eval(x33,y33,y2_33,n33,xq);
        thn[n2] = acos(xq) * 180.0 / M_PI;
        P11n[n2] = p11;
        P12n[n2] = p12;
        P33n[n2] = p33;
        n2++;
    }

    int rc = rt_aerosol_compute_vector_legendre_gauss(P11n, P12n, P33n, thn, n2, L_max, mie_n_mu,
                                                            betal_out, gammal_out, alphal_out, zetal_out);
    if (A_out) *A_out = A;

    free(x11); free(y11); free(x12); free(y12); free(x33); free(y33);
    free(y2_11); free(y2_12); free(y2_33); free(mu); free(w);
    free(p11p); free(p12p); free(p33p); free(p11t); free(P11n); free(P12n); free(P33n); free(thn);
    return rc;
}
