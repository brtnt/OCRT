/* rt_aerosol.h
 *
 * Aerosol-specific RT utilities for V3 Phase 3:
 *   - Legendre moment expansion of phase function (β_l of P11)
 *   - Forward-peak truncation: delta-fit (Hu-Stamnes 2000) and
 *     delta-M (Wiscombe 1977)
 *
 * The functions here operate on phase function tables P11/P12/P33 vs
 * scattering angle θ (typically loaded from .mie via shared/mie_io.h).
 * They produce inputs that are consumed by:
 *   - rt_kernel_phase_fourier(ws, m, betal)       — scalar
 *   - rt_kernel_phase_fourier_pol(ws, m, gammal)  — vector (Phase 3+)
 *
 * The scalar and vector moment projections are expressed directly in the
 * Legendre and real spin-2 bases used by OCRT.  Completion of the diagonal
 * vector coefficient families is delegated to rt_phase_moments.c, whose
 * implementation is derived from the generalized spherical-function
 * identities rather than from a reference-program source layout.
 *
 * Implementation in rt_aerosol.c.
 */
#ifndef OCRT_V2_RT_AEROSOL_H
#define OCRT_V2_RT_AEROSOL_H

/* ─── Scalar Legendre moments ─────────────────────────────────────────
 *
 * χ_l = (2l+1) × (1/2) × ∫_{-1}^{+1} P(μ) P_l(μ) dμ
 *
 * Normalization: if (1/2) ∫ P sin(θ) dθ = 1, then χ_0 = 1.
 * Implementation: PCHIP in θ, evaluate on dense uniform grid, trapezoid in μ.
 * Matches Python compute_legendre_moments() to ~1e-10.
 *
 * Inputs:
 *   P_values[N]        phase function values (typically P11) at θ_deg[i]
 *   theta_deg[N]       scattering angles in degrees, monotonic
 *                      (ascending or descending — handled internally)
 *   N                  number of angle samples
 *   n_max              compute χ_l for l = 0..n_max-1 (n_max moments)
 *   fine_grid_points   dense sampling for trapezoid (3601 typical)
 * Output:
 *   chi_out[n_max]     Legendre moments
 *
 * Returns 0 on success, -1 on alloc failure.
 */
int rt_aerosol_compute_legendre_moments(
    const double* P_values, const double* theta_deg, int N,
    int n_max, int fine_grid_points,
    double* chi_out
);

/* ─── delta-fit truncation (Hu-Stamnes 2000, roof-cap variant) ────────
 *
 * Removes the forward peak of P11 by capping the value at θ_cut, then
 * scales all three Mueller elements by 1/(1-f) to restore normalization:
 *
 *   f = (1/2) ∫_0^θ_cut (P11_orig − P11_cap) sin θ dθ
 *
 * Inputs/outputs operate on ascending or descending θ grid.
 *
 * Inputs:
 *   P11_in/P12_in/P33_in[N]   original phase matrix elements
 *   theta_deg[N]              scattering angles in degrees
 *   N                         number of angles
 *   theta_cut_deg             forward-peak cap angle (e.g. 5°)
 * Outputs:
 *   P11_out/P12_out/P33_out[N]   truncated phase matrix
 *   *f_out                       forward-peak fraction (0 if no truncation)
 *
 * Returns 0 on success, -1 on alloc failure.
 */
int rt_aerosol_delta_fit_truncate(
    const double* P11_in, const double* P12_in, const double* P33_in,
    const double* theta_deg, int N,
    double theta_cut_deg,
    double* P11_out, double* P12_out, double* P33_out,
    double* f_out
);

/* ─── Real-angle log10-linear forward-peak truncation ─────────────────
 *
 * Implements forward-peak truncation by replacing the steep forward
 * region of P11 with a log10(P11)-θ linear segment through two anchor
 * angles, then renormalizing P12, P22, P33 by the pointwise ratio
 * P11_truncated / P11_original.  The truncation coefficient
 *
 *   A = 2 (1 − β₀_truncated)
 *
 * is the fraction of phase-function mass removed and is reported
 * through *f_out.  When A is below `trunc_threshold` (default 0.1)
 * the operator is a no-op — useful for highly absorbing aerosols
 * whose phase functions lack a strong forward peak.
 *
 * Defaults: mu1 = 0.8 (θ ≈ 36.87°), mu2 = 0.94 (θ ≈ 19.95°),
 *           trunc_threshold = 0.1.
 *
 * Distinct from delta-M (Legendre-space) and delta-fit (roof cap):
 * operates entirely in θ space and produces a different forward-peak
 * mass measure A (range typically 0.5–1.5 for forward-peaked
 * aerosols), not the Wiscombe forward fraction f (range 0–1).
 *
 * The caller must apply the matching optical-thickness renormalization:
 *   τ_sca' = τ_sca · (1 − A/2)
 *   τ_ext' = τ_ext · (1 − ω · A/2)
 *   ω'     = (1 − A/2)·ω / (1 − ω·A/2)
 * See `rt_aerosol_runtime.c` for the wired-up path.
 *
 * References:
 *   Potter, J. F. (1970). The delta function approximation in radiative
 *     transfer theory. J. Atmos. Sci., 27, 943–949.
 *   Wiscombe, W. J. (1977). The delta-M method. J. Atmos. Sci., 34,
 *     1408–1422.
 *   Hu, Y.-X., and Stamnes, K. (2000). An accurate parameterization
 *     of the radiative properties of water clouds. J. Atmos. Sci., 57,
 *     1568–1591.
 *
 * Verified to ~1% against an independent vector RT implementation
 * (Chami et al. 2015, Opt. Express 23, 27829) on the FresnelOcean
 * aerosol benchmark.
 *
 * Returns 0 on success, -1 on alloc failure.
 */
int rt_aerosol_loglin_truncate(
    const double* P11_in, const double* P12_in, const double* P33_in,
    const double* theta_deg, int N,
    double mu1, double mu2, double trunc_threshold,
    double* P11_out, double* P12_out, double* P33_out,
    double* f_out
);

/* ─── delta-M truncation (Wiscombe 1977, P11 only; scalar 1/(1-f) for P12/P33)
 *
 *   1. χ_l of P11, l = 0..N_trunc
 *   2. f = χ_{N_trunc} / (2N_trunc + 1)
 *   3. χ_l_trunc = (χ_l − f·(2l+1)) / (1−f) for l < N_trunc
 *   4. P11_trunc(θ) = Σ_{l<N_trunc} χ_l_trunc · P_l(cos θ)
 *   5. P12_trunc = P12 / (1−f),   P33_trunc = P33 / (1−f)
 *
 * Returns 0 on success, -1 on alloc failure. *f_out = 0 if no truncation
 * (f computed ≤ 0 or ≥ 0.999).
 */
int rt_aerosol_delta_m_truncate(
    const double* P11_in, const double* P12_in, const double* P33_in,
    const double* theta_deg, int N, int N_trunc,
    double* P11_out, double* P12_out, double* P33_out,
    double* f_out
);

/* ─── Vector generalized-moment expansion ─────────────────────────────
 *
 * Decomposes the aerosol 3-element phase matrix (P11, P12, P33) into
 * generalized-Legendre coefficients used by the vector RT kernel:
 *
 *   β_l   = (2l+1)/2 × ∫ P11(μ) P_l(μ) dμ              (standard P_l)
 *   γ_l   = (2l+1)/2 × ∫ P12(μ) P²_l(μ) dμ             (Mishchenko P²_l)
 *   δ_l   = (2l+1)/2 × ∫ P33(μ) P_l(μ) dμ              (ordinary basis)
 *   α_l   = i(i-1)/((i+1)(i+2)) × β_l
 *           − 4(2i+1)/(i(i-1)(i+1)(i+2)) × Σ           (recursive sum)
 *   ζ_l   = i(i-1)/((i+1)(i+2)) × δ_l                  (≡ 0 for spherical)
 *           − 4(2i+1)/(i(i-1)(i+1)(i+2)) × Σ
 *
 * Final normalization: α/β/γ/ζ /= β_0  (so β_0 ≡ 1 after).
 *
 * For reciprocal spherical-particle tables, P22=P11.  P33 remains an
 * independent tabulated element and contributes to both α_l and ζ_l.
 *
 * Inputs:
 *   P11, P12, P33[N]   phase matrix elements at θ_deg[i]
 *   theta_deg[N]       monotonic scattering angles (deg)
 *   N                  number of angle samples
 *   L_max              max Legendre order (compute l = 0..L_max)
 *   fine_grid_points   trapezoid grid for integration (3601 typical)
 * Outputs (each size L_max+1, caller-allocated):
 *   betal_out  γ_out  alphal_out  zetal_out
 *
 * Returns 0 on success, -1 on alloc/bad-args failure.
 */
int rt_aerosol_compute_vector_legendre(
    const double* P11, const double* P12, const double* P33,
    const double* theta_deg, int N,
    int L_max, int fine_grid_points,
    double* betal_out, double* gammal_out,
    double* alphal_out, double* zetal_out
);

/* ─── Gauss-quadrature vector Legendre expansion for external hydrosol phase ──
 *
 * For external spherical-particle phase tables:
 *   - interpolate P11/P12/P33 as a clamped cubic spline in μ=cosθ
 *   - integrate on the full ±Gauss grid defined by mie_n_mu positive nodes
 *   - build β11, γ12, α, ζ using the generalized spherical-function identities
 *
 * This is intended for water .mie → hydrosol moment parity diagnostics and
 * production vector-RT validation.  The generic dense-θ PCHIP path
 * above is retained for atmospheric OPAC compatibility.
 */

/* External hydrosol forward-peak truncation + vector Legendre expansion.
 * Applies a log10-linear forward-peak replacement on the selected
 * Gauss grid, then computes alpha/beta/gamma/zeta moments. A_out is twice
 * the removed forward fraction (A=2f). */
int rt_aerosol_compute_vector_legendre_gauss_truncated(
    const double* P11, const double* P12, const double* P33,
    const double* theta_deg, int N,
    int L_max, int mie_n_mu,
    double mu1_tronca, double mu2_tronca, double trunc_threshold,
    double* betal_out, double* gammal_out,
    double* alphal_out, double* zetal_out,
    double* A_out
);
int rt_aerosol_compute_vector_legendre_gauss(
    const double* P11, const double* P12, const double* P33,
    const double* theta_deg, int N,
    int L_max, int mie_n_mu,
    double* betal_out, double* gammal_out,
    double* alphal_out, double* zetal_out
);

#endif /* OCRT_V2_RT_AEROSOL_H */
