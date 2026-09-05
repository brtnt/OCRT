#ifndef OCRT_V2_RT_SOLVER_H
#define OCRT_V2_RT_SOLVER_H

#include "rt_types.h"
#include "rt_kernel.h"

/* Distinct failure code for an external Fourier bottom-source whose declared
 * shape does not match the atmospheric angular ring.  Public callers normally
 * see this only through an internal ocean-coupling pass, but keeping the code
 * explicit prevents a silent fallback to out-of-bounds mode slicing. */
enum { RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE = -7 };

/* Solve a single case end-to-end (Phase 1: scalar Rayleigh-only).
 *
 * Drives the Step 4-8 chain in one call:
 *   1. tau_R = rt_rayleigh_tau_model(lambda, opts->rayleigh_model)
 *   2. rt_atm_alloc + rt_atm_build_rayleigh
 *   3. rt_legendre_workspace_alloc once (l_max = max(2, m_max))
 *   4. for m = 0..opts->fourier_m_max:
 *        rt_legendre_compute -> rt_kernel_phase_fourier(NULL = pure Rayleigh)
 *        -> rt_solver_primary_source -> rt_solver_integrate(LINEAR)
 *        -> rt_solver_sos -> sample I_total^m at view direction (linear
 *           interpolation between two bracketing positive mu nodes)
 *   5. rt_solver_reconstruct_phi + rt_solver_reflectance_from_intensity
 *   6. Free.
 *
 * Inputs
 *   cs     Geometry + spectrum + surface kind (Phase 1 expects
 *          surface=BLACK and aerosol_on=0; rayleigh_on must be 1).
 *          sza_deg, vza_deg in [0, 90); raa_deg arbitrary; wavelength
 *          in nm, > 0.
 *   opts   Pass NULL to use rt_options_default(); otherwise caller
 *          provides.  Reads n_mu, n_layers, fourier_m_max,
 *          rayleigh_model, integration_method, sos.*.
 *
 * Output
 *   out    rho_I (Phase 1; rho_Q/rho_U set to 0); plus diagnostics
 *          (I_TOA, tau_R_total, sos_orders_per_m[m], sos_residual_per_m[m],
 *          n_orders_used = max over m, converged = AND over m).
 *
 * Returns 0 on success, -1 on bad arguments (NULL, sza/vza >= 90,
 * lambda <= 0, opts->fourier_m_max not in [0, 2] for Phase 1), -2 on
 * downstream non-finite (NaN propagation from primary / SOS).
 *
 * Thread-safety: pure function in inputs; allocates internally and
 * frees before return on every path.  Multiple threads OK with
 * separate (cs, opts, out) triples.
 */
/* rt_solve_case (스칼라 I-전용): 2026-07-14 삭제 (M1, Jae 지시). */

/* Vector RT extension of rt_solve_case (Step 6 — closes the I/Q/U
 * end-to-end pipeline introduced in Steps 1-5).
 *
 * Drives the same chain as rt_solve_case but with vector RT primitives:
 *   1. tau_R + atm setup (identical to scalar path).
 *   2. Per-m loop:
 *        rt_legendre_compute_pol  (scalar plm + vector rrl, rtl, atm->xrl/xtl)
 *        rt_kernel_phase_fourier_pol(NULL γ_l → Rayleigh-only Phase 2)
 *        rt_solver_primary_source     → primary I source
 *        rt_solver_primary_source_pol → primary Q, U source
 *        rt_solver_integrate × 3      → primary I, Q, U fields
 *        rt_solver_sos_pol            → total I, Q, U fields (per m)
 *        Sample (I^m, Q^m, U^m) at mu_view (linear interp or
 *        view_as_node, same option as scalar path).
 *   3. Fourier reconstruction:
 *        I_TOA = Σ_m  cos(m(Δφ+π)) × I^m   (rt_solver_reconstruct_phi)
 *        Q_TOA = Σ_m  cos(m(Δφ+π)) × Q^m   (rt_solver_reconstruct_phi)
 *        U_TOA = Σ_m  sin(m(Δφ+π)) × U^m   (rt_solver_reconstruct_phi_sin)
 *   4. rho_X = X_TOA / mu_solar for X ∈ {I, Q, U}.
 *
 * Numerical equivalence with scalar path:
 *   For Phase-2 Rayleigh-only at the same (cs, opts), rho_I from this
 *   function MUST match rho_I from rt_solve_case to within floating-
 *   point round-off.  This is the primary regression criterion (the
 *   I path uses the same primary source and per-m integration; only the
 *   SOS iteration is replaced by the vector-coupled rt_solver_sos_pol,
 *   which at m=0 has identically-zero U coupling and reduces to the
 *   scalar SOS for I plus an independent Q channel).
 *
 * Inputs / Output / Return codes / Thread-safety:
 *   Same contract as rt_solve_case.  Output additionally fills rho_Q,
 *   rho_U, Q_TOA, U_TOA; Phase-3 aerosol path NOT yet wired (NULL γ_l).
 */
int rt_solve_case_pol(const rt_case_t *cs, const rt_options_t *opts,
                      rt_result_t *out);

/* ─── Phase 3: Aerosol+Rayleigh single-case solve ─────────────────────
 *
 * Caller pre-computes the aerosol Legendre expansion coefficients (via
 * rt_aerosol_compute_vector_legendre) and passes them through this
 * struct.  The vector solver path is invoked, populating out->rho_I,
 * out->rho_Q, out->rho_U at the requested geometry.
 */
typedef struct {
    double tau_a;                      /* column AOD at λ */
    double ssa_a;                      /* aerosol single scattering albedo */
    int    L_max;                      /* max Legendre order (≥ 2) */
    const double *betal_aer;           /* size L_max+1 (β_0 = 1) */
    const double *gammal_aer;          /* size L_max+1 */
    const double *alphal_aer;          /* size L_max+1 */
    const double *zetal_aer;           /* size L_max+1 */
    /* Angle-space value-kernel path (Gibbs-free; mirrors the in-water
     * fixed_bulk_direct_phase_fourier).  When use_value_kernel != 0 the SOS
     * builds ws->phase_fourier_m by direct azimuth integration of these
     * (loglin-capped) phase VALUES instead of from the Legendre moments,
     * which the moment kernel cannot represent at the backscatter glory.
     * theta_phase/P11_phase/P12_phase/P33_phase have size n_ang_phase. */
    int    use_value_kernel;           /* 0 = moment kernel (legacy); 1 = value kernel */
    int    n_ang_phase;                /* number of phase angles (e.g. 361) */
    const double *theta_phase;         /* scattering angles (deg), ascending */
    const double *P11_phase;           /* capped P11 values */
    const double *P12_phase;           /* capped P12 values (vector) */
    const double *P33_phase;           /* capped P33 values (vector) */
} rt_aerosol_input_t;

/* Same logic as rt_solve_case_pol but with aerosol input.  When aer is
 * NULL or aer->tau_a == 0, behavior is bit-level identical to
 * rt_solve_case_pol (Rayleigh-only). */
int rt_solve_case_pol_aerosol(const rt_case_t *cs, const rt_options_t *opts,
                               const rt_aerosol_input_t *aer,
                               rt_result_t *out);

/* ─── Phase B (in-water RT) — RT_SURFACE_OCEAN solve ──────────────────
 *
 * Drives the in-water SOS solver (rt_water_rt_sos_pure) for ocean surface
 * cases. Stage 1 (B.4): in-water solver only, no atmospheric coupling
 * (assumes atmosphere-free; uses cs->F_sun as just-above-surface BOA solar
 * irradiance directly). Stage 2 (B.4): atmospheric SOS → T_aw → in-water
 * boundary full coupling.
 *
 * Caller pre-loads the water IOP LUTs (rt_water_iop_lut_t,
 * rt_water_iop_psi_T_lut_t — declared in rt_water_iop.h) and passes
 * pointers. Required LUTs:
 *   aw_lut    — pure water a_w, b_w spectral LUT (NASA z09)
 *   psi_T_lut — Röttgers 2014 T-correction ψ_T LUT (may be NULL → no T-corr)
 *
 * Populates ocean-specific fields in rt_result_t. Also populates legacy
 * rho_I, rho_Q, rho_U with above-water Stokes ρ = π·L/(F_sun·μ_sun_air)
 * for atmospheric-pipeline output compatibility.
 *
 * Returns 0 on success, negative on failure.
 */
#include "rt_water_iop.h"

int rt_solve_case_ocean(const rt_case_t *cs, const rt_options_t *opts,
                         const rt_water_iop_lut_t *aw_lut,
                         const rt_water_iop_psi_T_lut_t *psi_T_lut,
                         const rt_aerosol_input_t *aer,
                         rt_result_t *out);

/* Native coupled ocean-atmosphere angular LUT.  The atmospheric and in-water
 * SOS fields are shared by the complete VZA x RAA grid.  Water target
 * projection is performed once per VZA and all RAAs are Fourier reconstructed;
 * exact atmospheric samples are added only outside the Gauss interpolation
 * domain.  Results are row-major: results[iv*n_raa + ir]. */
typedef struct {
    int n_vza;
    const double *vza_deg;
    int n_raa;
    const double *raa_deg;
    rt_result_t *results;
    int coupled_calls;       /* diagnostic: expensive orchestration entries */
    int water_cold_solves;   /* diagnostic: expected 1 */
    int water_view_calls;    /* target postprocess calls: expected n_vza, not n_vza*n_raa */
    int near_nadir_exact_rows;
} rt_ocean_lut_grid_out_t;

/* Reset all worker-private atmospheric/all-view replay caches.
 * Call once at the beginning of every independent full-grid case so a worker
 * may process changing wavelengths/conditions without inheriting mutable
 * state from a preceding batch row. */
void rt_solver_all_view_cache_reset(void);

int rt_solve_case_ocean_lut(const rt_case_t *base_cs,
                             const rt_options_t *opts,
                             const rt_water_iop_lut_t *aw_lut,
                             const rt_water_iop_psi_T_lut_t *psi_T_lut,
                             const rt_aerosol_input_t *aer,
                             rt_ocean_lut_grid_out_t *grid);

/* ─── Phase B.4 Stage 2b — atmospheric BOA Stokes export ──────────────
 *
 * rt_atm_boa_export_alloc / _free are convenience helpers for the export
 * struct lifecycle. After allocation, the struct can be passed to
 * rt_solve_case_pol_for_ocean to capture the diffuse BOA downward I^m,
 * Q^m, U^m field. The captured data is consumed by the ocean coupling
 * (Snell refraction + T_aw Mueller → in-water boundary condition).
 */
int  rt_atm_boa_export_alloc(rt_atm_boa_export_t *be, int m_max, int n_mu);
void rt_atm_boa_export_free (rt_atm_boa_export_t *be);

int rt_solve_case_pol_for_ocean(const rt_case_t *cs, const rt_options_t *opts,
                                 const rt_aerosol_input_t *aer,
                                 rt_result_t *out,
                                 rt_atm_boa_export_t *boa_export);

/* v1.01 LUT path.  Same algorithm and inputs as rt_solve_case_pol_aerosol
 * (aer may be NULL for Rayleigh-only), but additionally reconstructs the
 * radiation field at every (vza, raa) point in lut_out's grid arrays.
 * SOS is executed ONCE; grid is filled from the same quadrature solution.
 * The single-point rt_result_t *out is also filled using cs->vza_deg /
 * cs->raa_deg as the representative direction.
 *
 * Caller must allocate lut_out->rho_I/Q/U as n_vza * n_raa doubles, and
 * populate vza_deg[0..n_vza-1] in [0, 90) and raa_deg[0..n_raa-1] in
 * [0, 360). */
int rt_solve_case_pol_lut(const rt_case_t *cs, const rt_options_t *opts,
                          const rt_aerosol_input_t *aer,
                          rt_result_t *out,
                          rt_lut_grid_out_t *lut_out);

/* Allocate / free a per-Fourier-order field. */
int  rt_field_alloc(rt_field_t *f, int n_levels, int n_mu);
void rt_field_free(rt_field_t *f);

/* ----------------------------------------------------------------------------
 * rt_solver_primary_source
 *
 *   Compute the m-th azimuthal Fourier component of the 1st-order
 *   (single-scattering) source J₁_m(τ_k, μ_j) on every (level k,
 *   direction j) pair.
 *
 *   Per-layer formula (matches OS.f L348-353, term-for-term):
 *
 *       source_first[k][j] = atm.ch[k] · ( atm.xdel[k] · sa_aerosol(j)
 *                                        + atm.ydel[k] · sa_rayleigh(j) )
 *
 *   where
 *       sa_rayleigh(j) = β₀^m + atm.beta2 · atm.xpl[j] · atm.xpl[0]
 *       β₀^m           = atm.beta0 if m == 0, else 0
 *                        (the isotropic Rayleigh term contributes only to m=0)
 *       sa_aerosol(j)  = ws.phase_fourier_m[0][j]
 *
 *   For m ≥ 3 in pure-Rayleigh setups, atm.xpl[j] = Plm[2][j] = 0
 *   (rt_legendre_compute zeros the Plm[l<m] band) so sa_rayleigh(j) = 0
 *   automatically. With NULL betal, sa_aerosol(j) = 0 too, and the entire
 *   source vanishes.
 *
 * Convention B — the (½) factor
 *
 *   The source includes a (½) prefactor via atm.ch[k] = exp(-τ_k/μ_s)/2.
 *   The matching × 0.5 lives in the linear-in-τ vertical integration step
 *   (Step 6). Their product (½)(½) = ¼ recovers the textbook prefactor
 *   ω₀ F_⊙ / (4π) for ω₀ = 1 (Rayleigh) when F_⊙ = π — see the long
 *   comment block at the top of rt_types.h (RT_F_SOLAR_PI).
 *
 *   Convention A (alternative, NOT used here): would put the entire ¼
 *   in the source and use unit weight in the integration. Mathematically
 *   equivalent. We use Convention B because it (i) keeps each routine a
 *   pure structural multiply with no free constants, (ii) lets us cross-
 *   check against 6SV at the i2-level with bit-equality, and (iii) reuses
 *   atm.ch[k] without re-multiplying.
 *
 * Solar slot poisoning
 *
 *   source_first[k][0] is set to NaN, NOT to its formal value. Index
 *   j = 0 is the solar slot (rm[0] = -μ_solar) and never participates in
 *   diffuse-field integration in any downstream loop. NaN-poisoning makes
 *   accidental reads propagate immediately so bugs surface early. Step 6
 *   (vertical integration) and Step 7 (SOS) will assert finite on the
 *   slots they read.
 *
 * Inputs
 *   atm    Already built (Step 3). atm.ch, ydel, xdel, beta0, beta2 must
 *          all be populated; atm.xpl must correspond to the requested m
 *          — caller invokes rt_legendre_compute(ws, atm, m) first.
 *   m      Fourier order (0 .. ws->l_max).
 *   ws     Workspace whose phase_fourier_m[0][j] holds the aerosol
 *          contribution; for Phase 1 (pure Rayleigh) call
 *          rt_kernel_phase_fourier(ws, m, NULL) first to zero it.
 *
 * Output
 *   source_first   Caller-allocated, flat row-major layout
 *                  source_first[k * (2*n_mu+1) + (j + n_mu)]
 *                  for k = 0..atm->n_layers, j = -n_mu..+n_mu.
 *                  Slot at j = 0 is poisoned to NaN (see above).
 *
 * Returns 0 on success, -1 on bad args.
 *
 * Thread-safety: pure function in its inputs; no globals, no allocation.
 *   Multiple threads may call concurrently with separate (atm, ws,
 *   source_first) triples.
 */
int rt_solver_primary_source(const rt_atm_t *atm, int m,
                             const rt_legendre_workspace_t *ws,
                             double *source_first);

/* ----------------------------------------------------------------------------
 * rt_solver_primary_source_pol  (Phase 2: vector Q, U)
 *
 *   Same as rt_solver_primary_source, but produces Q and U Stokes source
 *   buffers in addition to (or instead of) the I source.  Reference:
 *   6SV OSPOL.f line 320-348:
 *
 *     For m ≤ 2 (Rayleigh contributes):
 *       sb1 = γ_2 × xrl(j) × xpl(0)        (Rayleigh Q, vector)
 *       sb2 = gr(0, j)                     (Aerosol Q, from gr_pol)
 *       sc1 = γ_2 × xtl(j) × xpl(0)        (Rayleigh U, vector)
 *       sc2 = gt(0, j)                     (Aerosol U, from gt_pol)
 *     For m > 2 (Rayleigh has no l>2):
 *       sb1 = sc1 = 0;  sb2, sc2 from aerosol only.
 *
 *     Source at level k (matching OSPOL.f line 344-346):
 *       q_source[k, j] = + ch(k) × (xdel × sb2 + ydel × sb1)
 *       u_source[k, j] = - ch(k) × (xdel × sc2 + ydel × sc1)
 *
 *   The U-source carries a leading minus sign (6SV convention; matches
 *   Stokes [I, Q, U] convention with U sign tied to azimuth orientation).
 *
 *   Pre-conditions:
 *     - rt_legendre_compute_pol(ws, atm, m) must have been called
 *       (atm->xrl, xtl filled for current m; ws->rrl, rtl filled).
 *     - rt_kernel_phase_fourier_pol(ws, m, gammal_aer) must have been
 *       called (ws->gr_pol, gt_pol filled; pass NULL gammal for Phase 2
 *       Rayleigh-only pure case → zeros).
 *
 *   Layouts identical to rt_solver_primary_source's source_first.
 *   Solar slot j=0 is poisoned to NaN.
 *
 * Returns 0 on success, -1 on bad args.
 */
int rt_solver_primary_source_pol(const rt_atm_t *atm, int m,
                                 const rt_legendre_workspace_t *ws,
                                 double *source_q,
                                 double *source_u);

/* ----------------------------------------------------------------------------
 * rt_solver_integrate
 *
 *   Integrate the m-th Fourier source J_m(τ_k, μ_j) (Step 5 output) along
 *   each direction μ_j to produce the diffuse intensity I_m(τ_k, μ_j).
 *
 *   For each direction j_signed ∈ {-n_mu, …, -1, +1, …, +n_mu}
 *   (j_signed = 0 is the solar slot — skipped; the input source value
 *   there is the NaN sentinel from Step 5):
 *
 *     UPWARD (j_signed > 0)        BC I(nt, j) = 0 (black surface);
 *                                  k = nt-1 → 0
 *     DOWNWARD (j_signed < 0)      BC I(0, j) = 0 (no diffuse at TOA);
 *                                  k = 1 → nt
 *
 *   Per-layer update (LINEAR mode, default):
 *
 *     a   = (J(k_next) - J(k_now)) / Δτ
 *     b   = J(now) - a · h(now)
 *     c   = exp(-Δτ / |μ|)
 *     I_new = c · I_old +
 *             ((1 - c)(b + a·μ_signed) + a·(h(k_now) - h(k_next)·c)) · 0.5
 *
 *   The (×0.5) is Convention B's match to atm.ch[k]'s (½) prefactor
 *   in Step 5 — see RT_F_SOLAR_PI in rt_types.h.
 *
 *   In CONSTANT mode (v1 reproduction; do NOT use for production):
 *     a is forced to zero, b = J(k_now); the slope term is dropped.
 *     This intentionally reproduces v1's grid-divergence pathology so
 *     it can be visualized as a comparison figure.
 *
 *   Inputs:
 *     atm           Step-3 atmosphere; reads atm.h[k] and atm.rm[j]
 *                   (and atm.n_layers / atm.n_mu).
 *     m             Fourier order (used only for diagnostics — the
 *                   integration is m-independent).
 *     source_first  Step-5 output, layout (n_layers+1) × (2 n_mu + 1),
 *                   source_first[k * dirs + (j + n_mu)].  source_first
 *                   at j_signed = 0 must be NaN (Step 5 sentinel);
 *                   integration loops never read it.
 *     method        LINEAR (production) or CONSTANT (v1 reproduction).
 *
 *   Output:
 *     radiation     Caller-allocated, same shape as source_first.
 *                   radiation[k * dirs + (j + n_mu)] = I_m(τ_k, μ_j).
 *                   j_signed = 0 slot is set to NaN (defensive).
 *
 *   Returns 0 on success, -1 on bad args, -2 if a non-finite source
 *   is read at j_signed ≠ 0 (real bug — source corruption).
 *
 * Thread-safety: pure function in inputs; no globals, no allocation.
 */
int rt_solver_integrate(const rt_atm_t *atm, int m,
                        const double *source_first,
                        rt_integration_method_t method,
                        double *radiation);

/* ----------------------------------------------------------------------------
 * rt_solver_integrate_with_surface_bc
 *
 *   Same as rt_solver_integrate but with non-zero surface boundary condition.
 *   For each upward direction j > 0, the BC at the surface (k = nt) is taken
 *   from `surface_bc[j]` (size n_mu, indexed 1..n_mu by Fortran convention).
 *
 *   The BC array layout:
 *     surface_bc[j-1] = upward radiance at (τ_max, μ_q[j])  for j = 1..n_mu
 *
 *   When surface_bc = NULL, behavior is identical to rt_solver_integrate
 *   (BC = 0, black ocean).
 *
 *   Used by rt_solver_sos_pol_with_surface for atm-surface MS coupling.
 *   Phase 4: flat Fresnel and Cox-Munk Fresnel surfaces.
 *
 * Returns 0 on success, -1 on bad inputs, -2 on non-finite source.
 */
int rt_solver_integrate_with_surface_bc(const rt_atm_t *atm, int m,
                                         const double *source_first,
                                         const double *surface_bc,
                                         rt_integration_method_t method,
                                         double *radiation);

/* rt_solver_integrate_bcs — like rt_solver_integrate_with_surface_bc but ALSO
 * accepts a TOP downward boundary condition (level k=0, the water surface from
 * below). For each downward direction j<0, the BC at k=0 is top_down_bc[|j|-1]
 * (size n_mu, |j|=1..n_mu). top_down_bc=NULL → BC=0 (current behavior).
 * Used by the P0-A surface internal-reflection feedback prototype. */
int rt_solver_integrate_bcs(const rt_atm_t *atm, int m,
                            const double *source_first,
                            const double *surface_bc,
                            const double *top_down_bc,
                            rt_integration_method_t method,
                            double *radiation);

/* ----------------------------------------------------------------------------
 * rt_solver_sos
 *
 *   Successive Orders of Scattering iteration.  Given the primary
 *   intensity field I_1 (Step 6 output), iteratively generate higher
 *   orders I_2, I_3, ... using
 *
 *       J_n^m(τ_k, μ_j) = Σ_{j' ≠ 0}
 *           [ ydel(k) · P_R^m(j, j')  +  xdel(k) · P_aer^m(j, j') ]
 *           · I_{n-1}^m(τ_k, μ_{j'})
 *           · w(|j'|)
 *
 *   followed by vertical integration via rt_solver_integrate() (LINEAR
 *   mode hard-coded — CONSTANT mode is for v1 reproduction at the
 *   primary level only and would compound errors across SOS).
 *
 *   The phase pieces P_R^m and P_aer^m are layer-independent; the
 *   layer mix (ydel, xdel) is applied inside the source assembly.
 *
 *   Convention B: SOS source has NO explicit prefactor (matches OS.f
 *   L444 `i2 = Σ ...`); the (½) lives in rt_solver_integrate's ×0.5.
 *
 *   Inputs
 *     atm           Step-3 atmosphere; reads h, rm, gb, xdel, ydel,
 *                   beta0, beta2, xpl.  atm.xpl must correspond to m
 *                   (caller invokes rt_legendre_compute first).
 *     m             Fourier order.
 *     ws            Workspace whose phase_fourier_m[j][k] is populated
 *                   from rt_kernel_phase_fourier(ws, m, betal).
 *     primary_field Step-6 output I_1; layout (n_layers+1) × (2 n_mu+1).
 *                   Must have NaN at j=0 slot (Step 5 sentinel).
 *     opts          NULL → defaults (max_iter=20, tol=1e-7,
 *                   acceleration=PLAIN, save_orders=0).
 *
 *   Outputs
 *     total_field   Caller-allocated, same shape as primary_field.
 *                   = Σ_{n=1..N} I_n.  j=0 NaN-poisoned.
 *     I_per_order   Optional. If non-NULL, layout
 *                   (max_iterations+1) × (n_layers+1) × (2 n_mu+1).
 *                   I_per_order[n] = I_{n+1}, with index 0 = primary.
 *                   Caller allocates max_iterations+1 buffer rows.
 *     result        Caller-allocated; receives n_orders_used,
 *                   converged flag, final_residual.
 *
 *   Returns 0 on success, -1 on bad args, -2 on NaN propagation.
 *
 *   Thread-safety: pure function in inputs; allocates two scratch
 *   buffers internally (J_curr + I_curr), freed before return.
 *   Multiple threads may call concurrently with separate output
 *   buffers.
 */
int rt_solver_sos(const rt_atm_t *atm, int m,
                  const rt_legendre_workspace_t *ws,
                  const double *primary_field,
                  const rt_solver_sos_options_t *opts,
                  double *total_field,
                  double *I_per_order,                       /* may be NULL */
                  rt_solver_sos_result_t *result);

/* Vector RT extension of rt_solver_sos (Phase 2).  Iterates Successive
 * Orders of Scattering for I, Q, U Stokes components simultaneously.
 * Each iteration:
 *   J_I, J_Q, J_U  ← rt_sos_operator_apply_vector(atm, m, ws, I_prev, Q_prev, U_prev)
 *   I_curr/Q_curr/U_curr ← rt_solver_integrate × 3
 *   total_I/Q/U   += I_curr / Q_curr / U_curr
 * Convergence residual: max(||I_curr||, ||Q_curr||, ||U_curr||) / max-totals.
 * I_per_order_* omitted (downstream Fourier reconstruction only needs totals).
 * Returns 0 on success; same error codes as rt_solver_sos.
 */
int rt_solver_sos_pol(const rt_atm_t *atm, int m,
                      const rt_legendre_workspace_t *ws,
                      const double *primary_I,
                      const double *primary_Q,
                      const double *primary_U,
                      const rt_solver_sos_options_t *opts,
                      double *total_I,
                      double *total_Q,
                      double *total_U,
                      rt_solver_sos_result_t *result);

/* rt_solver_sos_pol_intrefl — P0-A prototype. Same as rt_solver_sos_pol, but at
 * each SOS order injects the water-side surface internal reflection as a TOP
 * downward boundary condition via the FULL Fresnel Mueller matrix (2026-06-02):
 *   [I,Q,U]_dn^m(k=0,-mu_j) = (-1)^m * M_Rww(mu_j) * [I,Q,U]_up^m(k=0,+mu_j),
 * using the upwelling at the surface from the previous order. Summed over
 * orders this converges to the self-consistent reflecting-boundary field
 * (flat surface; TIR for theta_w>theta_c via M_Rww=identity). The (-1)^m and
 * block structure [[M00,M01,0],[M10,M11,0],[0,0,M22]] mirror the validated
 * rt_solver_sos_pol_with_surface specular BC.
 * Rww_M (size n_mu*9) = per-node row-major 3x3 internal-reflection Mueller
 * M_Rww(mu_q[j]) from rt_air_water_R_ww. Env-gated upstream. */
int rt_solver_sos_pol_intrefl(const rt_atm_t *atm, int m,
                              const rt_legendre_workspace_t *ws,
                              const double *primary_I,
                              const double *primary_Q,
                              const double *primary_U,
                              const double *Rww_M,
                              const rt_solver_sos_options_t *opts,
                              double *total_I,
                              double *total_Q,
                              double *total_U,
                              rt_solver_sos_result_t *result);
/* rt_solver_sos_pol_intrefl_rough - B2 (2026-06-03). Rough-surface (Cox-Munk)
 * water-side internal reflection. Same as rt_solver_sos_pol_intrefl but the top
 * downward BC reflects the previous-order upwelling through the angle-coupling
 * m-mode kernel Rww_K (surface_R_ww_coxmunk_fourier_kernel), not the per-node
 * specular Mueller. Contraction: sum_k (2pi mu_k w_k) Rww_K[j,k] up_k (no
 * (-1)^m; kernel carries azimuth). Rww_K size n_mu*n_mu*9, m-mode. wind > 0. */
int rt_solver_sos_pol_intrefl_rough(const rt_atm_t *atm, int m,
                              const rt_legendre_workspace_t *ws,
                              const double *primary_I,
                              const double *primary_Q,
                              const double *primary_U,
                              const double *Rww_K,
                              const rt_solver_sos_options_t *opts,
                              double *total_I,
                              double *total_Q,
                              double *total_U,
                              rt_solver_sos_result_t *result);

/* ----------------------------------------------------------------------------
 * rt_solver_sos_pol_with_surface
 *
 *   Same as rt_solver_sos_pol, but with vector surface boundary coupling.
 *   At each iteration, the upward radiation BC at the surface is computed
 *   from the previous-order downward radiation via the surface Mueller
 *   matrix M^flat (Phase A) or the Cox-Munk Fourier kernel (Phase B).
 *
 *   For flat Fresnel (wind_speed=0):
 *     [I_up^m]              [M[0][0] M[0][1]    0  ] [I_dn^m]
 *     [Q_up^m] = (-1)^m  ·  [M[1][0] M[1][1]    0  ] [Q_dn^m]
 *     [U_up^m]              [   0       0    M[2][2]] [U_dn^m]
 *
 *   M^flat = surface_flat_fresnel_matrix_raw(μ_o, n_water, q_convention).
 *
 *   Cox-Munk (wind_speed > 0): kernel_R_m_per_mu computed via numerical
 *   φ-quadrature integration of surface_R_coxmunk_trig (Phase B; not yet
 *   implemented — currently treats wind_speed > 0 as flat).
 *
 *   Returns 0 on success.
 */
int rt_solver_sos_pol_with_surface(const rt_atm_t *atm, int m,
                                    const rt_legendre_workspace_t *ws,
                                    const double *primary_I,
                                    const double *primary_Q,
                                    const double *primary_U,
                                    const rt_solver_sos_options_t *opts,
                                    /* Surface params */
                                    rt_surface_kind_t surface_kind,
                                    double n_water,
                                    double wind_speed,
                                    int    sigma_type,
                                    int    q_convention,
                                    /* Outputs */
                                    double *total_I,
                                    double *total_Q,
                                    double *total_U,
                                    rt_solver_sos_result_t *result);

/* Reconstruct I/Q-like cosine modes at a public relative azimuth.
 *
 * OCRT stores only non-negative azimuthal modes.  Reality of the Stokes field
 * folds the negative modes into a factor two for m>0.  The public relative-
 * azimuth origin differs from the stored propagation azimuth by pi, giving
 *
 *   X(phi) = X_0 + 2 sum_{m=1}^M X_m cos[m(phi+pi)].
 */
double rt_solver_reconstruct_phi(const double *I_total_per_m_at_view,
                                 int m_max,
                                 double delta_phi_rad);

/* Reconstruct the sine-series U component under the same convention:
 *
 *   U(phi) = 2 sum_{m=1}^M U_m sin[m(phi+pi)].
 *
 * Unlike the cosine series, reversing the sign of phi reverses U.
 */
double rt_solver_reconstruct_phi_sin(const double *U_total_per_m_at_view,
                                     int m_max,
                                     double delta_phi_rad);

/* ----------------------------------------------------------------------------
 * rt_solver_reflectance_from_intensity
 *
 *   Convert TOA upwelling intensity to TOA reflectance.
 *
 *       ρ = I_TOA / μ_solar         (Convention 2: F_⊙ = π absorbed in source)
 *
 *   The π factor is *not* present in the conversion: under Convention 2
 *   (used throughout v2 -- see RT_F_SOLAR_PI in rt_types.h), the F_⊙ = π
 *   normalization is folded into the source's (¼) prefactor (ω₀F_⊙/(4π) =
 *   ¼).  The matching reflectance definition is therefore
 *   `ρ = π·L/(μ_s·F_⊙) = π·L/(μ_s·π) = L/μ_s`.
 *
 *   This matches OS.f's own commented reflectance formula at L670:
 *   `c  write(6,*) 'reflectance ', xl(mu,1)/xmus`.
 *
 *   v1 used Convention 1 (F_⊙ = 1, source = (1/4π)·P·exp, so
 *   ρ = π·I/μ_s) -- different normalization but same final ρ.
 *
 *   Returns 0 if μ_solar ≤ 0 (defensive against ill-formed inputs).
 *
 * Thread-safety: pure function (static inline, header-only).
 */
static inline double rt_solver_reflectance_from_intensity(double I_TOA,
                                                          double mu_solar) {
    return (mu_solar > 0.0) ? (I_TOA / mu_solar) : 0.0;
}

#endif /* OCRT_V2_RT_SOLVER_H */
