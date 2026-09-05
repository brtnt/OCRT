/* ============================================================================
 * rt_water_rt.h — In-water vector SOS solver (Phase B.3)
 *
 * Phase B.3 implementation (2026-05-22).
 *
 * Strategy: REUSE existing atmospheric SOS engine (rt_solver_sos_pol) by
 * populating rt_atm_t with in-water values. This guarantees numerical
 * consistency between atmospheric and in-water paths, with zero code
 * duplication.
 *
 * Mapping (atmospheric → in-water):
 *   tau_total      ← (a_w + b_w) × z_max         (extinction OD over depth)
 *   mu_sun         ← Snell-refracted in-water solar mu
 *   depol          ← δ_w = 0.039 (Zhang 2009 const)
 *   beta2, ...     ← Rayleigh-like (δ_w 대입, atmospheric form 동일)
 *   xdel[k]        ← ω_w = b_w / (a_w + b_w)     (SSA via aerosol slot)
 *   ydel[k]        ← 0                            (Rayleigh OFF, water is aerosol-like)
 *   betal_aer      ← pure water Rayleigh-like Legendre with δ_w
 *   L_max          ← 2 (Rayleigh-like, no l>2)
 *   aerosol_active ← 1 (use aerosol slot for water phase)
 *
 * Result extraction:
 *   total_I/Q/U[k=0, j_pos] = just-below-surface upwelling Stokes for mode m
 *   Sum over m → directional radiance at view
 *   Hemispheric integral → Ed(0-), Eu(0-)
 *   Kd, Ku from log-derivative of Ed, Eu near surface
 *
 * Phase A integration:
 *   rt_water_iop_aw, _bw, _aw_T (NASA z09 + Röttgers ψ_T) already available.
 *   bb_w = 0.5 × b_w (Zhang 2009, mathematical exact).
 *   δ_w = 0.039 (Zhang 2009).
 *
 * Reference:
 *   - Atmospheric SOS: rt_solver_sos_pol (rt_solver.c L1412+)
 *   - Pure water IOP: rt_water_iop.c (Phase A, v1.02)
 *   - Mobley 1994, "Light and Water", Ch. 8-9
 * ============================================================================ */

#ifndef OCRT_RT_WATER_RT_H
#define OCRT_RT_WATER_RT_H

#include <stddef.h>
#include "rt_water_iop.h"  /* rt_water_iop_lut_t, rt_water_iop_psi_T_lut_t */

#ifdef __cplusplus
extern "C" {
#endif

/* v1.09 #19: reset the (threadprivate) full-grid replay cache. */
void rt_water_grid_cache_reset(void);

/* v1.09 #21: one-shot env snapshot (call before OpenMP regions; idempotent). */
void rt_water_env_init(void);

/* (forward declarations removed; included from rt_water_iop.h directly) */

/* ===========================================================================
 * Solver options and result
 * =========================================================================== */

typedef struct {
    /* SOS iteration control */
    int    max_iterations;          /* default 20 (matches atmospheric SOS) */
    double tolerance;               /* default 1e-7 */

    /* Vertical discretization */
    int    n_layers_water;          /* <=0: auto from final tau_max/layer_dtau_target; >0: debug-only override. */
    double tau_max_target;          /* <=0: auto start/base optical depth; >0: debug-only base override. */
    double max_tau_max_target;      /* production maximum optical-depth cap; debug-only override. */
    double max_z_max_m;             /* production maximum physical water depth [m]; debug-only override. */
    double depth_bottom_tol;        /* adaptive depth target for exp(-transport_ext*z). */
    double layer_dtau_target;       /* auto vertical grid target Δτ per layer; debug-only override. */

    /* Direction quadrature */
    int    n_mu_water;              /* default 48 (raised from 24, 2026-06-30); >= 48
                                     * recommended for in-water RT incl. LUT mode. n_mu=32
                                     * under-converged ~1.5% for forward-peaked high-omega
                                     * particles. Floor rises with omega. See rt_water_rt.c. */
    /* v1.09 commit #20 (S-007 root fix): air-side vza list appended as
     * zero-weight in-water nodes (OSOAA UserAngFile-equivalent).  The per-mode
     * view extraction then reads these nodes DIRECTLY instead of linearly
     * interpolating between Gauss nodes.  NULL/0 = off (legacy behavior). */
    const double *view_vza_deg_list;
    int    n_view_vza;
    int    view_as_node;            /* default 0: continuous reconstruction from GL nodes.
                                     * 1: add a zero-weight output node at the requested
                                     *    water-side view direction when it lies outside the
                                     *    GL range. Diagnostic only; CCRR v0.5.0 reference
                                     *    uses continuous_reconstruction. */

    /* Fourier modes */
    int    m_max_water;             /* default 2 (Rayleigh-like has no l>2 contribution) */

    /* Q convention */
    int    q_convention;            /* default 1 (Mishchenko, matches OCRT atm) */

    /* CDOM 2-parameter absorption (B.5, 2026-05-23):
     *   a_CDOM(λ) = a_CDOM_440 · exp(-S_CDOM · (λ - λ_ref))
     * CDOM 은 흡광만 추가하며 b_w, Mueller matrix 에는 영향이 없다.
     * Default: a_CDOM_440=0 (off). */
    double a_cdom_440_m_inv;        /* a_CDOM(λ_ref) [m⁻¹], default 0 */
    double S_cdom_nm_inv;           /* spectral slope [nm⁻¹], default 0.014 */
    double cdom_ref_lambda_nm;      /* reference λ [nm], default 440 */

    /* Public branch provenance; rt_water_input_mode_t.  The common RT solver
     * still consumes the derived ccrr_mode/fixed_bulk_iop_mode flags below. */
    int    water_input_mode;

    /* Water constituent mode (pre-solver only):
     * CHL/TSM/aDOM -> OCRT RT IOP and vector-phase inputs. */
    int    ccrr_mode;               /* derived: 0 = native pure water, 1 = constituent adapter active */
    int    water_constituent_model; /* rt_water_constituent_model_t; OCRT default */
    double ccrr_chl_mg_m3;          /* chlorophyll concentration [mg m^-3] */
    /* 0 = phytoplankton absorption only.  Species phase generation remains
     * available as an offline/API facility, but its RT scattering path is
     * disabled until forward-peak truncation is implemented consistently. */
    int    organic_phyto_scattering;
    double ccrr_min_g_m3;           /* TSM dry-weight concentration [g m^-3] */
    int    tsm_species;               /* ahn_species_t; 0=red_clay default */
    int    organic_phyto_group;       /* organic_phyto_group_t; micro default */
    double detritus_a440_m_inv;       /* optional a_d(440) [m^-1] */
    double detritus_slope_nm_inv;     /* default 0.0109 nm^-1 */
    const char *ccrr_phase_moments_path; /* optional IOCCG21/CCRR particle moment CSV; pure water remains OCRT vector */
    const char *ccrr_particle_phase_lut_path; /* positive high-resolution particle P11 LUT; replaces finite moments */
    const char *ccrr_particle_phase_case_id;  /* optional row selector for particle P11 LUT */
    double ccrr_particle_phase_wavelength_nm; /* optional row selector for particle P11 LUT */
    int    ccrr_particle_phase_lmax;          /* retained order for formal delta-M, default 10 */
    int    ccrr_particle_phase_nphi;          /* azimuth quadrature for direct Fourier kernel, default 720 */

    /* Fixed-total-IOP bulk-phase validation mode for HydroLight Mroot/native
     * nadir tables.  It uses a_total,b_total,bb_total directly and treats the
     * full scattering budget as one scalar bulk phase matched to bb/b. */
    int    fixed_bulk_iop_mode;
    double fixed_a_total_m_inv;
    double fixed_b_total_m_inv;
    double fixed_bb_total_m_inv;
    int    fixed_bulk_lmax;
    int    fixed_bulk_phase_model; /* 0=delta2bb isotropic, 1=HG direct Fourier, 2=tabulated P11 direct Fourier, 3=tabulated P11 with formal delta-M, 4=analytic scalar FF with formal delta-M */
    int    fixed_bulk_phase_nphi;
    double fixed_bulk_ff_n;        /* FF relative refractive index, default 1.18 */
    double fixed_bulk_ff_mu;       /* optional FF Junge slope; <=0 solves from bb/b */
    const char *fixed_bulk_phase_lut_path;
    const char *water_mie_phase_path; /* #0: optional .mie vector phase (P11/P12/P33) for bulk-water particle; bulk IOP from fixed-bulk */
    int water_mie_moment_mode;       /* 0=dense theta/PCHIP, 1=Gauss/μ-cubic (default) */
    int water_mie_truncation_mode;   /* branch-prefixed --ocrt/--iop-mie-truncation; 0=off, 1=delta truncation */
    int water_mie_ss_mode;           /* branch-prefixed --ocrt/--iop-mie-ss-mode: 0=none, 1=IMS, 2=NT-TMS */
    int water_mie_moment_n_mu;       /* Gauss projection nodes; default 400 */
    const char *fixed_bulk_phase_case_id;
    double fixed_bulk_phase_wavelength_nm;

    /* Cox-Munk water->air transmission T_wa branch (option A: full facet
     * integration over the in-water radiance field).
     *   wind_speed <= 0 : flat Fresnel T_wa (single refracted view direction)
     *   wind_speed  > 0 : surface_T_coxmunk_trig BTDF integrated over in-water
     *                     mu-nodes x azimuth (mixes multiple in-water dirs). */
    double wind_speed;          /* [m/s], default 0 (flat) */
    int    cox_munk_sigma_type; /* slope-variance model (matches atm side), default 1 */

    /* Particle (FF) phase-function truncation / kernel method:
     *   0 = value  (path B, default): OSOAA-style angular forward-peak cap on
     *               the analytic P11 VALUES + value-based azimuth Fourier kernel.
     *               Positive radiance at any nphi (no Legendre Gibbs).
     *   1 = moment (path A, SHELL): Wiscombe delta-M scaled Legendre moments +
     *               addition-theorem kernel. Backscatter-negative at low OS_NB;
     *               placeholder only (TODO: cap -> moments). */
    int    water_phase_kernel;

    /* [FIX-SKY-EDLU 2026-06-28] Sub-cone skylight equivalent-beam support.
     * When > 0, the in-water solar-beam cosine is taken DIRECTLY as this value
     * (the in-water zenith cosine), bypassing the air->water Snell refraction of
     * sza_deg_air.  This lets the diffuse-skylight Lu superposition drive equivalent
     * beams for in-water directions below the flat critical angle (mu_w < mu_crit),
     * which the rough sea surface transmits into (sub-cone leak) but which have no
     * real air incidence angle.  F_sun is then the in-water beam flux directly.
     * Default 0 (disabled -> standard air-refraction path). */
    double mu_sun_water_override;
    /* v1.10 B-0c.2b: per-m diffuse TOP incident (F_sun=pi convention),
     * sampled on ext_top_mu[] (the coupling GL grid).  The water solver
     * maps its ring columns by |dmu|<1e-12 (specials -> 0).  m-range:
     * 0..ext_top_m_max.  NULL => feature off (equivalent beams path). */
    const double *ext_top_I;
    const double *ext_top_Q;
    const double *ext_top_U;
    const double *ext_top_mu;
    int           ext_top_n;
    int           ext_top_m_max;
    /* v1.10 D3 (2026-07-10): the full-grid cache KEY does NOT include the
     * ext_top injection content (verified by key-field inventory), so
     * operator-assembly driver solves - which vary ONLY ext_top - would all
     * collide onto one entry (observed: 153 basis columns replayed from one
     * solve in 0.1 s, an all-zero operator).  The env switch cannot help at
     * run time (init-once snapshot).  Driver solves set this to 1 to skip
     * BOTH the cache lookup and the cache store; production leaves 0 -
     * behavior unchanged. */
    int           bypass_grid_cache;
} rt_water_rt_options_t;

/* Return options with defaults */
rt_water_rt_options_t rt_water_rt_options_default(void);

typedef struct {
    /* Iteration meta (across all m modes) */
    int    max_orders_used;
    int    all_converged;            /* 1 if all m modes converged */
    double max_residual;             /* worst across m modes */

    /* IOP used (for diagnostics, at requested λ and T) */
    double lambda_nm_used;
    double T_water_C_used;
    double a_w_used;                 /* m⁻¹ */
    double b_w_used;                 /* m⁻¹ */
    double bb_w_used;                /* m⁻¹, = 0.5 × b_w */
    double a_cdom_used;              /* m⁻¹ — a_CDOM(λ) at solved λ (B.5) */
    double a_pig_used;               /* m⁻¹ — CCRR pigment absorption */
    double a_chl_used;               /* m⁻¹ — phytoplankton absorption only */
    double b_pig_used;               /* m⁻¹ — CCRR pigment scattering */
    double bb_pig_used;              /* m⁻¹ — CCRR pigment backscatter */
    double a_min_used;               /* m⁻¹ — Ahn mineral/TSM absorption */
    double b_min_used;               /* m⁻¹ — Ahn mineral/TSM scattering */
    double bb_min_used;              /* m⁻¹ — Ahn mineral/TSM backscatter */
    double a_total_used;             /* total absorption */
    double b_total_used;             /* total scattering */
    double bb_total_used;            /* total backscatter */
    double omega_used;               /* SSA = b_total / (a_total + b_total) */
    double tau_max_used;             /* dimensionless */
    double z_max_used;               /* m (= tau_max / (a_w + b_w)) */

    /* Snell refraction */
    double mu_sun_air;               /* cos(SZA_air) */
    double mu_sun_water;             /* in-water solar cos, post-Snell */

    /* Just-below-surface upwelling Stokes at view direction (μ_v_water).
     * mu_v_water = Snell-refracted in-water view cos.
     * If view direction TIR-trapped (impossible for atm-side input view since
     * mu_v_air > 0 always refracts into water), value is 0 by convention.
     * These are summed over m=0..m_max_water at view azimuth raa. */
    double I_0minus_view;
    double Q_0minus_view;
    double U_0minus_view;

    /* Just-ABOVE-surface upwelling Stokes at view direction (μ_v_air).
     * Obtained by applying water→air radiance Mueller transmission
     * (rt_air_water_T_wa) to the just-below-surface Stokes.
     * Includes the 1/n_w² radiance reduction (n² law inverse). */
    double I_0plus_view;
    double Q_0plus_view;
    double U_0plus_view;

    /* Hemispheric integrals at z=0- (m=0 mode only — physically azimuth-avg)
     * Ed_0minus_water = ∫_{μ<0} I^{m=0}(μ) · |μ| · 2π dμ              [W·m⁻²·nm⁻¹]
     * Eu_0minus_water = ∫_{μ>0} I^{m=0}(μ) ·  μ  · 2π dμ              [W·m⁻²·nm⁻¹]
     */
    double Ed_0minus_water;
    double Eu_0minus_water;

    /* Just-above-surface downwelling irradiance.
     * Simple no-atmosphere case (B.3): Ed_0plus_air = F_sun × μ_sun_air
     * Full atmosphere-ocean coupling (B.4): includes atmospheric Ed at BOA.
     */
    double Ed_0plus_air;

    /* Remote-sensing reflectance ratios (lowercase r_rs = underwater,
     * uppercase R_rs = above-water; user convention). */
    double r_rs_0minus;     /* = I_0minus_view / Ed_0minus_water */
    double r_rs_0minus_Q;   /* = Q_0minus_view / Ed_0minus_water */
    double r_rs_0minus_U;   /* = U_0minus_view / Ed_0minus_water */
    double R_rs_0plus;      /* = I_0plus_view  / Ed_0plus_air    */
    double R_rs_0plus_Q;    /* = Q_0plus_view  / Ed_0plus_air    */
    double R_rs_0plus_U;    /* = U_0plus_view  / Ed_0plus_air    */

    /* Diffuse attenuation at z=0- (m⁻¹) */
    double Kd_0minus;
    double Ku_0minus;

    /* Internal near-surface diagnostics used by the atmosphere-ocean
     * orchestrator to reconstruct Kd(0-) with the unscattered transmitted
     * skylight included at both z=0 and the first water level.  These values
     * intentionally exclude that unscattered skylight term; the ext_top
     * scattering response is already present in Ed_level1_water. */
    double Ed_level1_water;
    double tau_level1_used;
    double z_level1_used;

    /* ---- Optional output (B.5 up-coupling): z=0⁻ upwelling Fourier field ----
     * Per-mode upwelling radiance at every positive water μ-node at z=0⁻, for
     * the caller's water→air up-transmission (build the 0⁺ field at all atm
     * μ-nodes, then inject into the atm SOS bottom-source). This is the SAME
     * field the in-water solver already builds for the Cox-Munk T_wa BTDF
     * (I_m_node), here exposed to the caller.
     *
     * Memory: CALLER pre-allocates I/Q/U_up_per_m (each size
     * (m_max_water+1)*n_mu_water) and mu_water_pos (size n_mu_water), and frees
     * them after use. If any pointer is NULL the solver skips the copy entirely
     * (no behaviour change → backward compatible with `= {0}` callers).
     * Indexing: [m*n_mu_water + (k-1)] for mode m=0..m_max_filled, node
     * k=1..n_mu_water_filled. */
    double *I_up_per_m;        /* may be NULL */
    double *Q_up_per_m;        /* may be NULL */
    double *U_up_per_m;        /* may be NULL */
    double *mu_water_pos;      /* may be NULL; positive water μ nodes, [k-1] */
    double *w_water_pos;       /* v1.10 B-0b.2: ring quadrature weights aligned
                                * with mu_water_pos (0 for zero-weight specials) */
    int     n_mu_water_filled; /* set by solver: number of positive μ nodes */
    int     m_max_filled;      /* set by solver: highest mode index (modes 0..m_max_filled) */

    /* Optional Fourier samples at the requested view direction.  The caller
     * allocates each array with at least (m_max_water+1) entries.  The below-
     * surface arrays are copied from the already solved water field; the
     * above-surface arrays are produced by the same linear T_wa coupling used
     * by the production output.  They are used only by the native coupled LUT
     * orchestration and add no SOS iteration. */
    double *I_view_per_m;
    double *Q_view_per_m;
    double *U_view_per_m;
    double *I_air_view_per_m;
    double *Q_air_view_per_m;
    double *U_air_view_per_m;
    int     view_m_max_capacity;
    int     view_m_max_filled;
    int     view_modes_exact;
} rt_water_rt_result_t;

/* ===========================================================================
 * Main entry — pure water in-water SOS
 * ===========================================================================
 *
 * Inputs:
 *   sza_deg_air    : solar zenith angle in air (degrees)
 *   vza_deg_air    : view zenith angle in air (degrees, observed direction)
 *   raa_deg        : relative azimuth (degrees)
 *   lambda_nm      : wavelength
 *   T_water_C      : water temperature (°C); Röttgers ψ_T auto-applied if ≠ 20
 *   S_water_gkg    : water salinity (g/kg) — informational only in B.3
 *   n_water        : water refractive index (e.g., 1.34)
 *   F_sun          : solar irradiance at top of water column (= T_aw_radiance · π
 *                    × F_TOA / mu_sun_air, but unit-convention-wise, F_sun=1
 *                    gives normalized output; caller multiplies for SI units)
 *   aw_lut         : NASA z09 a_w LUT (must be loaded; not auto-loaded here)
 *   psi_T_lut      : Röttgers ψ_T LUT (NULL → ψ_T=0 fallback)
 *   opts           : solver options (NULL → defaults)
 *
 * Outputs:
 *   result         : filled by callee
 *
 * Returns 0 on success, negative on error:
 *   -1 invalid input
 *   -2 IOP computation failure (a_w or b_w not retrievable)
 *   -3 SOS divergence in any m mode
 *   -4 memory allocation failure
 *
 * NOTE: This is the *in-water-only* solver. Air-water coupling (T_aw, T_wa)
 * is NOT applied here. Caller is responsible for:
 *   1. Computing in-water boundary condition (BOA atmospheric field → T_aw)
 *   2. Passing the resulting F_sun (or full Stokes BC)
 *   3. Applying T_wa to result fields to obtain above-water Lu_wl(0+)
 *
 * This separation keeps the in-water solver testable in isolation.
 * Full coupled atmosphere-ocean entry will live in rt_solver.c (B.4).
 */
int rt_water_rt_sos_pure(double sza_deg_air, double vza_deg_air, double raa_deg,
                          double lambda_nm,
                          double T_water_C, double S_water_gkg,
                          double n_water,
                          double F_sun,
                          const rt_water_iop_lut_t* aw_lut,
                          const rt_water_iop_psi_T_lut_t* psi_T_lut,
                          const rt_water_rt_options_t* opts,
                          rt_water_rt_result_t* result);

/* ===========================================================================
 * Diagnostic / single-scattering analytic reference
 * =========================================================================== */

/* Single-scattering analytic Lu(0-) at view direction.
 *   For pure water Rayleigh-like phase, SS approximation:
 *     Lu(0-, μ_v_water) = (ω_w · F_sun / (4π)) · P11(scattering_angle)
 *                       × (1 - exp(-τ_max × (1/μ_v_w + 1/μ_sun_w)))
 *                       / (1 + μ_v_w/μ_sun_w)
 *   Equivalent to the integrated source over [0, τ_max] in direction μ_v.
 *
 *   Useful for verifying SOS solver in the single-scattering limit
 *   (where multiple scattering contribution is negligible — strong absorption
 *   regime, e.g. red/NIR where a_w is large and ω is small).
 *
 *   theta_scat_deg : scattering angle between sun and view directions in water
 *                    (caller responsibility — compute from refracted angles)
 */
double rt_water_rt_single_scatter_analytic(double mu_sun_water,
                                            double mu_view_water,
                                            double theta_scat_deg,
                                            double tau_max,
                                            double omega_w,
                                            double F_sun);

#ifdef __cplusplus
}
#endif

#endif /* OCRT_RT_WATER_RT_H */
