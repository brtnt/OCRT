#ifndef OCRT_V2_RT_TYPES_H
#define OCRT_V2_RT_TYPES_H

#include <stddef.h>

/* ============================================================================
 *  OC-RT solar normalization convention  (Convention 2 / "Convention B")
 * ============================================================================
 *
 *  The whole pipeline uses F_solar = π for the TOA solar irradiance.
 *  The standard single-scattering source
 *
 *      J₁ = (ω₀ F_⊙ / 4π) · P(cos Θ_sv) · exp(-τ / μ_solar)
 *
 *  is implemented as the product of two ½-factors:
 *
 *      (a) atm.ch[k]  = exp(-τ_k / μ_solar) / 2   (Step 3, rt_atm)
 *      (b) × 0.5      in the linear-in-τ vertical integration  (Step 6)
 *
 *  Their product (½)(½) = ¼ recovers ω₀ F_⊙ / (4π) when ω₀ = 1 (Rayleigh)
 *  and F_⊙ = π. The split keeps each step a pure structural multiply, no
 *  free constants floating around.
 *
 *  Consequence — TOA reflectance under Convention 2 is
 *
 *      ρ(μ_solar, μ_view, φ) = I_TOA(μ_view, φ) / μ_solar
 *
 *  i.e. NO factor of π in the conversion.  This follows from the
 *  textbook bidirectional reflectance definition
 *      ρ = π · L / (μ_s · F_⊙)
 *  with F_⊙ = π absorbed already in the source's (¼) prefactor.
 *  This matches OS.f's own commented reference formula at L670:
 *      c   write(6,*) 'reflectance ', xl(mu,1)/xmus
 *
 *  Compare with Convention 1 (used by v1 baseline, and many textbooks):
 *      F_⊙ = 1; source = (1/4π)·P·exp(...);  ρ = π·I/μ_s.
 *  Both conventions yield the same final ρ for the same physical
 *  scenario; OCRT uses Convention 2 so the source and formal-solution
 *  operators each carry one explicit factor of one half.
 *
 *  Changing the convention (e.g. F_⊙ = 1) requires editing (a), (b),
 *  AND the reflectance formula in lockstep.  Don't change one without
 *  the others.
 *
 *  RT_F_SOLAR_PI is a documentation anchor — never multiplied directly into
 *  any expression (the constant is folded into ch[k] and the integration
 *  weight).  Code that derives reflectance from intensity should reference
 *  this comment and use rt_solver_reflectance_from_intensity().
 */
#define RT_F_SOLAR_PI  3.14159265358979323846

/* Phase 1: scalar Rayleigh-only.
 * Phase 2 will add Stokes Q,U fields here.
 * Phase 3 will add aerosol layer description.
 */

typedef enum {
    RT_SURFACE_BLACK = 0,
    RT_SURFACE_LAMBERT = 1,    /* Phase 4 */
    RT_SURFACE_BLACK_FRESNEL_OCEAN = 2, /* rough Fresnel interface over a black ocean; no water-leaving radiance */
    RT_SURFACE_FLAT    = 3,    /* AF1982 frs.f exact flat Fresnel */
    RT_SURFACE_OCEAN   = 4     /* Phase A: black ocean (flat Fresnel) + in-water RT */
} rt_surface_kind_t;

/* Deprecated source-compatibility alias. New code must use the physical
 * boundary-condition name RT_SURFACE_BLACK_FRESNEL_OCEAN. */
#define RT_SURFACE_COXMUNK RT_SURFACE_BLACK_FRESNEL_OCEAN

/* SOS (Successive Orders of Scattering) acceleration method.
 *
 *   PLAIN (default):
 *     Plain accumulation.  I_total = I_1 + I_2 + ... + I_N.  Stop
 *     when ‖I_n‖_∞ / ‖I_total‖_∞ < tolerance, or n == max_iterations.
 *     Used by Phase 1 verification for clarity.
 *
 *   GEOMETRIC (optional, future):
 *     Geometric-series tail extrapolation of the Neumann-series tail.
 *     Detects geometric decay I_n ≈ r · I_{n-1} and adds the
 *     extrapolated tail Σ_{n>N} I_n ≈ I_N · r/(1-r) to the running
 *     sum.  Saves 2-4 iterations for typical Rayleigh.  Reserved for
 *     production after Phase 1 verification establishes the plain-mode
 *     baseline.
 */
typedef enum {
    RT_SOS_ACCELERATION_PLAIN     = 0,
    RT_SOS_ACCELERATION_GEOMETRIC = 1
} rt_sos_acceleration_t;

/* Iteration control for rt_solver_sos().  Pass NULL to use defaults. */
typedef struct {
    int                   max_iterations;   /* default 20 */
    double                tolerance;        /* default 1e-7 (Phase 1 PASS = 5e-4) */
    rt_sos_acceleration_t acceleration;     /* default PLAIN */
    int                   save_orders;      /* 0 = discard per-order;
                                             * 1 = keep all I_n (caller buffer) */
    /* C2 (water->air up-coupling): external upward bottom-source = the
     * water-leaving 0+ Fourier field for the CURRENT mode, at the solver's
     * atm mu-nodes (indexed by mu_pos[k_o], k_o=0..n_mu-1). When all three are
     * non-NULL the solver attenuates this upward BOA field to every level and
     * adds it to BOTH total (direct upward transmission) and the order-1
     * iteration seed, so the SOS scatters it -> full multiple-scattering
     * propagation (mirror of the sunglint surface-reflection block). NULL =>
     * no-op (backward compatible; designated/positional initializers zero-fill). */
    const double         *ext_bottom_I;     /* may be NULL */
    const double         *ext_bottom_Q;     /* may be NULL */
    const double         *ext_bottom_U;     /* may be NULL */
    /* When non-zero, skip the solar source blocks (FRESNEL-diff, sunglint,
     * lower-boundary primary reflection) so ONLY ext_bottom drives the field
     * (water-leaving propagation pass; caller also passes primary=0). The SOS
     * iteration machinery (scattering + surface BC) still runs. */
    int                   bottom_source_only;
} rt_solver_sos_options_t;

/* Diagnostic output from rt_solver_sos(). */
typedef struct {
    int    n_orders_used;      /* total orders n=1..N, including primary */
    int    converged;          /* 1 if residual < tolerance; 0 if max hit */
    double final_residual;     /* ‖I_n_last‖_∞ / ‖I_total‖_∞ at last step */
} rt_solver_sos_result_t;

static inline rt_solver_sos_options_t rt_solver_sos_options_default(void) {
    rt_solver_sos_options_t o = { .max_iterations = 20,
                                  .tolerance      = 1.0e-7,
                                  .acceleration   = RT_SOS_ACCELERATION_PLAIN,
                                  .save_orders    = 0 };
    return o;
}

/* Vertical-integration source-interpolation method.
 *
 *   LINEAR (default, production):
 *     Per-layer source approximated as J(τ) = a τ + b (linear in τ).
 *     Integrand integrated analytically across each layer.  Result is
 *     grid-converged: nt = 10, 20, 50, 100 give the same answer to
 *     machine precision.  The analytical integral follows directly from
 *     the integrating-factor solution of the linear-source VRTE.
 *
 *   CONSTANT (v1 reproduction; do NOT use for production):
 *     Per-layer source held constant at the upper boundary value
 *     (J(τ) ≡ J(level k)).  This was v1's pathology — the missing
 *     slope term gave O(Δτ) bias per layer that did NOT vanish
 *     with refinement.  Exposed as a flag purely so we can
 *     visually demonstrate the v1 grid-divergence bug for
 *     paper-figure / talk material.
 */
typedef enum {
    RT_INTEGRATION_METHOD_LINEAR   = 0,
    RT_INTEGRATION_METHOD_CONSTANT = 1
} rt_integration_method_t;

/* Rayleigh optical-depth model.
 *
 * OCRT production uses the first-principles Bodhaine et al. (1999)
 * formulation with wavelength-dependent King correction, pressure scaling,
 * and WGS84 gravity.  Historical reference-code comparison models were
 * removed from the production core; callers retain this one-value enum only
 * so the public option structures remain source compatible during the v1.2
 * transition.
 */
typedef enum {
    RT_RAYLEIGH_MODEL_BODHAINE_1999 = 0
} rt_rayleigh_model_t;

/* Source of per-case Rayleigh optical depth.
 *
 *   AUTO (default, production):
 *       tau_R = rt_rayleigh_tau_model(wavelength_nm, opts.rayleigh_model)
 *       i.e. computed internally from the selected rayleigh_model.
 *
 *   FROM_INPUT (validation / algorithm-pure):
 *       tau_R = cs->tau_R_input.  Used for controlled comparisons in
 *       which the RT solvers must receive an identical column optical depth.
 *
 * The opts.rayleigh_model field is still recorded in atm metadata for
 * diagnostics, but the actual tau_R used is governed by opts.tau_r_source.
 */
typedef enum {
    RT_TAU_R_AUTO       = 0,
    RT_TAU_R_FROM_INPUT = 1
} rt_tau_r_source_t;

/* Atmosphere description.
 *
 * Direction array convention:
 *   rm[0]   = -mu_sun           (solar slot, negative)
 *   rm[+j]  = +mu_quad[j-1]     for j = 1..n_mu  (upward)
 *   rm[-j]  = -mu_quad[j-1]     for j = 1..n_mu  (downward, mirror)
 *   gb[0]   = 0                 (unused solar slot)
 *   gb[±j]  =  w_quad[|j|-1]    Gauss weight, sum_{j=1..n_mu} = 1
 *   xpl[k]  = plm[2][k] for the CURRENT m (renormalized Ptil_2^m;
 *             classical P_2 only at m = 0; even parity holds per mode)
 *
 * The rm/gb/xpl members are offset pointers into rm_storage/etc., so that
 * rm[-n_mu .. +n_mu] indexing is legal. Storage pointers are owned and
 * freed by rt_atm_free(); never free the offset pointer.
 */
typedef struct {
    int    n_layers;        /* nt: number of layers (nt+1 levels k=0..nt) */
    int    n_mu;            /* mu: number of positive Gauss nodes */
    double tau_total;       /* total Rayleigh optical depth (column) */
    double depol;            /* raw depolarization factor (~0.0279) */
    double ssa;              /* single-scattering albedo (=1 for Rayleigh) */
    double mu_sun;          /* cos(SZA), positive */
    rt_rayleigh_model_t rayleigh_model;  /* model that produced tau_total
                                          * (metadata for diagnostics; actual
                                          * tau dispatch happens in caller). */

    /* Rayleigh phase-function Legendre coefficients (m=0):
     *   aaaa = δ / (2 - δ)                     (depol convention transform)
     *   ron  = (1 - aaaa) / (1 + 2*aaaa)
     *   β_0  = 1
     *   β_2  = 0.5 × ron                       (Rayleigh I)
     *   γ_2  = -ron × √(1.5)                   (Rayleigh Q, vector RT)
     *   α_2  = 3 × ron                          (Rayleigh U, vector RT)
     */
    double beta0;
    double beta2;
    double gamma2;          /* vector — Rayleigh Q phase coefficient */
    double alpha2;          /* vector — Rayleigh U phase coefficient */

    /* Layer arrays [n_layers + 1] (level-wise, k = 0..nt). */
    double *h;              /* cumulative tau at level k from TOA (incl. gas absorption if applied) */
    double *ch;             /* exp(-h[k]/mu_sun) / 2  (note /2 factor) */
    double *xdel;           /* aerosol fraction (=0 in Phase 1) */
    double *ydel;           /* Rayleigh fraction (=1 in Phase 1) */

    /* v1.02 atmospheric absorption (SOS-integrated) extensions.
     * Populated by rt_atm_apply_gas_absorption() when --use-absorption is active.
     * NULL when absorption is disabled (bit-exact baseline regression). */
    double *z_km_level;     /* altitude (km) at each level k = 0..nt; nt = TOA, 0 = ground */
    double *tau_abs_layer;  /* per-layer τ_abs (incremental, between levels k and k+1); size nt */
    double  tau_abs_total;  /* sum(tau_abs_layer) */

    /* v1.11 (2026-09-05): the legacy average-secant Chapman PSSA fields
     * (the average-secant Chapman slant-path caches) were REMOVED together
     * with the legacy source files.  IPSS never modifies rt_atm_t: it
     * rescales the finished plane-parallel Stokes vector (see rt_ipss.h). */

    /* Direction arrays — owned storage + offset pointers. */
    double *rm_storage;     /* size 2*n_mu+1 */
    double *gb_storage;
    double *xpl_storage;
    double *xrl_storage;    /* vector — R_2(rm[k]) storage */
    double *xtl_storage;    /* vector — T_2(rm[k]) storage */
    double *rm;             /* &rm_storage[n_mu]  → indexable [-n_mu..+n_mu] */
    double *gb;
    double *xpl;
    double *xrl;            /* vector — R_2 at quadrature directions */
    double *xtl;            /* vector — T_2 at quadrature directions */

    /* ─── Phase 3: Aerosol expansion coefficients (m=0 phase) ─────────
     *
     * When aerosol_active = 0: the four arrays are NULL and the solver
     * passes NULL into rt_kernel_phase_fourier{,_pol}() (pure Rayleigh
     * path -- bit-level identical to Phase 2).
     *
     * When aerosol_active = 1: arrays of length L_max+1 hold the
     * aerosol Legendre coefficients computed by
     * rt_aerosol_compute_vector_legendre().
     *
     *   betal_aer[l]  = scalar β_l of P11
     *   gammal_aer[l] = vector γ_l of P12 (Mishchenko generalized P²_l)
     *   alphal_aer[l] = vector α_l (recursive from β, δ)
     *   zetal_aer[l]  = vector ζ_l (recursive from β, δ)
     *
     * Note: For spherical Mie aerosols (.mie format used by V3), the
     * P34 element vanishes so δ_l = 0; α_l/ζ_l simplify but the same
     * generalized spherical-function completion is used for all orders.
     */
    int    aerosol_active;
    int    L_max;             /* max Legendre order (truncated) */
    double *betal_aer;         /* size L_max+1, owned */
    double *gammal_aer;        /* size L_max+1, owned */
    double *alphal_aer;        /* size L_max+1, owned */
    double *zetal_aer;         /* size L_max+1, owned */

    /* Angle-space value-kernel path (Gibbs-free). When aer_use_value_kernel != 0
     * the SOS builds ws->phase_fourier_m by direct azimuth integration of these
     * (loglin-capped) phase VALUES rather than from the Legendre moments above,
     * which cannot represent the backscatter glory. These pointers are BORROWED
     * from the rt_aerosol_input_t (NOT owned/freed by rt_atm_free). */
    int    aer_use_value_kernel;
    int    aer_n_ang_phase;
    const double *aer_theta_phase;   /* ascending degrees, size aer_n_ang_phase */
    const double *aer_P11_phase;
    const double *aer_P12_phase;
    const double *aer_P33_phase;

    /* Immutable formal-solution transport cache.  Once h[] and rm[] are
     * finalized, transmission(direction, layer) is invariant across Stokes
     * components and all scattering orders.  The cache stores exact exp()
     * results plus byte snapshots of h/rm; rt_solver_integrate_bcs rebuilds it
     * only when the medium geometry changes.  Owned/freed by rt_atm_free(). */
    double *transport_T;          /* [signed direction 0..2*n_mu][layer] */
    double *transport_h_key;      /* snapshot [n_layers+1] */
    double *transport_rm_key;     /* snapshot [2*n_mu+1], signed ring */
    int     transport_cache_ready;
    unsigned long long transport_cache_builds;
    unsigned long long transport_cache_hits;

    /* v1.059: incident-beam polarization. Q/I of the solar beam after air-water
     * Fresnel transmission into the medium (= M_T_aw[1,0]/M_T_aw[0,0]). The
     * pre-v1.059 primary source treated the in-water beam as unpolarized
     * [I,0,0]; with beam_q != 0 the polarized primary source adds the 2nd-column
     * (1,2)/(2,2)/(3,2) phase coupling of the beam Q. Set only by the vector
     * water path (rt_water_rt_sos_pure); 0 in the scalar path and atmosphere,
     * so scalar regressions are bit-invariant. */
    double beam_q;
} rt_atm_t;

/* Radiation field per Fourier order m (Phase 1: scalar I only). */
typedef struct {
    int     n_levels;       /* n_layers+1 */
    int     n_mu;           /* same as atm */
    double *i1;             /* [(n_levels) * (2*n_mu+1)] field per direction */
    double *i2;             /* same shape: source function */
    int     n_orders;       /* SOS orders consumed */
} rt_field_t;

/* Per-case input geometry/spectrum. */
typedef struct {
    double sza_deg;         /* solar zenith angle */
    double vza_deg;         /* view zenith */
    double raa_deg;         /* relative azimuth */
    double wavelength_nm;
    rt_surface_kind_t surface;
    int    rayleigh_on;     /* DERIVED (not a user toggle): set = (pressure_hpa > 0) after arg parse. */
    int    aerosol_on;      /* DERIVED (not a user toggle): set = (aerosol AOD > 0) after arg parse. */
    double pressure_hpa;    /* surface pressure (hPa); scales Rayleigh tau (tau_R ∝ P/1013.25).
                             * 0 => no Rayleigh.  Default 1013.25 (standard atmosphere). */
    /* Pre-supplied Rayleigh optical depth for opts.tau_r_source =
     * RT_TAU_R_FROM_INPUT.  Used directly as tau_total in
     * rt_atm_build_rayleigh.  Set to a positive value to enable;
     * any value <= 0 is treated as "unset" and the AUTO path is
     * required.  Default 0 (caller must explicitly set). */
    double tau_R_input;

    /* ====================================================================
     * Phase 4 (rough Fresnel BRDF over a black ocean) — used when surface != RT_SURFACE_BLACK
     * ==================================================================== */
    /* Refractive index of water (real).  Default 1.34 -- 2026-07-15 Jae(G4)
     * 확정(Step57 의 Quan-Fry(λ) 자동 기본은 폐지; OSOAA/AF1982 1.34 정합).
     * Used for both flat-Fresnel (RT_SURFACE_FLAT, ws=0) and rough Fresnel surfaces. */
    double n_water;
    /* Wind speed (m/s) for the rough Fresnel surface slope distribution.
     * 0.0 → flat Fresnel (regardless of surface kind).
     * >0  → wave-slope-weighted rough Fresnel BRDF. */
    double wind_speed;
    /* Slope variance model:
     *   0 = Nakajima-Tanaka (sigma = 0.0731 sqrt(ws), default)
     *   1 = OCRT floor law (sigma^2 = 0.003 + 0.00512 max(0.01,ws); CM54-derived) */
    int    sigma_type;
    /* Q kernel sign convention (see V1 sos_solver_t::q_convention notes).
     *   0 = legacy (Rs-Rp)/2 — matches PSTAR/SeaDAS reference at <2%
     *   1 = Mishchenko (Rp-Rs)/2 — matches PSTAR principal-plane bit-level */
    int    q_convention;
    /* 1 → exclude direct sunglint single-bounce term from target output
     *     (post-hoc subtraction; AF1982 reference is decoupled) */
    int    decouple_sunglint;

    /* ====================================================================
     * Phase B (in-water RT — RT_SURFACE_OCEAN) — used when surface == OCEAN
     * Defaults are safe and unused for non-ocean surface kinds.
     * ==================================================================== */
    double T_water_C;       /* water temperature [°C], default 20.0 (Röttgers 2014 ψ_T LUT range -2…40) */
    double S_water_gkg;     /* salinity [g/kg]. default 38.4 (B.3 ignores; psi_S deferred to B.5+) */
    double F_sun;           /* incident solar irradiance just-above-surface [W/m²/nm or unit-normalized].
                             * default π → OCRT internal F_sun=π convention. */

    /* CDOM 2-parameter absorption model (B.5, 2026-05-23):
     *   a_CDOM(λ) = a_CDOM(λ_ref) · exp(-S_CDOM · (λ - λ_ref))
     * CDOM 은 흡광만 가지며 산란/위상함수에 기여하지 않는다 (Mueller matrix
     * 변경 없음). 단 ω₀ = b_w/(a_w+a_CDOM+b_w) 가 감소하여 다중산란 결과
     * (특히 polarization Q,U) 가 변한다 — 이는 RT 솔버가 자동 처리한다.
     *
     * Defaults: a_CDOM(440)=0 (즉 CDOM off), S=0.014 nm⁻¹ (Bricaud 1981
     * classical value, open-ocean reference), λ_ref=440 nm. */
    double a_cdom_440_m_inv;    /* CDOM absorption at reference λ [m⁻¹], default 0 (off) */
    double S_cdom_nm_inv;       /* CDOM spectral slope [nm⁻¹], default 0.014 */
    double cdom_ref_lambda_nm;  /* reference wavelength [nm], default 440 */

    /* Public water-input branch.  Exactly one of OCRT/CCRR/IOP must be
     * selected for RT_SURFACE_OCEAN; see rt_water_input_mode_t.  The legacy
     * ccrr_mode/fixed_bulk_iop_mode fields below are derived execution flags
     * retained for the existing water solver internals. */
    int    water_input_mode;

    /* Constituent-to-IOP adapter.  This does not port the CCRR RT solver; it
     * only maps CHL/TSM/aDOM into OCRT water-RT inputs. */
    int    ccrr_mode;
    /* EAP phytoplankton scattering gate.  The 17-species phase generator is
     * integrated, but constituent-model EAP scattering is disabled until a
     * validated forward-peak treatment is available for the L=200 water
     * moment representation.  0 = absorption only; 1 = internal diagnostic
     * re-enable (not exposed by the production CLI). */
    int    organic_phyto_scattering;
    int    water_constituent_model; /* rt_water_constituent_model_t; OCRT default */
    double ccrr_chl_mg_m3;
    double ccrr_min_g_m3;
    int    tsm_species;       /* ahn_species_t; 0=red_clay default */
    int    organic_phyto_group;      /* organic_phyto_group_t; micro default */
    double detritus_a440_m_inv;      /* optional organic-detritus a_d(440) */
    double detritus_slope_nm_inv;    /* Bricaud-Stramski S_d; default 0.0109 */
    const char *ccrr_phase_moments_path; /* optional CSV: phase_name,l,chi_l[/betal_l]; particle moments only */
    const char *ccrr_particle_phase_lut_path; /* positive high-res particle P11 LUT */
    const char *ccrr_particle_phase_case_id;  /* optional row selector */
    double ccrr_particle_phase_wavelength_nm; /* optional row selector */
    int    ccrr_particle_phase_lmax;          /* formal delta-M retained order */
    int    ccrr_particle_phase_nphi;          /* azimuth quadrature */
    int    water_phase_kernel;                /* 0=value (path B, default), 1=moment (path A shell) */

    /* HydroLight fixed-total-IOP bulk-phase validation mode.  This bypasses
     * constituent conversion and injects a_total,b_total,bb_total directly. */
    int    fixed_bulk_iop_mode;
    double fixed_a_total_m_inv;
    double fixed_b_total_m_inv;
    double fixed_bb_total_m_inv;
    int    fixed_bulk_lmax;
    int    fixed_bulk_phase_model;
    int    fixed_bulk_phase_nphi;
    double fixed_bulk_ff_n;
    double fixed_bulk_ff_mu;
    const char *fixed_bulk_phase_lut_path;
    const char *water_mie_phase_path; /* #0: optional .mie vector phase (P11/P12/P33) for bulk-water particle */
    int water_mie_moment_mode;       /* 0=dense theta/PCHIP, 1=Gauss/μ-cubic */
    int water_mie_truncation_mode;   /* branch-prefixed --ocrt/--iop-mie-truncation; 0=off, 1=delta truncation */
    int water_mie_ss_mode;           /* branch-prefixed --ocrt/--iop-mie-ss-mode: 0=none, 1=IMS, 2=NT-TMS */
    int water_mie_moment_n_mu;       /* default 400 */
    const char *fixed_bulk_phase_case_id;
    double fixed_bulk_phase_wavelength_nm;
    /* Scalar Fournier-Forand fixed-total-IOP validation mode.
     * Reads a_total,b_total from an external ASCII/CSV table and uses
     * bb_total=b_total*scalar_ff_bb_over_b with analytic positive FF P11. */
    const char *scalar_ff_iop_path;
    double scalar_ff_bb_over_b;
    double scalar_ff_refractive_index;
    double scalar_ff_mu;

    /* Debug-only in-water vertical domain overrides. Production uses automatic
     * deep-water defaults so ordinary users cannot under-resolve the slab.
     * Populated only by --debug-water-* options gated by OCRT_DEBUG=1. */
    int    n_layers_water_override;
    double tau_max_target_override;       /* debug base optical depth */
    double water_max_tau_override;        /* debug maximum optical depth cap */
    double water_max_z_max_m_override;    /* debug physical depth ceiling [m] */
    double water_depth_bottom_tol_override; /* debug adaptive bottom attenuation tolerance */
    double water_layer_dtau_target_override;
    int    n_mu_water_override;     /* 0 → use default 24; CLI --n-mu-water */
    int    n_mu_water_sky;         /* 0 => same as n_mu_water; CLI --n-mu-water-sky.
                                    * Coarse in-water SOS grid for the atmosphere sky-light
                                    * coupling beams ONLY.  Skylight is a small (~2%) additive
                                    * correction to rrs(0-), so a coarse grid there is acceptable
                                    * while the direct solar beam keeps the full (accurate)
                                    * n_mu_water.  Cuts the per-sky-beam SOS cost ~ (full/sky)^2. */
    int    water_m_max_override;    /* 0 → auto; CLI --water-m-max */

    /* Debug-only in-water RT SOS max-iteration override. Production uses
     * automatic caps and exits early on convergence. */
    int    water_max_orders_override;
    /* Debug-only: auto minimum SOS cap for particle-rich water media.
     * Production default is 100 and the solver may still raise it higher
     * for fixed-bulk / CCRR / OSOAA-PM paths. This value never lowers a
     * larger auto-selected cap. */
    int    debug_high_particle_water_max_orders;

    /* In-water view extraction policy.
     * 0 (default): continuous reconstruction/interpolation from the quadrature
     *    directions. This matches CCRR v0.5.0 direct-only reference bookkeeping.
     * 1: add requested view direction as a zero-weight output node when needed
     *    (near/exact nadir diagnostic; can reduce azimuthal artefacts but changes
     *    the discrete ordinates solution relative to continuous reconstruction). */
    /* v1.09 #20 (S-007 root fix): vza list (deg, air side) appended as
     * zero-weight in-water view nodes for direct grid extraction. */
    const double *water_view_vza_list;
    int           n_water_view_vza;
    /* v1.09 #22 (B-0: OSOAA-parity shared angular grid): when set, the water
     * quadrature uses the SAME node count as the atmosphere (identical
     * generator -> bit-identical mu set), dissolving n_mu_water as an
     * independent parameter.  OSOAA uses one grid for air and sea. */
    int           water_shared_grid;
    int    water_view_as_node;
} rt_case_t;

/* Per-case output reflectances at TOA.
 *
 * Phase 1: scalar (rho_I only); rho_Q/rho_U are zero placeholders.
 * Phase 2 will populate Q/U.
 *
 * Diagnostics (I_TOA, tau_R_total, sos_*_per_m) are filled in Step 9
 * by rt_solve_case().  External callers never need them and may
 * ignore the diagnostic block; it's kept for paper figures, batch
 * runs, and debugging.  sos_orders_per_m / sos_residual_per_m are
 * sized for fourier_m_max <= 2 (Phase 1); Phase 3 will replace the
 * fixed [3] arrays with dynamic allocation when m_max grows.
 */
typedef struct {
    double rho_I;
    double rho_Q;           /* zero in Phase 1 */
    double rho_U;           /* zero in Phase 1 */
    int    n_orders_used;   /* max over m */
    int    converged;       /* AND over m */

    /* Diagnostics (Step 9) */
    double I_TOA;                 /* upwelling intensity, pre rho conversion */
    double Q_TOA;                 /* Phase 2 vector — 0 if rt_solve_case
                                   * (scalar) was used; populated by
                                   * rt_solve_case_pol */
    double U_TOA;                 /* Phase 2 vector — same caveat as Q_TOA */
    double tau_R_total;           /* total Rayleigh column from Bodhaine model */
    int    sos_orders_per_m[33];  /* m = 0..32 (Phase 3 aerosol) */
    double sos_residual_per_m[33]; /* m = 0..32 (Phase 3 aerosol) */

    /* v1.01 Atmospheric transmittance decomposition.
     * F_TOA = π convention used internally (RT_F_SOLAR_PI).
     *   T_dir_dn       = exp(-τ_tot / μ_sun)          (Beer-Lambert direct beam)
     *   T_diff_dn_hemi = (2π Σ_quad μ·w·I_dn^{m=0}(μ)) / (π·μ_sun)
     *                                                 (hemispheric diffuse down)
     *   T_diff_dn_dir  = 2π·μ_v·I_dn(μ_v, φ) / (π·μ_sun)   (directional, at vza/raa)
     *   T_sg_up_dir    = ρ_glint(μ_v, φ) · μ_sun      (legacy flux-ratio diagnostic)
     *   T_total_up_dir = ρ_TOA(μ_v, φ) · μ_sun        (legacy flux-ratio diagnostic)
     *
     * These last two are NOT atmospheric transmission coefficients.  Public
     * output labels identify them as diagnostics. */
    double T_dir_dn;
    double T_diff_dn_hemi;
    double T_diff_dn_dir;
    double T_sg_up_dir;
    double T_total_up_dir;

    /* Explicit radiometric decomposition for scientific output.
     * All rho_* fields follow OCRT Convention B and are dimensionless BRF-like
     * directional reflectances at TOA.  For ocean runs the components sum to
     * rho_{I,Q,U} (glint-decoupled), while rho_*_glint_direct is the additional
     * direct solar sunglint term. */
    double rho_atm_path_I, rho_atm_path_Q, rho_atm_path_U;
    double rho_water_direct_I, rho_water_direct_Q, rho_water_direct_U;
    double rho_water_sky_I, rho_water_sky_Q, rho_water_sky_U;
    double rho_water_total_I, rho_water_total_Q, rho_water_total_U;
    double rho_glint_direct_I, rho_glint_direct_Q, rho_glint_direct_U;

    /* Atmospheric transmission diagnostics.
     *
     * Downward quantities are obtained from the solved BOA atmosphere field:
     *   T_dir_dn         exact unscattered beam attenuation;
     *   T_diff_dn_hemi   hemispheric diffuse BOA irradiance / incident irradiance;
     *   T_total_dn_hemi  their exact sum.
     *
     * Upward quantities are defined for the ACTUAL ocean water-leaving field
     * at the requested view direction, so atmospheric angular redistribution
     * and the water BRDF are retained:
     *
     *   T_total_up_view = (I_TOA_total - I_atm_path) / Lu(0+)
     *   T_diff_up_view  = T_total_up_view - T_dir_up_view
     *
     * T_dir_up_view is the exact unscattered Beer-Lambert component.
     * T_up_rt_valid is one only when the upward water signal was propagated by
     * the rigorous atmospheric bottom-source SOS pass (or no atmosphere is
     * present).  Invalid/undefined values are NaN; no reciprocity or
     * hemispheric approximation is exported. */
    double T_total_dn_hemi;
    double T_dir_up_view;
    double T_diff_up_view;
    double T_total_up_view;
    double I_TOA_water_signal; /* exact numerator: I_TOA_total - I_atm_path */
    int    T_up_rt_valid;

    /* Aerosol spectral-input metadata. AOD is always supplied at a fixed
     * reference wavelength and converted internally using .mie extinction. */
    double aerosol_aod_ref;
    double aerosol_aod_ref_nm;
    double aerosol_aod_band;
    double aerosol_extinction_ratio;

    /* ====================================================================
     * Phase B (in-water RT — RT_SURFACE_OCEAN) — populated when ocean solver
     * is used. Zero for non-ocean surface kinds.
     * ====================================================================
     * Notation (per user convention):
     *   r_rs = lowercase r → underwater (0⁻)  remote-sensing reflectance
     *   R_rs = uppercase R → above-water (0⁺) remote-sensing reflectance
     */
    double Lu_0minus_view;  /* in-water upwelling radiance at view direction [W·m⁻²·sr⁻¹·nm⁻¹] */
    double Lu_0plus_view;   /* above-water (air-side) upwelling radiance at view direction */
    double Qu_0minus_view;  /* in-water Stokes Q at view direction */
    double Qu_0plus_view;   /* above-water Stokes Q (after T_wa Mueller) */
    double Uu_0minus_view;  /* in-water Stokes U */
    double Uu_0plus_view;   /* above-water Stokes U */
    double Ed_0minus_water; /* just-below-surface downward irradiance [W·m⁻²·nm⁻¹] */
    double Ed_0plus_air;    /* just-above-surface downward irradiance */
    double Eu_0minus_water; /* just-below-surface upward irradiance */

    /* Direct/diffuse irradiance decomposition at the sea surface.
     *
     * Definition used by these fields is the standard angular radiometric
     * decomposition, NOT a photon-source provenance decomposition:
     *   direct  = the unscattered collimated solar beam at that interface,
     *   diffuse = the remaining hemispheric radiance integral.
     *
     * Therefore, with the current deep-water / black-bottom model there is no
     * upward collimated beam at 0-, and Eu_0minus_direct is identically zero;
     * all upwelling irradiance is diffuse.  The totals above are preserved and
     * the following identities hold (to round-off):
     *   Ed_0plus_air    = Ed_0plus_direct  + Ed_0plus_diffuse
     *   Ed_0minus_water = Ed_0minus_direct + Ed_0minus_diffuse
     *   Eu_0minus_water = Eu_0minus_direct + Eu_0minus_diffuse
     *
     * These six diagnostics are assembled from quantities already computed by
     * the coupled solve; no additional SOS pass or angular integration is
     * performed. */
    double Ed_0plus_direct;
    double Ed_0plus_diffuse;
    double Ed_0minus_direct;
    double Ed_0minus_diffuse;
    double Eu_0minus_direct;
    double Eu_0minus_diffuse;

    double r_rs_0minus;     /* underwater I/Ed(0⁻) [sr⁻¹] */
    double r_rs_0minus_Q;   /* underwater Q/Ed(0⁻) [sr⁻¹] */
    double r_rs_0minus_U;   /* underwater U/Ed(0⁻) [sr⁻¹] */
    double R_rs_0plus;      /* above-water I/Ed(0⁺) [sr⁻¹] */
    double R_rs_0plus_Q;    /* above-water Q/Ed(0⁺) [sr⁻¹] */
    double R_rs_0plus_U;    /* above-water U/Ed(0⁺) [sr⁻¹] */
    double Kd_0minus;       /* downward diffuse attenuation at z=0⁻ [m⁻¹] */
    double Ku_0minus;       /* upward diffuse attenuation at z=0⁻ [m⁻¹] */

    /* IOP metadata (record IOPs actually used — useful for diagnostics & post-hoc QAA) */
    double a_w_used;        /* pure water absorption [m⁻¹] */
    double b_w_used;        /* pure water scattering [m⁻¹] */
    double bb_w_used;       /* pure water back-scattering [m⁻¹] */
    double a_cdom_used;
    double a_pig_used;
    double a_chl_used;       /* phytoplankton absorption only [m⁻¹] */
    double b_pig_used;
    double bb_pig_used;
    double a_min_used;
    double b_min_used;
    double bb_min_used;
    double a_total_used;
    double b_total_used;
    double bb_total_used;
    double omega_water;     /* single-scattering albedo b_total/(a_total+b_total) */
    double tau_max_water;   /* water-column optical depth used (infinite-deep proxy) */

    /* Dual sunglint output: rho_I/Q/U above are glint-DECOUPLED (AF1982 ref);
     * these carry the SAME TOA Stokes WITH the direct solar sunglint added.
     * A single ocean run emits both glint-off (rho_I) and glint-on
     * (rho_I_glint) — no option toggles between them. */
    double rho_I_glint;
    double rho_Q_glint;
    double rho_U_glint;
} rt_result_t;

/* Solver options. */
typedef struct {
    int    n_mu;            /* Gauss-Legendre positive nodes (default 24) */
    int    n_layers;        /* atmosphere subdivision (default 40) */
    int    max_orders;      /* SOS truncation (default 20) */
    int    l_max;           /* max Legendre order; -1 = auto (Phase 1 → 2,
                             * Phase 3 → 2*n_mu - 2) */
    rt_rayleigh_model_t rayleigh_model;  /* default = BODHAINE_1999 */
    rt_integration_method_t integration_method;  /* default = LINEAR */
    rt_solver_sos_options_t sos;        /* SOS iteration control */
    int    fourier_m_max;   /* max Fourier order m (default 2 for Phase 1
                             * Rayleigh; Phase 3 aerosol may need 8-16) */
    rt_tau_r_source_t tau_r_source;     /* AUTO (default) or FROM_INPUT */
    /* 1 = spherical-shell correction (IPSS, Zhai & Hu 2022); default 0.
     * v1.11 (2026-09-05): plain on/off.  There is no algorithm selector any
     * more — the legacy average-secant Chapman PSSA was deleted, so the SOS
     * solve is ALWAYS plane-parallel and IPSS only rescales its output. */
    int    pssa;
    int    view_as_node;    /* 1 (Step162 production default) = append mu_view as
                             * an extra zero-weight output direction, preserving
                             * all original GL nodes/weights and reading I/Q/U
                             * directly at that node.  The old linear
                             * interpolation route is diagnostic only and is not
                             * used by the physical-strict water route. */
    int    vector_mode;     /* 1 (default, 2026-07-14 Jae 지시) = full vector
                             * I/Q/U path via rt_solve_case_pol.  The scalar
                             * I-only solver rt_solve_case was DELETED the same
                             * date (no future validation planned); this flag is
                             * retained only because LUT/batch output shaping
                             * still reads it.  CLI --vector is a no-op kept for
                             * script compatibility. */

    /* Aerosol runtime configuration shared by single-case and batch modes.
     * These fields are populated by main.c from --mie / --aod-555|--aod-865 / --aer-l-max
     * / truncation CLI options.  Batch mode uses aerosol_aod only when the
     * input CSV has no per-row aod_555/aod_865 column. */
    const char *aerosol_mie_path;
    double aerosol_aod;          /* AOD at aerosol_aod_ref_nm */
    double aerosol_aod_ref_nm;   /* fixed reference wavelength, default 555 nm */
    int    aerosol_l_max;
    double aerosol_theta_cut_deg;
    int    aerosol_delta_m_N;
    int    aerosol_apply_nt_tau;

    /* Output verbosity for batch CSV writer (v1.01).
     * 0 = SIMPLE (default): case geometry + rho_I/Q/U + tau_R only
     * 1 = DEBUG: all v0.9 columns (ref, delta, n_orders, walltime)
     * Note: single-case stdout is unchanged by this flag. */
    int    output_mode_debug;

    /* v1.01 LUT batch options (used when rt_io_run_batch is invoked with
     * lut_enable=1; ignored otherwise). */
    int          lut_enable;        /* 0 = no LUT output, 1 = per-case LUT */
    const char  *lut_output_dir;    /* output directory for case_<id>.csv */
    double       lut_vza_step;      /* degrees (default 2.5) */
    double       lut_raa_step;      /* degrees (default 5.0) */
    double       lut_vza_max;       /* degrees (default 85.0) */

    /* v1.02: SOS-integrated gas absorption.
     * If non-NULL, rt_atm_apply_gas_absorption is called after build,
     * which modifies h[k], ch[k], xdel[k], ydel[k] per AFGL profile + xsec.
     * NULL = no absorption (bit-exact baseline regression).
     *
     * Type-erased to keep rt_types.h independent of rt_absorption.h
     * (cyclic dependency avoidance). Caller passes a rt_absorption_t*.
     */
    const void  *abs_state;
    double aer_h_km;        /* OCRT aerosol profile scale height [km].
                             * Production model: exp(-z/2 km).
                             * Positive advanced overrides are allowed; zero
                             * and negative values are invalid. */

    double conv_tol;        /* relative tolerance (legacy; see opts.sos.tolerance) */
    int    verbose;

    /* C3 (water->air up-coupling) pass-2 driver: per-mode water-leaving 0+
     * Fourier field at the atm mu-nodes, indexed [m*ext_bottom_n_mu + k].
     *
     * The three pointers form one shaped object.  ext_bottom_n_mu is the
     * positive-direction stride and ext_bottom_m_max is the highest stored
     * Fourier mode (inclusive).  These dimensions are mandatory whenever the
     * pointers are non-NULL; the atmospheric solver validates them before it
     * creates a mode slice.  Modes above ext_bottom_m_max are a physically zero
     * external source and are never dereferenced.
     *
     * When bottom_source_only is set, rt_solve_case_pol_impl zeros the solar
     * primary and feeds this object as the atm SOS bottom source, yielding the
     * TOA water-leaving contribution.  Set by the ocean driver for pass 2 ONLY;
     * NULL/0/-1 for every normal (pass-1) solve. */
    const double *ext_bottom_per_m_I;  /* may be NULL */
    const double *ext_bottom_per_m_Q;  /* may be NULL */
    const double *ext_bottom_per_m_U;  /* may be NULL */
    int           ext_bottom_n_mu;     /* source stride; 0 when unused */
    int           ext_bottom_m_max;    /* inclusive mode bound; -1 when unused */
    int           bottom_source_only;  /* 0 = normal solar solve */
    /* v1.10 S17 (2026-07-12, bit-safe speed): OPTIONAL per-m view-sample
     * export from the vector solve - the caller can redo the phi
     * reconstruction itself (same function, same inputs => bit-identical)
     * to fold the raa axis of grid batches.  NULL (default) = inert. */
    double *toa_view_per_m_I; double *toa_view_per_m_Q; double *toa_view_per_m_U;

} rt_options_t;

/* Defaults — keep in sync with main.c argument parsing.
 * Phase 1: n_layers bumped to 40 (was 20 in scaffold) for safer
 * accuracy at large τ_R; cheap.
 */
static inline rt_options_t rt_options_default(void) {
    rt_options_t o = { .n_mu = 24, .n_layers = 40, .max_orders = 20,
                       .l_max = -1,
                       .rayleigh_model = RT_RAYLEIGH_MODEL_BODHAINE_1999,
                       .integration_method = RT_INTEGRATION_METHOD_LINEAR,
                       .sos = { .max_iterations = 20,
                                .tolerance      = 1.0e-7,
                                .acceleration   = RT_SOS_ACCELERATION_PLAIN,
                                .save_orders    = 0 },
                       .fourier_m_max = 2,
                       .tau_r_source = RT_TAU_R_AUTO,
                       .pssa = 0,
                       .view_as_node = 1,
                       .vector_mode  = 1,
                       .aerosol_mie_path = NULL,
                       .aerosol_aod = 0.0,
                       .aerosol_aod_ref_nm = 555.0,
                       .aerosol_l_max = 80,
                       .aerosol_theta_cut_deg = 0.0,
                       .aerosol_delta_m_N = 0,
                       .aerosol_apply_nt_tau = 0,
                       .output_mode_debug = 0,
                       .lut_enable = 0,
                       .lut_output_dir = NULL,
                       .lut_vza_step = 2.5,
                       .lut_raa_step = 5.0,
                       .lut_vza_max = 85.0,
                       .abs_state = NULL,
                       .aer_h_km = 2.0,
                       .ext_bottom_per_m_I = NULL,
                       .ext_bottom_per_m_Q = NULL,
                       .ext_bottom_per_m_U = NULL,
                       .ext_bottom_n_mu = 0,
                       .ext_bottom_m_max = -1,
                       .bottom_source_only = 0,
                       .conv_tol = 1.0e-6, .verbose = 0 };
    return o;
}

/* ─── LUT grid output (v1.01).  When provided to rt_solve_case_pol_lut(),
 *     the solver runs SOS ONCE and reconstructs the radiation field at
 *     every (vza, raa) grid point — no extra SOS iterations.
 *
 *     Conventions:
 *       - vza_deg[i] in [0, 89]  (90° is horizon; numerically unstable).
 *       - raa_deg[i] in [0, 360) interpreted in the same convention as
 *         rt_case_t.raa_deg (OSOAA Phi internal: 0 = forward, 180 = back).
 *       - rho_I / rho_Q / rho_U have row-major layout rho[iv * n_raa + ir].
 *       - All arrays are caller-allocated.
 *       - Sun-glint direct contribution (when decouple_sunglint=0 and
 *         surface=BLACK_FRESNEL_OCEAN) is added to rho_I/Q/U at each grid point. */
typedef struct {
    int n_vza;
    const double *vza_deg;     /* size n_vza */
    int n_raa;
    const double *raa_deg;     /* size n_raa */
    double *rho_I;             /* size n_vza * n_raa */
    double *rho_Q;             /* size n_vza * n_raa */
    double *rho_U;             /* size n_vza * n_raa */

    /* v1.10 S7 (atm grid cache): optional RAW TOA intensity grids — the
     * pure atmospheric-path per-direction reconstruction BEFORE the rho
     * division and BEFORE any sunglint addition.  Bit-exact source for
     * replaying atm_res.I/Q/U_TOA per row.  NULL to skip (legacy callers). */
    double *raw_I; double *raw_Q; double *raw_U;

    /* Optional per-view Fourier samples before the RAA reconstruction.
     * Layout: view_per_m_X[m * n_vza + iv], m=0..view_m_max.
     * These arrays let the coupled ocean LUT path reconstruct the complete
     * RAA axis without re-entering any atmospheric solver.  NULL keeps the
     * historical behavior. */
    double *view_per_m_I;
    double *view_per_m_Q;
    double *view_per_m_U;
    int     view_m_max;          /* caller capacity, inclusive; -1 if unused */
    int     view_m_max_filled;   /* solver output, inclusive */

    /* v1.01 transmittances at each grid point (caller-allocated; may be NULL
     * to skip). Layout matches rho_*: row-major [iv * n_raa + ir]. */
    double *T_diff_dn_dir;     /* downward diffuse transmittance per direction */
    double *T_sg_up_dir;       /* legacy diagnostic: rho_glint * mu_sun, not atm T */
    double *T_total_up_dir;    /* legacy diagnostic: rho_TOA * mu_sun, not atm T */

    /* Case-level (scalar) transmittances filled by solver (not per direction). */
    double T_dir_dn;           /* exp(-τ/μ_sun), independent of view */
    double T_diff_dn_hemi;     /* hemispheric flux ratio, independent of view */
} rt_lut_grid_out_t;

/* ====================================================================
 * Atmospheric BOA downward Stokes field export (Phase B.4 Stage 2b)
 *
 * Captures the full atmospheric SOS solution at z = BOA (level k = nt)
 * in the downward direction, per Fourier m mode, per quadrature node.
 * Used by the ocean coupling path to construct the in-water boundary
 * condition: each (μ_j, m) atmospheric BOA radiance is Snell-refracted
 * into water and T_aw-Mueller-transmitted to form the in-water inbound
 * field at z = 0⁻ for the in-water SOS source.
 *
 * Caller responsibilities:
 *   1. Allocate the struct fields (or call rt_atm_boa_export_alloc).
 *      I_per_m, Q_per_m, U_per_m each [(m_max+1) × n_mu] row-major:
 *      field[m * n_mu + (j-1)] = X^m(μ_j) where μ_j > 0 is the j-th
 *      positive Gauss-Legendre node (atm convention: positive μ in this
 *      export struct represents the *downward direction* with cos(zenith)=μ_j).
 *   2. mu_quad_pos [n_mu]: copy of the positive-μ quadrature grid.
 *   3. After call, free with rt_atm_boa_export_free (or manual free).
 *
 * Conventions:
 *   - Fourier reconstruction: I(μ_j, φ) = Σ_m (2-δ_m,0) X^m(μ_j) cos(m·Δφ)
 *     (same as TOA convention; user does reconstruction externally).
 *   - These are *diffuse* Stokes only (direct beam is computed analytically
 *     by the caller as F_sun_TOA × exp(-τ_atm/μ_sun_air) at direction -μ_sun_air).
 *   - F_sun=π normalized (OCRT internal convention).
 * ==================================================================== */
typedef struct {
    int     m_max;             /* max Fourier mode (= opts.fourier_m_max) */
    int     n_mu;              /* number of positive quadrature nodes */
    double *mu_quad_pos;       /* [n_mu] positive μ grid (atm side) */
    double *I_per_m;           /* [(m_max+1) × n_mu] diffuse I^m at BOA, downward */
    double *Q_per_m;           /* same for Q */
    double *U_per_m;           /* same for U */
    /* Atm metadata */
    double  tau_atm_total;     /* total atm optical depth (Rayleigh + aerosol + abs) */
    double  mu_sun_air;        /* cos(SZA) on air side */
    int     allocated;         /* 1 if rt_atm_boa_export_alloc was used (for safe free) */
} rt_atm_boa_export_t;

#endif /* OCRT_V2_RT_TYPES_H */
