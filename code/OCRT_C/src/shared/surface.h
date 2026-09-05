/* ============================================================================
 * shared/surface.h — Vector Fresnel + Cox-Munk surface BRDF kernels
 *
 * Phase 4 (Path 3 graft from V1 vrt_solver.c L1400-3076).
 *
 * Three primitives:
 *   1. surface_flat_fresnel_matrix_raw  — Mueller matrix for flat water
 *      (3x3 in meridian basis, no rotation, no slope weighting).
 *   2. surface_R_coxmunk_trig           — full Cox-Munk Fresnel BRDF kernel
 *      Returns 3x3 Mueller R such that I_up = R · I_dn for one (μ_o, μ_i, φ)
 *      direction pair, including:
 *        - Cox-Munk Gaussian wave slope (isotropic)
 *        - Sancer bistatic shadowing
 *        - Fresnel Mueller at microfacet incidence angle
 *        - Stokes basis rotation (incoming/outgoing meridian → scatter plane)
 *   3. surface_exact_1st_order_scatter  — exact specular peak handling for
 *      1st-order direct sunglint (replaces coarse-grid quadrature; needed
 *      because specular peak is sub-grid).
 *
 * Conventions:
 *   - mu_in_dn > 0  is downward incidence cosine (negate-sign meaning).
 *   - q_convention = 1 (Mishchenko) recommended for vector RT solver
 *     (matches V3/OSPOL.f Q sign).
 *   - sigma_type   = 0 Nakajima-Tanaka, 1 OCRT floor law (default).
 *
 * Dependencies: shared/mat3.h (mat3_zero, mat3_mul, build_rotation_L).
 * ============================================================================ */
#ifndef SHARED_SURFACE_H
#define SHARED_SURFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Compute Cox-Munk slope variance.
 * sigma_type 0: Nakajima-Tanaka (sigma = 0.0731 sqrt(ws), sigma_sq floor 1e-10)
 * sigma_type 1: OCRT floor law (CM54-derived coefficients with explicit low-wind floor): sigma_sq = 0.003 + 0.00512 max(0.01,ws) */
double surface_slope_variance(double wind_speed_ms, int sigma_type);

/* Direct-beam air->water IRRADIANCE transmittance for a Cox-Munk rough surface
 * (slope-integrated facet Fresnel, no shadowing). Returns 1-R(theta_sun) in the
 * flat limit (sigma2->0). Use for wind>0; for flat/wind<=0 use flat-Fresnel. */
double surface_T_aw_coxmunk_direct(double mu_sun_air, double n_water,
                                   double wind_speed, int sigma_type);

/* surface_T_aw_coxmunk_btdf_scalar — scalar (intensity) air->water transmission
 * BTDF for a Cox-Munk rough surface (Walter et al. 2007 microfacet transmission,
 * no shadowing).  f_t(i->o) for incident air dir mu_i and transmitted water dir
 * mu_o separated by azimuth dphi.  Validated: ∫ f_t cos(theta_o) dOmega_o over
 * the water hemisphere == surface_T_aw_coxmunk_direct(mu_i).  Includes the n^2
 * radiance enhancement.  sigma2<=0 (flat) returns 0 (flat is a delta — handle
 * separately via Snell + flat Fresnel). */
double surface_T_aw_coxmunk_btdf_scalar(double mu_i, double mu_o, double dphi,
                                        double n_water, double wind_speed,
                                        int sigma_type);

/* surface_R_ww_coxmunk_direct — water->water (internal) irradiance reflectance
 * for an upwelling beam at mu_up on a Cox-Munk rough surface (slope-integrated
 * water->air Fresnel reflectance incl. TIR; no shadowing).  sigma2->0 gives the
 * flat internal reflectance (=1 for mu_up<mu_crit=TIR).  Use only for wind>0. */
double surface_R_ww_coxmunk_direct(double mu_up_water, double n_water,
                                   double wind_speed, int sigma_type);

/* Flat-water Fresnel Mueller matrix (3x3 in meridian basis).
 *   mu       : abs cosine of incidence angle (>= 0).
 *   n_water  : water refractive index.
 *   q_convention: 0 legacy, 1 Mishchenko.
 *   MF[9]    : output, row-major.
 *
 * Total internal reflection branch returns identity.
 *
 * Backward-compatible wrapper around surface_flat_fresnel_R_matrix_general
 * with (n1, n2) = (1, n_water). Behavior bit-exact to the pre-Phase-B.2
 * implementation.
 */
void surface_flat_fresnel_matrix_raw(double mu, double n_water,
                                      int q_convention, double *MF);

/* ----------------------------------------------------------------------------
 * Generalized flat-interface Fresnel Mueller matrices (Phase B.2 additions)
 *
 * These accept (n1, n2) explicitly so the same routine handles:
 *   air -> water reflection at flat surface : (n1, n2) = (1, n_water)
 *   water-side internal reflection          : (n1, n2) = (n_water, 1)
 *   air -> water radiance transmission      : (n1, n2) = (1, n_water)
 *   water -> air radiance transmission      : (n1, n2) = (n_water, 1)
 *
 * Fresnel core (Snell + s/p amplitudes) is computed once internally and shared
 * between the reflection and transmission outputs. This eliminates the prior
 * inline duplication of Fresnel formulas in surface_flat_fresnel_matrix_raw
 * and surface_R_coxmunk_trig, and guarantees the reflection and transmission
 * branches use identical numerics.
 *
 * Mueller form (meridian basis, [I, Q, U] 3-vector):
 *   M[0][0] = (X_s + X_p) / 2,     X in {R, T_radiance}
 *   M[0][1] = M[1][0] = Q_kernel  (sign per q_convention)
 *   M[1][1] = (X_s + X_p) / 2
 *   M[2][2] = (UU element: amplitude product)
 *
 * Radiance transmittance (M_T) includes the n^2 radiance law:
 *   T_X_radiance = (n2/n1)^2 * T_X_flux,   T_X_flux = (n2*mu_t)/(n1*mu_i) * t_X^2
 * Energy conservation (validated by surface_flat_fresnel_check_energy):
 *   R_X(flux) + T_X(flux) = 1  per polarization (X in {s, p})
 * ---------------------------------------------------------------------------- */

/* Generalized flat-interface reflection Mueller matrix.
 *   mu_i         : abs cos of incidence angle in medium 1 (>= 0)
 *   n1, n2       : refractive index of incidence / transmission medium
 *   q_convention : 0 legacy, 1 Mishchenko
 *   MR[9]        : output, row-major, 3x3
 *
 * TIR (n1 > n2 and sin(theta_i) > n2/n1): returns identity matrix.
 */
void surface_flat_fresnel_R_matrix_general(double mu_i, double n1, double n2,
                                            int q_convention, double *MR);

/* Generalized flat-interface radiance transmission Mueller matrix.
 *   mu_i         : abs cos of incidence angle in medium 1 (>= 0)
 *   n1, n2       : refractive index of incidence / transmission medium
 *   q_convention : 0 legacy, 1 Mishchenko
 *   MT[9]        : output, row-major, 3x3
 *
 * TIR branch: returns zero matrix (no radiation transmitted).
 *
 * Output is in radiance form (n^2 law applied), suitable for direct use as
 * a boundary coupling kernel on Stokes [I, Q, U] vectors in the OCRT SOS
 * solver.
 */
void surface_flat_fresnel_T_matrix(double mu_i, double n1, double n2,
                                    int q_convention, double *MT);

/* Energy conservation diagnostic.
 *   For (mu_i, n1, n2), returns the residual
 *     | 1 - (R_X_flux + T_X_flux) |   summed over X in {s, p}
 *   Should be < 1e-12 for non-TIR cases. Returns 0.0 for TIR (no transmission,
 *   R=1 by definition).
 */
double surface_flat_fresnel_check_energy(double mu_i, double n1, double n2);

/* Cox-Munk + Fresnel BRDF kernel (3x3 Mueller, scatter-plane → meridian rotated).
 *   mu_out, mu_in_dn  : view+incidence cosines (incidence > 0 = downward).
 *   cos_phi, sin_phi  : relative azimuth (phi_view - phi_in).
 *   ws                : wind speed m/s.
 *   sigma_type_int    : 0 Nakajima-Tanaka, 1 OCRT floor law.
 *   n_water           : water refractive index.
 *   q_convention      : 0 legacy, 1 Mishchenko.
 *   R[9]              : output, row-major.
 *
 * R is the BRDF, dimensionless (units of [sr^-1] absorbed via the
 * I_TOA = π · I / (μ_sun F_sun) convention). To use as boundary kernel:
 *   I_up(mu_o) = ∫ R(mu_o, mu_i, φ) · I_dn(mu_i, φ_in) · 2 mu_i dmu_i dφ_in
 */
void surface_R_coxmunk_trig(double mu_out, double mu_in_dn,
                             double cos_phi, double sin_phi,
                             double ws, int sigma_type_int, double n_water,
                             int q_convention, double *R);

/* Cox-Munk + Fresnel BTDF kernel for water -> air TRANSMISSION (3x3 Mueller T).
 *   I_air(mu_out_air) = T · I_water(mu_in_water)   at one (mu2, mu1, φ) triple.
 * Transmission analog of surface_R_coxmunk_trig: refraction-facet geometry
 * (N = n_w p1 - p2), per-facet Fresnel transmission (fresnel core), Walter 2007
 * radiance BTDF with the same NDF (P_slope/cos^4β), Sancer shadowing, and
 * Hovenier rotations. Returns zero matrix on TIR (omega_w > critical angle). */
void surface_T_coxmunk_trig(double mu_out_air, double mu_in_water,
                            double cos_phi, double sin_phi,
                            double ws, int sigma_type_int, double n_water,
                            int q_convention, double *T);

/* Single direct sunglint Mueller kernel for view direction (μ_v, φ_v) and
 * solar direction (μ_0, φ_0). Returns the *direct* surface contribution
 * (1st-order Cox-Munk reflection of solar beam, attenuated by atmospheric
 * transmittance through tau_total).
 *
 * If decouple_sunglint=1 in the caller, this term is what gets subtracted
 * post-hoc from the target-direction radiance.
 *
 *   Output: I_glint[3] (I, Q, U); contribution to TOA reflectance ρ.
 *
 * Note: F_sun is normalized to π in the rho convention used throughout V3.
 *       This function returns ρ_glint such that the contribution to TOA ρ
 *       is exactly I_glint[k] (already premultiplied by atmospheric
 *       transmittance and 4π/μ_v factor).
 */
void surface_direct_sunglint_rho(double mu_v, double phi_v,
                                  double mu_0, double phi_0,
                                  double ws, int sigma_type, double n_water,
                                  double tau_total, int q_convention,
                                  double *rho_glint /* [3] = I,Q,U */);

/* ============================================================================
 * Phase B: Cox-Munk Fourier-mode surface kernel
 *
 *   Compute m-mode Fourier coefficients R^m_kl(μ_o, μ_i) of the Cox-Munk
 *   Fresnel BRDF kernel by numerical φ-quadrature integration over [0, 2π).
 *
 *   For Stokes basis with I,Q ~ cos(mφ) and U ~ sin(mφ) (V3/OSPOL.f conv):
 *
 *     R^m_kl = (1/π) ∫₀^{2π} R_kl(μ_o, μ_i, φ) · cos(mφ) dφ   (k,l ∈ {I,Q})
 *     R^m_UU = (1/π) ∫₀^{2π} R_UU(μ_o, μ_i, φ) · cos(mφ) dφ
 *     R^m_kl = (1/π) ∫₀^{2π} R_kl(μ_o, μ_i, φ) · sin(mφ) dφ   (one of k,l = U)
 *
 *   For m=0, the (1/π) prefactor becomes (1/2π) due to the Fourier
 *   normalization conventions used in V3 (consistent with the (2-δ_{m,0})
 *   reconstruction in rt_solver_reconstruct_phi).
 *
 *   The output stores all 9 Mueller entries per (μ_o, μ_i) pair, packed
 *   row-major: kernel[(j_o*n_mu_pairs + j_i)*9 + (k*3 + l)].
 *
 *   Inputs:
 *     mu_o[n_o]    : upward μ values (positive, > 0)
 *     mu_i[n_i]    : downward μ values, *positive value of |μ|* (kernel
 *                    routine takes mu_in_dn directly)
 *     m            : Fourier mode (m >= 0)
 *     n_phi_quad   : number of φ-quadrature points (uniform on [0, 2π))
 *     ws, sigma_type, n_water, q_convention : Cox-Munk parameters
 *
 *   Output:
 *     R_m[n_o * n_i * 9] : packed Mueller kernel m-mode coefficients
 *
 *   Returns 0 on success.
 */
int surface_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                    const double *mu_i, int n_i,
                                    int m, int n_phi_quad,
                                    double ws, int sigma_type,
                                    double n_water, int q_convention,
                                    double *R_m);

/* v1.11-speed S1 (2026-08-25): declare the true Fourier-mode upper bound of
 * the NEXT kernel-building m loop so the all-m batched kernel cache builds
 * only modes 0..m_max instead of the full compile ceiling.  Thread-local;
 * values outside the valid range restore the uncapped default.  Purely a
 * performance hint - built modes are bit-identical to the per-m evaluation. */
void surface_fkc_set_build_mmax(int m_max);

/* WATER-SIDE internal-reflection analogues (B2, 2026-06-03). Same signatures
 * and conventions as the air-side R_coxmunk_trig / coxmunk_fourier_kernel, with
 * the facet Fresnel changed to water(n_water)->air(1.0) (internal, TIR).
 * mu_out/mu_o = reflected DOWNWARD |mu| (into water); mu_in_dn/mu_i = incident
 * UPWARD |mu| (in-water upwelling). Used by the SOS internal-reflection
 * feedback at wind > 0 (rough surface). */
void surface_R_ww_coxmunk_trig(double mu_out, double mu_in_dn,
                               double cos_phi, double sin_phi,
                               double ws, int sigma_type_int, double n_water,
                               int q_convention, double *R);
int surface_R_ww_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                        const double *mu_i, int n_i,
                                        int m, int n_phi_quad,
                                        double ws, int sigma_type,
                                        double n_water, int q_convention,
                                        double *R_m);

/* WATER->AIR TRANSMISSION Fourier-mode Mueller kernel (Milestone 2b, 2026-06-04).
 * Transmission analog of surface_R_ww_coxmunk_fourier_kernel (same normalization
 * + cos/sin layout), built by phi-integrating surface_T_coxmunk_trig. mu_o = AIR
 * outgoing |mu|; mu_i = in-water upwelling |mu|. Packed T_m[(j_o*n_i+j_i)*9+k*3+l].
 * Used by the rough (wind>0) reverse coupling rt_air_water_couple_water_to_atm. */
/* v1.10 B-0c.1: air->water rough-surface transmission Fourier matrix.
 * mu_o = IN-WATER outgoing nodes, mu_i = AIR incident nodes.  Intensity
 * block is anchored to the TESTED scalar surface_T_aw_coxmunk_btdf_scalar;
 * the polarization structure is inherited from the wa Mueller by
 * reciprocity (D T_wa^T D, D=diag(1,1,-1) on IQU).  Flat (wind<=0) is a
 * delta lobe and returns -2 (caller keeps its flat path). */
int surface_T_aw_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                        const double *mu_i, int n_i,
                                        int m, int n_phi_quad,
                                        double ws, int sigma_type,
                                        double n_water, int q_convention,
                                        double *T_m);

int surface_T_wa_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                        const double *mu_i, int n_i,
                                        int m, int n_phi_quad,
                                        double ws, int sigma_type,
                                        double n_water, int q_convention,
                                        double *T_m);

#ifdef __cplusplus
}
#endif

#endif /* SHARED_SURFACE_H */
void surface_T_aw_coxmunk_trig_test(double,double,double,double,double,int,double,int,double*);

/* Multi-encounter water->air cascade (surface_multibounce.c).  Exact energy
 * closure T_h = T_up + W_ret within the uncorrelated-facet model; used for
 * the option-2 interface closure (2026-07-10).  Deterministic fixed seed. */
int surface_wa_multibounce_cascade(double n_water, double wind_speed,
                                   double mu_w, long n_rays,
                                   double *T_h, double *T_up, double *W_ret,
                                   double *esc_hist, double *escd_hist,
                                   double *ret_hist, int nb);
