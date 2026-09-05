#ifndef OCRT_RT_PHASE_MOMENTS_H
#define OCRT_RT_PHASE_MOMENTS_H

/* Mathematical support for the generalized spherical-function expansion of
 * a reciprocal spherical-particle 3-Stokes scattering matrix.
 *
 * These routines are derived from the published generalized Legendre/Wigner
 * basis identities and the de Haan-Hovenier vector phase expansion.  They do
 * not depend on any external radiative-transfer source implementation.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* m=0 spin-2 basis used to project P12.  out[0..l_max] is filled. */
void rt_phase_spin2_m0(int l_max, double mu, double *out);

/* Complete the alpha/zeta coefficient families for spherical particles from
 * ordinary moments of P11/P22/P33.  beta22 may alias beta11.  alpha and zeta
 * are initialized by this function; orders 0 and 1 are zero. */
int rt_phase_complete_spherical(int l_max,
                                const double *beta11,
                                const double *beta22,
                                const double *delta33,
                                double *alpha,
                                double *zeta);

/* Normalize all four coefficient families by beta[0]. */
int rt_phase_normalize_vector(int l_max,
                              double *beta,
                              double *gamma,
                              double *alpha,
                              double *zeta);

#ifdef __cplusplus
}
#endif

#endif /* OCRT_RT_PHASE_MOMENTS_H */
