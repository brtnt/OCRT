#ifndef OCRT_RT_MOLECULAR_PHASE_H
#define OCRT_RT_MOLECULAR_PHASE_H

/* Depolarized molecular (Rayleigh) phase representation used by the OCRT
 * 3-Stokes Fourier solver.  The formulas follow directly from the normalized
 * molecular Mueller matrix and its l=0,2 generalized spherical-function
 * expansion; no reference-program source layout is assumed. */

typedef struct {
    double beta0;
    double beta2;
    double gamma2;
    double alpha2;
} rt_molecular_phase_coeffs_t;

/* Convert depolarization ratio delta in [0,1) to the nonzero phase moments. */
int rt_molecular_phase_coefficients(double depolarization,
                                    rt_molecular_phase_coeffs_t *out);

/* Evaluate the m=0, l=2 scalar and real spin-2 basis values at signed cosine
 * mu.  The third basis is zero by symmetry for m=0. */
void rt_molecular_phase_l2_m0(double mu,
                              double *scalar_p2,
                              double *linear_r2,
                              double *linear_t2);

#endif /* OCRT_RT_MOLECULAR_PHASE_H */
