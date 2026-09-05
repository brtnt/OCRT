#include "rt_molecular_phase.h"

#include <math.h>

int rt_molecular_phase_coefficients(double depolarization,
                                    rt_molecular_phase_coeffs_t *out)
{
    if (!out || !(depolarization >= 0.0 && depolarization < 1.0)) return -1;

    /* For the normalized depolarized molecular Mueller matrix, the l=2
     * anisotropy factor is 2(1-delta)/(2+delta). */
    const double anisotropy =
        2.0 * (1.0 - depolarization) / (2.0 + depolarization);

    out->beta0  = 1.0;
    out->beta2  = 0.5 * anisotropy;
    out->gamma2 = -anisotropy * sqrt(1.5);
    out->alpha2 = 3.0 * anisotropy;
    return 0;
}

void rt_molecular_phase_l2_m0(double mu,
                              double *scalar_p2,
                              double *linear_r2,
                              double *linear_t2)
{
    const double mu2 = mu * mu;
    const double transverse2 = 1.0 - mu2;

    if (scalar_p2) *scalar_p2 = 0.5 * (3.0 * mu2 - 1.0);
    if (linear_r2) *linear_r2 =
        3.0 * transverse2 / (2.0 * sqrt(6.0));
    if (linear_t2) *linear_t2 = 0.0;
}
