#include "rt_phase_moments.h"

#include <math.h>
#include <stddef.h>

void rt_phase_spin2_m0(int l_max, double mu, double *out)
{
    if (out == NULL || l_max < 0) return;
    out[0] = 0.0;
    if (l_max == 0) return;
    out[1] = 0.0;
    if (l_max == 1) return;

    out[2] = 3.0 * (1.0 - mu * mu) / (2.0 * sqrt(6.0));
    for (int ell = 2; ell < l_max; ++ell) {
        const double scale = (2.0 * ell + 1.0) /
                             sqrt((double)(ell - 1) * (ell + 3));
        const double back = sqrt((double)(ell - 2) * (ell + 2)) /
                            (2.0 * ell + 1.0);
        out[ell + 1] = scale * (mu * out[ell] - back * out[ell - 1]);
    }
}

int rt_phase_complete_spherical(int l_max,
                                const double *beta11,
                                const double *beta22,
                                const double *delta33,
                                double *alpha,
                                double *zeta)
{
    if (l_max < 0 || beta11 == NULL || beta22 == NULL || delta33 == NULL ||
        alpha == NULL || zeta == NULL) {
        return -1;
    }

    for (int ell = 0; ell <= l_max; ++ell) {
        alpha[ell] = 0.0;
        zeta[ell] = 0.0;
    }

    for (int ell = 2; ell <= l_max; ++ell) {
        const double e = (double)ell;
        const double lower_scale = e * (e - 1.0) /
                                   ((e + 1.0) * (e + 2.0));
        const double coupling = 4.0 * (2.0 * e + 1.0) /
                                (e * (e - 1.0) * (e + 1.0) * (e + 2.0));
        const double base = (e - 1.0) * (e - 1.0);

        double same_beta = 0.0;
        double same_delta = 0.0;
        double opposite_beta = 0.0;
        double opposite_delta = 0.0;

        for (int order = ell - 1; order >= 0; --order) {
            const int gap = ell - order;
            const double weight = base -
                1.5 * (double)(gap - 1) * (double)(ell + order);
            if ((gap & 1) == 0) {
                same_beta += weight * beta11[order];
                same_delta += weight * delta33[order];
            } else {
                opposite_beta += weight * beta22[order];
                opposite_delta += weight * delta33[order];
            }
        }

        zeta[ell] = lower_scale * delta33[ell] -
                    coupling * (same_delta - opposite_beta);
        alpha[ell] = lower_scale * beta11[ell] -
                     coupling * (same_beta - opposite_delta);
    }
    return 0;
}

int rt_phase_normalize_vector(int l_max,
                              double *beta,
                              double *gamma,
                              double *alpha,
                              double *zeta)
{
    if (l_max < 0 || beta == NULL || gamma == NULL || alpha == NULL ||
        zeta == NULL) {
        return -1;
    }
    const double norm = beta[0];
    if (norm == 0.0 || !isfinite(norm)) return -1;
    for (int ell = 0; ell <= l_max; ++ell) {
        alpha[ell] /= norm;
        beta[ell] /= norm;
        gamma[ell] /= norm;
        zeta[ell] /= norm;
    }
    return 0;
}
