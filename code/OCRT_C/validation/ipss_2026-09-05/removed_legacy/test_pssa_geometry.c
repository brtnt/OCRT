#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rt_atm.h"
#include "rt_pssa.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int close_rel(double a, double b, double rtol, double atol)
{
    const double d = fabs(a - b);
    const double s = fmax(fabs(a), fabs(b));
    return d <= atol + rtol * s;
}

int main(void)
{
    const int n_layers = 400;
    const double sza_deg = 85.0;
    const double theta0 = sza_deg * M_PI / 180.0;
    const double mu0 = cos(theta0);
    rt_atm_t atm = {0};

    if (rt_atm_alloc(&atm, n_layers, 1) != 0) {
        fprintf(stderr, "rt_atm_alloc failed\n");
        return 2;
    }
    if (rt_atm_build_rayleigh(&atm, 0.318, 0.0279, mu0,
                              RT_RAYLEIGH_MODEL_BODHAINE_1999) != 0) {
        fprintf(stderr, "rt_atm_build_rayleigh failed\n");
        rt_atm_free(&atm);
        return 2;
    }
    if (rt_pssa_apply(&atm) != 0) {
        fprintf(stderr, "rt_pssa_apply failed\n");
        rt_atm_free(&atm);
        return 2;
    }

    double max_eq8_ratio = 0.0;
    double max_eq8_angle = 0.0;
    double max_ch = 0.0;
    int fail = 0;
    for (int k = 0; k <= n_layers; ++k) {
        const double alpha = atm.pssa_alpha[k];
        const double beta = atm.pssa_beta[k];
        const double z = atm.z_km_level[k];
        const double lhs = sin(alpha) / fmax(sin(beta), 1.0e-300);
        const double rhs = (RT_PSSA_EARTH_RADIUS_KM + z) /
                           RT_PSSA_EARTH_RADIUS_KM;
        const double ratio_err = fabs(lhs - rhs) / rhs;
        const double angle_err = fabs(2.0 * alpha - beta - theta0);
        const double ch_expected = 0.5 * exp(-atm.pssa_xi_dn[k]);
        const double ch_err = fabs(atm.ch[k] - ch_expected);
        if (ratio_err > max_eq8_ratio) max_eq8_ratio = ratio_err;
        if (angle_err > max_eq8_angle) max_eq8_angle = angle_err;
        if (ch_err > max_ch) max_ch = ch_err;
        if (!(beta >= 0.0 && beta <= 0.5 * M_PI)) fail = 1;
        if (k > 0 && atm.pssa_xi_dn[k] < atm.pssa_xi_dn[k - 1]) fail = 1;
    }

    if (!close_rel(atm.pssa_alpha[n_layers], theta0, 0.0, 2.0e-13) ||
        !close_rel(atm.pssa_beta[n_layers], theta0, 0.0, 4.0e-13) ||
        !close_rel(atm.pssa_xi_refl[n_layers], atm.pssa_xi_dn[n_layers],
                   2.0e-13, 2.0e-14) ||
        max_eq8_ratio > 3.0e-12 || max_eq8_angle > 3.0e-13 ||
        max_ch > 2.0e-16) {
        fail = 1;
    }

    printf("PSSA_GEOMETRY max_eq8_ratio=%.17g max_eq8_angle=%.17g "
           "max_ch_abs=%.17g xi_surface=%.17g\n",
           max_eq8_ratio, max_eq8_angle, max_ch,
           atm.pssa_xi_dn[n_layers]);

    rt_atm_free(&atm);
    if (fail) {
        fprintf(stderr, "PSSA geometry invariant failure\n");
        return 1;
    }
    puts("PASS: He et al. Eq. (8), surface and attenuation invariants");
    return 0;
}
