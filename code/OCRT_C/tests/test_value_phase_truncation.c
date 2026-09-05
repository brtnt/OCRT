#include "rt_phase_fr631.h"
#include "rt_value_phase.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int nearly(double a, double b, double atol, double rtol)
{
    const double scale = fmax(fabs(a), fabs(b));
    return fabs(a - b) <= atol + rtol * scale;
}

int main(void)
{
    double theta[RT_FR631_N_ANGLE];
    double p11[RT_FR631_N_ANGLE];
    double p12[RT_FR631_N_ANGLE];
    double p33[RT_FR631_N_ANGLE];
    rt_fr631_fill_theta(theta);

    for (int i = 0; i < RT_FR631_N_ANGLE; ++i) {
        const double x = theta[i] / 1.2;
        p11[i] = 1.0 + 4.0e4 * exp(-x * x);
        p12[i] = -0.28 * p11[i] * sin(theta[i] * M_PI / 180.0);
        p33[i] =  0.72 * p11[i];
    }

    rt_value_phase_interp_t raw = {0}, trunc = {0};
    if (rt_value_phase_interp_build(&raw, theta, p11, p12, p33,
                                    RT_FR631_N_ANGLE, 0) != 0) {
        fprintf(stderr, "raw phase build failed\n");
        return 1;
    }
    double A = 0.0;
    if (rt_value_phase_interp_loglinear_truncate(
            &raw, 0.85, 0.92, 0.1, &trunc, &A) != 0) {
        fprintf(stderr, "truncation failed\n");
        rt_value_phase_interp_free(&raw);
        return 1;
    }
    const double f = 0.5 * A;
    if (!(A > 0.1 && A < 2.0)) {
        fprintf(stderr, "unexpected A %.17g\n", A);
        return 1;
    }
    if (!nearly(trunc.norm_before, 1.0 - f, 2e-12, 2e-12)) {
        fprintf(stderr, "mass mismatch norm0=%.17g expected=%.17g\n",
                trunc.norm_before, 1.0 - f);
        return 1;
    }
    if (!nearly(raw.bb_b_ratio, (1.0 - f) * trunc.bb_b_ratio,
                3e-12, 3e-12)) {
        fprintf(stderr, "backscatter mismatch raw=%.17g residual=%.17g\n",
                raw.bb_b_ratio, (1.0 - f) * trunc.bb_b_ratio);
        return 1;
    }
    for (int i = 0; i < trunc.n; ++i) {
        if (!(trunc.p11[i] >= 0.0) || fabs(trunc.p12[i]) > trunc.p11[i] + 1e-13 ||
            fabs(trunc.p33[i]) > trunc.p11[i] + 1e-13) {
            fprintf(stderr, "physical bound failed at %d\n", i);
            return 1;
        }
    }

    double iso11[RT_FR631_N_ANGLE], iso12[RT_FR631_N_ANGLE], iso33[RT_FR631_N_ANGLE];
    for (int i = 0; i < RT_FR631_N_ANGLE; ++i) {
        iso11[i] = 1.0; iso12[i] = 0.0; iso33[i] = 1.0;
    }
    rt_value_phase_interp_t iso = {0}, iso_t = {0};
    double A_iso = -1.0;
    if (rt_value_phase_interp_build(&iso, theta, iso11, iso12, iso33,
                                    RT_FR631_N_ANGLE, 0) != 0 ||
        rt_value_phase_interp_loglinear_truncate(
            &iso, 0.85, 0.92, 0.1, &iso_t, &A_iso) != 0) {
        fprintf(stderr, "isotropic no-op test failed\n");
        return 1;
    }
    const size_t nbytes = (size_t)iso.n * sizeof(double);
    if (A_iso != 0.0 || iso.n != iso_t.n || iso.is_fr631 != iso_t.is_fr631 ||
        iso.norm_before != iso_t.norm_before || iso.g_asym != iso_t.g_asym ||
        iso.bb_b_ratio != iso_t.bb_b_ratio ||
        memcmp(iso.theta_deg, iso_t.theta_deg, nbytes) != 0 ||
        memcmp(iso.p11, iso_t.p11, nbytes) != 0 ||
        memcmp(iso.p12, iso_t.p12, nbytes) != 0 ||
        memcmp(iso.p33, iso_t.p33, nbytes) != 0) {
        fprintf(stderr, "isotropic cap was not an exact no-op A=%.17g\n", A_iso);
        return 1;
    }

    printf("PASS A=%.15g f=%.15g raw_bb_b=%.15g residual_bb_b=%.15g\n",
           A, f, raw.bb_b_ratio, trunc.bb_b_ratio);
    rt_value_phase_interp_free(&raw);
    rt_value_phase_interp_free(&trunc);
    rt_value_phase_interp_free(&iso);
    rt_value_phase_interp_free(&iso_t);
    return 0;
}
