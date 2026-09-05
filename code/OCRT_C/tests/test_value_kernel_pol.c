#include "rt_atm.h"
#include "rt_value_phase.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static size_t kidx(int m, int j, int k, int n_mu)
{
    const int kw = 2 * n_mu + 1;
    return ((size_t)m * (size_t)(n_mu + 1) + (size_t)j) *
           (size_t)kw + (size_t)(k + n_mu);
}

int main(void)
{
    enum { N = 361, NMU = 6, MCOUNT = 5, NPHI = 512 };
    double th[N], p11[N], p12[N], p33[N];
    for (int i = 0; i < N; ++i) {
        th[i] = 0.5 * (double)i;
        const double mu = cos(th[i] * M_PI / 180.0);
        p11[i] = 0.75 * (1.0 + mu * mu);
        p12[i] = -0.75 * (1.0 - mu * mu);
        p33[i] = 1.5 * mu;
    }

    rt_value_phase_interp_t phase;
    memset(&phase, 0, sizeof phase);
    if (rt_value_phase_interp_build(&phase, th, p11, p12, p33, N, 0) != 0) {
        fputs("FAIL build theta-linear phase\n", stderr);
        return 2;
    }
    if (!phase.is_fr631 || phase.n != 631) {
        fprintf(stderr, "FAIL legacy 361 grid was not canonicalized: n=%d fr631=%d\n",
                phase.n, phase.is_fr631);
        return 3;
    }

    double max_node_error = 0.0;
    for (int i = 0; i < N; ++i) {
        double a, b, c;
        rt_value_phase_interp_eval_theta(&phase, th[i], &a, &b, &c);
        const double expected[3] = {
            p11[i] / phase.norm_before,
            p12[i] / phase.norm_before,
            p33[i] / phase.norm_before
        };
        const double got[3] = {a, b, c};
        for (int q = 0; q < 3; ++q) {
            const double e = fabs(got[q] - expected[q]);
            if (e > max_node_error) max_node_error = e;
        }
    }

    for (int i = 0; i < phase.n; ++i) {
        if (!(phase.p11[i] >= 0.0) || fabs(phase.p12[i]) > phase.p11[i] + 1e-13 ||
            fabs(phase.p33[i]) > phase.p11[i] + 1e-13) {
            fprintf(stderr, "FAIL physical bound at node %d\n", i);
            return 4;
        }
    }

    rt_atm_t atm;
    memset(&atm, 0, sizeof atm);
    if (rt_atm_alloc(&atm, 1, NMU) != 0) return 5;
    if (rt_atm_build_rayleigh(&atm, 0.1, 0.0279,
                              cos(30.0 * M_PI / 180.0),
                              RT_RAYLEIGH_MODEL_BODHAINE_1999) != 0) return 6;

    const int kw = 2 * NMU + 1;
    const size_t plane = (size_t)(NMU + 1) * (size_t)kw;
    const size_t all = (size_t)MCOUNT * plane;
    double *scalar = (double *)calloc(all, sizeof(double));
    double *buf = (double *)calloc(6u * all, sizeof(double));
    double *buf2 = (double *)calloc(6u * all, sizeof(double));
    if (!scalar || !buf || !buf2) return 7;
    double *pf = buf;
    double *gr = pf + all;
    double *gt = gr + all;
    double *arr = gt + all;
    double *art = arr + all;
    double *att = art + all;

    if (rt_value_phase_fourier_scalar_allm(&atm, MCOUNT, &phase, NPHI, scalar) != 0 ||
        rt_aerosol_value_phase_fourier_pol_allm(
            &atm, MCOUNT, &phase, NPHI, pf, gr, gt, arr, art, att) != 0 ||
        rt_aerosol_value_phase_fourier_pol_allm(
            &atm, MCOUNT, &phase, NPHI,
            buf2, buf2 + all, buf2 + 2u * all, buf2 + 3u * all,
            buf2 + 4u * all, buf2 + 5u * all) != 0) {
        fputs("FAIL direct Fourier build\n", stderr);
        return 8;
    }

    double max_scalar_vector = 0.0;
    double max_repeat = 0.0;
    for (int m = 0; m < MCOUNT; ++m) {
        for (int j = 0; j <= NMU; ++j) {
            for (int k = -NMU; k <= NMU; ++k) {
                const size_t z = kidx(m, j, k, NMU);
                const double e = fabs(scalar[z] - pf[z]);
                if (e > max_scalar_vector) max_scalar_vector = e;
            }
        }
    }
    for (size_t i = 0; i < 6u * all; ++i) {
        const double e = fabs(buf[i] - buf2[i]);
        if (e > max_repeat) max_repeat = e;
        if (!isfinite(buf[i])) {
            fprintf(stderr, "FAIL nonfinite kernel at %zu\n", i);
            return 9;
        }
    }

    printf("n=%d fr631=%d norm_before=%.17g bb_b=%.17g g=%.17g "
           "max_node_error=%.3e max_scalar_vector=%.3e max_repeat=%.3e\n",
           phase.n, phase.is_fr631, phase.norm_before, phase.bb_b_ratio,
           phase.g_asym, max_node_error, max_scalar_vector, max_repeat);

    free(scalar); free(buf); free(buf2);
    rt_atm_free(&atm);
    rt_value_phase_interp_free(&phase);
    return (max_node_error <= 2e-14 && max_scalar_vector <= 2e-14 &&
            max_repeat == 0.0) ? 0 : 10;
}
