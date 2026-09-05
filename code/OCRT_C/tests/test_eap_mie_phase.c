#include "rt_eap_mie_phase.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int check_norm(const double *p11, const double *theta, size_t n)
{
    double s = 0.0;
    for (size_t j = 0; j + 1 < n; ++j) {
        double t0 = theta[j] * M_PI / 180.0;
        double t1 = theta[j + 1] * M_PI / 180.0;
        s += 0.5 * (t1 - t0) * (p11[j] * sin(t0) + p11[j + 1] * sin(t1));
    }
    s *= 0.5;
    if (fabs(s - 1.0) > 2.0e-12) {
        fprintf(stderr, "normalization failure: %.17g\n", s);
        return 1;
    }
    return 0;
}

int main(void)
{
    enum { NTH = 361, NW = 2 };
    double theta[NTH], wl[NW] = {412.0, 443.0};
    double *p11 = calloc(NW * NTH, sizeof *p11);
    double *p12 = calloc(NW * NTH, sizeof *p12);
    double *p33 = calloc(NW * NTH, sizeof *p33);
    double *q11 = calloc(NW * NTH, sizeof *q11);
    double *q12 = calloc(NW * NTH, sizeof *q12);
    double *q33 = calloc(NW * NTH, sizeof *q33);
    if (!p11 || !p12 || !p33 || !q11 || !q12 || !q33) return 2;
    for (int j = 0; j < NTH; ++j) theta[j] = 0.5 * (double)j;

    int rc = rt_eap_mie_phase_compute(11, wl, NW, theta, NTH, p11, p12, p33);
    if (rc != RT_EAP_PHASE_OK) {
        fprintf(stderr, "species 11 compute failed: %d\n", rc);
        return 3;
    }
    rc = rt_eap_mie_phase_compute(11, wl, NW, theta, NTH, q11, q12, q33);
    if (rc != RT_EAP_PHASE_OK) return 4;
    if (memcmp(p11, q11, NW * NTH * sizeof *p11) ||
        memcmp(p12, q12, NW * NTH * sizeof *p12) ||
        memcmp(p33, q33, NW * NTH * sizeof *p33)) {
        fprintf(stderr, "repeat call is not bit-identical\n");
        return 5;
    }
    for (int iw = 0; iw < NW; ++iw) {
        if (check_norm(p11 + iw * NTH, theta, NTH)) return 6;
        for (int j = 0; j < NTH; ++j) {
            size_t k = (size_t)iw * NTH + (size_t)j;
            if (!(p11[k] >= 0.0) || fabs(p12[k]) > p11[k] + 1e-12 * p11[0] ||
                fabs(p33[k]) > p11[k] + 1e-12 * p11[0]) {
                fprintf(stderr, "physicality failure iw=%d j=%d\n", iw, j);
                return 7;
            }
        }
    }
    if (rt_eap_mie_phase_compute(17, wl, NW, theta, NTH, q11, q12, q33) !=
        RT_EAP_PHASE_ESPECIES) return 8;
    double bad_theta[NTH]; memcpy(bad_theta, theta, sizeof bad_theta); bad_theta[NTH-1] = 179.5;
    if (rt_eap_mie_phase_compute(11, wl, NW, bad_theta, NTH, q11, q12, q33) !=
        RT_EAP_PHASE_ERANGE) return 9;

    printf("PASS: EAP phase API normalization, physicality, determinism and errors\n");
    free(p11); free(p12); free(p33); free(q11); free(q12); free(q33);
    return 0;
}
