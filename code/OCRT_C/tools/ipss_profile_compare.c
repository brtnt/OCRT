/* ============================================================================
 * tools/ipss_profile_compare.c
 * Single-scatter PSS-vs-IPSS error over the paper's geometry grid, for two
 * vertical profiles:
 *   (a) UNIFORM density 0-100 km  — the configuration Zhai & Hu (2022) Sec. 3
 *       validated against the Korkin et al. (2020) benchmark;
 *   (b) EXPONENTIAL (scale height 8 km) — what a real atmosphere, and hence an
 *       OCRT production run on the US62 grid, actually looks like.
 * "PSS" here is the literature convention: spherical line of sight, solar beam
 * evaluated on the nadir column.  (OCRT's legacy Chapman PSSA is weaker still:
 * it also leaves the viewing path plane-parallel.)
 * Usage: ipss_profile_compare <tau_total> <sza_deg>
 * ========================================================================== */
#include "rt_ipss.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void build(rt_ipss_t *G, double tau_t, double sza, int nlay, int exp_prof)
{
    memset(G, 0, sizeof(*G));
    G->nlay = nlay; G->Re = RT_IPSS_EARTH_RADIUS_KM; G->H = 100.0;
    G->mu_sun = cos(sza * M_PI / 180.0); G->n_quad = 400;
    G->h_edge = (double *)malloc((size_t)(nlay + 1) * sizeof(double));
    G->beta   = (double *)malloc((size_t)nlay * sizeof(double));
    G->omega  = (double *)malloc((size_t)nlay * sizeof(double));
    for (int k = 0; k <= nlay; ++k) G->h_edge[k] = 100.0 * (double)k / nlay;
    double sum = 0.0;
    for (int j = 0; j < nlay; ++j) {
        const double hm = 0.5 * (G->h_edge[j] + G->h_edge[j + 1]);
        G->beta[j] = exp_prof ? exp(-hm / 8.0) : 1.0;
        sum += G->beta[j] * (G->h_edge[j + 1] - G->h_edge[j]);
        G->omega[j] = 1.0;
    }
    for (int j = 0; j < nlay; ++j) G->beta[j] *= tau_t / sum;
    rt_ipss_finalize(G);
}

int main(int argc, char **argv)
{
    const double tau = (argc > 1) ? atof(argv[1]) : 0.25;
    const double sza = (argc > 2) ? atof(argv[2]) : 84.26;
    const double raa[3] = { 180.0, 90.0, 0.0 };   /* paper phi_v = 0, 90, 180 */
    printf("tau=%.2f  SZA=%.2f  (PSS-IPSS)/IPSS [%%] at single-scatter level\n",
           tau, sza);
    printf("%-12s %10s %10s %10s | %10s %10s %10s\n", "",
           "UNIF phi0", "UNIF phi90", "UNIF phi180",
           "EXP8 phi0", "EXP8 phi90", "EXP8 phi180");
    double umin = 1e9, umax = -1e9, emin = 1e9, emax = -1e9;
    for (double vza = 0.0; vza <= 70.001; vza += 10.0) {
        double u[3], e[3];
        rt_ipss_t G;
        build(&G, tau, sza, 100, 0);
        for (int i = 0; i < 3; ++i) {
            const double a = rt_ipss_single_scatter(&G, vza, raa[i], RT_IPSS_SUN_NADIR);
            const double b = rt_ipss_single_scatter(&G, vza, raa[i], RT_IPSS_SUN_SPHERICAL);
            u[i] = (b > 0.0) ? 100.0 * (a / b - 1.0) : 0.0;
            if (u[i] < umin) umin = u[i];
            if (u[i] > umax) umax = u[i];
        }
        rt_ipss_free(&G);
        build(&G, tau, sza, 100, 1);
        for (int i = 0; i < 3; ++i) {
            const double a = rt_ipss_single_scatter(&G, vza, raa[i], RT_IPSS_SUN_NADIR);
            const double b = rt_ipss_single_scatter(&G, vza, raa[i], RT_IPSS_SUN_SPHERICAL);
            e[i] = (b > 0.0) ? 100.0 * (a / b - 1.0) : 0.0;
            if (e[i] < emin) emin = e[i];
            if (e[i] > emax) emax = e[i];
        }
        rt_ipss_free(&G);
        printf("VZA %5.1f    %10.2f %10.2f %10.2f | %10.2f %10.2f %10.2f\n",
               vza, u[0], u[1], u[2], e[0], e[1], e[2]);
    }
    printf("range        UNIFORM %+.1f .. %+.1f %%     EXPONENTIAL %+.1f .. %+.1f %%\n",
           umin, umax, emin, emax);
    return 0;
}
