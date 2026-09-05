/* OCRT rt_ipss κ on an externally supplied layer profile (profile.csv from
 * make_case.py), for the RTSOS cross-check.  Prints kappa for both anchors. */
#include "rt_ipss.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s profile.csv\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r");
    if (!f) return 1;
    char line[512];
    double hb[512], ht[512], tau[512], om[512];
    int n = 0;
    fgets(line, sizeof line, f);                       /* header */
    while (fgets(line, sizeof line, f) && n < 512) {
        if (sscanf(line, "%lf,%lf,%lf,%lf", &hb[n], &ht[n], &tau[n], &om[n]) == 4) ++n;
    }
    fclose(f);
    /* profile.csv is top-down; rt_ipss wants ascending shells (index 0 = bottom) */
    const double sza_list[6] = { 0, 40, 60, 70, 80, 85 };
    const double vza_list[5] = { 0, 20, 40, 55, 70 };
    const double raa_list[3] = { 0, 90, 180 };
    printf("sza,vza,raa,kappa_surf,kappa_toa\n");
    for (int is = 0; is < 6; ++is) {
        rt_ipss_t G; memset(&G, 0, sizeof G);
        G.nlay = n; G.Re = RT_IPSS_EARTH_RADIUS_KM; G.H = ht[0];
        G.mu_sun = cos(sza_list[is] * M_PI / 180.0); G.n_quad = 4;
        G.h_edge = malloc((n + 1) * sizeof(double));
        G.beta = malloc(n * sizeof(double)); G.omega = malloc(n * sizeof(double));
        for (int j = 0; j < n; ++j) {                 /* j-th from bottom = row n-1-j */
            const int r = n - 1 - j;
            G.h_edge[j] = hb[r];
            G.beta[j] = tau[r] / (ht[r] - hb[r]);
            G.omega[j] = om[r];
        }
        G.h_edge[n] = ht[0];
        rt_ipss_finalize(&G);
        for (int iv = 0; iv < 5; ++iv)
            for (int ir = 0; ir < 3; ++ir) {
                double ks = 1, i1s, i1p; int g;
                rt_ipss_kappa(&G, vza_list[iv], raa_list[ir], &ks, &i1s, &i1p, &g);
                const double a_toa = rt_ipss_single_scatter(&G, vza_list[iv], raa_list[ir],
                                                            RT_IPSS_SUN_SPHERICAL);
                const double b_pp = rt_ipss_single_scatter_pp(&G, vza_list[iv]);
                const double kt = (b_pp > 0 && a_toa > 0) ? a_toa / b_pp : NAN;
                printf("%.1f,%.1f,%.1f,%.9f,%.9f\n", sza_list[is], vza_list[iv],
                       raa_list[ir], ks, kt);
            }
        rt_ipss_free(&G);
    }
    return 0;
}
