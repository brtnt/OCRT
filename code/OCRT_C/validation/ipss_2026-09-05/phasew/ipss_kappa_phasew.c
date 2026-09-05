/* κ on the solver's own layer profile (OCRT_IPSS_PROFILE_DUMP), three source
 * weightings, same attenuation (total extinction incl. gas):
 *   A  OCRT current : w_j = (dtau - tau_abs_gas)/dtau            (phase-free)
 *   B  paper Eq.(7) : w_j = ydel_j*P_R(Theta) + xdel_j*P_A(Theta) (phase-weighted)
 *   C  scattering   : w_j = ydel_j + xdel_j                        (phase-free, aerosol
 *                                                                   absorption excluded)
 * P_A from the aerosol P11 Legendre moments in the dump, P_R depolarized Rayleigh.
 * Theta is constant along a straight line of sight under a parallel beam, so the
 * phase functions enter only through the height-dependent mixture ratio. */
#include "rt_ipss.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double legendre_sum(const double *b, int L, double x)
{
    double p0 = 1.0, p1 = x, s = b[0];
    if (L >= 1) s += b[1] * x;
    for (int l = 2; l <= L; ++l) {
        const double p2 = ((2.0 * l - 1.0) * x * p1 - (l - 1.0) * p0) / l;
        s += b[l] * p2; p0 = p1; p1 = p2;
    }
    return s;
}

int main(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: %s prof.csv vza raa [sza...]\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) return 1;
    const double vza = atof(argv[2]), raa = atof(argv[3]);
    char line[65536]; double depol = 0.0; int L = 0; double *bl = NULL;
    double zt[2048], zb[2048], dt[2048], xd[2048], yd[2048], ta[2048]; int n = 0;
    while (fgets(line, sizeof line, f)) {
        if (!strncmp(line, "# depol=", 8)) { sscanf(line, "# depol=%lf L_max=%d", &depol, &L); continue; }
        if (!strncmp(line, "# betal_aer:", 12)) {
            bl = calloc((size_t)L + 1, sizeof(double)); char *p = line + 12;
            for (int l = 0; l <= L; ++l) { bl[l] = strtod(p, &p); }
            continue;
        }
        if (line[0] == 'k') continue;
        int k; if (sscanf(line, "%d,%lf,%lf,%lf,%lf,%lf,%lf", &k, &zt[n], &zb[n], &dt[n], &xd[n], &yd[n], &ta[n]) == 7) ++n;
    }
    fclose(f);
    const double g = depol / (2.0 - depol);
    printf("# n=%d depol=%g L_max=%d  vza=%g raa=%g\n", n, depol, L, vza, raa);
    printf("sza,Theta_deg,P_R,P_A,kappa_A_ocrt,kappa_B_phase,kappa_C_scat\n");
    for (int ia = 4; ia < argc; ++ia) {
        const double sza = atof(argv[ia]);
        const double mv = cos(vza * M_PI / 180), ms = cos(sza * M_PI / 180);
        const double fp = (180.0 - raa) * M_PI / 180;
        const double cth = -(sin(sza * M_PI / 180) * sin(vza * M_PI / 180) * cos(fp) + ms * mv);
        const double PR = 3.0 / (4.0 * (1 + 2 * g)) * ((1 + 3 * g) + (1 - g) * cth * cth);
        const double PA = bl ? legendre_sum(bl, L, cth) : 0.0;
        double kap[3];
        for (int mode = 0; mode < 3; ++mode) {
            rt_ipss_t G; memset(&G, 0, sizeof G);
            G.nlay = n; G.Re = RT_IPSS_EARTH_RADIUS_KM; G.mu_sun = ms; G.n_quad = 4;
            G.h_edge = malloc((n + 1) * sizeof(double));
            G.beta = malloc(n * sizeof(double)); G.omega = malloc(n * sizeof(double));
            for (int j = 0; j < n; ++j) {          /* ascending j <- dump row n-1-j */
                const int r = n - 1 - j;
                G.h_edge[j] = zb[r];
                G.beta[j] = dt[r] / (zt[r] - zb[r]);
                if (mode == 0)      G.omega[j] = (dt[r] - ta[r]) / dt[r];
                else if (mode == 1) G.omega[j] = yd[r] * PR + xd[r] * PA;
                else                G.omega[j] = yd[r] + xd[r];
            }
            G.h_edge[n] = zt[0]; G.H = zt[0];
            rt_ipss_finalize(&G);
            double ks = 1, i1s, i1p; int gd;
            rt_ipss_kappa(&G, vza, raa, &ks, &i1s, &i1p, &gd);
            kap[mode] = ks;
            rt_ipss_free(&G);
        }
        printf("%.1f,%.3f,%.5f,%.5f,%.9f,%.9f,%.9f\n", sza, acos(cth) * 180 / M_PI, PR, PA, kap[0], kap[1], kap[2]);
    }
    free(bl);
    return 0;
}
