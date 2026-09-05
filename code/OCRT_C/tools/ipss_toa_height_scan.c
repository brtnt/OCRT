/* 모델 TOA 높이 H 의존성: 지상 앵커 vs TOA 앵커.
 * 물리적 답은 H(모델 절단 고도)에 둔감해야 한다 — 30 km 위에는 소광이 없다. */
#include "rt_ipss.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
static void build(rt_ipss_t *G, double tau_t, double sza, int nlay, double Hs, double H)
{
    memset(G, 0, sizeof(*G));
    G->nlay = nlay; G->Re = RT_IPSS_EARTH_RADIUS_KM; G->H = H;
    G->mu_sun = cos(sza * M_PI / 180.0); G->n_quad = 64;
    G->h_edge = malloc((size_t)(nlay + 1) * sizeof(double));
    G->beta = malloc((size_t)nlay * sizeof(double));
    G->omega = malloc((size_t)nlay * sizeof(double));
    for (int k = 0; k <= nlay; ++k) G->h_edge[k] = H * (double)k / nlay;
    double sum = 0.0;
    for (int j = 0; j < nlay; ++j) {
        const double hm = 0.5 * (G->h_edge[j] + G->h_edge[j + 1]);
        G->beta[j] = (Hs > 0.0) ? exp(-hm / Hs) : 1.0;
        sum += G->beta[j] * (G->h_edge[j + 1] - G->h_edge[j]);
        G->omega[j] = 1.0;
    }
    /* 연직 tau 를 H 와 무관하게 고정: 실제 대기를 더 높이까지 잘라도 질량은 같다 */
    for (int j = 0; j < nlay; ++j) G->beta[j] *= tau_t / sum;
    rt_ipss_finalize(G);
}
int main(int argc, char **argv)
{
    const double tau = (argc > 1) ? atof(argv[1]) : 0.0935;
    const double sza = (argc > 2) ? atof(argv[2]) : 0.0;
    const double Hs  = (argc > 3) ? atof(argv[3]) : 8.0;
    const double vza = (argc > 4) ? atof(argv[4]) : 55.0;
    const double Re = RT_IPSS_EARTH_RADIUS_KM;
    printf("H_km,dtheta_deg,dev_surf_pct,dev_toa_pct\n");
    for (double H = 20.0; H <= 200.001; H += 5.0) {
        rt_ipss_t G; build(&G, tau, sza, 400, Hs, H);
        const double dth = asin(Re / (Re + H) * sin(vza * M_PI / 180.0)) * 180.0 / M_PI - vza;
        const double b_pp = rt_ipss_single_scatter_pp(&G, vza);
        const double a_s = rt_ipss_single_scatter_surf(&G, vza, 90.0, RT_IPSS_SUN_SPHERICAL);
        const double a_t = rt_ipss_single_scatter(&G, vza, 90.0, RT_IPSS_SUN_SPHERICAL);
        printf("%.1f,%.6f,%.6f,%.6f\n", H, dth,
               100.0 * (a_s / b_pp - 1.0), 100.0 * (a_t / b_pp - 1.0));
        rt_ipss_free(&G);
    }
    return 0;
}
