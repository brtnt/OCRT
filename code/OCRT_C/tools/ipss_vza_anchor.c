/* ============================================================================
 * tools/ipss_vza_anchor.c
 * VZA 앵커 규약 민감도 진단.
 *
 * 평면평행 해에는 "이 μ_v 가 어느 고도의 천정각인가"라는 개념이 없다.
 * 구면화하는 순간 그 앵커를 정해야 하고, Eq.(7) 의 κ = I1_sph/I1_pp 는
 * 앵커에 따라 값과 부호가 모두 달라진다:
 *
 *   TOA 앵커 :  κ_toa(θ)   = I1_sph(θ_TOA = θ) / I1_pp(θ)
 *                구면 광선의 국소 천정각은 아래로 갈수록 커지므로 경로가
 *                길어지고 κ > 1.
 *   지상 앵커 :  κ_gnd(θ_g) = I1_sph(θ_TOA = asin(Re/(Re+H) sinθ_g)) / I1_pp(θ_g)
 *                같은 지상 천정각을 평면평행이 전 고도에 유지하므로 평면평행이
 *                경로를 과대평가하고 κ < 1.
 *
 * 위성 L1B 의 sensor zenith 는 통상 화소(지상) 기준이다.
 * Usage: ipss_vza_anchor [tau_total] [sza_deg] [scale_height_km | 0=uniform]
 * ========================================================================== */
#include "rt_ipss.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void build(rt_ipss_t *G, double tau_t, double sza, int nlay, double Hs)
{
    memset(G, 0, sizeof(*G));
    G->nlay = nlay; G->Re = RT_IPSS_EARTH_RADIUS_KM; G->H = 100.0;
    G->mu_sun = cos(sza * M_PI / 180.0); G->n_quad = 64;
    G->h_edge = (double *)malloc((size_t)(nlay + 1) * sizeof(double));
    G->beta   = (double *)malloc((size_t)nlay * sizeof(double));
    G->omega  = (double *)malloc((size_t)nlay * sizeof(double));
    for (int k = 0; k <= nlay; ++k) G->h_edge[k] = 100.0 * (double)k / nlay;
    double sum = 0.0;
    for (int j = 0; j < nlay; ++j) {
        const double hm = 0.5 * (G->h_edge[j] + G->h_edge[j + 1]);
        G->beta[j] = (Hs > 0.0) ? exp(-hm / Hs) : 1.0;
        sum += G->beta[j] * (G->h_edge[j + 1] - G->h_edge[j]);
        G->omega[j] = 1.0;
    }
    for (int j = 0; j < nlay; ++j) G->beta[j] *= tau_t / sum;
    rt_ipss_finalize(G);
}

int main(int argc, char **argv)
{
    const double tau = (argc > 1) ? atof(argv[1]) : 0.0935;
    const double sza = (argc > 2) ? atof(argv[2]) : 0.0;
    const double Hs  = (argc > 3) ? atof(argv[3]) : 8.0;
    const double Re  = RT_IPSS_EARTH_RADIUS_KM, Hatm = 100.0;
    rt_ipss_t G;
    build(&G, tau, sza, 400, Hs);

    printf("tau_total=%.4f  SZA=%.2f  profile=%s  RAA=90\n", tau, sza,
           (Hs > 0.0) ? "exponential" : "uniform");
    printf("%8s %10s %12s %12s | %10s %12s %12s\n",
           "VZA", "th_gnd", "kappa_TOA", "dev_TOA[%]",
           "th_TOA", "kappa_srf", "dev_srf[%]");
    for (double v = 0.0; v <= 80.001; v += 5.0) {
        /* --- TOA 앵커: 입력 VZA 가 TOA 천정각 (논문 §3 규약) */
        const double a_toa = rt_ipss_single_scatter(&G, v, 90.0, RT_IPSS_SUN_SPHERICAL);
        const double b_pp  = rt_ipss_single_scatter_pp(&G, v);
        const double k_toa = (b_pp > 0.0 && a_toa > 0.0) ? a_toa / b_pp : 1.0;
        const double sg = fmin(1.0, (Re + Hatm) / Re * sin(v * M_PI / 180.0));
        const double th_g_of_toa = asin(sg) * 180.0 / M_PI;
        /* --- 지상(화소) 앵커: 입력 VZA 가 지상 천정각 = L1B senz (생산 기본) */
        const double a_srf = rt_ipss_single_scatter_surf(&G, v, 90.0,
                                                         RT_IPSS_SUN_SPHERICAL);
        const double k_srf = (b_pp > 0.0 && a_srf > 0.0) ? a_srf / b_pp : 1.0;
        const double th_toa_of_g = asin(Re / (Re + Hatm)
                                        * sin(v * M_PI / 180.0)) * 180.0 / M_PI;
        printf("%8.2f %10.3f %12.6f %+12.3f | %10.3f %12.6f %+12.3f\n",
               v, th_g_of_toa, k_toa, 100.0 * (1.0 / k_toa - 1.0),
               th_toa_of_g, k_srf, 100.0 * (1.0 / k_srf - 1.0));
    }
    printf("dev = (plane-parallel / IPSS - 1) [%%] = 1/kappa - 1\n");
    rt_ipss_free(&G);
    return 0;
}
