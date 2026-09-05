/* ============================================================================
 * tests/test_ipss_gates.c — IPSS gates G-I0 / G-I1 / G-I2 / G-I3 / G-I4
 *   G-I0  geometry vs the Python oracle (tau_view, tau_sun, shadow flag)
 *   G-I1  single-scatter PSS-vs-IPSS anchors (reference package T9 table)
 *   G-I2  plane-parallel limit (Re x 1e4 -> kappa == 1) and nadir behaviour
 *
 * Build:  gcc -std=c11 -O2 -Isrc tests/test_ipss_gates.c src/rt_ipss.c \
 *              src/rt_atm.c ... -lm      (see scripts/run_ipss_gates.sh)
 * The test builds its own rt_ipss_t directly so it needs no solver context.
 * Output: one line per check, "GATE:" summary, exit 0 only when all pass.
 * ========================================================================== */
#include "rt_ipss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_fail = 0;

/* Build the paper's homogeneous validation atmosphere directly:
 * 0-100 km, n uniform layers, conservative scattering, total tau = tau_t. */
static void build_homog(rt_ipss_t *G, double tau_t, double sza_deg, int nlay,
                        double Re)
{
    memset(G, 0, sizeof(*G));
    G->nlay   = nlay;
    G->Re     = Re;
    G->H      = 100.0;
    G->mu_sun = cos(sza_deg * M_PI / 180.0);
    G->n_quad = 2;    /* nodes PER LAYER PANEL (segment-wise rule) */
    G->h_edge = (double *)malloc((size_t)(nlay + 1) * sizeof(double));
    G->beta   = (double *)malloc((size_t)nlay * sizeof(double));
    G->omega  = (double *)malloc((size_t)nlay * sizeof(double));
    for (int k = 0; k <= nlay; ++k) G->h_edge[k] = 100.0 * (double)k / nlay;
    for (int j = 0; j < nlay; ++j) { G->beta[j] = tau_t / 100.0; G->omega[j] = 1.0; }
    rt_ipss_finalize(G);
}

/* Exponential (scale height 8 km) profile with vertical tau = tau_t; used for
 * the Eq.(5c) radius check where the literal and exact forms differ. */
static void build_exp(rt_ipss_t *G, double tau_t, double sza_deg, int nlay)
{
    build_homog(G, tau_t, sza_deg, nlay, RT_IPSS_EARTH_RADIUS_KM);
    double sum = 0.0;
    for (int j = 0; j < nlay; ++j) {
        const double hm = 0.5 * (G->h_edge[j] + G->h_edge[j + 1]);
        G->beta[j] = exp(-hm / 8.0);
        sum += G->beta[j] * (G->h_edge[j + 1] - G->h_edge[j]);
    }
    for (int j = 0; j < nlay; ++j) G->beta[j] *= tau_t / sum;
    rt_ipss_finalize(G);
}

/* Exponential profile with an arbitrary scale height; used by G-I3 to collapse
 * the atmosphere onto the surface. */
static void build_exp_hs(rt_ipss_t *G, double tau_t, double sza_deg, int nlay,
                         double Hs)
{
    build_homog(G, tau_t, sza_deg, nlay, RT_IPSS_EARTH_RADIUS_KM);
    double sum = 0.0;
    for (int j = 0; j < nlay; ++j) {
        const double hm = 0.5 * (G->h_edge[j] + G->h_edge[j + 1]);
        G->beta[j] = exp(-hm / Hs);
        sum += G->beta[j] * (G->h_edge[j + 1] - G->h_edge[j]);
    }
    for (int j = 0; j < nlay; ++j) G->beta[j] *= tau_t / sum;
    rt_ipss_finalize(G);
}

int main(void)
{
    rt_ipss_t G;

    /* ---- G-I0a/G-I0c: geometry vs the Python oracle fixture ------------
     * tests/fixtures/ipss_oracle_geometry.csv is emitted by the reference
     * implementation (ipss_geometry_reference.py) at full double precision:
     *   tau_view,vza,0,0,0,value        line-of-sight slant optical depth
     *   tau_sun,rBx,rBy,rBz,sza,value   solar slant depth (-1 = shadowed)
     * Gate: relative difference < 1e-9 and identical shadow classification. */
    {
        build_homog(&G, 0.25, 0.0, 100, RT_IPSS_EARTH_RADIUS_KM);
        FILE *fp = fopen("tests/fixtures/ipss_oracle_geometry.csv", "r");
        if (!fp) { printf("[FAIL] G-I0a oracle fixture missing\n"); ++g_fail; }
        else {
            char line[512];
            double worst_v = 0.0, worst_s = 0.0;
            int n_v = 0, n_s = 0, n_shadow = 0, shadow_mismatch = 0;
            if (!fgets(line, sizeof line, fp)) line[0] = 0;   /* header */
            while (fgets(line, sizeof line, fp)) {
                char kind[32];
                double a, b, c, d, val;
                if (sscanf(line, "%31[^,],%lf,%lf,%lf,%lf,%lf",
                           kind, &a, &b, &c, &d, &val) != 6) continue;
                if (!strcmp(kind, "tau_view")) {
                    const double t = rt_ipss_tau_view_toa(&G, a);
                    const double rel = fabs(t - val) / (val > 0 ? val : 1.0);
                    if (rel > worst_v) worst_v = rel;
                    ++n_v;
                } else if (!strcmp(kind, "tau_sun")) {
                    G.mu_sun = cos(d * M_PI / 180.0);
                    const double rB[3] = { a, b, c };
                    const double t = rt_ipss_tau_sun_at(&G, rB);
                    const int shadow_ref = (val < 0.0);
                    const int shadow_got = !isfinite(t);
                    if (shadow_ref) ++n_shadow;
                    if (shadow_ref != shadow_got) ++shadow_mismatch;
                    else if (!shadow_ref) {
                        const double rel = fabs(t - val) / (val > 1e-12 ? val : 1.0);
                        if (rel > worst_s) worst_s = rel;
                    }
                    ++n_s;
                }
            }
            fclose(fp);
            printf("       tau_view n=%d worst rel %.2e | tau_sun n=%d (shadowed %d) worst rel %.2e\n",
                   n_v, worst_v, n_s, n_shadow, worst_s);
            printf("[%s] G-I0a geometry vs oracle (gate rel 1e-9)\n",
                   (worst_v < 1e-9 && worst_s < 1e-9) ? "PASS" : "FAIL");
            if (!(worst_v < 1e-9 && worst_s < 1e-9)) ++g_fail;
            printf("[%s] G-I0c shadow classification identical (%d mismatch of %d, %d shadowed)\n",
                   shadow_mismatch == 0 ? "PASS" : "FAIL", shadow_mismatch, n_s, n_shadow);
            if (shadow_mismatch != 0) ++g_fail;
        }
        rt_ipss_free(&G);
    }

    /* ---- G-I0b: Eq.(5c) radius — exact form must differ ~1.3% from the
     * literal Re form on a non-uniform profile.  Here we only assert that the
     * exact (ray-traced) tau_view differs from the literal closed form by the
     * documented magnitude, i.e. that we are NOT using the literal form. ---- */
    {
        build_exp(&G, 0.25, 0.0, 100);
        const double tv = 40.0;
        const double t_exact = rt_ipss_tau_view_toa(&G, tv);
        /* literal (5c): d_i = sin(theta_c)/sin(theta_v) * Re */
        const double tvr = tv * M_PI / 180.0;
        double t_lit = 0.0;
        double d_prev = 0.0;
        for (int k = G.nlay; k >= 0; --k) {
            const double s = (G.Re + G.H) / (G.Re + G.h_edge[k]) * sin(tvr);
            const double th = asin(s > 1.0 ? 1.0 : s);
            const double d = sin(th - tvr) / sin(tvr) * G.Re;
            if (k < G.nlay) t_lit += (d - d_prev) * G.beta[k];
            d_prev = d;
        }
        const double rel = fabs(t_lit - t_exact) / t_exact * 100.0;
        printf("       exp-atm tau_view exact %.6f  literal-(5c) %.6f  diff %.4f%%\n",
               t_exact, t_lit, rel);
        printf("[%s] G-I0b Eq.(5c) exact radius in use (diff vs literal %.2f%%, expect ~1.3%%)\n",
               (rel > 0.8 && rel < 2.0) ? "PASS" : "FAIL", rel);
        if (!(rel > 0.8 && rel < 2.0)) ++g_fail;
        rt_ipss_free(&G);
    }

    /* ---- G-I0d: vertical-sun limit, tau_sun on the axis = residual tau ---- */
    {
        build_homog(&G, 0.25, 0.0, 100, RT_IPSS_EARTH_RADIUS_KM);
        const double rB[3] = { 0.0, 0.0, RT_IPSS_EARTH_RADIUS_KM + 30.0 };
        const double tau = rt_ipss_tau_sun_at(&G, rB);
        const double want = 0.25 * 70.0 / 100.0;
        printf("       tau_sun(30 km, SZA 0) = %.9f  expect %.9f\n", tau, want);
        printf("[%s] G-I0d vertical sun limit\n", fabs(tau - want) < 1e-9 ? "PASS" : "FAIL");
        if (!(fabs(tau - want) < 1e-9)) ++g_fail;
        rt_ipss_free(&G);
    }

    /* ---- G-I1: reference-package T9 anchors, PSS vs IPSS single scatter ---
     * tau_t = 0.25, theta_s = 84.26 deg, paper phi_v = 0 (OCRT raa = 180). --- */
    {
        build_homog(&G, 0.25, 84.26, 100, RT_IPSS_EARTH_RADIUS_KM);
        const double vza[5] = { 0.0, 20.0, 40.0, 60.0, 70.0 };
        const double ref[5] = { 0.000, -1.151, -2.572, -4.970, -7.365 };
        double worst = 0.0;
        printf("       vza    (I1_pss/I1_ipss-1)%%   oracle%%    diff pp\n");
        for (int i = 0; i < 5; ++i) {
            const double ip = rt_ipss_single_scatter(&G, vza[i], 180.0,
                                                     RT_IPSS_SUN_SPHERICAL);
            const double ps = rt_ipss_single_scatter(&G, vza[i], 180.0,
                                                     RT_IPSS_SUN_NADIR);
            const double d = (ps / ip - 1.0) * 100.0;
            const double diff = fabs(d - ref[i]);
            if (diff > worst) worst = diff;
            printf("       %4.0f   %18.3f   %8.3f   %6.3f\n", vza[i], d, ref[i], diff);
        }
        printf("[%s] G-I1 T9 anchors (worst %.4f pp, gate 0.1 pp)\n",
               worst < 0.1 ? "PASS" : "FAIL", worst);
        if (!(worst < 0.1)) ++g_fail;

        /* quadrature convergence: refining the per-panel rule must not move it */
        const double a = rt_ipss_single_scatter(&G, 70.0, 180.0, RT_IPSS_SUN_SPHERICAL);
        G.n_quad = 8;
        const double b = rt_ipss_single_scatter(&G, 70.0, 180.0, RT_IPSS_SUN_SPHERICAL);
        const double relq = fabs(b / a - 1.0);
        printf("[%s] G-I1b per-panel quadrature 2->8 nodes rel %.2e (gate 1e-8)\n",
               relq < 1e-8 ? "PASS" : "FAIL", relq);
        if (!(relq < 1e-8)) ++g_fail;
        rt_ipss_free(&G);
    }

    /* ---- G-I2a: plane-parallel limit Re x 1e4 -> kappa == 1 -------------- */
    {
        build_homog(&G, 0.25, 70.47, 100, RT_IPSS_EARTH_RADIUS_KM * 1.0e4);
        double kap = 0.0; int guard = 0;
        double worst = 0.0;
        for (double vza = 0.0; vza <= 70.0; vza += 10.0) {
            rt_ipss_kappa(&G, vza, 180.0, &kap, NULL, NULL, &guard);
            const double d = fabs(kap - 1.0);
            if (d > worst) worst = d;
        }
        printf("[%s] G-I2a plane-parallel limit |kappa-1| max %.2e (gate 1e-4)\n",
               worst < 1e-4 ? "PASS" : "FAIL", worst);
        if (!(worst < 1e-4)) ++g_fail;
        rt_ipss_free(&G);
    }

    /* ---- G-I2b: azimuth independence at nadir, azimuth dependence at VZA -- */
    {
        build_homog(&G, 0.25, 84.26, 100, RT_IPSS_EARTH_RADIUS_KM);
        double k0 = 0, k180 = 0, k70a = 0, k70b = 0;
        rt_ipss_kappa(&G, 0.0,  0.0,   &k0,   NULL, NULL, NULL);
        rt_ipss_kappa(&G, 0.0,  180.0, &k180, NULL, NULL, NULL);
        rt_ipss_kappa(&G, 70.0, 0.0,   &k70a, NULL, NULL, NULL);
        rt_ipss_kappa(&G, 70.0, 180.0, &k70b, NULL, NULL, NULL);
        const double dn = fabs(k0 / k180 - 1.0);
        const double dv = fabs(k70a / k70b - 1.0);
        printf("       kappa nadir %.6f/%.6f ; vza70 raa0 %.6f raa180 %.6f\n",
               k0, k180, k70a, k70b);
        printf("[%s] G-I2b nadir azimuth-invariant (%.1e) and slant azimuth-dependent (%.3f)\n",
               (dn < 1e-12 && dv > 1e-3) ? "PASS" : "FAIL", dn, dv);
        if (!(dn < 1e-12 && dv > 1e-3)) ++g_fail;
        rt_ipss_free(&G);
    }

    /* ---- G-I2c: guard on vanishing I_1,pp -------------------------------- */
    {
        build_homog(&G, 0.0, 40.0, 10, RT_IPSS_EARTH_RADIUS_KM);   /* tau = 0 */
        double kap = 0.0; int guard = 0;
        rt_ipss_kappa(&G, 40.0, 90.0, &kap, NULL, NULL, &guard);
        printf("[%s] G-I2c I1_pp guard (kappa=%.3f, guard=%d)\n",
               (guard == 1 && kap == 1.0) ? "PASS" : "FAIL", kap, guard);
        if (!(guard == 1 && kap == 1.0)) ++g_fail;
        rt_ipss_free(&G);
    }

    /* ---- G-I3: VZA anchor.  Collapsing the atmosphere onto the surface must
     * drive the correction to zero, because all the scattering then happens at
     * the level where the spherical ray and mu_v share the same zenith angle.
     * Only the SURFACE anchor satisfies this; the TOA anchor keeps a ~3 %
     * offset forever, which is how the anchor bug was found (see rt_ipss.h). */
    {
        const double hs[4] = { 8.0, 2.0, 0.5, 0.1 };
        double dev_surf[4], dev_toa[4];
        for (int i = 0; i < 4; ++i) {
            build_exp_hs(&G, 0.0935, 0.0, 400, hs[i]);
            double kap = 1.0;
            rt_ipss_kappa(&G, 55.0, 90.0, &kap, NULL, NULL, NULL);
            dev_surf[i] = 100.0 * (kap - 1.0);
            const double a_toa = rt_ipss_single_scatter(&G, 55.0, 90.0,
                                                        RT_IPSS_SUN_SPHERICAL);
            const double b_pp  = rt_ipss_single_scatter_pp(&G, 55.0);
            dev_toa[i] = 100.0 * (a_toa / b_pp - 1.0);
            rt_ipss_free(&G);
        }
        printf("       Hs[km]     8.0     2.0     0.5     0.1  (VZA 55, RAA 90, SZA 0)\n");
        printf("       surface %+7.3f %+7.3f %+7.3f %+7.3f  %% <- must -> 0\n",
               dev_surf[0], dev_surf[1], dev_surf[2], dev_surf[3]);
        printf("       TOA     %+7.3f %+7.3f %+7.3f %+7.3f  %% <- does not\n",
               dev_toa[0], dev_toa[1], dev_toa[2], dev_toa[3]);
        const int mono = fabs(dev_surf[3]) < fabs(dev_surf[2])
                      && fabs(dev_surf[2]) < fabs(dev_surf[1])
                      && fabs(dev_surf[1]) < fabs(dev_surf[0]);
        const int small = fabs(dev_surf[3]) < 0.02 && fabs(dev_surf[0]) < 0.5;
        const int toa_bad = fabs(dev_toa[3]) > 2.0;   /* documents the contrast */
        printf("[%s] G-I3 surface-collapse limit (kappa->1 only with the pixel anchor)\n",
               (mono && small && toa_bad) ? "PASS" : "FAIL");
        if (!(mono && small && toa_bad)) ++g_fail;
    }

    /* ---- G-I4: phase-weighted source (v1.11.1).  Eq. (7) divides the REAL
     * single-scattering radiances, so with two scatterers whose mixing ratio
     * changes with height the phase functions do not cancel.
     *   G-I4a  height-independent mixture -> kappa identical to phase-free
     *          (exact cancellation, 1e-12).
     *   G-I4b  Rayleigh (8 km, tau 0.0935, depol 0.0279) over aerosol (2 km,
     *          tau 0.1, omega 0.95, HG g = 0.7): kappa must match the
     *          independent scipy quadrature (indep_two_component_quad.py,
     *          continuous profiles, no shared code) and must differ from the
     *          phase-free value where that reference says so. */
    {
        const double HS_R = 8.0, TR = 0.0935, HS_A = 2.0, TA = 0.1, OMA = 0.95, GHG = 0.7;
        const int NL = 400, LHG = 120;
        /* fixture from indep_two_component_quad.py (2026-09-05) */
        const double fx[7][5] = {   /* vza, sza, raa, kappa_phase, kappa_free */
            { 55,  0,  90, 0.997782012, 0.998564708 },
            { 55, 40,  90, 0.997864991, 0.998606864 },
            { 55, 80,  90, 1.007490360, 1.009912232 },
            { 70, 80,   0, 1.001467078, 1.007529564 },
            { 70, 80, 180, 1.006466295, 1.001975469 },
            { 40, 60,  45, 0.999593314, 0.999985608 },
            {  0, 40,   0, 1.000058665, 1.000075029 } };
        double worst_a = 0.0, worst_b = 0.0, worst_free = 0.0, min_sep = 1e9;
        for (int c = 0; c < 7; ++c) {
            const double vza = fx[c][0], sza = fx[c][1], raa = fx[c][2];
            for (int mode = 0; mode < 3; ++mode) {   /* 0 uniform-mix phase, 1 uniform-mix free,
                                                      * 2 separated phase */
                build_homog(&G, 1.0, sza, NL, RT_IPSS_EARTH_RADIUS_KM);
                const double nR = TR / (HS_R * (1.0 - exp(-100.0 / HS_R)));
                const double nA = TA / (HS_A * (1.0 - exp(-100.0 / HS_A)));
                G.w_ray = (double *)malloc((size_t)NL * sizeof(double));
                G.w_aer = (double *)malloc((size_t)NL * sizeof(double));
                for (int j = 0; j < NL; ++j) {
                    const double hm = 0.5 * (G.h_edge[j] + G.h_edge[j + 1]);
                    double bR = nR * exp(-hm / HS_R), bA = nA * exp(-hm / HS_A);
                    if (mode < 2) { bA = bR; }           /* same profile: ratio fixed */
                    G.beta[j]  = bR + bA;
                    G.w_ray[j] = bR / G.beta[j];
                    G.w_aer[j] = OMA * bA / G.beta[j];
                    G.omega[j] = G.w_ray[j] + G.w_aer[j];
                }
                if (mode == 1) { free(G.w_ray); free(G.w_aer); G.w_ray = G.w_aer = NULL; }
                G.depol = 0.0279; G.L_aer = LHG;
                G.betal_aer = (double *)malloc((size_t)(LHG + 1) * sizeof(double));
                for (int l = 0; l <= LHG; ++l) G.betal_aer[l] = (2.0 * l + 1.0) * pow(GHG, l);
                rt_ipss_finalize(&G);
                double kap = 1.0;
                rt_ipss_kappa(&G, vza, raa, &kap, NULL, NULL, NULL);
                static double k_uni_phase = 0.0;
                if (mode == 0) k_uni_phase = kap;
                else if (mode == 1) {
                    const double d = fabs(kap / k_uni_phase - 1.0);
                    if (d > worst_a) worst_a = d;
                } else {
                    const double d = fabs(kap - fx[c][3]);
                    if (d > worst_b) worst_b = d;
                    const double sep = fabs(fx[c][3] - fx[c][4]);
                    if (sep < min_sep) min_sep = sep;
                    /* phase-free value of the same separated profile, for the table */
                    double kf = 1.0;
                    free(G.w_src); G.w_src = NULL;
                    rt_ipss_kappa(&G, vza, raa, &kf, NULL, NULL, NULL);
                    const double df = fabs(kf - fx[c][4]);
                    if (df > worst_free) worst_free = df;
                    printf("       vza %3.0f sza %2.0f raa %3.0f  kappa phase %.6f (py %.6f)  free %.6f (py %.6f)\n",
                           vza, sza, raa, kap, fx[c][3], kf, fx[c][4]);
                }
                rt_ipss_free(&G);
            }
        }
        const int ok_a = worst_a < 1e-12, ok_b = worst_b < 3e-5 && worst_free < 3e-5;
        printf("[%s] G-I4a uniform mixture: phase-weighted == phase-free (rel %.1e, gate 1e-12)\n",
               ok_a ? "PASS" : "FAIL", worst_a);
        printf("[%s] G-I4b separated mixture vs independent quadrature: |dkappa| phase %.1e, free %.1e (gate 3e-5; definitions differ by >= %.1e)\n",
               ok_b ? "PASS" : "FAIL", worst_b, worst_free, min_sep);
        if (!ok_a || !ok_b) ++g_fail;
    }

    printf("\nIPSS GATES: %s\n", g_fail == 0 ? "ALL PASS" : "FAILURES PRESENT");
    return g_fail == 0 ? 0 : 1;
}
