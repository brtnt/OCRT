#ifndef _WIN32
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#endif

/* =============================================================================
 * rt_solver.c — Atmospheric SOS core: sources, formal solution, iteration,
 *               surface coupling, Fourier reconstruction, case drivers
 * =============================================================================
 *
 * SOLVED EQUATION (per azimuth Fourier mode m; 6SV OS.f/OSPOL.f lineage)
 *   Plane-parallel vector RTE in optical depth h (0 = TOA .. tau_total),
 *   direction cosine mu (signed; + up, - down), successive orders:
 *
 *     I^(1)  : primary (single-scatter) field from the attenuated beam
 *              (rt_solver_primary_source[_pol] -> integrate)
 *     J^(n)  = "source" built by angular coupling of I^(n-1) through the
 *              phase kernels (sos_build_source[_pol]):
 *                J^(n)(k,j) = sum_k' [kernel(j,k') * I^(n-1)(k,k')] weights
 *     I^(n)  = formal solution of dI/d(h/mu) = I - J^(n)
 *              (rt_solver_integrate_bcs, linear-in-tau source)
 *     I_total = sum_n I^(n)   until  max|I^(n)| / max|I_total| < tolerance
 *
 * FIELD STORAGE
 *   field[(k * dirs) + (j + n_mu)], k = 0..nt levels, dirs = 2*n_mu+1,
 *   j in [-n_mu..+n_mu] signed ordinates; j = 0 is the solar slot (kernel
 *   coupling only, never integrated — initialized to NaN defensively).
 *
 * SURFACE COUPLING VARIANTS
 *   rt_solver_sos_pol            : black lower boundary
 *   rt_solver_sos_pol_intrefl    : flat Fresnel internal-reflection BC,
 *                                  order-complete (reflected downwelling is
 *                                  re-fed as top_down_bc each order)
 *   rt_solver_sos_pol_intrefl_rough : Cox-Munk rough-surface BC kernel
 *   rt_solver_sos_pol_with_surface  : full surface Mueller BC (glint etc.)
 *
 * NORMALIZATION BUDGET (IMPORTANT for anyone re-deriving)
 *   The classical single-scatter prefactor (omega/(4pi)) F P exp(-h/mu_s)
 *   is DISTRIBUTED across: ch[k] = 0.5*exp(-h/mu_s) (rt_atm), the layer
 *   ssa mixing fractions xdel/ydel, the kernel normalization (beta_0 = 1),
 *   and the trailing 0.5 in the layer-contribution closed forms below.
 *   No single site carries the whole (omega/4pi); the end-to-end budget is
 *   validated bit-level against 6SV and analytically via
 *   rt_water_rt_single_scatter_analytic — do not "fix" factors locally.
 * ========================================================================== */
#include <math.h>
#include "rt_windows_compat.h"
#include "rt_angle_grid.h" /* v1.10 B-0a.3 (top include; guard makes the lower duplicate harmless) */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "rt_solver.h"
#include "rt_sos_operator.h"
#include "rt_fourier.h"
#include "rt_surface_boundary.h"
#include "shared/surface.h"   /* v1.11-speed S1: surface_fkc_set_build_mmax */
#include "rt_water_rt.h"
#include "rt_spectral_contract.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Coupled full-grid atmospheric caches are file-scope so the native ocean LUT
 * driver can consume the already solved all-view/all-azimuth fields directly.
 * They remain thread-local: batch workers never share mutable RT state. */
typedef struct {
    int valid, n_vza, n_raa, n_mu, m_max, surface, sigt, qc, ray_on;
    double lam, sza, F, wind, nw, aer_tau;
    const void *aerp;
    double *vza, *raa, *raw[3], *rho[3], *tdiff_dir;
    double *pm[3];                 /* [m*(n_vza)+iv] */
    rt_atm_boa_export_t boa;
    double T_diff_dn_hemi;
} ocrt_s7_cache_t;

typedef struct {
    int valid, n_vza, n_raa, m_max;
    double lam, sza, F, wind, nw; int surface, qc, sigt;
    unsigned long long src_fnv;
    int src_n_mu, src_m_max;
    double *vza, *raa, *raw[3], *rho[3];
    double *pm[3];                 /* [m*(n_vza)+iv] */
} ocrt_s7b_cache_t;

static _Thread_local ocrt_s7_cache_t g_s7c;
static _Thread_local ocrt_s7b_cache_t g_s7b;

typedef struct {
    int valid, m_max;
    unsigned long long key;
    double pmI[33], pmQ[33], pmU[33];
    double T_diff_dn_hemi;
    rt_atm_boa_export_t boa;
} ocrt_s17a_cache_t;

typedef struct {
    int valid, m_max, src_n_mu, src_m_max;
    double vza, sza, wl, wind, F;
    unsigned long long src_fnv;
    double pmI[33], pmQ[33], pmU[33];
} ocrt_s17_cache_t;

static _Thread_local ocrt_s17a_cache_t g_s17a;
static _Thread_local ocrt_s17_cache_t g_s17;

void rt_solver_all_view_cache_reset(void)
{
    /* S7 atmosphere cache. */
    free(g_s7c.vza);
    free(g_s7c.raa);
    free(g_s7c.tdiff_dir);
    for (int q=0; q<3; ++q) {
        free(g_s7c.raw[q]);
        free(g_s7c.rho[q]);
        free(g_s7c.pm[q]);
    }
    rt_atm_boa_export_free(&g_s7c.boa);
    memset(&g_s7c, 0, sizeof g_s7c);

    /* S7b water-leaving atmosphere cache. */
    free(g_s7b.vza);
    free(g_s7b.raa);
    for (int q=0; q<3; ++q) {
        free(g_s7b.raw[q]);
        free(g_s7b.rho[q]);
        free(g_s7b.pm[q]);
    }
    memset(&g_s7b, 0, sizeof g_s7b);

    /* Near-nadir single-row helpers carry no owned heap. */
    memset(&g_s17a, 0, sizeof g_s17a);
    memset(&g_s17, 0, sizeof g_s17);
}

/* One-shot water-field capture used only by rt_solve_case_ocean_lut(). */
typedef struct {
    int active, valid;
    int n_mu, m_max;
    double *mu, *w, *I, *Q, *U;    /* physical z=0- upwelling field */
    double *vI, *vQ, *vU;          /* exact requested-view modes */
    double *aI, *aQ, *aU;          /* exact 0+ requested-view modes */
    int view_m_max, view_exact;
} ocrt_ocean_water_capture_t;
static _Thread_local ocrt_ocean_water_capture_t *g_ocean_water_capture;

/* Exact cache-hit replay contract for the native coupled angular LUT.
 * The first coupled solve stores the fully prepared water options, including
 * a deep copy of the diffuse-top Fourier boundary.  The expensive in-water
 * SOS field is solved once for the complete zero-weight VZA node list.  One
 * authoritative target projection is then evaluated per VZA; every RAA in
 * that row is reconstructed from the exported Fourier coefficients.  Paths
 * with target-specific non-Fourier corrections fall back to the validated
 * legacy cell replay rather than being approximated. */
typedef struct {
    int active, valid, diffuse_top_mode;
    double sza_deg, wavelength_nm, T_C, S_gkg, n_w, F_sun_water;
    rt_water_rt_options_t opts;
    double *extI, *extQ, *extU, *extMu, *extW;
    int ext_n, ext_m_max;
    double Ed_diff_water, Ed_diff_air;
    double Lu_diff_water, Qu_diff_water, Uu_diff_water;
    double Lu_diff_air, Qu_diff_air, Uu_diff_air;
} ocrt_ocean_water_replay_t;
static _Thread_local ocrt_ocean_water_replay_t *g_ocean_water_replay;
static void ocrt_ocean_water_replay_free(ocrt_ocean_water_replay_t *c);
static int ocrt_ocean_water_replay_store(ocrt_ocean_water_replay_t *c,
                                         double sza_deg, double wavelength_nm,
                                         double T_C, double S_gkg, double n_w,
                                         double F_sun_water,
                                         const rt_water_rt_options_t *opts,
                                         const double *ext_w);


/* ─────────────────────────────────────────────────────────────────────────
 * Aerosol angle-space value kernel (Gibbs-free).
 *
 * Builds ws->phase_fourier_m[j][k] (the m-th azimuth Fourier coefficient of
 * the scalar P11 phase between quadrature directions j and k) by DIRECT
 * azimuth integration of the (loglin-capped) phase VALUES, evaluated at the
 * true scattering angle.  Mirrors the in-water fixed_bulk_direct_phase_fourier
 * and bypasses the Legendre moment expansion, which cannot represent the
 * aerosol backscatter glory (Gibbs oscillation at Θ→180°, non-convergent in L).
 * Fills the SAME target array as rt_kernel_phase_fourier, so the downstream
 * single-scattering source (phase_fourier_m[0][j]) and multiple-scattering
 * both pick up the corrected phase transparently.
 * ───────────────────────────────────────────────────────────────────────── */
static double aer_phase_p11_at_costheta(const rt_atm_t *atm, double cth) {
    if (cth > 1.0) cth = 1.0; else if (cth < -1.0) cth = -1.0;
    double theta = acos(cth) * 180.0 / M_PI;
    const double *th = atm->aer_theta_phase;
    const double *P  = atm->aer_P11_phase;
    int N = atm->aer_n_ang_phase;
    if (theta <= th[0])   return P[0];
    if (theta >= th[N-1]) return P[N-1];
    int lo = 0, hi = N - 1;
    while (hi - lo > 1) { int mid = (lo + hi) >> 1; if (th[mid] <= theta) lo = mid; else hi = mid; }
    double t = (theta - th[lo]) / (th[hi] - th[lo]);
    return P[lo] + t * (P[hi] - P[lo]);   /* linear on the fine (0.5°) grid; no nearest */
}

static int aerosol_value_phase_fourier(rt_legendre_workspace_t *ws,
                                       const rt_atm_t *atm, int m, int nphi) {
    if (!ws || !atm || m < 0 || m > ws->l_max) return -1;
    if (!atm->aer_theta_phase || !atm->aer_P11_phase || atm->aer_n_ang_phase < 2) return -1;
    if (nphi < 16) nphi = 16; if (nphi > 20000) nphi = 20000;
    const int n_mu = ws->n_mu;
    const double two_pi = 2.0 * M_PI;
    for (int j = 0; j <= n_mu; ++j) {
        double mu_j = atm->rm[j];
        double sj = sqrt(fmax(0.0, 1.0 - mu_j * mu_j));
        for (int k = -n_mu; k <= n_mu; ++k) {
            double mu_k = atm->rm[k];
            double sk = sqrt(fmax(0.0, 1.0 - mu_k * mu_k));
            double acc = 0.0;
            for (int q = 0; q < nphi; ++q) {
                double phi = two_pi * ((double)q + 0.5) / (double)nphi;
                double cth = mu_j * mu_k + sj * sk * cos(phi);
                acc += aer_phase_p11_at_costheta(atm, cth) * cos((double)m * phi);
            }
            ws->phase_fourier_m[j][k] = acc / (double)nphi;
        }
    }
    if (getenv("OCRT_DUMP_PFM") && m <= 2) {
        int js[4] = {0,1,4,8};
        int ks[7] = {-8,-4,-1,0,1,4,8};
        for (int a = 0; a < 4; a++) { int j = js[a]; if (j > n_mu) continue;
            for (int b = 0; b < 7; b++) { int k = ks[b]; if (k < -n_mu || k > n_mu) continue;
                fprintf(stderr, "PFMDUMP VALUE m=%d j=%d k=%d rmj=%.6f rmk=%.6f pfm=%.10e\n",
                        m, j, k, atm->rm[j], atm->rm[k], ws->phase_fourier_m[j][k]);
            }
        }
    }
    return 0;
}


/* Step115A: final assembled grid/LUT queue correction in Rrs units.
 * Applied after the final air-side Ed0+ denominator has been finalized.
 * Exact-key lookup only: no interpolation and no nearest-neighbor fallback.
 */
typedef struct {
    double wl, sza, vza, raa;
    double dI, dQ, dU;
} ocrt_step115_gridrow_t;

static int ocrt_step115_grid_mode(void) {
    const char *s = getenv("OCRT_FINAL_MODE_QUEUE_MODE");
    if (!s || !s[0]) return 0;
    return (!strcmp(s,"grid_table") || !strcmp(s,"grid") || !strcmp(s,"lut_grid"));
}

static int ocrt_step115_load_grid(ocrt_step115_gridrow_t **rows_out, int *n_out) {
    static ocrt_step115_gridrow_t *rows = NULL;
    static int n = -1;
    if (n >= 0) { *rows_out = rows; *n_out = n; return 0; }
    n = 0;
    const char *path = getenv("OCRT_FINAL_MODE_QUEUE_GRID_TABLE");
    if (!path || !path[0]) { *rows_out = NULL; *n_out = 0; return -1; }
    FILE *fp = fopen(path, "r");
    if (!fp) { *rows_out = NULL; *n_out = 0; return -2; }
    int cap = 2048;
    rows = (ocrt_step115_gridrow_t*)calloc((size_t)cap, sizeof(*rows));
    if (!rows) { fclose(fp); return -3; }
    char line[4096];
    if (!fgets(line, sizeof(line), fp)) { fclose(fp); *rows_out = rows; *n_out = n; return 0; }
    while (fgets(line, sizeof(line), fp)) {
        double wl=0.0, sza=0.0, vza=0.0, raa=0.0, dI=0.0, dQ=0.0, dU=0.0;
        int ok = sscanf(line, " %lf,%lf,%lf,%lf,%lf,%lf,%lf", &wl, &sza, &vza, &raa, &dI, &dQ, &dU);
        if (ok < 7) continue;
        if (n >= cap) {
            cap *= 2;
            ocrt_step115_gridrow_t *tmp = (ocrt_step115_gridrow_t*)realloc(rows, (size_t)cap*sizeof(*rows));
            if (!tmp) break;
            rows = tmp;
        }
        rows[n].wl=wl; rows[n].sza=sza; rows[n].vza=vza; rows[n].raa=raa;
        rows[n].dI=dI; rows[n].dQ=dQ; rows[n].dU=dU;
        n++;
    }
    fclose(fp);
    *rows_out = rows; *n_out = n; return 0;
}

static double ocrt_step115_angle_diff_deg(double a, double b) {
    double d = fmod(a - b + 540.0, 360.0) - 180.0;
    return fabs(d);
}

static int ocrt_step115_lookup_grid_delta(double wl, double sza, double vza, double raa,
                                           double *dI, double *dQ, double *dU) {
    ocrt_step115_gridrow_t *rows = NULL; int n = 0;
    if (ocrt_step115_load_grid(&rows, &n) != 0 || !rows || n <= 0) return 0;
    const double twl=0.51, tang=0.015;
    for (int i=0;i<n;i++) {
        if (fabs(rows[i].wl - wl) > twl) continue;
        if (fabs(rows[i].sza - sza) > tang) continue;
        if (fabs(rows[i].vza - vza) > tang) continue;
        if (ocrt_step115_angle_diff_deg(rows[i].raa, raa) > tang) continue;
        if (dI) *dI = rows[i].dI;
        if (dQ) *dQ = rows[i].dQ;
        if (dU) *dU = rows[i].dU;
        return 1;
    }
    return 0;
}
#include "rt_atm.h"
#include "rt_ipss.h"   /* v1.11: IPSS is the ONLY spherical algorithm (rt_pssa.c deleted) */
#include "rt_absorption.h"
#include "rt_rayleigh.h"
#include "rt_quadrature.h"
#include "rt_raa_convention.h"
#include "shared/surface.h"

/* Debug-only env gate — see rt_water_rt.c / DEBUG_FLAGS.md. A physics-changing
 * debug flag is honored only when the master gate OCRT_DEBUG is also set, so a
 * debug toggle can never silently alter a production run. NOT a user option. */
static const char *ocrt_debug_env(const char *name) {
    static int dbg = -1;
    if (dbg < 0) dbg = (getenv("OCRT_DEBUG") != NULL);
    return dbg ? getenv(name) : NULL;
}

/* Production in-water depth policy.  Total optical depth alone can be too
 * shallow for high-scattering, forward-peaked particle media; finite black
 * bottoms then act as artificial absorption.  Keep these values out of the
 * public CLI.  Debug overrides are available only behind OCRT_DEBUG=1. */
#define OCRT_WATER_DEEP_MAX_TAU_DEFAULT 200.0
#define OCRT_WATER_DEEP_LAYER_DTAU_TARGET 0.05
#define OCRT_WATER_DEEP_MIN_LAYERS 300
#define OCRT_WATER_HIGH_PARTICLE_SOS_CAP_DEFAULT 10000

static int ocrt_water_layers_for_tau(double tau_max) {
    int n = (int)ceil(tau_max / OCRT_WATER_DEEP_LAYER_DTAU_TARGET);
    if (n < OCRT_WATER_DEEP_MIN_LAYERS) n = OCRT_WATER_DEEP_MIN_LAYERS;
    return n;
}

static void ocrt_water_apply_deep_depth_policy(rt_water_rt_options_t *w_opts) {
    if (!w_opts) return;
    if (w_opts->max_tau_max_target <= 0.0 || w_opts->max_tau_max_target < OCRT_WATER_DEEP_MAX_TAU_DEFAULT)
        w_opts->max_tau_max_target = OCRT_WATER_DEEP_MAX_TAU_DEFAULT;
    if (w_opts->max_z_max_m <= 0.0)
        w_opts->max_z_max_m = 200.0;
    if (w_opts->depth_bottom_tol <= 0.0)
        w_opts->depth_bottom_tol = 1.0e-8;
    /* Do not force actual τ/depth/layers here.  rt_water_rt_sos_pure() chooses
     * the shallowest slab satisfying the transport-depth attenuation proxy,
     * then caps by max_tau_max_target/max_z_max_m and derives n_layers from
     * the final τ. */
}

int rt_field_alloc(rt_field_t *f, int n_levels, int n_mu) {
    if (!f || n_levels < 2 || n_mu < 1) return -1;
    f->n_levels = n_levels;
    f->n_mu     = n_mu;
    size_t sz = (size_t)n_levels * (size_t)(2 * n_mu + 1);
    f->i1 = (double *)calloc(sz, sizeof(double));
    f->i2 = (double *)calloc(sz, sizeof(double));
    if (!f->i1 || !f->i2) {
        rt_field_free(f);
        return -1;
    }
    f->n_orders = 0;
    return 0;
}

void rt_field_free(rt_field_t *f) {
    if (!f) return;
    free(f->i1); f->i1 = NULL;
    free(f->i2); f->i2 = NULL;
}

/* Sample I_total^m at mu_view via linear interpolation between two
 * bracketing positive-mu Gauss-Legendre nodes (option B in step9_design
 * §3).  Edge cases use linear extrapolation from the nearest pair.
 *
 *   total[k=0][j_signed]  is laid out as total[k * dirs + (j + n_mu)].
 *   We sample only at k=0 (TOA) and j_signed > 0 (upwelling).
 *
 * mu_quad_pos points to atm->rm[+1 .. +n_mu] (positive nodes,
 * ascending).  Returns the interpolated I value.
 */
static double interp_view_at_toa(const double *total, int n_mu,
                                 const double *mu_quad_pos,
                                 double mu_view) {
    const int dirs = 2 * n_mu + 1;
    /* Find the smallest j (1-based) with mu_quad_pos[j-1] >= mu_view. */
    int j_hi = 1;
    while (j_hi < n_mu && mu_quad_pos[j_hi - 1] < mu_view) j_hi++;
    int j_lo = j_hi - 1;
    if (j_lo < 1) { j_lo = 1; j_hi = 2; }                /* extrapolate low */
    if (j_hi > n_mu) { j_hi = n_mu; j_lo = n_mu - 1; }   /* extrapolate high */

    const double mu_lo = mu_quad_pos[j_lo - 1];
    const double mu_hi = mu_quad_pos[j_hi - 1];
    const double t     = (mu_view - mu_lo) / (mu_hi - mu_lo);
    const double I_lo  = total[(size_t)0 * (size_t)dirs + (size_t)(+j_lo + n_mu)];
    const double I_hi  = total[(size_t)0 * (size_t)dirs + (size_t)(+j_hi + n_mu)];
    return (1.0 - t) * I_lo + t * I_hi;
}

/* 2026-07-13 LUT-PCHIP (Jae 승인): monotone cubic (PCHIP, Fritsch-Carlson
 * derivatives, Butland/scipy weighted-harmonic form) sampling of the TOA
 * per-m field at an off-node mu_view.  Replaces linear interpolation in the
 * LUT grid fill ONLY (single-view non-node path at L515/L860 and the LUT
 * downward T_diff grid remain linear — scope per instruction).
 * Verified vs external scipy PchipInterpolator emulation: linear max err
 * 0.138%(of I) -> PCHIP max 0.006% against zero-weight exact-node truth
 * (555/W2/SZA40 and 865/W10/SZA80, 63 pts each).  Outside the node range the
 * behaviour falls back to the SAME clamped linear extrapolation as
 * interp_view_at_toa (PCHIP is undefined for extrapolation; near-nadir
 * mu_view > last node keeps legacy behaviour). Cost: O(1) local 4-point
 * stencil per call (~2 ms per LUT run, 0.14% of runtime). */
static double pchip_derivative_interior(double h0, double h1,
                                        double s0, double s1) {
    if (s0 == 0.0 || s1 == 0.0 || (s0 > 0.0) != (s1 > 0.0)) return 0.0;
    const double w1 = 2.0 * h1 + h0;
    const double w2 = h1 + 2.0 * h0;
    return (w1 + w2) / (w1 / s0 + w2 / s1);
}
static double pchip_derivative_edge(double h0, double h1,
                                    double s0, double s1) {
    /* scipy _edge_case: one-sided three-point estimate + monotonicity clamp */
    double d = ((2.0 * h0 + h1) * s0 - h0 * s1) / (h0 + h1);
    if (d == 0.0 || s0 == 0.0 || (d > 0.0) != (s0 > 0.0)) {
        if (s0 == 0.0 || (d > 0.0) != (s0 > 0.0)) d = 0.0;
    } else if (((s0 > 0.0) != (s1 > 0.0)) && fabs(d) > 3.0 * fabs(s0)) {
        d = 3.0 * s0;
    }
    return d;
}
static double pchip_view_at_toa(const double *total, int n_mu,
                                const double *mu_quad_pos,
                                double mu_view) {
    const int dirs = 2 * n_mu + 1;
    /* Bracket exactly as interp_view_at_toa. */
    int j_hi = 1;
    while (j_hi < n_mu && mu_quad_pos[j_hi - 1] < mu_view) j_hi++;
    int j_lo = j_hi - 1;
    /* Out of range -> legacy clamped linear extrapolation (bit-parity). */
    if (j_lo < 1 || j_hi > n_mu || n_mu < 4)
        return interp_view_at_toa(total, n_mu, mu_quad_pos, mu_view);

    #define PCHIP_Y(j) total[(size_t)0 * (size_t)dirs + (size_t)(+(j) + n_mu)]
    const double x1 = mu_quad_pos[j_lo - 1];
    const double x2 = mu_quad_pos[j_hi - 1];
    const double y1 = PCHIP_Y(j_lo);
    const double y2 = PCHIP_Y(j_hi);
    const double h  = x2 - x1;
    const double s  = (y2 - y1) / h;

    double d1, d2;
    if (j_lo == 1) {
        const double x3 = mu_quad_pos[j_hi];      /* node j_hi+1 */
        const double s2 = (PCHIP_Y(j_hi + 1) - y2) / (x3 - x2);
        d1 = pchip_derivative_edge(h, x3 - x2, s, s2);
        d2 = pchip_derivative_interior(h, x3 - x2, s, s2);
    } else if (j_hi == n_mu) {
        const double x0 = mu_quad_pos[j_lo - 2];  /* node j_lo-1 */
        const double s0 = (y1 - PCHIP_Y(j_lo - 1)) / (x1 - x0);
        d1 = pchip_derivative_interior(x1 - x0, h, s0, s);
        d2 = pchip_derivative_edge(h, x1 - x0, s, s0);
    } else {
        const double x0 = mu_quad_pos[j_lo - 2];
        const double x3 = mu_quad_pos[j_hi];
        const double s0 = (y1 - PCHIP_Y(j_lo - 1)) / (x1 - x0);
        const double s2 = (PCHIP_Y(j_hi + 1) - y2) / (x3 - x2);
        d1 = pchip_derivative_interior(x1 - x0, h, s0, s);
        d2 = pchip_derivative_interior(h, x3 - x2, s, s2);
    }
    #undef PCHIP_Y

    /* Cubic Hermite on [x1, x2]. */
    const double t  = (mu_view - x1) / h;
    const double t2 = t * t, t3 = t2 * t;
    return (2.0 * t3 - 3.0 * t2 + 1.0) * y1 + (t3 - 2.0 * t2 + t) * h * d1
         + (-2.0 * t3 + 3.0 * t2) * y2 + (t3 - t2) * h * d2;
}


/* Replace the quadrature geometry in an already-built atmosphere with
 * `n_mu_gl` true Gauss-Legendre positive nodes plus one appended zero-weight
 * target-view node at index ±(n_mu_gl+1).
 *
 * This is the output-only target-node scheme: all original quadrature nodes and
 * weights remain untouched for the source integration, while the extra target
 * node receives a transported radiance solution but contributes zero incoming
 * quadrature weight to every source integral.  This avoids the older
 * --view-as-node behaviour that replaced the last GL node and degraded the
 * quadrature itself.
 */
/* v1.10 B-0a.3: fill the atm ring from the unified OSOAA-parity table
 * (sorted, deduped).  Replaces rt_solver_append_view_node at both call
 * sites; the old function is kept above for reference. */
static int rt_solver_fill_ring_from_uangles(rt_atm_t *atm,
                                            const rt_uangles_t *u) {
    if (!atm || !u || u->n_total < 2) return -1;
    if (atm->n_mu != u->n_total) return -1;
    for (int j = 1; j <= u->n_total; ++j) {
        atm->rm[+j] = +u->mu[j-1];
        atm->rm[-j] = -u->mu[j-1];
        atm->gb[+j] =  u->w [j-1];
        atm->gb[-j] =  u->w [j-1];
    }
    for (int j = 1; j <= u->n_total; ++j) {
        const double mu = fabs(atm->rm[+j]);
        const double xx = 1.0 - mu * mu;
        const double p2 = 0.5 * (3.0 * mu * mu - 1.0);
        const double r2 = 3.0 * xx / (2.0 * sqrt(6.0));
        atm->xpl[+j] = p2; atm->xpl[-j] = p2;
        atm->xrl[+j] = r2; atm->xrl[-j] = r2;
        atm->xtl[+j] = 0.0; atm->xtl[-j] = 0.0;
    }

    return 0;
}

static int rt_solver_append_view_node(rt_atm_t *atm, int n_mu_gl,
                                      double mu_view) {
    if (!atm || n_mu_gl < 2) return -1;
    const int n_mu_total = atm->n_mu;
    if (n_mu_total != n_mu_gl + 1) return -1;
    if (!(mu_view > 0.0 && mu_view <= 1.0)) return -1;
    if (n_mu_gl > 256) return -1;

    double mu_quad[256], w_quad[256];
    if (rt_quadrature_gauss_legendre_pos(n_mu_gl, mu_quad, w_quad) != 0)
        return -2;

    for (int j = 1; j <= n_mu_gl; ++j) {
        const double mu = mu_quad[j - 1];
        const double wt = w_quad[j - 1];
        atm->rm[+j] = +mu;
        atm->rm[-j] = -mu;
        atm->gb[+j] = wt;
        atm->gb[-j] = wt;
    }

    const int jt = n_mu_total;
    atm->rm[+jt] = +mu_view;
    atm->rm[-jt] = -mu_view;
    atm->gb[+jt] = 0.0;
    atm->gb[-jt] = 0.0;

    /* Initialize m=0 auxiliary values.  These are overwritten per Fourier mode
     * by rt_legendre_compute{,_pol}, but keeping them consistent avoids stale
     * GL(n_mu_gl+1) values in diagnostic paths that inspect the atmosphere
     * immediately after build. */
    for (int j = 1; j <= n_mu_total; ++j) {
        const double mu = fabs(atm->rm[+j]);
        const double xx = 1.0 - mu * mu;
        const double p2 = 0.5 * (3.0 * mu * mu - 1.0);
        const double r2 = 3.0 * xx / (2.0 * sqrt(6.0));
        atm->xpl[+j] = p2; atm->xpl[-j] = p2;
        atm->xrl[+j] = r2; atm->xrl[-j] = r2;
        atm->xtl[+j] = 0.0; atm->xtl[-j] = 0.0;
    }

    return 0;
}

/* Attach a shaped external Fourier bottom source to a solver option object.
 * Keeping the pointer triplet and its dimensions together at every call site
 * is the central invariant of OCRT-RT-C3D-EXTBOTTOM-001. */
static void rt_options_set_ext_bottom_source(rt_options_t *o,
                                              const double *I,
                                              const double *Q,
                                              const double *U,
                                              int n_mu,
                                              int m_max) {
    if (!o) return;
    o->ext_bottom_per_m_I = I;
    o->ext_bottom_per_m_Q = Q;
    o->ext_bottom_per_m_U = U;
    o->ext_bottom_n_mu = n_mu;
    o->ext_bottom_m_max = m_max;
}

/* 2026-07-14 M1 (Jae 지시): 스칼라 I-전용 솔버 rt_solve_case 를 완전
 * 삭제했다.  향후 스칼라 경로 검증 계획이 없고, 모든 실행은 벡터
 * rt_solve_case_pol 로 통일된다.  삭제 전 본체는 v1.11 패키지
 * (wpchip 이전) 소스에 보존되어 있다. */

/* Step 6: vector RT end-to-end (rt_solve_case_pol).
 *
 * Mirrors rt_solve_case but threads I, Q, U through the per-m loop:
 *   - Each m step calls the *_pol Legendre + kernel + primary-source
 *     vector functions (Steps 2b, 2c, 3) and rt_solver_sos_pol (Step 5c).
 *   - Three-buffer field storage (src/primary/total per Stokes component);
 *     I path is independent of the Q,U path's SOS coupling at m=0
 *     (vector kernel has gt=art=att=0 there) and therefore I rho_I match
 *     to the scalar path bit-for-bit at m=0; m>=1 differs only at the
 *     final reflectance via vector vs scalar SOS source coupling, but
 *     the scalar I source produced from rt_solver_primary_source is
 *     identical to the I component of the vector primary source by
 *     construction (see rt_solver_primary_source_pol doc).
 *
 * View sampling identical to the scalar path: linear interpolation
 * between bracketing positive-mu GL nodes (default) or the appended
 * zero-weight output-only target node (--view-as-node).
 *
 * Reconstruction:
 *   I, Q  → rt_solver_reconstruct_phi      (cos)
 *   U     → rt_solver_reconstruct_phi_sin  (sin)
 *
 * No new error code beyond rt_solve_case: -1 on bad args, -2 on NaN.
 */
/* Internal implementation: rt_solve_case_pol with optional aerosol input.
 * Public wrappers below dispatch to this function.
 *
 * boa_export (optional, may be NULL):
 *   When non-NULL and pre-allocated with consistent (m_max, n_mu), receives
 *   the diffuse BOA downward I^m, Q^m, U^m field at each positive μ_j node.
 *   No effect on the main TOA computation — bit-exact preservation when NULL. */
static int rt_solve_case_pol_impl(const rt_case_t *cs, const rt_options_t *opts,
                                    const rt_aerosol_input_t *aer,
                                    rt_result_t *out,
                                    rt_lut_grid_out_t *lut_out,
                                    rt_atm_boa_export_t *boa_export) {
    if (!cs || !out) return -1;
    rt_options_t opts_default = rt_options_default();
    if (!opts) opts = &opts_default;

    if (!rt_spectral_wavelength_supported(cs->wavelength_nm)) return -1;
    if (cs->sza_deg < 0.0 || cs->sza_deg >= 90.0) return -1;
    if (cs->vza_deg < 0.0 || cs->vza_deg >= 90.0) return -1;
    if (opts->fourier_m_max < 0 || opts->fourier_m_max > 32) return -1;
    if (opts->n_mu < 2 || opts->n_layers < 2)               return -1;
    if (aer) {
        if (aer->tau_a < 0.0)                  return -1;
        if (aer->ssa_a < 0.0 || aer->ssa_a > 1.0) return -1;
        if (aer->L_max < 2)                    return -1;
        if (!aer->betal_aer || !aer->gammal_aer
            || !aer->alphal_aer || !aer->zetal_aer) return -1;
    }

    memset(out, 0, sizeof(*out));

    const double pi    = RT_F_SOLAR_PI;
    const double mu_solar = cos(cs->sza_deg * pi / 180.0);
    const double mu_view  = cos(cs->vza_deg * pi / 180.0);
    const double dphi     = rt_raa_to_atm_fourier_phi(cs->raa_deg);
    const double depol    = 0.0279;

    double tau_R;
    if (cs->rayleigh_on == 0) {
        /* Phase 3 aerosol-only mode: force tau_R = 0 (matches V1 --no-rayleigh,
         * OSOAA AP.MOT=0). Without this override, build_aerosol_rayleigh always
         * mixes in Rayleigh from rt_rayleigh_tau_model. */
        tau_R = 0.0;
    } else if (opts->tau_r_source == RT_TAU_R_FROM_INPUT) {
        if (cs->tau_R_input <= 0.0) return -1;
        tau_R = cs->tau_R_input;
    } else {
        tau_R = rt_rayleigh_tau_model_full(cs->wavelength_nm, cs->pressure_hpa, 45.0, 0.0, opts->rayleigh_model);
    }
    out->tau_R_total = tau_R;

    const int    n_mu_requested = opts->n_mu;
    const int    append_view_node = opts->view_as_node ? 1 : 0;
    /* v1.10 B-0a.3: atm ring via the unified table.  The view angle is
     * inserted at its SORTED position and DEDUPED against the GL nodes
     * (full-grid rows: mu_view == a GL node to <1e-15, so the ring keeps
     * n_gauss nodes and the view index aliases that GL node - the OSOAA
     * IMUS rule).  n_mu is the TRUE ring size. */
    rt_uangles_t ocrt_uatm; int ocrt_view_ring_idx = 0;
    int n_mu = n_mu_requested;
    if (append_view_node) {
        int ocrt_vidx0 = -1;
        if (rt_uangles_init(&ocrt_uatm, n_mu_requested) != 0) return -1;
        if (rt_uangles_add(&ocrt_uatm, mu_view, &ocrt_vidx0) != 0) return -1;
        n_mu = ocrt_uatm.n_total;
        ocrt_view_ring_idx = ocrt_vidx0 + 1;
    }


    /* External bottom-source shape contract.
     *
     * Before this guard, the atmosphere loop sliced the arrays at m*n_mu
     * without knowing how many modes the producer had allocated.  A coupled
     * ocean run with water_m_max < atmospheric m_max therefore read beyond the
     * water-leaving array and could return finite-looking or enormous corrupt
     * TOA Stokes values while reporting convergence. */
    const int ext_any = (opts->ext_bottom_per_m_I != NULL) ||
                        (opts->ext_bottom_per_m_Q != NULL) ||
                        (opts->ext_bottom_per_m_U != NULL);
    const int ext_all = (opts->ext_bottom_per_m_I != NULL) &&
                        (opts->ext_bottom_per_m_Q != NULL) &&
                        (opts->ext_bottom_per_m_U != NULL);
    if (ext_any && !ext_all) {
        fprintf(stderr,
                "rt_solver: external bottom source requires I/Q/U together\n");
        return RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE;
    }
    if (ext_all &&
        (opts->ext_bottom_n_mu != n_mu || opts->ext_bottom_m_max < 0)) {
        fprintf(stderr,
                "rt_solver: external bottom source shape mismatch "
                "(source n_mu=%d m_max=%d, solver n_mu=%d m_max=%d)\n",
                opts->ext_bottom_n_mu, opts->ext_bottom_m_max,
                n_mu, opts->fourier_m_max);
        return RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE;
    }
    if(getenv("DBG_BOA"))fprintf(stderr,"NMUDET nreq=%d nmu=%d view=%d\n",n_mu_requested,n_mu,opts->view_as_node);
    const int    nt    = opts->n_layers;
    const int    dirs  = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;
    /* Phase 3: when aerosol active, l_max must reach aer->L_max so that
     * rt_kernel_phase_fourier{,_pol} can sum to that order. */
    int l_max = (opts->fourier_m_max > 2) ? opts->fourier_m_max : 2;
    if (aer && aer->tau_a > 0.0 && aer->L_max > l_max) l_max = aer->L_max;
    const int    m_max = opts->fourier_m_max;

    rt_atm_t                  atm  = {0};
    rt_legendre_workspace_t   ws   = {0};
    /* Three Stokes components × {src, primary, total} = 9 buffers. */
    double *src_i    = NULL, *src_q    = NULL, *src_u    = NULL;
    double *prim_i   = NULL, *prim_q   = NULL, *prim_u   = NULL;
    double *total_i  = NULL, *total_q  = NULL, *total_u  = NULL;
    int     rc       = 0;

    if (rt_atm_alloc(&atm, nt, n_mu) != 0) { rc = -1; goto cleanup; }
    /* Phase 3: aerosol-aware atmosphere build. */
    if (aer && aer->tau_a > 0.0) {
        if (rt_atm_build_aerosol_rayleigh(&atm, tau_R, depol,
                                          aer->tau_a, aer->ssa_a, aer->L_max,
                                          aer->betal_aer, aer->gammal_aer,
                                          aer->alphal_aer, aer->zetal_aer,
                                          mu_solar, opts->rayleigh_model,
                                          opts->aer_h_km) != 0) {
            rc = -1; goto cleanup;
        }
        /* Borrow the (loglin-capped) angle-space phase for the value kernel.
         * Pointers are owned by aer; rt_atm_free must NOT free them. */
        atm.aer_use_value_kernel = aer->use_value_kernel;
        atm.aer_n_ang_phase      = aer->n_ang_phase;
        atm.aer_theta_phase      = aer->theta_phase;
        atm.aer_P11_phase        = aer->P11_phase;
        atm.aer_P12_phase        = aer->P12_phase;
        atm.aer_P33_phase        = aer->P33_phase;
    } else {
        if (getenv("OCRT_DUMP_TAUR"))
        fprintf(stderr, "[TAUR] wl=%g tau_R=%.10f depol=%.8f\n",
                cs->wavelength_nm, tau_R, depol);
    if (rt_atm_build_rayleigh(&atm, tau_R, depol, mu_solar,
                                  opts->rayleigh_model) != 0) {
            rc = -1; goto cleanup;
        }
    }
    /* v1.02: apply SOS-integrated gas absorption if requested. */
    if (opts->abs_state) {
        if (rt_atm_apply_gas_absorption(&atm,
                (const struct rt_absorption_state *)opts->abs_state,
                cs->wavelength_nm) != 0) {
            rc = -1; goto cleanup;
        }
    }

    /* PSSA 2026-07-18: apply after gas absorption finalizes h[].
     * v1.11 (2026-09-05): --pssa selects IPSS (Zhai & Hu 2022); the legacy
     * average-secant Chapman path was DELETED from the tree.  The SOS solve
     * therefore stays plane-parallel unconditionally, which is what Eq. (7)
     * requires because it rescales I_pp (see rt_ipss.h). */
    rt_ipss_t ipss_geo;
    int ipss_ready = 0;
    double ipss_tau_dir_dn = -1.0;
    memset(&ipss_geo, 0, sizeof(ipss_geo));
    if (opts->pssa) {
        const int irc = rt_ipss_init_from_atm(&ipss_geo, &atm,
                                              RT_IPSS_EARTH_RADIUS_KM);
        if (irc != 0) {
            fprintf(stderr,
                    "error: IPSS geometry setup failed (rc=%d).  A strictly "
                    "decreasing altitude grid is required.\n", irc);
            rc = -1; goto cleanup;
        }
        ipss_ready = 1;
        ipss_tau_dir_dn = rt_ipss_tau_dir_dn_boa(&ipss_geo);
        {   /* Diagnostic (env-gated, inert by default): dump the layer profile
             * that feeds kappa -- altitude edges, total extinction, aerosol /
             * Rayleigh scattering fractions, gas absorption -- plus the aerosol
             * P11 Legendre moments and depol, so that kappa can be recomputed
             * outside the solver (e.g. with phase-function-weighted sources). */
            const char *pdmp = getenv("OCRT_IPSS_PROFILE_DUMP");
            if (pdmp && pdmp[0]) {
                FILE *pf = fopen(pdmp, "w");
                if (pf) {
                    const int nt = atm.n_layers;
                    fprintf(pf, "# depol=%.10g L_max=%d aerosol_active=%d\n",
                            atm.depol, atm.L_max, atm.aerosol_active);
                    fprintf(pf, "# betal_aer:");
                    if (atm.betal_aer)
                        for (int l = 0; l <= atm.L_max; ++l)
                            fprintf(pf, " %.12g", atm.betal_aer[l]);
                    fprintf(pf, "\n");
                    /* tabulated P11 (the production value kernel), if present */
                    fprintf(pf, "# p11_table n=%d use_value_kernel=%d:",
                            atm.aer_n_ang_phase, atm.aer_use_value_kernel);
                    if (atm.aer_theta_phase && atm.aer_P11_phase)
                        for (int a = 0; a < atm.aer_n_ang_phase; ++a)
                            fprintf(pf, " %.6g:%.9g", atm.aer_theta_phase[a],
                                    atm.aer_P11_phase[a]);
                    fprintf(pf, "\n");
                    fprintf(pf, "# ipss_weights=%s P_R/P_A at (vza,raa) are set per geometry\n",
                            ipss_geo.w_src ? "phase-weighted" : "phase-free");
                    fprintf(pf, "k,z_top_km,z_bot_km,dtau,xdel,ydel,tau_abs\n");
                    for (int k = 0; k < nt; ++k)
                        fprintf(pf, "%d,%.9f,%.9f,%.12e,%.12e,%.12e,%.12e\n", k,
                                atm.z_km_level[k], atm.z_km_level[k + 1],
                                atm.h[k + 1] - atm.h[k],
                                atm.xdel ? atm.xdel[k + 1] : 0.0,
                                atm.ydel ? atm.ydel[k + 1] : 1.0,
                                atm.tau_abs_layer ? atm.tau_abs_layer[k] : 0.0);
                    fclose(pf);
                }
            }
        }
        {   /* traceability: campaigns record stderr, and the bit gates compare
             * CSV outputs only, so this banner is regression-inert. */
            static int ipss_banner = 0;
            if (!ipss_banner) {
                ipss_banner = 1;
                fprintf(stderr,
                        "# spherical=ipss (Zhai & Hu 2022 JQSRT 282,108132; "
                        "Re=%.1f km, H=%.1f km, %d layers, %d nodes/panel)\n",
                        ipss_geo.Re, ipss_geo.H, ipss_geo.nlay, ipss_geo.n_quad);
            }
        }
    }

    if (append_view_node) {
        if (rt_solver_fill_ring_from_uangles(&atm, &ocrt_uatm) != 0) {
            rc = -1; goto cleanup;
        }
    }

    if (rt_legendre_workspace_alloc(&ws, n_mu, l_max) != 0) {
        rc = -1; goto cleanup;
    }

    src_i   = (double *)calloc(buf_size, sizeof(double));
    src_q   = (double *)calloc(buf_size, sizeof(double));
    src_u   = (double *)calloc(buf_size, sizeof(double));
    prim_i  = (double *)calloc(buf_size, sizeof(double));
    prim_q  = (double *)calloc(buf_size, sizeof(double));
    prim_u  = (double *)calloc(buf_size, sizeof(double));
    total_i = (double *)calloc(buf_size, sizeof(double));
    total_q = (double *)calloc(buf_size, sizeof(double));
    total_u = (double *)calloc(buf_size, sizeof(double));
    if (!src_i || !src_q || !src_u || !prim_i || !prim_q || !prim_u
            || !total_i || !total_q || !total_u) {
        rc = -1; goto cleanup;
    }

    double I_total_per_m[33] = {0};
    double Q_total_per_m[33] = {0};
    double U_total_per_m[33] = {0};

    /* v1.01: capture downward-at-BOA intensity per Fourier mode m and per
     * quadrature node (positive μ index). Used for atmospheric diffuse
     * downward transmittance calculation. Layout: [m * n_mu + (j_pos-1)]. */
    double *I_dn_BOA_per_m = (double*)calloc((size_t)(m_max + 1) * n_mu, sizeof(double));
    if (!I_dn_BOA_per_m) { rc = -1; goto cleanup; }

    int    max_orders_seen = 0;
    int    all_converged   = 1;

    /* v1.01 LUT grid buffers: per-m sampled values at each vza grid point.
     * Layout: [m * n_vza + iv]. Allocated only when lut_out != NULL. */
    double *I_per_m_grid = NULL;
    double *Q_per_m_grid = NULL;
    double *U_per_m_grid = NULL;
    double *mu_v_grid    = NULL;  /* cos(vza_grid) precomputed */
    if (lut_out) {
        const int n_vza_g = lut_out->n_vza;
        if (n_vza_g <= 0 || !lut_out->vza_deg || lut_out->n_raa <= 0 || !lut_out->raa_deg ||
            !lut_out->rho_I || !lut_out->rho_Q || !lut_out->rho_U) {
            rc = -1; goto cleanup;
        }
        I_per_m_grid = (double*)calloc((size_t)(m_max + 1) * n_vza_g, sizeof(double));
        Q_per_m_grid = (double*)calloc((size_t)(m_max + 1) * n_vza_g, sizeof(double));
        U_per_m_grid = (double*)calloc((size_t)(m_max + 1) * n_vza_g, sizeof(double));
        mu_v_grid    = (double*)calloc((size_t)n_vza_g, sizeof(double));
        if (!I_per_m_grid || !Q_per_m_grid || !U_per_m_grid || !mu_v_grid) {
            free(I_per_m_grid); free(Q_per_m_grid); free(U_per_m_grid); free(mu_v_grid);
            I_per_m_grid = Q_per_m_grid = U_per_m_grid = mu_v_grid = NULL;
            rc = -1; goto cleanup;
        }
        for (int iv = 0; iv < n_vza_g; ++iv) {
            double v = lut_out->vza_deg[iv];
            if (v < 0.0 || v >= 90.0) {
                rc = -1; goto cleanup;
            }
            mu_v_grid[iv] = cos(v * pi / 180.0);
        }
    }

    surface_fkc_set_build_mmax(m_max);   /* v1.11-speed S1: cap air-kernel all-m build */
    for (int m = 0; m <= m_max; ++m) {
        /* Phase 3: pass atm->{gammal,betal}_aer when aerosol active; NULL
         * preserves bit-level Phase 2 Rayleigh-only behavior. */
        const double *gammal_aer = atm.aerosol_active ? atm.gammal_aer : NULL;
        const double *betal_aer  = atm.aerosol_active ? atm.betal_aer  : NULL;
        const double *alphal_aer = atm.aerosol_active ? atm.alphal_aer : NULL;
        const double *zetal_aer  = atm.aerosol_active ? atm.zetal_aer  : NULL;
        if (rt_legendre_compute_pol(&ws, &atm, m) != 0)              { rc = -1; goto cleanup; }
        if (rt_kernel_phase_fourier_pol(&ws, m, gammal_aer) != 0)    { rc = -1; goto cleanup; }
        /* Phase 3: aerosol arr/art/att kernels (KERNELPOL.f L160-178).
         * Zero-init when aerosol inactive — bit-level Phase 2 preserved. */
        if (rt_kernel_phase_fourier_aerosol_full(&ws, m, alphal_aer,
                                                  zetal_aer) != 0) {
            rc = -1; goto cleanup;
        }
        /* rt_kernel_phase_fourier (scalar) is also called inside
         * rt_legendre_compute_pol → rt_legendre_compute, but the I-channel
         * phase function must be re-summed here because the legacy scalar
         * primary needs ws->phase_fourier_m.  Calling _pol after _scalar
         * is harmless (idempotent on rrl/rtl, refills phase_fourier_m). */
        if (atm.aerosol_active && atm.aer_use_value_kernel && atm.aer_P11_phase) {
            /* Gibbs-free angle-space value kernel for the scalar I phase
             * (fixes the backscatter glory the moment kernel cannot represent). */
            if (aerosol_value_phase_fourier(&ws, &atm, m, 720) != 0)    { rc = -1; goto cleanup; }
        } else {
            if (rt_kernel_phase_fourier(&ws, m, betal_aer) != 0)        { rc = -1; goto cleanup; }
        }
        if (rt_sos_operator_prepare(&atm, m, &ws) != 0) {
            rc = -1; goto cleanup;
        }

        /* Primary sources (I, Q, U). */
        if (rt_solver_primary_source(&atm, m, &ws, src_i) != 0) {
            rc = -1; goto cleanup;
        }
        if (rt_solver_primary_source_pol(&atm, m, &ws, src_q, src_u) != 0) {
            rc = -1; goto cleanup;
        }

        /* Vertical integration × 3 (source-agnostic — Step 4 verified). */
        if (rt_solver_integrate(&atm, m, src_i, opts->integration_method,
                                prim_i) != 0) { rc = -2; goto cleanup; }
        if (rt_solver_integrate(&atm, m, src_q, opts->integration_method,
                                prim_q) != 0) { rc = -2; goto cleanup; }
        if (rt_solver_integrate(&atm, m, src_u, opts->integration_method,
                                prim_u) != 0) { rc = -2; goto cleanup; }
        /* C3 pass-2: zero the solar primary so ONLY the injected bottom source
         * (ext_bottom, applied inside the atm SOS) drives the field. */
        if (opts->bottom_source_only) {
            memset(prim_i, 0, buf_size * sizeof(double));
            memset(prim_q, 0, buf_size * sizeof(double));
            memset(prim_u, 0, buf_size * sizeof(double));
        }

        /* Vector SOS iteration (Step 5c).
         * Surface-aware path: COXMUNK 또는 FLAT 둘 다 surface-aware solver.
         * - RT_SURFACE_FLAT: 항상 flat Fresnel (wind 무관)
         * - RT_SURFACE_BLACK_FRESNEL_OCEAN: wind 분기 (W<=0 → flat, W>0 → Cox-Munk)
         * - RT_SURFACE_BLACK / LAMBERT: atm-only solver */
        rt_solver_sos_options_t sos_opt = opts->sos;
        rt_solver_sos_result_t  sos_res = {0};
        int sos_rc;
        /* C3 pass-2: slice the per-mode water-leaving 0+ field into the SOS
         * bottom source for this mode (no-op when NULL / normal solar solve). */
        sos_opt.bottom_source_only = opts->bottom_source_only;
        /* Never inherit a stale slice from opts->sos.  A source producer may
         * legitimately provide fewer Fourier modes than the atmosphere; the
         * missing modes are exactly zero and therefore remain NULL here. */
        sos_opt.ext_bottom_I = NULL;
        sos_opt.ext_bottom_Q = NULL;
        sos_opt.ext_bottom_U = NULL;
        if (ext_all && m <= opts->ext_bottom_m_max) {
            const size_t off = (size_t)m * (size_t)opts->ext_bottom_n_mu;
            sos_opt.ext_bottom_I = opts->ext_bottom_per_m_I + off;
            sos_opt.ext_bottom_Q = opts->ext_bottom_per_m_Q + off;
            sos_opt.ext_bottom_U = opts->ext_bottom_per_m_U + off;
        }
        if (cs->surface == RT_SURFACE_BLACK_FRESNEL_OCEAN ||
            cs->surface == RT_SURFACE_FLAT) {
            sos_rc = rt_solver_sos_pol_with_surface(
                &atm, m, &ws, prim_i, prim_q, prim_u, &sos_opt,
                cs->surface,
                cs->n_water, cs->wind_speed, cs->sigma_type, cs->q_convention,
                total_i, total_q, total_u, &sos_res);
        } else {
            sos_rc = rt_solver_sos_pol(&atm, m, &ws, prim_i, prim_q, prim_u,
                                        &sos_opt, total_i, total_q, total_u,
                                        &sos_res);
        }
        if (sos_rc != 0) {
            rc = -2; goto cleanup;
        }
        out->sos_orders_per_m[m]   = sos_res.n_orders_used;
        out->sos_residual_per_m[m] = sos_res.final_residual;
        if (sos_res.n_orders_used > max_orders_seen)
            max_orders_seen = sos_res.n_orders_used;
        if (!sos_res.converged) all_converged = 0;

        if (ocrt_debug_env("OCRT_DUMP_PASS2TOA") && opts->bottom_source_only && m==0) {
            fprintf(stderr,"[PASS2TOA] m=0 n_mu=%d total_i at TOA upward nodes:", n_mu);
            for (int k=0;k<n_mu;k++){ size_t ix=(size_t)0*(size_t)dirs+(size_t)(n_mu+1+k);
                fprintf(stderr," [mu=%.3f]%.4e", atm.rm[+(k+1)], total_i[ix]); }
            fprintf(stderr,"\n");
        }

        /* Sample TOA upwelling at mu_view — same two-path (output-only
         * appended node vs linear interp) as scalar rt_solve_case, applied
         * to all three Stokes components. */
        if (opts->view_as_node) {
            const size_t idx = (size_t)0 * (size_t)dirs +
                               (size_t)(+n_mu + ocrt_view_ring_idx);
            I_total_per_m[m] = total_i[idx];
            Q_total_per_m[m] = total_q[idx];
            U_total_per_m[m] = total_u[idx];
        } else {
            const double *mu_quad_pos = &atm.rm[+1];
            I_total_per_m[m] = interp_view_at_toa(total_i, n_mu,
                                                  mu_quad_pos, mu_view);
            Q_total_per_m[m] = interp_view_at_toa(total_q, n_mu,
                                                  mu_quad_pos, mu_view);
            U_total_per_m[m] = interp_view_at_toa(total_u, n_mu,
                                                  mu_quad_pos, mu_view);
        }

        /* v1.01 LUT: sample TOA upwelling at every vza grid point.
         * 2026-07-13 LUT-PCHIP (Jae 승인): linear -> monotone cubic (PCHIP)
         * per-m sampling; see pchip_view_at_toa.  Off-node max error drops
         * 0.138% -> 0.006% (of I) vs exact-node truth; out-of-range falls
         * back to the legacy linear extrapolation. */
        if (lut_out) {
            const double *mu_quad_pos = &atm.rm[+1];
            const int n_vza_g = lut_out->n_vza;
            for (int iv = 0; iv < n_vza_g; ++iv) {
                I_per_m_grid[m * n_vza_g + iv] =
                    pchip_view_at_toa(total_i, n_mu, mu_quad_pos, mu_v_grid[iv]);
                Q_per_m_grid[m * n_vza_g + iv] =
                    pchip_view_at_toa(total_q, n_mu, mu_quad_pos, mu_v_grid[iv]);
                U_per_m_grid[m * n_vza_g + iv] =
                    pchip_view_at_toa(total_u, n_mu, mu_quad_pos, mu_v_grid[iv]);
            }
        }

        /* v1.01: capture I_dn^m at BOA (level k = nt) at each positive
         * quadrature node. Index: total_i[k=nt][j=-j_pos + n_mu]
         * (negative-j slot stores downward intensity by convention).
         *
         * Phase B.4 Stage 2b: also capture Q_dn^m, U_dn^m and (if boa_export
         * non-NULL) copy I/Q/U^m into the export struct for ocean coupling. */
        {
            const int nt = atm.n_layers;
            const int dirs_ = 2 * n_mu + 1;
            for (int jp = 1; jp <= n_mu; ++jp) {
                const size_t idx_dn = (size_t)nt * (size_t)dirs_ +
                                       (size_t)(-jp + n_mu);
                I_dn_BOA_per_m[m * n_mu + (jp - 1)] = total_i[idx_dn];
                if (getenv("DBG_BOA") && jp <= 3)
                    fprintf(stderr, "DBGBOA m=%d jp=%d idx_dn=%zu I_dn=%.6e Q_dn=%.6e\n",
                            m, jp, idx_dn, total_i[idx_dn], total_q[idx_dn]);

                if (getenv("DBG_BOA") && jp == 1)
                    fprintf(stderr, "BOACOND m=%d be=%p I_per_m=%p n_mu=%d be_n_mu=%d be_mmax=%d\n",
                            m, (void*)boa_export,
                            boa_export ? (void*)boa_export->I_per_m : NULL,
                            n_mu, boa_export ? boa_export->n_mu : -1,
                            boa_export ? boa_export->m_max : -1);
                if (boa_export &&
                    boa_export->I_per_m && boa_export->Q_per_m && boa_export->U_per_m &&
                    boa_export->n_mu == n_mu &&
                    m <= boa_export->m_max) {
                    boa_export->I_per_m[m * n_mu + (jp - 1)] = total_i[idx_dn];
                    boa_export->Q_per_m[m * n_mu + (jp - 1)] = total_q[idx_dn];
                    boa_export->U_per_m[m * n_mu + (jp - 1)] = total_u[idx_dn];
                }
            }
        }
    }

    /* Phase B.4 Stage 2b: fill boa_export metadata after m loop completes.
     * Copies positive μ grid and τ_atm metadata for the caller's coupling step. */
    if (boa_export && boa_export->mu_quad_pos && boa_export->n_mu == n_mu) {
        for (int jp = 1; jp <= n_mu; ++jp) {
            boa_export->mu_quad_pos[jp - 1] = atm.rm[+jp];
        }
        boa_export->tau_atm_total = atm.tau_total + atm.tau_abs_total;
        boa_export->mu_sun_air    = mu_solar;
    }

    /* Export the per-view Fourier samples used by the LUT reconstruction.
     * This is a pure copy of arrays already built above; it adds no SOS work. */
    if (lut_out && lut_out->view_per_m_I && lut_out->view_per_m_Q &&
        lut_out->view_per_m_U && lut_out->view_m_max >= m_max) {
        const size_t npm = (size_t)(m_max + 1) * (size_t)lut_out->n_vza;
        memcpy(lut_out->view_per_m_I, I_per_m_grid, npm * sizeof(double));
        memcpy(lut_out->view_per_m_Q, Q_per_m_grid, npm * sizeof(double));
        memcpy(lut_out->view_per_m_U, U_per_m_grid, npm * sizeof(double));
        lut_out->view_m_max_filled = m_max;
    } else if (lut_out) {
        lut_out->view_m_max_filled = -1;
    }

    /* Fourier reconstruction at view azimuth: I/Q use cosine modes and U
     * uses sine modes under the OCRT meridian-basis convention. */
    if (opts->toa_view_per_m_I) {   /* v1.10 S17 export (bit-inert; see rt_types.h) */
        memcpy(opts->toa_view_per_m_I, I_total_per_m, (size_t)(m_max+1)*sizeof(double));
        memcpy(opts->toa_view_per_m_Q, Q_total_per_m, (size_t)(m_max+1)*sizeof(double));
        memcpy(opts->toa_view_per_m_U, U_total_per_m, (size_t)(m_max+1)*sizeof(double));
    }
    const double I_TOA = rt_solver_reconstruct_phi    (I_total_per_m, m_max, dphi);
    const double Q_TOA = rt_solver_reconstruct_phi    (Q_total_per_m, m_max, dphi);
    /* OCRT reports U in its fixed meridian-basis sine convention. */
    double U_TOA = rt_solver_reconstruct_phi_sin(U_total_per_m, m_max, dphi);

    out->I_TOA         = I_TOA;
    out->Q_TOA         = Q_TOA;
    out->U_TOA         = U_TOA;
    out->rho_I         = rt_solver_reflectance_from_intensity(I_TOA, mu_solar);
    out->rho_Q         = rt_solver_reflectance_from_intensity(Q_TOA, mu_solar);
    out->rho_U         = rt_solver_reflectance_from_intensity(U_TOA, mu_solar);

    /* v1.11 IPSS Eq. (7): rescale the plane-parallel atmospheric-path Stokes
     * vector by kappa = I_1,ss / I_1,pp for this exact (VZA, RAA).  Applied in
     * ANGLE space (not per Fourier mode) because tau_sun along the line of
     * sight depends on azimuth.  The direct sunglint term below keeps its own
     * spherical slant factor. */
    if (ipss_ready) {
        double kappa = 1.0; int kguard = 0;
        rt_ipss_kappa(&ipss_geo, cs->vza_deg, cs->raa_deg,
                      &kappa, NULL, NULL, &kguard);
        out->rho_I *= kappa;
        out->rho_Q *= kappa;
        out->rho_U *= kappa;
        out->I_TOA *= kappa;
        out->Q_TOA *= kappa;
        out->U_TOA *= kappa;
        if (getenv("OCRT_IPSS_TRACE"))
            fprintf(stderr, "[ipss] vza=%.3f raa=%.3f kappa=%.9f guard=%d\n",
                    cs->vza_deg, cs->raa_deg, kappa, kguard);
    }

    /* Direct sunglint contribution (Phase 4): when surface is Cox-Munk and
     * decouple_sunglint = 0, ADD the direct surface specular bounce of the
     * solar beam to the TOA reflectance.
     *
     * For Cox-Munk wind > 0: slope-weighted Fresnel kernel.
     * For wind = 0: caller should use surface = FLAT or accept that
     *               sigma_type=1 gives σ²=0.003 finite-width "almost flat"
     *               specular peak. AF1982 reference uses true flat Fresnel
     *               (delta function at exact specular) when wind=0, which
     *               cannot be represented by Cox-Munk integration here.
     *               Our with-glint output for wind=0 will diverge sharply
     *               at exact specular (sza=vza, raa=180 in V3 conv = 0 in
     *               AF1982 conv). */
    if (cs->surface == RT_SURFACE_BLACK_FRESNEL_OCEAN && cs->decouple_sunglint == 0) {
        /* Convert input raa back to AF convention for sunglint geometry.
         * V3 input raa = (raa_AF + 180) mod 360. Specular condition in V3
         * conv requires cos(Δφ) = -1 (Δφ = π) but the surface kernel
         * surface_R_coxmunk_trig uses scatter cosine cT including
         * mu_v·(-mu_0) (downward sign convention), which yields specular
         * at cos_phi = +1 (i.e. φ_v = φ_0). We pass phi_v as the AF
         * convention φ to align: phi_v = (raa_input - 180) deg. */
        const double phi_v = rt_raa_to_surface_view_phi(cs->raa_deg);
        const double phi_0 = 0.0;       /* solar at φ=0 reference */
        const double tau_total = atm.h[atm.n_layers]; /* τ at surface */
        const double mu_v = cs->vza_deg > 0 ?
                            cos(cs->vza_deg * RT_F_SOLAR_PI / 180.0) : 1.0;

        double rho_glint[3];
        surface_direct_sunglint_rho(mu_v, phi_v, mu_solar, phi_0,
                                     cs->wind_speed, cs->sigma_type,
                                     cs->n_water, tau_total,
                                     cs->q_convention, rho_glint);
        if (ipss_ready && ipss_tau_dir_dn >= 0.0) {
            /* exact ray-traced slant depth of the direct beam */
            const double c_glint = exp(-(ipss_tau_dir_dn - tau_total / mu_solar));
            rho_glint[0] *= c_glint; rho_glint[1] *= c_glint; rho_glint[2] *= c_glint;
        }
        out->rho_I += rho_glint[0];
        out->rho_Q += rho_glint[1];
        out->rho_U += -rho_glint[2];   /* match U sign convention */
    }
    out->n_orders_used = max_orders_seen;
    out->converged     = all_converged;

    /* v1.01: Atmospheric transmittance decomposition.
     *
     *   T_dir_dn       = exp(-τ_tot/μ_sun)
     *   F_TOA_dn_solar = π·μ_sun     (F_solar=π convention)
     *   T_diff_dn_hemi = (2π Σ_quad μ·w·I_dn^{m=0}(μ)) / (π·μ_sun)
     *
     * Note: I^m at BOA carries the (2-δ_m0) Fourier-mode normalization
     * embedded in I_dn(μ,φ) = Σ_m (2-δ_m0)·I_dn^m(μ)·cos(m·φ). For
     * hemispheric integration over φ, only m=0 survives, so we use the
     * m=0 slice directly. */
    {
        const double tau_tot = atm.h[atm.n_layers];
        out->T_dir_dn = (ipss_ready && ipss_tau_dir_dn >= 0.0)
                          ? exp(-ipss_tau_dir_dn)
                          : exp(-tau_tot / mu_solar);

        double F_dn_diff = 0.0;
        for (int jp = 1; jp <= n_mu; ++jp) {
            const double mu_q = atm.rm[+jp];
            const double w_q  = atm.gb[+jp];
            const double I0   = I_dn_BOA_per_m[0 * n_mu + (jp - 1)];
            F_dn_diff += mu_q * w_q * I0;
        }
        F_dn_diff *= 2.0 * pi;
        out->T_diff_dn_hemi = F_dn_diff / (pi * mu_solar);
        out->T_total_dn_hemi = out->T_dir_dn + out->T_diff_dn_hemi;

        /* Directional T_diff_dn at single-view (μ_view, raa). */
        const double *mu_quad_pos = &atm.rm[+1];
        double I_dn_per_m_view[33] = {0};
        for (int m = 0; m <= m_max; ++m) {
            /* Build a "BOA-only" total_i-shaped slice for interp helper:
             * since interp_view_at_toa reads index [k=0 * dirs + ...], we
             * need a small inline equivalent that reads from a flat
             * positive-quadrature array at level k=nt. Implement directly. */
            const double *Im = &I_dn_BOA_per_m[m * n_mu];
            int j_hi = 1;
            while (j_hi < n_mu && mu_quad_pos[j_hi - 1] < mu_view) j_hi++;
            int j_lo = j_hi - 1;
            if (j_lo < 1) { j_lo = 1; j_hi = 2; }
            if (j_hi > n_mu) { j_hi = n_mu; j_lo = n_mu - 1; }
            const double mu_lo = mu_quad_pos[j_lo - 1];
            const double mu_hi = mu_quad_pos[j_hi - 1];
            const double t     = (mu_view - mu_lo) / (mu_hi - mu_lo);
            I_dn_per_m_view[m] = (1.0 - t) * Im[j_lo - 1] + t * Im[j_hi - 1];
        }
        const double I_dn_view = rt_solver_reconstruct_phi(I_dn_per_m_view, m_max, dphi);
        out->T_diff_dn_dir = (2.0 * pi * mu_view * I_dn_view) / (pi * mu_solar);

        /* Upward total at single-view: T_total_up_dir = ρ_TOA(μ_v,φ) · μ_sun. */
        out->T_total_up_dir = out->rho_I * mu_solar;
        /* Standalone atmosphere/surface path: the solver does not run a
         * companion black-surface decomposition.  Record the solved total as
         * atmospheric/surface path and expose the direct sunglint separately. */
        out->rho_atm_path_I = out->rho_I;
        out->rho_atm_path_Q = out->rho_Q;
        out->rho_atm_path_U = out->rho_U;
        out->T_dir_up_view = exp(-tau_tot / ((mu_view > 0.0) ? mu_view : 1.0));
        /* A solved atmosphere+surface case has no prescribed upward source at
         * BOA, so an upward atmospheric transmission cannot be inferred from
         * rho_TOA.  The former reciprocity shortcut T_diff_up=T_diff_dn was not
         * an RT solution and is deliberately removed. */
        out->T_diff_up_view = NAN;
        out->T_total_up_view = NAN;
        out->I_TOA_water_signal = NAN;
        out->T_up_rt_valid = 0;

        /* Direct sunglint upward at single-view: T_sg_up_dir = ρ_glint · μ_sun.
         * Cox-Munk only; flat-Fresnel direct sunglint is part of the SOS BC. */
        out->T_sg_up_dir = 0.0;
        if (cs->surface == RT_SURFACE_BLACK_FRESNEL_OCEAN && cs->wind_speed > 0.0) {
            const double phi_v = rt_raa_to_surface_view_phi(cs->raa_deg);
            const double phi_0 = 0.0;
            const double mu_v_use = (cs->vza_deg > 0.0) ? mu_view : 1.0;
            double rho_glint[3];
            surface_direct_sunglint_rho(mu_v_use, phi_v, mu_solar, phi_0,
                                         cs->wind_speed, cs->sigma_type,
                                         cs->n_water, tau_tot,
                                         cs->q_convention, rho_glint);
            if (ipss_ready && ipss_tau_dir_dn >= 0.0) {
                /* keep the exported diagnostic on the same footing as the
                 * glint actually folded into rho_I above */
                const double c_glint = exp(-(ipss_tau_dir_dn - tau_tot / mu_solar));
                rho_glint[0] *= c_glint; rho_glint[1] *= c_glint; rho_glint[2] *= c_glint;
            }
            out->T_sg_up_dir = rho_glint[0] * mu_solar;
            out->rho_glint_direct_I = rho_glint[0];
            out->rho_glint_direct_Q = rho_glint[1];
            out->rho_glint_direct_U = rho_glint[2];
        }
    }

    /* v1.01: post-SOS LUT grid reconstruction (Fourier sum at each
     * (vza, raa) grid point, using the per-m grid samples cached above).
     * Sun-glint direct addition applies identically per grid point. */
    if (lut_out) {
        const int n_vza_g = lut_out->n_vza;
        const int n_raa_g = lut_out->n_raa;
        const double tau_total_g = atm.h[atm.n_layers];
        const double ipss_cglint =
            (ipss_ready && ipss_tau_dir_dn >= 0.0)
                ? exp(-(ipss_tau_dir_dn - tau_total_g / mu_solar)) : 1.0;
        double I_pm[33], Q_pm[33], U_pm[33];

        /* Case-level scalars (identical to single-case path). */
        lut_out->T_dir_dn = out->T_dir_dn;
        lut_out->T_diff_dn_hemi = out->T_diff_dn_hemi;

        /* Pre-compute I_dn_per_m at each grid vza for directional T_diff_dn. */
        const double *mu_quad_pos = &atm.rm[+1];
        double *I_dn_per_m_grid_v = NULL;
        if (lut_out->T_diff_dn_dir) {
            I_dn_per_m_grid_v = (double*)calloc((size_t)(m_max + 1) * n_vza_g, sizeof(double));
            if (!I_dn_per_m_grid_v) { rc = -1; goto cleanup; }
            for (int m = 0; m <= m_max; ++m) {
                const double *Im = &I_dn_BOA_per_m[m * n_mu];
                for (int iv = 0; iv < n_vza_g; ++iv) {
                    const double mu_g = mu_v_grid[iv];
                    int j_hi = 1;
                    while (j_hi < n_mu && mu_quad_pos[j_hi - 1] < mu_g) j_hi++;
                    int j_lo = j_hi - 1;
                    if (j_lo < 1) { j_lo = 1; j_hi = 2; }
                    if (j_hi > n_mu) { j_hi = n_mu; j_lo = n_mu - 1; }
                    const double mu_lo = mu_quad_pos[j_lo - 1];
                    const double mu_hi = mu_quad_pos[j_hi - 1];
                    const double t = (mu_g - mu_lo) / (mu_hi - mu_lo);
                    I_dn_per_m_grid_v[m * n_vza_g + iv] =
                        (1.0 - t) * Im[j_lo - 1] + t * Im[j_hi - 1];
                }
            }
        }

        for (int iv = 0; iv < n_vza_g; ++iv) {
            for (int m = 0; m <= m_max; ++m) {
                I_pm[m] = I_per_m_grid[m * n_vza_g + iv];
                Q_pm[m] = Q_per_m_grid[m * n_vza_g + iv];
                U_pm[m] = U_per_m_grid[m * n_vza_g + iv];
            }
            const double mu_v_g = mu_v_grid[iv];

            /* Per-vza I_dn^m extraction for directional T_diff_dn. */
            double I_dn_pm[33] = {0};
            if (I_dn_per_m_grid_v) {
                for (int m = 0; m <= m_max; ++m) {
                    I_dn_pm[m] = I_dn_per_m_grid_v[m * n_vza_g + iv];
                }
            }

            for (int ir = 0; ir < n_raa_g; ++ir) {
                const double raa_g     = lut_out->raa_deg[ir];
                const double dphi_g    = rt_raa_to_atm_fourier_phi(raa_g);
                const double I_TOA_g   = rt_solver_reconstruct_phi    (I_pm, m_max, dphi_g);
                const double Q_TOA_g   = rt_solver_reconstruct_phi    (Q_pm, m_max, dphi_g);
                const double U_TOA_g   = rt_solver_reconstruct_phi_sin(U_pm, m_max, dphi_g);

                if (lut_out->raw_I) {   /* v1.10 S7: pre-rho, pre-glint */
                    const size_t idx_raw = (size_t)iv * (size_t)n_raa_g + (size_t)ir;
                    lut_out->raw_I[idx_raw] = I_TOA_g;
                    lut_out->raw_Q[idx_raw] = Q_TOA_g;
                    lut_out->raw_U[idx_raw] = U_TOA_g;
                }
                double rho_I_g = rt_solver_reflectance_from_intensity(I_TOA_g, mu_solar);
                double rho_Q_g = rt_solver_reflectance_from_intensity(Q_TOA_g, mu_solar);
                double rho_U_g = rt_solver_reflectance_from_intensity(U_TOA_g, mu_solar);

                /* v1.11 IPSS Eq. (7) per grid direction (see rt_ipss.h). */
                if (ipss_ready) {
                    double kappa_g = 1.0;
                    rt_ipss_kappa(&ipss_geo, lut_out->vza_deg[iv], raa_g,
                                  &kappa_g, NULL, NULL, NULL);
                    rho_I_g *= kappa_g; rho_Q_g *= kappa_g; rho_U_g *= kappa_g;
                    /* Diagnostic (env-gated, inert by default): dump kappa for
                     * IPSS and for the literature "PSS" convention (spherical
                     * line of sight, solar beam evaluated on the nadir column),
                     * which is what Zhai & Hu (2022) Sec. 3 reports errors for.
                     * OCRT's legacy Chapman PSSA is weaker still: it leaves the
                     * viewing path plane-parallel. */
                    const char *dmp = getenv("OCRT_IPSS_DUMP");
                    if (dmp && dmp[0]) {
                        const double i_sph = rt_ipss_single_scatter(
                            &ipss_geo, lut_out->vza_deg[iv], raa_g,
                            RT_IPSS_SUN_SPHERICAL);
                        const double i_nad = rt_ipss_single_scatter(
                            &ipss_geo, lut_out->vza_deg[iv], raa_g,
                            RT_IPSS_SUN_NADIR);
                        const double i_pp = rt_ipss_single_scatter_pp(
                            &ipss_geo, lut_out->vza_deg[iv]);
                        FILE *fd = fopen(dmp, "a");
                        if (fd) {
                            fprintf(fd, "%.6f,%.6f,%.12e,%.12e,%.12e,%.12e\n",
                                    lut_out->vza_deg[iv], raa_g, kappa_g,
                                    (i_pp > 0.0) ? i_nad / i_pp : 1.0,
                                    i_sph, i_pp);
                            fclose(fd);
                        }
                    }
                }

                /* Direct sunglint contribution + capture T_sg_up_dir. */
                double rho_sg_I_g = 0.0;
                if (cs->surface == RT_SURFACE_BLACK_FRESNEL_OCEAN && cs->wind_speed > 0.0) {
                    const double phi_v_g = rt_raa_to_surface_view_phi(raa_g);
                    const double phi_0   = 0.0;
                    const double mu_v_use = (lut_out->vza_deg[iv] > 0.0) ? mu_v_g : 1.0;
                    double rho_glint_g[3];
                    surface_direct_sunglint_rho(mu_v_use, phi_v_g, mu_solar, phi_0,
                                                 cs->wind_speed, cs->sigma_type,
                                                 cs->n_water, tau_total_g,
                                                 cs->q_convention, rho_glint_g);
                    rho_glint_g[0] *= ipss_cglint;
                    rho_glint_g[1] *= ipss_cglint;
                    rho_glint_g[2] *= ipss_cglint;
                    if (cs->decouple_sunglint == 0) {
                        rho_I_g += rho_glint_g[0];
                        rho_Q_g += rho_glint_g[1];
                        rho_U_g += -rho_glint_g[2];
                    }
                    rho_sg_I_g = rho_glint_g[0];
                }

                const int idx = iv * n_raa_g + ir;
                lut_out->rho_I[idx] = rho_I_g;
                lut_out->rho_Q[idx] = rho_Q_g;
                lut_out->rho_U[idx] = rho_U_g;

                /* Transmittances per direction. */
                if (lut_out->T_diff_dn_dir) {
                    const double I_dn_dir = rt_solver_reconstruct_phi(I_dn_pm, m_max, dphi_g);
                    lut_out->T_diff_dn_dir[idx] = (2.0 * pi * mu_v_g * I_dn_dir) / (pi * mu_solar);
                }
                if (lut_out->T_sg_up_dir) {
                    lut_out->T_sg_up_dir[idx] = rho_sg_I_g * mu_solar;
                }
                if (lut_out->T_total_up_dir) {
                    lut_out->T_total_up_dir[idx] = rho_I_g * mu_solar;
                }
            }
        }
        free(I_dn_per_m_grid_v);
    }

cleanup:
    rt_ipss_free(&ipss_geo);
    free(I_dn_BOA_per_m);
    free(I_per_m_grid); free(Q_per_m_grid); free(U_per_m_grid); free(mu_v_grid);
    free(total_u); free(total_q); free(total_i);
    free(prim_u);  free(prim_q);  free(prim_i);
    free(src_u);   free(src_q);   free(src_i);
    rt_legendre_workspace_free(&ws);
    rt_atm_free(&atm);
    return rc;
}

/* Public wrapper: Rayleigh-only path (legacy signature, bit-level Phase 2). */
int rt_solve_case_pol(const rt_case_t *cs, const rt_options_t *opts,
                      rt_result_t *out) {
    return rt_solve_case_pol_impl(cs, opts, NULL, out, NULL, NULL);
}

/* Public wrapper: aerosol-aware path (Phase 3). */
int rt_solve_case_pol_aerosol(const rt_case_t *cs, const rt_options_t *opts,
                               const rt_aerosol_input_t *aer,
                               rt_result_t *out) {
    return rt_solve_case_pol_impl(cs, opts, aer, out, NULL, NULL);
}

/* Public wrapper: LUT grid path (v1.01).  Solves SOS once and reconstructs
 * the radiation field at every (vza, raa) point in lut_out's grid arrays.
 * The single-point rt_result_t out is also filled, using cs->vza_deg /
 * cs->raa_deg as the representative direction (same as the non-LUT call). */
int rt_solve_case_pol_lut(const rt_case_t *cs, const rt_options_t *opts,
                          const rt_aerosol_input_t *aer,
                          rt_result_t *out,
                          rt_lut_grid_out_t *lut_out) {
    return rt_solve_case_pol_impl(cs, opts, aer, out, lut_out, NULL);
}

/* Phase B.4 Stage 2b: atmospheric SOS with BOA downward Stokes field export.
 * Same logic as rt_solve_case_pol_aerosol (aer may be NULL for Rayleigh-only),
 * but additionally captures the diffuse BOA downward I^m/Q^m/U^m field into
 * boa_export. boa_export->mu_quad_pos, I/Q/U_per_m arrays must be pre-
 * allocated with size [(m_max+1) × n_mu] matching opts->fourier_m_max and
 * opts->n_mu. */
int rt_solve_case_pol_for_ocean(const rt_case_t *cs, const rt_options_t *opts,
                                 const rt_aerosol_input_t *aer,
                                 rt_result_t *out,
                                 rt_atm_boa_export_t *boa_export) {
    return rt_solve_case_pol_impl(cs, opts, aer, out, NULL, boa_export);
}

/* BOA export struct alloc/free helpers (Phase B.4 Stage 2b). */
int rt_atm_boa_export_alloc(rt_atm_boa_export_t *be, int m_max, int n_mu) {
    if (!be || m_max < 0 || n_mu < 1) return -1;
    memset(be, 0, sizeof(*be));
    be->m_max = m_max;
    be->n_mu  = n_mu;
    size_t per_m_size = (size_t)(m_max + 1) * (size_t)n_mu;
    be->mu_quad_pos = (double*)calloc((size_t)n_mu, sizeof(double));
    be->I_per_m     = (double*)calloc(per_m_size, sizeof(double));
    be->Q_per_m     = (double*)calloc(per_m_size, sizeof(double));
    be->U_per_m     = (double*)calloc(per_m_size, sizeof(double));
    if (!be->mu_quad_pos || !be->I_per_m || !be->Q_per_m || !be->U_per_m) {
        rt_atm_boa_export_free(be);
        return -1;
    }
    be->allocated = 1;
    return 0;
}

void rt_atm_boa_export_free(rt_atm_boa_export_t *be) {
    if (!be) return;
    if (be->mu_quad_pos) free(be->mu_quad_pos);
    if (be->I_per_m)     free(be->I_per_m);
    if (be->Q_per_m)     free(be->Q_per_m);
    if (be->U_per_m)     free(be->U_per_m);
    memset(be, 0, sizeof(*be));
}

/* Diagnostic-only water/ocean order trace. When OCRT_DUMP_WATER_ORDER_TRACE
 * is set, append the per-Fourier-mode scattering-order contribution at the
 * top boundary (k=0) for all signed ordinates. This is side-effect free and
 * is intended only for OSOAA/OCRT order-budget audits. */
static void dump_water_order_trace_if_requested(const char *solver_name,
                                                int m, int order,
                                                const rt_atm_t *atm,
                                                const double *I,
                                                const double *Q,
                                                const double *U) {
    const char *path = getenv("OCRT_DUMP_WATER_ORDER_TRACE");
    if (!path || !path[0] || !atm || !I || !Q || !U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    static int wrote_header = 0;
    if (!wrote_header) {
        fprintf(fp, "solver,m,order,k,j_signed,mu,I,Q,U\n");
        wrote_header = 1;
    }
    for (int j = -n_mu; j <= n_mu; ++j) {
        if (j == 0) continue;
        size_t idx = (size_t)0 * (size_t)dirs + (size_t)(j + n_mu);
        fprintf(fp, "%s,%d,%d,0,%d,%.17g,%.17e,%.17e,%.17e\n",
                solver_name ? solver_name : "sos", m, order, j, atm->rm[j],
                I[idx], Q[idx], U[idx]);
    }
    fclose(fp);
}


/* Diagnostic-only decomposition of internally-reflected ocean SOS orders.
 * When OCRT_DUMP_WATER_ORDER_COMPONENT_TRACE is set, append top-boundary
 * fields split by origin:
 *   source_direct  : volume source driven by the direct-source lineage;
 *   source_feedback: volume source driven by previously internally-reflected lineage;
 *   bc_feedback    : current-order downward field injected by Rww top BC;
 *   total_recon    : source_direct + source_feedback + bc_feedback.
 * This is a bookkeeping audit only; it does not alter production radiance. */
static void dump_water_order_component_trace_if_requested(const char *solver_name,
                                                          const char *component,
                                                          int m, int order,
                                                          const rt_atm_t *atm,
                                                          const double *I,
                                                          const double *Q,
                                                          const double *U) {
    const char *path = getenv("OCRT_DUMP_WATER_ORDER_COMPONENT_TRACE");
    if (!path || !path[0] || !atm || !I || !Q || !U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    if (file_len == 0L) {
        fprintf(fp, "solver,component,m,order,k,j_signed,mu,I,Q,U\n");
    }
    for (int j = -n_mu; j <= n_mu; ++j) {
        if (j == 0) continue;
        size_t idx = (size_t)0 * (size_t)dirs + (size_t)(j + n_mu);
        fprintf(fp, "%s,%s,%d,%d,0,%d,%.17g,%.17e,%.17e,%.17e\n",
                solver_name ? solver_name : "sos",
                component ? component : "component",
                m, order, j, atm->rm[j], I[idx], Q[idx], U[idx]);
    }
    fclose(fp);
}


/* Diagnostic-only split of the transported order field into the field obtained
 * by propagating the volume source alone and the residual.  In the flat
 * water-side internal-reflection solver, the residual is the same-order top
 * boundary/internal-reflection feedback contribution.  For the standard delayed
 * top-BC formulation the residual is expected to live in the downwelling branch
 * and to be nearly zero in the top-up numerator; this diagnostic verifies that
 * expectation explicitly.
 * Enabled by OCRT_DUMP_WATER_ORDER_COMPONENT_TRACE=/path/to.csv. */
static void dump_water_order_source_feedback_trace_if_requested(const char *solver_name,
                                                                int m, int order,
                                                                const rt_atm_t *atm,
                                                                const double *J_I,
                                                                const double *J_Q,
                                                                const double *J_U,
                                                                const double *T_I,
                                                                const double *T_Q,
                                                                const double *T_U) {
    const char *path = getenv("OCRT_DUMP_WATER_ORDER_COMPONENT_TRACE");
    if (!path || !path[0] || !atm || !J_I || !J_Q || !J_U || !T_I || !T_Q || !T_U) return;
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;

    double *S_I = (double *)calloc(buf_size, sizeof(double));
    double *S_Q = (double *)calloc(buf_size, sizeof(double));
    double *S_U = (double *)calloc(buf_size, sizeof(double));
    if (!S_I || !S_Q || !S_U) { free(S_I); free(S_Q); free(S_U); return; }

    if (rt_solver_integrate_bcs(atm, m, J_I, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, S_I) != 0 ||
        rt_solver_integrate_bcs(atm, m, J_Q, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, S_Q) != 0 ||
        rt_solver_integrate_bcs(atm, m, J_U, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, S_U) != 0) {
        free(S_I); free(S_Q); free(S_U); return;
    }

    FILE *fp = fopen(path, "a+");
    if (!fp) { free(S_I); free(S_Q); free(S_U); return; }
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    if (file_len == 0L) {
        fprintf(fp, "solver,component,m,order,k,j_signed,mu,I,Q,U\n");
    }
    for (int j = -n_mu; j <= n_mu; ++j) {
        if (j == 0) continue;
        const size_t idx = (size_t)0 * (size_t)dirs + (size_t)(j + n_mu);
        const double FI = T_I[idx] - S_I[idx];
        const double FQ = T_Q[idx] - S_Q[idx];
        const double FU = T_U[idx] - S_U[idx];
        fprintf(fp, "%s,source_only,%d,%d,0,%d,%.17g,%.17e,%.17e,%.17e\n",
                solver_name ? solver_name : "sos", m, order, j, atm->rm[j], S_I[idx], S_Q[idx], S_U[idx]);
        fprintf(fp, "%s,feedback_residual,%d,%d,0,%d,%.17g,%.17e,%.17e,%.17e\n",
                solver_name ? solver_name : "sos", m, order, j, atm->rm[j], FI, FQ, FU);
        fprintf(fp, "%s,total,%d,%d,0,%d,%.17g,%.17e,%.17e,%.17e\n",
                solver_name ? solver_name : "sos", m, order, j, atm->rm[j], T_I[idx], T_Q[idx], T_U[idx]);
    }
    fclose(fp);
    free(S_I); free(S_Q); free(S_U);
}


/* Diagnostic-only water/ocean source trace. When OCRT_DUMP_WATER_SOURCE_TRACE
 * is set, append the Fourier-mode source J_n^m(tau_k, mu_j) before vertical
 * transport.  This is intended to be compared with OSOAA's order-by-order
 * sea-source dumps; it does not affect production physics. */
static void dump_water_source_trace_if_requested(const char *solver_name,
                                                 const char *source_stage,
                                                 int m, int order,
                                                 const rt_atm_t *atm,
                                                 const double *J_I,
                                                 const double *J_Q,
                                                 const double *J_U) {
    const char *path = getenv("OCRT_DUMP_WATER_SOURCE_TRACE");
    if (!path || !path[0] || !atm || !J_I || !J_Q || !J_U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    if (file_len == 0L) {
        fprintf(fp, "solver,stage,m,order,k,j_signed,mu,J_I,J_Q,J_U\n");
    }
    for (int ksel_i = 0; ksel_i < 3; ++ksel_i) {
        int k = (ksel_i == 0) ? 0 : ((ksel_i == 1) ? (nt / 2) : nt);
        if (k < 0) k = 0;
        if (k > nt) k = nt;
        for (int j = -n_mu; j <= n_mu; ++j) {
            if (j == 0) continue;
            size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
            fprintf(fp, "%s,%s,%d,%d,%d,%d,%.17g,%.17e,%.17e,%.17e\n",
                    solver_name ? solver_name : "sos",
                    source_stage ? source_stage : "source",
                    m, order, k, j, atm->rm[j], J_I[idx], J_Q[idx], J_U[idx]);
        }
    }
    fclose(fp);
}




/* Step22 diagnostic-only feedback-stage transport trace. This records the
 * internal-reflection feedback chain explicitly:
 *   B_rww_downfield : downwelling field generated by Rww × previous upwelling,
 *                     then transported through the column with zero volume source;
 *   C_volume_from_B : volume-source propagated field generated only from that
 *                     reflected downfield.
 * It does not alter the solver state. */
static void dump_water_feedback_stage_trace_if_requested(const char *solver_name,
                                                         const char *stage,
                                                         int m, int order,
                                                         const rt_atm_t *atm,
                                                         const double *I,
                                                         const double *Q,
                                                         const double *U) {
    const char *path = getenv("OCRT_DUMP_WATER_FEEDBACK_STAGE_TRACE");
    if (!path || !path[0] || !atm || !I || !Q || !U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    if (len == 0L) {
        fprintf(fp, "solver,stage,wavelength_nm,m,order,k,j_signed,mu,I,Q,U\n");
    }
    double wl = NAN;
    const char *wls = getenv("OCRT_TRACE_WAVELENGTH_NM");
    if (wls && wls[0]) wl = atof(wls);
    for (int kk_sel = 0; kk_sel < 3; ++kk_sel) {
        int k = (kk_sel == 0) ? 0 : ((kk_sel == 1) ? nt / 2 : nt);
        if (k < 0) k = 0;
        if (k > nt) k = nt;
        for (int j = -n_mu; j <= n_mu; ++j) {
            if (j == 0) continue;
            const size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
            fprintf(fp, "%s,%s,%.17g,%d,%d,%d,%d,%.17g,%.17e,%.17e,%.17e\n",
                    solver_name ? solver_name : "water",
                    stage ? stage : "stage", wl, m, order, k, j, atm->rm[j],
                    I[idx], Q[idx], U[idx]);
        }
    }
    fclose(fp);
}

static void dump_water_feedback_chain_if_requested(const char *solver_name,
                                                   int m, int order,
                                                   const rt_atm_t *atm,
                                                   const rt_legendre_workspace_t *ws,
                                                   const double *zero_src,
                                                   const double *top_bc_I,
                                                   const double *top_bc_Q,
                                                   const double *top_bc_U) {
    const char *path = getenv("OCRT_DUMP_WATER_FEEDBACK_STAGE_TRACE");
    if (!path || !path[0] || !atm || !ws || !zero_src || !top_bc_I || !top_bc_Q || !top_bc_U) return;
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;
    double *B_I = (double*)calloc(buf_size, sizeof(double));
    double *B_Q = (double*)calloc(buf_size, sizeof(double));
    double *B_U = (double*)calloc(buf_size, sizeof(double));
    double *JB_I = (double*)calloc(buf_size, sizeof(double));
    double *JB_Q = (double*)calloc(buf_size, sizeof(double));
    double *JB_U = (double*)calloc(buf_size, sizeof(double));
    double *C_I = (double*)calloc(buf_size, sizeof(double));
    double *C_Q = (double*)calloc(buf_size, sizeof(double));
    double *C_U = (double*)calloc(buf_size, sizeof(double));
    if (!B_I || !B_Q || !B_U || !JB_I || !JB_Q || !JB_U || !C_I || !C_Q || !C_U) {
        free(B_I); free(B_Q); free(B_U); free(JB_I); free(JB_Q); free(JB_U); free(C_I); free(C_Q); free(C_U); return;
    }
    if (rt_solver_integrate_bcs(atm, m, zero_src, NULL, top_bc_I,
                                RT_INTEGRATION_METHOD_LINEAR, B_I) != 0 ||
        rt_solver_integrate_bcs(atm, m, zero_src, NULL, top_bc_Q,
                                RT_INTEGRATION_METHOD_LINEAR, B_Q) != 0 ||
        rt_solver_integrate_bcs(atm, m, zero_src, NULL, top_bc_U,
                                RT_INTEGRATION_METHOD_LINEAR, B_U) != 0) {
        free(B_I); free(B_Q); free(B_U); free(JB_I); free(JB_Q); free(JB_U); free(C_I); free(C_Q); free(C_U); return;
    }
    dump_water_feedback_stage_trace_if_requested(solver_name, "B_rww_downfield", m, order, atm, B_I, B_Q, B_U);
    if (rt_sos_operator_apply_vector(atm, m, ws, B_I, B_Q, B_U, JB_I, JB_Q, JB_U) == 0 &&
        rt_solver_integrate_bcs(atm, m, JB_I, NULL, NULL,
                                RT_INTEGRATION_METHOD_LINEAR, C_I) == 0 &&
        rt_solver_integrate_bcs(atm, m, JB_Q, NULL, NULL,
                                RT_INTEGRATION_METHOD_LINEAR, C_Q) == 0 &&
        rt_solver_integrate_bcs(atm, m, JB_U, NULL, NULL,
                                RT_INTEGRATION_METHOD_LINEAR, C_U) == 0) {
        dump_water_feedback_stage_trace_if_requested(solver_name, "C_volume_from_B", m, order, atm, C_I, C_Q, C_U);
    }
    free(B_I); free(B_Q); free(B_U); free(JB_I); free(JB_Q); free(JB_U); free(C_I); free(C_Q); free(C_U);
}

/* ----------------------------------------------------------------------------
 * SOS multiple-scattering iteration
 * --------------------------------------------------------------------------*/

/* Build the m-th SOS source J_n^m(τ_k, μ_j) from the previous-order
 * field I_{n-1}.  Per OS.f L422-468 (with m≤2 / m>2 branches unified —
 * atm.xpl[j] is automatically zero for l<m so the Rayleigh term
 * vanishes naturally).
 *
 * For each (kk_layer, k_dir > 0), the loop sums over j' > 0 only and
 * folds in the negative-j' contribution via the explicit
 *   I[+j'] · phase(+j', +k)  +  I[-j'] · phase(+j', -k)
 * structure (matches 6SV's storage convention; phase_fourier_m stored
 * with first index ≥ 0).
 *
 * Output J_curr has the same shape as I_prev; the j=0 (solar slot)
 * column is set to NaN (defensive — same convention as Step 5/6).
 */

/* Diagnostic-only water/ocean source prefactor trace.  This complements
 * OCRT_DUMP_WATER_SOURCE_TRACE by recording the scalar/source factors before
 * vertical transport.  It is intentionally environment-variable gated and has
 * no production side effects.
 *
 * CSV columns:
 * solver,stage,m,order,k,j_signed,mu,ch,xdel,ydel,beam_q,phase_I,phase_Q,phase_U,J_I,J_Q,J_U,wavelength_nm
 *
 * For primary_source, phase_* are the direct-beam Fourier coefficients that
 * multiply ch*xdel (hydrosol-only) plus the beam-Q terms where present.  For
 * ms_source, phase_* are reported as an effective source prefactor J/xdel
 * (not a pure phase coefficient) because the MS source is already a previous-
 * order radiance weighted integral over incoming directions.
 */
static void dump_water_source_pref_trace_if_requested(const char *solver_name,
                                                      const char *stage,
                                                      int m, int order,
                                                      const rt_atm_t *atm,
                                                      const double *phase_I,
                                                      const double *phase_Q,
                                                      const double *phase_U,
                                                      const double *J_I,
                                                      const double *J_Q,
                                                      const double *J_U) {
    const char *path = getenv("OCRT_DUMP_WATER_SOURCE_PREF_TRACE");
    if (!path || !path[0] || !atm || !J_I || !J_Q || !J_U) return;
    FILE *fp = fopen(path, "a+");
    if (!fp) return;
    fseek(fp, 0, SEEK_END);
    long file_len = ftell(fp);
    if (file_len == 0L) {
        fprintf(fp, "solver,stage,m,order,k,j_signed,mu,ch,xdel,ydel,beam_q,phase_I,phase_Q,phase_U,J_I,J_Q,J_U,wavelength_nm\n");
    }
    double wl = NAN;
    const char *wls = getenv("OCRT_TRACE_WAVELENGTH_NM");
    if (wls && wls[0]) wl = atof(wls);
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    for (int ksel_i = 0; ksel_i < 3; ++ksel_i) {
        int kk = (ksel_i == 0) ? 0 : ((ksel_i == 1) ? (nt / 2) : nt);
        if (kk < 0) kk = 0;
        if (kk > nt) kk = nt;
        const double ch = atm->ch ? atm->ch[kk] : NAN;
        const double x  = atm->xdel ? atm->xdel[kk] : NAN;
        const double y  = atm->ydel ? atm->ydel[kk] : NAN;
        for (int j = -n_mu; j <= n_mu; ++j) {
            if (j == 0) continue;
            const size_t idx = (size_t)kk * (size_t)dirs + (size_t)(j + n_mu);
            const double pI = phase_I ? phase_I[j + n_mu] : ((fabs(x) > 1e-300) ? J_I[idx] / x : NAN);
            const double pQ = phase_Q ? phase_Q[j + n_mu] : ((fabs(x) > 1e-300) ? J_Q[idx] / x : NAN);
            const double pU = phase_U ? phase_U[j + n_mu] : ((fabs(x) > 1e-300) ? J_U[idx] / x : NAN);
            fprintf(fp, "%s,%s,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17g\n",
                    solver_name ? solver_name : "water",
                    stage ? stage : "source_pref",
                    m, order, kk, j, atm->rm[j], ch, x, y, atm->beam_q,
                    pI, pQ, pU, J_I[idx], J_Q[idx], J_U[idx], wl);
        }
    }
    fclose(fp);
}



/* ‖field‖_∞ over j_signed ≠ 0 (skips solar slot). */


static double field_max_abs(const double *field, int n_layers, int n_mu) {
    const int dirs = 2 * n_mu + 1;
    double mx = 0.0;
    for (int k = 0; k <= n_layers; k++) {
        for (int j = -n_mu; j <= n_mu; j++) {
            if (j == 0) continue;
            double v = field[(size_t)k * (size_t)dirs + (size_t)(j + n_mu)];
            double a = fabs(v);
            if (a > mx) mx = a;
        }
    }
    return mx;
}

/* Successive-orders iteration (scalar I; the _pol variant mirrors this for
 * {I,Q,U} with the polarized kernels):
 *
 *   I_total = I^(1) + sum_{n>=2} I^(n)
 *   I^(n)   = FormalSolution( J^(n) ),  J^(n) = Kernel . I^(n-1)
 *
 * Convergence criterion (relative sup-norm of the latest order):
 *   residual = max|I^(n)| / max|I_total|  <  tolerance
 * i.e. geometric-tail truncation of the Neumann series; the number of
 * orders needed grows with the effective single-scattering albedo (the
 * printed `orders=` value; high-omega Csed50 cases run 400+ orders).
 * Physical validity requires spectral radius of the transport operator
 * < 1, guaranteed for omega < 1 (some absorption) or a lossy boundary.
 * Integration inside the loop is hard-coded LINEAR (the CONSTANT method
 * would compound its O(dtau) error across hundreds of orders).
 * S-003 RESOLVED (v1.09 commit #4): RT_SOS_ACCELERATION_GEOMETRIC was
 * accepted but never implemented (silent no-op).  The CLI now REJECTS
 * "geometric"; the enum value is retained for ABI stability and this
 * validation still tolerates it for direct-API callers, but nothing
 * implements it — do not add new users. */
int rt_solver_sos(const rt_atm_t *atm, int m,
                  const rt_legendre_workspace_t *ws,
                  const double *primary_field,
                  const rt_solver_sos_options_t *opts,
                  double *total_field,
                  double *I_per_order,
                  rt_solver_sos_result_t *result) {
    if (!atm || !ws || !primary_field || !total_field || !result) return -1;
    if (m < 0)                                                     return -1;
    if (atm->n_mu != ws->n_mu)                                     return -1;

    rt_solver_sos_options_t opt_default = rt_solver_sos_options_default();
    if (!opts) opts = &opt_default;
    if (opts->max_iterations < 1)                                  return -1;
    if (!(opts->tolerance > 0.0))                                  return -1;
    if (opts->acceleration != RT_SOS_ACCELERATION_PLAIN &&
        opts->acceleration != RT_SOS_ACCELERATION_GEOMETRIC)       return -1;

    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;

    /* total_field starts as a copy of primary_field. */
    memcpy(total_field, primary_field, buf_size * sizeof(double));

    /* Quick-out: max_iterations == 1 means primary only. */
    if (opts->max_iterations == 1) {
        result->n_orders_used  = 1;
        result->converged      = 1;
        result->final_residual = 0.0;
        if (I_per_order)
            memcpy(I_per_order, primary_field, buf_size * sizeof(double));
        return 0;
    }

    double *I_prev = (double *)calloc(buf_size, sizeof(double));
    double *J_curr = (double *)calloc(buf_size, sizeof(double));
    double *I_curr = (double *)calloc(buf_size, sizeof(double));
    if (!I_prev || !J_curr || !I_curr) {
        free(I_prev); free(J_curr); free(I_curr);
        return -1;
    }
    memcpy(I_prev, primary_field, buf_size * sizeof(double));

    if (I_per_order)
        memcpy(I_per_order, primary_field, buf_size * sizeof(double));

    int n = 1;
    double residual = 0.0;
    int converged = 0;

    for (int it = 2; it <= opts->max_iterations; it++) {
        int rc = rt_sos_operator_apply_scalar(atm, m, ws, I_prev, J_curr);
        if (rc != 0) {
            free(I_prev); free(J_curr); free(I_curr);
            return rc;
        }
        /* LINEAR mode hard-coded — CONSTANT is for v1 primary-only demo
         * and would compound errors across SOS iterations. */
        if (rt_solver_integrate(atm, m, J_curr,
                                RT_INTEGRATION_METHOD_LINEAR, I_curr) != 0) {
            free(I_prev); free(J_curr); free(I_curr);
            return -2;
        }
        for (int k = 0; k <= nt; k++) {
            for (int j = -n_mu; j <= n_mu; j++) {
                if (j == 0) continue;
                size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
                total_field[idx] += I_curr[idx];
            }
        }
        n = it;

        if (I_per_order) {
            memcpy(I_per_order + (size_t)(it - 1) * buf_size,
                   I_curr, buf_size * sizeof(double));
        }

        double max_in    = field_max_abs(I_curr, nt, n_mu);
        double max_total = field_max_abs(total_field, nt, n_mu);
        residual = (max_total > 0.0) ? max_in / max_total : 0.0;
        if (residual < opts->tolerance) {
            converged = 1;
            break;
        }

        memcpy(I_prev, I_curr, buf_size * sizeof(double));
    }

    free(I_prev); free(J_curr); free(I_curr);

    result->n_orders_used  = n;
    result->converged      = converged;
    result->final_residual = residual;
    return 0;
}

int rt_solver_sos_pol(const rt_atm_t *atm, int m,
                      const rt_legendre_workspace_t *ws,
                      const double *primary_I,
                      const double *primary_Q,
                      const double *primary_U,
                      const rt_solver_sos_options_t *opts,
                      double *total_I,
                      double *total_Q,
                      double *total_U,
                      rt_solver_sos_result_t *result) {
    if (!atm || !ws || !primary_I || !primary_Q || !primary_U ||
        !total_I || !total_Q || !total_U || !result) return -1;
    if (m < 0)                                                     return -1;
    if (atm->n_mu != ws->n_mu)                                     return -1;

    rt_solver_sos_options_t opt_default = rt_solver_sos_options_default();
    if (!opts) opts = &opt_default;
    if (opts->max_iterations < 1)                                  return -1;
    if (!(opts->tolerance > 0.0))                                  return -1;

    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;

    /* total_X starts as a copy of primary_X. */
    memcpy(total_I, primary_I, buf_size * sizeof(double));
    memcpy(total_Q, primary_Q, buf_size * sizeof(double));
    memcpy(total_U, primary_U, buf_size * sizeof(double));

    dump_water_order_trace_if_requested("plain", m, 1, atm, primary_I, primary_Q, primary_U);

    /* Quick-out: max_iterations == 1 means primary only. */
    if (opts->max_iterations == 1) {
        result->n_orders_used  = 1;
        result->converged      = 1;
        result->final_residual = 0.0;
        return 0;
    }

    double *I_prev = calloc(buf_size, sizeof(double));
    double *Q_prev = calloc(buf_size, sizeof(double));
    double *U_prev = calloc(buf_size, sizeof(double));
    double *J_I    = calloc(buf_size, sizeof(double));
    double *J_Q    = calloc(buf_size, sizeof(double));
    double *J_U    = calloc(buf_size, sizeof(double));
    double *I_curr = calloc(buf_size, sizeof(double));
    double *Q_curr = calloc(buf_size, sizeof(double));
    double *U_curr = calloc(buf_size, sizeof(double));
    if (!I_prev || !Q_prev || !U_prev || !J_I || !J_Q || !J_U ||
        !I_curr || !Q_curr || !U_curr) {
        free(I_prev); free(Q_prev); free(U_prev);
        free(J_I); free(J_Q); free(J_U);
        free(I_curr); free(Q_curr); free(U_curr);
        return -1;
    }
    memcpy(I_prev, primary_I, buf_size * sizeof(double));
    memcpy(Q_prev, primary_Q, buf_size * sizeof(double));
    memcpy(U_prev, primary_U, buf_size * sizeof(double));

    int n = 1;
    double residual = 0.0;
    int converged = 0;

    for (int it = 2; it <= opts->max_iterations; it++) {
        int rc = rt_sos_operator_apply_vector(atm, m, ws, I_prev, Q_prev, U_prev,
                                       J_I, J_Q, J_U);
        dump_water_source_trace_if_requested("plain", "ms_source", m, it,
                                             atm, J_I, J_Q, J_U);
        dump_water_source_pref_trace_if_requested("plain", "ms_source", m, it,
                                                  atm, NULL, NULL, NULL,
                                                  J_I, J_Q, J_U);
        if (rc != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr);
            return rc;
        }
        if (rt_solver_integrate(atm, m, J_I,
                                RT_INTEGRATION_METHOD_LINEAR, I_curr) != 0 ||
            rt_solver_integrate(atm, m, J_Q,
                                RT_INTEGRATION_METHOD_LINEAR, Q_curr) != 0 ||
            rt_solver_integrate(atm, m, J_U,
                                RT_INTEGRATION_METHOD_LINEAR, U_curr) != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr);
            return -2;
        }

        /* Step87 diagnostic-only isolated previous-Stokes propagation dump.
         * This answers a narrower question than Step82 source splitting:
         * if only the previous-order I, Q, or U field is allowed to build
         * the next MS source, how much source and propagated 0- field does
         * it generate?  No production output is changed. */
        {
            const char *cp = getenv("OCRT_DUMP_WATER_MS_CHANNEL_PROP_TRACE");
            if (cp && cp[0]) {
                int ok_m = 1;
                const char *sm = getenv("OCRT_MS_CHANNEL_PROP_M");
                if (sm && sm[0] && m != atoi(sm)) ok_m = 0;
                int ok_o = 1;
                const char *so = getenv("OCRT_MS_CHANNEL_PROP_ORDER");
                if (so && so[0] && it != atoi(so)) ok_o = 0;
                if (ok_m && ok_o) {
                    double *Z  = calloc(buf_size, sizeof(double));
                    double *cJI = calloc(buf_size, sizeof(double));
                    double *cJQ = calloc(buf_size, sizeof(double));
                    double *cJU = calloc(buf_size, sizeof(double));
                    double *cFI = calloc(buf_size, sizeof(double));
                    double *cFQ = calloc(buf_size, sizeof(double));
                    double *cFU = calloc(buf_size, sizeof(double));
                    if (Z && cJI && cJQ && cJU && cFI && cFQ && cFU) {
                        FILE *pf = fopen(cp, "a+");
                        if (pf) {
                            fseek(pf, 0, SEEK_END);
                            if (ftell(pf) == 0L) {
                                fprintf(pf,
                                        "solver,wavelength_nm,m,order,prev_channel,out_component,source_abs,source_sum,field_all_abs,field_all_sum,top_up_abs,top_up_sum,bottom_down_abs,bottom_down_sum,source_to_top_abs_eff,source_to_all_abs_eff\n");
                            }
                            const char *pnames[3] = {"prev_I", "prev_Q", "prev_U"};
                            for (int pc = 0; pc < 3; ++pc) {
                                const double *inI = (pc == 0) ? I_prev : Z;
                                const double *inQ = (pc == 1) ? Q_prev : Z;
                                const double *inU = (pc == 2) ? U_prev : Z;
                                for (size_t bi = 0; bi < buf_size; ++bi) {
                                    cJI[bi] = cJQ[bi] = cJU[bi] = 0.0;
                                    cFI[bi] = cFQ[bi] = cFU[bi] = 0.0;
                                }
                                int rc_iso = rt_sos_operator_apply_vector(atm, m, ws, inI, inQ, inU, cJI, cJQ, cJU);
                                if (rc_iso == 0) {
                                    int ri = rt_solver_integrate(atm, m, cJI, RT_INTEGRATION_METHOD_LINEAR, cFI);
                                    int rq = rt_solver_integrate(atm, m, cJQ, RT_INTEGRATION_METHOD_LINEAR, cFQ);
                                    int ru = rt_solver_integrate(atm, m, cJU, RT_INTEGRATION_METHOD_LINEAR, cFU);
                                    if (ri == 0 && rq == 0 && ru == 0) {
                                        const double *srcs[3] = {cJI, cJQ, cJU};
                                        const double *flds[3] = {cFI, cFQ, cFU};
                                        const char *cnames[3] = {"I", "Q", "U"};
                                        for (int cc = 0; cc < 3; ++cc) {
                                            double s_abs = 0.0, s_sum = 0.0;
                                            double f_abs = 0.0, f_sum = 0.0;
                                            double tu_abs = 0.0, tu_sum = 0.0;
                                            double bd_abs = 0.0, bd_sum = 0.0;
                                            for (int kk2 = 0; kk2 <= nt; ++kk2) {
                                                for (int jj = -n_mu; jj <= n_mu; ++jj) {
                                                    if (jj == 0) continue;
                                                    size_t idx = (size_t)kk2 * (size_t)dirs + (size_t)(jj + n_mu);
                                                    s_abs += fabs(srcs[cc][idx]); s_sum += srcs[cc][idx];
                                                    f_abs += fabs(flds[cc][idx]); f_sum += flds[cc][idx];
                                                }
                                            }
                                            for (int jj = 1; jj <= n_mu; ++jj) {
                                                size_t top_idx = (size_t)0 * (size_t)dirs + (size_t)(jj + n_mu);
                                                tu_abs += fabs(flds[cc][top_idx]); tu_sum += flds[cc][top_idx];
                                                size_t bot_idx = (size_t)nt * (size_t)dirs + (size_t)(-jj + n_mu);
                                                bd_abs += fabs(flds[cc][bot_idx]); bd_sum += flds[cc][bot_idx];
                                            }
                                            const double eff_top = (s_abs > 0.0) ? tu_abs / s_abs : 0.0;
                                            const double eff_all = (s_abs > 0.0) ? f_abs / s_abs : 0.0;
                                            double wl = 0.0;
                                            const char *wls = getenv("OCRT_TRACE_WAVELENGTH_NM");
                                            if (wls && wls[0]) wl = atof(wls);
                                            fprintf(pf,
                                                    "OCRT,%.17g,%d,%d,%s,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                                                    wl, m, it, pnames[pc], cnames[cc],
                                                    s_abs, s_sum, f_abs, f_sum, tu_abs, tu_sum, bd_abs, bd_sum, eff_top, eff_all);
                                        }
                                    }
                                }
                            }
                            fclose(pf);
                        }
                    }
                    free(Z); free(cJI); free(cJQ); free(cJU); free(cFI); free(cFQ); free(cFU);
                }
            }
        }

        dump_water_order_trace_if_requested("plain", m, it, atm, I_curr, Q_curr, U_curr);

        for (int k = 0; k <= nt; k++) {
            for (int j = -n_mu; j <= n_mu; j++) {
                if (j == 0) continue;
                size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
                total_I[idx] += I_curr[idx];
                total_Q[idx] += Q_curr[idx];
                total_U[idx] += U_curr[idx];
            }
        }
        n = it;

        /* Convergence: max of (max-of-each-component) over all (k, j) ≠ solar. */
        double max_in_I = field_max_abs(I_curr, nt, n_mu);
        double max_in_Q = field_max_abs(Q_curr, nt, n_mu);
        double max_in_U = field_max_abs(U_curr, nt, n_mu);
        double max_in   = max_in_I;
        if (max_in_Q > max_in) max_in = max_in_Q;
        if (max_in_U > max_in) max_in = max_in_U;

        double max_t_I = field_max_abs(total_I, nt, n_mu);
        double max_t_Q = field_max_abs(total_Q, nt, n_mu);
        double max_t_U = field_max_abs(total_U, nt, n_mu);
        double max_t   = max_t_I;
        if (max_t_Q > max_t) max_t = max_t_Q;
        if (max_t_U > max_t) max_t = max_t_U;

        residual = (max_t > 0.0) ? max_in / max_t : 0.0;
        if (residual < opts->tolerance) {
            converged = 1;
            break;
        }

        memcpy(I_prev, I_curr, buf_size * sizeof(double));
        memcpy(Q_prev, Q_curr, buf_size * sizeof(double));
        memcpy(U_prev, U_curr, buf_size * sizeof(double));
    }

    free(I_prev); free(Q_prev); free(U_prev);
    free(J_I); free(J_Q); free(J_U);
    free(I_curr); free(Q_curr); free(U_curr);

    result->n_orders_used  = n;
    result->converged      = converged;
    result->final_residual = residual;
    return 0;
}


/* Order-complete water-side flat internal-reflection down-field injection.
 *
 * Per positive node mu_j (jp = 1..n_mu), the top-level (k=0) upwelling
 * Stokes vector of the CURRENT order is specularly reflected back down:
 *
 *   [I,Q,U]_dn^m(k=0, -mu_j) = M_Rww(mu_j) . [I,Q,U]_up^m(k=0, +mu_j)
 *
 *   M_Rww(mu_j) : 3x3 water-side Fresnel internal-reflection Mueller for
 *       node jp, row-major at Rww_M[(jp-1)*9]; TIR nodes carry the identity
 *       (total reflection).  Flat surface => specular: node +j maps to -j
 *       with no cross-node coupling.
 *   No (-1)^m factor is applied: water-side specular reflection flips only
 *       the vertical propagation component and preserves the physical azimuth.
 *
 * The reflected field is then propagated source-free through the stack
 * (rt_solver_integrate_bcs with zero_src and top_down_bc = reflected
 * values) and ADDED to the order's field.  Summed over SOS orders this
 * converges to the self-consistent reflecting upper boundary WITHOUT
 * counting the reflection as a volume scattering event (it enters the
 * next order only through the field it adds, hence "order-complete").
 * Buffers zero_src, B_I/B_Q/B_U, top_bc are caller-provided scratch
 * (hot loop; avoids per-order allocation). */
static int add_flat_intrefl_downfield_order_complete(const rt_atm_t *atm, int m,
                                                     const double *Rww_M,
                                                     double *I, double *Q, double *U,
                                                     double *zero_src,
                                                     double *B_I, double *B_Q, double *B_U,
                                                     double *top_bc) {
    if (!atm || !Rww_M || !I || !Q || !U || !zero_src || !B_I || !B_Q || !B_U || !top_bc)
        return -1;
    const int nt = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;
    const double msign = 1.0;  /* specular water reflection preserves azimuth; no (-1)^m */
    double *tb_I = top_bc;
    double *tb_Q = top_bc + n_mu;
    double *tb_U = top_bc + 2 * n_mu;
    memset(top_bc, 0, (size_t)(3 * n_mu) * sizeof(double));
    for (int jp = 1; jp <= n_mu; ++jp) {
        const size_t up = (size_t)0 * (size_t)dirs + (size_t)(+jp + n_mu);
        const double I_up = I[up], Q_up = Q[up], U_up = U[up];
        const double *M = Rww_M + (size_t)(jp - 1) * 9;
        tb_I[jp - 1] = msign * (M[0] * I_up + M[1] * Q_up + M[2] * U_up);
        tb_Q[jp - 1] = msign * (M[3] * I_up + M[4] * Q_up + M[5] * U_up);
        tb_U[jp - 1] = msign * (M[6] * I_up + M[7] * Q_up + M[8] * U_up);
    }
    memset(B_I, 0, buf_size * sizeof(double));
    memset(B_Q, 0, buf_size * sizeof(double));
    memset(B_U, 0, buf_size * sizeof(double));
    if (rt_solver_integrate_bcs(atm, m, zero_src, NULL, tb_I,
                                RT_INTEGRATION_METHOD_LINEAR, B_I) != 0 ||
        rt_solver_integrate_bcs(atm, m, zero_src, NULL, tb_Q,
                                RT_INTEGRATION_METHOD_LINEAR, B_Q) != 0 ||
        rt_solver_integrate_bcs(atm, m, zero_src, NULL, tb_U,
                                RT_INTEGRATION_METHOD_LINEAR, B_U) != 0) {
        return -2;
    }
    for (size_t idx = 0; idx < buf_size; ++idx) {
        I[idx] += B_I[idx];
        Q[idx] += B_Q[idx];
        U[idx] += B_U[idx];
    }
    return 0;
}

/* ========================================================================
 * rt_solver_sos_pol_intrefl — P0-A prototype: SOS with water-side surface
 * internal-reflection feedback as a TOP downward boundary condition.
 * Copy of rt_solver_sos_pol; the only difference is that at each order the
 * downward field at k=0 is seeded with M_Rww(mu_j) * [I,Q,U]_up^m(k=0,+j)
 * from the previous order (full Fresnel Mueller; flat surface; TIR via
 * M_Rww=identity). Summed over orders the top downward field -> M_Rww applied
 * to the total upwelling, i.e. the self-consistent reflecting boundary.
 * ====================================================================== */
int rt_solver_sos_pol_intrefl(const rt_atm_t *atm, int m,
                              const rt_legendre_workspace_t *ws,
                              const double *primary_I,
                              const double *primary_Q,
                              const double *primary_U,
                              const double *Rww_M,
                              const rt_solver_sos_options_t *opts,
                              double *total_I,
                              double *total_Q,
                              double *total_U,
                              rt_solver_sos_result_t *result) {
    if (!atm || !ws || !primary_I || !primary_Q || !primary_U || !Rww_M ||
        !total_I || !total_Q || !total_U || !result) return -1;
    if (m < 0)                                                     return -1;
    if (atm->n_mu != ws->n_mu)                                     return -1;

    rt_solver_sos_options_t opt_default = rt_solver_sos_options_default();
    if (!opts) opts = &opt_default;
    if (opts->max_iterations < 1)                                  return -1;
    if (!(opts->tolerance > 0.0))                                  return -1;

    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;
    const double msign = 1.0;  /* specular water reflection preserves azimuth; no (-1)^m */

    memcpy(total_I, primary_I, buf_size * sizeof(double));
    memcpy(total_Q, primary_Q, buf_size * sizeof(double));
    memcpy(total_U, primary_U, buf_size * sizeof(double));

    dump_water_order_trace_if_requested("intrefl", m, 1, atm, primary_I, primary_Q, primary_U);

    if (opts->max_iterations == 1) {
        result->n_orders_used  = 1;
        result->converged      = 1;
        result->final_residual = 0.0;
        return 0;
    }

    double *I_prev = calloc(buf_size, sizeof(double));
    double *Q_prev = calloc(buf_size, sizeof(double));
    double *U_prev = calloc(buf_size, sizeof(double));
    double *J_I    = calloc(buf_size, sizeof(double));
    double *J_Q    = calloc(buf_size, sizeof(double));
    double *J_U    = calloc(buf_size, sizeof(double));
    double *I_curr = calloc(buf_size, sizeof(double));
    double *Q_curr = calloc(buf_size, sizeof(double));
    double *U_curr = calloc(buf_size, sizeof(double));
    double *top_bc = calloc((size_t)(3 * n_mu), sizeof(double)); /* [I|Q|U] reflected-down BC */
    double *zero_src = calloc(buf_size, sizeof(double));
    double *B_I = calloc(buf_size, sizeof(double));
    double *B_Q = calloc(buf_size, sizeof(double));
    double *B_U = calloc(buf_size, sizeof(double));
    if (!I_prev || !Q_prev || !U_prev || !J_I || !J_Q || !J_U ||
        !I_curr || !Q_curr || !U_curr || !top_bc || !zero_src || !B_I || !B_Q || !B_U) {
        free(I_prev); free(Q_prev); free(U_prev);
        free(J_I); free(J_Q); free(J_U);
        free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
        free(zero_src); free(B_I); free(B_Q); free(B_U);
        return -1;
    }
    memcpy(I_prev, primary_I, buf_size * sizeof(double));
    memcpy(Q_prev, primary_Q, buf_size * sizeof(double));
    memcpy(U_prev, primary_U, buf_size * sizeof(double));

    const int boundary_complete_intrefl = getenv("OCRT_INTREFL_BOUNDARY_COMPLETE") ? 1 : 0;
    if (boundary_complete_intrefl) {
        if (add_flat_intrefl_downfield_order_complete(atm, m, Rww_M,
                                                      I_prev, Q_prev, U_prev,
                                                      zero_src, B_I, B_Q, B_U, top_bc) != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
            free(zero_src); free(B_I); free(B_Q); free(B_U);
            return -2;
        }
    }

    int n = 1;
    double residual = 0.0;
    int converged = 0;

    for (int it = 2; it <= opts->max_iterations; it++) {
        /* Diagnostic/proposed ordering fix: the previous-order radiance field
         * used to build J_n should already include water-side specular internal
         * reflection at the top boundary.  Without this, scattering of the
         * reflected part is delayed by one SOS order.  Keep this env-gated
         * until OSOAA parity and energy accounting are fully checked. */
        const int preseed_transport_only = getenv("OCRT_INTREFL_SOURCE_PRESEED_TRANSPORT_ONLY") ? 1 : 0;
        const int preseed_intrefl_transport = (getenv("OCRT_INTREFL_SOURCE_PRESEED_TRANSPORT") || preseed_transport_only) ? 1 : 0;
        const int preseed_intrefl_source = (getenv("OCRT_INTREFL_SOURCE_PRESEED") || getenv("OCRT_INTREFL_SOURCE_PRESEED_ONLY")) ? 1 : 0;
        const int preseed_source_only = getenv("OCRT_INTREFL_SOURCE_PRESEED_ONLY") ? 1 : 0;
        double *tb_I = top_bc, *tb_Q = top_bc + n_mu, *tb_U = top_bc + 2 * n_mu;

        /* Step15 diagnostic: boundary-complete the previous-order field used
         * for source construction. Step14 inserted only the reflected top slot
         * I(0,-mu).  Here Rww*I_{n-1}(0,+mu) is propagated down through the
         * column with zero volume source, and the resulting downwelling branch
         * is added to a temporary source-build copy.  I_prev is left unchanged
         * unless the older top-slot preseed mode is explicitly requested. */
        const double *src_I_for_J = I_prev;
        const double *src_Q_for_J = Q_prev;
        const double *src_U_for_J = U_prev;
        double *I_src = NULL, *Q_src = NULL, *U_src = NULL;
        if (preseed_intrefl_transport) {
            I_src = (double *)malloc(buf_size * sizeof(double));
            Q_src = (double *)malloc(buf_size * sizeof(double));
            U_src = (double *)malloc(buf_size * sizeof(double));
            double *Z_src = (double *)calloc(buf_size, sizeof(double));
            double *R_I = (double *)malloc(buf_size * sizeof(double));
            double *R_Q = (double *)malloc(buf_size * sizeof(double));
            double *R_U = (double *)malloc(buf_size * sizeof(double));
            if (!I_src || !Q_src || !U_src || !Z_src || !R_I || !R_Q || !R_U) {
                free(I_src); free(Q_src); free(U_src); free(Z_src);
                free(R_I); free(R_Q); free(R_U);
                free(I_prev); free(Q_prev); free(U_prev);
                free(J_I); free(J_Q); free(J_U);
                free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
                return -1;
            }
            memcpy(I_src, I_prev, buf_size * sizeof(double));
            memcpy(Q_src, Q_prev, buf_size * sizeof(double));
            memcpy(U_src, U_prev, buf_size * sizeof(double));
            for (int jp = 1; jp <= n_mu; ++jp) {
                size_t up = (size_t)0 * (size_t)dirs + (size_t)(+jp + n_mu);
                double I_up = I_prev[up], Q_up = Q_prev[up], U_up = U_prev[up];
                const double *M = Rww_M + (size_t)(jp - 1) * 9;
                tb_I[jp - 1] = msign * (M[0] * I_up + M[1] * Q_up + M[2] * U_up);
                tb_Q[jp - 1] = msign * (M[3] * I_up + M[4] * Q_up + M[5] * U_up);
                tb_U[jp - 1] = msign * (M[6] * I_up + M[7] * Q_up + M[8] * U_up);
            }
            if (rt_solver_integrate_bcs(atm, m, Z_src, NULL, tb_I,
                                        RT_INTEGRATION_METHOD_LINEAR, R_I) != 0 ||
                rt_solver_integrate_bcs(atm, m, Z_src, NULL, tb_Q,
                                        RT_INTEGRATION_METHOD_LINEAR, R_Q) != 0 ||
                rt_solver_integrate_bcs(atm, m, Z_src, NULL, tb_U,
                                        RT_INTEGRATION_METHOD_LINEAR, R_U) != 0) {
                free(I_src); free(Q_src); free(U_src); free(Z_src);
                free(R_I); free(R_Q); free(R_U);
                free(I_prev); free(Q_prev); free(U_prev);
                free(J_I); free(J_Q); free(J_U);
                free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
                return -2;
            }
            for (int k = 0; k <= nt; ++k) {
                for (int jp = 1; jp <= n_mu; ++jp) {
                    size_t dn = (size_t)k * (size_t)dirs + (size_t)(-jp + n_mu);
                    I_src[dn] += R_I[dn];
                    Q_src[dn] += R_Q[dn];
                    U_src[dn] += R_U[dn];
                }
            }
            free(Z_src); free(R_I); free(R_Q); free(R_U);
            src_I_for_J = I_src;
            src_Q_for_J = Q_src;
            src_U_for_J = U_src;
        } else if (preseed_intrefl_source) {
            for (int jp = 1; jp <= n_mu; ++jp) {
                size_t up = (size_t)0 * (size_t)dirs + (size_t)(+jp + n_mu);
                size_t dn = (size_t)0 * (size_t)dirs + (size_t)(-jp + n_mu);
                double I_up = I_prev[up], Q_up = Q_prev[up], U_up = U_prev[up];
                const double *M = Rww_M + (size_t)(jp - 1) * 9;
                tb_I[jp - 1] = msign * (M[0] * I_up + M[1] * Q_up + M[2] * U_up);
                tb_Q[jp - 1] = msign * (M[3] * I_up + M[4] * Q_up + M[5] * U_up);
                tb_U[jp - 1] = msign * (M[6] * I_up + M[7] * Q_up + M[8] * U_up);
                I_prev[dn] = tb_I[jp - 1];
                Q_prev[dn] = tb_Q[jp - 1];
                U_prev[dn] = tb_U[jp - 1];
            }
        }

        int rc = rt_sos_operator_apply_vector(atm, m, ws, src_I_for_J, src_Q_for_J, src_U_for_J,
                                       J_I, J_Q, J_U);
        dump_water_source_trace_if_requested("intrefl", "ms_source", m, it,
                                             atm, J_I, J_Q, J_U);
        dump_water_source_pref_trace_if_requested("intrefl", "ms_source", m, it,
                                                  atm, NULL, NULL, NULL,
                                                  J_I, J_Q, J_U);
        if (rc != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
            return rc;
        }

        /* Step19 diagnostic-only: decompose the propagated order-n numerator
         * by the sign of the previous-order incoming field used in the source
         * integral.  At the water top boundary, the -mu previous field is the
         * internal-reflection feedback lineage.  We construct two source fields:
         * one with previous +mu streams only and one with previous -mu streams
         * only, then integrate each with no top BC. */
        if (getenv("OCRT_DUMP_WATER_ORDER_ORIGIN_TRACE")) {
            double *Ip_up = (double *)calloc(buf_size, sizeof(double));
            double *Qp_up = (double *)calloc(buf_size, sizeof(double));
            double *Up_up = (double *)calloc(buf_size, sizeof(double));
            double *Ip_dn = (double *)calloc(buf_size, sizeof(double));
            double *Qp_dn = (double *)calloc(buf_size, sizeof(double));
            double *Up_dn = (double *)calloc(buf_size, sizeof(double));
            double *Ju_I = (double *)calloc(buf_size, sizeof(double));
            double *Ju_Q = (double *)calloc(buf_size, sizeof(double));
            double *Ju_U = (double *)calloc(buf_size, sizeof(double));
            double *Jd_I = (double *)calloc(buf_size, sizeof(double));
            double *Jd_Q = (double *)calloc(buf_size, sizeof(double));
            double *Jd_U = (double *)calloc(buf_size, sizeof(double));
            double *Fu_I = (double *)calloc(buf_size, sizeof(double));
            double *Fu_Q = (double *)calloc(buf_size, sizeof(double));
            double *Fu_U = (double *)calloc(buf_size, sizeof(double));
            double *Fd_I = (double *)calloc(buf_size, sizeof(double));
            double *Fd_Q = (double *)calloc(buf_size, sizeof(double));
            double *Fd_U = (double *)calloc(buf_size, sizeof(double));
            if (Ip_up && Qp_up && Up_up && Ip_dn && Qp_dn && Up_dn &&
                Ju_I && Ju_Q && Ju_U && Jd_I && Jd_Q && Jd_U &&
                Fu_I && Fu_Q && Fu_U && Fd_I && Fd_Q && Fd_U) {
                for (int kk = 0; kk <= nt; ++kk) {
                    for (int jj = 1; jj <= n_mu; ++jj) {
                        size_t pidx = (size_t)kk * (size_t)dirs + (size_t)(+jj + n_mu);
                        size_t didx = (size_t)kk * (size_t)dirs + (size_t)(-jj + n_mu);
                        Ip_up[pidx] = src_I_for_J[pidx]; Qp_up[pidx] = src_Q_for_J[pidx]; Up_up[pidx] = src_U_for_J[pidx];
                        Ip_dn[didx] = src_I_for_J[didx]; Qp_dn[didx] = src_Q_for_J[didx]; Up_dn[didx] = src_U_for_J[didx];
                    }
                }
                if (rt_sos_operator_apply_vector(atm, m, ws, Ip_up, Qp_up, Up_up, Ju_I, Ju_Q, Ju_U) == 0 &&
                    rt_sos_operator_apply_vector(atm, m, ws, Ip_dn, Qp_dn, Up_dn, Jd_I, Jd_Q, Jd_U) == 0 &&
                    rt_solver_integrate_bcs(atm, m, Ju_I, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, Fu_I) == 0 &&
                    rt_solver_integrate_bcs(atm, m, Ju_Q, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, Fu_Q) == 0 &&
                    rt_solver_integrate_bcs(atm, m, Ju_U, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, Fu_U) == 0 &&
                    rt_solver_integrate_bcs(atm, m, Jd_I, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, Fd_I) == 0 &&
                    rt_solver_integrate_bcs(atm, m, Jd_Q, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, Fd_Q) == 0 &&
                    rt_solver_integrate_bcs(atm, m, Jd_U, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, Fd_U) == 0) {
                    dump_water_order_component_trace_if_requested("intrefl", "src_from_prev_up", m, it, atm, Fu_I, Fu_Q, Fu_U);
                    dump_water_order_component_trace_if_requested("intrefl", "src_from_prev_down", m, it, atm, Fd_I, Fd_Q, Fd_U);
                }
                    }
            free(Ip_up); free(Qp_up); free(Up_up); free(Ip_dn); free(Qp_dn); free(Up_dn);
            free(Ju_I); free(Ju_Q); free(Ju_U); free(Jd_I); free(Jd_Q); free(Jd_U);
            free(Fu_I); free(Fu_Q); free(Fu_U); free(Fd_I); free(Fd_Q); free(Fd_U);
        }

        /* Step87 diagnostic-only isolated previous-Stokes propagation dump for
         * internal-reflection water SOS path.  Use the actual source-build field
         * src_*_for_J, because it may include boundary-completed/preseeded
         * reflected downfield.  Integrate isolated sources with no delayed BC
         * to measure volume source -> propagated 0- transfer. */
        {
            const char *cp = getenv("OCRT_DUMP_WATER_MS_CHANNEL_PROP_TRACE");
            if (cp && cp[0]) {
                int ok_m = 1;
                const char *sm = getenv("OCRT_MS_CHANNEL_PROP_M");
                if (sm && sm[0] && m != atoi(sm)) ok_m = 0;
                int ok_o = 1;
                const char *so = getenv("OCRT_MS_CHANNEL_PROP_ORDER");
                if (so && so[0] && it != atoi(so)) ok_o = 0;
                if (ok_m && ok_o) {
                    double *Z  = calloc(buf_size, sizeof(double));
                    double *cJI = calloc(buf_size, sizeof(double));
                    double *cJQ = calloc(buf_size, sizeof(double));
                    double *cJU = calloc(buf_size, sizeof(double));
                    double *cFI = calloc(buf_size, sizeof(double));
                    double *cFQ = calloc(buf_size, sizeof(double));
                    double *cFU = calloc(buf_size, sizeof(double));
                    if (Z && cJI && cJQ && cJU && cFI && cFQ && cFU) {
                        FILE *pf = fopen(cp, "a+");
                        if (pf) {
                            fseek(pf, 0, SEEK_END);
                            if (ftell(pf) == 0L) {
                                fprintf(pf,
                                        "solver,wavelength_nm,route,m,order,prev_channel,out_component,source_abs,source_sum,field_all_abs,field_all_sum,top_up_abs,top_up_sum,bottom_down_abs,bottom_down_sum,source_to_top_abs_eff,source_to_all_abs_eff\n");
                            }
                            const char *pnames[3] = {"prev_I", "prev_Q", "prev_U"};
                            for (int pc = 0; pc < 3; ++pc) {
                                const double *inI = (pc == 0) ? src_I_for_J : Z;
                                const double *inQ = (pc == 1) ? src_Q_for_J : Z;
                                const double *inU = (pc == 2) ? src_U_for_J : Z;
                                for (size_t bi = 0; bi < buf_size; ++bi) {
                                    cJI[bi] = cJQ[bi] = cJU[bi] = 0.0;
                                    cFI[bi] = cFQ[bi] = cFU[bi] = 0.0;
                                }
                                int rc_iso = rt_sos_operator_apply_vector(atm, m, ws, inI, inQ, inU, cJI, cJQ, cJU);
                                if (rc_iso == 0) {
                                    int ri = rt_solver_integrate_bcs(atm, m, cJI, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, cFI);
                                    int rq = rt_solver_integrate_bcs(atm, m, cJQ, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, cFQ);
                                    int ru = rt_solver_integrate_bcs(atm, m, cJU, NULL, NULL, RT_INTEGRATION_METHOD_LINEAR, cFU);
                                    if (ri == 0 && rq == 0 && ru == 0) {
                                        const double *srcs[3] = {cJI, cJQ, cJU};
                                        const double *flds[3] = {cFI, cFQ, cFU};
                                        const char *cnames[3] = {"I", "Q", "U"};
                                        for (int cc = 0; cc < 3; ++cc) {
                                            double s_abs = 0.0, s_sum = 0.0;
                                            double f_abs = 0.0, f_sum = 0.0;
                                            double tu_abs = 0.0, tu_sum = 0.0;
                                            double bd_abs = 0.0, bd_sum = 0.0;
                                            for (int kk2 = 0; kk2 <= nt; ++kk2) {
                                                for (int jj = -n_mu; jj <= n_mu; ++jj) {
                                                    if (jj == 0) continue;
                                                    size_t idx = (size_t)kk2 * (size_t)dirs + (size_t)(jj + n_mu);
                                                    s_abs += fabs(srcs[cc][idx]); s_sum += srcs[cc][idx];
                                                    f_abs += fabs(flds[cc][idx]); f_sum += flds[cc][idx];
                                                }
                                            }
                                            for (int jj = 1; jj <= n_mu; ++jj) {
                                                size_t top_idx = (size_t)0 * (size_t)dirs + (size_t)(jj + n_mu);
                                                tu_abs += fabs(flds[cc][top_idx]); tu_sum += flds[cc][top_idx];
                                                size_t bot_idx = (size_t)nt * (size_t)dirs + (size_t)(-jj + n_mu);
                                                bd_abs += fabs(flds[cc][bot_idx]); bd_sum += flds[cc][bot_idx];
                                            }
                                            const double eff_top = (s_abs > 0.0) ? tu_abs / s_abs : 0.0;
                                            const double eff_all = (s_abs > 0.0) ? f_abs / s_abs : 0.0;
                                            double wl = 0.0;
                                            const char *wls = getenv("OCRT_TRACE_WAVELENGTH_NM");
                                            if (wls && wls[0]) wl = atof(wls);
                                            fprintf(pf,
                                                    "OCRT,%.17g,intrefl,%d,%d,%s,%s,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                                                    wl, m, it, pnames[pc], cnames[cc],
                                                    s_abs, s_sum, f_abs, f_sum, tu_abs, tu_sum, bd_abs, bd_sum, eff_top, eff_all);
                                        }
                                    }
                                }
                            }
                            fclose(pf);
                        }
                    }
                    free(Z); free(cJI); free(cJQ); free(cJU); free(cFI); free(cFQ); free(cFU);
                }
            }
        }
        free(I_src); free(Q_src); free(U_src);

        if (!boundary_complete_intrefl) {
            /* P0-A legacy: reflect previous-order upwelling at surface (k=0)
             * into the downward top BC of the next integrated field. */
            tb_I = top_bc; tb_Q = top_bc + n_mu; tb_U = top_bc + 2 * n_mu;
            for (int jp = 1; jp <= n_mu; ++jp) {
                size_t up = (size_t)0 * (size_t)dirs + (size_t)(+jp + n_mu);
                double I_up = I_prev[up], Q_up = Q_prev[up], U_up = U_prev[up];
                const double *M = Rww_M + (size_t)(jp - 1) * 9;
                tb_I[jp - 1] = msign * (M[0] * I_up + M[1] * Q_up + M[2] * U_up);
                tb_Q[jp - 1] = msign * (M[3] * I_up + M[4] * Q_up + M[5] * U_up);
                tb_U[jp - 1] = msign * (M[6] * I_up + M[7] * Q_up + M[8] * U_up);
            }

            dump_water_feedback_chain_if_requested("intrefl", m, it, atm, ws, zero_src, tb_I, tb_Q, tb_U);

            const double *bc_I_for_integrate = (preseed_source_only || preseed_transport_only) ? NULL : tb_I;
            const double *bc_Q_for_integrate = (preseed_source_only || preseed_transport_only) ? NULL : tb_Q;
            const double *bc_U_for_integrate = (preseed_source_only || preseed_transport_only) ? NULL : tb_U;
            if (rt_solver_integrate_bcs(atm, m, J_I, NULL, bc_I_for_integrate,
                                        RT_INTEGRATION_METHOD_LINEAR, I_curr) != 0 ||
                rt_solver_integrate_bcs(atm, m, J_Q, NULL, bc_Q_for_integrate,
                                        RT_INTEGRATION_METHOD_LINEAR, Q_curr) != 0 ||
                rt_solver_integrate_bcs(atm, m, J_U, NULL, bc_U_for_integrate,
                                        RT_INTEGRATION_METHOD_LINEAR, U_curr) != 0) {
                free(I_prev); free(Q_prev); free(U_prev);
                free(J_I); free(J_Q); free(J_U);
                free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
                free(zero_src); free(B_I); free(B_Q); free(B_U);
                return -2;
            }
        } else {
            /* Order-complete diagnostic: J_n is built from boundary-completed
             * I_{n-1}; integrate J_n with no delayed top BC, then complete the
             * resulting I_n by adding its same-order reflected downward field. */
            if (rt_solver_integrate_bcs(atm, m, J_I, NULL, NULL,
                                        RT_INTEGRATION_METHOD_LINEAR, I_curr) != 0 ||
                rt_solver_integrate_bcs(atm, m, J_Q, NULL, NULL,
                                        RT_INTEGRATION_METHOD_LINEAR, Q_curr) != 0 ||
                rt_solver_integrate_bcs(atm, m, J_U, NULL, NULL,
                                        RT_INTEGRATION_METHOD_LINEAR, U_curr) != 0) {
                free(I_prev); free(Q_prev); free(U_prev);
                free(J_I); free(J_Q); free(J_U);
                free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
                free(zero_src); free(B_I); free(B_Q); free(B_U);
                return -2;
            }
            if (add_flat_intrefl_downfield_order_complete(atm, m, Rww_M,
                                                          I_curr, Q_curr, U_curr,
                                                          zero_src, B_I, B_Q, B_U, top_bc) != 0) {
                free(I_prev); free(Q_prev); free(U_prev);
                free(J_I); free(J_Q); free(J_U);
                free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
                free(zero_src); free(B_I); free(B_Q); free(B_U);
                return -2;
            }
        }

        dump_water_order_trace_if_requested("intrefl", m, it, atm, I_curr, Q_curr, U_curr);
        dump_water_order_source_feedback_trace_if_requested("intrefl", m, it, atm,
                                                            J_I, J_Q, J_U,
                                                            I_curr, Q_curr, U_curr);

        for (int k = 0; k <= nt; k++) {
            for (int j = -n_mu; j <= n_mu; j++) {
                if (j == 0) continue;
                size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
                total_I[idx] += I_curr[idx];
                total_Q[idx] += Q_curr[idx];
                total_U[idx] += U_curr[idx];
            }
        }
        n = it;

        double max_in_I = field_max_abs(I_curr, nt, n_mu);
        double max_in_Q = field_max_abs(Q_curr, nt, n_mu);
        double max_in_U = field_max_abs(U_curr, nt, n_mu);
        double max_in   = max_in_I;
        if (max_in_Q > max_in) max_in = max_in_Q;
        if (max_in_U > max_in) max_in = max_in_U;

        double max_t_I = field_max_abs(total_I, nt, n_mu);
        double max_t_Q = field_max_abs(total_Q, nt, n_mu);
        double max_t_U = field_max_abs(total_U, nt, n_mu);
        double max_t   = max_t_I;
        if (max_t_Q > max_t) max_t = max_t_Q;
        if (max_t_U > max_t) max_t = max_t_U;

        residual = (max_t > 0.0) ? max_in / max_t : 0.0;
        if (residual < opts->tolerance) {
            converged = 1;
            break;
        }

        memcpy(I_prev, I_curr, buf_size * sizeof(double));
        memcpy(Q_prev, Q_curr, buf_size * sizeof(double));
        memcpy(U_prev, U_curr, buf_size * sizeof(double));
    }

    free(I_prev); free(Q_prev); free(U_prev);
    free(J_I); free(J_Q); free(J_U);
    free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
    free(zero_src); free(B_I); free(B_Q); free(B_U);

    result->n_orders_used  = n;
    result->converged      = converged;
    result->final_residual = residual;
    return 0;
}

/* ========================================================================
 * rt_solver_sos_pol_intrefl_rough — B2: rough-surface (Cox-Munk) water-side
 * internal-reflection feedback. Identical to rt_solver_sos_pol_intrefl but the
 * top downward BC seed reflects the previous-order upwelling through the
 * angle-coupling m-mode kernel Rww_K (surface_R_ww_coxmunk_fourier_kernel)
 * instead of the per-node specular Mueller M_Rww. Used at wind > 0.
 * ======================================================================== */
int rt_solver_sos_pol_intrefl_rough(const rt_atm_t *atm, int m,
                              const rt_legendre_workspace_t *ws,
                              const double *primary_I,
                              const double *primary_Q,
                              const double *primary_U,
                              const double *Rww_K,   /* m-mode Cox-Munk water-side kernel [n_mu*n_mu*9] */
                              const rt_solver_sos_options_t *opts,
                              double *total_I,
                              double *total_Q,
                              double *total_U,
                              rt_solver_sos_result_t *result) {
    if (!atm || !ws || !primary_I || !primary_Q || !primary_U || !Rww_K ||
        !total_I || !total_Q || !total_U || !result) return -1;
    if (m < 0)                                                     return -1;
    if (atm->n_mu != ws->n_mu)                                     return -1;

    rt_solver_sos_options_t opt_default = rt_solver_sos_options_default();
    if (!opts) opts = &opt_default;
    if (opts->max_iterations < 1)                                  return -1;
    if (!(opts->tolerance > 0.0))                                  return -1;

    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;
    const double msign = 1.0;  /* specular water reflection preserves azimuth; no (-1)^m */

    memcpy(total_I, primary_I, buf_size * sizeof(double));
    memcpy(total_Q, primary_Q, buf_size * sizeof(double));
    memcpy(total_U, primary_U, buf_size * sizeof(double));

    dump_water_order_trace_if_requested("intrefl_rough", m, 1, atm, primary_I, primary_Q, primary_U);

    if (opts->max_iterations == 1) {
        result->n_orders_used  = 1;
        result->converged      = 1;
        result->final_residual = 0.0;
        return 0;
    }

    double *I_prev = calloc(buf_size, sizeof(double));
    double *Q_prev = calloc(buf_size, sizeof(double));
    double *U_prev = calloc(buf_size, sizeof(double));
    double *J_I    = calloc(buf_size, sizeof(double));
    double *J_Q    = calloc(buf_size, sizeof(double));
    double *J_U    = calloc(buf_size, sizeof(double));
    double *I_curr = calloc(buf_size, sizeof(double));
    double *Q_curr = calloc(buf_size, sizeof(double));
    double *U_curr = calloc(buf_size, sizeof(double));
    double *top_bc = calloc((size_t)(3 * n_mu), sizeof(double)); /* [I|Q|U] reflected-down BC */
    if (!I_prev || !Q_prev || !U_prev || !J_I || !J_Q || !J_U ||
        !I_curr || !Q_curr || !U_curr || !top_bc) {
        free(I_prev); free(Q_prev); free(U_prev);
        free(J_I); free(J_Q); free(J_U);
        free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
        return -1;
    }
    memcpy(I_prev, primary_I, buf_size * sizeof(double));
    memcpy(Q_prev, primary_Q, buf_size * sizeof(double));
    memcpy(U_prev, primary_U, buf_size * sizeof(double));

    int n = 1;
    double residual = 0.0;
    int converged = 0;

    for (int it = 2; it <= opts->max_iterations; it++) {
        int rc = rt_sos_operator_apply_vector(atm, m, ws, I_prev, Q_prev, U_prev,
                                       J_I, J_Q, J_U);
        dump_water_source_trace_if_requested("intrefl_rough", "ms_source", m, it,
                                             atm, J_I, J_Q, J_U);
        dump_water_source_pref_trace_if_requested("intrefl_rough", "ms_source", m, it,
                                                  atm, NULL, NULL, NULL,
                                                  J_I, J_Q, J_U);
        if (rc != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
            return rc;
        }

        /* B2 rough: reflect previous-order upwelling at the surface (k=0) into
         * the downward top BC via the Cox-Munk WATER-SIDE internal-reflection
         * m-mode kernel K^m (angle-coupling, not per-node specular). Contraction
         * follows the same Fourier reconstruction convention as the rough
         * atmosphere-surface BC: C_0=2*pi, C_{m>0}=pi, because radiance fields are
         * reconstructed as X(phi)=X^0+2*sum_{m>0}X^m cos/sin(m phi). As wind->0
         * the slope PDF collapses to the specular delta and this reduces to the
         * flat per-node M_Rww limit.
         *   [I,Q,U]_dn^m(mu_j) = sum_k (C_m mu_k w_k) K^m[j,k] [I,Q,U]_up^m(mu_k) */
        double *tb_I = top_bc, *tb_Q = top_bc + n_mu, *tb_U = top_bc + 2 * n_mu;
        for (int jp = 1; jp <= n_mu; ++jp) {
            double aI = 0.0, aQ = 0.0, aU = 0.0;
            for (int kp = 1; kp <= n_mu; ++kp) {
                size_t up = (size_t)0 * (size_t)dirs + (size_t)(+kp + n_mu);
                double I_up = I_prev[up], Q_up = Q_prev[up], U_up = U_prev[up];
                const double *K = Rww_K + ((size_t)(jp - 1) * (size_t)n_mu + (size_t)(kp - 1)) * 9;
                const double az_factor = (m == 0) ? (2.0 * RT_F_SOLAR_PI) : RT_F_SOLAR_PI;
                double wk = az_factor * atm->rm[+kp] * atm->gb[+kp];
                aI += wk * (K[0] * I_up + K[1] * Q_up + K[2] * U_up);
                aQ += wk * (K[3] * I_up + K[4] * Q_up + K[5] * U_up);
                aU += wk * (K[6] * I_up + K[7] * Q_up + K[8] * U_up);
            }
            tb_I[jp - 1] = msign * aI;
            tb_Q[jp - 1] = msign * aQ;
            tb_U[jp - 1] = msign * aU;
        }

        if (rt_solver_integrate_bcs(atm, m, J_I, NULL, tb_I,
                                    RT_INTEGRATION_METHOD_LINEAR, I_curr) != 0 ||
            rt_solver_integrate_bcs(atm, m, J_Q, NULL, tb_Q,
                                    RT_INTEGRATION_METHOD_LINEAR, Q_curr) != 0 ||
            rt_solver_integrate_bcs(atm, m, J_U, NULL, tb_U,
                                    RT_INTEGRATION_METHOD_LINEAR, U_curr) != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr); free(top_bc);
            return -2;
        }

        dump_water_order_trace_if_requested("intrefl_rough", m, it, atm, I_curr, Q_curr, U_curr);

        for (int k = 0; k <= nt; k++) {
            for (int j = -n_mu; j <= n_mu; j++) {
                if (j == 0) continue;
                size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
                total_I[idx] += I_curr[idx];
                total_Q[idx] += Q_curr[idx];
                total_U[idx] += U_curr[idx];
            }
        }
        n = it;

        double max_in_I = field_max_abs(I_curr, nt, n_mu);
        double max_in_Q = field_max_abs(Q_curr, nt, n_mu);
        double max_in_U = field_max_abs(U_curr, nt, n_mu);
        double max_in   = max_in_I;
        if (max_in_Q > max_in) max_in = max_in_Q;
        if (max_in_U > max_in) max_in = max_in_U;

        double max_t_I = field_max_abs(total_I, nt, n_mu);
        double max_t_Q = field_max_abs(total_Q, nt, n_mu);
        double max_t_U = field_max_abs(total_U, nt, n_mu);
        double max_t   = max_t_I;
        if (max_t_Q > max_t) max_t = max_t_Q;
        if (max_t_U > max_t) max_t = max_t_U;

        residual = (max_t > 0.0) ? max_in / max_t : 0.0;
        if (residual < opts->tolerance) {
            converged = 1;
            break;
        }

        memcpy(I_prev, I_curr, buf_size * sizeof(double));
        memcpy(Q_prev, Q_curr, buf_size * sizeof(double));
        memcpy(U_prev, U_curr, buf_size * sizeof(double));
    }

    free(I_prev); free(Q_prev); free(U_prev);
    free(J_I); free(J_Q); free(J_U);
    free(I_curr); free(Q_curr); free(U_curr); free(top_bc);

    result->n_orders_used  = n;
    result->converged      = converged;
    result->final_residual = residual;
    return 0;
}


/* ========================================================================
 * rt_solver_sos_pol_with_surface — SOS with vector flat/Cox-Munk Fresnel BC
 *
 * Phase A (Path 3 graft) implementation: flat Fresnel surface coupling.
 *
 * For flat ocean (specular only):
 *   I_up^m(τ_max, μ_o) = (-1)^m · M^flat(μ_o) · [I_dn^m(τ_max, μ_o)]
 *
 * The factor (-1)^m comes from the φ → φ+π specular relation in Fourier
 * expansion (cos(m(φ+π)) = (-1)^m cos(mφ), same for sin).
 *
 * Cox-Munk (wind_speed > 0): Phase B will replace the specular relation
 * with numerical φ-quadrature evaluation of surface_R_coxmunk_trig over
 * all (μ_o, μ_in) pairs.
 *
 * IMPORTANT: This function adds atm-surface MS coupling to total_X but
 * does NOT add the *direct solar specular bounce* (sunglint). That term
 * is added separately in the driver (currently always decoupled when
 * decouple_sunglint=1 in the rt_case).
 * ====================================================================== */

int rt_solver_sos_pol_with_surface(const rt_atm_t *atm, int m,
                                    const rt_legendre_workspace_t *ws,
                                    const double *primary_I,
                                    const double *primary_Q,
                                    const double *primary_U,
                                    const rt_solver_sos_options_t *opts,
                                    rt_surface_kind_t surface_kind,
                                    double n_water,
                                    double wind_speed,
                                    int    sigma_type,
                                    int    q_convention,
                                    double *total_I,
                                    double *total_Q,
                                    double *total_U,
                                    rt_solver_sos_result_t *result) {
    if (!atm || !ws || !primary_I || !primary_Q || !primary_U ||
        !total_I || !total_Q || !total_U || !result) return -1;
    if (m < 0)                                                     return -1;
    if (atm->n_mu != ws->n_mu)                                     return -1;
    if (n_water <= 0.0)                                            return -1;

    rt_solver_sos_options_t opt_default = rt_solver_sos_options_default();
    if (!opts) opts = &opt_default;
    if (opts->max_iterations < 1)                                  return -1;
    if (!(opts->tolerance > 0.0))                                  return -1;

    const int nt   = atm->n_layers;
    const int n_mu = atm->n_mu;
    const int dirs = 2 * n_mu + 1;
    const size_t buf_size = (size_t)(nt + 1) * (size_t)dirs;
    rt_surface_boundary_t surface_boundary;
    int surface_status = rt_surface_boundary_init(
        &surface_boundary, atm, m, surface_kind, n_water, wind_speed,
        sigma_type, q_convention);
    if (surface_status != 0) return surface_status;

    const int use_surf_first_order =
        (surface_kind == RT_SURFACE_BLACK_FRESNEL_OCEAN && wind_speed > 0.0);

    /* total_X starts as a copy of primary_X. */
    memcpy(total_I, primary_I, buf_size * sizeof(double));
    memcpy(total_Q, primary_Q, buf_size * sizeof(double));
    memcpy(total_U, primary_U, buf_size * sizeof(double));

    /* First-order lower-boundary fields that seed higher scattering orders. */
    double *order1_I = (double *)calloc(buf_size, sizeof(double));
    double *order1_Q = (double *)calloc(buf_size, sizeof(double));
    double *order1_U = (double *)calloc(buf_size, sizeof(double));
    double *surf_seed_I = use_surf_first_order
        ? (double *)calloc(buf_size, sizeof(double)) : NULL;
    double *surf_seed_Q = use_surf_first_order
        ? (double *)calloc(buf_size, sizeof(double)) : NULL;
    double *surf_seed_U = use_surf_first_order
        ? (double *)calloc(buf_size, sizeof(double)) : NULL;
    if (!order1_I || !order1_Q || !order1_U ||
        (use_surf_first_order &&
         (!surf_seed_I || !surf_seed_Q || !surf_seed_U))) {
        free(order1_I); free(order1_Q); free(order1_U);
        free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U);
        rt_surface_boundary_free(&surface_boundary);
        return -1;
    }

    surface_status = rt_surface_build_initial_sources(
        &surface_boundary, atm, ws, m, opts, surface_kind, n_water,
        wind_speed, sigma_type, q_convention,
        total_I, total_Q, total_U,
        order1_I, order1_Q, order1_U,
        surf_seed_I, surf_seed_Q, surf_seed_U);
    if (surface_status != 0) {
        free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U);
        free(order1_I); free(order1_Q); free(order1_U);
        rt_surface_boundary_free(&surface_boundary);
        return surface_status;
    }

    if (opts->max_iterations == 1) {
        if (use_surf_first_order) {
            rt_surface_remove_direct_seed(atm, surf_seed_I, surf_seed_Q,
                                          surf_seed_U,
                                          total_I, total_Q, total_U);
        }
        result->n_orders_used = 1;
        result->converged = 1;
        result->final_residual = 0.0;
        free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U);
        free(order1_I); free(order1_Q); free(order1_U);
        rt_surface_boundary_free(&surface_boundary);
        return 0;
    }

    double *I_prev = calloc(buf_size, sizeof(double));
    double *Q_prev = calloc(buf_size, sizeof(double));
    double *U_prev = calloc(buf_size, sizeof(double));
    double *J_I    = calloc(buf_size, sizeof(double));
    double *J_Q    = calloc(buf_size, sizeof(double));
    double *J_U    = calloc(buf_size, sizeof(double));
    double *I_curr = calloc(buf_size, sizeof(double));
    double *Q_curr = calloc(buf_size, sizeof(double));
    double *U_curr = calloc(buf_size, sizeof(double));
    double *bc_I   = calloc((size_t)n_mu, sizeof(double));
    double *bc_Q   = calloc((size_t)n_mu, sizeof(double));
    double *bc_U   = calloc((size_t)n_mu, sizeof(double));
    if (!I_prev || !Q_prev || !U_prev || !J_I || !J_Q || !J_U ||
        !I_curr || !Q_curr || !U_curr || !bc_I || !bc_Q || !bc_U) {
        free(I_prev); free(Q_prev); free(U_prev);
        free(J_I); free(J_Q); free(J_U);
        free(I_curr); free(Q_curr); free(U_curr);
        free(bc_I); free(bc_Q); free(bc_U);
        free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U); free(order1_I); free(order1_Q); free(order1_U);
        rt_surface_boundary_free(&surface_boundary);
        return -1;
    }
    /* Phase B (T0-1, 2026-05-10): I_prev seed = primary + order1_eff.
     * order1_eff buffer 가 FRESNEL_DIFF (flat path) 또는 RII/RQQ/RUU (rough sea
     * mer agitée) 의 surface first-order contribution 누적.  iteration 의
     * source J = K · I_prev 가 *이 first-order field* 를 사용하면 *higher-order
     * coupling (sun → wave/Fresnel reflection → atm-scatter → atm-scatter → ...)*
     * 정확히 처리됨.  이전 V3: I_prev = primary 만 (atm-only) 이라 *surface
     * higher-order* 누락. */
    for (size_t i = 0; i < buf_size; i++) {
        I_prev[i] = primary_I[i] + order1_I[i];
        Q_prev[i] = primary_Q[i] + order1_Q[i];
        U_prev[i] = primary_U[i] + order1_U[i];
    }

    int n = 1;
    double residual = 0.0;
    int converged = 0;

    for (int it = 2; it <= opts->max_iterations; it++) {
        int rc = rt_sos_operator_apply_vector(atm, m, ws, I_prev, Q_prev, U_prev,
                                       J_I, J_Q, J_U);
        if (rc != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr);
            free(bc_I); free(bc_Q); free(bc_U);
            free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U); free(order1_I); free(order1_Q); free(order1_U);
            rt_surface_boundary_free(&surface_boundary);
            return rc;
        }

        rc = rt_surface_boundary_apply(&surface_boundary, atm, n_water,
                                       q_convention,
                                       I_prev, Q_prev, U_prev,
                                       bc_I, bc_Q, bc_U);
        if (rc != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr);
            free(bc_I); free(bc_Q); free(bc_U);
            free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U);
            free(order1_I); free(order1_Q); free(order1_U);
            rt_surface_boundary_free(&surface_boundary);
            return rc;
        }

        if (rt_solver_integrate_with_surface_bc(atm, m, J_I, bc_I,
                RT_INTEGRATION_METHOD_LINEAR, I_curr) != 0 ||
            rt_solver_integrate_with_surface_bc(atm, m, J_Q, bc_Q,
                RT_INTEGRATION_METHOD_LINEAR, Q_curr) != 0 ||
            rt_solver_integrate_with_surface_bc(atm, m, J_U, bc_U,
                RT_INTEGRATION_METHOD_LINEAR, U_curr) != 0) {
            free(I_prev); free(Q_prev); free(U_prev);
            free(J_I); free(J_Q); free(J_U);
            free(I_curr); free(Q_curr); free(U_curr);
            free(bc_I); free(bc_Q); free(bc_U);
            free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U); free(order1_I); free(order1_Q); free(order1_U);
            rt_surface_boundary_free(&surface_boundary);
            return -2;
        }

        for (int k = 0; k <= nt; k++) {
            for (int j = -n_mu; j <= n_mu; j++) {
                if (j == 0) continue;
                size_t idx = (size_t)k * (size_t)dirs + (size_t)(j + n_mu);
                total_I[idx] += I_curr[idx];
                total_Q[idx] += Q_curr[idx];
                total_U[idx] += U_curr[idx];
            }
        }
        n = it;

        double max_in_I = field_max_abs(I_curr, nt, n_mu);
        double max_in_Q = field_max_abs(Q_curr, nt, n_mu);
        double max_in_U = field_max_abs(U_curr, nt, n_mu);
        double max_in   = max_in_I;
        if (max_in_Q > max_in) max_in = max_in_Q;
        if (max_in_U > max_in) max_in = max_in_U;

        double max_t_I = field_max_abs(total_I, nt, n_mu);
        double max_t_Q = field_max_abs(total_Q, nt, n_mu);
        double max_t_U = field_max_abs(total_U, nt, n_mu);
        double max_t   = max_t_I;
        if (max_t_Q > max_t) max_t = max_t_Q;
        if (max_t_U > max_t) max_t = max_t_U;

        residual = (max_t > 0.0) ? max_in / max_t : 0.0;
        if (residual < opts->tolerance) {
            converged = 1;
            break;
        }

        memcpy(I_prev, I_curr, buf_size * sizeof(double));
        memcpy(Q_prev, Q_curr, buf_size * sizeof(double));
        memcpy(U_prev, U_curr, buf_size * sizeof(double));
    }

    if (use_surf_first_order) {
        rt_surface_remove_direct_seed(atm, surf_seed_I, surf_seed_Q,
                                      surf_seed_U,
                                      total_I, total_Q, total_U);
    }

    free(surf_seed_I); free(surf_seed_Q); free(surf_seed_U);
    free(I_prev); free(Q_prev); free(U_prev);
    free(J_I); free(J_Q); free(J_U);
    free(I_curr); free(Q_curr); free(U_curr);
    free(bc_I); free(bc_Q); free(bc_U);
    free(order1_I); free(order1_Q); free(order1_U);
    rt_surface_boundary_free(&surface_boundary);

    result->n_orders_used  = n;
    result->converged      = converged;
    result->final_residual = residual;
    return 0;
}

/* Public compatibility wrappers around the OCRT Fourier-convention module. */
double rt_solver_reconstruct_phi(const double *I_total_per_m_at_view,
                                 int m_max,
                                 double delta_phi_rad)
{
    return rt_fourier_reconstruct_cos(I_total_per_m_at_view, m_max,
                                      delta_phi_rad);
}

double rt_solver_reconstruct_phi_sin(const double *U_total_per_m_at_view,
                                     int m_max,
                                     double delta_phi_rad)
{
    return rt_fourier_reconstruct_sin(U_total_per_m_at_view, m_max,
                                      delta_phi_rad);
}

/* ============================================================================
 * Phase B.4 (Stage 1) — rt_solve_case_ocean
 *
 * Drives the in-water SOS solver for RT_SURFACE_OCEAN cases.
 * Stage 1 (this implementation): in-water solver only, no atmospheric coupling.
 *   - F_sun (cs->F_sun, default π) is treated as just-above-surface BOA solar
 *     irradiance directly (T_aw=1 assumption).
 *   - Ed_0plus_air = F_sun · μ_sun_air  (no atmospheric absorption/scattering)
 *   - In-water source receives the same F_sun normalization via rt_water_rt_sos_pure.
 *
 * Stage 2 (next turn): atmospheric SOS BOA Stokes → T_aw Mueller → in-water
 * boundary condition full coupling.
 * ============================================================================ */

#include "rt_water_rt.h"
#include "rt_air_water.h"
#include "rt_air_water_coupling.h"
#include "rt_quadrature.h"
#include "rt_raa_convention.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Scalar-FF fixed-bulk IOP table reader.
 * Accepts either whitespace/CSV with no header:
 *   wavelength_nm  a_total  b_total
 * or a header containing wavelength_nm/lambda_nm/wl_nm, a/a_total, b/b_total.
 * Interpolation is strict inside the table range; no endpoint extrapolation. */
typedef struct { double wl, a, b; } scalar_ff_iop_row_t;
static int scalar_ff_iop_cmp(const void *a, const void *b) {
    const scalar_ff_iop_row_t *ra=(const scalar_ff_iop_row_t*)a, *rb=(const scalar_ff_iop_row_t*)b;
    return (ra->wl < rb->wl) ? -1 : (ra->wl > rb->wl);
}
static int scalar_ff_split_line(char *line, char **fields, int maxf) {
    int n=0; char *p=line;
    while (*p && n<maxf) {
        while (*p==' '||*p=='\t'||*p==','||*p=='\r'||*p=='\n') p++;
        if (!*p || *p=='#') break;
        fields[n++]=p;
        while (*p && *p!=' '&&*p!='\t'&&*p!=','&&*p!='\r'&&*p!='\n') p++;
        if (*p) *p++='\0';
    }
    return n;
}
static void scalar_ff_lower(char *s) { for(; s && *s; ++s) *s=(char)tolower((unsigned char)*s); }
static int scalar_ff_col_any(char **fields, int n, const char **names) {
    for(int i=0;i<n;i++){ char tmp[128]; snprintf(tmp,sizeof tmp,"%s", fields[i]?fields[i]:""); scalar_ff_lower(tmp); for(int k=0; names[k]; ++k) if(!strcmp(tmp,names[k])) return i; }
    return -1;
}
static int scalar_ff_is_numeric(const char *s) { if(!s||!*s) return 0; char *e=NULL; (void)strtod(s,&e); return e && *e=='\0'; }
static int scalar_ff_iop_table_lookup(const char *path, double wavelength_nm, double *a_out, double *b_out) {
    if(!path||!*path||!a_out||!b_out) return -1;
    FILE *fp=fopen(path,"r"); if(!fp) return -2;
    scalar_ff_iop_row_t *rows=NULL; int n=0, cap=0; int iw=0, ia=1, ib=2, header_done=0;
    char line[8192]; const char *wl_names[]={"wavelength_nm","lambda_nm","wl_nm","wavelength","lambda","wl",NULL};
    const char *a_names[]={"a_total_m-1","a_total_m_inv","a_total","a_m-1","a_m_inv","a",NULL};
    const char *b_names[]={"b_total_m-1","b_total_m_inv","b_total","b_m-1","b_m_inv","b",NULL};
    while(fgets(line,sizeof line,fp)){
        char *p=line; while(*p==' '||*p=='\t') p++; if(!*p||*p=='#'||*p=='\n'||*p=='\r') continue;
        char *fields[64]; int nf=scalar_ff_split_line(p,fields,64); if(nf<3) continue;
        if(!header_done && !scalar_ff_is_numeric(fields[0])){
            iw=scalar_ff_col_any(fields,nf,wl_names); ia=scalar_ff_col_any(fields,nf,a_names); ib=scalar_ff_col_any(fields,nf,b_names);
            if(iw<0||ia<0||ib<0){ fclose(fp); free(rows); return -3; }
            header_done=1; continue;
        }
        header_done=1; if(nf<=iw||nf<=ia||nf<=ib) continue;
        double wl=atof(fields[iw]), aa=atof(fields[ia]), bb=atof(fields[ib]);
        if(!(wl>0.0)||!(aa>=0.0)||!(bb>=0.0)) continue;
        if(n>=cap){ cap=cap?cap*2:64; scalar_ff_iop_row_t *tmp=(scalar_ff_iop_row_t*)realloc(rows,(size_t)cap*sizeof(*rows)); if(!tmp){ fclose(fp); free(rows); return -4; } rows=tmp; }
        rows[n++]=(scalar_ff_iop_row_t){wl,aa,bb};
    }
    fclose(fp); if(n<1){ free(rows); return -5; }
    qsort(rows,(size_t)n,sizeof(*rows),scalar_ff_iop_cmp);
    if(wavelength_nm < rows[0].wl-1e-9 || wavelength_nm > rows[n-1].wl+1e-9){ free(rows); return -6; }
    if(n==1 || fabs(wavelength_nm-rows[0].wl)<1e-9){ *a_out=rows[0].a; *b_out=rows[0].b; free(rows); return 0; }
    for(int i=0;i<n-1;i++){
        if(fabs(wavelength_nm-rows[i].wl)<1e-9){ *a_out=rows[i].a; *b_out=rows[i].b; free(rows); return 0; }
        if(rows[i].wl <= wavelength_nm && wavelength_nm <= rows[i+1].wl){ double f=(wavelength_nm-rows[i].wl)/(rows[i+1].wl-rows[i].wl); *a_out=rows[i].a*(1.0-f)+rows[i+1].a*f; *b_out=rows[i].b*(1.0-f)+rows[i+1].b*f; free(rows); return 0; }
    }
    free(rows); return -7;
}

#include "rt_angle_grid.h" /* v1.10 B-0a.3 */
/* v-diag S6 trace (env OCRT_S6_TRACE=1): per-row block wall-time budget. */
#include <time.h>
static double ocrt_s6_now(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }
static int ocrt_s6_trace_on(void){ static int v=-1; if(v<0){const char*e=getenv("OCRT_S6_TRACE"); v=(e&&*e&&strcmp(e,"0"))?1:0;} return v; }
#define OCRT_S6T(...) do{ if(ocrt_s6_trace_on()) fprintf(stderr, __VA_ARGS__); }while(0)

/* 2026-07-14 S7-D FIX (Jae 승인): S7/S7b 재생은 노드장 재구성이므로 뷰
 * μ가 최대 GL 노드를 넘는 행(연직 부근; nmw=96 기준 vza < 약 1.4도)은
 * 범위 밖 선형 외삽이 되어 m>=1 성분이 소거되지 않는다(실측: 연직 TOA_I가
 * raa 0/180 에서 ±1.7% 갈라짐 — 물리적으로 raa 불변이어야 함).  해당 행만
 * 캐시를 우회해 행별 표준 계산(뷰 0-가중 노드 정확 평가)으로 폴백한다.
 * 비용: 격자당 대기 solve 1회 수준.  GL 최대 노드는 n별 TLS 캐시. */
static int ocrt_s7_view_within_nodes(int n_mu, double mu_view) {
    static _Thread_local int    c_n = -1;
    static _Thread_local double c_mu_max = 0.0;
    if (c_n != n_mu) {
        double mu[512], wg[512];
        if (n_mu < 1 || n_mu > 512 ||
            rt_quadrature_gauss_legendre_pos(n_mu, mu, wg) != 0) return 0;
        c_mu_max = (mu[n_mu - 1] > mu[0]) ? mu[n_mu - 1] : mu[0];
        c_n = n_mu;
    }
    return mu_view <= c_mu_max;
}


/* v1.2 LUT-PERF S17a (2026-07-21): exact near-nadir RAA folding.
 *
 * S7 deliberately refuses to extrapolate the Gauss-node LUT beyond its
 * largest mu node.  That protects exact/near-nadir rows, but the old fallback
 * then repeated the complete atmosphere pass for every RAA.  The atmosphere
 * solve is RAA-independent up to the final Fourier reconstruction, so retain
 * the exact view-as-node per-m samples and the BOA Fourier field from the first
 * RAA and replay only the reconstruction for the remaining RAA values.
 *
 * This is not an angular interpolation and therefore preserves the S7-D
 * accuracy rule.  OCRT_S17A_OFF=1 disables the optimization. */
static unsigned long long ocrt_fnv1a_bytes(unsigned long long h,
                                            const void *ptr, size_t n)
{
    const unsigned char *p = (const unsigned char*)ptr;
    for (size_t i = 0; i < n; ++i) h = (h ^ (unsigned long long)p[i]) * 1099511628211ULL;
    return h;
}

static unsigned long long ocrt_s17a_key(const rt_case_t *cs,
                                         const rt_options_t *opts,
                                         const rt_aerosol_input_t *aer)
{
    unsigned long long h = 1469598103934665603ULL;
    rt_case_t kc = *cs;
    rt_options_t ko = *opts;
    /* RAA is the folded dimension.  Export/result pointers are not physics. */
    kc.raa_deg = 0.0;
    ko.ext_bottom_per_m_I = ko.ext_bottom_per_m_Q = ko.ext_bottom_per_m_U = NULL;
    ko.ext_bottom_n_mu = 0; ko.ext_bottom_m_max = -1; ko.bottom_source_only = 0;
    ko.toa_view_per_m_I = ko.toa_view_per_m_Q = ko.toa_view_per_m_U = NULL;
    h = ocrt_fnv1a_bytes(h, &kc, sizeof(kc));
    h = ocrt_fnv1a_bytes(h, &ko, sizeof(ko));
    if (!aer) {
        const unsigned long long z = 0;
        return ocrt_fnv1a_bytes(h, &z, sizeof(z));
    }
    h = ocrt_fnv1a_bytes(h, &aer->tau_a, sizeof(aer->tau_a));
    h = ocrt_fnv1a_bytes(h, &aer->ssa_a, sizeof(aer->ssa_a));
    h = ocrt_fnv1a_bytes(h, &aer->L_max, sizeof(aer->L_max));
    const size_t nl = (size_t)(aer->L_max + 1);
    h = ocrt_fnv1a_bytes(h, aer->betal_aer,  nl * sizeof(double));
    h = ocrt_fnv1a_bytes(h, aer->gammal_aer, nl * sizeof(double));
    h = ocrt_fnv1a_bytes(h, aer->alphal_aer, nl * sizeof(double));
    h = ocrt_fnv1a_bytes(h, aer->zetal_aer,  nl * sizeof(double));
    h = ocrt_fnv1a_bytes(h, &aer->use_value_kernel, sizeof(aer->use_value_kernel));
    h = ocrt_fnv1a_bytes(h, &aer->n_ang_phase, sizeof(aer->n_ang_phase));
    if (aer->use_value_kernel && aer->n_ang_phase > 0) {
        const size_t na = (size_t)aer->n_ang_phase;
        h = ocrt_fnv1a_bytes(h, aer->theta_phase, na * sizeof(double));
        h = ocrt_fnv1a_bytes(h, aer->P11_phase,   na * sizeof(double));
        h = ocrt_fnv1a_bytes(h, aer->P12_phase,   na * sizeof(double));
        h = ocrt_fnv1a_bytes(h, aer->P33_phase,   na * sizeof(double));
    }
    return h;
}

static int ocrt_boa_copy(rt_atm_boa_export_t *dst,
                          const rt_atm_boa_export_t *src)
{
    if (!dst || !src || src->n_mu < 1 || src->m_max < 0) return -1;
    if (!dst->allocated || dst->n_mu != src->n_mu || dst->m_max != src->m_max) {
        rt_atm_boa_export_free(dst);
        if (rt_atm_boa_export_alloc(dst, src->m_max, src->n_mu) != 0) return -1;
    }
    const size_t n = (size_t)(src->m_max + 1) * (size_t)src->n_mu;
    memcpy(dst->mu_quad_pos, src->mu_quad_pos, (size_t)src->n_mu * sizeof(double));
    memcpy(dst->I_per_m, src->I_per_m, n * sizeof(double));
    memcpy(dst->Q_per_m, src->Q_per_m, n * sizeof(double));
    memcpy(dst->U_per_m, src->U_per_m, n * sizeof(double));
    dst->tau_atm_total = src->tau_atm_total;
    dst->mu_sun_air = src->mu_sun_air;
    return 0;
}

static double ocrt_boa_T_diff_dn_dir(const rt_atm_boa_export_t *be,
                                      double mu_view, double raa_deg,
                                      double mu_solar)
{
    if (!be || be->n_mu < 2 || be->m_max < 0 || be->m_max > 32 ||
        !(mu_solar > 0.0)) return 0.0;
    double pm[33] = {0};
    for (int m = 0; m <= be->m_max; ++m) {
        const double *Im = &be->I_per_m[(size_t)m * (size_t)be->n_mu];
        int j_hi = 1;
        while (j_hi < be->n_mu && be->mu_quad_pos[j_hi - 1] < mu_view) ++j_hi;
        int j_lo = j_hi - 1;
        if (j_lo < 1) { j_lo = 1; j_hi = 2; }
        if (j_hi > be->n_mu) { j_hi = be->n_mu; j_lo = be->n_mu - 1; }
        const double x0 = be->mu_quad_pos[j_lo - 1];
        const double x1 = be->mu_quad_pos[j_hi - 1];
        const double t = (mu_view - x0) / (x1 - x0);
        pm[m] = (1.0 - t) * Im[j_lo - 1] + t * Im[j_hi - 1];
    }
    const double phi = rt_raa_to_atm_fourier_phi(raa_deg);
    const double I = rt_solver_reconstruct_phi(pm, be->m_max, phi);
    return (2.0 * M_PI * mu_view * I) / (M_PI * mu_solar);
}

/* v1.11 (2026-09-05) coupled-ocean scope guard for --pssa.
 *
 * IPSS Eq. (7) is validated for the atmosphere-only / black-ocean TOA path
 * (Phase I).  Applying it to a coupled ocean-atmosphere solve is an OPEN
 * design decision — "single scattering" and "total field" both need new
 * definitions once the water body is in the solve — so the combination is
 * refused loudly.  Before v1.11 this function silently fell back to the
 * legacy average-secant Chapman beam, which is exactly the wrong-entry
 * behaviour the IPSS replacement was meant to remove; that fallback and the
 * legacy algorithm itself are gone.
 *
 * Returns -1.0 (refusal) whenever --pssa is on; the caller aborts the solve. */
static double ocrt_pssa_T_dir_dn_for_ocean(const rt_options_t *opts)
{
    if (opts && opts->pssa) {
        static int warned_once = 0;
        if (!warned_once) {
            warned_once = 1;
            fprintf(stderr,
                "error: --pssa (IPSS) is validated for the atmosphere-only / "
                "black-ocean path only (Phase I).\n"
                "       Drop --pssa for a coupled --water-model run.  The "
                "legacy Chapman fallback was removed in v1.11.\n");
        }
        return -1.0;
    }
    return 1.0;
}

/* Recompute Kd(0-) when a diffuse-top water boundary is active.
 *
 * The water SOS includes scattering generated by the atmospheric skylight,
 * but the unscattered diffuse-top radiance is a boundary term and is not part
 * of the internal SOS field.  The historical finite-difference Kd therefore
 * omitted this term at both z=0 and the first water level.  Re-add it from the
 * m=0 angular boundary field before evaluating the logarithmic slope.
 *
 * sky_I_m0 is on the common F_sun_water-normalized scale used by the water
 * solve (the orchestrator pre-scales the atmospheric snapshot accordingly).
 * The 2*pi*(F_sun_water/pi) factor converts the quadrature sum to physical
 * irradiance. */
static double ocrt_kd_with_unscattered_sky(
        const rt_water_rt_result_t *wr,
        const double *sky_I_m0,
        const double *mu_water,
        const double *w_water,
        int n_mu_water,
        double F_sun_water,
        double public_Ed_sky0,
        double fallback_kd) {
    if (!wr || !sky_I_m0 || !mu_water || !w_water ||
        n_mu_water <= 0 || !(F_sun_water > 0.0) ||
        !(wr->Ed_level1_water > 0.0) ||
        !(wr->tau_level1_used > 0.0) || !(wr->z_level1_used > 0.0))
        return fallback_kd;

    double q0 = 0.0, q1 = 0.0;
    for (int c = 0; c < n_mu_water; ++c) {
        const double mu = mu_water[c];
        const double wt = w_water[c];
        const double I0 = sky_I_m0[c];
        if (!(mu > 0.0) || !(wt > 0.0) || !isfinite(I0)) continue;
        const double term = I0 * mu * wt;
        q0 += term;
        q1 += term * exp(-wr->tau_level1_used / mu);
    }
    if (!(q0 > 0.0) || !(q1 > 0.0)) return fallback_kd;

    const double scale = 2.0 * M_PI * (F_sun_water / M_PI);
    const double Ed_sky0 = scale * q0;
    const double Ed_sky1 = scale * q1;
    const double Ed0 = wr->Ed_0minus_water + Ed_sky0;
    const double Ed1 = wr->Ed_level1_water + Ed_sky1;
    if (!(Ed0 > 0.0) || !(Ed1 > 0.0)) return fallback_kd;

    const double kd = -log(Ed1 / Ed0) / wr->z_level1_used;
    if (!isfinite(kd)) return fallback_kd;

    if (getenv("OCRT_DUMP_KD_SKY_FIX")) {
        fprintf(stderr,
            "KD_SKY_FIX Ed0_field=%.12e Ed0_sky_quad=%.12e "
            "Ed0_sky_public=%.12e Ed1_field=%.12e Ed1_sky=%.12e "
            "z1=%.12e Kd_old=%.12e Kd_new=%.12e\n",
            wr->Ed_0minus_water, Ed_sky0, public_Ed_sky0,
            wr->Ed_level1_water, Ed_sky1, wr->z_level1_used,
            fallback_kd, kd);
    }
    return kd;
}

int rt_solve_case_ocean(const rt_case_t *cs, const rt_options_t *opts,
                         const rt_water_iop_lut_t *aw_lut,
                         const rt_water_iop_psi_T_lut_t *psi_T_lut,
                         const rt_aerosol_input_t *aer,
                         rt_result_t *out) {
    const double ocrt_s6_t0 = ocrt_s6_now();
    OCRT_S6T("[S6T] row vza=%.3f raa=%.1f enter abs=%.3f\n", cs->vza_deg, cs->raa_deg, ocrt_s6_now());
    int ocrt_s6_solved = 0, ocrt_s6_skipped = 0; (void)ocrt_s6_solved; (void)ocrt_s6_skipped;
    if (!cs || !opts || !aw_lut || !out) return -1;
    if (cs->surface != RT_SURFACE_OCEAN) {
        fprintf(stderr, "[B.4] rt_solve_case_ocean: surface != OCEAN (got %d)\n",
                (int)cs->surface);
        return -1;
    }

    memset(out, 0, sizeof(*out));

    /* Configure in-water SOS options. The in-water options have their own
     * field names (suffix _water) distinct from atmospheric rt_options_t.
     * Stage 1 uses defaults; future expansion of rt_options_t can pipe through
     * water-specific knobs (n_layers_water, tau_max_target, etc.). */
    rt_water_rt_options_t w_opts = rt_water_rt_options_default();
    w_opts.water_input_mode = cs->water_input_mode;
    w_opts.q_convention = cs->q_convention;
    w_opts.view_as_node = cs->water_view_as_node;
    w_opts.view_vza_deg_list = cs->water_view_vza_list;   /* #20 */
    w_opts.n_view_vza        = cs->n_water_view_vza;
    w_opts.water_phase_kernel = cs->water_phase_kernel;  /* common: value(0)/moment(1) for all in-water paths (fixed-bulk + CCRR) */
    /* Cox-Munk water->air transmission branch (option A). We are inside the
     * ocean solver, so always forward the wind and slope-variance model; the
     * water RT picks flat T_wa for wind_speed<=0 and surface_T_coxmunk_trig
     * facet integration for wind_speed>0 (matches the atmosphere-side branch,
     * which also goes Cox-Munk whenever wind_speed>0 regardless of OCEAN vs
     * COXMUNK surface enum). */
    w_opts.wind_speed = cs->wind_speed;
    w_opts.cox_munk_sigma_type = cs->sigma_type;
    /* In-water depth/grid controls are production-selected below. Debug-only
     * overrides are applied after mode-specific depth policy. */
    if (cs->n_mu_water_override > 0) {
        w_opts.n_mu_water = cs->n_mu_water_override;
    }
    /* B.5 (2026-05-23): CDOM 2-parameter absorption.
     *   cs->a_cdom_440_m_inv <= 0  → CDOM 비활성 (pre-B.5 baseline 유지).
     *   cs->S_cdom_nm_inv    <= 0  → option default (0.014 nm⁻¹) 사용.
     *   cs->cdom_ref_lambda_nm <= 0 → option default (440 nm) 사용. */
    if (cs->a_cdom_440_m_inv > 0.0) {
        w_opts.a_cdom_440_m_inv = cs->a_cdom_440_m_inv;
        if (cs->S_cdom_nm_inv > 0.0)
            w_opts.S_cdom_nm_inv = cs->S_cdom_nm_inv;
        if (cs->cdom_ref_lambda_nm > 0.0)
            w_opts.cdom_ref_lambda_nm = cs->cdom_ref_lambda_nm;
    }
    if (cs->ccrr_mode) {
        w_opts.ccrr_mode = 1;
        w_opts.water_constituent_model = cs->water_constituent_model;
        w_opts.ccrr_chl_mg_m3 = cs->ccrr_chl_mg_m3;
        w_opts.organic_phyto_scattering = cs->organic_phyto_scattering;
        w_opts.ccrr_min_g_m3 = cs->ccrr_min_g_m3;
        w_opts.tsm_species = cs->tsm_species;
        w_opts.organic_phyto_group = cs->organic_phyto_group;
        w_opts.detritus_a440_m_inv = cs->detritus_a440_m_inv;
        w_opts.detritus_slope_nm_inv = cs->detritus_slope_nm_inv;
        w_opts.ccrr_phase_moments_path = cs->ccrr_phase_moments_path;
        w_opts.ccrr_particle_phase_lut_path = cs->ccrr_particle_phase_lut_path;
        w_opts.ccrr_particle_phase_case_id = cs->ccrr_particle_phase_case_id;
        w_opts.ccrr_particle_phase_wavelength_nm = cs->ccrr_particle_phase_wavelength_nm;
        w_opts.ccrr_particle_phase_lmax = cs->ccrr_particle_phase_lmax;
        w_opts.ccrr_particle_phase_nphi = cs->ccrr_particle_phase_nphi;
        w_opts.water_phase_kernel = cs->water_phase_kernel;
        w_opts.water_mie_phase_path = cs->water_mie_phase_path;
        w_opts.water_mie_moment_mode = cs->water_mie_moment_mode;
        w_opts.water_mie_truncation_mode = cs->water_mie_truncation_mode;
        w_opts.water_mie_ss_mode = cs->water_mie_ss_mode;
        w_opts.water_mie_moment_n_mu = cs->water_mie_moment_n_mu;
        w_opts.fixed_bulk_lmax = cs->fixed_bulk_lmax;

        /* CCRR constituent mode is intentionally pre-solver only:
         * it maps CHL/MIN/aDOM to OCRT water IOPs and scalar particle phase
         * moments.  The CCRR v0.5.0 reference bookkeeping is direct-only with
         * a black lower boundary at optical depth 15 and continuous view
         * reconstruction.  We therefore set only the boundary/bookkeeping
         * knobs needed for row-by-row comparison; the RT solution is still the
         * OCRT in-water SOS solver, not the CCRR solver.  Public depth/order
         * overrides are disabled; debug-only overrides remain behind OCRT_DEBUG. */
        w_opts.tau_max_target = 15.0;  /* pure/reference baseline; particle-rich policy may raise */
        w_opts.n_layers_water = 300;   /* Δτ≈0.05 for τ=15; particle-rich policy may raise */
        if ((cs->ccrr_chl_mg_m3 > 0.0 || cs->ccrr_min_g_m3 > 0.0) &&
            w_opts.m_max_water < 24) w_opts.m_max_water = 24;
        if (w_opts.max_iterations < 140) w_opts.max_iterations = 140;
        /* CCRR/IOCCG21 reference bookkeeping uses continuous reconstruction.
         * Do not force target-node extraction in CCRR mode; users may still
         * request the diagnostic exact-view node explicitly with
         * --water-view-as-node. */
        w_opts.view_as_node = cs->water_view_as_node;
    w_opts.view_vza_deg_list = cs->water_view_vza_list;   /* #20 */
    w_opts.n_view_vza        = cs->n_water_view_vza;
    }
    if (cs->fixed_bulk_iop_mode) {
        w_opts.fixed_bulk_iop_mode = 1;

        if (cs->scalar_ff_iop_path && *cs->scalar_ff_iop_path) {
            double a_tab = 0.0, b_tab = 0.0;
            int lrc = scalar_ff_iop_table_lookup(cs->scalar_ff_iop_path,
                                                 cs->wavelength_nm,
                                                 &a_tab, &b_tab);
            if (lrc != 0) {
                fprintf(stderr,
                    "[scalar-ff] failed to read/interpolate IOP table '%s' "
                    "at wavelength %.10g nm (rc=%d). Table must cover the "
                    "requested wavelength without endpoint extrapolation.\n",
                    cs->scalar_ff_iop_path, cs->wavelength_nm, lrc);
                return -1;
            }
            if (!(cs->scalar_ff_bb_over_b > 0.0 && cs->scalar_ff_bb_over_b < 0.5)) {
                fprintf(stderr,
                    "[scalar-ff] --scalar-ff-bbfrac must be in (0,0.5); got %.12g\n",
                    cs->scalar_ff_bb_over_b);
                return -1;
            }
            w_opts.fixed_a_total_m_inv  = a_tab;
            w_opts.fixed_b_total_m_inv  = b_tab;
            w_opts.fixed_bb_total_m_inv = b_tab * cs->scalar_ff_bb_over_b;
            w_opts.fixed_bulk_phase_model = 4; /* analytic scalar FF + formal delta-M */
            /* L=10 is the current physically formal retained-order candidate
             * for the FF fixed-bulk validation path. Users may override with
             * --scalar-ff-lmax / --fixed-bulk-lmax. */
            w_opts.fixed_bulk_lmax = (cs->fixed_bulk_lmax > 0 && cs->fixed_bulk_lmax != 30) ?
                                      cs->fixed_bulk_lmax : 10;
            w_opts.fixed_bulk_ff_n = (cs->scalar_ff_refractive_index > 1.0) ?
                                      cs->scalar_ff_refractive_index : 1.18;
            w_opts.fixed_bulk_ff_mu = cs->scalar_ff_mu;
        } else {
            w_opts.fixed_a_total_m_inv = cs->fixed_a_total_m_inv;
            w_opts.fixed_b_total_m_inv = cs->fixed_b_total_m_inv;
            w_opts.fixed_bb_total_m_inv = cs->fixed_bb_total_m_inv;
            w_opts.fixed_bulk_lmax = cs->fixed_bulk_lmax;
            w_opts.fixed_bulk_phase_model = cs->fixed_bulk_phase_model;
            w_opts.fixed_bulk_ff_n = cs->fixed_bulk_ff_n;
            w_opts.fixed_bulk_ff_mu = cs->fixed_bulk_ff_mu;
        }

        w_opts.fixed_bulk_phase_nphi = cs->fixed_bulk_phase_nphi;
        w_opts.fixed_bulk_phase_lut_path = cs->fixed_bulk_phase_lut_path;
        w_opts.water_mie_phase_path = cs->water_mie_phase_path;
        w_opts.water_mie_moment_mode = cs->water_mie_moment_mode;
        w_opts.water_mie_truncation_mode = cs->water_mie_truncation_mode;
        w_opts.water_mie_ss_mode = cs->water_mie_ss_mode;
        w_opts.water_mie_moment_n_mu = cs->water_mie_moment_n_mu;
        w_opts.fixed_bulk_phase_case_id = cs->fixed_bulk_phase_case_id;
        w_opts.fixed_bulk_phase_wavelength_nm = cs->fixed_bulk_phase_wavelength_nm;
        w_opts.tau_max_target = 15.0;  /* baseline; particle-rich depth policy may raise */
        w_opts.n_layers_water = 300;   /* Δτ≈0.05 for τ=15; policy may raise */
        w_opts.m_max_water = (w_opts.fixed_bulk_lmax > 0) ? w_opts.fixed_bulk_lmax : 30;
        if (w_opts.max_iterations < 160) w_opts.max_iterations = 160;
        w_opts.view_as_node = cs->water_view_as_node;
    w_opts.view_vza_deg_list = cs->water_view_vza_list;   /* #20 */
    w_opts.n_view_vza        = cs->n_water_view_vza;
    }
    /* Particle-rich water media need safe SOS and depth defaults.  The SOS cap
     * is a minimum only and never lowers stronger mode-specific caps such as
     * fixed-bulk 160.  The vertical slab is also
     * deepened: total τ=15 can be too shallow in forward-peaked high-TSM media
     * when judged by transport/diffusion depth.
     *
     * Trigger policy:
     *   - constituent mode: any Chl or mineral particle concentration > 0
     *   - fixed-bulk mode: total scattering/backscattering above pure-water scale
     */
    {
        int high_particle_medium = 0;
        int high_particle_cap = (cs->debug_high_particle_water_max_orders > 0)
                                ? cs->debug_high_particle_water_max_orders
                                : OCRT_WATER_HIGH_PARTICLE_SOS_CAP_DEFAULT;
        if (cs->ccrr_mode && (cs->ccrr_chl_mg_m3 > 0.0 || cs->ccrr_min_g_m3 > 0.0))
            high_particle_medium = 1;
        if (cs->fixed_bulk_iop_mode &&
            (cs->fixed_b_total_m_inv > 2.0e-2 || cs->fixed_bb_total_m_inv > 5.0e-4))
            high_particle_medium = 1;
        if (high_particle_medium) {
            if (w_opts.max_iterations < high_particle_cap)
                w_opts.max_iterations = high_particle_cap;
            ocrt_water_apply_deep_depth_policy(&w_opts);
        }
    }

    if (cs->water_m_max_override > 0) {
        w_opts.m_max_water = cs->water_m_max_override;
    }

    /* opts->n_mu / opts->n_layers refer to ATMOSPHERIC quadrature/layers; the
     * in-water grid is independent.  Production depth/grid controls are not
     * public CLI.  Debug overrides require OCRT_DEBUG=1 and may intentionally
     * lower accuracy for sensitivity tests. */
    if (cs->tau_max_target_override > 0.0) {
        w_opts.tau_max_target = cs->tau_max_target_override;
        if (cs->n_layers_water_override <= 0)
            w_opts.n_layers_water = ocrt_water_layers_for_tau(w_opts.tau_max_target);
    }
    if (cs->water_max_tau_override > 0.0) {
        w_opts.max_tau_max_target = cs->water_max_tau_override;
    }
    if (cs->water_max_z_max_m_override > 0.0) {
        w_opts.max_z_max_m = cs->water_max_z_max_m_override;
    }
    if (cs->water_depth_bottom_tol_override > 0.0) {
        w_opts.depth_bottom_tol = cs->water_depth_bottom_tol_override;
    }
    if (cs->water_layer_dtau_target_override > 0.0) {
        w_opts.layer_dtau_target = cs->water_layer_dtau_target_override;
    }
    if (cs->n_layers_water_override > 0) {
        w_opts.n_layers_water = cs->n_layers_water_override;
        /* Exact n-layer debug override: do not auto-raise via layer_dtau_target. */
        w_opts.layer_dtau_target = 0.0;
    }

    (void)opts;

    /* Debug-only water SOS max-iteration override. Production auto-selects a
     * cap and the SOS loop exits early on convergence. */
    {
        if (cs->water_max_orders_override > 0) {
            /* [2026-06-04] Guard: the fixed-bulk / CCRR branches above auto-raise
             * the in-water SOS cap to 160 / 140 so high-omega0 cases converge
             * (an omega0=0.9 fixed-bulk case needs ~90 orders). A CLI override
             * BELOW that auto-cap silently under-converges the in-water field
             * (lowers Lu / Rrs). Warn so a low override is used only for an
             * intentional single-scatter / diagnostic isolation, never for a
             * production / validation run. */
            if (cs->water_max_orders_override < w_opts.max_iterations) {
                fprintf(stderr,
                    "warning: --debug-water-max-orders %d is below the auto-selected "
                    "in-water SOS cap %d for this mode; high-omega0 in-water RT "
                    "may not converge (Lu/Rrs under-estimated). Intended only for "
                    "single-scatter / diagnostic isolation.\n",
                    cs->water_max_orders_override, w_opts.max_iterations);
            }
            w_opts.max_iterations = cs->water_max_orders_override;
        } else {
            const char *env_morders = ocrt_debug_env("OCRT_FORCE_WATER_MAX_ORDERS");
            if (env_morders && *env_morders) {
                int n = atoi(env_morders);
                if (n > 0) {
                    w_opts.max_iterations = n;
                    fprintf(stderr,
                        "warning: OCRT_FORCE_WATER_MAX_ORDERS is deprecated; "
                        "use --debug-water-max-orders %d instead\n", n);
                }
            }
        }
    }

    /* Resolve cs->F_sun with a safe default */
    double F_sun_TOA = (cs->F_sun > 0.0) ? cs->F_sun : M_PI;
    double T_C    = cs->T_water_C;       /* 0.0 valid (cold-water case); user-supplied via CLI */
    double S_gkg  = cs->S_water_gkg;
    double n_w    = (cs->n_water > 0.0) ? cs->n_water : 1.34;

    /* ====================================================================
     * Stage 2a — atmospheric direct beam attenuation
     *
     *   F_sun_at_BOA_direct = F_sun_TOA × exp(-τ_atm / μ_sun_air)
     *
     * Stage 2b step3a (this implementation): additionally couple atmospheric
     * diffuse skylight at BOA into in-water Ed via Snell + T_aw Mueller.
     *   Ed_water = Ed_direct_in_water + Ed_diffuse_in_water_from_atm
     *
     * NOTE (limitation): The diffuse coupling here only adds the *Ed
     * contribution*. The in-water Lu remains the *direct-beam-only* response.
     * Full Stage 2b (step3b — next turn) will pass the diffuse boundary field
     * into the in-water SOS source so that diffuse skylight also produces
     * upward Lu via scattering. This first-cut variant captures the dominant
     * Rrs reduction from added skylight Ed but underestimates Lu slightly.
     * ==================================================================== */
    double mu_sun_air = cos(cs->sza_deg * M_PI / 180.0);
    double tau_atm = 0.0;
    if (cs->rayleigh_on) {
        double tau_R = rt_rayleigh_tau_model_full(cs->wavelength_nm, cs->pressure_hpa, 45.0, 0.0, opts->rayleigh_model);
        tau_atm += tau_R;
    }
    /* Stage 2c — aerosol extinction added to atm optical depth.
     * aer may be NULL (no aerosol) or aer->tau_a == 0 (Rayleigh-only). */
    if (aer && aer->tau_a > 0.0) {
        tau_atm += aer->tau_a;
    }
    /* Gas absorption — direct-beam attenuation through atm column.
     * (Diffuse / multi-scattered Ed is treated via the layered τ_abs path
     * inside rt_atm_apply_gas_absorption, which is applied during atm SOS.) */
    double tau_abs_col = 0.0;
    if (opts->abs_state) {
        tau_abs_col = rt_absorption_tau_total(
            (const rt_absorption_t *)opts->abs_state, cs->wavelength_nm);
        tau_atm += tau_abs_col;
    }

    double F_sun_BOA = F_sun_TOA;
    if (mu_sun_air > 0.0 && tau_atm > 0.0) {
        F_sun_BOA = F_sun_TOA * exp(-tau_atm / mu_sun_air);
    }
    double ocrt_T_dir_beam = (mu_sun_air > 0 && tau_atm > 0)
        ? exp(-tau_atm / mu_sun_air) : 1.0;
    if (opts->pssa) {
        ocrt_T_dir_beam = ocrt_pssa_T_dir_dn_for_ocean(opts);
        /* negative = the IPSS scope guard refused this combination (see the
         * function body); abort instead of continuing with a bogus beam. */
        if (!(ocrt_T_dir_beam >= 0.0)) return -1;
        F_sun_BOA = F_sun_TOA * ocrt_T_dir_beam;
    }

    /* Stage 2b step3a/b: atm SOS + boundary coupling. Stage 2d (this turn):
     * also captures the atm path radiance I_TOA from the same atm SOS run
     * for use in TOA composition.
     *
     * atm_path_I/Q/U: atm path radiance at view direction (sun → atm scatter
     * → TOA, WITHOUT surface contribution since we run with surface=BLACK).
     * These are valid only if the atm SOS branch executes. */
    double Ed_diff_water_from_atm = 0.0;
    double Ed_diff_BOA_air = 0.0;        /* hemispheric diffuse Ed at atm BOA (air-side) */
    double Lu_diff_water_from_atm = 0.0;
    /* v1.10 B-0c.2b FINAL (Jae, 2026-07-10): diffuse-top is the DEFAULT.
     * Adjudicated METHODOLOGICALLY: OSOAA_SOS_CORE.F treats the interface as
     * per-mode matrix x diffuse-field (line 2485: XI1*TAW11+XQ1*TAW12+XU1*TAW13)
     * with only the solar direct as a beam (line 1747) - identical in form to
     * our diffuse top source.  Equivalent beams do not exist in OSOAA and are
     * retired; OCRT_EQUIV_BEAMS=1 restores them for archaeology only.
     * OCRT_DIFFUSE_TOP=m0 keeps the m0-only diagnostic; =0 also restores beams. */
    int ocrt_dt_mode = 1;
    { const char *e_dt = getenv("OCRT_DIFFUSE_TOP");
      const char *e_eb = getenv("OCRT_EQUIV_BEAMS");
      if (e_eb && e_eb[0] && e_eb[0] != '0') ocrt_dt_mode = 0;
      else if (e_dt && e_dt[0]) {
          if      (strcmp(e_dt, "m0") == 0) ocrt_dt_mode = 2;
          else if (strcmp(e_dt, "0")  == 0) ocrt_dt_mode = 0;
          else                              ocrt_dt_mode = 1;
      } }
    double *ocrt_dt_I = NULL, *ocrt_dt_Q = NULL, *ocrt_dt_U = NULL;
    /* The diffuse-top source stores ext_top_mu as a borrowed pointer and
     * consumes it later in the same solve. Keep the quadrature arrays at
     * function scope so their lifetime covers the complete water solve. */
    double mu_water_grid[256], w_water_grid[256];
    double Qu_diff_water_from_atm = 0.0;
    double Uu_diff_water_from_atm = 0.0;
    /* Cox-Munk T_wa (option A): above-water (0+) atm-sky-light contribution.
     * Each per-node beam call applies its own water->air transmission (flat or
     * Cox-Munk facet integration) inside rt_water_rt_sos_pure; we sum the
     * already-transmitted above-water Stokes here. By linearity of T_wa this
     * equals transmitting the summed in-water field (so flat/wind=0 is identical
     * to the previous single-Mueller-on-the-sum path). */
    double Lu_above_diff_from_atm = 0.0;
    double Qu_above_diff_from_atm = 0.0;
    double Uu_above_diff_from_atm = 0.0;
    double atm_path_I = 0.0, atm_path_Q = 0.0, atm_path_U = 0.0;
    double atm_T_diff_dn_hemi_cache = 0.0;  /* Stage 2d step2 */
    double atm_T_diff_dn_dir_cache = 0.0;
    int    atm_branch_executed = 0;
    /* C3-full: skylight(확산광)-유발 수중 0- 상향장의 m=0 mode (방위 등방).
     * BOA skylight m=0 은 방위 등방 → 등방 입사에 대한 수중 상향 응답도 m=0 only
     * (단일 빔 응답을 모든 방위에 대해 적분하면 m=0 만 잔존; R_iso=R^{m=0}).
     * 아래 equivalent-beam 루프에서 각 빔의 m=0 upwelling field 만 누적하고, 직달빔
     * (C3d)과 동일하게 역결합(water→atm) + 별도 atm SOS bottom-source pass 로 엄밀
     * 전파한다 → first-cut hemispheric T_atm_up 대체.  대기 RTE 선형성에 의해
     * TOA_wl_direct + TOA_wl_sky 가 전체 water-leaving 의 엄밀 TOA.  default ON;
     * OCRT_DEBUG=1 OCRT_C3FULL_OFF=1 로 first-cut 복귀(A/B). */
    const int c3full_on = 1; /* rigorous water-leaving: always on (physical) */
    double *sky_up_I_m0 = NULL, *sky_up_Q_m0 = NULL, *sky_up_U_m0 = NULL;
    double *sky_mu_pos  = NULL;
    double *sky_w_pos   = NULL;   /* v1.10 B-0b.2 */
    int     sky_n_mu_m0 = 0;
    /* C3-full 진단: field 기반 first-cut 비교용 — rigorous 와 동일한 m=0 above-water
     * 복사를 view 방향에서 보간한 값(transport 만 hemispheric T_atm_up 으로 대체).
     * rho_I 계산엔 미사용(순수 진단; kr0only 검증 런과도 호환). */
    double  sky_above_m0_at_view = 0.0;
    /* Stage 2b step3a/b + Stage 2c: atm SOS branch executes when atmosphere
     * is present — Rayleigh OR aerosol (OR both). */
    int atm_present = (cs->rayleigh_on || (aer && aer->tau_a > 0.0)) && tau_atm > 0.0;
    if (atm_present) {
        rt_case_t cs_atm = *cs;
        /* B.4 Stage 2 (bb=0 sanity check 응답): atm SOS surface BC를 FLAT
         * (또는 wind_speed > 0이면 COXMUNK)로 변경. 이로써 atm SOS가 *Fresnel
         * surface reflection contribution*을 *atm path radiance에 포함*시킨다.
         * 이전 RT_SURFACE_BLACK은 BOA downward radiance를 surface에서 *전부
         * 흡수*해서 *5-6% TOA under-prediction*을 만들었다.
         *
         * 이 변경으로:
         *   atm_path_I = atm scatter + Fresnel reflection of BOA downward
         *   I_TOA_total = atm_path_I + Lu_above × T_atm_up
         *                 ↑ flat Fresnel ocean과 일치        ↑ in-water signal */
        if (cs->wind_speed > 0.0) {
            cs_atm.surface = RT_SURFACE_BLACK_FRESNEL_OCEAN;
        } else {
            cs_atm.surface = RT_SURFACE_FLAT;
        }
        cs_atm.F_sun    = F_sun_TOA;            /* atmospheric uses TOA F_sun */
        cs_atm.n_water  = n_w;                  /* ensure surface uses our n */

        rt_atm_boa_export_t boa;
        /* FIX [diffuse-coupling]: the ocean-branch atmosphere solve runs with
         * view_as_node forced on (n_mu_atm = opts->n_mu + 1).  The BOA diffuse
         * export buffer must be allocated with that SAME n_mu, otherwise the
         * populate guard (boa.n_mu == n_mu_atm) fails silently and the entire
         * transmitted-skylight field stays zero (atmosphere then enters the
         * in-water RT only as direct-beam attenuation). */
        /* ---- v1.10 S7 atm-solve grid cache (gate) ------------------------
         * Batch grid rows re-solve an IDENTICAL atmosphere every row; the
         * heavy SOS mode fields are view-invariant.  When the #20 grid vza
         * list is present and OCRT_ATM_GRID_CACHE is set (batch driver
         * setenv; runtime-set like #22's water cache -> live getenv, NOT
         * snapshot-eligible), solve ONCE per key via the v1.01 LUT path
         * (SOS once, every (vza,raa) reconstructed from the same quadrature
         * solution; "bit-exact at grid points" per its header) with
         * view_as_node OFF (grid vza ARE the shared-N Gauss nodes; the
         * uangles dedup measured boa.n_mu 25->24 confirms coincidence), and
         * replay rows from the cached raw grids + boa export.  Pristine
         * per-row path is untouched when the gate is off (base rows,
         * standalone CLI, strict semantics). */
        int ocrt_s7_on = 0;
#ifdef OCRT_FAST_KERNELS
        /* FASTK-only (sprint rule): replay differs by ~1 ulp in near-zero U
         * (node-read vs ring-read reconstruction ordering) - inside the
         * 1e-10 fastk gate, but the STRICT build stays pristine bit-wise. */
        { const char *s7e = getenv("OCRT_ATM_GRID_CACHE");
          ocrt_s7_on = (s7e && s7e[0] && s7e[0] != '0') &&
                       cs->water_view_vza_list && cs->n_water_view_vza >= 2 &&
                       getenv("OCRT_S7_RAA_LIST") != NULL &&
                       /* S7-D: 뷰 μ가 노드 범위 밖(연직 부근)이면 행별 폴백 */
                       ocrt_s7_view_within_nodes(w_opts.n_mu_water,
                           cos(cs->vza_deg * M_PI / 180.0));
          /* 2026-07-14 M2 (Jae 지시): 기존 전제(격자 vza == 공유-N Gauss 노드,
           * 즉 grid n_mu == n_mu_water)를 제거한다.  M2 통일 격자는 균일 vza
           * (기본 2.5도)라 노드와 일치하지 않는다.  LUT 재구성은 노드 일치점
           * bit-exact, 비노드는 PCHIP(LUT-PCHIP, 실측 최대 0.005%)으로
           * 채워지므로 종전 "선형 보간 최대 292%(0-근방 상대오차 지표 발산)"
           * 우려는 해소되었다.  수중부는 #20 0-가중 노드로 여전히 정확 추출. */
        }
#endif
        int n_mu_atm_bc = w_opts.n_mu_water +
                          ((!ocrt_s7_on && opts->view_as_node) ? 1 : 0);
        if (!ocrt_s7_on && opts->view_as_node) {
            rt_uangles_t ocrt_ub; int ocrt_vd = -1;
            if (rt_uangles_init(&ocrt_ub, w_opts.n_mu_water) == 0 &&
                rt_uangles_add(&ocrt_ub, cos(cs->vza_deg * M_PI / 180.0),
                               &ocrt_vd) == 0)
                n_mu_atm_bc = ocrt_ub.n_total;   /* dedup-matched capacity */
        }
        ocrt_s7_cache_t *s7cp = &g_s7c;
        if (rt_atm_boa_export_alloc(&boa, opts->fourier_m_max, n_mu_atm_bc) == 0) {
            rt_result_t atm_res = {0};
            rt_options_t opts_atm = *opts;
            opts_atm.n_mu = w_opts.n_mu_water;   /* v1.10 B-0b.1: shared N
                * (OSOAA single-table rule): the atm SOS runs on the SAME
                * Gauss count as the water solve for this row.  --n-mu no
                * longer sizes the coupled-path atmosphere. */
            int rc_atm = -1;
            int s7_key_ok = ocrt_s7_on && s7cp->valid &&
                s7cp->lam == cs->wavelength_nm && s7cp->sza == cs->sza_deg &&
                s7cp->F == F_sun_TOA && s7cp->wind == cs->wind_speed &&
                s7cp->nw == n_w && s7cp->surface == (int)cs_atm.surface &&
                s7cp->sigt == cs->sigma_type && s7cp->qc == cs->q_convention &&
                s7cp->ray_on == cs->rayleigh_on && s7cp->aerp == (const void*)aer &&
                s7cp->aer_tau == (aer ? aer->tau_a : 0.0) &&
                s7cp->n_mu == w_opts.n_mu_water && s7cp->m_max == opts->fourier_m_max &&
                s7cp->n_vza == cs->n_water_view_vza &&
                s7cp->boa.n_mu == n_mu_atm_bc;
            OCRT_S6T("[S6T]  prep_end t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            int s7_iv = -1, s7_ir = -1;
            if (s7_key_ok) {
                for (int q = 0; q < s7cp->n_vza; ++q)
                    if (fabs(s7cp->vza[q] - cs->vza_deg) < 1e-9) { s7_iv = q; break; }
                for (int q = 0; q < s7cp->n_raa; ++q)
                    if (fabs(s7cp->raa[q] - cs->raa_deg) < 1e-9) { s7_ir = q; break; }
            }
            if (s7_iv >= 0 && s7_ir >= 0) {
                /* -------- S7 CACHE HIT: replay, no atm solve -------- */
                OCRT_S6T("[S6T]  atm1 s7-hit iv=%d ir=%d t=%.4f\n",
                         s7_iv, s7_ir, ocrt_s6_now()-ocrt_s6_t0);
                size_t pm = (size_t)(s7cp->boa.m_max + 1) * (size_t)s7cp->boa.n_mu;
                memcpy(boa.mu_quad_pos, s7cp->boa.mu_quad_pos,
                       (size_t)s7cp->boa.n_mu * sizeof(double));
                memcpy(boa.I_per_m, s7cp->boa.I_per_m, pm * sizeof(double));
                memcpy(boa.Q_per_m, s7cp->boa.Q_per_m, pm * sizeof(double));
                memcpy(boa.U_per_m, s7cp->boa.U_per_m, pm * sizeof(double));
                boa.tau_atm_total = s7cp->boa.tau_atm_total;
                boa.mu_sun_air    = s7cp->boa.mu_sun_air;
                size_t g = (size_t)s7_iv * (size_t)s7cp->n_raa + (size_t)s7_ir;
                atm_res.I_TOA = s7cp->raw[0][g];
                atm_res.Q_TOA = s7cp->raw[1][g];
                atm_res.U_TOA = s7cp->raw[2][g];
                atm_res.T_diff_dn_hemi = s7cp->T_diff_dn_hemi;
                rc_atm = 0;
            } else if (ocrt_s7_on) {
                /* -------- S7 MISS: one LUT-path solve, fill the cache ---- */
                OCRT_S6T("[S6T]  atm1 begin (s7-fill) t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
                /* (re)build the cache buffers for this key */
                int nv = cs->n_water_view_vza, nr = 0;
                { const char *rl = getenv("OCRT_S7_RAA_LIST");
                  char tmp[4096]; strncpy(tmp, rl, sizeof tmp - 1); tmp[sizeof tmp-1]=0;   /* M2: raa 2.5도 격자 수용 */
                  free(s7cp->raa); s7cp->raa = (double*)calloc(512, sizeof(double));
                  for (char *p = strtok(tmp, ","); p && nr < 512; p = strtok(NULL, ","))
                      s7cp->raa[nr++] = atof(p); }
                free(s7cp->vza); s7cp->vza = (double*)calloc((size_t)nv, sizeof(double));
                memcpy(s7cp->vza, cs->water_view_vza_list, (size_t)nv * sizeof(double));
                for (int q3 = 0; q3 < 3; ++q3) {
                    free(s7cp->raw[q3]); free(s7cp->rho[q3]); free(s7cp->pm[q3]);
                    s7cp->raw[q3] = (double*)calloc((size_t)nv * (size_t)nr, sizeof(double));
                    s7cp->rho[q3] = (double*)calloc((size_t)nv * (size_t)nr, sizeof(double));
                    s7cp->pm[q3] = (double*)calloc((size_t)(opts->fourier_m_max + 1) * (size_t)nv, sizeof(double));
                }
                free(s7cp->tdiff_dir);
                s7cp->tdiff_dir = (double*)calloc((size_t)nv * (size_t)nr, sizeof(double));
                rt_atm_boa_export_free(&s7cp->boa);
                rt_lut_grid_out_t s7lut; memset(&s7lut, 0, sizeof s7lut);
                s7lut.n_vza = nv; s7lut.vza_deg = s7cp->vza;
                s7lut.n_raa = nr; s7lut.raa_deg = s7cp->raa;
                s7lut.rho_I = s7cp->rho[0]; s7lut.rho_Q = s7cp->rho[1]; s7lut.rho_U = s7cp->rho[2];
                s7lut.raw_I = s7cp->raw[0]; s7lut.raw_Q = s7cp->raw[1]; s7lut.raw_U = s7cp->raw[2];
                s7lut.view_per_m_I = s7cp->pm[0]; s7lut.view_per_m_Q = s7cp->pm[1]; s7lut.view_per_m_U = s7cp->pm[2];
                s7lut.view_m_max = opts->fourier_m_max;
                s7lut.T_diff_dn_dir = s7cp->tdiff_dir;
                rt_options_t opts_s7 = opts_atm;
                opts_s7.view_as_node = 0;   /* grid vza ARE nodes; interp path
                                             * is bit-exact at coinciding nodes */
                int rc_fill = -1;
                if (rt_atm_boa_export_alloc(&s7cp->boa, opts->fourier_m_max,
                                            n_mu_atm_bc) == 0)
                    rc_fill = rt_solve_case_pol_impl(&cs_atm, &opts_s7, aer,
                                                     &atm_res, &s7lut, &s7cp->boa);
                OCRT_S6T("[S6T]  atm1 end (s7-fill) rc=%d t=%.4f\n", rc_fill, ocrt_s6_now()-ocrt_s6_t0);
                if (rc_fill == 0) {
                    s7cp->valid = 1; s7cp->n_vza = nv; s7cp->n_raa = nr;
                    s7cp->n_mu = w_opts.n_mu_water; s7cp->m_max = opts->fourier_m_max;
                    s7cp->lam = cs->wavelength_nm; s7cp->sza = cs->sza_deg;
                    s7cp->F = F_sun_TOA; s7cp->wind = cs->wind_speed; s7cp->nw = n_w;
                    s7cp->surface = (int)cs_atm.surface; s7cp->sigt = cs->sigma_type;
                    s7cp->qc = cs->q_convention; s7cp->ray_on = cs->rayleigh_on;
                    s7cp->aerp = (const void*)aer; s7cp->aer_tau = aer ? aer->tau_a : 0.0;
                    s7cp->T_diff_dn_hemi = s7lut.T_diff_dn_hemi;
                    /* replay THIS row from the grids (same source as later rows) */
                    int iv2 = -1, ir2 = -1;
                    for (int q = 0; q < nv; ++q)
                        if (fabs(s7cp->vza[q] - cs->vza_deg) < 1e-9) { iv2 = q; break; }
                    for (int q = 0; q < nr; ++q)
                        if (fabs(s7cp->raa[q] - cs->raa_deg) < 1e-9) { ir2 = q; break; }
                    if (iv2 >= 0 && ir2 >= 0) {
                        size_t g2 = (size_t)iv2 * (size_t)nr + (size_t)ir2;
                        atm_res.I_TOA = s7cp->raw[0][g2];
                        atm_res.Q_TOA = s7cp->raw[1][g2];
                        atm_res.U_TOA = s7cp->raw[2][g2];
                        atm_res.T_diff_dn_hemi = s7cp->T_diff_dn_hemi;
                    }
                    memcpy(boa.mu_quad_pos, s7cp->boa.mu_quad_pos,
                           (size_t)s7cp->boa.n_mu * sizeof(double));
                    size_t pm2 = (size_t)(s7cp->boa.m_max + 1) * (size_t)s7cp->boa.n_mu;
                    memcpy(boa.I_per_m, s7cp->boa.I_per_m, pm2 * sizeof(double));
                    memcpy(boa.Q_per_m, s7cp->boa.Q_per_m, pm2 * sizeof(double));
                    memcpy(boa.U_per_m, s7cp->boa.U_per_m, pm2 * sizeof(double));
                    boa.tau_atm_total = s7cp->boa.tau_atm_total;
                    boa.mu_sun_air    = s7cp->boa.mu_sun_air;
                    rc_atm = 0;
                } else { s7cp->valid = 0; rc_atm = rc_fill; }
            } else {
                /* -------- exact near-nadir fallback + S17a RAA fold --------
                 * S7 is off here because the requested view is outside the
                 * Gauss-node interpolation range (or because the grid cache is
                 * disabled).  In an actual full-grid run, keep the exact
                 * view-as-node solution and fold only RAA. */
                ocrt_s17a_cache_t *s17ap = &g_s17a;
                const char *s17a_grid = getenv("OCRT_ATM_GRID_CACHE");
                const int s17a_on = !getenv("OCRT_S17A_OFF") &&
                    s17a_grid && s17a_grid[0] && s17a_grid[0] != '0' &&
                    getenv("OCRT_S7_RAA_LIST") != NULL &&
                    opts_atm.fourier_m_max >= 0 && opts_atm.fourier_m_max < 33;
                const unsigned long long k17 = s17a_on ?
                    ocrt_s17a_key(&cs_atm, &opts_atm, aer) : 0ULL;
                if (s17a_on && s17ap->valid && s17ap->key == k17 &&
                    s17ap->m_max == opts_atm.fourier_m_max &&
                    ocrt_boa_copy(&boa, &s17ap->boa) == 0) {
                    const double phi = rt_raa_to_atm_fourier_phi(cs->raa_deg);
                    atm_res.I_TOA = rt_solver_reconstruct_phi(s17ap->pmI, s17ap->m_max, phi);
                    atm_res.Q_TOA = rt_solver_reconstruct_phi(s17ap->pmQ, s17ap->m_max, phi);
                    atm_res.U_TOA = rt_solver_reconstruct_phi_sin(s17ap->pmU, s17ap->m_max, phi);
                    atm_res.T_diff_dn_hemi = s17ap->T_diff_dn_hemi;
                    atm_res.T_diff_dn_dir = ocrt_boa_T_diff_dn_dir(
                        &s17ap->boa, cos(cs->vza_deg * M_PI / 180.0),
                        cs->raa_deg, cos(cs->sza_deg * M_PI / 180.0));
                    rc_atm = 0;
                    OCRT_S6T("[S6T]  atm1 s17a-hit t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
                } else {
                    if (s17a_on) {
                        opts_atm.toa_view_per_m_I = s17ap->pmI;
                        opts_atm.toa_view_per_m_Q = s17ap->pmQ;
                        opts_atm.toa_view_per_m_U = s17ap->pmU;
                    }
                    OCRT_S6T("[S6T]  atm1 begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
                    rc_atm = rt_solve_case_pol_for_ocean(&cs_atm, &opts_atm, aer,
                                                         &atm_res, &boa);
                    OCRT_S6T("[S6T]  atm1 end rc=%d t=%.4f\n", rc_atm, ocrt_s6_now()-ocrt_s6_t0);
                    if (rc_atm == 0 && s17a_on &&
                        ocrt_boa_copy(&s17ap->boa, &boa) == 0) {
                        s17ap->valid = 1;
                        s17ap->key = k17;
                        s17ap->m_max = opts_atm.fourier_m_max;
                        s17ap->T_diff_dn_hemi = atm_res.T_diff_dn_hemi;
                    } else if (s17a_on) {
                        s17ap->valid = 0;
                    }
                }
            }
            if (rc_atm == 0) {
                if (getenv("OCRT_S6_TRACE")) {
                    double ocrt_bmax = 0.0;
                    for (int jj = 0; jj < n_mu_atm_bc; ++jj) {
                        double v = boa.I_per_m[0 * n_mu_atm_bc + jj];
                        if (v > ocrt_bmax) ocrt_bmax = v;
                    }
                    fprintf(stderr, "[S6B] boa.n_mu=%d expect=%d Imax_m0=%.3e\n",
                            boa.n_mu, n_mu_atm_bc, ocrt_bmax);
                }
                /* Stage 2d: capture atm path radiance for TOA composition */
                atm_path_I = atm_res.I_TOA;
                atm_path_Q = atm_res.Q_TOA;
                atm_path_U = atm_res.U_TOA;
                atm_T_diff_dn_hemi_cache = atm_res.T_diff_dn_hemi;  /* step2d-2 */
                atm_T_diff_dn_dir_cache = atm_res.T_diff_dn_dir;
                atm_branch_executed = 1;
                const int n_mu_water = w_opts.n_mu_water;
                /* [2026-05-31] BUG-WRT-003 fix.  The transmitted diffuse skylight
                 * populates only the Snell refraction cone mu_w in [mu_crit, 1];
                 * below mu_crit the in-water downwelling boundary source is
                 * identically zero (TIR).  Representing this kinked integrand on a
                 * [0,1] Gauss-Legendre grid mishandles the cone-edge transition and
                 * over-estimates the transmitted diffuse irradiance by ~7.5% at the
                 * default n_mu_water=24 (effective skylight transmittance ~0.94 vs
                 * the correct ~0.86; converges only slowly/oscillatorily with
                 * n_mu_water: 24->0.94, 48->0.92, 96->0.85).  Fix: build the
                 * quadrature that represents the transmitted skylight (used below
                 * for BOTH Ed(0-) and the equivalent-beam Lu superposition) on the
                 * cone [mu_crit, 1] by a linear map of the [0,1] Gauss rule.  The
                 * integrand is smooth there, so M~16 already gives 0.864.  The
                 * in-water SOS keeps its own independent [0,1] grid (generated
                 * inside rt_water_rt_sos_pure), so multiple-scattering is unaffected. */
                {
                    const double mu_crit_cone = sqrt(1.0 - 1.0/(n_w*n_w));
                    double xg[256], wg[256];
                    rt_quadrature_gauss_legendre_pos(n_mu_water, xg, wg);
                    if (cs->wind_speed > 0.0) {
                        /* [FIX-SKY-EDLU 2026-06-28] rough sea surface transmits skylight
                         * BELOW the flat critical angle (sub-cone leak, mu_w<mu_crit),
                         * which the air-side Ed(0-) integral includes.  Sample the FULL
                         * hemisphere [0,1] so the air->water Cox-Munk BTDF coupling
                         * populates the sub-cone in-water field too, and the Lu
                         * equivalent-beam superposition (below) can drive sub-cone beams
                         * -> numerator and denominator carry the same full transmitted
                         * field.  The rough-surface field is smooth across mu_crit (no
                         * flat-Snell cone-edge kink), so a plain [0,1] GL rule is fine. */
                        for (int k = 0; k < n_mu_water; ++k) {
                            mu_water_grid[k] = xg[k];
                            w_water_grid[k]  = wg[k];
                        }
                    } else {
                        /* flat (wind=0): transmitted field is identically zero below
                         * mu_crit (TIR), a kinked integrand on [0,1] -> cone-map onto
                         * [mu_crit,1] (BUG-WRT-003 fix). */
                        for (int k = 0; k < n_mu_water; ++k) {
                            mu_water_grid[k] = mu_crit_cone + (1.0 - mu_crit_cone) * xg[k];
                            w_water_grid[k]  = (1.0 - mu_crit_cone) * wg[k];
                        }
                    }
                }

                rt_aw_coupled_field_t coupled;
                if (rt_aw_coupled_field_alloc(&coupled, opts->fourier_m_max, n_mu_water) == 0) {
                    /* v1.10 B-0c.1: rebuild the atm ring table (identical
                     * construction to the boa ring) and hand its weights to
                     * the forward coupling for the full T_aw matrix path;
                     * size-mismatch -> NULL -> legacy hybrid. */
                    rt_uangles_t ocrt_ufw; int ocrt_vfw = -1;
                    const double *w_atm_fw = NULL;
                    if (rt_uangles_init(&ocrt_ufw, w_opts.n_mu_water) == 0) {
                        int okf = 1;
                        /* 2026-07-14 S7-boa FIX (Jae 승인): 링의 관측노드 추가
                         * 조건을 boa 구성 조건(n_mu_atm_bc 산정식)과 동일하게
                         * 맞춘다.  S7 재생 행의 boa는 노드 순수(뷰 없음)인데
                         * 종전 코드는 무조건 뷰를 추가해 크기 25 != boa 24 가
                         * 되었고, 아래 크기 검사 실패로 w_atm_fw=NULL -> 결합이
                         * 레거시 하이브리드(평면 Fresnel + 선형 보간)로 조용히
                         * 강등되었다.  실측: 격자 rrs0-가 단일 기하 대비
                         * 0.11~0.30% 이탈(n_mu_water 반비례 감소).  뷰 노드는
                         * 가중 0이라 결합 합(wj==0 skip)에 어차피 불참하므로,
                         * 조건 일치 후 두 모드의 결합은 수학적으로 동일하다. */
                        if (!ocrt_s7_on && opts->view_as_node)
                            okf = (rt_uangles_add(&ocrt_ufw,
                                      cos(cs->vza_deg * M_PI / 180.0),
                                      &ocrt_vfw) == 0);
                        if (okf && ocrt_ufw.n_total == boa.n_mu)
                            w_atm_fw = ocrt_ufw.w;
                    }
                    if (rt_air_water_couple_atm_to_water(&boa, n_w,
                                                          cs->q_convention,
                                                          mu_water_grid, n_mu_water,
                                                          w_atm_fw,
                                                          cs->wind_speed,
                                                          cs->sigma_type,
                                                          &coupled) == 0) {
                        /* v1.054 NORM FIX: coupled.{I,Q,U}_inwater_per_m inherit the
                         * F_sun=pi-normalized convention of boa.I_per_m (the coupling
                         * applies only the T_aw Mueller transform). Convert to SI by
                         * × F_sun_TOA/pi so the in-water diffuse Ed and the step3c
                         * equivalent-beam flux F_sun_eq (→ Lu_diff) are consistent with
                         * the SI direct beam (w_res) and with Ed_diff_BOA_air (already
                         * SI). Prior bug: in-water diffuse was pi× too large → unphysical
                         * Ed(0-)>Ed(0+) and ~+43% inflated Rrs(0+). */
                        {
                            if (ocrt_dt_mode) {
                                /* B-0c.2b: snapshot the pi-normalized coupled field
                                 * BEFORE SI conversion and hand it to the water solver
                                 * as the diffuse top source (beams then bypass). */
                                size_t ndt = (size_t)(opts->fourier_m_max + 1)
                                           * (size_t)n_mu_water;
                                ocrt_dt_I = (double*)malloc(ndt * sizeof(double));
                                ocrt_dt_Q = (double*)malloc(ndt * sizeof(double));
                                ocrt_dt_U = (double*)malloc(ndt * sizeof(double)); /* v1.10 A-FIX 2026-07-12: wire the U-incident (was "next increment") */
                                if (ocrt_dt_I && ocrt_dt_Q && ocrt_dt_U) {
                                    memcpy(ocrt_dt_I, coupled.I_inwater_per_m, ndt * sizeof(double));
                                    memcpy(ocrt_dt_Q, coupled.Q_inwater_per_m, ndt * sizeof(double));
                                    memcpy(ocrt_dt_U, coupled.U_inwater_per_m, ndt * sizeof(double));
                                    /* DTPSIGN (2026-08-01, parity v8, Jae approved): the boa/coupled
                                     * per-m Fourier fields are in the atmosphere azimuth base
                                     * (phi + pi); the water SOS source basis is cos/sin(m*phi).
                                     * Convert odd modes by (-1)^m at this hand-off, so the
                                     * injector's source/operator closure identity stays intact
                                     * (test_stage2_diffuse_top_uclosure).  Verified against OSOAA
                                     * (pseudo-pure 443/555, SZA 20/40/60, wind 3/7): m=1 gap
                                     * -2.40% -> +0.11% of c0; even modes bit-identical; full-azimuth
                                     * Rrs(0+) reaches I 0.106% / Q 0.120% / U 0.056% / DoLP 0.110pp.
                                     * Runs once per case outside hot loops: zero runtime cost. */
                                    for (int mm_ = 1; mm_ <= opts->fourier_m_max; mm_ += 2)
                                        for (int j_ = 0; j_ < n_mu_water; ++j_) {
                                            size_t z_ = (size_t)mm_ * (size_t)n_mu_water + (size_t)j_;
                                            ocrt_dt_I[z_] = -ocrt_dt_I[z_];
                                            ocrt_dt_Q[z_] = -ocrt_dt_Q[z_];
                                            ocrt_dt_U[z_] = -ocrt_dt_U[z_];
                                        }
                                    { const char *oc_ = getenv("OCRT_DTP_ONLYC");  /* v1.10 probe: keep one incident column */
                                      if (oc_) { int c_ = atoi(oc_);
                                        for (int mm_ = 0; mm_ <= opts->fourier_m_max; ++mm_)
                                          for (int j_ = 0; j_ < n_mu_water; ++j_)
                                            if (j_ != c_) { ocrt_dt_I[(size_t)mm_*n_mu_water+j_]=0.0;
                                                            ocrt_dt_Q[(size_t)mm_*n_mu_water+j_]=0.0;
                                                            ocrt_dt_U[(size_t)mm_*n_mu_water+j_]=0.0; } } }
                                    if (ocrt_dt_mode == 2) {
                                        for (size_t z = (size_t)n_mu_water; z < ndt; ++z) {
                                            ocrt_dt_I[z] = 0.0; ocrt_dt_Q[z] = 0.0;
                                        }
                                    }
                                    w_opts.ext_top_I     = ocrt_dt_I;
                                    w_opts.ext_top_Q     = ocrt_dt_Q;
                                    w_opts.ext_top_U     = ocrt_dt_U;   /* v1.10 A-FIX 2026-07-12 */
                                    w_opts.ext_top_mu    = mu_water_grid;
                                    w_opts.ext_top_n     = n_mu_water;
                                    w_opts.ext_top_m_max = opts->fourier_m_max;
                                }
                            }
                            double si_conv = F_sun_TOA / M_PI;
                            int n_el = (opts->fourier_m_max + 1) * n_mu_water;
                            for (int e = 0; e < n_el; ++e) {
                                coupled.I_inwater_per_m[e] *= si_conv;
                                coupled.Q_inwater_per_m[e] *= si_conv;
                                coupled.U_inwater_per_m[e] *= si_conv;
                            }
                        }
                        /* [2026-05-31] FIX-B: diffuse skylight transmitted
                         * irradiance Ed(0-) via the EXACT air-side rough-surface
                         * (Cox-Munk) transmittance integral.  For each air
                         * direction surface_T_aw_coxmunk_direct gives the slope-
                         * integrated facet transmittance (= 1-R(θ) flat at
                         * wind=0; OSOAA-validated).  By the identity
                         *   ∫_0^1 μ_w K^0(μ_w,μ_air) dμ_w = T_aw_direct(μ_air)
                         * this equals the hemispheric integral of the transmitted
                         * in-water field, but is computed on the smooth air grid
                         * (no Snell-cone-edge kink — supersedes the FIX-A cone-quad)
                         * and is exact for all wind, including the sub-cone leak.
                         * Same air grid and F_sun=π→SI normalization as the air-
                         * side incident Ed_diff_BOA_air below. */
                        /* v1.10 B-0b.1 FIX: this integral must ride the SAME
                         * unified ring as the boa export (shared N + view).
                         * Stale GL(opts->n_mu) weights paired with 24-node boa
                         * radiances read the six GRAZING-most sky directions
                         * with full-hemisphere weights, inflating Ed by ~9% -
                         * caught by OSOAA flux arbitration (Ed(0+): OSOAA
                         * 2.676, legacy-6 2.639, broken-mix 2.884). */
                        rt_uangles_t ocrt_ued; int ocrt_ved = -1;
                        int ocrt_ued_ok = (rt_uangles_init(&ocrt_ued,
                                               w_opts.n_mu_water) == 0);
                        /* 2026-07-14 S7-boa FIX #2 (Jae 승인): 5886 지점과 동일
                         * 결함.  S7 재생 행의 boa 는 노드 순수(24)인데 링이
                         * 무조건 뷰를 추가해 25 가 되면, 아래 Ed 적분과 S9
                         * 블록이 boa.I_per_m 을 n_mu_atm=25 보폭으로 읽어
                         * m=0 행 경계를 넘는 오배열 접근이 된다(실측: 격자
                         * rrs0- 잔여 0.104% 이탈의 원인).  boa 구성 조건
                         * (n_mu_atm_bc 산정식)과 동일하게 게이트한다. */
                        if (ocrt_ued_ok && !ocrt_s7_on && opts->view_as_node)
                            ocrt_ued_ok = (rt_uangles_add(&ocrt_ued,
                                cos(cs->vza_deg * M_PI / 180.0),
                                &ocrt_ved) == 0);
                        const int n_mu_atm = ocrt_ued_ok ? ocrt_ued.n_total : 0;
                        const double *mu_atm = ocrt_ued.mu;
                        const double *w_atm  = ocrt_ued.w;
                        /* [FIX-SKY-EDLU 2026-06-28] transmitted diffuse Ed(0-) via the
                         * EXACT air-side rough-surface (Cox-Munk) transmittance integral
                         * surface_T_aw_coxmunk_direct.  This is the physically correct
                         * transmitted irradiance: it includes the rough-surface sub-cone
                         * leak (skylight transmitted into in-water mu_w<mu_crit through
                         * tilted facets), which grows with wind and with grazing sun.
                         * The Lu(0-) numerator is made consistent by ALSO driving
                         * equivalent beams for the sub-cone in-water directions (see the
                         * full-hemisphere grid for wind>0 and the sub-cone branch in the
                         * Lu superposition loop below), so numerator and denominator are
                         * the water response to one identical full transmitted field. */
                        OCRT_S6T("[S6T]  edlu begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
#ifdef OCRT_FAST_KERNELS
                        /* v1.10 S9 (2026-07-10): row-invariant Td hoist.
                         * surface_T_aw_coxmunk_direct(mu_j, n_w, wind, sigma)
                         * depends ONLY on the fixed air-node set and the
                         * (n_w, wind, sigma) triple - identical for every row
                         * of a batch grid (S6T probe: this loop = 5.9 ms/row,
                         * the wind>0 path is a 101x101 slope-Fresnel integral
                         * per node).  Cache the Td[] table thread-locally,
                         * keyed by count + FNV of node mus + the triple;
                         * same-input double replay, bit-identical by design.
                         * FASTK-only by S2/S7 convention; strict pristine. */
                        {
                            static _Thread_local double s9_td[64];
                            static _Thread_local double s9_nw, s9_ws;
                            static _Thread_local unsigned long long s9_h;
                            static _Thread_local int s9_n = -1, s9_st;
                            unsigned long long h_ = 1469598103934665603ULL;
                            for (int j = 0; j < n_mu_atm; ++j) {
                                union { double d; unsigned long long u; } cv_;
                                cv_.d = mu_atm[j];
                                h_ ^= cv_.u; h_ *= 1099511628211ULL;
                            }
                            if (!(s9_n == n_mu_atm && s9_h == h_ &&
                                  s9_nw == n_w && s9_ws == cs->wind_speed &&
                                  s9_st == cs->sigma_type) ) {
                                for (int j = 0; j < n_mu_atm && j < 64; ++j)
                                    s9_td[j] = surface_T_aw_coxmunk_direct(
                                                   mu_atm[j], n_w,
                                                   cs->wind_speed,
                                                   cs->sigma_type);
                                s9_n = n_mu_atm; s9_h = h_; s9_nw = n_w;
                                s9_ws = cs->wind_speed; s9_st = cs->sigma_type;
                            }
                            for (int j = 0; j < n_mu_atm; ++j) {
                                double Td = (j < 64) ? s9_td[j]
                                          : surface_T_aw_coxmunk_direct(
                                                mu_atm[j], n_w, cs->wind_speed,
                                                cs->sigma_type);
                                double muL = w_atm[j] * mu_atm[j] *
                                             boa.I_per_m[0 * n_mu_atm + j];
                                Ed_diff_water_from_atm += muL * Td;
                                Ed_diff_BOA_air        += muL;
                            }
                        }
#else
                        for (int j = 0; j < n_mu_atm; ++j) {
                            double Td = surface_T_aw_coxmunk_direct(
                                            mu_atm[j], n_w, cs->wind_speed,
                                            cs->sigma_type);
                            double muL = w_atm[j] * mu_atm[j] *
                                         boa.I_per_m[0 * n_mu_atm + j];
                            Ed_diff_water_from_atm += muL * Td;   /* transmitted */
                            Ed_diff_BOA_air        += muL;        /* incident */
                        }
#endif
                        OCRT_S6T("[S6T]  edlu end t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
                        Ed_diff_water_from_atm *= 2.0 * M_PI * (F_sun_TOA / M_PI);
                        if (getenv("OCRT_DTP_FLUX") && ocrt_dt_I) {   /* v1.10 probe 2026-07-12 */
                            double fx = 0.0;
                            for (int c_ = 0; c_ < n_mu_water; ++c_)
                                fx += 2.0 * M_PI * mu_water_grid[c_] * w_water_grid[c_] * ocrt_dt_I[c_];
                            fx *= (F_sun_TOA / M_PI);   /* snapshot is pi-normalized */
                            if (getenv("OCRT_DTP_FLUX")[0]=='2')
                                for (int c_ = 0; c_ < n_mu_water; ++c_)
                                    fprintf(stderr, "[DTPdist] c=%d mu_w=%.4f w=%.4f I0=%.4e\n",
                                            c_, mu_water_grid[c_], w_water_grid[c_], ocrt_dt_I[c_]);
                            fprintf(stderr, "[DTPflux] inj_m0_flux=%.6e  Ed_diff_analytic=%.6e  ratio=%.4f\n",
                                    fx, Ed_diff_water_from_atm,
                                    (Ed_diff_water_from_atm != 0.0) ? fx / Ed_diff_water_from_atm : -1.0);
                        }
                        Ed_diff_BOA_air        *= 2.0 * M_PI * (F_sun_TOA / M_PI);

                        if (getenv("OCRT_DUMP_SKY")) {
                            for (int jj = 0; jj < n_mu_atm; ++jj)
                                fprintf(stderr, "SKYAIR %.6f %.10e\n",
                                        mu_atm[jj], boa.I_per_m[0 * n_mu_atm + jj]);
                            for (int jj = 0; jj < n_mu_water; ++jj)
                                fprintf(stderr, "SKYWAT %.6f %.10e %.10e\n",
                                        mu_water_grid[jj],
                                        coupled.I_inwater_per_m[0 * n_mu_water + jj],
                                        w_water_grid[jj]);
                            fprintf(stderr, "SKYED %.10e %.10e\n",
                                    Ed_diff_water_from_atm, Ed_diff_BOA_air);
                        }

                        /* step3b: Lu contribution via multi-call superposition.
                         * step3c (this turn): correct the implementation by
                         * averaging each j'-call over N=4 raa-shifted samples.
                         *
                         * Mathematical correction:
                         *   True m=0 inbound boundary condition gives an
                         *   azimuthally-isotropic incoming radiance. The in-water
                         *   response is then
                         *     Lu_true = Σ_j' (1/2π) ∫_0^{2π} Lu_single(raa) d(raa)
                         *             = Σ_j' Lu^{m=0}_single
                         *   step3b used a single raa-fixed call, which conflates
                         *   m=0 and m≥1 single-call modes (raa-dependent error
                         *   ~1.7% measured). raa-quadrature with N=4 cancels
                         *   m=1, m=2, m=3 modes exactly (atm m_max=2 → sufficient).
                         *
                         * Limitation remaining: m≥1 atm BOA amplitude contributions
                         * (Rayleigh polarization angular pattern) still treated as
                         * zero in this step3c-1 — full step3c-2 deferred.
                         */
                        const double mu_crit = sqrt(1.0 - 1.0/(n_w*n_w));
                        /* C3-full rho_I(ON) 은 skylight m=0 field(kr=0, 방위 무관)만 쓴다.
                         * OCRT_DEBUG=1 OCRT_C3FULL_KR0ONLY=1 이면 N_raa=1 로 줄여 skylight
                         * equivalent-beam 루프를 4× 단축 → rho_I(ON) 정확(불변), view radiance
                         * (0+ 출력·first-cut 진단)만 단일-azimuth 근사가 된다.  검증 런타임
                         * 제약(단일코어·tool 한도)용; 정상 런은 N_raa=4(불변). */
                        const int    N_raa_quad =
                            (ocrt_debug_env("OCRT_C3FULL_KR0ONLY") != NULL) ? 1 : 4;
                        OCRT_S6T("[S6T]  beams begin t=%.2f n_mu_water=%d mu0=%.4f mu_last=%.4f\n", ocrt_s6_now()-ocrt_s6_t0, n_mu_water, mu_water_grid[0], mu_water_grid[n_mu_water-1]);
                        for (int j = 0; j < n_mu_water; ++j) {
                            if (ocrt_dt_mode) break;   /* B-0c.2b: diffuse-top supersedes equivalent beams */
                            double mu_w = mu_water_grid[j];
                            double w_j  = w_water_grid[j];

                            double I_in = coupled.I_inwater_per_m[0 * n_mu_water + j];
                            if (I_in <= 0.0) continue;

                            /* [FIX-SKY-EDLU 2026-06-28] equivalent-beam direction.
                             * cone (mu_w>mu_crit): a real air incidence angle exists,
                             *   drive the in-water beam via its air equivalent (Snell).
                             * sub-cone (mu_w<=mu_crit, wind>0 only): rough-surface leak;
                             *   no air angle -> drive the in-water beam DIRECTLY at mu_w
                             *   (mu_sun_water_override).  Both carry the same physical
                             *   transmitted-skylight field that the air-side Ed integral
                             *   sums, so rrs(0-) numerator/denominator are consistent. */
                            double sza_air_eq_d;
                            double mu_sun_w_override = 0.0;  /* 0 => air-refraction path */
                            if (mu_w > mu_crit) {
                                double sin2_w   = 1.0 - mu_w*mu_w;
                                double sin2_air = n_w*n_w * sin2_w;
                                if (sin2_air >= 1.0) continue;
                                double mu_air_eq    = sqrt(1.0 - sin2_air);
                                sza_air_eq_d = acos(mu_air_eq) * 180.0 / M_PI;
                            } else {
                                if (cs->wind_speed <= 0.0) continue; /* flat: no sub-cone */
                                sza_air_eq_d     = 0.0;   /* placeholder; mu_sun_air=1>0 */
                                mu_sun_w_override = mu_w;  /* in-water beam cosine direct */
                            }
                            /* FIX-D (2026-05-31 KST 2단계): step3b diffuse 상향 정규화 — 과잉 /mu_w 제거.
                               diffuse 방향 j의 수중 horizontal irradiance = I_in*mu_w*(2pi*w_j);
                               이를 대체하는 등가 collimated 빔의 in-water normal flux
                               = horizontal/mu_w = I_in*(2pi*w_j) (mu_w 상쇄). 이전의 추가 /mu_w 는
                               각 diffuse 방향 상향 기여를 1/mu_w배 과대화(cone edge=grazing skylight에서 최대)
                               시켜 grazing rrs(0-) 반등을 유발했다. 직달 항 F_sun_water(.../mu_sun_water_corr)
                               와 동일 normal-flux convention 으로 통일하는 물리 기반 수정(임의 튜닝 아님). */
                            { const char *oj_ = getenv("OCRT_BEAM_ONLYJ");   /* v1.10 probe: keep one beam */
                              if (oj_ && j != atoi(oj_)) continue; }
                            double F_sun_eq    = I_in * 2.0 * M_PI * w_j;
#ifdef OCRT_FAST_KERNELS
                            /* v1.09-opt S6a: numerically-null skylight beams.
                             * Sub-cone Cox-Munk leak directions can carry
                             * denormal-level flux (measured 1e-45..1e-52);
                             * their contribution is far below one ulp of any
                             * accumulated output, yet each costs a full
                             * multi-mode multi-order water solve.  Skip. */
                            if (!isfinite(F_sun_eq) || F_sun_eq < 1e-30) { ocrt_s6_skipped++; continue; }
#endif


                            double Lu_j_avg = 0.0, Q_j_avg = 0.0, U_j_avg = 0.0;
                            double La_j_avg = 0.0, Qa_j_avg = 0.0, Ua_j_avg = 0.0;
                            for (int kr = 0; kr < N_raa_quad; ++kr) {
                                double raa_shift_deg = (360.0 * kr) / N_raa_quad;
                                double raa_call = cs->raa_deg + raa_shift_deg;
                                /* keep raa_call in [0, 360) — solver handles cyclic */
                                while (raa_call >= 360.0) raa_call -= 360.0;
                                while (raa_call <    0.0) raa_call += 360.0;

                                rt_water_rt_options_t w_opts_local = w_opts;
                                w_opts_local.view_vza_deg_list = NULL;  /* #20: sky beams need no view nodes */
                                w_opts_local.n_view_vza = 0;
                                /* [FIX-SKY-EDLU] sub-cone beams drive the in-water solar
                                 * cosine directly (no air-refraction); 0 for cone beams. */
                                w_opts_local.mu_sun_water_override = mu_sun_w_override;
                                /* sky-light coupling speedup (CLI --n-mu-water-sky): solve these
                                 * small (~2%) skylight beams on a COARSE in-water grid.  The
                                 * direct solar beam (w_res, above) keeps the full n_mu_water, so
                                 * the dominant signal stays accurate; the downstream water->atm
                                 * coupling (rt_air_water_couple_water_to_atm) interpolates from
                                 * sky_mu_pos, so a coarse grid here is self-consistent and safe.
                                 * Cuts the per-beam SOS cost ~ (n_mu_water / n_mu_water_sky)^2. */
                                if (cs->n_mu_water_sky > 0 &&
                                    cs->n_mu_water_sky < w_opts_local.n_mu_water)
                                    w_opts_local.n_mu_water = cs->n_mu_water_sky;
                                rt_water_rt_result_t r_k = {0};
                                /* C3-full: on the first azimuth sample only, capture the
                                 * in-water 0- upwelling m=0 field (azimuth-independent).
                                 * Allocating these buffers makes rt_water_rt_sos_pure
                                 * expose the per-mode 0- field (no-op otherwise). */
                                const int sky_capture = (c3full_on && kr == 0);
                                if (sky_capture) {
                                    int cM = w_opts_local.m_max_water + 1;
                                    /* FIX(2026-06-13): same capacity bug as the
                                     * C3d site — solver appends up to 4 nodes. */
                                    int cN = w_opts_local.n_mu_water + 4;   /* view list off in sky-beam locals */
                                    r_k.I_up_per_m   = (double*)calloc((size_t)cM*cN, sizeof(double));
                                    r_k.Q_up_per_m   = (double*)calloc((size_t)cM*cN, sizeof(double));
                                    r_k.U_up_per_m   = (double*)calloc((size_t)cM*cN, sizeof(double));
                                    r_k.mu_water_pos = (double*)calloc((size_t)cN, sizeof(double));
                                    r_k.w_water_pos  = (double*)calloc((size_t)cN, sizeof(double));
                                }
                                ocrt_s6_solved++;
                                OCRT_S6T("[S6T]   beam j=%d kr=%d F=%.3e mu_w=%.4f t=%.4f\n", j, kr, F_sun_eq, mu_w, ocrt_s6_now()-ocrt_s6_t0);
#ifdef OCRT_FAST_KERNELS
                                /* v1.09-opt S7a: solve sky beams at UNIT flux and
                                 * scale outputs by F_sun_eq (solver is linear in
                                 * F_sun).  F_sun_eq inherits bit-level jitter from
                                 * the per-row atm solve (view node in the atm
                                 * grid), which perturbed the cache key every row;
                                 * with F pinned to 1.0 each beam solves ONCE per
                                 * grid and replays everywhere else. */
                                const double ocrt_beam_scale = F_sun_eq;
#else
                                const double ocrt_beam_scale = 1.0;
#endif
                                int rc_k = rt_water_rt_sos_pure(sza_air_eq_d,
                                                                  cs->vza_deg, raa_call,
                                                                  cs->wavelength_nm,
                                                                  T_C, S_gkg, n_w,
#ifdef OCRT_FAST_KERNELS
                                                                  1.0,
#else
                                                                  F_sun_eq,
#endif
                                                                  aw_lut, psi_T_lut,
                                                                  &w_opts_local, &r_k);
                                if (rc_k != 0) {
                                    if (sky_capture) {
                                        free(r_k.I_up_per_m); free(r_k.Q_up_per_m);
                                        free(r_k.U_up_per_m); free(r_k.mu_water_pos); free(r_k.w_water_pos);
                                    }
                                    continue;
                                }
#ifdef OCRT_FAST_KERNELS
                                r_k.I_0minus_view *= ocrt_beam_scale;
                                r_k.Q_0minus_view *= ocrt_beam_scale;
                                r_k.U_0minus_view *= ocrt_beam_scale;
                                r_k.I_0plus_view  *= ocrt_beam_scale;
                                r_k.Q_0plus_view  *= ocrt_beam_scale;
                                r_k.U_0plus_view  *= ocrt_beam_scale;
#else
                                (void)ocrt_beam_scale;
#endif
                                Lu_j_avg += r_k.I_0minus_view / N_raa_quad;
                                Q_j_avg  += r_k.Q_0minus_view / N_raa_quad;
                                U_j_avg  += r_k.U_0minus_view / N_raa_quad;
                                /* above-water (already T_wa-transmitted per call) */
                                La_j_avg += r_k.I_0plus_view / N_raa_quad;
                                Qa_j_avg += r_k.Q_0plus_view / N_raa_quad;
                                Ua_j_avg += r_k.U_0plus_view / N_raa_quad;
                                /* C3-full: accumulate this beam's m=0 0- upwelling field.
                                 * The field is already F_sun_eq/π-scaled (physical), so the
                                 * sum over beams = total skylight-induced 0- field.  m=0
                                 * slice = indices [0..nf-1] (stride n_mu = n_mu_water_filled). */
                                if (sky_capture && r_k.I_up_per_m && r_k.n_mu_water_filled > 0) {
                                    int nf = r_k.n_mu_water_filled;
                                    if (sky_n_mu_m0 == 0) {
                                        sky_up_I_m0 = (double*)calloc((size_t)nf, sizeof(double));
                                        sky_up_Q_m0 = (double*)calloc((size_t)nf, sizeof(double));
                                        sky_up_U_m0 = (double*)calloc((size_t)nf, sizeof(double));
                                        sky_mu_pos  = (double*)calloc((size_t)nf, sizeof(double));
                                        sky_w_pos   = (double*)calloc((size_t)nf, sizeof(double));
                                        if (sky_up_I_m0 && sky_up_Q_m0 && sky_up_U_m0 && sky_mu_pos) {
                                            sky_n_mu_m0 = nf;
                                            for (int kk = 0; kk < nf; ++kk) {
                                                sky_mu_pos[kk] = r_k.mu_water_pos[kk];
                                                if (sky_w_pos && r_k.w_water_pos) sky_w_pos[kk] = r_k.w_water_pos[kk];
                                            }
                                        }
                                    }
                                    if (sky_n_mu_m0 == nf) {
                                        for (int kk = 0; kk < nf; ++kk) {
                                            sky_up_I_m0[kk] += r_k.I_up_per_m[kk] * ocrt_beam_scale; /* m=0 */
                                            sky_up_Q_m0[kk] += r_k.Q_up_per_m[kk] * ocrt_beam_scale;
                                            sky_up_U_m0[kk] += r_k.U_up_per_m[kk] * ocrt_beam_scale;
                                        }
                                    }
                                }
                                if (sky_capture) {
                                    free(r_k.I_up_per_m); free(r_k.Q_up_per_m);
                                    free(r_k.U_up_per_m); free(r_k.mu_water_pos); free(r_k.w_water_pos);
                                }
                            }
                            Lu_diff_water_from_atm += Lu_j_avg;
                            Qu_diff_water_from_atm += Q_j_avg;
                            Uu_diff_water_from_atm += U_j_avg;
                            Lu_above_diff_from_atm += La_j_avg;
                            Qu_above_diff_from_atm += Qa_j_avg;
                            Uu_above_diff_from_atm += Ua_j_avg;
                        }
                    }
                    OCRT_S6T("[S6T]  beams end solved=%d skipped=%d t=%.4f\n", ocrt_s6_solved, ocrt_s6_skipped, ocrt_s6_now()-ocrt_s6_t0);
                    rt_aw_coupled_field_free(&coupled);
                }
            }
            rt_atm_boa_export_free(&boa);
        }
    }

    /* [v1.04+xsec_norm_cli+2, 2026-05-24] BUG-WRT-001 fix (옵션 A):
     *   F_sun_water = F_sun_BOA · μ_air · T_aw_flux / μ_water
     * water solver 에 air-water Fresnel transmission 이 반영된 F_sun 전달.
     * 호출 후 Ed_0plus_air 는 air-side BOA value 로 재정정 (denominator).
     * 상세: CHANGELOG_v1_04_2026-05-24_KST0240_water_rt_fixes.md §1.1 */
    /* PSSA scope boundary (2026-07-18): --pssa corrects the atmospheric
     * TOA->0+ direct-beam path only.  After refraction the in-water beam is
     * deliberately propagated with the existing plane-parallel water-column
     * geometry.  No underwater spherical-shell correction is applied here.
     * This is an explicit model-scope decision, documented in rt_ipss.h and
     * docs/PSSA_SCOPE_AND_RADIOMETRY_SPLIT_2026-07-18.md. */
    double F_sun_water = F_sun_BOA;
    {
        double mu_sun_air_corr = cos(cs->sza_deg * M_PI / 180.0);
        double mu_sun_water_corr = rt_air_water_mu_refracted_down(mu_sun_air_corr, n_w);
        if (mu_sun_air_corr > 0.0 && mu_sun_water_corr > 0.0) {
            double T_aw;
            if ((cs->surface == RT_SURFACE_OCEAN ||
                 cs->surface == RT_SURFACE_BLACK_FRESNEL_OCEAN) && cs->wind_speed > 0.0) {
                /* [2026-05-31] BUG-WRT-002 fix: rough-surface (Cox-Munk) direct-beam
                 * transmittance via slope-integrated facet Fresnel. Previously this
                 * path always used flat-Fresnel 1-R(SZA) regardless of wind, omitting
                 * the grazing transmission enhancement (-22% vs OSOAA at SZA=85,W=1).
                 * Validated vs OSOAA: exact <=80 deg, wind-dependence matched at 85.
                 * Consistent with the water->air Cox-Munk BTDF already used for wind>0. */
                T_aw = surface_T_aw_coxmunk_direct(mu_sun_air_corr, n_w,
                                                   cs->wind_speed, cs->sigma_type);
            } else {
                /* flat-Fresnel (FLAT surface, or wind<=0 almost-flat sea) */
                double M_R_aa[9];
                rt_air_water_R_aa(mu_sun_air_corr, n_w, cs->q_convention, M_R_aa);
                T_aw = 1.0 - M_R_aa[0];
            }
            if (T_aw < 0.0) T_aw = 0.0;
            if (T_aw > 1.0) T_aw = 1.0;
            F_sun_water = F_sun_BOA * mu_sun_air_corr * T_aw / mu_sun_water_corr;
            /* v1.10 A-FIX-2 (2026-07-12): the diffuse-top snapshot (ocrt_dt_*,
             * pi-normalized sky field) is consumed inside the water solve,
             * whose entire source/output scale is F_sun_water (the DIRECT-
             * beam converted flux above, incl. exp(-tau/mu0)*mu_air*T_aw/
             * mu_w).  The injected sky must be scaled by F_sun_TOA instead;
             * pre-multiply the snapshot by (F_sun_TOA / F_sun_water) so the
             * solve's common F_sun_water scaling lands it on F_sun_TOA.
             * Probe-proven: a single global rescale reproduces the
             * equivalent-beam Lu to 6 digits at sza40 AND sza80 (ledger
             * session 8).  Beams path unaffected (dt arrays bypassed). */
            if (ocrt_dt_I && F_sun_water > 0.0) {
                const double dtsc_ = F_sun_TOA / F_sun_water;
                size_t ndt_ = (size_t)(opts->fourier_m_max + 1) * (size_t)w_opts.n_mu_water;
                for (size_t z_ = 0; z_ < ndt_; ++z_) {
                    ocrt_dt_I[z_] *= dtsc_;
                    ocrt_dt_Q[z_] *= dtsc_;
                    if (ocrt_dt_U) ocrt_dt_U[z_] *= dtsc_;
                }
            }
        }
    }

    /* Call in-water SOS solver */
    OCRT_S6T("[S6T]  direct begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
    rt_water_rt_result_t w_res = {0};
    /* C3d: rigorous direct-beam water-leaving up-transmission (atm present only).
     * Expose the in-water 0- field so we can reverse-couple it onto the atm grid
     * and propagate it via a 2nd atm SOS bottom-source pass. */
    double TOA_wl_direct_I = 0.0, TOA_wl_direct_Q = 0.0, TOA_wl_direct_U = 0.0;
    int wl_rigorous_done = 0;
    if (atm_present) {
        /* FIX(2026-06-13): in-water solver can append up to 1 view node +
         * 3 zero-weight OSOAA slots (zero_slot_solve_mode==1 unconditionally
         * in step162) => internal n_mu can reach n_mu_water+4. The previous
         * "+2" under-allocated and rt_water_rt.c:3824 wrote past the end. */
        int _cM = w_opts.m_max_water + 1, _cN = w_opts.n_mu_water + 4 + ((w_opts.n_view_vza > 0) ? w_opts.n_view_vza : 0);   /* #20 */
        w_res.I_up_per_m  = (double*)calloc((size_t)_cM*_cN, sizeof(double));
        w_res.Q_up_per_m  = (double*)calloc((size_t)_cM*_cN, sizeof(double));
        w_res.U_up_per_m  = (double*)calloc((size_t)_cM*_cN, sizeof(double));
        w_res.mu_water_pos= (double*)calloc((size_t)_cN, sizeof(double));
        w_res.w_water_pos = (double*)calloc((size_t)_cN, sizeof(double));
        if (g_ocean_water_capture && g_ocean_water_capture->active) {
            w_res.I_view_per_m = (double*)calloc((size_t)_cM, sizeof(double));
            w_res.Q_view_per_m = (double*)calloc((size_t)_cM, sizeof(double));
            w_res.U_view_per_m = (double*)calloc((size_t)_cM, sizeof(double));
            w_res.I_air_view_per_m = (double*)calloc((size_t)_cM, sizeof(double));
            w_res.Q_air_view_per_m = (double*)calloc((size_t)_cM, sizeof(double));
            w_res.U_air_view_per_m = (double*)calloc((size_t)_cM, sizeof(double));
            w_res.view_m_max_capacity = _cM - 1;
        }
    }
    /* v1.10 D3-PROD (2026-07-10): OPT-IN water-stage swap.
     * env OCRT_D3_PROD=path to an operator file written by OCRT_D3_RB_OUT.
     * If the file's basis grid matches this run's ext_top injection grid
     * (node count + mu values to 1e-9) the atmospheric-diffuse injection is
     * REMOVED from the solve (beam-only water solve) and its contribution is
     * reconstructed afterwards by the operator application:
     *   up[m][k] += sum_cb sum_jb RB[cb][jb][cout][m][k] * d[cb][m][jb]
     *   Ed_0minus += sum_cb sum_jb EdRB[cb][jb] * d[cb][m=0][jb]
     * (field-linear contract, SURGERY_DESIGNS §D3-PROD; every quantity the
     * grid consumes is covered - certified by the D3-0/1/2 + Ed gates).
     * On ANY mismatch the swap is skipped and the classical path runs
     * unchanged.  Default (env unset): dead code. */
    int    d3p_on = 0, d3p_n = 0, d3p_Mc = 0;
    double *d3p_mu = NULL, *d3p_w = NULL, *d3p_Ed = NULL, *d3p_RB = NULL, *d3p_Vw = NULL;
    const double *d3p_dI = NULL, *d3p_dQ = NULL, *d3p_dU = NULL;
    int d3p_map[64];
    int    d3p_dn = 0, d3p_dm = 0;
    { const char *pp = getenv("OCRT_D3_PROD");
      if (pp && pp[0] && w_opts.ext_top_I && w_opts.ext_top_mu && w_opts.ext_top_n > 0) {
        FILE *fi = fopen(pp, "rb");
        if (fi) {
            int hdr[4] = {0,0,0,0};
            if (fread(hdr, sizeof(int), 4, fi) == 4 && hdr[0] == 0x44335242 && hdr[1] == 2) {
                d3p_n = hdr[2]; d3p_Mc = hdr[3];
                d3p_mu = (double*)malloc((size_t)d3p_n*sizeof(double));
                d3p_w  = (double*)malloc((size_t)d3p_n*sizeof(double));
                d3p_Ed = (double*)malloc((size_t)3*d3p_n*sizeof(double));
                d3p_Vw = (double*)malloc((size_t)3*d3p_n*3*sizeof(double));
                size_t rbsz = (size_t)3*d3p_n*3*(size_t)d3p_Mc*d3p_n;
                d3p_RB = (double*)malloc(rbsz*sizeof(double));
                int ok = d3p_mu && d3p_w && d3p_Ed && d3p_Vw && d3p_RB &&
                         fread(d3p_mu, sizeof(double), (size_t)d3p_n, fi) == (size_t)d3p_n &&
                         fread(d3p_w,  sizeof(double), (size_t)d3p_n, fi) == (size_t)d3p_n &&
                         fread(d3p_Ed, sizeof(double), (size_t)3*d3p_n, fi) == (size_t)3*d3p_n &&
                         fread(d3p_Vw, sizeof(double), (size_t)3*d3p_n*3, fi) == (size_t)3*d3p_n*3 &&
                         fread(d3p_RB, sizeof(double), rbsz, fi) == rbsz;
                if (ok && w_opts.ext_top_n <= d3p_n) {
                    /* injection grid may be a SUBSET of the assembly grid
                     * (coupling GL nodes vs internal grid with zero-weight
                     * slots): build the index map by mu matching. */
                    int match = 1;
                    for (int j = 0; j < w_opts.ext_top_n; ++j) {
                        int hit = -1;
                        for (int k = 0; k < d3p_n; ++k)
                            if (fabs(w_opts.ext_top_mu[j] - d3p_mu[k]) <= 1e-9) { hit = k; break; }
                        if (hit < 0) { match = 0; break; }
                        d3p_map[j] = hit;
                    }
                    if (match) {
                        d3p_dI = w_opts.ext_top_I; d3p_dQ = w_opts.ext_top_Q; d3p_dU = w_opts.ext_top_U;
                        d3p_dn = w_opts.ext_top_n; d3p_dm = w_opts.ext_top_m_max;
                        w_opts.ext_top_I = NULL; w_opts.ext_top_Q = NULL; w_opts.ext_top_U = NULL;
                        w_opts.ext_top_mu = NULL; w_opts.ext_top_n = 0; w_opts.ext_top_m_max = -1;
                        d3p_on = 1;
                    }
                }
            }
            fclose(fi);
        }
        if (!d3p_on) {
            fprintf(stderr, "[D3-PROD] mismatch: file n=%d Mc=%d vs ext_top_n=%d m_max=%d; ",
                    d3p_n, d3p_Mc, w_opts.ext_top_n, w_opts.ext_top_m_max);
            if (d3p_mu && w_opts.ext_top_mu && w_opts.ext_top_n == d3p_n) {
                for (int k = 0; k < d3p_n; ++k)
                    if (fabs(w_opts.ext_top_mu[k] - d3p_mu[k]) > 1e-9) {
                        fprintf(stderr, "first mu diff at k=%d: run=%.9f file=%.9f", k, w_opts.ext_top_mu[k], d3p_mu[k]);
                        break;
                    }
            }
            fprintf(stderr, " - classical path kept\n");
        }
      }
    }
    if (g_ocean_water_replay && g_ocean_water_replay->active) {
        if (ocrt_ocean_water_replay_store(g_ocean_water_replay,
                cs->sza_deg, cs->wavelength_nm, T_C, S_gkg, n_w, F_sun_water,
                &w_opts, w_water_grid) == 0) {
            g_ocean_water_replay->diffuse_top_mode = ocrt_dt_mode ? 1 : 0;
        }
    }
    int rc = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                    cs->wavelength_nm,
                                    T_C, S_gkg, n_w, F_sun_water,
                                    aw_lut, psi_T_lut, &w_opts, &w_res);
    double Kd_0minus_total = w_res.Kd_0minus;
    OCRT_S6T("[S6T]  direct end rc=%d t=%.4f\n", rc, ocrt_s6_now()-ocrt_s6_t0);
    if (d3p_on && rc == 0) {
        /* reconstruct the diffuse contribution from the operator */
        int McA = d3p_dm + 1; if (McA > d3p_Mc) McA = d3p_Mc;
        int McR = w_res.m_max_filled + 1; if (McA > McR) McA = McR;
        int nA = (w_res.n_mu_water_filled == d3p_n) ? d3p_n : 0;
        size_t plf = (size_t)d3p_Mc * (size_t)d3p_n;
        double upmax = 0.0, edadd = 0.0;
        /* v1.10 D3-PROD wiring-3 (2026-07-10): VIEW-SCALAR reconstruction.
         * Batch warm rows consume w_res.{I,Q,U}_0minus_view (computed inside
         * the solver from the cached BEAM-ONLY fields), so the mode-field
         * addition alone leaves every row's Lu without the diffuse share
         * (observed 17%).  Add the operator response AT THE VIEW NODE,
         * phi-synthesized with the solver's own exported convention
         * (rt_solver_reconstruct_phi/_sin, Dphi via rt_raa_to_atm_fourier_phi;
         * view node located by the standard air->water Snell cosine - the
         * batch grid carries the view directions as grid nodes). */
        int d3p_kv = -1;
        double d3p_vI[8] = {0}, d3p_vQ[8] = {0}, d3p_vU[8] = {0};
        {
            double mu_va = cos(cs->vza_deg * M_PI / 180.0);
            double mu_vw = sqrt(fmax(0.0, 1.0 - (1.0 - mu_va*mu_va) / (n_w * n_w)));
            for (int k = 0; k < w_res.n_mu_water_filled; ++k)
                if (fabs(w_res.mu_water_pos[k] - mu_vw) < 1e-9) { d3p_kv = k; break; }
        }
        if (nA > 0) {
            for (int cb = 0; cb < 3; ++cb) {
                const double *dsrc = (cb==0)? d3p_dI : (cb==1)? d3p_dQ : d3p_dU;
                if (!dsrc) continue;
                for (int jb = 0; jb < d3p_dn; ++jb) {
                    const int jf = d3p_map[jb];   /* injection node -> assembly column */
                    const double *col = d3p_RB + (((size_t)cb*d3p_n + jf)*3)*plf;
                    edadd += d3p_Ed[(size_t)cb*d3p_n + jf] * dsrc[(size_t)0*d3p_dn + jb];
                    for (int mm = 0; mm < McA; ++mm) {
                        double c0 = dsrc[(size_t)mm*d3p_dn + jb];
                        if (c0 == 0.0) continue;
                        for (int k = 0; k < nA; ++k) {
                            size_t tf = (size_t)mm*d3p_n + k;
                            size_t tr = (size_t)mm*nA + k;
                            w_res.I_up_per_m[tr] += c0 * col[0*plf + tf];
                            w_res.Q_up_per_m[tr] += c0 * col[1*plf + tf];
                            w_res.U_up_per_m[tr] += c0 * col[2*plf + tf];
                            if (k == d3p_kv && mm < 8) {
                                d3p_vI[mm] += c0 * col[0*plf + tf];
                                d3p_vQ[mm] += c0 * col[1*plf + tf];
                                d3p_vU[mm] += c0 * col[2*plf + tf];
                            }
                            double v = fabs(c0 * col[0*plf + tf]);
                            if (v > upmax) upmax = v;
                        }
                    }
                }
            }
            w_res.Ed_0minus_water += edadd;
            { static int tr1 = 0; if (!tr1++) { fprintf(stderr, "[D3tr] beamEd=%.9e dEd=%.9e afterEd=%.9e beamUp00=%.9e dUp00=%.9e\n",
                    w_res.Ed_0minus_water - edadd, edadd, w_res.Ed_0minus_water,
                    w_res.I_up_per_m[0] - 0.0, 0.0); fflush(stderr); } }
            if (d3p_kv >= 0) {
                double dphi = rt_raa_to_atm_fourier_phi(cs->raa_deg);
                double dvI = rt_solver_reconstruct_phi    (d3p_vI, McA - 1, dphi);
                double dvQ = rt_solver_reconstruct_phi    (d3p_vQ, McA - 1, dphi);
                double dvU =  rt_solver_reconstruct_phi_sin(d3p_vU, McA - 1, dphi);
                /* Canonical public reconstruction: I/Q use cosine modes and
                 * U uses the positive sine modes, identical to the main water
                 * view synthesis and the atmospheric output convention. */
                w_res.I_0minus_view += dvI;
                w_res.Q_0minus_view += dvQ;
                w_res.U_0minus_view += dvU;
                /* v1.10 D3-PROD (2026-07-10): 0+ view scalars.  The solver's
                 * 0- -> 0+ transform is LINEAR in the 0- view Stokes (both
                 * default branches apply the SAME exported 3x3,
                 * rt_air_water_T_wa at the refracted view direction; the
                 * spline branch composes linear ops).  Apply the pointwise
                 * Mueller to the reconstruction delta.  The 30-column gate
                 * adjudicates whether the run's active branch matches;
                 * spline-branch runs would show a residual here and are then
                 * explicitly unsupported rather than approximated. */
                { double mu_va2 = cos(cs->vza_deg * M_PI / 180.0);
                  double mu_vw2 = sqrt(fmax(0.0, 1.0 - (1.0 - mu_va2*mu_va2) / (n_w * n_w)));
                  double Mtw[9];
                  rt_air_water_T_wa(mu_vw2, n_w, cs->q_convention, Mtw);
                  w_res.I_0plus_view += Mtw[0]*dvI + Mtw[1]*dvQ + Mtw[2]*dvU;
                  w_res.Q_0plus_view += Mtw[3]*dvI + Mtw[4]*dvQ + Mtw[5]*dvU;
                  w_res.U_0plus_view += Mtw[6]*dvI + Mtw[7]*dvQ + Mtw[8]*dvU;
                }
            }
        }
        fprintf(stderr, "[D3-PROD] swap active: beam-only solve + operator reconstruction (modes %d, max|dup|=%.3e, dEd=%.6e, gridmatch=%d filled=%d file_n=%d kv=%d)\n",
                McA, upmax, edadd, nA > 0, w_res.n_mu_water_filled, d3p_n, d3p_kv);
    }
    free(d3p_mu); free(d3p_w); free(d3p_Ed); free(d3p_RB);
    if (rc != 0) {
        fprintf(stderr, "[B.4] rt_water_rt_sos_pure FAIL rc=%d\n", rc);
        return rc - 1000;  /* prefix for traceability */
    }

    /* v1.10 D3-0 (2026-07-10): R_B assembly driver - INSTRUMENTATION ONLY.
     * env OCRT_D3_ASSEMBLE=1 -> after the production solve, assemble the
     * water-body reflection-operator columns by unit-basis injections
     * through the ext_top channel (B-0c plumbing).  Bases: I and Q only;
     * U injection is NOT plumbed in ocrt_add_diffuse_top_primary yet -
     * explicitly deferred to D3-0b.  Basis response = solve(basis) -
     * solve(no-injection): the subtraction removes the direct-beam part
     * exactly (linearity), avoiding the f_scale=F_sun/pi=0 trap that
     * F_sun=0 would cause in the B.5 exposure scaling.  Responses are the
     * z=0- upwelling per-mode fields (I/Q/U_up_per_m).  LINEARITY GATE:
     * the response to a smooth test field must equal the basis
     * recomposition (reported to stderr).  Writes /tmp/d3_RB.bin.
     * Production outputs untouched (runs after them, files/stderr only). */
    if (getenv("OCRT_D3_ASSEMBLE") && rc == 0) {
        const double d3_t0 = ocrt_s6_now();
        const int Mc3 = w_opts.m_max_water + 1;
        const int cap3 = w_opts.n_mu_water + 4 +
                         ((w_opts.n_view_vza > 0) ? w_opts.n_view_vza : 0);
        #define D3_ALLOC(rr) do { memset(&(rr),0,sizeof(rr)); \
            (rr).I_up_per_m=(double*)calloc((size_t)Mc3*cap3,sizeof(double)); \
            (rr).Q_up_per_m=(double*)calloc((size_t)Mc3*cap3,sizeof(double)); \
            (rr).U_up_per_m=(double*)calloc((size_t)Mc3*cap3,sizeof(double)); \
            (rr).mu_water_pos=(double*)calloc((size_t)cap3,sizeof(double)); \
            (rr).w_water_pos=(double*)calloc((size_t)cap3,sizeof(double)); } while(0)
        #define D3_FREE(rr) do { free((rr).I_up_per_m); free((rr).Q_up_per_m); \
            free((rr).U_up_per_m); free((rr).mu_water_pos); free((rr).w_water_pos); } while(0)
        setenv("OCRT_EXTTOP_KERNEL", "table", 1);   /* D3-0c: driver solves use the table kernel */
        rt_water_rt_result_t d3_r0; D3_ALLOC(d3_r0);
        rt_water_rt_options_t d3_o = w_opts;
        /* operator assembly is a one-off precompute: force high SOS
         * precision so the operator (and its linearity residual) is not
         * limited by the production tolerance.  Gate scales with this. */
        d3_o.tolerance = 1e-11;
        d3_o.bypass_grid_cache = 1;
        { const char *te = getenv("OCRT_D3_TOL");   /* verification runs may loosen: gate = O(tol) proven */
          if (te && atof(te) > 0.0) d3_o.tolerance = atof(te); }
        d3_o.ext_top_I = NULL; d3_o.ext_top_Q = NULL; d3_o.ext_top_U = NULL;
        d3_o.ext_top_mu = NULL; d3_o.ext_top_n = 0; d3_o.ext_top_m_max = -1;
        int d3_rc0 = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                          cs->wavelength_nm, T_C, S_gkg, n_w,
                                          F_sun_water, aw_lut, psi_T_lut,
                                          &d3_o, &d3_r0);
        const int n3 = d3_r0.n_mu_water_filled;
        if (d3_rc0 == 0 && n3 > 0 && n3 <= cap3) {
            const size_t plane = (size_t)Mc3 * (size_t)n3;
            /* columns: [comp_in(0=I,1=Q)][node_in jb] -> response plane
             * per (m, k_out, comp_out in {I,Q,U}) */
            double *RB = (double*)calloc((size_t)3 * n3 * 3 * plane, sizeof(double));
            double *extI = (double*)calloc((size_t)Mc3 * n3, sizeof(double));
            double *extQ = (double*)calloc((size_t)Mc3 * n3, sizeof(double));
            double *extU = (double*)calloc((size_t)Mc3 * n3, sizeof(double));
            double *EdRB = (double*)calloc((size_t)3 * n3, sizeof(double));   /* D3-PROD contract: Ed response row */
            int d3_ok = (RB && extI && extQ && extU);
            int d3_solves = 0;
            int d3_nb = n3;   /* diag: cap basis nodes via env */
            { const char *e = getenv("OCRT_D3_NB"); if (e && atoi(e) > 0 && atoi(e) < n3) d3_nb = atoi(e); }
            for (int cb = 0; cb < 3 && d3_ok; ++cb)   /* D3-0b: I,Q,U */
              for (int jb = 0; jb < d3_nb && d3_ok; ++jb) {
                memset(extI, 0, (size_t)Mc3*n3*sizeof(double));
                memset(extQ, 0, (size_t)Mc3*n3*sizeof(double));
                memset(extU, 0, (size_t)Mc3*n3*sizeof(double));
                for (int mm = 0; mm < Mc3; ++mm)
                    (cb == 0 ? extI : cb == 1 ? extQ : extU)[(size_t)mm*n3 + jb] = 1.0;
                rt_water_rt_result_t d3_rb; D3_ALLOC(d3_rb);
                d3_o.ext_top_I = extI; d3_o.ext_top_Q = extQ;
                d3_o.ext_top_U = extU;   /* D3-0b */
                d3_o.ext_top_mu = d3_r0.mu_water_pos;
                d3_o.ext_top_n = n3; d3_o.ext_top_m_max = Mc3 - 1;
                int rcb = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg,
                              cs->raa_deg, cs->wavelength_nm, T_C, S_gkg,
                              n_w, F_sun_water, aw_lut, psi_T_lut,
                              &d3_o, &d3_rb);
                if (getenv("OCRT_D3_DIAG")) {
                    double q0=0, u0=0;
                    for (int mm2 = 0; mm2 < Mc3; ++mm2)
                      for (int k2 = 0; k2 < n3; ++k2) {
                        double vq = fabs(d3_rb.Q_up_per_m[(size_t)mm2*n3+k2] - d3_r0.Q_up_per_m[(size_t)mm2*n3+k2]);
                        double vu = fabs(d3_rb.U_up_per_m[(size_t)mm2*n3+k2] - d3_r0.U_up_per_m[(size_t)mm2*n3+k2]);
                        if (vq>q0) q0=vq; if (vu>u0) u0=vu;
                      }
                    fprintf(stderr, "[D3diagQU] cb=%d jb=%d dQ=%.3e dU=%.3e y0Qmax=%.3e\n", cb, jb, q0, u0,
                            fabs(d3_r0.Q_up_per_m[(size_t)1*n3+3]));
                    double m0=0, m1=0, y0m0=0, y0m1=0;
                    for (int mm2 = 0; mm2 < Mc3; ++mm2)
                      for (int k2 = 0; k2 < n3; ++k2) {
                        double v = fabs(d3_rb.I_up_per_m[(size_t)mm2*n3+k2]);
                        double w = fabs(d3_r0.I_up_per_m[(size_t)mm2*n3+k2]);
                        if (mm2==0) { if (v>m0) m0=v; if (w>y0m0) y0m0=w; }
                        else        { if (v>m1) m1=v; if (w>y0m1) y0m1=w; }
                      }
                    fprintf(stderr, "[D3diag] cb=%d jb=%d raw m0=%.3e m1+=%.3e | y0 m0=%.3e m1+=%.3e\n",
                            cb, jb, m0, m1, y0m0, y0m1);
                }
                if (rcb == 0 && d3_rb.n_mu_water_filled == n3 && EdRB)
                    EdRB[(size_t)cb*n3 + jb] = d3_rb.Ed_0minus_water - d3_r0.Ed_0minus_water;   /* D3-PROD: scalar Ed response per column */
                if (rcb == 0 && d3_rb.n_mu_water_filled == n3) {
                    double *col = RB + (((size_t)cb*n3 + jb) * 3) * plane;
                    for (size_t t = 0; t < plane; ++t) {
                        col[0*plane + t] = d3_rb.I_up_per_m[t] - d3_r0.I_up_per_m[t];
                        col[1*plane + t] = d3_rb.Q_up_per_m[t] - d3_r0.Q_up_per_m[t];
                        col[2*plane + t] = d3_rb.U_up_per_m[t] - d3_r0.U_up_per_m[t];
                    }
                    ++d3_solves;
                } else d3_ok = 0;
                D3_FREE(d3_rb);
              }
            /* LINEARITY GATE: smooth I-only test field d[j] */
            double d3_gate = -1.0;
            if (d3_ok) {
                double *d3_d = (double*)calloc((size_t)Mc3*n3, sizeof(double));
                for (int mm = 0; mm < Mc3; ++mm)
                    for (int j = 0; j < n3; ++j)
                        d3_d[(size_t)mm*n3+j] = 0.10 + 0.02*j;
                rt_water_rt_result_t d3_rt; D3_ALLOC(d3_rt);
                d3_o.ext_top_I = d3_d; d3_o.ext_top_Q = NULL; d3_o.ext_top_U = NULL;
                int rct = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg,
                              cs->raa_deg, cs->wavelength_nm, T_C, S_gkg,
                              n_w, F_sun_water, aw_lut, psi_T_lut,
                              &d3_o, &d3_rt);
                if (rct == 0) {
                    double mx = 0.0, ref = 0.0;
                    for (size_t t = 0; t < plane; ++t) {
                        double y = d3_rt.I_up_per_m[t] - d3_r0.I_up_per_m[t];
                        double yc = 0.0;
                        /* d is m-uniform: coefficient of node jb = d[jb] (comp I) */
                        for (int jb = 0; jb < n3; ++jb)
                            yc += (0.10 + 0.02*jb) *
                                  RB[(((size_t)0*n3 + jb)*3 + 0)*plane + t];
                        double e = fabs(y - yc);
                        if (fabs(y) > ref) ref = fabs(y);
                        if (e > mx) mx = e;
                    }
                    d3_gate = (ref > 0.0) ? mx / ref : 0.0;
                }
                D3_FREE(d3_rt); free(d3_d);
            }
            /* D3-0b GATES: (a) m=0 U-basis response must vanish (t_l m=0 is 0);
             * (b) reciprocity: with W_j = mu_j*w_j, the matrices
             * A(k,j) = I-response(k) to U-basis(j) and
             * B(k,j) = U-response(k) to I-basis(j) must satisfy
             * W_k A(k,j) = -W_j B(j,k) up to one global constant; we report
             * the max relative spread of the ratio over well-conditioned
             * entries (m>=1 planes). */
            double d3_m0u = -1.0, d3_rec = -1.0;
            if (d3_ok) {
                double mx0 = 0.0, ref0 = 0.0;
                for (int jb = 0; jb < n3; ++jb) {
                    const double *colU = RB + (((size_t)2*n3 + jb)*3)*plane;
                    const double *colI = RB + (((size_t)0*n3 + jb)*3)*plane;
                    for (int k = 0; k < n3; ++k) {
                        double v = fabs(colU[0*plane + (size_t)0*n3 + k]);
                        if (v > mx0) mx0 = v;
                        double w = fabs(colI[0*plane + (size_t)0*n3 + k]);
                        if (w > ref0) ref0 = w;
                    }
                }
                d3_m0u = (ref0 > 0.0) ? mx0 / ref0 : mx0;
                double rmin = 0.0, rmax = 0.0; int rn = 0;
                for (int mm = 1; mm < Mc3; ++mm)
                  for (int jb = 0; jb < n3; ++jb)
                    for (int k = 0; k < n3; ++k) {
                        double Wk = d3_r0.mu_water_pos[k] *
                                    (d3_r0.w_water_pos ? d3_r0.w_water_pos[k] : 1.0);
                        double Wj = d3_r0.mu_water_pos[jb] *
                                    (d3_r0.w_water_pos ? d3_r0.w_water_pos[jb] : 1.0);
                        double Aa = RB[(((size_t)2*n3 + jb)*3 + 0)*plane
                                       + (size_t)mm*n3 + k];
                        double Bb = RB[(((size_t)0*n3 + (size_t)k)*3 + 2)*plane
                                       + (size_t)mm*n3 + jb];
                        double num = Wk * Aa, den = -Wj * Bb;
                        if (fabs(den) < 1e-14 || fabs(num) < 1e-14) continue;
                        double r = num / den;
                        if (!rn) { rmin = r; rmax = r; }
                        else { if (r < rmin) rmin = r; if (r > rmax) rmax = r; }
                        ++rn;
                    }
                d3_rec = (rn > 0 && rmax != 0.0) ? (rmax - rmin) / fabs(rmax) : -2.0;
            }
            if (d3_ok) {
                FILE *fb = fopen("/tmp/d3_RB.bin", "wb");
                if (fb) {
                    int hdr[4] = { Mc3, n3, 3, 3 };
                    fwrite(hdr, sizeof(int), 4, fb);
                    fwrite(d3_r0.mu_water_pos, sizeof(double), (size_t)n3, fb);
                    fwrite(d3_r0.w_water_pos, sizeof(double), (size_t)n3, fb);   /* D3-0c: weights for gate scans */
                    fwrite(RB, sizeof(double), (size_t)3*n3*3*plane, fb);
                    fclose(fb);
                }
            }
            fprintf(stderr, "[D3-0] solves=%d n=%d Mc=%d lin=%.3e m0U=%.3e rec_spread=%.3e t=%.1fs ok=%d\n",
                    d3_solves, n3, Mc3, d3_gate, d3_m0u, d3_rec, ocrt_s6_now()-d3_t0, d3_ok);
            if (d3_ok && EdRB) {
                FILE *fe = fopen("/tmp/d3_EdRB.bin", "wb");
                if (fe) { int hdr[2]={3,n3}; fwrite(hdr,sizeof(int),2,fe);
                          fwrite(EdRB,sizeof(double),(size_t)3*n3,fe); fclose(fe); }
                double emx=0; for (size_t t2=0;t2<(size_t)3*n3;++t2) if (fabs(EdRB[t2])>emx) emx=fabs(EdRB[t2]);
                fprintf(stderr, "[D3-Ed] Ed response row recorded: max|dEd|=%.3e (y0 Ed=%.6e)\n", emx, d3_r0.Ed_0minus_water);
            }
            free(RB); free(extI); free(extQ); free(extU); free(EdRB);
        } else {
            fprintf(stderr, "[D3-0] baseline solve failed rc=%d n=%d\n", d3_rc0, n3);
        }
        unsetenv("OCRT_EXTTOP_KERNEL");
        D3_FREE(d3_r0);
        #undef D3_ALLOC
        #undef D3_FREE
    }

    /* C3d: rigorous direct-beam water-leaving up-transmission.  The in-water
     * 0- upwelling field (per Fourier mode, on the water mu grid) is reverse-
     * coupled (Snell + T_wa Mueller, FLAT) onto the atm mu grid, then injected
     * as the bottom source of a 2nd atm SOS pass (solar source OFF) so the atm
     * multiple scattering + attenuation act on it.  Skylight stays first-cut. */
    if (atm_present && w_res.I_up_per_m && w_res.n_mu_water_filled > 0) {
        int n_mu_gl  = w_opts.n_mu_water;   /* v1.10 B-0b.1: shared N */
        rt_uangles_t ocrt_u2; int ocrt_vd2 = -1;
        int ocrt_u2_ok = (rt_uangles_init(&ocrt_u2, n_mu_gl) == 0);
        if (ocrt_u2_ok && opts->view_as_node)
            ocrt_u2_ok = (rt_uangles_add(&ocrt_u2,
                              cos(cs->vza_deg * M_PI / 180.0), &ocrt_vd2) == 0);
        int n_mu_atm = ocrt_u2_ok ? ocrt_u2.n_total
                                  : n_mu_gl + (opts->view_as_node ? 1 : 0);
        /* v1.10 B-0a.3: pass2 bottom-source grid = the unified table
         * (must match the pass2 SOS internal ring index-for-index). */
        const int atm_mmax = opts->fourier_m_max;
        int water_src_mmax = atm_mmax;
        if (water_src_mmax > w_res.m_max_filled)
            water_src_mmax = w_res.m_max_filled;
        /* Historical diagnostics below use `mmax` for the modes actually
         * produced by the water solve.  The storage passed to the atmosphere,
         * however, is deliberately sized to atm_mmax and zero-padded above
         * water_src_mmax. */
        const int mmax = water_src_mmax;
        double *mu_atm = (double*)calloc((size_t)n_mu_atm, sizeof(double));
        double *w_atm  = (double*)calloc((size_t)n_mu_atm, sizeof(double));
        size_t na = (size_t)(atm_mmax + 1) * (size_t)n_mu_atm;
        double *wlI = (double*)calloc(na, sizeof(double));
        double *wlQ = (double*)calloc(na, sizeof(double));
        double *wlU = (double*)calloc(na, sizeof(double));
        int _wl_grid_ok = (mu_atm && w_atm && wlI && wlQ && wlU && ocrt_u2_ok);
        if (_wl_grid_ok) {
            for (int kU = 0; kU < n_mu_atm; ++kU) {
                mu_atm[kU] = ocrt_u2.mu[kU];
                w_atm [kU] = ocrt_u2.w [kU];
            }
        }
        /* view node lives at its sorted table position (B-0a.3). */
        OCRT_S6T("[S6T]  cpl begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
        if (_wl_grid_ok &&
            rt_air_water_couple_water_to_atm(
                w_res.I_up_per_m, w_res.Q_up_per_m, w_res.U_up_per_m,
                w_res.mu_water_pos, w_res.w_water_pos,
                w_res.n_mu_water_filled, water_src_mmax,
                n_w, cs->q_convention, cs->wind_speed, cs->sigma_type,
                mu_atm, n_mu_atm, w_atm, wlI, wlQ, wlU) == 0) {
            if (ocrt_debug_env("OCRT_DUMP_WLBOT")) {
                double bs_max=0, bs_sum=0; int kmax=-1;
                for (int k=0;k<n_mu_atm;k++){ double v=wlI[k]; if(fabs(v)>fabs(bs_max)){bs_max=v;kmax=k;} bs_sum+=v; }
                double iw_max=0;
                for (int k=0;k<w_res.n_mu_water_filled;k++){ double v=w_res.I_up_per_m[k]; if(fabs(v)>fabs(iw_max))iw_max=v; }
                fprintf(stderr,"[WLBOT] m=0 bottom_src(atm,n=%d): max=%.6e@mu=%.4f sum=%.6e | in-water 0- max(m=0)=%.6e | n_water_filled=%d\n",
                        n_mu_atm, bs_max, (kmax>=0?mu_atm[kmax]:-1.0), bs_sum, iw_max, w_res.n_mu_water_filled);
            }
            /* pass 2: atm SOS driven ONLY by this water-leaving bottom source */
            OCRT_S6T("[S6T]  cpl end t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            OCRT_S6T("[S6T]  wl2 begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            rt_case_t cs_atm2 = *cs;
            cs_atm2.surface = (cs->wind_speed > 0.0) ? RT_SURFACE_BLACK_FRESNEL_OCEAN : RT_SURFACE_FLAT;
            cs_atm2.F_sun   = F_sun_TOA;
            cs_atm2.n_water = n_w;
            rt_options_t opts_p2 = *opts;
            opts_p2.n_mu = n_mu_gl;   /* B-0b.1 shared N */
            OCRT_S6T("[S6T]  atm2(C3d) begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            opts_p2.bottom_source_only  = 1;
            rt_options_set_ext_bottom_source(&opts_p2, wlI, wlQ, wlU,
                                             n_mu_atm, atm_mmax);
            rt_result_t atm_res2 = {0};
            /* ---- v1.10 S7b: atm2(C3d) grid cache -------------------------
             * Same pattern as S7 (see atm1 site): the bottom source wlI/Q/U
             * comes from the #16 grid-cached water field via the reverse
             * coupling — row-invariant (verified by the FNV fingerprint in
             * the key).  One LUT-path solve per key (view_as_node OFF; grid
             * vza ARE the shared-N nodes under the S7 precondition), rows
             * replay TOA_wl_direct from the raw grids.  Memory-for-speed per
             * Jae's 2026-07-10 policy (grids + source snapshot kept in TLS).
             * FASTK-only; pristine per-row path when the gate is off. */
            int s7b_on = 0;
#ifdef OCRT_FAST_KERNELS
            { const char *s7e = getenv("OCRT_ATM_GRID_CACHE");
              s7b_on = (s7e && s7e[0] && s7e[0] != '0') &&
                       cs->water_view_vza_list && cs->n_water_view_vza >= 2 &&
                       getenv("OCRT_S7_RAA_LIST") != NULL &&
                       /* S7-D: 연직 부근 행별 폴백 (S7 게이트와 동일 조건) */
                       ocrt_s7_view_within_nodes(w_opts.n_mu_water,
                           cos(cs->vza_deg * M_PI / 180.0));
              /* 2026-07-14 M2: 노드 일치 전제 제거 (S7 atm1 지점 주석 참조). */ }
#endif
            /* S7b solves the atmosphere on the pure Gauss integration ring
             * (view_as_node=0).  The normal row path may carry one additional
             * zero-weight view node, so its external bottom source has a
             * different stride.  Build a second, Gauss-only projection for
             * the cache instead of pretending the row-shaped object matches
             * the LUT solver.  This also makes the fingerprint independent of
             * whichever VZA happened to trigger the cache fill. */
            double *s7b_owned_I = NULL, *s7b_owned_Q = NULL, *s7b_owned_U = NULL;
            const double *s7b_src_I = wlI, *s7b_src_Q = wlQ, *s7b_src_U = wlU;
            int s7b_src_n_mu = opts_p2.ext_bottom_n_mu;
            int s7b_src_m_max = opts_p2.ext_bottom_m_max;
            if (s7b_on && s7b_src_n_mu != n_mu_gl) {
                rt_uangles_t s7b_gl;
                const size_t n7 = (size_t)(atm_mmax + 1) * (size_t)n_mu_gl;
                s7b_owned_I = (double*)calloc(n7, sizeof(double));
                s7b_owned_Q = (double*)calloc(n7, sizeof(double));
                s7b_owned_U = (double*)calloc(n7, sizeof(double));
                if (rt_uangles_init(&s7b_gl, n_mu_gl) != 0 ||
                    !s7b_owned_I || !s7b_owned_Q || !s7b_owned_U ||
                    rt_air_water_couple_water_to_atm(
                        w_res.I_up_per_m, w_res.Q_up_per_m, w_res.U_up_per_m,
                        w_res.mu_water_pos, w_res.w_water_pos,
                        w_res.n_mu_water_filled, water_src_mmax,
                        n_w, cs->q_convention, cs->wind_speed, cs->sigma_type,
                        s7b_gl.mu, n_mu_gl, s7b_gl.w,
                        s7b_owned_I, s7b_owned_Q, s7b_owned_U) != 0) {
                    free(s7b_owned_I); free(s7b_owned_Q); free(s7b_owned_U);
                    s7b_owned_I = s7b_owned_Q = s7b_owned_U = NULL;
                    s7b_on = 0;  /* safe row-shaped fallback below */
                } else {
                    s7b_src_I = s7b_owned_I;
                    s7b_src_Q = s7b_owned_Q;
                    s7b_src_U = s7b_owned_U;
                    s7b_src_n_mu = n_mu_gl;
                    s7b_src_m_max = atm_mmax;
                }
            }

            ocrt_s7b_cache_t *s7bp = &g_s7b;
            unsigned long long s7b_h = 1469598103934665603ULL;
            if (s7b_on) {
                /* Hash exactly the shaped object consumed by the cache solve.
                 * The high atmospheric modes are explicit zero padding when
                 * the water solve has a smaller Fourier cap. */
                size_t nsrc = (size_t)(s7b_src_m_max + 1) *
                              (size_t)s7b_src_n_mu;
                const double *srcs[3] = { s7b_src_I, s7b_src_Q, s7b_src_U };
                s7b_h = (s7b_h ^ (unsigned long long)s7b_src_n_mu) *
                        1099511628211ULL;
                s7b_h = (s7b_h ^ (unsigned long long)(s7b_src_m_max + 1)) *
                        1099511628211ULL;
                for (int q3 = 0; q3 < 3; ++q3)
                    for (size_t k = 0; k < nsrc; ++k) {
                        unsigned long long b;
                        memcpy(&b, &srcs[q3][k], 8);
                        s7b_h = (s7b_h ^ b) * 1099511628211ULL;
                    }
            }
            int s7b_iv = -1, s7b_ir = -1;
            if (s7b_on && s7bp->valid && s7bp->src_fnv == s7b_h &&
                s7bp->src_n_mu == s7b_src_n_mu &&
                s7bp->src_m_max == s7b_src_m_max &&
                s7bp->lam == cs->wavelength_nm && s7bp->sza == cs->sza_deg &&
                s7bp->F == F_sun_TOA && s7bp->wind == cs->wind_speed &&
                s7bp->nw == n_w && s7bp->surface == (int)cs_atm2.surface &&
                s7bp->qc == cs->q_convention && s7bp->sigt == cs->sigma_type &&
                s7bp->n_vza == cs->n_water_view_vza) {
                for (int q = 0; q < s7bp->n_vza; ++q)
                    if (fabs(s7bp->vza[q] - cs->vza_deg) < 1e-9) { s7b_iv = q; break; }
                for (int q = 0; q < s7bp->n_raa; ++q)
                    if (fabs(s7bp->raa[q] - cs->raa_deg) < 1e-9) { s7b_ir = q; break; }
            }
            if (s7b_iv >= 0 && s7b_ir >= 0) {
                size_t g = (size_t)s7b_iv * (size_t)s7bp->n_raa + (size_t)s7b_ir;
                TOA_wl_direct_I = s7bp->raw[0][g];
                TOA_wl_direct_Q = s7bp->raw[1][g];
                TOA_wl_direct_U = s7bp->raw[2][g];
                wl_rigorous_done = 1;
                OCRT_S6T("[S6T]  atm2(C3d) s7b-hit t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            } else if (s7b_on) {
                int nv = cs->n_water_view_vza, nr = 0;
                { const char *rl = getenv("OCRT_S7_RAA_LIST");
                  char tmp[4096]; strncpy(tmp, rl, sizeof tmp - 1); tmp[sizeof tmp-1]=0;   /* M2: raa 2.5도 격자 수용 */
                  free(s7bp->raa); s7bp->raa = (double*)calloc(512, sizeof(double));
                  for (char *p = strtok(tmp, ","); p && nr < 512; p = strtok(NULL, ","))
                      s7bp->raa[nr++] = atof(p); }
                free(s7bp->vza); s7bp->vza = (double*)calloc((size_t)nv, sizeof(double));
                memcpy(s7bp->vza, cs->water_view_vza_list, (size_t)nv * sizeof(double));
                for (int q3 = 0; q3 < 3; ++q3) {
                    free(s7bp->raw[q3]); free(s7bp->rho[q3]); free(s7bp->pm[q3]);
                    s7bp->raw[q3] = (double*)calloc((size_t)nv * (size_t)nr, sizeof(double));
                    s7bp->rho[q3] = (double*)calloc((size_t)nv * (size_t)nr, sizeof(double));
                    s7bp->pm[q3] = (double*)calloc((size_t)(opts_p2.fourier_m_max + 1) * (size_t)nv, sizeof(double));
                }
                rt_lut_grid_out_t l2; memset(&l2, 0, sizeof l2);
                l2.n_vza = nv; l2.vza_deg = s7bp->vza;
                l2.n_raa = nr; l2.raa_deg = s7bp->raa;
                l2.rho_I = s7bp->rho[0]; l2.rho_Q = s7bp->rho[1]; l2.rho_U = s7bp->rho[2];
                l2.raw_I = s7bp->raw[0]; l2.raw_Q = s7bp->raw[1]; l2.raw_U = s7bp->raw[2];
                l2.view_per_m_I = s7bp->pm[0]; l2.view_per_m_Q = s7bp->pm[1]; l2.view_per_m_U = s7bp->pm[2];
                l2.view_m_max = opts_p2.fourier_m_max;
                rt_options_t op2c = opts_p2;
                op2c.view_as_node = 0;
                rt_options_set_ext_bottom_source(&op2c,
                                                 s7b_src_I, s7b_src_Q, s7b_src_U,
                                                 s7b_src_n_mu, s7b_src_m_max);
                if (rt_solve_case_pol_impl(&cs_atm2, &op2c, aer, &atm_res2,
                                           &l2, NULL) == 0) {
                    s7bp->valid = 1; s7bp->n_vza = nv; s7bp->n_raa = nr; s7bp->m_max = opts_p2.fourier_m_max;
                    s7bp->src_fnv = s7b_h;
                    s7bp->src_n_mu = s7b_src_n_mu;
                    s7bp->src_m_max = s7b_src_m_max;
                    s7bp->lam = cs->wavelength_nm; s7bp->sza = cs->sza_deg;
                    s7bp->F = F_sun_TOA; s7bp->wind = cs->wind_speed; s7bp->nw = n_w;
                    s7bp->surface = (int)cs_atm2.surface;
                    s7bp->qc = cs->q_convention; s7bp->sigt = cs->sigma_type;
                    int iv2 = -1, ir2 = -1;
                    for (int q = 0; q < nv; ++q)
                        if (fabs(s7bp->vza[q] - cs->vza_deg) < 1e-9) { iv2 = q; break; }
                    for (int q = 0; q < nr; ++q)
                        if (fabs(s7bp->raa[q] - cs->raa_deg) < 1e-9) { ir2 = q; break; }
                    if (iv2 >= 0 && ir2 >= 0) {
                        size_t g2 = (size_t)iv2 * (size_t)nr + (size_t)ir2;
                        TOA_wl_direct_I = s7bp->raw[0][g2];
                        TOA_wl_direct_Q = s7bp->raw[1][g2];
                        TOA_wl_direct_U = s7bp->raw[2][g2];
                        wl_rigorous_done = 1;
                    }
                } else s7bp->valid = 0;
            } else {
                /* v1.10 S17 (2026-07-12, bit-safe speed): fold the raa axis of
                 * the atm2 (C3d) pass.  In lut-grid batches rows iterate raa
                 * inside a fixed vza; the atm2 solve depends on raa ONLY via
                 * the final phi reconstruction.  Memoize the per-m view
                 * samples for the last (vza | case scalars) and, on a hit,
                 * redo ONLY the reconstruction with the SAME functions and
                 * inputs the impl uses => bit-identical outputs (gated:
                 * production batch byte-compare, ledger S17).  Miss path is
                 * the untouched original call.  Disable: OCRT_S17_OFF=1. */
                {
                    ocrt_s17_cache_t *s17p = &g_s17;
                    /* The old S17 key omitted the external bottom-source
                     * field.  A discarded base solve could therefore seed the
                     * cache with a different water-leaving field and corrupt
                     * the first full-grid row.  Hash the complete shaped
                     * source, as S7b does. */
                    unsigned long long s17_h = 1469598103934665603ULL;
                    s17_h = ocrt_fnv1a_bytes(s17_h, &opts_p2.ext_bottom_n_mu,
                                             sizeof(opts_p2.ext_bottom_n_mu));
                    s17_h = ocrt_fnv1a_bytes(s17_h, &opts_p2.ext_bottom_m_max,
                                             sizeof(opts_p2.ext_bottom_m_max));
                    if (opts_p2.ext_bottom_n_mu > 0 && opts_p2.ext_bottom_m_max >= 0) {
                        const size_t ns17 = (size_t)opts_p2.ext_bottom_n_mu *
                                             (size_t)(opts_p2.ext_bottom_m_max + 1);
                        if (opts_p2.ext_bottom_per_m_I)
                            s17_h = ocrt_fnv1a_bytes(s17_h, opts_p2.ext_bottom_per_m_I,
                                                    ns17 * sizeof(double));
                        if (opts_p2.ext_bottom_per_m_Q)
                            s17_h = ocrt_fnv1a_bytes(s17_h, opts_p2.ext_bottom_per_m_Q,
                                                    ns17 * sizeof(double));
                        if (opts_p2.ext_bottom_per_m_U)
                            s17_h = ocrt_fnv1a_bytes(s17_h, opts_p2.ext_bottom_per_m_U,
                                                    ns17 * sizeof(double));
                    }
                    const int s17_off = (getenv("OCRT_S17_OFF") != NULL);
                    const int s17_hit = (!s17_off && s17p->valid &&
                        s17p->vza == cs->vza_deg && s17p->sza == cs->sza_deg &&
                        s17p->wl == cs->wavelength_nm && s17p->wind == cs->wind_speed &&
                        s17p->F == F_sun_TOA && s17p->m_max == opts_p2.fourier_m_max &&
                        s17p->src_n_mu == opts_p2.ext_bottom_n_mu &&
                        s17p->src_m_max == opts_p2.ext_bottom_m_max &&
                        s17p->src_fnv == s17_h &&
                        opts_p2.fourier_m_max < 33);
                    if (s17_hit) {
                        const double s17_dphi = cs->raa_deg * M_PI / 180.0;
                        TOA_wl_direct_I = rt_solver_reconstruct_phi    (s17p->pmI, s17p->m_max, s17_dphi);
                        TOA_wl_direct_Q = rt_solver_reconstruct_phi    (s17p->pmQ, s17p->m_max, s17_dphi);
                        TOA_wl_direct_U = rt_solver_reconstruct_phi_sin(s17p->pmU, s17p->m_max, s17_dphi);
                        wl_rigorous_done = 1;
                        goto s17_done;
                    }
                    if (!s17_off && opts_p2.fourier_m_max < 33) {
                        opts_p2.toa_view_per_m_I = s17p->pmI;
                        opts_p2.toa_view_per_m_Q = s17p->pmQ;
                        opts_p2.toa_view_per_m_U = s17p->pmU;
                    }
                    if (rt_solve_case_pol_for_ocean(&cs_atm2, &opts_p2, aer, &atm_res2, NULL) == 0) {
                        if (opts_p2.toa_view_per_m_I) {
                            s17p->valid = 1; s17p->m_max = opts_p2.fourier_m_max;
                            s17p->vza = cs->vza_deg; s17p->sza = cs->sza_deg;
                            s17p->wl = cs->wavelength_nm; s17p->wind = cs->wind_speed;
                            s17p->F = F_sun_TOA;
                            s17p->src_n_mu = opts_p2.ext_bottom_n_mu;
                            s17p->src_m_max = opts_p2.ext_bottom_m_max;
                            s17p->src_fnv = s17_h;
                        }
                    TOA_wl_direct_I = atm_res2.I_TOA;
                    TOA_wl_direct_Q = atm_res2.Q_TOA;
                    TOA_wl_direct_U = atm_res2.U_TOA;
                    wl_rigorous_done = 1;
                    }
                    s17_done: ;
                }
                OCRT_S6T("[S6T]  wl2 end t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            }
            /* A cache-fill failure must never drop the water-leaving source.
             * Re-run the exact row-shaped path as a correctness fallback. */
            if (s7b_on && !wl_rigorous_done) {
                rt_result_t atm_res2_fb = {0};
                if (rt_solve_case_pol_for_ocean(&cs_atm2, &opts_p2, aer,
                                                 &atm_res2_fb, NULL) == 0) {
                    TOA_wl_direct_I = atm_res2_fb.I_TOA;
                    TOA_wl_direct_Q = atm_res2_fb.Q_TOA;
                    TOA_wl_direct_U = atm_res2_fb.U_TOA;
                    wl_rigorous_done = 1;
                }
            }
            free(s7b_owned_I); free(s7b_owned_Q); free(s7b_owned_U);
                /* v1.10 D3-PROD wiring-2 (2026-07-10): the D3-1a/PP/FAC
                 * instrumentation blocks were hosted INSIDE the direct
                 * (for_ocean) branch above and never fired on batch
                 * (lut-raa) runs, which take the S7b cached branch.
                 * Relocated here - after the branch join, still inside
                 * the cs_atm2/opts_p2/wl* scope (before their frees) -
                 * so BOTH run types can assemble.  Fire-once guard stops
                 * per-row re-assembly in batches (rows share one
                 * water/atm configuration; the first row's operator is
                 * valid for all rows).  Inner blocks remain env-gated,
                 * default-inert. */
                { static int d3_blocks_fired = 0;
                  if (!d3_blocks_fired && (getenv("OCRT_D3_ATMOP") || getenv("OCRT_D3_PP"))) {
                    d3_blocks_fired = 1;
                    /* D3 driver solves must NOT hit the water grid cache
                     * (batch runs replay the production solve into every
                     * basis column - observed: 153 cols in 0.1 s, all-zero
                     * operator).  Save and disable for the whole block;
                     * restored below. */
                    const char *d3_gc_save = getenv("OCRT_WATER_GRID_CACHE");
                    setenv("OCRT_WATER_GRID_CACHE", "0", 1);
                    /* v1.10 D3-1a (2026-07-10): R_A* assembly driver -
                     * INSTRUMENTATION ONLY (env OCRT_D3_ATMOP=1).  Unit
                     * upward-radiance bases through the SAME ext_bottom
                     * channel this C3d pass uses (bottom_source_only=1:
                     * no solar source, so columns need no subtraction; a
                     * zero-injection purity solve is still gated).  BOA
                     * downward responses recorded via rt_atm_boa_export_t.
                     * The solver's INTERNAL node count may exceed n_mu_gl
                     * (inserted view node on single-point runs) and the
                     * export fill requires an exact n_mu match, so the
                     * driver detects it by trying n_mu_gl..+4 with the
                     * actual wl injection until the export fills.
                     * Production outputs untouched. */
                    if (getenv("OCRT_D3_ATMOP")) {
                        const int Mc_a = mmax + 1;   /* CLAMPED mode count = the wl/ext_bottom grid convention */
                        const int n_in = n_mu_atm;             /* injection grid = pass-2 bottom-source grid (GL + inserted view) */
                        const size_t pl_in = (size_t)Mc_a * (size_t)n_in;
                        int n_det = -1;                        /* detected internal n_mu */
                        const size_t pl_max = (size_t)Mc_a * (size_t)(n_in + 4);
                        double *exI = (double*)calloc(pl_in, sizeof(double));
                        double *exQ = (double*)calloc(pl_in, sizeof(double));
                        double *exU = (double*)calloc(pl_in, sizeof(double));
                        rt_atm_boa_export_t ex; memset(&ex, 0, sizeof ex);
                        ex.m_max = Mc_a - 1;
                        ex.mu_quad_pos = (double*)calloc((size_t)(n_in + 4), sizeof(double));
                        ex.I_per_m = (double*)calloc(pl_max, sizeof(double));
                        ex.Q_per_m = (double*)calloc(pl_max, sizeof(double));
                        ex.U_per_m = (double*)calloc(pl_max, sizeof(double));
                        int okA = (exI && exQ && exU && ex.mu_quad_pos &&
                                   ex.I_per_m && ex.Q_per_m && ex.U_per_m);
                        if (okA) {
                            for (int tryn = n_in; tryn <= n_in + 4; ++tryn) {
                                ex.n_mu = tryn;
                                memset(ex.I_per_m, 0, pl_max * sizeof(double));
                                rt_options_t ob = opts_p2;
                                rt_options_set_ext_bottom_source(&ob, wlI, wlQ, wlU,
                                                                 n_in, atm_mmax);
                                rt_result_t rbp = {0};
                                if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rbp, &ex) != 0) { okA = 0; break; }
                                double v = 0.0;
                                for (size_t t2 = 0; t2 < (size_t)Mc_a * (size_t)tryn; ++t2) v += fabs(ex.I_per_m[t2]);
                                if (v > 0.0) { n_det = tryn; break; }
                            }
                            if (n_det < 0) okA = 0;
                        }
                        const size_t pl_det = (okA) ? (size_t)Mc_a * (size_t)n_det : 0;
                        double pur = -1.0;
                        if (okA) {
                            memset(exI,0,pl_in*sizeof(double)); memset(exQ,0,pl_in*sizeof(double)); memset(exU,0,pl_in*sizeof(double));
                            rt_options_t ob = opts_p2;
                            rt_options_set_ext_bottom_source(&ob, exI, exQ, exU,
                                                             n_in, mmax);
                            rt_result_t rb0 = {0};
                            memset(ex.I_per_m,0,pl_max*sizeof(double)); memset(ex.Q_per_m,0,pl_max*sizeof(double)); memset(ex.U_per_m,0,pl_max*sizeof(double));
                            if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rb0, &ex) == 0) {
                                pur = 0.0;
                                for (size_t t2 = 0; t2 < pl_det; ++t2) {
                                    double v = fabs(ex.I_per_m[t2]) + fabs(ex.Q_per_m[t2]) + fabs(ex.U_per_m[t2]);
                                    if (v > pur) pur = v;
                                }
                            }
                        }
                        double *RA = okA ? (double*)calloc((size_t)3*n_in*3*pl_det, sizeof(double)) : NULL;
                        double *y0I = okA ? (double*)calloc(pl_det, sizeof(double)) : NULL;
                        double *y0Q = okA ? (double*)calloc(pl_det, sizeof(double)) : NULL;
                        double *y0U = okA ? (double*)calloc(pl_det, sizeof(double)) : NULL;
                        if (okA && (!RA || !y0I || !y0Q || !y0U)) okA = 0;
                        if (okA) { memcpy(y0I, ex.I_per_m, pl_det*sizeof(double));
                                   memcpy(y0Q, ex.Q_per_m, pl_det*sizeof(double));
                                   memcpy(y0U, ex.U_per_m, pl_det*sizeof(double)); }
                        int nsv = 0;
                        for (int cb = 0; cb < 3 && okA; ++cb)
                          for (int jb = 0; jb < n_in && okA; ++jb) {
                            memset(exI,0,pl_in*sizeof(double)); memset(exQ,0,pl_in*sizeof(double)); memset(exU,0,pl_in*sizeof(double));
                            double *tgt = (cb==0)?exI:(cb==1)?exQ:exU;
                            for (int mm = 0; mm < Mc_a; ++mm) tgt[(size_t)mm*n_in + jb] = 1.0;
                            rt_options_t ob = opts_p2;
                            rt_options_set_ext_bottom_source(&ob, exI, exQ, exU,
                                                             n_in, mmax);
                            rt_result_t rb = {0};
                            memset(ex.I_per_m,0,pl_max*sizeof(double)); memset(ex.Q_per_m,0,pl_max*sizeof(double)); memset(ex.U_per_m,0,pl_max*sizeof(double));
                            if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rb, &ex) == 0) {
                                double *col = RA + (((size_t)cb*n_in + jb)*3)*pl_det;
                                for (size_t t2 = 0; t2 < pl_det; ++t2) {
                                    col[0*pl_det+t2] = ex.I_per_m[t2] - y0I[t2];
                                    col[1*pl_det+t2] = ex.Q_per_m[t2] - y0Q[t2];
                                    col[2*pl_det+t2] = ex.U_per_m[t2] - y0U[t2];
                                }
                                ++nsv;
                            } else okA = 0;
                          }
                        double linA = -1.0;
                        if (okA) {
                            rt_options_t ob = opts_p2;
                            rt_options_set_ext_bottom_source(&ob, wlI, wlQ, wlU,
                                                             n_in, atm_mmax);
                            rt_result_t rb = {0};
                            memset(ex.I_per_m,0,pl_max*sizeof(double)); memset(ex.Q_per_m,0,pl_max*sizeof(double)); memset(ex.U_per_m,0,pl_max*sizeof(double));
                            if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rb, &ex) == 0) {
                                double mx=0.0, rf=0.0;
                                for (size_t t2 = 0; t2 < pl_det; ++t2) {
                                    double yI=0, yQ=0, yU=0;
                                    int mm2 = (int)(t2 / (size_t)n_det);
                                    for (int jb = 0; jb < n_in; ++jb) {
                                        double cI = wlI[(size_t)mm2*n_in+jb], cQ = wlQ[(size_t)mm2*n_in+jb], cU = wlU[(size_t)mm2*n_in+jb];
                                        yI += cI*RA[(((size_t)0*n_in+jb)*3+0)*pl_det+t2] + cQ*RA[(((size_t)1*n_in+jb)*3+0)*pl_det+t2] + cU*RA[(((size_t)2*n_in+jb)*3+0)*pl_det+t2];
                                        yQ += cI*RA[(((size_t)0*n_in+jb)*3+1)*pl_det+t2] + cQ*RA[(((size_t)1*n_in+jb)*3+1)*pl_det+t2] + cU*RA[(((size_t)2*n_in+jb)*3+1)*pl_det+t2];
                                        yU += cI*RA[(((size_t)0*n_in+jb)*3+2)*pl_det+t2] + cQ*RA[(((size_t)1*n_in+jb)*3+2)*pl_det+t2] + cU*RA[(((size_t)2*n_in+jb)*3+2)*pl_det+t2];
                                    }
                                    double eI=fabs(ex.I_per_m[t2]-y0I[t2]-yI), eQ=fabs(ex.Q_per_m[t2]-y0Q[t2]-yQ), eU=fabs(ex.U_per_m[t2]-y0U[t2]-yU);
                                    double v=fabs(ex.I_per_m[t2])+fabs(ex.Q_per_m[t2])+fabs(ex.U_per_m[t2]);
                                    if (v>rf) rf=v;
                                    if (eI>mx) mx=eI; if (eQ>mx) mx=eQ; if (eU>mx) mx=eU;
                                }
                                linA = (rf>0.0)? mx/rf : mx;
                            }
                            FILE *fa = fopen("/tmp/d3_RA.bin","wb");
                            if (fa) { int hdr[4]={Mc_a,n_det,3,n_in}; fwrite(hdr,sizeof(int),4,fa);
                                      fwrite(RA,sizeof(double),(size_t)3*n_in*3*pl_det,fa); fclose(fa); }
                        }
                        fprintf(stderr, "[D3-1a] RA solves=%d n_in=%d n_det=%d Mc=%d purity=%.3e linA=%.3e ok=%d\n",
                                nsv, n_in, n_det, Mc_a, pur, linA, okA);
                        free(exI); free(exQ); free(exU); free(RA); free(y0I); free(y0Q); free(y0U);
                        free(ex.mu_quad_pos); free(ex.I_per_m); free(ex.Q_per_m); free(ex.U_per_m);
                    }
                    /* v1.10 D3-1b (2026-07-10): PING-PONG exact multibounce
                     * reference - INSTRUMENTATION ONLY (env OCRT_D3_PP=1).
                     * Repeats the production coupling chain to convergence:
                     *   a_dn^b   = atm2(ext_bottom = wl^{b-1}) BOA export
                     *   d^b      = rt_air_water_couple_atm_to_water(a_dn^b)
                     *   up^b     = water solve(ext_top = d^b) - water solve(0)
                     *   wl^b     = rt_air_water_couple_water_to_atm(up^b)
                     *   dTOA^b   = atm2(ext_bottom = wl^b).I_TOA
                     * Every stage is the SAME production transform/solver;
                     * scales are chain-consistent by construction.  Water
                     * solves run with the D3-0c table injection kernel
                     * (setenv/unsetenv around, #22-style live getenv inside).
                     * Reported: per-bounce TOA increments relative to the
                     * production water-leaving TOA term.  Production outputs
                     * untouched. */
                    if (getenv("OCRT_D3_PP")) {
                        /* v1.10 D3-PROD assembly-cost (2026-07-10): rbonly
                         * file-export runs need only p0 (subtraction base);
                         * the 3 ping-pong bounces are gate instrumentation -
                         * skip them to cut per-assembly fixed cost. */
                        const char *ppfm = getenv("OCRT_D3_FAC");
                        const int PPB = (ppfm && !strcmp(ppfm, "rbonly")) ? 0 : 3;
                        const int Mc_p = mmax + 1;
                        const size_t pla = (size_t)Mc_p * (size_t)n_mu_atm;
                        rt_atm_boa_export_t px; memset(&px, 0, sizeof px);
                        px.m_max = Mc_p - 1; px.n_mu = n_mu_atm;      /* pass-2 grid already includes the view node (D3-1a: n_det==n_mu_atm) */
                        const size_t plx = (size_t)Mc_p * (size_t)(n_mu_atm + 4);
                        px.mu_quad_pos = (double*)calloc((size_t)(n_mu_atm + 4), sizeof(double));
                        px.I_per_m = (double*)calloc(plx, sizeof(double));
                        px.Q_per_m = (double*)calloc(plx, sizeof(double));
                        px.U_per_m = (double*)calloc(plx, sizeof(double));
                        double *wlb_I = (double*)calloc(pla, sizeof(double));
                        double *wlb_Q = (double*)calloc(pla, sizeof(double));
                        double *wlb_U = (double*)calloc(pla, sizeof(double));
                        memcpy(wlb_I, wlI, pla*sizeof(double));
                        memcpy(wlb_Q, wlQ, pla*sizeof(double));
                        memcpy(wlb_U, wlU, pla*sizeof(double));
                        double *wl0_I = (double*)calloc(pla, sizeof(double));
                        double *wl0_Q = (double*)calloc(pla, sizeof(double));
                        double *wl0_U = (double*)calloc(pla, sizeof(double));
                        double *b1_I  = (double*)calloc(pla, sizeof(double));
                        double *b1_Q  = (double*)calloc(pla, sizeof(double));
                        double *b1_U  = (double*)calloc(pla, sizeof(double));
                        if (wl0_I) memcpy(wl0_I, wlI, pla*sizeof(double));
                        if (wl0_Q) memcpy(wl0_Q, wlQ, pla*sizeof(double));
                        if (wl0_U) memcpy(wl0_U, wlU, pla*sizeof(double));
                        double *WopKeep = NULL;   /* D3-3: whole-chain W handed out of the OP block for the FAC gate */
                        rt_aw_coupled_field_t cf; memset(&cf, 0, sizeof cf);
                        int okP = (px.mu_quad_pos && px.I_per_m && px.Q_per_m && px.U_per_m &&
                                   wlb_I && wlb_Q && wlb_U &&
                                   rt_aw_coupled_field_alloc(&cf, Mc_p - 1, w_res.n_mu_water_filled) == 0);
                        /* water no-injection baseline (same everything, ext_top NULL) */
                        const int nw3 = w_res.n_mu_water_filled;
                        const int McW = w_res.m_max_filled + 1;
                        const size_t plw = (size_t)McW * (size_t)nw3;
                        rt_water_rt_result_t p0; memset(&p0, 0, sizeof p0);
                        p0.I_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                        p0.Q_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                        p0.U_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                        p0.mu_water_pos = (double*)calloc((size_t)(nw3+4), sizeof(double));
                        p0.w_water_pos  = (double*)calloc((size_t)(nw3+4), sizeof(double));
                        okP = okP && p0.I_up_per_m && p0.Q_up_per_m && p0.U_up_per_m &&
                              p0.mu_water_pos && p0.w_water_pos;
                        const size_t plw_pp = (size_t)Mc_p * (size_t)nw3;
                        double *d1_I = (double*)calloc(plw_pp, sizeof(double));   /* D3 Ed gate: bounce-1 in-water downwelling */
                        double *d1_Q = (double*)calloc(plw_pp, sizeof(double));
                        double *d1_U = (double*)calloc(plw_pp, sizeof(double));
                        double dEd1_actual = 0.0;
                        if (okP) {
                            fprintf(stderr, "[D3pp] p0 solve begin (nw3=%d McW=%d)\n", nw3, McW); fflush(stderr);
                            setenv("OCRT_EXTTOP_KERNEL", "table", 1);
                            rt_water_rt_options_t wo0 = w_opts;
                            { const char *te = getenv("OCRT_D3_TOL"); if (te && atof(te) > 0.0) wo0.tolerance = atof(te); }
                            wo0.bypass_grid_cache = 1;
                            wo0.m_max_water = Mc_p - 1;   /* injection carries only Mc_p modes; mode
                                                           * decoupling (proven by the linearity gates)
                                                           * makes higher water modes irrelevant here */
                            wo0.ext_top_I = NULL; wo0.ext_top_Q = NULL; wo0.ext_top_U = NULL;
                            wo0.ext_top_mu = NULL; wo0.ext_top_n = 0; wo0.ext_top_m_max = -1;
                            okP = (rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                        cs->wavelength_nm, T_C, S_gkg, n_w, F_sun_water,
                                        aw_lut, psi_T_lut, &wo0, &p0) == 0 &&
                                   p0.n_mu_water_filled == nw3);
                        }
                        double dTOA_I[4] = {0,0,0,0};
                        int bdone = 0;
                        for (int b = 1; b <= PPB && okP; ++b) {
                            /* (1) atm backscatter of current wl */
                            rt_options_t ob = opts_p2;
                            rt_options_set_ext_bottom_source(&ob, wlb_I, wlb_Q, wlb_U,
                                                             n_mu_atm, mmax);
                            rt_result_t rba = {0};
                            memset(px.I_per_m,0,plx*sizeof(double)); memset(px.Q_per_m,0,plx*sizeof(double)); memset(px.U_per_m,0,plx*sizeof(double));
                            if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rba, &px) != 0) { okP = 0; break; }
                            /* (2) transmit down into water (T_aw Mueller) */
                            if (rt_air_water_couple_atm_to_water(&px, n_w, cs->q_convention,
                                    p0.mu_water_pos, nw3, w_atm,
                                    cs->wind_speed, cs->sigma_type, &cf) != 0) { okP = 0; break; }
                            /* (3) water response to that downwelling only */
                            rt_water_rt_result_t pr; memset(&pr, 0, sizeof pr);
                            pr.I_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                            pr.Q_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                            pr.U_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                            pr.mu_water_pos = (double*)calloc((size_t)(nw3+4), sizeof(double));
                            pr.w_water_pos  = (double*)calloc((size_t)(nw3+4), sizeof(double));
                            rt_water_rt_options_t wob = w_opts;
                            { const char *te = getenv("OCRT_D3_TOL"); if (te && atof(te) > 0.0) wob.tolerance = atof(te); }
                            wob.bypass_grid_cache = 1;
                                wob.m_max_water = Mc_p - 1;
                                wob.ext_top_I = cf.I_inwater_per_m;
                            wob.ext_top_Q = cf.Q_inwater_per_m;
                            wob.ext_top_U = cf.U_inwater_per_m;
                            wob.ext_top_mu = p0.mu_water_pos;
                            wob.ext_top_n = nw3; wob.ext_top_m_max = McW - 1;
                            int rcw = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                          cs->wavelength_nm, T_C, S_gkg, n_w, F_sun_water,
                                          aw_lut, psi_T_lut, &wob, &pr);
                            if (rcw != 0 || pr.n_mu_water_filled != nw3) { okP = 0; }
                            if (okP) {
                                for (size_t t2 = 0; t2 < plw; ++t2) {
                                    pr.I_up_per_m[t2] -= p0.I_up_per_m[t2];
                                    pr.Q_up_per_m[t2] -= p0.Q_up_per_m[t2];
                                    pr.U_up_per_m[t2] -= p0.U_up_per_m[t2];
                                }
                                if (b == 1) {   /* D3 Ed gate captures */
                                    dEd1_actual = pr.Ed_0minus_water - p0.Ed_0minus_water;
                                    if (d1_I && cf.I_inwater_per_m) memcpy(d1_I, cf.I_inwater_per_m, plw_pp*sizeof(double));
                                    if (d1_Q && cf.Q_inwater_per_m) memcpy(d1_Q, cf.Q_inwater_per_m, plw_pp*sizeof(double));
                                    if (d1_U && cf.U_inwater_per_m) memcpy(d1_U, cf.U_inwater_per_m, plw_pp*sizeof(double));
                                }
                                /* (4) back up through the interface to 0+ */
                                if (rt_air_water_couple_water_to_atm(
                                        pr.I_up_per_m, pr.Q_up_per_m, pr.U_up_per_m,
                                        p0.mu_water_pos, p0.w_water_pos, nw3, mmax,
                                        n_w, cs->q_convention, cs->wind_speed, cs->sigma_type,
                                        mu_atm, n_mu_atm, w_atm,
                                        wlb_I, wlb_Q, wlb_U) != 0) okP = 0;
                            }
                            free(pr.I_up_per_m); free(pr.Q_up_per_m); free(pr.U_up_per_m);
                            free(pr.mu_water_pos); free(pr.w_water_pos);
                            if (!okP) break;
                            if (b == 1 && b1_I && b1_Q && b1_U) {
                                memcpy(b1_I, wlb_I, pla*sizeof(double));
                                memcpy(b1_Q, wlb_Q, pla*sizeof(double));
                                memcpy(b1_U, wlb_U, pla*sizeof(double));
                            }
                            /* (5) TOA increment of this bounce */
                            rt_options_t ob2 = opts_p2;
                            rt_options_set_ext_bottom_source(&ob2, wlb_I, wlb_Q, wlb_U,
                                                             n_mu_atm, mmax);
                            rt_result_t rb2 = {0};
                            if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob2, aer, &rb2, NULL) != 0) { okP = 0; break; }
                            dTOA_I[b] = rb2.I_TOA;
                            bdone = b;
                        }
                        /* v1.10 D3-1c (2026-07-10): OPERATOR CLOSED FORM -
                         * INSTRUMENTATION ONLY (env OCRT_D3_OP=1, requires
                         * OCRT_D3_PP=1 scope).  Assemble the whole-chain
                         * bounce operator W column-by-column (unit wl bases
                         * through the SAME four production stages as the
                         * ping-pong above), then solve the per-mode linear
                         * system (I - W_m) s_m = (W wl0)_m by Gaussian
                         * elimination with partial pivoting, and evaluate
                         * dTOA_op = atm2(ext_bottom = s).I_TOA.  Gates:
                         * (g1) field: W*wl0 vs a direct one-stage run;
                         * (g2) sum:  dTOA_op vs the ping-pong geometric sum. */
                        if (okP && getenv("OCRT_D3_OP")) {
                            const size_t colsz = (size_t)3 * pla;
                            const int NB = 3 * n_mu_atm;
                            double *Wop = (double*)calloc((size_t)NB * colsz, sizeof(double));
                            double *bx_I = (double*)calloc(pla, sizeof(double));
                            double *bx_Q = (double*)calloc(pla, sizeof(double));
                            double *bx_U = (double*)calloc(pla, sizeof(double));
                            int okO = (Wop && bx_I && bx_Q && bx_U);
                            int nsw = 0;
                            for (int cb = 0; cb < 3 && okO; ++cb)
                              for (int jb = 0; jb < n_mu_atm && okO; ++jb) {
                                memset(bx_I,0,pla*sizeof(double)); memset(bx_Q,0,pla*sizeof(double)); memset(bx_U,0,pla*sizeof(double));
                                double *tg = (cb==0)?bx_I:(cb==1)?bx_Q:bx_U;
                                for (int mm = 0; mm < Mc_p; ++mm) tg[(size_t)mm*n_mu_atm + jb] = 1.0;
                                /* stage (1): atm backscatter of the basis */
                                rt_options_t ob = opts_p2;
                                rt_options_set_ext_bottom_source(&ob, bx_I, bx_Q, bx_U,
                                                                 n_mu_atm, mmax);
                                rt_result_t rba = {0};
                                memset(px.I_per_m,0,plx*sizeof(double)); memset(px.Q_per_m,0,plx*sizeof(double)); memset(px.U_per_m,0,plx*sizeof(double));
                                if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rba, &px) != 0) { okO = 0; break; }
                                /* stage (2): down through the interface */
                                if (rt_air_water_couple_atm_to_water(&px, n_w, cs->q_convention,
                                        p0.mu_water_pos, nw3, w_atm,
                                        cs->wind_speed, cs->sigma_type, &cf) != 0) { okO = 0; break; }
                                /* stage (3): water response (minus baseline) */
                                rt_water_rt_result_t pr; memset(&pr, 0, sizeof pr);
                                pr.I_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                                pr.Q_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                                pr.U_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                                pr.mu_water_pos = (double*)calloc((size_t)(nw3+4), sizeof(double));
                                pr.w_water_pos  = (double*)calloc((size_t)(nw3+4), sizeof(double));
                                rt_water_rt_options_t wob = w_opts;
                                { const char *te = getenv("OCRT_D3_TOL"); if (te && atof(te) > 0.0) wob.tolerance = atof(te); }
                                wob.bypass_grid_cache = 1;
                                wob.m_max_water = Mc_p - 1;
                                wob.ext_top_I = cf.I_inwater_per_m;
                                wob.ext_top_Q = cf.Q_inwater_per_m;
                                wob.ext_top_U = cf.U_inwater_per_m;
                                wob.ext_top_mu = p0.mu_water_pos;
                                wob.ext_top_n = nw3; wob.ext_top_m_max = McW - 1;
                                int rcw = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                              cs->wavelength_nm, T_C, S_gkg, n_w, F_sun_water,
                                              aw_lut, psi_T_lut, &wob, &pr);
                                if (rcw == 0 && pr.n_mu_water_filled == nw3) {
                                    for (size_t t2 = 0; t2 < plw; ++t2) {
                                        pr.I_up_per_m[t2] -= p0.I_up_per_m[t2];
                                        pr.Q_up_per_m[t2] -= p0.Q_up_per_m[t2];
                                        pr.U_up_per_m[t2] -= p0.U_up_per_m[t2];
                                    }
                                    /* stage (4): back up to 0+ = W column */
                                    double *col = Wop + ((size_t)cb*n_mu_atm + jb) * colsz;
                                    if (rt_air_water_couple_water_to_atm(
                                            pr.I_up_per_m, pr.Q_up_per_m, pr.U_up_per_m,
                                            p0.mu_water_pos, p0.w_water_pos, nw3, mmax,
                                            n_w, cs->q_convention, cs->wind_speed, cs->sigma_type,
                                            mu_atm, n_mu_atm, w_atm,
                                            col + 0*pla, col + 1*pla, col + 2*pla) != 0) okO = 0;
                                    else ++nsw;
                                } else okO = 0;
                                free(pr.I_up_per_m); free(pr.Q_up_per_m); free(pr.U_up_per_m);
                                free(pr.mu_water_pos); free(pr.w_water_pos);
                              }
                            /* g1 field gate: W*wl0 vs pingpong wl1 (b1_*) */
                            double g1 = -1.0;
                            double *s_I = (double*)calloc(pla, sizeof(double));
                            double *s_Q = (double*)calloc(pla, sizeof(double));
                            double *s_U = (double*)calloc(pla, sizeof(double));
                            if (okO && s_I && s_Q && s_U) {
                                double mx = 0.0, rf = 0.0;
                                for (int mm = 0; mm < Mc_p; ++mm)
                                  for (int k = 0; k < n_mu_atm; ++k) {
                                    size_t t2 = (size_t)mm*n_mu_atm + k;
                                    double yI=0, yQ=0, yU=0;
                                    for (int cb = 0; cb < 3; ++cb)
                                      for (int jb = 0; jb < n_mu_atm; ++jb) {
                                        double c0 = (cb==0?wl0_I:(cb==1?wl0_Q:wl0_U))[(size_t)mm*n_mu_atm + jb];
                                        if (c0 == 0.0) continue;
                                        const double *col = Wop + ((size_t)cb*n_mu_atm + jb) * colsz;
                                        yI += c0 * col[0*pla + t2];
                                        yQ += c0 * col[1*pla + t2];
                                        yU += c0 * col[2*pla + t2];
                                      }
                                    s_I[t2]=yI; s_Q[t2]=yQ; s_U[t2]=yU;   /* temp: W*wl0 */
                                    double eI=fabs(yI-b1_I[t2]), eQ=fabs(yQ-b1_Q[t2]), eU=fabs(yU-b1_U[t2]);
                                    double v=fabs(b1_I[t2])+fabs(b1_Q[t2])+fabs(b1_U[t2]);
                                    if (v>rf) rf=v;
                                    if (eI>mx) mx=eI; if (eQ>mx) mx=eQ; if (eU>mx) mx=eU;
                                  }
                                g1 = (rf>0.0)? mx/rf : mx;
                            }
                            /* closed form: per mode solve (I - W_m) s_m = (W wl0)_m */
                            double dTOA_op = 0.0; int lu_ok = okO;
                            if (okO) {
                                int N3 = 3 * n_mu_atm;
                                double *A = (double*)malloc((size_t)N3*N3*sizeof(double));
                                double *rhs = (double*)malloc((size_t)N3*sizeof(double));
                                int *piv = (int*)malloc((size_t)N3*sizeof(int));
                                for (int mm = 0; mm < Mc_p && lu_ok && A && rhs && piv; ++mm) {
                                    for (int cc = 0; cc < 3; ++cc)
                                      for (int k = 0; k < n_mu_atm; ++k) {
                                        int row = cc*n_mu_atm + k;
                                        size_t t2 = (size_t)mm*n_mu_atm + k;
                                        rhs[row] = (cc==0? s_I:(cc==1? s_Q:s_U))[t2];
                                        for (int cb = 0; cb < 3; ++cb)
                                          for (int jb = 0; jb < n_mu_atm; ++jb) {
                                            int colj = cb*n_mu_atm + jb;
                                            const double *col = Wop + ((size_t)cb*n_mu_atm + jb) * colsz;
                                            double wv = col[(size_t)cc*pla + t2];
                                            A[(size_t)row*N3 + colj] = (row==colj ? 1.0 : 0.0) - wv;
                                          }
                                      }
                                    /* LU w/ partial pivot, solve in place */
                                    for (int c2 = 0; c2 < N3 && lu_ok; ++c2) {
                                        int p2 = c2; double mxp = fabs(A[(size_t)c2*N3+c2]);
                                        for (int r2 = c2+1; r2 < N3; ++r2) { double v=fabs(A[(size_t)r2*N3+c2]); if (v>mxp){mxp=v;p2=r2;} }
                                        if (mxp < 1e-14) { lu_ok = 0; break; }
                                        if (p2 != c2) { for (int q2=0;q2<N3;++q2){double t3=A[(size_t)c2*N3+q2];A[(size_t)c2*N3+q2]=A[(size_t)p2*N3+q2];A[(size_t)p2*N3+q2]=t3;} double t3=rhs[c2];rhs[c2]=rhs[p2];rhs[p2]=t3; }
                                        for (int r2 = c2+1; r2 < N3; ++r2) {
                                            double f2 = A[(size_t)r2*N3+c2]/A[(size_t)c2*N3+c2];
                                            if (f2 == 0.0) continue;
                                            for (int q2=c2;q2<N3;++q2) A[(size_t)r2*N3+q2] -= f2*A[(size_t)c2*N3+q2];
                                            rhs[r2] -= f2*rhs[c2];
                                        }
                                    }
                                    for (int r2 = N3-1; r2 >= 0 && lu_ok; --r2) {
                                        double acc2 = rhs[r2];
                                        for (int q2=r2+1;q2<N3;++q2) acc2 -= A[(size_t)r2*N3+q2]*rhs[q2];
                                        rhs[r2] = acc2 / A[(size_t)r2*N3+r2];
                                    }
                                    for (int cc = 0; cc < 3 && lu_ok; ++cc)
                                      for (int k = 0; k < n_mu_atm; ++k) {
                                        size_t t2 = (size_t)mm*n_mu_atm + k;
                                        (cc==0? s_I:(cc==1? s_Q:s_U))[t2] = rhs[cc*n_mu_atm + k];
                                      }
                                }
                                free(A); free(rhs); free(piv);
                                if (lu_ok) {
                                    rt_options_t ob2 = opts_p2;
                                    rt_options_set_ext_bottom_source(&ob2, s_I, s_Q, s_U,
                                                                     n_mu_atm, mmax);
                                    rt_result_t rb3 = {0};
                                    if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob2, aer, &rb3, NULL) == 0)
                                        dTOA_op = rb3.I_TOA;
                                    else lu_ok = 0;
                                }
                            }
                            double pp_sum = dTOA_I[1] + dTOA_I[2] + dTOA_I[3];
                            fprintf(stderr, "[D3-1c] Wcols=%d g1(field W*wl0 vs pp b1)=%.3e dTOA_op=%.6e pp_sum=%.6e reldiff=%.3e lu_ok=%d\n",
                                    nsw, g1, dTOA_op, pp_sum,
                                    (pp_sum != 0.0) ? fabs(dTOA_op - pp_sum)/fabs(pp_sum) : -1.0, lu_ok);
                            WopKeep = Wop; Wop = NULL;   /* transfer to PP scope */
                            free(bx_I); free(bx_Q); free(bx_U);
                            free(s_I); free(s_Q); free(s_U);
                        }
                        /* v1.10 D3-3 (2026-07-10): PRODUCTION-SCALE FACTORIZED
                         * OPERATOR - INSTRUMENTATION ONLY (env OCRT_D3_FAC=1,
                         * independent of the OP whole-chain section).
                         * Assembles the four reusable pieces at THIS run's
                         * discretization, with WALL TIMERS on the two solver
                         * pieces (the amortized costs of the D3 production
                         * architecture):
                         *   RAo (3*n_a atm solves)  - per atm-config cost
                         *   RBw (3*nw water solves) - per (case,wl,wind) cost
                         *   Cdn/Cup                  - function-call matrices
                         * Then builds W_fac = Cup.RBw.Cdn.RAo, optionally
                         * gates it against the whole-chain W (if the OP
                         * section ran), solves the closed form from the
                         * PIECES ALONE, and gates dTOA_fac against the
                         * ping-pong geometric reference. */
                        if (okP && getenv("OCRT_D3_FAC")) {
                            const size_t plaw = (size_t)Mc_p * (size_t)nw3;
                            const size_t colsz = (size_t)3 * pla;
                            const int NBa = 3 * n_mu_atm, NBw = 3 * nw3;
                            double *RAo = (double*)calloc((size_t)NBa * 3 * pla,  sizeof(double));
                            double *Cdn = (double*)calloc((size_t)NBa * 3 * plaw, sizeof(double));
                            double *RBw = (double*)calloc((size_t)NBw * 3 * plaw, sizeof(double));
                            double *Cup = (double*)calloc((size_t)NBw * 3 * pla,  sizeof(double));
                            double *Wfac = (double*)calloc((size_t)NBa * colsz,   sizeof(double));
                            double *bx_I = (double*)calloc(pla, sizeof(double));
                            double *bx_Q = (double*)calloc(pla, sizeof(double));
                            double *bx_U = (double*)calloc(pla, sizeof(double));
                            double *EdRBf = (double*)calloc((size_t)3*nw3, sizeof(double));   /* Ed response row */
                            double *VwRBf = (double*)calloc((size_t)3*nw3*3, sizeof(double)); /* view-scalar response rows (I/Q/U out) */
                            double *wt_I = (double*)calloc(plaw, sizeof(double));
                            double *wt_Q = (double*)calloc(plaw, sizeof(double));
                            double *wt_U = (double*)calloc(plaw, sizeof(double));
                            int okF = (RAo && Cdn && RBw && Cup && Wfac && bx_I && bx_Q && bx_U && wt_I && wt_Q && wt_U);
                            double tA = 0.0, tW = 0.0, t0f;
                            const char *facmode = getenv("OCRT_D3_FAC");
                            const int fac_rbonly = (facmode && !strcmp(facmode, "rbonly"));
                            /* (i) RAo (skipped in rbonly file-export mode) */
                            t0f = ocrt_s6_now();
                            for (int cb = 0; cb < 3 && okF && !fac_rbonly; ++cb)
                              for (int jb = 0; jb < n_mu_atm && okF; ++jb) {
                                memset(bx_I,0,pla*sizeof(double)); memset(bx_Q,0,pla*sizeof(double)); memset(bx_U,0,pla*sizeof(double));
                                double *tg = (cb==0)?bx_I:(cb==1)?bx_Q:bx_U;
                                for (int mm = 0; mm < Mc_p; ++mm) tg[(size_t)mm*n_mu_atm + jb] = 1.0;
                                rt_options_t ob = opts_p2;
                                rt_options_set_ext_bottom_source(&ob, bx_I, bx_Q, bx_U,
                                                                 n_mu_atm, mmax);
                                rt_result_t rr = {0};
                                memset(px.I_per_m,0,plx*sizeof(double)); memset(px.Q_per_m,0,plx*sizeof(double)); memset(px.U_per_m,0,plx*sizeof(double));
                                if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob, aer, &rr, &px) != 0) { okF = 0; break; }
                                double *col = RAo + ((size_t)cb*n_mu_atm + jb) * 3 * pla;
                                memcpy(col + 0*pla, px.I_per_m, pla*sizeof(double));
                                memcpy(col + 1*pla, px.Q_per_m, pla*sizeof(double));
                                memcpy(col + 2*pla, px.U_per_m, pla*sizeof(double));
                              }
                            tA = ocrt_s6_now() - t0f;
                            /* (ii) Cdn (skipped in rbonly) */
                            for (int cb = 0; cb < 3 && okF && !fac_rbonly; ++cb)
                              for (int jb = 0; jb < n_mu_atm && okF; ++jb) {
                                memset(px.I_per_m,0,plx*sizeof(double)); memset(px.Q_per_m,0,plx*sizeof(double)); memset(px.U_per_m,0,plx*sizeof(double));
                                double *tg = (cb==0)?px.I_per_m:(cb==1)?px.Q_per_m:px.U_per_m;
                                for (int mm = 0; mm < Mc_p; ++mm) tg[(size_t)mm*n_mu_atm + jb] = 1.0;
                                if (rt_air_water_couple_atm_to_water(&px, n_w, cs->q_convention,
                                        p0.mu_water_pos, nw3, w_atm,
                                        cs->wind_speed, cs->sigma_type, &cf) != 0) { okF = 0; break; }
                                double *col = Cdn + ((size_t)cb*n_mu_atm + jb) * 3 * plaw;
                                memcpy(col + 0*plaw, cf.I_inwater_per_m, plaw*sizeof(double));
                                memcpy(col + 1*plaw, cf.Q_inwater_per_m, plaw*sizeof(double));
                                memcpy(col + 2*plaw, cf.U_inwater_per_m, plaw*sizeof(double));
                              }
                            /* (iii) RBw - optional column partitioning for
                             * sandbox time limits: OCRT_D3_RB_PART="k/N"
                             * computes only global columns with idx%N==k;
                             * files from parts are merged offline (columns
                             * not computed stay zero and are filled by the
                             * other parts). */
                            int rb_pk = 0, rb_pn = 1;
                            { const char *pe = getenv("OCRT_D3_RB_PART");
                              if (pe && sscanf(pe, "%d/%d", &rb_pk, &rb_pn) == 2 && rb_pn > 0) {}
                              else { rb_pk = 0; rb_pn = 1; } }
                            t0f = ocrt_s6_now();
                            fprintf(stderr, "[D3rb] begin nw3=%d Mc=%d part=%d/%d\n", nw3, Mc_p, rb_pk, rb_pn); fflush(stderr);
                            for (int cb = 0; cb < 3 && okF; ++cb)
                              for (int jb = 0; jb < nw3 && okF; ++jb) {
                                if (((cb*nw3 + jb) % rb_pn) != rb_pk) continue;
                                /* zero-quadrature-weight nodes (view/slot
                                 * insertions) have identically zero injection
                                 * (AI = 2 g_b ext = 0) - measured: only the
                                 * 24 GL columns of 153 are nonzero.  Skip the
                                 * solve; the calloc'd response stays 0. */
                                if (p0.w_water_pos && p0.w_water_pos[jb] == 0.0) continue;
                                { static int d3rb_ct = 0;
                                  if ((d3rb_ct++ % 10) == 0) {
                                      fprintf(stderr, "[D3rb] col %d (cb=%d jb=%d) t=%.1fs\n",
                                              d3rb_ct, cb, jb, ocrt_s6_now()-t0f);
                                      fflush(stderr);
                                  } }
                                memset(wt_I,0,plaw*sizeof(double)); memset(wt_Q,0,plaw*sizeof(double)); memset(wt_U,0,plaw*sizeof(double));
                                double *tg = (cb==0)?wt_I:(cb==1)?wt_Q:wt_U;
                                for (int mm = 0; mm < Mc_p; ++mm) tg[(size_t)mm*nw3 + jb] = 1.0;
                                rt_water_rt_result_t pr; memset(&pr, 0, sizeof pr);
                                pr.I_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                                pr.Q_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                                pr.U_up_per_m = (double*)calloc((size_t)(McW)*(nw3+4), sizeof(double));
                                pr.mu_water_pos = (double*)calloc((size_t)(nw3+4), sizeof(double));
                                pr.w_water_pos  = (double*)calloc((size_t)(nw3+4), sizeof(double));
                                rt_water_rt_options_t wob = w_opts;
                                { const char *te = getenv("OCRT_D3_TOL"); if (te && atof(te) > 0.0) wob.tolerance = atof(te); }
                                wob.bypass_grid_cache = 1;
                                wob.m_max_water = Mc_p - 1;
                                    wob.ext_top_I = wt_I; wob.ext_top_Q = wt_Q; wob.ext_top_U = wt_U;
                                wob.ext_top_mu = p0.mu_water_pos;
                                wob.ext_top_n = nw3; wob.ext_top_m_max = Mc_p - 1;
                                int rcw = rt_water_rt_sos_pure(cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                              cs->wavelength_nm, T_C, S_gkg, n_w, F_sun_water,
                                              aw_lut, psi_T_lut, &wob, &pr);
                                if (rcw == 0 && pr.n_mu_water_filled == nw3 && EdRBf)
                                    EdRBf[(size_t)cb*nw3 + jb] = pr.Ed_0minus_water - p0.Ed_0minus_water;
                                if (rcw == 0 && pr.n_mu_water_filled == nw3 && VwRBf) {
                                    VwRBf[((size_t)cb*nw3 + jb)*3 + 0] = pr.I_0minus_view - p0.I_0minus_view;
                                    VwRBf[((size_t)cb*nw3 + jb)*3 + 1] = pr.Q_0minus_view - p0.Q_0minus_view;
                                    VwRBf[((size_t)cb*nw3 + jb)*3 + 2] = pr.U_0minus_view - p0.U_0minus_view;
                                }
                                if (rcw == 0 && pr.n_mu_water_filled == nw3) {
                                    double *col = RBw + ((size_t)cb*nw3 + jb) * 3 * plaw;
                                    for (size_t tw = 0; tw < plaw; ++tw) {
                                        col[0*plaw+tw] = pr.I_up_per_m[tw] - p0.I_up_per_m[tw];
                                        col[1*plaw+tw] = pr.Q_up_per_m[tw] - p0.Q_up_per_m[tw];
                                        col[2*plaw+tw] = pr.U_up_per_m[tw] - p0.U_up_per_m[tw];
                                    }
                                } else okF = 0;
                                free(pr.I_up_per_m); free(pr.Q_up_per_m); free(pr.U_up_per_m);
                                free(pr.mu_water_pos); free(pr.w_water_pos);
                              }
                            tW = ocrt_s6_now() - t0f;
                            /* (iv) Cup (skipped in rbonly) */
                            for (int cb = 0; cb < 3 && okF && !fac_rbonly; ++cb)
                              for (int jb = 0; jb < nw3 && okF; ++jb) {
                                memset(wt_I,0,plaw*sizeof(double)); memset(wt_Q,0,plaw*sizeof(double)); memset(wt_U,0,plaw*sizeof(double));
                                double *tg = (cb==0)?wt_I:(cb==1)?wt_Q:wt_U;
                                for (int mm = 0; mm < Mc_p; ++mm) tg[(size_t)mm*nw3 + jb] = 1.0;
                                double *col = Cup + ((size_t)cb*nw3 + jb) * 3 * pla;
                                if (rt_air_water_couple_water_to_atm(
                                        wt_I, wt_Q, wt_U,
                                        p0.mu_water_pos, p0.w_water_pos, nw3, mmax,
                                        n_w, cs->q_convention, cs->wind_speed, cs->sigma_type,
                                        mu_atm, n_mu_atm, w_atm,
                                        col + 0*pla, col + 1*pla, col + 2*pla) != 0) { okF = 0; break; }
                              }
                            /* W_fac columns = Cup.RBw.Cdn.RAo */
                            if (okF && !fac_rbonly) {
                                double *v_dn = (double*)malloc(3*plaw*sizeof(double));
                                double *v_up = (double*)malloc(3*plaw*sizeof(double));
                                for (int cb = 0; cb < 3 && v_dn && v_up; ++cb)
                                  for (int jb = 0; jb < n_mu_atm; ++jb) {
                                    const double *v_exp = RAo + ((size_t)cb*n_mu_atm + jb)*3*pla;
                                    memset(v_dn, 0, 3*plaw*sizeof(double));
                                    for (int ci = 0; ci < 3; ++ci)
                                      for (int ja = 0; ja < n_mu_atm; ++ja) {
                                        const double *cc2 = Cdn + ((size_t)ci*n_mu_atm + ja)*3*plaw;
                                        for (int mm = 0; mm < Mc_p; ++mm) {
                                          double c0 = v_exp[(size_t)ci*pla + (size_t)mm*n_mu_atm + ja];
                                          if (c0 == 0.0) continue;
                                          for (int co = 0; co < 3; ++co)
                                            for (int kw = 0; kw < nw3; ++kw)
                                              v_dn[(size_t)co*plaw + (size_t)mm*nw3 + kw] += c0 * cc2[(size_t)co*plaw + (size_t)mm*nw3 + kw];
                                        }
                                      }
                                    memset(v_up, 0, 3*plaw*sizeof(double));
                                    for (int ci = 0; ci < 3; ++ci)
                                      for (int jw = 0; jw < nw3; ++jw) {
                                        const double *cc2 = RBw + ((size_t)ci*nw3 + jw)*3*plaw;
                                        for (int mm = 0; mm < Mc_p; ++mm) {
                                          double c0 = v_dn[(size_t)ci*plaw + (size_t)mm*nw3 + jw];
                                          if (c0 == 0.0) continue;
                                          for (int co = 0; co < 3; ++co)
                                            for (int kw = 0; kw < nw3; ++kw)
                                              v_up[(size_t)co*plaw + (size_t)mm*nw3 + kw] += c0 * cc2[(size_t)co*plaw + (size_t)mm*nw3 + kw];
                                        }
                                      }
                                    double *wcol = Wfac + ((size_t)cb*n_mu_atm + jb) * colsz;
                                    for (int ci = 0; ci < 3; ++ci)
                                      for (int jw = 0; jw < nw3; ++jw) {
                                        const double *cc2 = Cup + ((size_t)ci*nw3 + jw)*3*pla;
                                        for (int mm = 0; mm < Mc_p; ++mm) {
                                          double c0 = v_up[(size_t)ci*plaw + (size_t)mm*nw3 + jw];
                                          if (c0 == 0.0) continue;
                                          for (int co = 0; co < 3; ++co)
                                            for (int ka = 0; ka < n_mu_atm; ++ka)
                                              wcol[(size_t)co*pla + (size_t)mm*n_mu_atm + ka] += c0 * cc2[(size_t)co*pla + (size_t)mm*n_mu_atm + ka];
                                        }
                                      }
                                  }
                                free(v_dn); free(v_up);
                            }
                            /* optional gate vs whole-chain W */
                            double gF = -1.0;
                            if (okF && WopKeep) {
                                double mx = 0.0, rf = 0.0;
                                for (size_t t2 = 0; t2 < (size_t)NBa*colsz; ++t2) {
                                    double e = fabs(Wfac[t2] - WopKeep[t2]);
                                    double v = fabs(WopKeep[t2]);
                                    if (v > rf) rf = v;
                                    if (e > mx) mx = e;
                                }
                                gF = (rf > 0.0) ? mx/rf : mx;
                            }
                            /* closed form FROM PIECES: (I - Wfac_m) s_m = (Wfac wl0)_m */
                            double dTOA_fac = 0.0; int lu2 = okF && !fac_rbonly;
                            double *sf_I = (double*)calloc(pla, sizeof(double));
                            double *sf_Q = (double*)calloc(pla, sizeof(double));
                            double *sf_U = (double*)calloc(pla, sizeof(double));
                            if (okF && sf_I && sf_Q && sf_U) {
                                for (int mm = 0; mm < Mc_p; ++mm)
                                  for (int k = 0; k < n_mu_atm; ++k) {
                                    size_t t2 = (size_t)mm*n_mu_atm + k;
                                    double yI=0, yQ=0, yU=0;
                                    for (int cb = 0; cb < 3; ++cb)
                                      for (int jb = 0; jb < n_mu_atm; ++jb) {
                                        double c0 = (cb==0?wl0_I:(cb==1?wl0_Q:wl0_U))[(size_t)mm*n_mu_atm + jb];
                                        if (c0 == 0.0) continue;
                                        const double *col = Wfac + ((size_t)cb*n_mu_atm + jb) * colsz;
                                        yI += c0 * col[0*pla + t2];
                                        yQ += c0 * col[1*pla + t2];
                                        yU += c0 * col[2*pla + t2];
                                      }
                                    sf_I[t2]=yI; sf_Q[t2]=yQ; sf_U[t2]=yU;
                                  }
                                int N3 = 3 * n_mu_atm;
                                double *A = (double*)malloc((size_t)N3*N3*sizeof(double));
                                double *rhs = (double*)malloc((size_t)N3*sizeof(double));
                                for (int mm = 0; mm < Mc_p && lu2 && A && rhs; ++mm) {
                                    for (int cc = 0; cc < 3; ++cc)
                                      for (int k = 0; k < n_mu_atm; ++k) {
                                        int row = cc*n_mu_atm + k;
                                        size_t t2 = (size_t)mm*n_mu_atm + k;
                                        rhs[row] = (cc==0? sf_I:(cc==1? sf_Q:sf_U))[t2];
                                        for (int cb = 0; cb < 3; ++cb)
                                          for (int jb = 0; jb < n_mu_atm; ++jb) {
                                            int colj = cb*n_mu_atm + jb;
                                            const double *col = Wfac + ((size_t)cb*n_mu_atm + jb) * colsz;
                                            A[(size_t)row*N3 + colj] = (row==colj ? 1.0 : 0.0) - col[(size_t)cc*pla + t2];
                                          }
                                      }
                                    for (int c2 = 0; c2 < N3 && lu2; ++c2) {
                                        int p2 = c2; double mxp = fabs(A[(size_t)c2*N3+c2]);
                                        for (int r2 = c2+1; r2 < N3; ++r2) { double v=fabs(A[(size_t)r2*N3+c2]); if (v>mxp){mxp=v;p2=r2;} }
                                        if (mxp < 1e-14) { lu2 = 0; break; }
                                        if (p2 != c2) { for (int q2=0;q2<N3;++q2){double t3=A[(size_t)c2*N3+q2];A[(size_t)c2*N3+q2]=A[(size_t)p2*N3+q2];A[(size_t)p2*N3+q2]=t3;} double t3=rhs[c2];rhs[c2]=rhs[p2];rhs[p2]=t3; }
                                        for (int r2 = c2+1; r2 < N3; ++r2) {
                                            double f2 = A[(size_t)r2*N3+c2]/A[(size_t)c2*N3+c2];
                                            if (f2 == 0.0) continue;
                                            for (int q2=c2;q2<N3;++q2) A[(size_t)r2*N3+q2] -= f2*A[(size_t)c2*N3+q2];
                                            rhs[r2] -= f2*rhs[c2];
                                        }
                                    }
                                    for (int r2 = N3-1; r2 >= 0 && lu2; --r2) {
                                        double acc2 = rhs[r2];
                                        for (int q2=r2+1;q2<N3;++q2) acc2 -= A[(size_t)r2*N3+q2]*rhs[q2];
                                        rhs[r2] = acc2 / A[(size_t)r2*N3+r2];
                                    }
                                    for (int cc = 0; cc < 3 && lu2; ++cc)
                                      for (int k = 0; k < n_mu_atm; ++k)
                                        (cc==0? sf_I:(cc==1? sf_Q:sf_U))[(size_t)mm*n_mu_atm + k] = rhs[cc*n_mu_atm + k];
                                }
                                free(A); free(rhs);
                                if (lu2) {
                                    rt_options_t ob2 = opts_p2;
                                    rt_options_set_ext_bottom_source(&ob2, sf_I, sf_Q, sf_U,
                                                                     n_mu_atm, mmax);
                                    rt_result_t rb3 = {0};
                                    if (rt_solve_case_pol_for_ocean(&cs_atm2, &ob2, aer, &rb3, NULL) == 0)
                                        dTOA_fac = rb3.I_TOA;
                                    else lu2 = 0;
                                }
                            }
                            /* v1.10 D3-PROD (2026-07-10): formal operator file.
                             * env OCRT_D3_RB_OUT=path -> write the assembled
                             * water operator for consumption by the opt-in
                             * production path: header {magic 'D3RB', ver, nw3,
                             * Mc_p}, mu[nw3], w[nw3], EdRB[3*nw3],
                             * RB[(3*nw3) cols x 3 comps x (Mc_p*nw3)]. */
                            { const char *rbout = getenv("OCRT_D3_RB_OUT");
                              if (okF && rbout && rbout[0]) {
                                FILE *fo = fopen(rbout, "wb");
                                if (fo) {
                                    int hdr[4] = { 0x44335242, 2, nw3, Mc_p };
                                    fwrite(hdr, sizeof(int), 4, fo);
                                    fwrite(p0.mu_water_pos, sizeof(double), (size_t)nw3, fo);
                                    fwrite(p0.w_water_pos,  sizeof(double), (size_t)nw3, fo);
                                    fwrite(EdRBf, sizeof(double), (size_t)3*nw3, fo);
                                    fwrite(VwRBf, sizeof(double), (size_t)3*nw3*3, fo);
                                    fwrite(RBw, sizeof(double), (size_t)NBw*3*plaw, fo);
                                    fclose(fo);
                                    fprintf(stderr, "[D3-PROD] operator written: %s (n=%d Mc=%d)\n", rbout, nw3, Mc_p);
                                }
                              }
                            }
                            double gEd = -1.0;
                            if (okF && !fac_rbonly && EdRBf && d1_I) {
                                double pred = 0.0;
                                for (int cb = 0; cb < 3; ++cb)
                                  for (int jb = 0; jb < nw3; ++jb)
                                    pred += EdRBf[(size_t)cb*nw3 + jb] *
                                            ((cb==0?d1_I:(cb==1?d1_Q:d1_U))[(size_t)0*nw3 + jb]);   /* m=0 plane only */
                                double den = fabs(dEd1_actual) > 1e-300 ? fabs(dEd1_actual) : 1.0;
                                gEd = fabs(pred - dEd1_actual) / den;
                            }
                            double pp_sum2 = dTOA_I[1] + dTOA_I[2] + dTOA_I[3];
                            fprintf(stderr, "[D3-3] prod-scale pieces: RAo %d cols t=%.1fs, RBw %d cols t=%.1fs; gF(vs whole-chain)=%.3e gEd=%.3e dTOA_fac=%.6e pp_sum=%.6e reldiff=%.3e lu=%d ok=%d\n",
                                    NBa, tA, NBw, tW, gF, gEd, dTOA_fac, pp_sum2,
                                    (pp_sum2 != 0.0)? fabs(dTOA_fac - pp_sum2)/fabs(pp_sum2) : -1.0, lu2, okF);
                            free(RAo); free(Cdn); free(RBw); free(Cup); free(Wfac);
                            free(bx_I); free(bx_Q); free(bx_U);
                            free(wt_I); free(wt_Q); free(wt_U);
                            free(sf_I); free(sf_Q); free(sf_U); free(EdRBf); free(VwRBf);
                        }

                        unsetenv("OCRT_EXTTOP_KERNEL");
                        free(WopKeep);
                        free(d1_I); free(d1_Q); free(d1_U);
                        double ref = fabs(TOA_wl_direct_I) > 0.0 ? fabs(TOA_wl_direct_I) : 1.0;
                        fprintf(stderr, "[D3-1b] pingpong bounces=%d ref(TOA_wl)=%.6e dTOA1=%.3e (%.3e rel) dTOA2=%.3e dTOA3=%.3e ratio21=%.3e ok=%d\n",
                                bdone, TOA_wl_direct_I,
                                dTOA_I[1], dTOA_I[1]/ref, dTOA_I[2], dTOA_I[3],
                                (fabs(dTOA_I[1])>0.0)? dTOA_I[2]/dTOA_I[1] : 0.0, okP);
                        rt_aw_coupled_field_free(&cf);
                        free(px.mu_quad_pos); free(px.I_per_m); free(px.Q_per_m); free(px.U_per_m);
                        free(wlb_I); free(wlb_Q); free(wlb_U);
                        free(wl0_I); free(wl0_Q); free(wl0_U);
                        free(b1_I); free(b1_Q); free(b1_U);
                        free(p0.I_up_per_m); free(p0.Q_up_per_m); free(p0.U_up_per_m);
                        free(p0.mu_water_pos); free(p0.w_water_pos);
                    }
                    if (d3_gc_save) setenv("OCRT_WATER_GRID_CACHE", d3_gc_save, 1);
                    else unsetenv("OCRT_WATER_GRID_CACHE");
                  }
                }
        }
        free(mu_atm); free(w_atm); free(wlI); free(wlQ); free(wlU);
    }
    if (atm_present && g_ocean_water_capture && g_ocean_water_capture->active &&
        w_res.I_up_per_m && w_res.Q_up_per_m && w_res.U_up_per_m &&
        w_res.mu_water_pos && w_res.n_mu_water_filled > 0 && w_res.m_max_filled >= 0) {
        ocrt_ocean_water_capture_t *cap = g_ocean_water_capture;
        const int ncap = w_res.n_mu_water_filled;
        const int mcap = w_res.m_max_filled;
        const size_t nf = (size_t)ncap * (size_t)(mcap + 1);
        free(cap->mu); free(cap->w); free(cap->I); free(cap->Q); free(cap->U);
        cap->mu = (double*)malloc((size_t)ncap * sizeof(double));
        cap->w  = (double*)malloc((size_t)ncap * sizeof(double));
        cap->I  = (double*)malloc(nf * sizeof(double));
        cap->Q  = (double*)malloc(nf * sizeof(double));
        cap->U  = (double*)malloc(nf * sizeof(double));
        if (cap->mu && cap->w && cap->I && cap->Q && cap->U) {
            memcpy(cap->mu, w_res.mu_water_pos, (size_t)ncap * sizeof(double));
            memcpy(cap->w,  w_res.w_water_pos,  (size_t)ncap * sizeof(double));
            memcpy(cap->I,  w_res.I_up_per_m, nf * sizeof(double));
            memcpy(cap->Q,  w_res.Q_up_per_m, nf * sizeof(double));
            memcpy(cap->U,  w_res.U_up_per_m, nf * sizeof(double));
            cap->n_mu = ncap; cap->m_max = mcap; cap->valid = 1;
            free(cap->vI); free(cap->vQ); free(cap->vU);
            free(cap->aI); free(cap->aQ); free(cap->aU);
            cap->vI=cap->vQ=cap->vU=cap->aI=cap->aQ=cap->aU=NULL;
            if (w_res.view_m_max_filled >= 0 && w_res.view_modes_exact) {
                const size_t nm=(size_t)(w_res.view_m_max_filled+1);
                cap->vI=(double*)malloc(nm*sizeof(double)); cap->vQ=(double*)malloc(nm*sizeof(double)); cap->vU=(double*)malloc(nm*sizeof(double));
                cap->aI=(double*)malloc(nm*sizeof(double)); cap->aQ=(double*)malloc(nm*sizeof(double)); cap->aU=(double*)malloc(nm*sizeof(double));
                if (cap->vI&&cap->vQ&&cap->vU&&cap->aI&&cap->aQ&&cap->aU) {
                    memcpy(cap->vI,w_res.I_view_per_m,nm*sizeof(double)); memcpy(cap->vQ,w_res.Q_view_per_m,nm*sizeof(double)); memcpy(cap->vU,w_res.U_view_per_m,nm*sizeof(double));
                    memcpy(cap->aI,w_res.I_air_view_per_m,nm*sizeof(double)); memcpy(cap->aQ,w_res.Q_air_view_per_m,nm*sizeof(double)); memcpy(cap->aU,w_res.U_air_view_per_m,nm*sizeof(double));
                    cap->view_m_max=w_res.view_m_max_filled; cap->view_exact=1;
                }
            }
        } else cap->valid = 0;
    }

    if (atm_present) {
        if (!d3p_on && ocrt_dt_I) {
            Kd_0minus_total = ocrt_kd_with_unscattered_sky(
                &w_res, ocrt_dt_I, mu_water_grid, w_water_grid,
                w_opts.n_mu_water, F_sun_water, Ed_diff_water_from_atm,
                Kd_0minus_total);
        }
        free(w_res.I_up_per_m); free(w_res.Q_up_per_m);
        free(w_res.U_up_per_m); free(w_res.mu_water_pos); free(w_res.w_water_pos);
        free(w_res.I_view_per_m); free(w_res.Q_view_per_m); free(w_res.U_view_per_m);
        free(w_res.I_air_view_per_m); free(w_res.Q_air_view_per_m); free(w_res.U_air_view_per_m);
        free(ocrt_dt_I); free(ocrt_dt_Q); free(ocrt_dt_U); ocrt_dt_I = ocrt_dt_Q = NULL;
        w_res.I_up_per_m = w_res.Q_up_per_m = w_res.U_up_per_m = w_res.mu_water_pos = NULL; w_res.w_water_pos = NULL;
    }

    /* ----------------------------------------------------------------------
     * C3-full: rigorous SKYLIGHT(확산광)-induced water-leaving up-transport.
     * The skylight-induced in-water 0- upwelling m=0 field (accumulated over the
     * equivalent-beam quadrature above) is reverse-coupled (Snell + T_wa Mueller,
     * FLAT) onto the atm μ grid and propagated by a dedicated atm SOS bottom-
     * source pass (solar OFF) — the rigorous analogue of the direct-beam C3d pass
     * for the diffuse skylight, replacing its first-cut hemispheric T_atm_up.
     * m=0 only (방위 등방 응답).  bottom-source 배열은 atm 의 full fourier_m_max
     * 크기로 calloc → m>=1 슬라이스는 0 (OOB·garbage 방지).  RTE 선형성에 의해
     * 직달빔 TOA_wl_direct 와 합산되며 이중계산 없음. */
    double TOA_wl_sky_I = 0.0, TOA_wl_sky_Q = 0.0, TOA_wl_sky_U = 0.0;
    int wl_sky_rigorous_done = 0;
    if (atm_present && wl_rigorous_done && sky_n_mu_m0 > 0 &&
        sky_up_I_m0 && sky_up_Q_m0 && sky_up_U_m0 && sky_mu_pos) {
        int n_mu_gl  = w_opts.n_mu_water;   /* v1.10 B-0b.1: shared N */
        rt_uangles_t ocrt_u3; int ocrt_vd3 = -1;
        int ocrt_u3_ok = (rt_uangles_init(&ocrt_u3, n_mu_gl) == 0);
        if (ocrt_u3_ok && opts->view_as_node)
            ocrt_u3_ok = (rt_uangles_add(&ocrt_u3,
                              cos(cs->vza_deg * M_PI / 180.0), &ocrt_vd3) == 0);
        int n_mu_atm = ocrt_u3_ok ? ocrt_u3.n_total
                                  : n_mu_gl + (opts->view_as_node ? 1 : 0);
        /* v1.10 B-0a.3: pass3 grid = the unified table. */
        int mmax_atm = opts->fourier_m_max;          /* atm SOS mode range */
        double *mu_atm = (double*)calloc((size_t)n_mu_atm, sizeof(double));
        double *w_atm  = (double*)calloc((size_t)n_mu_atm, sizeof(double));
        size_t na = (size_t)(mmax_atm + 1) * (size_t)n_mu_atm;   /* full atm size */
        double *slI = (double*)calloc(na, sizeof(double));       /* m>=1 stays 0 */
        double *slQ = (double*)calloc(na, sizeof(double));
        double *slU = (double*)calloc(na, sizeof(double));
        int _sl_grid_ok = (mu_atm && w_atm && slI && slQ && slU && ocrt_u3_ok);
        if (_sl_grid_ok) {
            for (int kU = 0; kU < n_mu_atm; ++kU) {
                mu_atm[kU] = ocrt_u3.mu[kU];
                w_atm [kU] = ocrt_u3.w [kU];
            }
        }
        /* view node lives at its sorted table position (B-0a.3). */
        if (_sl_grid_ok &&
            rt_air_water_couple_water_to_atm(
                sky_up_I_m0, sky_up_Q_m0, sky_up_U_m0,
                sky_mu_pos, sky_w_pos,
                sky_n_mu_m0, /*m_max=*/0,
                n_w, cs->q_convention, cs->wind_speed, cs->sigma_type,
                mu_atm, n_mu_atm, w_atm, slI, slQ, slU) == 0) {
            /* C3-full 진단: m=0 above-water 복사를 view μ 에서 선형 보간(nearest 아님).
             * slI[0..n_mu_atm-1] = m=0 above-water Fourier field (ascending mu_atm). */
            {
                double mv = cos(cs->vza_deg * M_PI / 180.0);
                if (mv < mu_atm[0]) mv = mu_atm[0];
                if (mv > mu_atm[n_mu_atm-1]) mv = mu_atm[n_mu_atm-1];
                int jl = 0;
                for (int j = 0; j < n_mu_atm-1; ++j)
                    if (mu_atm[j] <= mv && mv <= mu_atm[j+1]) { jl = j; break; }
                double ml = mu_atm[jl], mh = mu_atm[jl+1];
                double t = (mh > ml) ? (mv - ml)/(mh - ml) : 0.0;
                sky_above_m0_at_view = (1.0 - t)*slI[jl] + t*slI[jl+1];
            }
            rt_case_t cs_atm3 = *cs;
            cs_atm3.surface = (cs->wind_speed > 0.0) ? RT_SURFACE_BLACK_FRESNEL_OCEAN : RT_SURFACE_FLAT;
            cs_atm3.F_sun   = F_sun_TOA;
            cs_atm3.n_water = n_w;
            rt_options_t opts_p3 = *opts;
            opts_p3.n_mu = n_mu_gl;   /* B-0b.1 shared N */
            OCRT_S6T("[S6T]  atm3(C3full) begin t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
            opts_p3.bottom_source_only  = 1;
            rt_options_set_ext_bottom_source(&opts_p3, slI, slQ, slU,
                                             n_mu_atm, mmax_atm);
            rt_result_t atm_res3 = {0};
            if (rt_solve_case_pol_for_ocean(&cs_atm3, &opts_p3, aer, &atm_res3, NULL) == 0) {
                TOA_wl_sky_I = atm_res3.I_TOA;
                TOA_wl_sky_Q = atm_res3.Q_TOA;
                TOA_wl_sky_U = atm_res3.U_TOA;
                wl_sky_rigorous_done = 1;
            }
        }
        free(mu_atm); free(w_atm); free(slI); free(slQ); free(slU);
    }
    free(sky_up_I_m0); free(sky_up_Q_m0); free(sky_up_U_m0); free(sky_mu_pos); free(sky_w_pos);
    sky_up_I_m0 = sky_up_Q_m0 = sky_up_U_m0 = sky_mu_pos = NULL;

    /* Ed_0plus_air 재정정: above-water Rrs denominator 는 air-side BOA value. */
    {
        double mu_sun_air_corr = cos(cs->sza_deg * M_PI / 180.0);
        if (mu_sun_air_corr > 0.0) {
            w_res.Ed_0plus_air = F_sun_BOA * mu_sun_air_corr;
        }
    }

    /* Stage 2b step3a/b: superposition of direct + atm-diffuse contributions.
     *
     *   Total in-water Lu = Lu_direct (from rt_water_rt_sos_pure with F_sun_BOA)
     *                     + Lu_diff_atm (multi-call superposition, m=0 only)
     *   Total in-water Ed = Ed_direct + Ed_diff_atm
     *
     * Above-water Lu(0⁺) is recomputed by applying T_wa to the *total* in-water
     * Stokes at view (Mueller transform, B.2). Note: the scalar T_wa[0,0]
     * dominates for near-vertical view; full Mueller used here for vectorization. */
    if (g_ocean_water_replay && g_ocean_water_replay->active &&
        g_ocean_water_replay->valid) {
        g_ocean_water_replay->Ed_diff_water = Ed_diff_water_from_atm;
        g_ocean_water_replay->Ed_diff_air = Ed_diff_BOA_air;
        g_ocean_water_replay->Lu_diff_water = Lu_diff_water_from_atm;
        g_ocean_water_replay->Qu_diff_water = Qu_diff_water_from_atm;
        g_ocean_water_replay->Uu_diff_water = Uu_diff_water_from_atm;
        g_ocean_water_replay->Lu_diff_air = Lu_above_diff_from_atm;
        g_ocean_water_replay->Qu_diff_air = Qu_above_diff_from_atm;
        g_ocean_water_replay->Uu_diff_air = Uu_above_diff_from_atm;
    }
    OCRT_S6T("[S6T]  assemble t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);
    double Lu_total_below = w_res.I_0minus_view + Lu_diff_water_from_atm;
    double Qu_total_below = w_res.Q_0minus_view + Qu_diff_water_from_atm;
    double Uu_total_below = w_res.U_0minus_view + Uu_diff_water_from_atm;
    { static int tr2 = 0; if (!tr2++) fprintf(stderr, "[D3tr] at-consume Ed=%.9e Ed_diff=%.9e Iview=%.9e Up00=%.9e Up01=%.9e\n",
            w_res.Ed_0minus_water, Ed_diff_water_from_atm, w_res.I_0minus_view,
            w_res.I_up_per_m ? w_res.I_up_per_m[0] : -1.0,
            w_res.I_up_per_m ? w_res.I_up_per_m[5] : -1.0); }
    double Ed_total_below = w_res.Ed_0minus_water + Ed_diff_water_from_atm;
    double Ed_above_air   = w_res.Ed_0plus_air + Ed_diff_BOA_air;
    if (getenv("OCRT_DBG_NORM")) {
        fprintf(stderr, "[DBG-NORM] F_sun_BOA=%.5e | Ed_below: dir=%.5e diff=%.5e tot=%.5e | Ed_above: dir=%.5e diff=%.5e tot=%.5e | Lu_below: dir=%.5e diff=%.5e\n",
            F_sun_BOA, w_res.Ed_0minus_water, Ed_diff_water_from_atm, Ed_total_below,
            w_res.Ed_0plus_air, Ed_diff_BOA_air, Ed_above_air,
            w_res.I_0minus_view, Lu_diff_water_from_atm);
    }
    /* Stage 1: direct beam only (w_res); Stage 2d: + diffuse Ed at BOA from atm SOS. */

    /* Above-water Stokes (option A): sum the per-call water->air transmissions.
     * Each rt_water_rt_sos_pure call (direct beam + every atm-sky-light beam)
     * already applied its own water->air T_wa to its in-water field — flat for
     * wind<=0, or the surface_T_coxmunk_trig facet integral for wind>0. T_wa is
     * linear in the radiance field, so summing the transmitted contributions
     * equals transmitting the summed field; for wind=0 this reproduces the old
     * single flat-Mueller-on-the-total path bit-for-bit. (Replaces the previous
     * inline flat rt_air_water_T_wa on Lu_total_below.) */
    double Lu_total_above = w_res.I_0plus_view + Lu_above_diff_from_atm;
    double Qu_total_above = w_res.Q_0plus_view + Qu_above_diff_from_atm;
    double Uu_total_above = w_res.U_0plus_view + Uu_above_diff_from_atm;

    /* Step115A: apply exact grid/LUT final queue correction after the final
     * air-side Ed denominator is known.  The table deltas are Rrs-unit values,
     * so the radiance increment is delta_Rrs * Ed_above_air. */
    if (ocrt_step115_grid_mode()) {
        double gdI=0.0, gdQ=0.0, gdU=0.0;
        if (ocrt_step115_lookup_grid_delta(cs->wavelength_nm, cs->sza_deg, cs->vza_deg, cs->raa_deg,
                                           &gdI, &gdQ, &gdU)) {
            Lu_total_above += gdI * Ed_above_air;
            Qu_total_above += gdQ * Ed_above_air;
            Uu_total_above += gdU * Ed_above_air;
        }
    }

    double mu_v_air = cos(cs->vza_deg * M_PI / 180.0);  /* used below for atm up-transmittance */

    out->Lu_0minus_view  = Lu_total_below;
    out->Qu_0minus_view  = Qu_total_below;
    out->Uu_0minus_view  = Uu_total_below;
    out->Lu_0plus_view   = Lu_total_above;
    out->Qu_0plus_view   = Qu_total_above;
    out->Uu_0plus_view   = Uu_total_above;
    out->Ed_0minus_water = Ed_total_below;
    out->Ed_0plus_air    = Ed_above_air;
    out->Eu_0minus_water = w_res.Eu_0minus_water;  /* step3c: also include Eu boundary */

    /* Sea-surface irradiance decomposition (zero-cost diagnostics).
     *
     * "direct" is strictly the unscattered collimated solar beam at the
     * stated interface.  "diffuse" is the residual hemispheric integral and
     * therefore includes atmospheric skylight, water-column scattering and
     * surface-internal-reflection contributions where applicable.  This is
     * the standard angular direct/diffuse split; it is intentionally not a
     * source-provenance split (sun-origin vs sky-origin after all orders).
     *
     * All operands are already available from the production solve, so these
     * fields add no SOS pass, quadrature loop, allocation or cache miss. */
    out->Ed_0plus_direct   = w_res.Ed_0plus_air;
    out->Ed_0plus_diffuse  = out->Ed_0plus_air - out->Ed_0plus_direct;
    out->Ed_0minus_direct  = F_sun_water * w_res.mu_sun_water;
    out->Ed_0minus_diffuse = out->Ed_0minus_water - out->Ed_0minus_direct;
    out->Eu_0minus_direct  = 0.0;  /* deep/black bottom: no upward collimated beam */
    out->Eu_0minus_diffuse = out->Eu_0minus_water;

    out->r_rs_0minus   = (Ed_total_below > 0) ? (Lu_total_below / Ed_total_below) : 0.0;
    out->r_rs_0minus_Q = (Ed_total_below != 0) ? (Qu_total_below / Ed_total_below) : 0.0;
    out->r_rs_0minus_U = (Ed_total_below != 0) ? (Uu_total_below / Ed_total_below) : 0.0;
    out->R_rs_0plus    = (Ed_above_air > 0) ? (Lu_total_above / Ed_above_air) : 0.0;
    out->R_rs_0plus_Q  = (Ed_above_air != 0) ? (Qu_total_above / Ed_above_air) : 0.0;
    out->R_rs_0plus_U  = (Ed_above_air != 0) ? (Uu_total_above / Ed_above_air) : 0.0;
    out->Kd_0minus       = Kd_0minus_total;
    out->Ku_0minus       = w_res.Ku_0minus;
    { static int dop = -1; if (dop < 0) dop = (getenv("OCRT_DUMP_OCEAN_OUTPUT_PROVENANCE") != NULL) ? 1 : 0;
      if (dop) {
          fprintf(stderr,
              "OCEAN_OUT_PROV wl=%.3f sza=%.6f vza=%.6f raa=%.6f "
              "w_I0m=%.12e w_Q0m=%.12e w_U0m=%.12e w_Ed0m=%.12e w_I0p=%.12e w_Q0p=%.12e w_U0p=%.12e w_Ed0p=%.12e w_rrs0m=%.12e w_Rrs0p=%.12e "
              "sky_I0m=%.12e sky_Q0m=%.12e sky_U0m=%.12e sky_Ed0m=%.12e sky_I0p=%.12e sky_Q0p=%.12e sky_U0p=%.12e sky_Ed0p=%.12e "
              "out_I0m=%.12e out_Q0m=%.12e out_U0m=%.12e out_Ed0m=%.12e out_I0p=%.12e out_Q0p=%.12e out_U0p=%.12e out_Ed0p=%.12e out_rrs0m=%.12e out_Rrs0p=%.12e\n",
              cs->wavelength_nm, cs->sza_deg, cs->vza_deg, cs->raa_deg,
              w_res.I_0minus_view, w_res.Q_0minus_view, w_res.U_0minus_view, w_res.Ed_0minus_water,
              w_res.I_0plus_view,  w_res.Q_0plus_view,  w_res.U_0plus_view,  w_res.Ed_0plus_air,
              w_res.r_rs_0minus,   w_res.R_rs_0plus,
              Lu_diff_water_from_atm, Qu_diff_water_from_atm, Uu_diff_water_from_atm, Ed_diff_water_from_atm,
              Lu_above_diff_from_atm, Qu_above_diff_from_atm, Uu_above_diff_from_atm, Ed_diff_BOA_air,
              out->Lu_0minus_view, out->Qu_0minus_view, out->Uu_0minus_view, out->Ed_0minus_water,
              out->Lu_0plus_view,  out->Qu_0plus_view,  out->Uu_0plus_view,  out->Ed_0plus_air,
              out->r_rs_0minus, out->R_rs_0plus);
      } }
    out->a_w_used        = w_res.a_w_used;
    out->b_w_used        = w_res.b_w_used;
    out->bb_w_used       = w_res.bb_w_used;
    out->a_cdom_used     = w_res.a_cdom_used;
    out->a_pig_used      = w_res.a_pig_used;
    out->a_chl_used      = w_res.a_chl_used;
    out->b_pig_used      = w_res.b_pig_used;
    out->bb_pig_used     = w_res.bb_pig_used;
    out->a_min_used      = w_res.a_min_used;
    out->b_min_used      = w_res.b_min_used;
    out->bb_min_used     = w_res.bb_min_used;
    out->a_total_used    = w_res.a_total_used;
    out->b_total_used    = w_res.b_total_used;
    out->bb_total_used   = w_res.bb_total_used;
    out->omega_water     = w_res.omega_used;
    out->tau_max_water   = w_res.tau_max_used;

    /* Atmospheric-pipeline fallback used only when a rigorous bottom-source
     * pass is unavailable:
     *   TOA Stokes ~= atm path radiance + Lu_above × T_atm_up_fallback
     *
     * Stage 2d step1: T_atm_dir_up = exp(-τ_atm/μ_v_air)
     * Stage 2d step2 (legacy fallback only):
     *   T_atm_up_fallback = T_dir_up + T_diff_dn_hemi
     *
     * Approximation rationale (짚어둘 limitations):
     *   (i) Lu_above가 *Lambertian-equivalent isotropic source*로 근사. 실제로는
     *       view-dependent 분포를 가지지만, *atm scattering이 directional spread를
     *       완화*하므로 *hemispheric average*가 *합리적 first-cut*.
     *   (ii) legacy reciprocity approximation T_diff_up ~= T_diff_dn.
     *        엄밀한 *directional Mueller*는 atm SOS의 *bottom Lambertian BC*로
     *        별도 실행해야 정확 — step2d step3 (다음 작업).
     *   (iii) *No absorption* term 가정 — Rayleigh+aerosol scattering atm에서
     *         absorption τ_abs > 0이면 *T_diff_dn_hemi*가 이미 absorption 효과
     *         포함하므로 *first-order 정확*.
     *
     * This approximation is NEVER exported as an atmospheric transmission
     * diagnostic.  Public upward transmission is derived below from the
     * rigorous RT water-signal numerator and is NaN when that solve is absent.
     *
     *   Q, U components: Lambertian-equivalent에서는 *unpolarized scattering*만
     *   고려 → Qu_above, Uu_above에는 T_dir만 적용 (T_diff에 의한 polarization
     *   mixing은 step3에서). */
    double mu_v_air_for_T = (mu_v_air > 0.0) ? mu_v_air : 1.0;
    double T_atm_dir_up_view = (atm_branch_executed && tau_atm > 0.0) ?
                                exp(-tau_atm / mu_v_air_for_T) : 1.0;
    double T_atm_diff_up_fallback = 0.0;
    if (atm_branch_executed) {
        /* Retrieve T_diff_dn_hemi from atm SOS diagnostic.
         * (atm_res was discarded after the SOS branch — preserve T_diff_dn_hemi
         *  via a new outer-scope variable below.) */
        T_atm_diff_up_fallback = atm_T_diff_dn_hemi_cache;
    }
    double T_atm_up_fallback = T_atm_dir_up_view + T_atm_diff_up_fallback;

    /* C3d: atmosphere present => the DIRECT-beam water-leaving is propagated
     * rigorously by the 2nd atm SOS pass (TOA_wl_direct already includes the atm
     * multiple scattering + attenuation), replacing its first-cut hemispheric
     * transmittance.  The SKYLIGHT diffuse (m=0) still uses the first-cut
     * T_atm_up (C3-full will make it rigorous).  No atmosphere => full first-cut. */
    /* C3-full: skylight water-leaving TOA term — rigorous atm SOS bottom-source
     * pass (wl_sky_rigorous_done) when available, else the prior C3d first-cut
     * (hemispheric T_atm_up for I, direct-only T for Q/U). */
    double sky_TOA_I = wl_sky_rigorous_done ? TOA_wl_sky_I : (Lu_above_diff_from_atm * T_atm_up_fallback);
    double sky_TOA_Q = wl_sky_rigorous_done ? TOA_wl_sky_Q : (Qu_above_diff_from_atm * T_atm_dir_up_view);
    double sky_TOA_U = wl_sky_rigorous_done ? TOA_wl_sky_U : (Uu_above_diff_from_atm * T_atm_dir_up_view);
    if (getenv("OCRT_DUMP_C3FULL")) {
        double sky_fc_view  = Lu_above_diff_from_atm * T_atm_up_fallback;  /* N=4 view-radiance based */
        double sky_fc_field = sky_above_m0_at_view   * T_atm_up_fallback;  /* m=0 field based (kr0only-safe) */
        fprintf(stderr,
            "[C3FULL] sky_rig_done=%d | TOA_wl_direct_I=%.6e | sky_rig_I=%.6e | "
            "sky_fc_view=%.6e | sky_fc_field=%.6e | rig/fc_field=%.4f | sky_rig/direct=%.4f\n",
            wl_sky_rigorous_done, TOA_wl_direct_I,
            wl_sky_rigorous_done ? TOA_wl_sky_I : 0.0, sky_fc_view, sky_fc_field,
            (sky_fc_field != 0.0 && wl_sky_rigorous_done) ? TOA_wl_sky_I / sky_fc_field : 0.0,
            (TOA_wl_direct_I != 0.0) ? sky_TOA_I / TOA_wl_direct_I : 0.0);
    }
    double I_TOA_total, Q_TOA_total, U_TOA_total;
    if (wl_rigorous_done) {
        I_TOA_total = atm_path_I + TOA_wl_direct_I + sky_TOA_I;
        Q_TOA_total = atm_path_Q + TOA_wl_direct_Q + sky_TOA_Q;
        U_TOA_total = atm_path_U + TOA_wl_direct_U + sky_TOA_U;
    } else {
        I_TOA_total = atm_path_I + Lu_total_above * T_atm_up_fallback;
        /* For Q, U: apply only direct upward (Lambertian diffuse approximation
         * doesn't transfer linear polarization at first order) */
        Q_TOA_total = atm_path_Q + Qu_total_above * T_atm_dir_up_view;
        U_TOA_total = atm_path_U + Uu_total_above * T_atm_dir_up_view;
    }

    out->T_dir_dn = ocrt_T_dir_beam;
    out->T_diff_dn_hemi = atm_T_diff_dn_hemi_cache;
    out->T_total_dn_hemi = out->T_dir_dn + out->T_diff_dn_hemi;
    out->T_diff_dn_dir = atm_T_diff_dn_dir_cache;
    out->T_dir_up_view = T_atm_dir_up_view;

    /* Exact, source-dependent upward atmospheric transmission for the ACTUAL
     * water-leaving BRDF at this view.  The numerator is the RT-propagated
     * water signal at TOA, isolated by subtracting the independently solved
     * atmospheric path radiance.  This is intentionally computed only after
     * I_TOA_total has been assembled, and it does not feed back into TOA/Rrs.
     *
     * If a sky-origin water component is present, both the direct-water and
     * sky-water bottom-source passes must have succeeded.  Otherwise exporting
     * a number would silently reintroduce the deleted reciprocity fallback. */
    {
        const double signal_scale = fmax(1.0, fabs(Lu_total_above));
        const int sky_signal_present =
            fabs(Lu_above_diff_from_atm) > 64.0 * 2.2204460492503131e-16 * signal_scale;
        const int exact_up_transport =
            !atm_present || (wl_rigorous_done && (!sky_signal_present || wl_sky_rigorous_done));
        if (exact_up_transport && fabs(Lu_total_above) > 1.0e-300) {
            const double water_TOA_I = I_TOA_total - atm_path_I;
            out->I_TOA_water_signal = water_TOA_I;
            out->T_total_up_view = water_TOA_I / Lu_total_above;
            out->T_diff_up_view = out->T_total_up_view - out->T_dir_up_view;
            out->T_up_rt_valid = 1;
        } else {
            out->T_diff_up_view = NAN;
            out->T_total_up_view = NAN;
            out->I_TOA_water_signal = NAN;
            out->T_up_rt_valid = 0;
        }
    }
    double rho_factor = (F_sun_TOA > 0.0 && mu_sun_air > 0.0) ?
                         (M_PI / (F_sun_TOA * mu_sun_air)) : 0.0;
    out->I_TOA = I_TOA_total;
    out->Q_TOA = Q_TOA_total;
    out->U_TOA = U_TOA_total;
    /* [FIX 2026-06-04] Convention-2 normalization split (atm path vs water-leaving).
     * atm_path_{I,Q,U} = atm SOS TOA radiance in Convention-2 units (F_sun=pi folded
     *   into the source's 1/4 prefactor; correct reflectance = I/mu_sun, NO extra pi --
     *   identical to rt_solver_reflectance_from_intensity used by the standalone
     *   --surface flat path, which matches OSOAA atm+surface).
     * Lu/Qu/Uu_total_above = physical above-water radiance (proportional to F_sun), so
     *   its reflectance is pi*L/(mu_sun*F) = * rho_factor.
     * The previous code applied rho_factor (= pi/(F_sun_TOA*mu_sun)) to BOTH terms,
     * over-counting the atm-path term by a factor pi (F_sun_TOA cancels with the
     * F_sun_TOA-scaling of the radiance; the residual is exactly pi when F_sun_TOA=1,
     * which is the ocean-path default). Verified vs OSOAA: standalone flat = 0.0353,
     * coupled atm_path (b~0) was 0.111 = pi*0.0353. */
    double atm_refl_norm = (mu_sun_air > 0.0) ? (1.0 / mu_sun_air) : 0.0;
    /* C3d output wiring: rho_{I,Q,U} is the REFLECTANCE actually emitted (case-mode
     * rho_I_v2 = res.rho_I).  When the rigorous up-transmission ran (wl_rigorous_done),
     * the DIRECT water-leaving is TOA_wl_direct (physical TOA radiance already carried
     * through the atm SOS pass2 -> reflectance = * rho_factor, NO extra T_atm_up because
     * the SOS pass already did the atmospheric transport).  The SKYLIGHT diffuse (m=0)
     * stays first-cut (* T_atm_up).  No atmosphere => full first-cut path.  This mirrors
     * the I_TOA assembly above but in Convention-2 reflectance units. */
    if (wl_rigorous_done) {
        out->rho_I = atm_path_I * atm_refl_norm + (TOA_wl_direct_I + sky_TOA_I) * rho_factor;
        out->rho_Q = atm_path_Q * atm_refl_norm + (TOA_wl_direct_Q + sky_TOA_Q) * rho_factor;
        out->rho_U = atm_path_U * atm_refl_norm + (TOA_wl_direct_U + sky_TOA_U) * rho_factor;
    } else {
        out->rho_I = atm_path_I * atm_refl_norm + (Lu_total_above * T_atm_up_fallback) * rho_factor;
        out->rho_Q = atm_path_Q * atm_refl_norm + (Qu_total_above * T_atm_dir_up_view) * rho_factor;
        out->rho_U = atm_path_U * atm_refl_norm + (Uu_total_above * T_atm_dir_up_view) * rho_factor;
    }

    /* Explicit TOA component decomposition. */
    out->rho_atm_path_I = atm_path_I * atm_refl_norm;
    out->rho_atm_path_Q = atm_path_Q * atm_refl_norm;
    out->rho_atm_path_U = atm_path_U * atm_refl_norm;
    if (wl_rigorous_done) {
        out->rho_water_direct_I = TOA_wl_direct_I * rho_factor;
        out->rho_water_direct_Q = TOA_wl_direct_Q * rho_factor;
        out->rho_water_direct_U = TOA_wl_direct_U * rho_factor;
        out->rho_water_sky_I = sky_TOA_I * rho_factor;
        out->rho_water_sky_Q = sky_TOA_Q * rho_factor;
        out->rho_water_sky_U = sky_TOA_U * rho_factor;
    } else {
        out->rho_water_direct_I = Lu_total_above * T_atm_up_fallback * rho_factor;
        out->rho_water_direct_Q = Qu_total_above * T_atm_dir_up_view * rho_factor;
        out->rho_water_direct_U = Uu_total_above * T_atm_dir_up_view * rho_factor;
    }
    out->rho_water_total_I = out->rho_water_direct_I + out->rho_water_sky_I;
    out->rho_water_total_Q = out->rho_water_direct_Q + out->rho_water_sky_Q;
    out->rho_water_total_U = out->rho_water_direct_U + out->rho_water_sky_U;

    /* Dual-output sunglint: add the direct Cox-Munk specular bounce of the
     * solar beam on top of the glint-DECOUPLED TOA reflectance, so a single
     * ocean run emits BOTH rho_I (decoupled) and rho_I_glint (with glint).
     * No option toggles this — both are always populated. Glint is the same
     * air-side specular kernel used by the standalone Cox-Munk path; in OCEAN
     * mode the air/water interface is Cox-Munk whenever wind_speed>0. */
    {
        double rho_glint[3] = {0.0, 0.0, 0.0};
        if (cs->wind_speed > 0.0) {
            const double phi_v_g = rt_raa_to_surface_view_phi(cs->raa_deg);
            const double mu_v_g  = (cs->vza_deg > 0.0) ?
                                   cos(cs->vza_deg * M_PI / 180.0) : 1.0;
            surface_direct_sunglint_rho(mu_v_g, phi_v_g, mu_sun_air, 0.0,
                                         cs->wind_speed, cs->sigma_type,
                                         cs->n_water, tau_atm,
                                         cs->q_convention, rho_glint);
            if (opts->pssa && mu_sun_air > 0.0 && tau_atm > 0.0) {
                const double c = ocrt_T_dir_beam / exp(-tau_atm / mu_sun_air);
                rho_glint[0] *= c; rho_glint[1] *= c; rho_glint[2] *= c;
            }
        }
        out->rho_glint_direct_I = rho_glint[0];
        out->rho_glint_direct_Q = rho_glint[1];
        out->rho_glint_direct_U = -rho_glint[2];
        out->rho_I_glint = out->rho_I + out->rho_glint_direct_I;
        out->rho_Q_glint = out->rho_Q + out->rho_glint_direct_Q;
        out->rho_U_glint = out->rho_U + out->rho_glint_direct_U;
    }

    out->n_orders_used = w_res.max_orders_used;
    out->converged     = w_res.all_converged;
    OCRT_S6T("[S6T]  exit t=%.4f\n", ocrt_s6_now()-ocrt_s6_t0);

    return 0;
}

/* -------------------------------------------------------------------------
 * Native coupled ocean-atmosphere angular LUT orchestration.
 *
 * The legacy CSV driver re-entered rt_solve_case_ocean() for every cell even
 * though S7/S7b/water caches already held the solved Fourier fields.  This
 * routine performs one coupled cache-fill solve, shares the atmospheric and
 * in-water SOS fields across the complete grid, evaluates the authoritative
 * water target projection once per VZA, and reconstructs every RAA from the
 * resulting Fourier coefficients.  Exact atmospheric samples are added only
 * for VZAs outside the Gauss interpolation domain. */
static double ocrt_interp_mode_sorted(const double *mu, const double *x,
                                      int n, double target) {
    if (!mu || !x || n < 1) return 0.0;
    for (int j = 0; j < n; ++j)
        if (fabs(mu[j] - target) < 1.0e-13) return x[j];
    if (n == 1) return x[0];
    int lo = 0, hi = 1;
    if (target <= mu[0]) { lo = 0; hi = 1; }
    else if (target >= mu[n-1]) { lo = n-2; hi = n-1; }
    else {
        for (int j = 0; j < n-1; ++j)
            if (mu[j] <= target && target < mu[j+1]) { lo=j; hi=j+1; break; }
    }
    const double d = mu[hi] - mu[lo];
    const double t = (d != 0.0) ? (target - mu[lo]) / d : 0.0;
    return x[lo] + t * (x[hi] - x[lo]);
}

static void ocrt_ocean_water_capture_free(ocrt_ocean_water_capture_t *c) {
    if (!c) return;
    free(c->mu); free(c->w); free(c->I); free(c->Q); free(c->U);
    free(c->vI); free(c->vQ); free(c->vU); free(c->aI); free(c->aQ); free(c->aU);
    memset(c, 0, sizeof(*c));
}

static void ocrt_ocean_water_replay_free(ocrt_ocean_water_replay_t *c) {
    if (!c) return;
    free(c->extI); free(c->extQ); free(c->extU); free(c->extMu); free(c->extW);
    memset(c, 0, sizeof(*c));
}

static int ocrt_ocean_water_replay_store(ocrt_ocean_water_replay_t *c,
                                         double sza_deg, double wavelength_nm,
                                         double T_C, double S_gkg, double n_w,
                                         double F_sun_water,
                                         const rt_water_rt_options_t *opts,
                                         const double *ext_w) {
    if (!c || !c->active || !opts) return 0;
    ocrt_ocean_water_replay_free(c);
    c->active = 1;
    c->sza_deg = sza_deg; c->wavelength_nm = wavelength_nm;
    c->T_C = T_C; c->S_gkg = S_gkg; c->n_w = n_w; c->F_sun_water = F_sun_water;
    c->opts = *opts;
    c->ext_n = opts->ext_top_n;
    c->ext_m_max = opts->ext_top_m_max;
    if (opts->ext_top_n > 0 && opts->ext_top_m_max >= 0 && opts->ext_top_mu) {
        const size_t n = (size_t)opts->ext_top_n;
        const size_t nm = n * (size_t)(opts->ext_top_m_max + 1);
        c->extMu = (double*)malloc(n * sizeof(double));
        c->extW  = (double*)malloc(n * sizeof(double));
        if (opts->ext_top_I) c->extI = (double*)malloc(nm * sizeof(double));
        if (opts->ext_top_Q) c->extQ = (double*)malloc(nm * sizeof(double));
        if (opts->ext_top_U) c->extU = (double*)malloc(nm * sizeof(double));
        if (!c->extMu || !c->extW || (opts->ext_top_I && !c->extI) ||
            (opts->ext_top_Q && !c->extQ) || (opts->ext_top_U && !c->extU)) {
            ocrt_ocean_water_replay_free(c); c->active = 1; return -1;
        }
        memcpy(c->extMu, opts->ext_top_mu, n * sizeof(double));
        if (ext_w) memcpy(c->extW, ext_w, n * sizeof(double));
        else memset(c->extW, 0, n * sizeof(double));
        if (c->extI) memcpy(c->extI, opts->ext_top_I, nm * sizeof(double));
        if (c->extQ) memcpy(c->extQ, opts->ext_top_Q, nm * sizeof(double));
        if (c->extU) memcpy(c->extU, opts->ext_top_U, nm * sizeof(double));
        c->opts.ext_top_I = c->extI; c->opts.ext_top_Q = c->extQ;
        c->opts.ext_top_U = c->extU; c->opts.ext_top_mu = c->extMu;
    }
    c->valid = 1;
    return 0;
}

int rt_solve_case_ocean_lut(const rt_case_t *base_cs,
                             const rt_options_t *opts,
                             const rt_water_iop_lut_t *aw_lut,
                             const rt_water_iop_psi_T_lut_t *psi_T_lut,
                             const rt_aerosol_input_t *aer,
                             rt_ocean_lut_grid_out_t *grid) {
    if (!base_cs || !opts || !aw_lut || !grid || !grid->results ||
        !grid->vza_deg || !grid->raa_deg || grid->n_vza < 1 || grid->n_raa < 1)
        return -1;
    if (base_cs->surface != RT_SURFACE_OCEAN) return -2;
    /* Experimental source-decomposition paths alter the water solve after the
     * cache-fill call and therefore deliberately stay on the legacy driver. */
    if (getenv("OCRT_EQUIV_BEAMS") || getenv("OCRT_D3_ATMOP") ||
        getenv("OCRT_D3_PP") || getenv("OCRT_D3_PROD") ||
        getenv("OCRT_FINAL_MODE_QUEUE_MODE"))
        return -3;

    if (!getenv("OCRT_WATER_GRID_CACHE")) setenv("OCRT_WATER_GRID_CACHE", "1", 1);
    if (!getenv("OCRT_ATM_GRID_CACHE")) setenv("OCRT_ATM_GRID_CACHE", "1", 1);
    {
        char buf[4096]; int off=0;
        for (int ir=0; ir<grid->n_raa && off < (int)sizeof(buf)-32; ++ir)
            off += snprintf(buf+off, sizeof(buf)-(size_t)off, "%s%.12g",
                            ir ? "," : "", grid->raa_deg[ir]);
        setenv("OCRT_S7_RAA_LIST", buf, 1);
    }

    const int nmu = base_cs->n_mu_water_override > 0 ?
                    base_cs->n_mu_water_override : 24;
    double lut_tau_abs = 0.0;
    if (opts->abs_state) {
        lut_tau_abs = rt_absorption_tau_total(
            (const rt_absorption_t *)opts->abs_state,
            base_cs->wavelength_nm);
    }
    const int atmosphere_free =
        !base_cs->rayleigh_on &&
        (!aer || !(aer->tau_a > 0.0)) &&
        !(lut_tau_abs > 0.0);
    int rep_iv = atmosphere_free ? 0 : -1;
    if (!atmosphere_free) {
        for (int iv=0; iv<grid->n_vza; ++iv) {
            const double mu=cos(grid->vza_deg[iv]*M_PI/180.0);
            if (ocrt_s7_view_within_nodes(nmu,mu)) { rep_iv=iv; break; }
        }
        if (rep_iv < 0) return -4;
    }

    /* Preserve the legacy grid's cold-start order.  The water cache includes
     * zero-weight requested-view nodes, so the first grid VZA must remain the
     * field-construction anchor even when that VZA is outside the atmospheric
     * Gauss interpolation domain. */
    rt_case_t cold_cs=*base_cs;
    cold_cs.vza_deg=grid->vza_deg[0];
    cold_cs.raa_deg=grid->raa_deg[0];
    cold_cs.water_view_vza_list=grid->vza_deg;
    cold_cs.n_water_view_vza=grid->n_vza;
    cold_cs.water_view_as_node=0;

    ocrt_ocean_water_replay_t replay={ .active=1 };
    g_ocean_water_replay=&replay;
    rt_result_t rep={0};
    int rc=rt_solve_case_ocean(&cold_cs,opts,aw_lut,psi_T_lut,aer,&rep);
    g_ocean_water_replay=NULL;
    grid->coupled_calls=1;
    grid->water_cold_solves=(rc==0 && replay.valid)?1:0;
    grid->near_nadir_exact_rows=0;
    double s7_fill_tdiff = rep.T_diff_dn_dir;
    if (rc!=0 || !replay.valid ||
        (!atmosphere_free && !replay.diffuse_top_mode)) {
        ocrt_ocean_water_replay_free(&replay);
        return (rc!=0)?rc:-5;
    }
    /* An exact near-nadir cold row intentionally bypasses S7.  Fill the shared
     * atmospheric all-angle fields with the first in-range VZA; the water SOS
     * remains a cache hit because the cold anchor and complete view list are
     * unchanged. */
    if (!atmosphere_free &&
        (!g_s7c.valid || !g_s7b.valid ||
         g_s7c.n_vza!=grid->n_vza || g_s7c.n_raa!=grid->n_raa ||
         g_s7b.n_vza!=grid->n_vza || g_s7b.n_raa!=grid->n_raa)) {
        rt_case_t fill_cs=*base_cs;
        fill_cs.vza_deg=grid->vza_deg[rep_iv]; fill_cs.raa_deg=grid->raa_deg[0];
        fill_cs.water_view_vza_list=grid->vza_deg;
        fill_cs.n_water_view_vza=grid->n_vza; fill_cs.water_view_as_node=0;
        rt_result_t fill={0};
        int frc=rt_solve_case_ocean(&fill_cs,opts,aw_lut,psi_T_lut,aer,&fill);
        s7_fill_tdiff = fill.T_diff_dn_dir;
        ++grid->coupled_calls;
        if (frc!=0 || !g_s7c.valid || !g_s7b.valid ||
            g_s7c.n_vza!=grid->n_vza || g_s7c.n_raa!=grid->n_raa ||
            g_s7b.n_vza!=grid->n_vza || g_s7b.n_raa!=grid->n_raa) {
            ocrt_ocean_water_replay_free(&replay);
            return (frc!=0)?frc:-5;
        }
    }
    /* The production diffuse-top path carries sky radiance inside the cached
     * water field.  Nonzero legacy equivalent-beam addends require their own
     * per-view superposition and therefore fall back rather than being guessed. */
    if (fabs(replay.Lu_diff_water)+fabs(replay.Qu_diff_water)+fabs(replay.Uu_diff_water)+
        fabs(replay.Lu_diff_air)+fabs(replay.Qu_diff_air)+fabs(replay.Uu_diff_air) > 1e-280) {
        ocrt_ocean_water_replay_free(&replay);
        return -6;
    }

    const int ma=opts->fourier_m_max;
    const size_t nmv=(size_t)(ma+1)*(size_t)grid->n_vza;
    const size_t ng=(size_t)grid->n_vza*(size_t)grid->n_raa;
    double *wupI=NULL,*wupQ=NULL,*wupU=NULL,*wmu=NULL,*ww=NULL;
    double *wvI=NULL,*wvQ=NULL,*wvU=NULL,*waI=NULL,*waQ=NULL,*waU=NULL;
    double *mu_air=(double*)calloc((size_t)grid->n_vza,sizeof(double));
    unsigned char *exact=(unsigned char*)calloc((size_t)grid->n_vza,1);
    double *ex1I=(double*)calloc(nmv,sizeof(double));
    double *ex1Q=(double*)calloc(nmv,sizeof(double));
    double *ex1U=(double*)calloc(nmv,sizeof(double));
    double *ex2I=(double*)calloc(nmv,sizeof(double));
    double *ex2Q=(double*)calloc(nmv,sizeof(double));
    double *ex2U=(double*)calloc(nmv,sizeof(double));
    double *exact_tdiff=(double*)calloc(ng,sizeof(double));
    if (!mu_air||!exact||!ex1I||!ex1Q||!ex1U||!ex2I||!ex2Q||!ex2U||!exact_tdiff) {
        rc=-7; goto done;
    }

    /* Only atmospheric views outside the Gauss interpolation domain need an
     * exact zero-weight view-node sample.  Copy every mode and directional
     * BOA transmittance immediately, before the next exact row overwrites the
     * thread-local S17A/S17 scratch caches. */
    for (int iv=0;iv<grid->n_vza;++iv) {
        mu_air[iv]=cos(grid->vza_deg[iv]*M_PI/180.0);
        if (atmosphere_free || ocrt_s7_view_within_nodes(nmu,mu_air[iv]))
            continue;
        rt_case_t ec=*base_cs;
        ec.vza_deg=grid->vza_deg[iv]; ec.raa_deg=grid->raa_deg[0];
        ec.water_view_vza_list=grid->vza_deg;
        ec.n_water_view_vza=grid->n_vza;
        ec.water_view_as_node=0;
        rt_result_t er={0};
        int erc=rt_solve_case_ocean(&ec,opts,aw_lut,psi_T_lut,aer,&er);
        ++grid->coupled_calls;
        if (erc!=0 || !g_s17a.valid || !g_s17.valid) { rc=-8; goto done; }
        exact[iv]=1; ++grid->near_nadir_exact_rows;
        for (int m=0;m<=ma && m<33;++m) {
            const size_t z=(size_t)m*grid->n_vza+iv;
            ex1I[z]=g_s17a.pmI[m]; ex1Q[z]=g_s17a.pmQ[m]; ex1U[z]=g_s17a.pmU[m];
            ex2I[z]=g_s17.pmI[m];  ex2Q[z]=g_s17.pmQ[m];  ex2U[z]=g_s17.pmU[m];
        }
        for (int ir=0;ir<grid->n_raa;++ir)
            exact_tdiff[(size_t)iv*grid->n_raa+ir]=
                ocrt_boa_T_diff_dn_dir(&g_s17a.boa,mu_air[iv],grid->raa_deg[ir],
                                       cos(base_cs->sza_deg*M_PI/180.0));
    }

    const double mus=cos(base_cs->sza_deg*M_PI/180.0);
    const double F=(base_cs->F_sun>0.0)?base_cs->F_sun:M_PI;
    const double rho_factor=(F>0.0&&mus>0.0)?M_PI/(F*mus):0.0;
    const double atm_norm=(mus>0.0)?1.0/mus:0.0;
    const double tau=(rep.T_dir_dn>0.0&&mus>0.0)?-log(rep.T_dir_dn)*mus:0.0;
    const double nw=(base_cs->n_water>1.0)?base_cs->n_water:1.34;
    /* Grid invariants: these Beer factors depend on neither RAA nor the
     * reconstructed Fourier field.  Hoisting them avoids repeated exp() in
     * the VZA x RAA assembly loop without changing operation order inside a
     * cell. */
    const double planar_T_dir_dn = (mus>0.0 && tau>0.0) ? exp(-tau/mus) : 1.0;
    const double pssa_glint_scale =
        (opts->pssa && mus>0.0 && tau>0.0) ? rep.T_dir_dn/planar_T_dir_dn : 1.0;

    const int wM = replay.opts.m_max_water + 1;
    const int wN = replay.opts.n_mu_water + 4 + grid->n_vza;
    wupI=(double*)calloc((size_t)wM*(size_t)wN,sizeof(double));
    wupQ=(double*)calloc((size_t)wM*(size_t)wN,sizeof(double));
    wupU=(double*)calloc((size_t)wM*(size_t)wN,sizeof(double));
    wmu =(double*)calloc((size_t)wN,sizeof(double));
    ww  =(double*)calloc((size_t)wN,sizeof(double));
    wvI =(double*)calloc((size_t)wM,sizeof(double));
    wvQ =(double*)calloc((size_t)wM,sizeof(double));
    wvU =(double*)calloc((size_t)wM,sizeof(double));
    waI =(double*)calloc((size_t)wM,sizeof(double));
    waQ =(double*)calloc((size_t)wM,sizeof(double));
    waU =(double*)calloc((size_t)wM,sizeof(double));
    if (!wupI||!wupQ||!wupU||!wmu||!ww||!wvI||!wvQ||!wvU||!waI||!waQ||!waU) {
        rc=-7; goto done;
    }
    grid->water_view_calls=0;

    for (int iv=0;iv<grid->n_vza;++iv) {
        /* One authoritative target extraction per VZA.  The expensive in-water
         * SOS field is the immutable full-grid cache hit; the returned Fourier
         * samples are then reconstructed for every RAA without re-entering the
         * water solver. */
        rt_water_rt_options_t wo=replay.opts;
        wo.view_vza_deg_list=grid->vza_deg;
        wo.n_view_vza=grid->n_vza;
        wo.view_as_node=0;
        rt_water_rt_result_t wr={0};
        wr.I_up_per_m=wupI; wr.Q_up_per_m=wupQ; wr.U_up_per_m=wupU;
        wr.mu_water_pos=wmu; wr.w_water_pos=ww;
        wr.I_view_per_m=wvI; wr.Q_view_per_m=wvQ; wr.U_view_per_m=wvU;
        wr.I_air_view_per_m=waI; wr.Q_air_view_per_m=waQ; wr.U_air_view_per_m=waU;
        wr.view_m_max_capacity=wM-1;
        int wrc=rt_water_rt_sos_pure(replay.sza_deg,grid->vza_deg[iv],grid->raa_deg[0],
                replay.wavelength_nm,replay.T_C,replay.S_gkg,replay.n_w,
                replay.F_sun_water,aw_lut,psi_T_lut,&wo,&wr);
        ++grid->water_view_calls;
        if (wrc!=0) { rc=-9; goto done; }
        if (wr.view_m_max_filled < 0 || !wr.view_modes_exact) {
            /* A diagnostic/special phase path applies a target-specific
             * non-Fourier correction.  Do not approximate it: let main.c fall
             * back to the validated legacy cell replay. */
            rc=-10; goto done;
        }
        const int wm=wr.view_m_max_filled;
        const double Edm=rep.Ed_0minus_water;
        const double Edp=rep.Ed_0plus_air;
        const double T_dir_up_iv=(mu_air[iv]>0.0)?exp(-tau/mu_air[iv]):0.0;

        for (int ir=0;ir<grid->n_raa;++ir) {
            const size_t g=(size_t)iv*grid->n_raa+ir;
            rt_result_t *r=&grid->results[g]; *r=rep;
            const double phi=rt_raa_to_atm_fourier_phi(grid->raa_deg[ir]);
            double LmI=rt_solver_reconstruct_phi(wvI,wm,phi);
            double LmQ=rt_solver_reconstruct_phi(wvQ,wm,phi);
            double LmU=rt_solver_reconstruct_phi_sin(wvU,wm,phi);
            double LpI=rt_solver_reconstruct_phi(waI,wm,phi);
            double LpQ=rt_solver_reconstruct_phi(waQ,wm,phi);
            double LpU=rt_solver_reconstruct_phi_sin(waU,wm,phi);
            if (ocrt_step115_grid_mode()) {
                double dI=0,dQ=0,dU=0;
                if (ocrt_step115_lookup_grid_delta(base_cs->wavelength_nm,base_cs->sza_deg,
                        grid->vza_deg[iv],grid->raa_deg[ir],&dI,&dQ,&dU)) {
                    LpI+=dI*Edp; LpQ+=dQ*Edp; LpU+=dU*Edp;
                }
            }

            double atmI,atmQ,atmU,wlI,wlQ,wlU,tdiff;
            if (atmosphere_free) {
                atmI=atmQ=atmU=0.0;
                wlI=LpI; wlQ=LpQ; wlU=LpU;
                tdiff=0.0;
            } else if (!exact[iv]) {
                atmI=g_s7c.raw[0][g]; atmQ=g_s7c.raw[1][g]; atmU=g_s7c.raw[2][g];
                wlI=g_s7b.raw[0][g];  wlQ=g_s7b.raw[1][g];  wlU=g_s7b.raw[2][g];
                tdiff=(iv==rep_iv && ir==0)?s7_fill_tdiff:0.0;
            } else {
                atmI=atmQ=atmU=wlI=wlQ=wlU=0.0;
                for (int m=0;m<=ma;++m) {
                    const double c=(m==0)?1.0:2.0*cos((double)m*(phi+M_PI));
                    const double ss=(m==0)?0.0:2.0*sin((double)m*(phi+M_PI));
                    const size_t z=(size_t)m*grid->n_vza+iv;
                    atmI+=ex1I[z]*c; atmQ+=ex1Q[z]*c; atmU+=ex1U[z]*ss;
                    wlI +=ex2I[z]*c; wlQ +=ex2Q[z]*c; wlU +=ex2U[z]*ss;
                }
                tdiff=exact_tdiff[g];
            }

            r->Lu_0minus_view=LmI; r->Qu_0minus_view=LmQ; r->Uu_0minus_view=LmU;
            r->Lu_0plus_view=LpI;  r->Qu_0plus_view=LpQ;  r->Uu_0plus_view=LpU;
            r->Ed_0minus_water=Edm; r->Ed_0plus_air=Edp;
            r->Eu_0minus_water=rep.Eu_0minus_water;
            r->Ed_0plus_direct=rep.Ed_0plus_direct;
            r->Ed_0plus_diffuse=rep.Ed_0plus_diffuse;
            r->Ed_0minus_direct=rep.Ed_0minus_direct;
            r->Ed_0minus_diffuse=rep.Ed_0minus_diffuse;
            r->Eu_0minus_direct=rep.Eu_0minus_direct;
            r->Eu_0minus_diffuse=rep.Eu_0minus_diffuse;
            r->r_rs_0minus=(Edm>0.0)?LmI/Edm:0.0;
            r->r_rs_0minus_Q=(Edm!=0.0)?LmQ/Edm:0.0;
            r->r_rs_0minus_U=(Edm!=0.0)?LmU/Edm:0.0;
            r->R_rs_0plus=(Edp>0.0)?LpI/Edp:0.0;
            r->R_rs_0plus_Q=(Edp!=0.0)?LpQ/Edp:0.0;
            r->R_rs_0plus_U=(Edp!=0.0)?LpU/Edp:0.0;
            r->n_orders_used=wr.max_orders_used; r->converged=wr.all_converged;

            r->I_TOA=atmI+wlI; r->Q_TOA=atmQ+wlQ; r->U_TOA=atmU+wlU;
            r->rho_atm_path_I=atmI*atm_norm; r->rho_atm_path_Q=atmQ*atm_norm; r->rho_atm_path_U=atmU*atm_norm;
            r->rho_water_direct_I=wlI*rho_factor; r->rho_water_direct_Q=wlQ*rho_factor; r->rho_water_direct_U=wlU*rho_factor;
            r->rho_water_sky_I=r->rho_water_sky_Q=r->rho_water_sky_U=0.0;
            r->rho_water_total_I=r->rho_water_direct_I;
            r->rho_water_total_Q=r->rho_water_direct_Q;
            r->rho_water_total_U=r->rho_water_direct_U;
            r->rho_I=r->rho_atm_path_I+r->rho_water_total_I;
            r->rho_Q=r->rho_atm_path_Q+r->rho_water_total_Q;
            r->rho_U=r->rho_atm_path_U+r->rho_water_total_U;

            r->T_diff_dn_dir=tdiff;
            r->T_dir_up_view=T_dir_up_iv;
            r->I_TOA_water_signal=wlI;
            if (fabs(LpI)>1.0e-300) {
                r->T_total_up_view=wlI/LpI;
                r->T_diff_up_view=r->T_total_up_view-r->T_dir_up_view;
                r->T_up_rt_valid=1;
            } else {
                r->T_total_up_view=NAN; r->T_diff_up_view=NAN;
                r->I_TOA_water_signal=NAN; r->T_up_rt_valid=0;
            }

            double gl[3]={0,0,0};
            if (base_cs->wind_speed>0.0) {
                surface_direct_sunglint_rho(mu_air[iv],rt_raa_to_surface_view_phi(grid->raa_deg[ir]),
                    mus,0.0,base_cs->wind_speed,base_cs->sigma_type,nw,tau,
                    base_cs->q_convention,gl);
                if (opts->pssa && mus>0.0 && tau>0.0) {
                    gl[0]*=pssa_glint_scale; gl[1]*=pssa_glint_scale; gl[2]*=pssa_glint_scale;
                }
            }
            r->rho_glint_direct_I=gl[0]; r->rho_glint_direct_Q=gl[1]; r->rho_glint_direct_U=-gl[2];
            r->rho_I_glint=r->rho_I+r->rho_glint_direct_I;
            r->rho_Q_glint=r->rho_Q+r->rho_glint_direct_Q;
            r->rho_U_glint=r->rho_U+r->rho_glint_direct_U;
            r->T_sg_up_dir=r->rho_glint_direct_I*mus;
            r->T_total_up_dir=r->rho_I*mus;
        }
    }
    rc=0;
done:
    free(wupI); free(wupQ); free(wupU); free(wmu); free(ww);
    free(wvI); free(wvQ); free(wvU); free(waI); free(waQ); free(waU);
    free(mu_air); free(exact); free(ex1I); free(ex1Q); free(ex1U);
    free(ex2I); free(ex2Q); free(ex2U); free(exact_tdiff);
    ocrt_ocean_water_replay_free(&replay);
    return rc;
}

