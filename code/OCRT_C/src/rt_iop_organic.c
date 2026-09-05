#define _USE_MATH_DEFINES
#include "rt_iop_organic.h"
#include "shared/mie_io.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ORGANIC_MAX_STAR 1024
/* Canonical detritus Mie files contain 771 phase-wavelength columns.
 * Keep explicit headroom and do not couple Chl execution to EAP tables.
 * DOC-REF: OCRT_EAP_RUNTIME_DISABLED_AND_WINDOWS_RUN_POLICY_2026-08-22. */
#define ORGANIC_MAX_PHW  1024
#define ORGANIC_PATH_MAX 1024

typedef struct {
    int n_star;
    double wl_star[ORGANIC_MAX_STAR];
    double astar[ORGANIC_MAX_STAR];
    double bstar[ORGANIC_MAX_STAR];

    int n_phw;
    double wl_phw[ORGANIC_MAX_PHW];
    double bbratio[ORGANIC_MAX_PHW];

    char mie_path[ORGANIC_PATH_MAX];
    int loaded;
} organic_table_t;

static organic_table_t g_phyto[ORGANIC_PHYTO_GROUP_COUNT];
static organic_table_t g_detritus;
static int g_phyto_abs_n = 0;
static double g_phyto_abs_wl[ORGANIC_MAX_STAR];
static double g_phyto_abs_star[ORGANIC_MAX_STAR];
static int g_ready = 0;

/* 이름표와 파일표는 순서가 enum 과 같아야 한다.  --ocrt-phyto-group 인자는
 * 이름표를 그대로 받는다. */
static const char *const k_group_name[ORGANIC_PHYTO_GROUP_COUNT] = {
    "pico", "nano", "micro",
    "eap_diatoms_pennate", "eap_chlorophytes", "eap_diatoms_centric",
    "eap_cryptophytes", "eap_cyano_blue", "eap_cyano_red",
    "eap_dinoflagellates", "eap_eustigmatophytes", "eap_hapto_pavlovaceae",
    "eap_pelagophytes", "eap_prasinophytes", "eap_prochlorococcus",
    "eap_hapto_prymnesiaceae", "eap_raphidophytes", "eap_rhodophytes",
    "eap_synechococcus", "eap_microcystis"
};
/* EAP phase tables remain archival/generator products only.  They are not
 * runtime dependencies while project policy keeps EAP scattering disabled. */
static const char *const k_detritus_file = "Detritus_Stramski2001.mie";

/* 대형 여부.  1 이면 지금은 쓸 수 없다.  근거는 헤더 주석에 적었다.
 * 순서는 위 두 표와 같아야 한다.
 *
 * 측정값(412 nm, L=200 되살림, 30~150도 음수 개수 / |beta200/beta2|):
 *   pico            0개 / 6.0e-05      prochlorococcus 0.5um  0개 / 5.4e-05
 *   synechococcus 1.2um 0개 / 4.6e-05  pelagophytes 3um      42개 / 6.1e-02
 *   prymnesiaceae 4um  95개 / 3.0e-01  microcystis 5um       71개 / 1.0e+00
 *   centric 6um       118개 / 1.4e+00  dinoflagellates 24um 121개 / 3.4e+01
 * 경계는 1.2 um 와 3 um 사이에 있다. */
static const int k_group_large[ORGANIC_PHYTO_GROUP_COUNT] = {
    0, 1, 1,          /* pico(가능), nano, micro */
    1, 1, 1, 1, 1, 1, /* pennate, chlorophytes, centric, cryptophytes, cyano_blue, cyano_red */
    1, 1, 1, 1, 1,    /* dinoflagellates, eustigmatophytes, pavlovaceae, pelagophytes, prasinophytes */
    0,                /* prochlorococcus 0.5 um — 사용 가능 */
    1, 1, 1,          /* prymnesiaceae, raphidophytes, rhodophytes */
    0,                /* synechococcus 1.2 um — 사용 가능 */
    1                 /* microcystis */
};

static int group_valid(organic_phyto_group_t group)
{
    return group >= ORGANIC_PHYTO_PICO && group < ORGANIC_PHYTO_GROUP_COUNT;
}

static double interpolate(const double *x, const double *y, int n, double xq)
{
    if (!x || !y || n <= 0 || !isfinite(xq)) return 0.0;
    if (xq <= x[0]) return y[0];
    if (xq >= x[n - 1]) return y[n - 1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (x[mid] <= xq) lo = mid;
        else hi = mid;
    }
    double dx = x[hi] - x[lo];
    if (!(dx > 0.0)) return y[lo];
    double t = (xq - x[lo]) / dx;
    return y[lo] + t * (y[hi] - y[lo]);
}

static int table_monotone(const double *x, const double *y, int n, int positive)
{
    if (!x || !y || n < 2) return 0;
    for (int i = 0; i < n; ++i) {
        if (!isfinite(x[i]) || !isfinite(y[i]) || !(x[i] > 0.0)) return 0;
        if (i > 0 && !(x[i] > x[i - 1])) return 0;
        if (positive ? !(y[i] > 0.0) : (y[i] < 0.0)) return 0;
    }
    return 1;
}

static void sort_pair(double *x, double *y, int n)
{
    for (int i = 1; i < n; ++i) {
        int j = i;
        while (j > 0 && x[j] < x[j - 1]) {
            double t = x[j]; x[j] = x[j - 1]; x[j - 1] = t;
            t = y[j]; y[j] = y[j - 1]; y[j - 1] = t;
            --j;
        }
    }
}

static int summarize_mie(const char *path, int need_star, organic_table_t *out)
{
    if (!path || !out) return -1;
    mie_data_t mie = {0};
    if (read_mie_file(path, &mie) != 0) return -2;

    memset(out, 0, sizeof(*out));
    if (snprintf(out->mie_path, sizeof out->mie_path, "%s", path) >=
        (int)sizeof out->mie_path) {
        mie_data_free(&mie);
        return -3;
    }

    if (need_star) {
        if (mie.n_wl < 2 || mie.n_wl > ORGANIC_MAX_STAR) {
            mie_data_free(&mie);
            return -4;
        }
        out->n_star = mie.n_wl;
        for (int i = 0; i < mie.n_wl; ++i) {
            double ext = SPECTRAL_AT(&mie, i, 4);
            double sca = SPECTRAL_AT(&mie, i, 5);
            if (!isfinite(ext) || !isfinite(sca) || ext < 0.0 || sca <= 0.0 ||
                sca > ext + 1.0e-10) {
                mie_data_free(&mie);
                return -5;
            }
            out->wl_star[i] = 1000.0 * mie.wavelengths[i];
            out->astar[i] = fmax(0.0, ext - sca);
            out->bstar[i] = sca;
        }
        if (!table_monotone(out->wl_star, out->astar, out->n_star, 0) ||
            !table_monotone(out->wl_star, out->bstar, out->n_star, 1)) {
            mie_data_free(&mie);
            return -6;
        }
    }

    if (mie.n_phase_wl < 2 || mie.n_phase_wl > ORGANIC_MAX_PHW || mie.n_ang < 3) {
        mie_data_free(&mie);
        return -7;
    }
    out->n_phw = mie.n_phase_wl;
    for (int iw = 0; iw < mie.n_phase_wl; ++iw) {
        double full = 0.0, back = 0.0;
        for (int ia = 0; ia + 1 < mie.n_ang; ++ia) {
            double t0 = mie.angles[ia] * M_PI / 180.0;
            double t1 = mie.angles[ia + 1] * M_PI / 180.0;
            double dt = fabs(t1 - t0);
            double p0 = BLOCK_AT(mie.P11, mie.n_phase_wl, ia, iw);
            double p1 = BLOCK_AT(mie.P11, mie.n_phase_wl, ia + 1, iw);
            double seg = 0.5 * (p0 * sin(t0) + p1 * sin(t1)) * dt;
            if (!isfinite(seg)) {
                mie_data_free(&mie);
                return -8;
            }
            full += seg;
            if (0.5 * (t0 + t1) >= 0.5 * M_PI) back += seg;
        }
        double ratio = (full > 0.0) ? back / full : -1.0;
        if (!(ratio > 0.0 && ratio < 0.5) || !isfinite(ratio)) {
            mie_data_free(&mie);
            return -9;
        }
        out->wl_phw[iw] = 1000.0 * mie.phase_wavelengths[iw];
        out->bbratio[iw] = ratio;
    }
    sort_pair(out->wl_phw, out->bbratio, out->n_phw);
    if (!table_monotone(out->wl_phw, out->bbratio, out->n_phw, 1)) {
        mie_data_free(&mie);
        return -10;
    }

    out->loaded = 1;
    mie_data_free(&mie);
    return 0;
}

void rt_iop_organic_free(void)
{
    memset(g_phyto, 0, sizeof g_phyto);
    memset(&g_detritus, 0, sizeof g_detritus);
    g_phyto_abs_n = 0;
    memset(g_phyto_abs_wl, 0, sizeof g_phyto_abs_wl);
    memset(g_phyto_abs_star, 0, sizeof g_phyto_abs_star);
    g_ready = 0;
}

static int organic_load_phyto_absorption(const char *data_dir)
{
    char path[ORGANIC_PATH_MAX];
    int n = snprintf(path, sizeof path, "%s/%s", data_dir,
                     "phyto_absorption_default.csv");
    if (n < 0 || (size_t)n >= sizeof path) return -1;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[organic] failed to load absorption table %s\n", path);
        return -1;
    }
    char line[512];
    int count = 0;
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        double wl = 0.0, a = 0.0;
        if (sscanf(line, " %lf , %lf", &wl, &a) != 2) continue;
        if (count >= ORGANIC_MAX_STAR || !(wl > 0.0) || !(a >= 0.0) ||
            !isfinite(wl) || !isfinite(a) ||
            (count > 0 && !(wl > g_phyto_abs_wl[count - 1]))) {
            fclose(fp); return -2;
        }
        g_phyto_abs_wl[count] = wl;
        g_phyto_abs_star[count] = a;
        ++count;
    }
    fclose(fp);
    if (count < 2) return -3;
    /* Production contract: the data file itself must cover 330--1100 nm.
     * The approved longwave closure is encoded as explicit zero rows from
     * 800 nm onward; no wavelength-specific runtime branch is permitted. */
    if (!rt_spectral_grid_covers(g_phyto_abs_wl[0],
                                 g_phyto_abs_wl[count - 1]))
        return -4;
    int have_zero_anchor = 0;
    for (int i = 0; i < count; ++i) {
        if (fabs(g_phyto_abs_wl[i] - ORGANIC_PHYTO_ABS_ZERO_START_NM) < 1e-9)
            have_zero_anchor = 1;
        if (g_phyto_abs_wl[i] >= ORGANIC_PHYTO_ABS_ZERO_START_NM - 1e-9 &&
            g_phyto_abs_star[i] != 0.0)
            return -5;
    }
    if (!have_zero_anchor) return -6;
    g_phyto_abs_n = count;
    return 0;
}

static int organic_load_detritus_table(const char *data_dir)
{
    char path[ORGANIC_PATH_MAX];
    int n = snprintf(path, sizeof path, "%s/%s", data_dir, k_detritus_file);
    if (n < 0 || (size_t)n >= sizeof path ||
        summarize_mie(path, 0, &g_detritus) != 0) {
        fprintf(stderr, "[organic] failed to load detritus table %s\n", path);
        return -1;
    }
    if (g_detritus.n_phw < 2 ||
        !(g_detritus.wl_phw[g_detritus.n_phw - 1] > g_detritus.wl_phw[0])) {
        fprintf(stderr, "[organic] invalid detritus phase wavelength axis in %s\n", path);
        return -2;
    }
    /* Do not silently extrapolate the active detritus phase. Canonical
     * production data cover 330--1100 nm with 771 explicit phase-wavelength
     * columns; callers fail loudly outside the loaded file range. */
    return 0;
}

int rt_iop_organic_init_selected(const char *data_dir,
                                 organic_phyto_group_t selected_group,
                                 int load_all_phyto)
{
    /* EAP_RUNTIME_DISABLED: selected_group/load_all_phyto are retained in the
     * ABI only.  Chl production loads the common absorption table and the
     * Stramski detritus phase; no EAP .mie path is opened or required. */
    (void)selected_group;
    (void)load_all_phyto;
    if (!data_dir || !data_dir[0]) data_dir = ".";
    rt_iop_organic_free();

    if (organic_load_phyto_absorption(data_dir) != 0) {
        rt_iop_organic_free();
        return -1;
    }
    if (organic_load_detritus_table(data_dir) != 0) {
        rt_iop_organic_free();
        return -2;
    }
    g_ready = 1;
    return 0;
}

int rt_iop_organic_init(const char *data_dir)
{
    return rt_iop_organic_init_selected(
        data_dir, ORGANIC_PHYTO_MICRO, 0);
}

int rt_iop_organic_ready(void) { return g_ready; }

int rt_iop_organic_group_parse(const char *name, organic_phyto_group_t *out)
{
    if (!name || !name[0] || !out) return -1;
    char norm[96];
    int n = 0;
    for (const unsigned char *p = (const unsigned char *)name;
         *p && n < (int)sizeof norm - 1; ++p) {
        if (isalnum(*p)) norm[n++] = (char)tolower(*p);
    }
    norm[n] = 0;

    /* 1) 이름표와 정확히 맞는 항목을 먼저 찾는다.  이 단계가 앞에 와야
     *    한다.  아래 별칭 규칙이 부분 문자열로 판정하므로, 예를 들어
     *    "eap_synechococcus" 가 pico 로, "eap_diatoms_pennate" 가 micro 로
     *    잘못 빨려 들어가는 것을 막는다. */
    for (int g = 0; g < ORGANIC_PHYTO_GROUP_COUNT; ++g) {
        char tn[96];
        int m = 0;
        for (const unsigned char *p = (const unsigned char *)k_group_name[g];
             *p && m < (int)sizeof tn - 1; ++p) {
            if (isalnum(*p)) tn[m++] = (char)tolower(*p);
        }
        tn[m] = 0;
        if (!strcmp(norm, tn)) { *out = (organic_phyto_group_t)g; return 0; }
    }

    /* 2) 기존 세 그룹의 별칭.  뒤로 물린 이유는 위 주석과 같다. */
    if (strstr(norm, "synechococcus")) { *out = ORGANIC_PHYTO_PICO;  return 0; }
    if (strstr(norm, "haptophyte"))    { *out = ORGANIC_PHYTO_NANO;  return 0; }
    if (strstr(norm, "diatom"))        { *out = ORGANIC_PHYTO_MICRO; return 0; }
    return -1;
}

organic_phyto_group_t rt_iop_organic_group_from_name(const char *name)
{
    organic_phyto_group_t group = ORGANIC_PHYTO_MICRO;
    (void)rt_iop_organic_group_parse(name, &group);
    return group;
}

const char *rt_iop_organic_group_name(organic_phyto_group_t group)
{
    if (!group_valid(group)) group = ORGANIC_PHYTO_MICRO;
    return k_group_name[group];
}

int rt_iop_organic_group_is_large(organic_phyto_group_t group)
{
    if (!group_valid(group)) return 0;
    return k_group_large[group];
}

const char *rt_iop_organic_phyto_phase_path(organic_phyto_group_t group)
{
    /* EAP_RUNTIME_DISABLED: never expose a phytoplankton phase path to RT. */
    (void)group;
    return NULL;
}

const char *rt_iop_organic_detritus_phase_path(void)
{
    return (g_ready && g_detritus.loaded) ? g_detritus.mie_path : NULL;
}

int rt_iop_organic_detritus_wavelength_supported(double lambda_nm)
{
    return g_ready && g_detritus.loaded && g_detritus.n_phw >= 2 &&
           isfinite(lambda_nm) &&
           lambda_nm >= g_detritus.wl_phw[0] - RT_SPECTRAL_TOL_NM &&
           lambda_nm <= g_detritus.wl_phw[g_detritus.n_phw - 1] + RT_SPECTRAL_TOL_NM;
}

double rt_iop_organic_detritus_wavelength_min_nm(void)
{
    return (g_ready && g_detritus.loaded && g_detritus.n_phw > 0)
               ? g_detritus.wl_phw[0] : NAN;
}

double rt_iop_organic_detritus_wavelength_max_nm(void)
{
    return (g_ready && g_detritus.loaded && g_detritus.n_phw > 0)
               ? g_detritus.wl_phw[g_detritus.n_phw - 1] : NAN;
}

int rt_iop_organic_wavelength_supported(double lambda_nm)
{
    return rt_spectral_wavelength_supported(lambda_nm);
}

double rt_iop_huot_bbp(double lambda_nm, double chl_mg_m3)
{
    if (!(chl_mg_m3 > 0.0) || !isfinite(lambda_nm)) return 0.0;
    double alpha1 = 2.267e-3 - 5.058e-6 * (lambda_nm - 550.0);
    double beta1 = 0.565 + 0.000486 * (lambda_nm - 550.0);
    if (!(alpha1 > 0.0) || !isfinite(beta1)) return 0.0;
    return alpha1 * pow(chl_mg_m3, beta1);
}

double rt_iop_fph_fraction(double chl_mg_m3)
{
    if (!(chl_mg_m3 > 0.0)) return 0.020;
    return 0.035 + 0.015 * tanh(log10(chl_mg_m3));
}

double rt_iop_detritus_abs_shape(double lambda_nm, double slope_nm_inv)
{
    double slope = (slope_nm_inv > 0.0) ? slope_nm_inv : ORGANIC_DETRITUS_SLOPE_DEFAULT;
    return exp(-slope * (lambda_nm - 440.0));
}

int rt_iop_eap_phyto_absorption_eval(
    double lambda_nm, double chl_mg_m3, organic_phyto_group_t group,
    double *a_out_m_inv)
{
    if (!a_out_m_inv || !isfinite(lambda_nm) || !isfinite(chl_mg_m3)) return -1;
    *a_out_m_inv = 0.0;
    if (!(chl_mg_m3 > 0.0)) return 0;
    if (!g_ready) return -2;
    if (!rt_iop_organic_wavelength_supported(lambda_nm)) return -3;
    (void)group; /* species-specific production absorption is not enabled */
    if (g_phyto_abs_n < 2) return -4;
    /* The 800--1100 nm zero closure is represented directly in the table;
     * runtime evaluation is ordinary in-range interpolation only. */
    const double astar = interpolate(g_phyto_abs_wl, g_phyto_abs_star,
                                     g_phyto_abs_n, lambda_nm);
    if (!(astar >= 0.0) || !isfinite(astar)) return -4;
    *a_out_m_inv = astar * chl_mg_m3;
    return isfinite(*a_out_m_inv) ? 0 : -5;
}

int rt_iop_eap_phyto_eval_with_phase_ratios(
    double lambda_nm, double chl_mg_m3, organic_phyto_group_t group,
    double bb_ratio_lambda, double bb_ratio_550, rt_iop_t *out)
{
    (void)lambda_nm; (void)chl_mg_m3; (void)group;
    (void)bb_ratio_lambda; (void)bb_ratio_550;
    if (out) out->a = out->b = out->bb = 0.0;
    return -8; /* EAP_RUNTIME_DISABLED */
}

int rt_iop_eap_phyto_eval(double lambda_nm, double chl_mg_m3,
                          organic_phyto_group_t group, rt_iop_t *out)
{
    (void)lambda_nm; (void)chl_mg_m3; (void)group;
    if (out) out->a = out->b = out->bb = 0.0;
    return -8; /* EAP_RUNTIME_DISABLED */
}

int rt_iop_detritus_eval_with_phase_ratios(
    double lambda_nm, double chl_mg_m3,
    double a_d440_m_inv, double slope_nm_inv,
    double bb_ratio_lambda, double bb_ratio_550, rt_iop_t *out)
{
    if (!out || !isfinite(lambda_nm) || !isfinite(chl_mg_m3) ||
        !isfinite(a_d440_m_inv) || !isfinite(slope_nm_inv)) return -1;
    out->a = out->b = out->bb = 0.0;
    if (!(chl_mg_m3 > 0.0)) return 0;
    if (!g_ready) return -2;
    if (!rt_iop_organic_wavelength_supported(lambda_nm)) return -3;
    if (!rt_iop_organic_detritus_wavelength_supported(lambda_nm)) return -7;
    if (a_d440_m_inv < 0.0) return -4;

    if (!(bb_ratio_lambda > 0.0 && bb_ratio_lambda < 0.5) ||
        !(bb_ratio_550 > 0.0 && bb_ratio_550 < 0.5) ||
        !isfinite(bb_ratio_lambda) || !isfinite(bb_ratio_550))
        return -5;

    out->a = a_d440_m_inv * rt_iop_detritus_abs_shape(lambda_nm, slope_nm_inv);
    double shape = bb_ratio_lambda / bb_ratio_550;
    out->bb = (1.0 - rt_iop_fph_fraction(chl_mg_m3)) *
              rt_iop_huot_bbp(550.0, chl_mg_m3) * shape;
    out->b = out->bb / bb_ratio_lambda;
    return (isfinite(out->a) && isfinite(out->b) && isfinite(out->bb)) ? 0 : -6;
}

int rt_iop_detritus_eval(double lambda_nm, double chl_mg_m3,
                         double a_d440_m_inv, double slope_nm_inv,
                         rt_iop_t *out)
{
    if (!g_ready) {
        if (out) out->a = out->b = out->bb = 0.0;
        return -2;
    }
    if (!rt_iop_organic_detritus_wavelength_supported(lambda_nm)) {
        if (out) out->a = out->b = out->bb = 0.0;
        return -7;
    }
    const double ratio = interpolate(g_detritus.wl_phw, g_detritus.bbratio,
                                     g_detritus.n_phw, lambda_nm);
    const double ratio550 = interpolate(g_detritus.wl_phw, g_detritus.bbratio,
                                        g_detritus.n_phw, 550.0);
    return rt_iop_detritus_eval_with_phase_ratios(
        lambda_nm, chl_mg_m3, a_d440_m_inv, slope_nm_inv,
        ratio, ratio550, out);
}
