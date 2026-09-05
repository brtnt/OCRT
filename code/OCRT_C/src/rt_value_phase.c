#include "rt_value_phase.h"

#include "rt_phase_fr631.h"
#include "shared/mat3.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

typedef struct {
    double theta_deg;
    double p11;
    double p12;
    double p33;
} phase_row_t;

static int row_cmp_theta(const void *aa, const void *bb)
{
    const phase_row_t *a = (const phase_row_t *)aa;
    const phase_row_t *b = (const phase_row_t *)bb;
    return (a->theta_deg > b->theta_deg) - (a->theta_deg < b->theta_deg);
}

static double phase_lerp(double a, double b, double w)
{
#if defined(FP_FAST_FMA) || defined(__FMA__)
    return fma(w, b - a, a);
#else
    return a + w * (b - a);
#endif
}

static void phase_eval_sorted_rows(const phase_row_t *rows, int n,
                                   double theta_deg,
                                   double *p11, double *p12, double *p33)
{
    int lo = 0, hi = n - 1;
    if (theta_deg <= rows[0].theta_deg) {
        lo = 0;
        theta_deg = rows[0].theta_deg;
    } else if (theta_deg >= rows[n - 1].theta_deg) {
        lo = n - 2;
        theta_deg = rows[n - 1].theta_deg;
    } else {
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (rows[mid].theta_deg <= theta_deg) lo = mid;
            else hi = mid;
        }
    }
    const double h = rows[lo + 1].theta_deg - rows[lo].theta_deg;
    const double w = (h > 0.0) ?
        (theta_deg - rows[lo].theta_deg) / h : 0.0;
    if (p11) *p11 = phase_lerp(rows[lo].p11, rows[lo + 1].p11, w);
    if (p12) *p12 = phase_lerp(rows[lo].p12, rows[lo + 1].p12, w);
    if (p33) *p33 = phase_lerp(rows[lo].p33, rows[lo + 1].p33, w);
}

/* Exact integral of a piecewise-linear P(theta) times sin(theta), where theta
 * is in radians and p0/p1 are the endpoint values. */
static double segment_int_p_sin(double t0, double t1, double p0, double p1)
{
    const double h = t1 - t0;
    if (!(h > 0.0)) return 0.0;
    const double slope = (p1 - p0) / h;
    const double f0 = -p0 * cos(t0) + slope * sin(t0);
    const double f1 = -p1 * cos(t1) + slope * sin(t1);
    return f1 - f0;
}

/* Exact integral of a piecewise-linear P(theta) times cos(theta)sin(theta). */
static double segment_int_p_cos_sin(double t0, double t1,
                                    double p0, double p1)
{
    const double h = t1 - t0;
    if (!(h > 0.0)) return 0.0;
    const double slope = (p1 - p0) / h;
    const double a = p0 - slope * t0;
#define FG(T) (0.5 * a * sin(T) * sin(T) + \
               slope * (0.5 * (T) * sin(T) * sin(T) - 0.25 * (T) + 0.125 * sin(2.0 * (T))))
    const double result = FG(t1) - FG(t0);
#undef FG
    return result;
}

static void table_integrals(const double *theta_deg, const double *p11, int n,
                            double *norm, double *g_asym, double *bb_b_ratio)
{
    double total = 0.0;
    double back = 0.0;
    double gnum = 0.0;
    const double half_pi = 0.5 * M_PI;

    for (int i = 0; i < n - 1; ++i) {
        const double t0 = theta_deg[i] * DEG2RAD;
        const double t1 = theta_deg[i + 1] * DEG2RAD;
        const double p0 = p11[i];
        const double p1 = p11[i + 1];
        const double seg = segment_int_p_sin(t0, t1, p0, p1);
        total += seg;
        gnum += segment_int_p_cos_sin(t0, t1, p0, p1);

        if (t0 >= half_pi) {
            back += seg;
        } else if (t1 > half_pi) {
            const double w = (half_pi - t0) / (t1 - t0);
            const double pm = phase_lerp(p0, p1, w);
            back += segment_int_p_sin(half_pi, t1, pm, p1);
        }
    }

    if (norm) *norm = 0.5 * total;
    if (g_asym) *g_asym = (total > 0.0) ? (gnum / total) : 0.0;
    if (bb_b_ratio) *bb_b_ratio = (total > 0.0) ? (back / total) : 0.0;
}

void rt_value_phase_interp_free(rt_value_phase_interp_t *t)
{
    if (!t) return;
    free(t->theta_deg);
    free(t->p11);
    free(t->p12);
    free(t->p33);
    memset(t, 0, sizeof *t);
}

static int value_phase_copy(const rt_value_phase_interp_t *src,
                            rt_value_phase_interp_t *out)
{
    if (!src || !out || src->n < 2 || !src->theta_deg || !src->p11 ||
        !src->p12 || !src->p33) return -1;
    memset(out, 0, sizeof *out);
    const size_t bytes = (size_t)src->n * sizeof(double);
    out->theta_deg = (double *)malloc(bytes);
    out->p11 = (double *)malloc(bytes);
    out->p12 = (double *)malloc(bytes);
    out->p33 = (double *)malloc(bytes);
    if (!out->theta_deg || !out->p11 || !out->p12 || !out->p33) {
        rt_value_phase_interp_free(out);
        return -2;
    }
    memcpy(out->theta_deg, src->theta_deg, bytes);
    memcpy(out->p11, src->p11, bytes);
    memcpy(out->p12, src->p12, bytes);
    memcpy(out->p33, src->p33, bytes);
    out->n = src->n;
    out->is_fr631 = src->is_fr631;
    out->norm_before = src->norm_before;
    out->g_asym = src->g_asym;
    out->bb_b_ratio = src->bb_b_ratio;
    return 0;
}

int rt_value_phase_interp_loglinear_truncate(
    const rt_value_phase_interp_t *src,
    double mu1, double mu2, double threshold_A,
    rt_value_phase_interp_t *out,
    double *A_out)
{
    if (A_out) *A_out = 0.0;
    if (!src || !out || src->n < 2 || !src->theta_deg || !src->p11 ||
        !src->p12 || !src->p33 || !isfinite(mu1) || !isfinite(mu2) ||
        !isfinite(threshold_A) || !(mu1 > -1.0 && mu1 < mu2 && mu2 < 1.0) ||
        threshold_A < 0.0) return -1;
    memset(out, 0, sizeof *out);

    const double theta1_deg = acos(mu1) * RAD2DEG;
    const double theta2_deg = acos(mu2) * RAD2DEG;
    if (!(theta1_deg > theta2_deg && theta2_deg > 0.0))
        return value_phase_copy(src, out);

    double p11_t1 = 0.0, p11_t2 = 0.0;
    rt_value_phase_interp_eval_theta(src, theta1_deg, &p11_t1, NULL, NULL);
    rt_value_phase_interp_eval_theta(src, theta2_deg, &p11_t2, NULL, NULL);
    if (!(p11_t1 > 0.0) || !(p11_t2 > 0.0) ||
        !isfinite(p11_t1) || !isfinite(p11_t2))
        return value_phase_copy(src, out);

    const double theta1_rad = theta1_deg * DEG2RAD;
    const double theta2_rad = theta2_deg * DEG2RAD;
    const double slope = (log10(p11_t2) - log10(p11_t1)) /
                         (theta2_rad - theta1_rad);
    const double intercept = log10(p11_t2);

    const int n = src->n;
    double *p11 = (double *)malloc((size_t)n * sizeof(double));
    double *p12 = (double *)malloc((size_t)n * sizeof(double));
    double *p33 = (double *)malloc((size_t)n * sizeof(double));
    if (!p11 || !p12 || !p33) {
        free(p11); free(p12); free(p33);
        return -2;
    }

    for (int i = 0; i < n; ++i) {
        const double theta = src->theta_deg[i];
        double q11 = src->p11[i];
        double ratio = 1.0;
        if (theta < theta2_deg) {
            const double logp = intercept + slope * (theta * DEG2RAD - theta2_rad);
            const double capped = pow(10.0, logp);
            if (!(capped >= 0.0) || !isfinite(capped)) {
                free(p11); free(p12); free(p33);
                return -3;
            }
            ratio = (q11 > 0.0) ? capped / q11 : 1.0;
            q11 = capped;
        }
        p11[i] = q11;
        p12[i] = src->p12[i] * ratio;
        p33[i] = src->p33[i] * ratio;
    }

    double chi0 = 0.0;
    table_integrals(src->theta_deg, p11, n, &chi0, NULL, NULL);
    double A = 2.0 * (1.0 - chi0);
    if (!isfinite(A) || !(chi0 > 0.0)) {
        free(p11); free(p12); free(p33);
        return -4;
    }

    /* Negative A means that the logarithmic continuation adds rather than
     * removes forward mass.  OSOAA's threshold contract treats this as no
     * truncation; do the same and never increase transport scattering. */
    if (!(A > 0.0) || A < threshold_A) {
        free(p11); free(p12); free(p33);
        if (A_out) *A_out = 0.0;
        return value_phase_copy(src, out);
    }
    if (!(A < 2.0)) {
        free(p11); free(p12); free(p33);
        return -5;
    }

    const int rc = rt_value_phase_interp_build(
        out, src->theta_deg, p11, p12, p33, n, 0);
    free(p11); free(p12); free(p33);
    if (rc != 0) return -6;
    if (A_out) *A_out = A;
    return 0;
}

static int sanitize_node(double *p11, double *p12, double *p33)
{
    const double scale = fmax(1.0, fabs(*p11));
    const double tiny = 2.0e-12 * scale;
    const double bound_tol = 2.0e-8 * scale;

    if (!isfinite(*p11) || !isfinite(*p12) || !isfinite(*p33)) return 0;
    if (*p11 < -tiny) return 0;
    if (*p11 < 0.0) *p11 = 0.0;

    if (fabs(*p12) > *p11 + bound_tol || fabs(*p33) > *p11 + bound_tol)
        return 0;
    if (fabs(*p12) > *p11) *p12 = copysign(*p11, *p12);
    if (fabs(*p33) > *p11) *p33 = copysign(*p11, *p33);
    return 1;
}

int rt_value_phase_interp_build(rt_value_phase_interp_t *out,
                                const double *theta_deg,
                                const double *p11,
                                const double *p12,
                                const double *p33,
                                int n,
                                int norm_n_mu)
{
    (void)norm_n_mu;
    if (!out || !theta_deg || !p11 || n < 2) return -1;
    memset(out, 0, sizeof *out);

    phase_row_t *rows = (phase_row_t *)malloc((size_t)n * sizeof(*rows));
    if (!rows) return -2;
    for (int i = 0; i < n; ++i) {
        rows[i].theta_deg = theta_deg[i];
        rows[i].p11 = p11[i];
        rows[i].p12 = p12 ? p12[i] : 0.0;
        rows[i].p33 = p33 ? p33[i] : p11[i];
        if (!isfinite(rows[i].theta_deg) ||
            !sanitize_node(&rows[i].p11, &rows[i].p12, &rows[i].p33)) {
            free(rows);
            return -3;
        }
    }
    qsort(rows, (size_t)n, sizeof(*rows), row_cmp_theta);

    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (m > 0 && fabs(rows[i].theta_deg - rows[m - 1].theta_deg) <= 1.0e-12) {
            rows[m - 1] = rows[i];
        } else {
            rows[m++] = rows[i];
        }
    }
    if (m < 2 || rows[0].theta_deg > 1.0e-7 ||
        rows[m - 1].theta_deg < 180.0 - 1.0e-7) {
        free(rows);
        return -4;
    }
    rows[0].theta_deg = 0.0;
    rows[m - 1].theta_deg = 180.0;
    for (int i = 1; i < m; ++i) {
        if (!(rows[i].theta_deg > rows[i - 1].theta_deg)) {
            free(rows);
            return -4;
        }
    }

    /* Canonicalize a legacy table only when every original node is already a
     * node of FR631 (notably the 361-node 0.5-degree grid).  Adding nodes on
     * existing straight segments then preserves the theta-linear function
     * exactly; arbitrary non-nested grids remain on the generic fallback. */
    int nested_fr631 = 1;
    {
        double grid[RT_FR631_N_ANGLE];
        rt_fr631_fill_theta(grid);
        for (int i = 0; i < m; ++i) {
            const rt_fr631_bracket_t b = rt_fr631_bracket(rows[i].theta_deg);
            const int i0 = (int)b.lower;
            const int i1 = i0 + 1;
            if (fabs(rows[i].theta_deg - grid[i0]) > 5.0e-7 &&
                fabs(rows[i].theta_deg - grid[i1]) > 5.0e-7) {
                nested_fr631 = 0;
                break;
            }
        }
    }
    const int out_n = (nested_fr631 ? RT_FR631_N_ANGLE : m);
    out->n = out_n;
    out->theta_deg = (double *)malloc((size_t)out_n * sizeof(double));
    out->p11 = (double *)malloc((size_t)out_n * sizeof(double));
    out->p12 = (double *)malloc((size_t)out_n * sizeof(double));
    out->p33 = (double *)malloc((size_t)out_n * sizeof(double));
    if (!out->theta_deg || !out->p11 || !out->p12 || !out->p33) {
        free(rows);
        rt_value_phase_interp_free(out);
        return -2;
    }
    if (nested_fr631) {
        rt_fr631_fill_theta(out->theta_deg);
        for (int i = 0; i < out_n; ++i) {
            phase_eval_sorted_rows(rows, m, out->theta_deg[i],
                                   &out->p11[i], &out->p12[i], &out->p33[i]);
        }
    } else {
        for (int i = 0; i < m; ++i) {
            out->theta_deg[i] = rows[i].theta_deg;
            out->p11[i] = rows[i].p11;
            out->p12[i] = rows[i].p12;
            out->p33[i] = rows[i].p33;
        }
    }
    free(rows);

    int descending = 0;
    out->is_fr631 = rt_fr631_validate_theta(out->theta_deg, (size_t)out_n,
                                            5.0e-7, &descending) && !descending;

    double norm = 0.0, g = 0.0, bb = 0.0;
    table_integrals(out->theta_deg, out->p11, out_n, &norm, &g, &bb);
    if (!(norm > 0.0) || !isfinite(norm)) {
        rt_value_phase_interp_free(out);
        return -5;
    }
    out->norm_before = norm;
    const double inv = 1.0 / norm;
    for (int i = 0; i < out_n; ++i) {
        out->p11[i] *= inv;
        out->p12[i] *= inv;
        out->p33[i] *= inv;
    }
    table_integrals(out->theta_deg, out->p11, out_n, NULL, &g, &bb);
    out->g_asym = g;
    out->bb_b_ratio = bb;
    return 0;
}

void rt_value_phase_interp_eval_cached(const rt_value_phase_interp_t *t,
                                       int lower, double weight,
                                       double *p11, double *p12, double *p33)
{
    if (!t || t->n < 2) {
        if (p11) *p11 = 0.0;
        if (p12) *p12 = 0.0;
        if (p33) *p33 = 0.0;
        return;
    }
    if (lower < 0) { lower = 0; weight = 0.0; }
    if (lower >= t->n - 1) { lower = t->n - 2; weight = 1.0; }
    if (weight < 0.0) weight = 0.0;
    if (weight > 1.0) weight = 1.0;
    if (p11) *p11 = phase_lerp(t->p11[lower], t->p11[lower + 1], weight);
    if (p12) *p12 = phase_lerp(t->p12[lower], t->p12[lower + 1], weight);
    if (p33) *p33 = phase_lerp(t->p33[lower], t->p33[lower + 1], weight);
}

void rt_value_phase_interp_eval_theta(const rt_value_phase_interp_t *t,
                                      double theta_deg,
                                      double *p11, double *p12, double *p33)
{
    if (!t || t->n < 2) {
        if (p11) *p11 = 0.0;
        if (p12) *p12 = 0.0;
        if (p33) *p33 = 0.0;
        return;
    }
    if (!(theta_deg > 0.0)) theta_deg = 0.0;
    else if (theta_deg >= 180.0) theta_deg = 180.0;

    if (t->is_fr631) {
        const rt_fr631_bracket_t b = rt_fr631_bracket(theta_deg);
        rt_value_phase_interp_eval_cached(t, (int)b.lower, b.weight,
                                          p11, p12, p33);
        return;
    }

    int lo = 0, hi = t->n - 1;
    if (theta_deg <= t->theta_deg[0]) {
        lo = 0;
        theta_deg = t->theta_deg[0];
    } else if (theta_deg >= t->theta_deg[t->n - 1]) {
        lo = t->n - 2;
        theta_deg = t->theta_deg[t->n - 1];
    } else {
        while (hi - lo > 1) {
            const int mid = (lo + hi) >> 1;
            if (t->theta_deg[mid] <= theta_deg) lo = mid;
            else hi = mid;
        }
    }
    const double h = t->theta_deg[lo + 1] - t->theta_deg[lo];
    const double w = (h > 0.0) ? (theta_deg - t->theta_deg[lo]) / h : 0.0;
    rt_value_phase_interp_eval_cached(t, lo, w, p11, p12, p33);
}

void rt_value_phase_interp_eval(const rt_value_phase_interp_t *t, double mu,
                                double *p11, double *p12, double *p33)
{
    if (mu > 1.0) mu = 1.0;
    else if (mu < -1.0) mu = -1.0;
    rt_value_phase_interp_eval_theta(t, acos(mu) * RAD2DEG, p11, p12, p33);
}

static size_t vidx(int m, int j, int k, int n_mu)
{
    const int kw = 2*n_mu + 1;
    return ((size_t)m*(size_t)(n_mu+1)+(size_t)j)*(size_t)kw + (size_t)(k+n_mu);
}

typedef struct {
    int valid;
    int m_count;
    int nphi;
    double *phi;
    double *cm;
    double *sm;
    unsigned long long build_count;
    unsigned long long hit_count;
} value_fourier_trig_cache_t;

static value_fourier_trig_cache_t g_value_trig_cache = {0};
#ifdef _OPENMP
#pragma omp threadprivate(g_value_trig_cache)
#endif

static int get_fourier_trig(int m_count, int nphi,
                            const double **phi_out,
                            const double **cos_out,
                            const double **sin_out)
{
    value_fourier_trig_cache_t *c=&g_value_trig_cache;
    if (c->valid && c->m_count==m_count && c->nphi==nphi &&
        c->phi && c->cm && c->sm) {
        ++c->hit_count;
        *phi_out=c->phi; *cos_out=c->cm; *sin_out=c->sm;
        return 0;
    }
    const size_t nmodephi=(size_t)m_count*(size_t)nphi;
    double *phi=(double*)realloc(c->phi,(size_t)nphi*sizeof(double));
    if (!phi) return -1;
    c->phi=phi;
    double *cm=(double*)realloc(c->cm,nmodephi*sizeof(double));
    if (!cm) return -1;
    c->cm=cm;
    double *sm=(double*)realloc(c->sm,nmodephi*sizeof(double));
    if (!sm) return -1;
    c->sm=sm;
    const double two_pi=2.0*M_PI;
    for (int q=0; q<nphi; ++q) {
        c->phi[q]=two_pi*((double)q+0.5)/(double)nphi;
        for (int m=0; m<m_count; ++m) {
            const size_t z=(size_t)m*(size_t)nphi+(size_t)q;
            c->cm[z]=cos((double)m*c->phi[q]);
            c->sm[z]=sin((double)m*c->phi[q]);
        }
    }
    c->valid=1; c->m_count=m_count; c->nphi=nphi;
    ++c->build_count;
    *phi_out=c->phi; *cos_out=c->cm; *sin_out=c->sm;
    return 0;
}

/* Geometry-only FR631 lookup map.  The nonzero GL/view ring is invariant to
 * wavelength, particle type, and the current solar beam in rm[0].  Cache the
 * lower FR631 index and common linear weight for that stable block so repeat
 * kernels avoid acos() and interval selection.  Solar row/column entries are
 * intentionally evaluated on the fly. */
typedef struct {
    int valid;
    int n_mu;
    int nphi;
    unsigned long long rm_hash;
    size_t count;
    uint16_t *lower;
    double *weight;
    unsigned long long build_count;
    unsigned long long hit_count;
} value_geom_cache_t;

static value_geom_cache_t g_value_geom_cache = {0};
#ifdef _OPENMP
#pragma omp threadprivate(g_value_geom_cache)
#endif

static unsigned long long geom_hash_bytes(unsigned long long h,
                                          const void *data, size_t n)
{
    const unsigned char *p = (const unsigned char *)data;
    for (size_t i = 0; i < n; ++i) {
        h ^= (unsigned long long)p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static unsigned long long value_geom_rm_hash(const rt_atm_t *atm)
{
    unsigned long long h = 1469598103934665603ULL;
    const int n_mu = atm->n_mu;
    h = geom_hash_bytes(h, &n_mu, sizeof n_mu);
    for (int k = -n_mu; k <= n_mu; ++k) {
        if (k == 0) continue;
        const double mu = atm->rm[k];
        h = geom_hash_bytes(h, &mu, sizeof mu);
    }
    return h;
}

static size_t value_geom_cache_cap_bytes(void)
{
    const char *env = getenv("OCRT_WATER_VALUE_GEOM_CACHE_MAX_MB");
    double mb = env ? atof(env) : 64.0;
    if (!isfinite(mb) || mb < 0.0) mb = 64.0;
    if (mb > 4096.0) mb = 4096.0;
    return (size_t)(mb * 1024.0 * 1024.0);
}

static size_t value_geom_index(int jin, int kout, int q, int n_mu, int nphi)
{
    const int jslot = jin - 1;
    const int kslot = (kout < 0) ? (kout + n_mu) : (n_mu - 1 + kout);
    return ((size_t)jslot * (size_t)(2 * n_mu) + (size_t)kslot) *
           (size_t)nphi + (size_t)q;
}

static int value_geom_cache_get(const rt_atm_t *atm, int nphi,
                                const uint16_t **lower_out,
                                const double **weight_out)
{
    if (lower_out) *lower_out = NULL;
    if (weight_out) *weight_out = NULL;
    if (!atm || atm->n_mu < 1 || nphi < 1 || !lower_out || !weight_out) return 0;

    const int n_mu = atm->n_mu;
    const size_t count = (size_t)n_mu * (size_t)(2 * n_mu) * (size_t)nphi;
    if (count > SIZE_MAX / (sizeof(uint16_t) + sizeof(double))) return 0;
    const size_t bytes = count * (sizeof(uint16_t) + sizeof(double));
    if (bytes > value_geom_cache_cap_bytes()) return 0;

    const unsigned long long hash = value_geom_rm_hash(atm);
    value_geom_cache_t *c = &g_value_geom_cache;
    if (c->valid && c->n_mu == n_mu && c->nphi == nphi &&
        c->rm_hash == hash && c->count == count && c->lower && c->weight) {
        ++c->hit_count;
        *lower_out = c->lower;
        *weight_out = c->weight;
        const char *trace = getenv("OCRT_VALUE_GEOM_CACHE_TRACE");
        if (trace && trace[0] && strcmp(trace, "0") != 0)
            fprintf(stderr, "[VALUE-GEOM-CACHE] hit nmu=%d nphi=%d bytes=%zu builds=%llu hits=%llu\n",
                    n_mu, nphi, bytes, c->build_count, c->hit_count);
        return 1;
    }

    c->valid = 0;
    uint16_t *new_lower = (uint16_t *)realloc(c->lower, count * sizeof(uint16_t));
    if (!new_lower) return 0;
    c->lower = new_lower;
    double *new_weight = (double *)realloc(c->weight, count * sizeof(double));
    if (!new_weight) return 0;
    c->weight = new_weight;

    const double two_pi = 2.0 * M_PI;
    for (int jin = 1; jin <= n_mu; ++jin) {
        const double mui = atm->rm[jin];
        const double si = sqrt(fmax(0.0, 1.0 - mui * mui));
        for (int kout = -n_mu; kout <= n_mu; ++kout) {
            if (kout == 0) continue;
            const double muo = atm->rm[kout];
            const double so = sqrt(fmax(0.0, 1.0 - muo * muo));
            for (int q = 0; q < nphi; ++q) {
                const double phi = two_pi * ((double)q + 0.5) / (double)nphi;
                double cT = mui * muo + si * so * cos(phi);
                if (cT > 1.0) cT = 1.0;
                else if (cT < -1.0) cT = -1.0;
                const rt_fr631_bracket_t b =
                    rt_fr631_bracket(acos(cT) * RAD2DEG);
                const size_t z = value_geom_index(jin, kout, q, n_mu, nphi);
                c->lower[z] = b.lower;
                c->weight[z] = b.weight;
            }
        }
    }
    c->valid = 1;
    c->n_mu = n_mu;
    c->nphi = nphi;
    c->rm_hash = hash;
    c->count = count;
    ++c->build_count;
    *lower_out = c->lower;
    *weight_out = c->weight;
    const char *trace = getenv("OCRT_VALUE_GEOM_CACHE_TRACE");
    if (trace && trace[0] && strcmp(trace, "0") != 0)
        fprintf(stderr, "[VALUE-GEOM-CACHE] build nmu=%d nphi=%d bytes=%zu builds=%llu hits=%llu\n",
                n_mu, nphi, bytes, c->build_count, c->hit_count);
    return 1;
}

int rt_value_phase_fourier_scalar_allm(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, double *pfm)
{
    if (!atm || !phase || phase->n < 2 || m_count < 1 || !pfm) return -1;
    if (nphi < 16) nphi = 16;
    if (nphi > 20000) nphi = 20000;
    const int n_mu = atm->n_mu;
    const uint16_t *geom_lower = NULL;
    const double *geom_weight = NULL;
    const int have_geom_cache = phase->is_fr631 &&
        value_geom_cache_get(atm, nphi, &geom_lower, &geom_weight);
    double *acc = (double *)calloc((size_t)m_count, sizeof(double));
    const double *phis = NULL, *cm = NULL, *sm = NULL;
    if (!acc || get_fourier_trig(m_count, nphi, &phis, &cm, &sm) != 0) {
        free(acc); return -2;
    }
    (void)sm;
    for (int jin = 0; jin <= n_mu; ++jin) {
        const double mui = atm->rm[jin];
        const double si = sqrt(fmax(0.0, 1.0 - mui*mui));
        for (int kout = -n_mu; kout <= n_mu; ++kout) {
            const double muo = atm->rm[kout];
            const double so = sqrt(fmax(0.0, 1.0 - muo*muo));
            memset(acc, 0, (size_t)m_count * sizeof(double));
            for (int q = 0; q < nphi; ++q) {
                double cT = mui*muo + si*so*cos(phis[q]);
                if (cT > 1.0) cT = 1.0;
                else if (cT < -1.0) cT = -1.0;
                double p = 0.0;
                if (have_geom_cache && jin > 0 && kout != 0) {
                    const size_t gz = value_geom_index(jin, kout, q, n_mu, nphi);
                    rt_value_phase_interp_eval_cached(
                        phase, (int)geom_lower[gz], geom_weight[gz],
                        &p, NULL, NULL);
                } else {
                    rt_value_phase_interp_eval(phase, cT, &p, NULL, NULL);
                }
                for (int m = 0; m < m_count; ++m)
                    acc[m] += p * cm[(size_t)m*(size_t)nphi + (size_t)q];
            }
            for (int m = 0; m < m_count; ++m)
                pfm[vidx(m, jin, kout, n_mu)] = acc[m] / (double)nphi;
        }
    }
    free(acc);
    return 0;
}

int rt_aerosol_value_phase_fourier_pol_rows(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, int j_begin, int j_end,
    double *pfm, double *gr, double *gt, double *arr,
    double *art, double *att)
{
    if (!atm || !phase || phase->n < 2 || m_count < 1 ||
        !pfm || !gr || !gt || !arr || !art || !att) return -1;
    if (nphi < 16) nphi = 16;
    if (nphi > 20000) nphi = 20000;
    const int n_mu = atm->n_mu;
    const uint16_t *geom_lower = NULL;
    const double *geom_weight = NULL;
    const int have_geom_cache = phase->is_fr631 &&
        value_geom_cache_get(atm, nphi, &geom_lower, &geom_weight);
    if (j_begin < 0) j_begin = 0;
    if (j_end > n_mu) j_end = n_mu;
    if (j_begin > j_end) return 0;

    double *cacc = (double *)calloc((size_t)m_count * 4u, sizeof(double));
    double *sacc = (double *)calloc((size_t)m_count * 2u, sizeof(double));
    const double *phis = NULL, *cm = NULL, *sm = NULL;
    if (!cacc || !sacc || get_fourier_trig(m_count, nphi, &phis, &cm, &sm) != 0) {
        free(cacc); free(sacc); return -2;
    }

    for (int jin = j_begin; jin <= j_end; ++jin) {
        const double mui = atm->rm[jin];
        const double si = sqrt(fmax(0.0, 1.0 - mui*mui));
        for (int kout = -n_mu; kout <= n_mu; ++kout) {
            const double muo = atm->rm[kout];
            const double so = sqrt(fmax(0.0, 1.0 - muo*muo));
            memset(cacc, 0, (size_t)m_count * 4u * sizeof(double));
            memset(sacc, 0, (size_t)m_count * 2u * sizeof(double));
            for (int q = 0; q < nphi; ++q) {
                const double phi = phis[q];
                const double cp = cos(phi), sp = sin(phi);
                double cT = mui*muo + si*so*cp;
                if (cT > 1.0) cT = 1.0;
                else if (cT < -1.0) cT = -1.0;
                const double sT = sqrt(fmax(0.0, 1.0 - cT*cT));
                double p11, p12, p33;
                if (have_geom_cache && jin > 0 && kout != 0) {
                    const size_t gz = value_geom_index(jin, kout, q, n_mu, nphi);
                    rt_value_phase_interp_eval_cached(
                        phase, (int)geom_lower[gz], geom_weight[gz],
                        &p11, &p12, &p33);
                } else {
                    rt_value_phase_interp_eval(phase, cT, &p11, &p12, &p33);
                }
                const double P[9] = {p11,p12,0.0,p12,p11,0.0,0.0,0.0,p33};
                double ci1,si1,ci2,si2,L1[9],L2[9],tmp[9],Z[9];
                compute_hovenier_rotation_cp(muo,mui,so,si,cT,sT,cp,sp,
                                           &ci1,&si1,&ci2,&si2);
                build_rotation_L(ci1,si1,L1);
                build_rotation_L(ci2,si2,L2);
                mat3_mul(P,L1,tmp);
                mat3_mul(L2,tmp,Z);
                for (int m = 0; m < m_count; ++m) {
                    const double c = cm[(size_t)m*(size_t)nphi + (size_t)q];
                    const double s = sm[(size_t)m*(size_t)nphi + (size_t)q];
                    cacc[4*m+0] += Z[0]*c;
                    cacc[4*m+1] += Z[3]*c;
                    cacc[4*m+2] += Z[4]*c;
                    cacc[4*m+3] += Z[8]*c;
                    sacc[2*m+0] += (-Z[6])*s;
                    sacc[2*m+1] += Z[5]*s;
                }
            }
            for (int m = 0; m < m_count; ++m) {
                const size_t z = vidx(m,jin,kout,n_mu);
                pfm[z] = cacc[4*m+0] / (double)nphi;
                gr[z]  = cacc[4*m+1] / (double)nphi;
                gt[z]  = sacc[2*m+0] / (double)nphi;
                arr[z] = cacc[4*m+2] / (double)nphi;
                art[z] = sacc[2*m+1] / (double)nphi;
                att[z] = cacc[4*m+3] / (double)nphi;
            }
        }
    }
    free(cacc); free(sacc);
    return 0;
}

int rt_aerosol_value_phase_fourier_pol_allm(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, double *pfm, double *gr, double *gt, double *arr,
    double *art, double *att)
{
    if (!atm) return -1;
    return rt_aerosol_value_phase_fourier_pol_rows(
        atm, m_count, phase, nphi, 0, atm->n_mu,
        pfm, gr, gt, arr, art, att);
}

/* Same algebra and quadrature order as the full builder, restricted to the
 * two beam-dependent strips. */
int rt_aerosol_value_phase_fourier_pol_solar_allm(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, double *pfm, double *gr, double *gt, double *arr,
    double *art, double *att)
{
    if (!atm || !phase || phase->n < 2 || m_count < 1 ||
        !pfm || !gr || !gt || !arr || !art || !att) return -1;
    if (nphi < 16) nphi = 16;
    if (nphi > 20000) nphi = 20000;
    const int n_mu = atm->n_mu;
    double *cacc = (double *)calloc((size_t)m_count * 4u, sizeof(double));
    double *sacc = (double *)calloc((size_t)m_count * 2u, sizeof(double));
    const double *phis = NULL, *cm = NULL, *sm = NULL;
    if (!cacc || !sacc || get_fourier_trig(m_count, nphi, &phis, &cm, &sm) != 0) {
        free(cacc); free(sacc); return -2;
    }

#define EVAL_PAIR(JIN,KOUT) do { \
        const int jin_ = (JIN), kout_ = (KOUT); \
        const double mui = atm->rm[jin_]; \
        const double si = sqrt(fmax(0.0, 1.0 - mui*mui)); \
        const double muo = atm->rm[kout_]; \
        const double so = sqrt(fmax(0.0, 1.0 - muo*muo)); \
        memset(cacc, 0, (size_t)m_count * 4u * sizeof(double)); \
        memset(sacc, 0, (size_t)m_count * 2u * sizeof(double)); \
        for (int q = 0; q < nphi; ++q) { \
            const double phi = phis[q]; \
            const double cp = cos(phi), sp = sin(phi); \
            double cT = mui*muo + si*so*cp; \
            if (cT > 1.0) cT = 1.0; else if (cT < -1.0) cT = -1.0; \
            const double sT = sqrt(fmax(0.0, 1.0 - cT*cT)); \
            double p11, p12, p33; \
            rt_value_phase_interp_eval(phase, cT, &p11, &p12, &p33); \
            const double P[9] = {p11,p12,0.0,p12,p11,0.0,0.0,0.0,p33}; \
            double ci1,si1,ci2,si2,L1[9],L2[9],tmp[9],Z[9]; \
            compute_hovenier_rotation_cp(muo,mui,so,si,cT,sT,cp,sp, \
                                       &ci1,&si1,&ci2,&si2); \
            build_rotation_L(ci1,si1,L1); build_rotation_L(ci2,si2,L2); \
            mat3_mul(P,L1,tmp); mat3_mul(L2,tmp,Z); \
            for (int mm = 0; mm < m_count; ++mm) { \
                const double c = cm[(size_t)mm*(size_t)nphi + (size_t)q]; \
                const double s = sm[(size_t)mm*(size_t)nphi + (size_t)q]; \
                cacc[4*mm+0] += Z[0]*c; \
                cacc[4*mm+1] += Z[3]*c; \
                cacc[4*mm+2] += Z[4]*c; \
                cacc[4*mm+3] += Z[8]*c; \
                sacc[2*mm+0] += (-Z[6])*s; \
                sacc[2*mm+1] += Z[5]*s; \
            } \
        } \
        for (int mm = 0; mm < m_count; ++mm) { \
            const size_t z = vidx(mm,jin_,kout_,n_mu); \
            pfm[z]=cacc[4*mm+0]/(double)nphi; \
            gr[z] =cacc[4*mm+1]/(double)nphi; \
            gt[z] =sacc[2*mm+0]/(double)nphi; \
            arr[z]=cacc[4*mm+2]/(double)nphi; \
            art[z]=sacc[2*mm+1]/(double)nphi; \
            att[z]=cacc[4*mm+3]/(double)nphi; \
        } \
    } while (0)

    for (int kout=-n_mu; kout<=n_mu; ++kout) EVAL_PAIR(0,kout);
    for (int jin=1; jin<=n_mu; ++jin) EVAL_PAIR(jin,0);
#undef EVAL_PAIR
    free(cacc); free(sacc);
    return 0;
}

int rt_aerosol_value_phase_fourier_pol_beamq_allm(
    const rt_atm_t *atm, int m_count, const rt_value_phase_interp_t *phase,
    int nphi, double *i_from_q, double *q_from_q, double *u_from_q)
{
    if (!atm || !phase || phase->n < 2 || m_count < 1 ||
        !i_from_q || !q_from_q || !u_from_q) return -1;
    if (nphi < 16) nphi = 16;
    if (nphi > 20000) nphi = 20000;

    const int n_mu = atm->n_mu;
    const int directions = 2 * n_mu + 1;
    double *c_iq = (double *)calloc((size_t)m_count, sizeof(double));
    double *c_qq = (double *)calloc((size_t)m_count, sizeof(double));
    double *s_uq = (double *)calloc((size_t)m_count, sizeof(double));
    const double *phis = NULL, *cm = NULL, *sm = NULL;
    if (!c_iq || !c_qq || !s_uq ||
        get_fourier_trig(m_count, nphi, &phis, &cm, &sm) != 0) {
        free(c_iq); free(c_qq); free(s_uq);
        return -2;
    }

    const double mui = atm->rm[0];
    const double si = sqrt(fmax(0.0, 1.0 - mui * mui));
    for (int kout = -n_mu; kout <= n_mu; ++kout) {
        const int ko = kout + n_mu;
        if (kout == 0) {
            for (int m = 0; m < m_count; ++m) {
                const size_t z = (size_t)m * (size_t)directions + (size_t)ko;
                i_from_q[z] = q_from_q[z] = u_from_q[z] = 0.0;
            }
            continue;
        }
        const double muo = atm->rm[kout];
        const double so = sqrt(fmax(0.0, 1.0 - muo * muo));
        memset(c_iq, 0, (size_t)m_count * sizeof(double));
        memset(c_qq, 0, (size_t)m_count * sizeof(double));
        memset(s_uq, 0, (size_t)m_count * sizeof(double));
        for (int q = 0; q < nphi; ++q) {
            const double phi = phis[q];
            const double cp = cos(phi), sp = sin(phi);
            double cT = mui * muo + si * so * cp;
            if (cT > 1.0) cT = 1.0;
            else if (cT < -1.0) cT = -1.0;
            const double sT = sqrt(fmax(0.0, 1.0 - cT * cT));
            double p11, p12, p33;
            rt_value_phase_interp_eval(phase, cT, &p11, &p12, &p33);
            const double P[9] = {p11,p12,0.0,p12,p11,0.0,0.0,0.0,p33};
            double ci1,si1,ci2,si2,L1[9],L2[9],tmp[9],Z[9];
            compute_hovenier_rotation(muo,mui,so,si,cT,sT,sp,
                                       &ci1,&si1,&ci2,&si2);
            build_rotation_L(ci1,si1,L1);
            build_rotation_L(ci2,si2,L2);
            mat3_mul(P,L1,tmp);
            mat3_mul(L2,tmp,Z);
            for (int m = 0; m < m_count; ++m) {
                const double c = cm[(size_t)m * (size_t)nphi + (size_t)q];
                const double ss = sm[(size_t)m * (size_t)nphi + (size_t)q];
                c_iq[m] += Z[1] * c;
                c_qq[m] += Z[4] * c;
                /* Existing OCRT primary-U convention carries a minus sign. */
                s_uq[m] += (-Z[7]) * ss;
            }
        }
        for (int m = 0; m < m_count; ++m) {
            const size_t z = (size_t)m * (size_t)directions + (size_t)ko;
            i_from_q[z] = c_iq[m] / (double)nphi;
            q_from_q[z] = c_qq[m] / (double)nphi;
            u_from_q[z] = s_uq[m] / (double)nphi;
        }
    }
    free(c_iq); free(c_qq); free(s_uq);
    return 0;
}

int rt_aerosol_value_phase_fourier_pol(rt_legendre_workspace_t *ws,
    const rt_atm_t *atm, int m, const rt_value_phase_interp_t *phase, int nphi)
{
    if (!ws || !atm || m<0 || m>ws->l_max) return -1;
    const int n_mu=ws->n_mu, kw=2*n_mu+1;
    const size_t plane=(size_t)(n_mu+1)*(size_t)kw;
    double *buf=(double*)calloc((size_t)6*(size_t)(m+1)*plane,sizeof(double));
    if(!buf)return -2;
    double *p=buf;
    double *gr=p+(size_t)(m+1)*plane;
    double *gt=gr+(size_t)(m+1)*plane;
    double *arr=gt+(size_t)(m+1)*plane;
    double *art=arr+(size_t)(m+1)*plane;
    double *att=art+(size_t)(m+1)*plane;
    int rc=rt_aerosol_value_phase_fourier_pol_allm(atm,m+1,phase,nphi,p,gr,gt,arr,art,att);
    if(rc==0){
        const size_t off=(size_t)m*plane;
        memcpy(ws->pfm_storage,p+off,plane*sizeof(double));
        memcpy(ws->gr_storage,gr+off,plane*sizeof(double));
        memcpy(ws->gt_storage,gt+off,plane*sizeof(double));
        memcpy(ws->arr_storage,arr+off,plane*sizeof(double));
        memcpy(ws->art_storage,art+off,plane*sizeof(double));
        memcpy(ws->att_storage,att+off,plane*sizeof(double));
        ws->sos_kernel_pack_valid=0; ws->sos_kernel_pack_m=-1;
    }
    free(buf); return rc;
}
