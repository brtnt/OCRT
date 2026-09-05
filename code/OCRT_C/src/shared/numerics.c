/* shared/numerics.c
 *
 * Implementation of numerical utilities. Sourced from vrt_solver.c v4.0
 * (validated against PSTAR/SeaDAS REF1/REF2 Rayleigh references and
 * Python mie_bh.py bit-level).
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "numerics.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Gauss-Legendre quadrature on [-1, +1] (Newton iteration on Legendre roots). */
void gauss_legendre(int n, double* x, double* w) {
    const int m = (n + 1) / 2;   /* nodes are symmetric */
    for (int i = 0; i < m; ++i) {
        /* Initial guess from A&S 25.4.29 */
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));
        double z1, pp = 0.0;
        int iter;
        for (iter = 0; iter < 50; ++iter) {
            /* Evaluate P_n(z) by upward recurrence */
            double p1 = 1.0, p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / (double)j;
            }
            /* Derivative: P'_n(z) = n (z P_n - P_{n-1}) / (z^2 - 1) */
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
            if (fabs(z - z1) < 1.0e-14) break;
        }
        /* x ordered ascending: -|z| first, +|z| last */
        x[i]           = -z;
        x[n - 1 - i]   =  z;
        w[i]           = 2.0 / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i]   = w[i];
    }
}

double trapezoid_int(const double* x, const double* y, int n) {
    double s = 0.0;
    for (int i = 1; i < n; ++i) {
        s += 0.5 * (x[i] - x[i-1]) * (y[i] + y[i-1]);
    }
    return s;
}

double linterp(const double* x, const double* y, int n, double xq) {
    if (n == 0) return 0.0;
    if (n == 1) return y[0];
    if (xq <= x[0])   return y[0];
    if (xq >= x[n-1]) return y[n-1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (x[mid] <= xq) lo = mid; else hi = mid;
    }
    double t = (xq - x[lo]) / (x[hi] - x[lo]);
    return y[lo] + t * (y[hi] - y[lo]);
}

/* PCHIP node derivative. The formula follows the weighted harmonic mean
 * used by scipy.interpolate.PchipInterpolator and OSOAA conv_local. */
static double pchip_node_slope(const double* x, const double* y, int n, int k) {
    if (n == 2) return (y[1] - y[0]) / (x[1] - x[0]);

    if (k == 0) {
        const double h0 = x[1] - x[0];
        const double h1 = x[2] - x[1];
        const double d0 = (y[1] - y[0]) / h0;
        const double d1 = (y[2] - y[1]) / h1;
        double m = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (m * d0 <= 0.0) return 0.0;
        if (d0 * d1 <= 0.0 && fabs(m) > 3.0 * fabs(d0)) return 3.0 * d0;
        return m;
    }
    if (k == n - 1) {
        const double h0 = x[n - 1] - x[n - 2];
        const double h1 = x[n - 2] - x[n - 3];
        const double d0 = (y[n - 1] - y[n - 2]) / h0;
        const double d1 = (y[n - 2] - y[n - 3]) / h1;
        double m = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (m * d0 <= 0.0) return 0.0;
        if (d0 * d1 <= 0.0 && fabs(m) > 3.0 * fabs(d0)) return 3.0 * d0;
        return m;
    }

    const double hm = x[k] - x[k - 1];
    const double hp = x[k + 1] - x[k];
    const double dm = (y[k] - y[k - 1]) / hm;
    const double dp = (y[k + 1] - y[k]) / hp;
    if (dm == 0.0 || dp == 0.0 || dm * dp <= 0.0) return 0.0;

    const double w1 = 2.0 * hp + hm;
    const double w2 = hp + 2.0 * hm;
    return (w1 + w2) / (w1 / dm + w2 / dp);
}

double pchip_interp(const double* x, const double* y, int n, double xq) {
    if (!x || !y || n <= 0) return 0.0;
    if (n == 1) return y[0];
    if (xq <= x[0]) return y[0];
    if (xq >= x[n - 1]) return y[n - 1];

    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (x[mid] <= xq) lo = mid;
        else hi = mid;
    }

    const double h = x[hi] - x[lo];
    if (!(h > 0.0)) return linterp(x, y, n, xq);
    const double m0 = pchip_node_slope(x, y, n, lo);
    const double m1 = pchip_node_slope(x, y, n, hi);
    const double t = (xq - x[lo]) / h;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;
    return h00 * y[lo] + h10 * h * m0 + h01 * y[hi] + h11 * h * m1;
}

double legendre_P(int l, double x) {
    if (l == 0) return 1.0;
    if (l == 1) return x;
    double p0 = 1.0, p1 = x, p2 = 0.0;
    for (int k = 1; k < l; ++k) {
        p2 = ((2.0 * k + 1.0) * x * p1 - (double)k * p0) / (double)(k + 1);
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

void legendre_P_all(int n_max, double x, double* out) {
    if (n_max <= 0) return;
    out[0] = 1.0;
    if (n_max == 1) return;
    out[1] = x;
    for (int k = 1; k < n_max - 1; ++k) {
        out[k + 1] = ((2.0 * k + 1.0) * x * out[k] - (double)k * out[k - 1]) / (double)(k + 1);
    }
}

int pchip_build(pchip_t* p, const double* x, const double* y, int n) {
    if (n < 2) return -1;
    p->n = n;
    p->x = (double*)malloc(n * sizeof(double));
    p->y = (double*)malloc(n * sizeof(double));
    p->m = (double*)malloc(n * sizeof(double));
    if (!p->x || !p->y || !p->m) {
        free(p->x); free(p->y); free(p->m);
        return -1;
    }
    memcpy(p->x, x, n * sizeof(double));
    memcpy(p->y, y, n * sizeof(double));

    /* Secant slopes d_k = (y[k+1] - y[k]) / (x[k+1] - x[k]) for k = 0..n-2 */
    double* d = (double*)malloc((n - 1) * sizeof(double));
    double* h = (double*)malloc((n - 1) * sizeof(double));
    if (!d || !h) { free(d); free(h); free(p->x); free(p->y); free(p->m); return -1; }
    for (int k = 0; k < n - 1; ++k) {
        h[k] = x[k + 1] - x[k];
        d[k] = (y[k + 1] - y[k]) / h[k];
    }

    /* Interior slopes (Fritsch-Carlson 1980) */
    for (int k = 1; k < n - 1; ++k) {
        if (d[k - 1] * d[k] <= 0.0) {
            p->m[k] = 0.0;
        } else {
            double w1 = 2.0 * h[k] + h[k - 1];
            double w2 = h[k] + 2.0 * h[k - 1];
            p->m[k] = (w1 + w2) / (w1 / d[k - 1] + w2 / d[k]);
        }
    }

    /* Endpoint slopes (3-point edge formula with shape-preserving clamp) */
    p->m[0] = ((2.0 * h[0] + h[1]) * d[0] - h[0] * d[1]) / (h[0] + h[1]);
    if (p->m[0] * d[0] <= 0.0) p->m[0] = 0.0;
    else if ((d[0] * d[1] <= 0.0) && (fabs(p->m[0]) > 3.0 * fabs(d[0]))) p->m[0] = 3.0 * d[0];

    int nn = n - 1;
    p->m[nn] = ((2.0 * h[nn - 1] + h[nn - 2]) * d[nn - 1] - h[nn - 1] * d[nn - 2])
               / (h[nn - 1] + h[nn - 2]);
    if (p->m[nn] * d[nn - 1] <= 0.0) p->m[nn] = 0.0;
    else if ((d[nn - 1] * d[nn - 2] <= 0.0) && (fabs(p->m[nn]) > 3.0 * fabs(d[nn - 1])))
        p->m[nn] = 3.0 * d[nn - 1];

    free(d); free(h);
    return 0;
}

double pchip_eval(const pchip_t* p, double xq) {
    const int n = p->n;
    const double* x = p->x;
    const double* y = p->y;
    const double* m = p->m;

    /* Bracket interval via binary search */
    int lo = 0, hi = n - 1;
    if (xq <= x[0]) { lo = 0; hi = 1; }
    else if (xq >= x[n - 1]) { lo = n - 2; hi = n - 1; }
    else {
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x[mid] <= xq) lo = mid; else hi = mid;
        }
    }
    const double h = x[hi] - x[lo];
    const double t = (xq - x[lo]) / h;
    const double t2 = t * t, t3 = t2 * t;
    const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    const double h10 = t3 - 2.0 * t2 + t;
    const double h01 = -2.0 * t3 + 3.0 * t2;
    const double h11 = t3 - t2;
    return h00 * y[lo] + h10 * h * m[lo] + h01 * y[hi] + h11 * h * m[hi];
}

void pchip_free(pchip_t* p) {
    free(p->x); free(p->y); free(p->m);
    memset(p, 0, sizeof(*p));
}
