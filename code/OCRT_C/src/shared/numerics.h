/* shared/numerics.h
 *
 * Numerical utilities used by both vrt_solver and ocean_solver:
 *   - Gauss-Legendre quadrature
 *   - linear and PCHIP interpolation
 *   - Legendre polynomial recurrence
 *   - trapezoidal integration
 *
 * Implementation in shared/numerics.c.
 */
#ifndef OCRT_SHARED_NUMERICS_H
#define OCRT_SHARED_NUMERICS_H

#include <stdlib.h>
#include <string.h>

/* Gauss-Legendre nodes/weights on [-1, +1]. Fills x[n], w[n]. */
void gauss_legendre(int n, double* x, double* w);

/* Trapezoidal integration of y[] on monotone grid x[], length n. */
double trapezoid_int(const double* x, const double* y, int n);

/* Linear interpolation y(xq) given monotonically INCREASING x[n], y[n].
 * Returns y[0] if xq < x[0], y[n-1] if xq > x[n-1]. */
double linterp(const double* x, const double* y, int n, double xq);

/* One-shot shape-preserving cubic Hermite interpolation on an ascending grid.
 * This uses the same Fritsch-Carlson/PCHIP slopes as pchip_build(), but avoids
 * heap allocation when a table is evaluated only once (for example, one
 * phase-matrix angle across the .mie wavelength grid). Queries outside the
 * tabulated wavelength range are clamped to the nearest endpoint, matching
 * linterp() and OCRT's no-spectral-extrapolation policy. */
double pchip_interp(const double* x, const double* y, int n, double xq);

/* Legendre polynomial P_l(x) via recurrence. */
double legendre_P(int l, double x);

/* Fill P_l(x) for l = 0..n_max-1 into out[n_max]. */
void legendre_P_all(int n_max, double x, double* out);

/* PCHIP (Fritsch-Carlson 1980) shape-preserving cubic interpolation.
 * Build once from (x[n], y[n]) ascending, evaluate many times.
 * pchip_t owns copies of x, y, m -- must free via pchip_free. */
typedef struct {
    int n;
    double* x;
    double* y;
    double* m;
} pchip_t;

int pchip_build(pchip_t* p, const double* x, const double* y, int n);
double pchip_eval(const pchip_t* p, double xq);
void pchip_free(pchip_t* p);

#endif /* OCRT_SHARED_NUMERICS_H */
