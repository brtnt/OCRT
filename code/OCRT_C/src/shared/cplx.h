/* shared/cplx.h
 *
 * Local complex number type and operations.
 *
 * We avoid <complex.h> for portability (GCC vs MSVC incompatibility).
 * All operations are inline header-only and have no external dependencies
 * beyond <math.h>. Both vrt_solver and ocean_solver include this file.
 */
#ifndef OCRT_SHARED_CPLX_H
#define OCRT_SHARED_CPLX_H

#include <math.h>

typedef struct { double re, im; } cplx;

static inline cplx cplx_make(double r, double i) {
    cplx c; c.re = r; c.im = i; return c;
}

static inline cplx cplx_add(cplx a, cplx b) {
    cplx r; r.re = a.re + b.re; r.im = a.im + b.im; return r;
}

static inline cplx cplx_sub(cplx a, cplx b) {
    cplx r; r.re = a.re - b.re; r.im = a.im - b.im; return r;
}

static inline cplx cplx_mul(cplx a, cplx b) {
    cplx r;
    r.re = a.re * b.re - a.im * b.im;
    r.im = a.re * b.im + a.im * b.re;
    return r;
}

static inline cplx cplx_div(cplx a, cplx b) {
    double d = b.re * b.re + b.im * b.im;
    cplx r;
    r.re = (a.re * b.re + a.im * b.im) / d;
    r.im = (a.im * b.re - a.re * b.im) / d;
    return r;
}

static inline double cplx_abs(cplx a) {
    return sqrt(a.re * a.re + a.im * a.im);
}

static inline double cplx_abs2(cplx a) {
    return a.re * a.re + a.im * a.im;
}

#endif /* OCRT_SHARED_CPLX_H */
