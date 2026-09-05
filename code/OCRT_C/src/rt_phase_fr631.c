#include "rt_phase_fr631.h"

#include <math.h>
#include <stddef.h>

static double lerp(double a, double b, double w)
{
#if defined(FP_FAST_FMA) || defined(__FMA__)
    return fma(w, b - a, a);
#else
    return a + w * (b - a);
#endif
}

void rt_fr631_fill_theta(double theta_deg[RT_FR631_N_ANGLE])
{
    int i;
    for (i = 0; i <= 40; ++i) {
        theta_deg[i] = 0.005 * (double)i;
    }
    for (i = 41; i <= 80; ++i) {
        theta_deg[i] = 0.2 + 0.02 * (double)(i - 40);
    }
    for (i = 81; i <= 160; ++i) {
        theta_deg[i] = 1.0 + 0.05 * (double)(i - 80);
    }
    for (i = 161; i <= 310; ++i) {
        theta_deg[i] = 5.0 + 0.1 * (double)(i - 160);
    }
    for (i = 311; i <= 630; ++i) {
        theta_deg[i] = 20.0 + 0.5 * (double)(i - 310);
    }
    theta_deg[0] = 0.0;
    theta_deg[40] = 0.2;
    theta_deg[80] = 1.0;
    theta_deg[160] = 5.0;
    theta_deg[310] = 20.0;
    theta_deg[630] = 180.0;
}

int rt_fr631_validate_theta(const double *theta_deg,
                            size_t n,
                            double abs_tol_deg,
                            int *is_descending)
{
    double expected[RT_FR631_N_ANGLE];
    int descending;
    size_t i;

    if (theta_deg == NULL || n != RT_FR631_N_ANGLE || abs_tol_deg < 0.0) {
        return 0;
    }
    rt_fr631_fill_theta(expected);
    descending = theta_deg[0] > theta_deg[n - 1];

    for (i = 0; i < n; ++i) {
        const size_t j = descending ? (n - 1u - i) : i;
        if (!isfinite(theta_deg[j]) || fabs(theta_deg[j] - expected[i]) > abs_tol_deg) {
            return 0;
        }
    }
    if (is_descending != NULL) {
        *is_descending = descending;
    }
    return 1;
}

void rt_fr631_eval3(const double *p11,
                    const double *p12,
                    const double *p33,
                    int descending,
                    double theta_deg,
                    double *out_p11,
                    double *out_p12,
                    double *out_p33)
{
    const rt_fr631_bracket_t bracket = rt_fr631_bracket(theta_deg);
    size_t i0 = (size_t)bracket.lower;
    size_t i1 = i0 + 1u;

    if (descending) {
        i0 = (size_t)RT_FR631_LAST_INDEX - i0;
        i1 = i0 - 1u;
    }

    *out_p11 = lerp(p11[i0], p11[i1], bracket.weight);
    *out_p12 = lerp(p12[i0], p12[i1], bracket.weight);
    *out_p33 = lerp(p33[i0], p33[i1], bracket.weight);
}

static int theta_is_ascending(const double *theta_deg, size_t n)
{
    return n >= 2u && theta_deg[n - 1u] > theta_deg[0];
}

int rt_phase_theta_linear_eval3(const double *theta_deg,
                                const double *p11,
                                const double *p12,
                                const double *p33,
                                size_t n,
                                double theta_query_deg,
                                double *out_p11,
                                double *out_p12,
                                double *out_p33)
{
    size_t lo = 0u;
    size_t hi;
    int ascending;
    double x0;
    double x1;
    double w;

    if (theta_deg == NULL || p11 == NULL || p12 == NULL || p33 == NULL ||
        out_p11 == NULL || out_p12 == NULL || out_p33 == NULL || n < 2u) {
        return 0;
    }

    ascending = theta_is_ascending(theta_deg, n);
    if (ascending) {
        if (theta_query_deg <= theta_deg[0]) {
            lo = 0u;
        } else if (theta_query_deg >= theta_deg[n - 1u]) {
            lo = n - 2u;
            theta_query_deg = theta_deg[n - 1u];
        } else {
            hi = n - 1u;
            while (hi - lo > 1u) {
                const size_t mid = lo + (hi - lo) / 2u;
                if (theta_deg[mid] <= theta_query_deg) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        }
        x0 = theta_deg[lo];
        x1 = theta_deg[lo + 1u];
    } else {
        if (theta_query_deg >= theta_deg[0]) {
            lo = 0u;
            theta_query_deg = theta_deg[0];
        } else if (theta_query_deg <= theta_deg[n - 1u]) {
            lo = n - 2u;
            theta_query_deg = theta_deg[n - 1u];
        } else {
            hi = n - 1u;
            while (hi - lo > 1u) {
                const size_t mid = lo + (hi - lo) / 2u;
                if (theta_deg[mid] >= theta_query_deg) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        }
        x0 = theta_deg[lo];
        x1 = theta_deg[lo + 1u];
    }

    w = (theta_query_deg - x0) / (x1 - x0);
    *out_p11 = lerp(p11[lo], p11[lo + 1u], w);
    *out_p12 = lerp(p12[lo], p12[lo + 1u], w);
    *out_p33 = lerp(p33[lo], p33[lo + 1u], w);
    return 1;
}
