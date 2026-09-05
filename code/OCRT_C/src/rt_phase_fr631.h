#ifndef OCRT_RT_PHASE_FR631_H
#define OCRT_RT_PHASE_FR631_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RT_FR631_N_ANGLE 631
#define RT_FR631_LAST_INDEX 630
#define RT_FR631_GRID_ID "FR631-2026-08-15"

typedef struct {
    uint16_t lower;
    double weight;
} rt_fr631_bracket_t;

/* Fill the canonical ascending 0..180 degree FR631 grid. */
void rt_fr631_fill_theta(double theta_deg[RT_FR631_N_ANGLE]);

/* Validate an ascending or descending FR631 grid. Returns 1 on success. */
int rt_fr631_validate_theta(const double *theta_deg,
                            size_t n,
                            double abs_tol_deg,
                            int *is_descending);

/* O(1), branch-limited bracket on the canonical ascending FR631 grid. */
static inline rt_fr631_bracket_t rt_fr631_bracket(double theta_deg)
{
    rt_fr631_bracket_t result;
    double x;
    int base;
    int k;

    if (!(theta_deg > 0.0)) {
        result.lower = 0u;
        result.weight = 0.0;
        return result;
    }
    if (theta_deg >= 180.0) {
        result.lower = 629u;
        result.weight = 1.0;
        return result;
    }

    if (theta_deg < 0.2) {
        x = theta_deg * 200.0;
        base = 0;
    } else if (theta_deg < 1.0) {
        x = (theta_deg - 0.2) * 50.0;
        base = 40;
    } else if (theta_deg < 5.0) {
        x = (theta_deg - 1.0) * 20.0;
        base = 80;
    } else if (theta_deg < 20.0) {
        x = (theta_deg - 5.0) * 10.0;
        base = 160;
    } else {
        x = (theta_deg - 20.0) * 2.0;
        base = 310;
    }

    k = (int)x; /* x is nonnegative; truncation equals floor. */
    result.lower = (uint16_t)(base + k);
    result.weight = x - (double)k;
    return result;
}

/* Evaluate P11/P12/P33 with one common bracket and weight.
 * descending=0: arrays are theta=0..180.
 * descending=1: arrays are theta=180..0, as written by legacy .mie writers.
 */
void rt_fr631_eval3(const double *p11,
                    const double *p12,
                    const double *p33,
                    int descending,
                    double theta_deg,
                    double *out_p11,
                    double *out_p12,
                    double *out_p33);

/* Generic monotone-theta linear fallback for legacy non-FR631 files. */
int rt_phase_theta_linear_eval3(const double *theta_deg,
                                const double *p11,
                                const double *p12,
                                const double *p33,
                                size_t n,
                                double theta_query_deg,
                                double *out_p11,
                                double *out_p12,
                                double *out_p33);

#ifdef __cplusplus
}
#endif

#endif
