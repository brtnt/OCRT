#include <math.h>
#include "rt_quadrature.h"

/* Evaluate P_N(t) and P_{N-1}(t) via the three-term recurrence
 *   (n+1) P_{n+1} = (2n+1) t P_n - n P_{n-1},   P_0 = 1, P_1 = t.
 * On return: *pn = P_N(t), *pn_minus_1 = P_{N-1}(t).
 */
static void legendre_eval(int N, double t, double *pn, double *pn_minus_1) {
    double p0 = 1.0;
    double p1 = t;
    for (int n = 1; n < N; n++) {
        double p2 = ((2.0 * n + 1.0) * t * p1 - (double)n * p0) / (double)(n + 1);
        p0 = p1;
        p1 = p2;
    }
    *pn          = p1;
    *pn_minus_1  = p0;
}

int rt_quadrature_gauss_legendre_pos(int n_mu, double *mu, double *w) {
    if (n_mu < 1 || !mu || !w) return -1;

    const int    N        = 2 * n_mu;          /* total nodes on [-1, 1] */
    const double pi       = 3.14159265358979323846;
    const double tol      = 1.0e-15;
    const int    max_iter = 100;

    /* Tricomi initial guess gives roots in descending order:
     *   k = 1 → t near +1     (largest positive root)
     *   k = N → t near -1     (largest negative root)
     * Positive roots correspond to k = 1 .. n_mu.
     * Store ascending in mu[]: mu[0] = smallest positive, mu[n_mu-1] = largest.
     */
    for (int k = 1; k <= n_mu; k++) {
        double t = cos(pi * ((double)k - 0.25) / ((double)N + 0.5));

        double pn = 0.0, pn_m1 = 0.0, dp = 0.0;
        int converged = 0;
        for (int it = 0; it < max_iter; it++) {
            legendre_eval(N, t, &pn, &pn_m1);
            /* P'_N(t) = N (t P_N - P_{N-1}) / (t² - 1) */
            dp = (double)N * (t * pn - pn_m1) / (t * t - 1.0);
            double dt = -pn / dp;
            t += dt;
            if (fabs(dt) < tol * (fabs(t) + 1.0)) { converged = 1; break; }
        }
        if (!converged) return -2;

        /* Re-evaluate at the converged root for accurate weight. */
        legendre_eval(N, t, &pn, &pn_m1);
        dp = (double)N * (t * pn - pn_m1) / (t * t - 1.0);
        /* Standard weight: w = 2 / ((1 - t²) (P'_N(t))²) */
        double w_full = 2.0 / ((1.0 - t * t) * dp * dp);

        int idx = n_mu - k;     /* k=1 (largest) → idx=n_mu-1; k=n_mu (smallest) → idx=0 */
        mu[idx] = t;
        w[idx]  = w_full;
    }

    return 0;
}
