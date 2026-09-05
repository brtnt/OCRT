#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double rt_test_aw_interp_on_unsorted(const double *mu, const double *vals,
                                      int n, double mu_target);

static int run_case(int n, double target) {
    double *mu = (double *)malloc((size_t)n * sizeof(double));
    double *v = (double *)malloc((size_t)n * sizeof(double));
    if (!mu || !v) return 1;
    /* Deliberately non-monotone permutation; v=3+2*mu must be reproduced
       exactly by linear interpolation and extrapolation. */
    for (int i = 0; i < n; ++i) {
        int j = (37 * i + 11) % n;
        mu[i] = 0.01 + 0.98 * (double)j / (double)(n - 1);
        v[i] = 3.0 + 2.0 * mu[i];
    }
    double got = rt_test_aw_interp_on_unsorted(mu, v, n, target);
    double expect = 3.0 + 2.0 * target;
    free(mu); free(v);
    if (!isfinite(got) || fabs(got - expect) > 2e-13) {
        fprintf(stderr, "FAIL n=%d target=%.17g got=%.17g expected=%.17g\n",
                n, target, got, expect);
        return 1;
    }
    return 0;
}

int main(void) {
    const int ns[] = {48, 64, 96, 320, 321, 513};
    const double targets[] = {0.005, 0.123456789, 0.995};
    for (size_t i = 0; i < sizeof(ns)/sizeof(ns[0]); ++i)
        for (size_t j = 0; j < sizeof(targets)/sizeof(targets[0]); ++j)
            if (run_case(ns[i], targets[j])) return 1;
    puts("PASS: air-water interpolation workspace stack/heap paths");
    return 0;
}
