#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rt_molecular_profile.h"

static int nearly(double a, double b, double atol, double rtol)
{
    const double d = fabs(a - b);
    return d <= atol + rtol * fmax(fabs(a), fabs(b));
}

int main(void)
{
    static const double z[] = {
        0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852
    };
    static const double p[] = {
        101325.0, 22632.1, 5474.89, 868.019,
        110.906, 66.9389, 3.95642, 0.3734
    };
    int fail = 0;

    for (size_t i = 0; i < sizeof(z) / sizeof(z[0]); ++i) {
        const double got = rt_molecular_us62_pressure_pa(z[i]);
        if (!nearly(got, p[i], 1e-12, 2e-15)) {
            fprintf(stderr, "anchor mismatch i=%zu z=%.9g got=%.17g expected=%.17g\n",
                    i, z[i], got, p[i]);
            fail++;
        }
    }

    double prev = 1.0;
    for (int i = 0; i <= 1000; ++i) {
        const double f = (double)i / 1000.0;
        const double alt = rt_molecular_us62_altitude_from_grid_fraction(f);
        const double back = rt_molecular_us62_grid_fraction_above(alt);
        if (!nearly(f, back, 2e-14, 2e-14)) {
            fprintf(stderr, "inverse mismatch f=%.17g z=%.17g back=%.17g\n",
                    f, alt, back);
            fail++;
            break;
        }
        const double current = rt_molecular_us62_grid_fraction_above(
            rt_molecular_us62_top_km() * (double)i / 1000.0);
        if (current > prev + 2e-15) {
            fprintf(stderr, "non-monotone cumulative at i=%d prev=%.17g current=%.17g\n",
                    i, prev, current);
            fail++;
            break;
        }
        prev = current;
    }

    double sum = 0.0;
    const int n = 400;
    const double ztop = rt_molecular_us62_top_km();
    for (int i = 0; i < n; ++i) {
        const double z0 = ztop * (double)i / (double)n;
        const double z1 = ztop * (double)(i + 1) / (double)n;
        sum += rt_molecular_us62_layer_fraction(z0, z1);
    }
    if (!nearly(sum, 1.0, 3e-14, 3e-14)) {
        fprintf(stderr, "column closure mismatch sum=%.17g\n", sum);
        fail++;
    }

    printf("US62_PROFILE anchors=8 inverse=1001 layer_sum=%.17g fail=%d\n",
           sum, fail);
    return fail ? EXIT_FAILURE : EXIT_SUCCESS;
}
