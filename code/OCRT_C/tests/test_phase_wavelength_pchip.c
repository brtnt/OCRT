#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared/mie_io.h"
#include "shared/numerics.h"

static int failures = 0;

static void expect_close(const char *name, double got, double expected,
                         double atol, double rtol) {
    const double tol = atol + rtol * fabs(expected);
    const double err = fabs(got - expected);
    if (!(err <= tol)) {
        fprintf(stderr,
                "FAIL %-36s got=%.17g expected=%.17g abs_err=%.3e tol=%.3e\n",
                name, got, expected, err, tol);
        ++failures;
    }
}

static void test_scipy_reference_values(void) {
    /* References generated with scipy.interpolate.PchipInterpolator 1.17.0.
     * Values outside the grid use OCRT's pre-existing endpoint clamp policy. */
    static const double x1[] = {0.40, 0.55, 0.67, 0.87, 1.02};
    static const double y1[] = {1.00, 1.80, 1.30, 2.40, 2.20};
    static const double q1[] = {0.40, 0.475, 0.61, 0.77, 0.945, 1.02};
    static const double r1[] = {
        1.0,
        1.5989583333333333,
        1.5500000000000003,
        1.8500000000000001,
        2.375,
        2.2000000000000002
    };

    static const double x2[] = {0.44, 0.67, 0.87, 1.02};
    static const double y2[] = {0.10, 0.50, 0.90, 1.40};
    static const double q2[] = {0.44, 0.555, 0.77, 0.945, 1.02};
    static const double r2[] = {
        0.10000000000000001,
        0.29241307415152812,
        0.68333420355310726,
        1.1242254733218586,
        1.3999999999999999
    };

    for (int i = 0; i < (int)(sizeof q1 / sizeof q1[0]); ++i) {
        char name[64];
        snprintf(name, sizeof name, "SciPy nonmonotone q=%.3f", q1[i]);
        expect_close(name, pchip_interp(x1, y1, 5, q1[i]), r1[i], 2e-15, 2e-14);
    }
    for (int i = 0; i < (int)(sizeof q2 / sizeof q2[0]); ++i) {
        char name[64];
        snprintf(name, sizeof name, "SciPy monotone q=%.3f", q2[i]);
        expect_close(name, pchip_interp(x2, y2, 4, q2[i]), r2[i], 2e-15, 2e-14);
    }
}

static void test_exact_nodes_and_clamp(void) {
    static const double x[] = {0.41, 0.49, 0.55, 0.67, 0.865};
    static const double y[] = {4.0, -1.0, 3.5, 2.0, 8.0};
    for (int i = 0; i < 5; ++i) {
        char name[64];
        snprintf(name, sizeof name, "exact phase wavelength node %d", i);
        expect_close(name, pchip_interp(x, y, 5, x[i]), y[i], 0.0, 0.0);
    }
    expect_close("lower endpoint clamp", pchip_interp(x, y, 5, 0.30), y[0], 0.0, 0.0);
    expect_close("upper endpoint clamp", pchip_interp(x, y, 5, 1.20), y[4], 0.0, 0.0);

    {
        static const double x2[] = {0.4, 0.5};
        static const double y2[] = {2.0, 4.0};
        expect_close("two-node cubic degenerates linear",
                     pchip_interp(x2, y2, 2, 0.45), 3.0, 0.0, 1e-15);
    }
    expect_close("single-node value", pchip_interp(x, y, 1, 9.0), y[0], 0.0, 0.0);
    expect_close("empty input value", pchip_interp(x, y, 0, 0.5), 0.0, 0.0, 0.0);
}

static void fill_phase_row(double *block, int n_wl, int ia,
                           const double *row, double scale) {
    for (int iw = 0; iw < n_wl; ++iw) {
        BLOCK_AT(block, n_wl, ia, iw) = scale * row[iw];
    }
}

static void test_mie_phase_pipeline(void) {
    enum { NANG = 3, NWL = 5 };
    static const double wl[] = {0.40, 0.55, 0.67, 0.87, 1.02};
    static const double p11_row[] = {1.00, 1.80, 1.30, 2.40, 2.20};
    static const double p12_row[] = {-0.10, -0.50, -0.90, -1.40, -1.60};
    static const double p33_row[] = {0.20, 0.40, 0.30, 0.80, 0.70};
    const double query_wl = 0.61;

    mie_data_t mie;
    aerosol_phase_interp_t phase;
    memset(&mie, 0, sizeof mie);
    memset(&phase, 0, sizeof phase);
    mie_data_init(&mie, NANG, NWL, NWL);
    if (!mie.wavelengths || !mie.phase_wavelengths || !mie.angles ||
        !mie.spectral || !mie.P11 || !mie.P12 || !mie.P33) {
        fprintf(stderr, "FAIL mie_data_init allocation\n");
        ++failures;
        mie_data_free(&mie);
        return;
    }

    memcpy(mie.wavelengths, wl, sizeof wl);
    memcpy(mie.phase_wavelengths, wl, sizeof wl);
    /* Deliberately descending to exercise mie_io angle sorting. */
    mie.angles[0] = 180.0;
    mie.angles[1] = 90.0;
    mie.angles[2] = 0.0;
    for (int iw = 0; iw < NWL; ++iw) {
        SPECTRAL_AT(&mie, iw, 0) = 1.0 + 0.1 * iw;
        SPECTRAL_AT(&mie, iw, 1) = 0.9 + 0.1 * iw;
        SPECTRAL_AT(&mie, iw, 2) = 0.90 + 0.01 * iw;
        SPECTRAL_AT(&mie, iw, 3) = 0.60 + 0.02 * iw;
        SPECTRAL_AT(&mie, iw, 4) = 1.0 + 0.1 * iw;
        SPECTRAL_AT(&mie, iw, 5) = 0.9 + 0.1 * iw;
    }
    for (int ia = 0; ia < NANG; ++ia) {
        const double scale = (double)(ia + 1);
        fill_phase_row(mie.P11, NWL, ia, p11_row, scale);
        fill_phase_row(mie.P12, NWL, ia, p12_row, scale);
        fill_phase_row(mie.P33, NWL, ia, p33_row, scale);
    }

    if (build_aerosol_interpolators(&mie, query_wl, &phase) != 0) {
        fprintf(stderr, "FAIL build_aerosol_interpolators\n");
        ++failures;
        mie_data_free(&mie);
        return;
    }

    for (int ia = 0; ia < NANG; ++ia) {
        const double theta = mie.angles[ia];
        const double scale = (double)(ia + 1);
        double p11 = 0.0, p12 = 0.0, p33 = 0.0;
        char name[80];
        eval_aerosol_phase(&phase, theta, &p11, &p12, &p33);

        snprintf(name, sizeof name, "mie pipeline P11 theta=%.0f", theta);
        expect_close(name, p11, scale * 1.55, 4e-14, 4e-14);
        snprintf(name, sizeof name, "mie pipeline P12 theta=%.0f", theta);
        expect_close(name, p12,
                     scale * pchip_interp(wl, p12_row, NWL, query_wl),
                     4e-14, 4e-14);
        snprintf(name, sizeof name, "mie pipeline P33 theta=%.0f", theta);
        expect_close(name, p33,
                     scale * pchip_interp(wl, p33_row, NWL, query_wl),
                     4e-14, 4e-14);
    }

    /* Stage-3B direct particle path: a single wavelength bracket/weight is
     * shared by P11/P12/P33.  At q=0.61 the bracket is 0.55..0.67 and w=0.5. */
    {
        double l11[NANG], l12[NANG], l33[NANG];
        if (mie_phase_nodes_at_wavelength_linear(
                &mie, query_wl, l11, l12, l33) != 0) {
            fprintf(stderr, "FAIL mie_phase_nodes_at_wavelength_linear\n");
            ++failures;
        } else {
            const double e11 = 0.5 * (p11_row[1] + p11_row[2]);
            const double e12 = 0.5 * (p12_row[1] + p12_row[2]);
            const double e33 = 0.5 * (p33_row[1] + p33_row[2]);
            for (int ia = 0; ia < NANG; ++ia) {
                const double scale = (double)(ia + 1);
                expect_close("mie common-linear P11", l11[ia], scale * e11, 0.0, 2e-15);
                expect_close("mie common-linear P12", l12[ia], scale * e12, 0.0, 2e-15);
                expect_close("mie common-linear P33", l33[ia], scale * e33, 0.0, 2e-15);
                if (l11[ia] < 0.0 || fabs(l12[ia]) > l11[ia] ||
                    fabs(l33[ia]) > l11[ia]) {
                    fprintf(stderr, "FAIL common-linear physical bound ia=%d\n", ia);
                    ++failures;
                }
            }
        }
    }

    /* At an exact table wavelength the new spectral interpolation must be exact. */
    aerosol_phase_interp_free(&phase);
    memset(&phase, 0, sizeof phase);
    if (build_aerosol_interpolators(&mie, wl[2], &phase) != 0) {
        fprintf(stderr, "FAIL exact-node build_aerosol_interpolators\n");
        ++failures;
    } else {
        for (int ia = 0; ia < NANG; ++ia) {
            double p11 = 0.0, p12 = 0.0, p33 = 0.0;
            eval_aerosol_phase(&phase, mie.angles[ia], &p11, &p12, &p33);
            expect_close("mie exact-node P11", p11,
                         BLOCK_AT(mie.P11, NWL, ia, 2), 4e-14, 4e-14);
            expect_close("mie exact-node P12", p12,
                         BLOCK_AT(mie.P12, NWL, ia, 2), 4e-14, 4e-14);
            expect_close("mie exact-node P33", p33,
                         BLOCK_AT(mie.P33, NWL, ia, 2), 4e-14, 4e-14);
        }
    }

    aerosol_phase_interp_free(&phase);
    mie_data_free(&mie);
}

int main(void) {
    test_scipy_reference_values();
    test_exact_nodes_and_clamp();
    test_mie_phase_pipeline();

    if (failures != 0) {
        fprintf(stderr, "FAIL phase-wavelength PCHIP tests=%d\n", failures);
        return 1;
    }
    puts("PASS phase-wavelength PCHIP: SciPy references, endpoint clamp, exact nodes, Mie pipeline");
    return 0;
}
