#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void surface_T_aw_coxmunk_trig_test(double mu_i_air, double mu_o_water,
                                    double cphi, double sphi,
                                    double ws, int sigma_type,
                                    double n_water, int q_convention,
                                    double *T);
int surface_T_aw_coxmunk_fourier_kernel(const double *mu_o, int n_o,
                                        const double *mu_i, int n_i,
                                        int m, int n_phi_quad,
                                        double ws, int sigma_type,
                                        double n_water, int q_convention,
                                        double *T_m);

static int close_value(const char *name, double got, double expected, double tol)
{
    const double err = fabs(got - expected);
    if (err > tol) {
        fprintf(stderr,
                "FAIL %s got=%.17g expected=%.17g |d|=%.3e tol=%.3e\n",
                name, got, expected, err, tol);
        return 1;
    }
    return 0;
}

static int test_polar_limit(void)
{
    const double n = 1.34;
    const double mu_i = cos(40.0 * M_PI / 180.0);
    const double theta_o = 1.0e-4 * M_PI / 180.0;
    const double mu_o = cos(theta_o);
    const double phi = M_PI + 30.0 * M_PI / 180.0;
    const double cphi = cos(phi), sphi = sin(phi);
    const double si = sqrt(fmax(0.0, 1.0 - mu_i * mu_i));
    const double so = sqrt(fmax(0.0, 1.0 - mu_o * mu_o));

    /* Independent microfacet incidence cosine.  These vectors are the
     * away-from-interface geometry used by the facet construction only. */
    const double ix = si, iy = 0.0, iz = mu_i;
    const double ox = so * cphi, oy = so * sphi, oz = -mu_o;
    double hx = -(ix + n * ox);
    double hy = -(iy + n * oy);
    double hz = -(iz + n * oz);
    const double hn = sqrt(hx*hx + hy*hy + hz*hz);
    hx /= hn; hy /= hn; hz /= hn;
    if (hz < 0.0) { hx=-hx; hy=-hy; hz=-hz; }
    const double cti = ix*hx + iy*hy + iz*hz;
    const double sti = sqrt(fmax(0.0, 1.0 - cti*cti));
    const double stt = sti / n;
    const double ctt = sqrt(fmax(0.0, 1.0 - stt*stt));
    const double rs = (cti - n*ctt) / (cti + n*ctt);
    const double rp = (n*cti - ctt) / (n*cti + ctt);
    const double Ts = 1.0 - rs*rs;
    const double Tp = 1.0 - rp*rp;
    const double A = 0.5 * (Ts + Tp);
    const double B = 0.5 * (Tp - Ts);
    const double C = sqrt(Ts * Tp);

    /* In the outgoing-pole limit sigma1->0 and sigma2->psi, with
     * psi = phi-pi = 30 deg for the physical propagation directions. */
    const double psi = 30.0 * M_PI / 180.0;
    const double c2 = cos(2.0 * psi), s2 = sin(2.0 * psi);
    const double expected[6] = {
        (B/A) * c2, c2, -(C/A) * s2,
        (B/A) * s2, s2,  (C/A) * c2
    };
    const int index[6] = {3,4,5,6,7,8};

    double T[9] = {0};
    surface_T_aw_coxmunk_trig_test(mu_i, mu_o, cphi, sphi,
                                    3.0, 1, n, 1, T);
    if (!(T[0] > 0.0)) {
        fprintf(stderr, "FAIL polar limit: T11 is non-positive\n");
        return 1;
    }
    int failed = 0;
    for (int k=0; k<6; ++k) {
        char name[48];
        snprintf(name, sizeof(name), "polar M%d%d/M11",
                 index[k]/3+1, index[k]%3+1);
        failed |= close_value(name, T[index[k]]/T[0], expected[k], 2.0e-6);
    }
    return failed;
}

static int test_fourier_u_column(void)
{
    const double mu_i[1] = {cos(40.0 * M_PI / 180.0)};
    const double mu_o[1] = {0.82};
    const int mode = 1;
    const int nphi = 4096;
    const double dphi = 2.0 * M_PI / (double)nphi;
    const int is_cos[9] = {1,1,0, 1,1,0, 0,0,1};
    double raw[9] = {0};
    for (int p=0; p<nphi; ++p) {
        const double phi = (p + 0.5) * dphi;
        const double cm = cos(mode * phi);
        const double sm = sin(mode * phi);
        double T[9];
        surface_T_aw_coxmunk_trig_test(mu_i[0], mu_o[0], cos(phi), sin(phi),
                                        3.0, 1, 1.34, 1, T);
        for (int k=0; k<9; ++k) raw[k] += T[k] * (is_cos[k] ? cm : sm);
    }
    for (int k=0; k<9; ++k) raw[k] *= dphi / M_PI;

    double stored[9] = {0};
    if (surface_T_aw_coxmunk_fourier_kernel(mu_o, 1, mu_i, 1,
                                             mode, nphi, 3.0, 1,
                                             1.34, 1, stored) != 0) {
        fprintf(stderr, "FAIL Fourier kernel call\n");
        return 1;
    }
    int failed = 0;
    for (int k=0; k<9; ++k) {
        const double expected = (k == 2 || k == 5) ? -raw[k] : raw[k];
        char name[48];
        snprintf(name, sizeof(name), "Fourier M%d%d", k/3+1, k%3+1);
        failed |= close_value(name, stored[k], expected,
                              2.0e-12 * fmax(1.0, fabs(expected)));
    }
    return failed;
}

int main(void)
{
    int failed = 0;
    failed |= test_polar_limit();
    failed |= test_fourier_u_column();
    if (failed) return EXIT_FAILURE;
    puts("PASS: air-to-water signed-rotation and Fourier U-column regression");
    return EXIT_SUCCESS;
}
