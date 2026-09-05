#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void surface_R_coxmunk_trig(double,double,double,double,double,int,double,int,double*);
void surface_T_coxmunk_trig(double,double,double,double,double,int,double,int,double*);
void surface_R_ww_coxmunk_trig(double,double,double,double,double,int,double,int,double*);
void surface_T_aw_coxmunk_trig_test(double,double,double,double,double,int,double,int,double*);

typedef void (*kernel_fn)(double,double,double,double,double,int,double,int,double*);

typedef struct {
    const char *name;
    kernel_fn fn;
    double mu1;
    double mu2;
    int pole_arg; /* 1 or 2 */
} pole_case_t;

static void eval(kernel_fn fn, double mu1, double mu2, double phi, double out[9]) {
    fn(mu1, mu2, cos(phi), sin(phi), 3.0, 1, 1.34, 1, out);
}

static double max_abs9(const double x[9]) {
    double m = 0.0;
    for (int i=0; i<9; ++i) if (fabs(x[i]) > m) m = fabs(x[i]);
    return m;
}

static int test_exact_pole_continuity(void) {
    const double mu_near = cos(5.0e-7); /* sin(theta) > 1e-8: original formula */
    const double phi = 30.0 * M_PI / 180.0;
    const pole_case_t cases[] = {
        {"R_air_out_pole", surface_R_coxmunk_trig,
         1.0, cos(40.0*M_PI/180.0), 1},
        {"R_air_in_pole", surface_R_coxmunk_trig,
         cos(40.0*M_PI/180.0), 1.0, 2},
        {"T_wa_out_pole", surface_T_coxmunk_trig,
         1.0, cos(5.0*M_PI/180.0), 1},
        {"T_wa_in_pole", surface_T_coxmunk_trig,
         cos(5.0*M_PI/180.0), 1.0, 2},
        {"R_ww_out_pole", surface_R_ww_coxmunk_trig,
         1.0, cos(25.0*M_PI/180.0), 1},
        {"R_ww_in_pole", surface_R_ww_coxmunk_trig,
         cos(25.0*M_PI/180.0), 1.0, 2},
        {"T_aw_in_pole", surface_T_aw_coxmunk_trig_test,
         1.0, cos(5.0*M_PI/180.0), 1},
        {"T_aw_out_pole", surface_T_aw_coxmunk_trig_test,
         cos(5.0*M_PI/180.0), 1.0, 2},
    };
    int failed = 0;
    for (size_t c=0; c<sizeof(cases)/sizeof(cases[0]); ++c) {
        double exact[9], near[9], d[9];
        double a = cases[c].mu1, b = cases[c].mu2;
        eval(cases[c].fn, a, b, phi, exact);
        if (cases[c].pole_arg == 1) a = mu_near; else b = mu_near;
        eval(cases[c].fn, a, b, phi, near);
        for (int k=0; k<9; ++k) d[k] = exact[k] - near[k];
        const double scale = fmax(max_abs9(near), 1.0e-300);
        const double rel = max_abs9(d) / scale;
        printf("pole-continuity %-16s rel=%.6e M11=%.9e\n",
               cases[c].name, rel, exact[0]);
        if (!(rel < 1.0e-4)) {
            fprintf(stderr, "FAIL pole continuity %s rel=%.17g\n",
                    cases[c].name, rel);
            failed = 1;
        }
    }
    return failed;
}

typedef struct {
    const char *name;
    kernel_fn fn;
    double mu1;
    double mu2;
} spin_case_t;

static int test_vertical_spin2_covariance(void) {
    const spin_case_t cases[] = {
        {"R_air", surface_R_coxmunk_trig,
         1.0, cos(40.0*M_PI/180.0)},
        {"T_wa", surface_T_coxmunk_trig,
         1.0, cos(5.0*M_PI/180.0)},
        {"R_ww", surface_R_ww_coxmunk_trig,
         1.0, cos(25.0*M_PI/180.0)},
        /* For T_aw the outgoing water direction is the second argument. */
        {"T_aw", surface_T_aw_coxmunk_trig_test,
         cos(5.0*M_PI/180.0), 1.0},
    };
    enum { NPHI = 4096 };
    const double dphi = 2.0*M_PI/(double)NPHI;
    int failed = 0;
    for (size_t c=0; c<sizeof(cases)/sizeof(cases[0]); ++c) {
        double q0=0.0,u0=0.0,qc2=0.0,qs2=0.0,uc2=0.0,us2=0.0;
        for (int p=0; p<NPHI; ++p) {
            const double phi = (p+0.5)*dphi;
            double M[9];
            eval(cases[c].fn, cases[c].mu1, cases[c].mu2, phi, M);
            const double q = M[3]; /* unpolarized input -> Q output */
            const double u = M[6]; /* unpolarized input -> U output */
            q0 += q; u0 += u;
            qc2 += q*cos(2.0*phi); qs2 += q*sin(2.0*phi);
            uc2 += u*cos(2.0*phi); us2 += u*sin(2.0*phi);
        }
        q0 /= NPHI; u0 /= NPHI;
        qc2 *= 2.0/NPHI; qs2 *= 2.0/NPHI;
        uc2 *= 2.0/NPHI; us2 *= 2.0/NPHI;
        const double amp = fmax(fabs(qc2), fabs(us2));
        const double leakage = fmax(fmax(fabs(q0),fabs(u0)),
                                    fmax(fabs(qs2),fabs(uc2)));
        const double pair_err = fabs(fabs(qc2)-fabs(us2));
        printf("spin2 %-5s amp=%.9e leakage/amp=%.3e pair/amp=%.3e\n",
               cases[c].name, amp, leakage/fmax(amp,1e-300),
               pair_err/fmax(amp,1e-300));
        if (!(amp > 1.0e-12 && leakage/amp < 1.0e-10 &&
              pair_err/amp < 1.0e-10)) {
            fprintf(stderr, "FAIL spin-2 pole covariance %s\n", cases[c].name);
            failed = 1;
        }
    }
    return failed;
}

int main(void) {
    int failed = 0;
    failed |= test_exact_pole_continuity();
    failed |= test_vertical_spin2_covariance();
    if (failed) return EXIT_FAILURE;
    puts("PASS: exact-pole meridian rotation limits and vertical spin-2 covariance");
    return EXIT_SUCCESS;
}
