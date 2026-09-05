#include <math.h>
#include <stdio.h>
#include <string.h>

#include "rt_solver.h"

static int all_exact_zero(const double *x, int lo, int hi) {
    for (int i = lo; i <= hi; ++i) {
        if (x[i] != 0.0) return 0;
    }
    return 1;
}

int main(void) {
    enum { N_MU = 8, SRC_M_MAX = 1, ATM_M_MAX = 4 };
    double srcI[(SRC_M_MAX + 1) * N_MU];
    double srcQ[(SRC_M_MAX + 1) * N_MU];
    double srcU[(SRC_M_MAX + 1) * N_MU];
    double pmI[ATM_M_MAX + 1], pmQ[ATM_M_MAX + 1], pmU[ATM_M_MAX + 1];
    memset(srcI, 0, sizeof srcI);
    memset(srcQ, 0, sizeof srcQ);
    memset(srcU, 0, sizeof srcU);
    memset(pmI, 0, sizeof pmI);
    memset(pmQ, 0, sizeof pmQ);
    memset(pmU, 0, sizeof pmU);

    for (int k = 0; k < N_MU; ++k) {
        srcI[k] = 1.0e-4 * (double)(k + 1);
        srcI[N_MU + k] = 2.0e-5 * (double)(k + 1);
        srcQ[N_MU + k] = -5.0e-6 * (double)(k + 1);
        srcU[N_MU + k] = 3.0e-6 * (double)(k + 1);
    }

    rt_case_t cs;
    memset(&cs, 0, sizeof cs);
    cs.sza_deg = 30.0;
    cs.vza_deg = 20.0;
    cs.raa_deg = 90.0;
    cs.wavelength_nm = 443.0;
    cs.surface = RT_SURFACE_FLAT;
    cs.rayleigh_on = 1;
    cs.pressure_hpa = 1013.25;
    cs.n_water = 1.34;
    cs.F_sun = RT_F_SOLAR_PI;

    rt_options_t opts = rt_options_default();
    opts.n_mu = N_MU;
    opts.n_layers = 12;
    opts.fourier_m_max = ATM_M_MAX;
    opts.view_as_node = 0;
    opts.bottom_source_only = 1;
    opts.ext_bottom_per_m_I = srcI;
    opts.ext_bottom_per_m_Q = srcQ;
    opts.ext_bottom_per_m_U = srcU;
    opts.ext_bottom_n_mu = N_MU;
    opts.ext_bottom_m_max = SRC_M_MAX;
    opts.toa_view_per_m_I = pmI;
    opts.toa_view_per_m_Q = pmQ;
    opts.toa_view_per_m_U = pmU;

    rt_result_t out;
    int rc = rt_solve_case_pol(&cs, &opts, &out);
    if (rc != 0) {
        fprintf(stderr, "valid short source failed rc=%d\n", rc);
        return 1;
    }
    if (!(fabs(pmI[0]) > 0.0 || fabs(pmQ[0]) > 0.0 || fabs(pmU[0]) > 0.0) ||
        !(fabs(pmI[1]) > 0.0 || fabs(pmQ[1]) > 0.0 || fabs(pmU[1]) > 0.0)) {
        fprintf(stderr, "provided modes were not transported\n");
        return 2;
    }
    if (!all_exact_zero(pmI, 2, ATM_M_MAX) ||
        !all_exact_zero(pmQ, 2, ATM_M_MAX) ||
        !all_exact_zero(pmU, 2, ATM_M_MAX)) {
        fprintf(stderr, "unprovided modes are not exact zero\n");
        return 3;
    }

    opts.ext_bottom_n_mu = N_MU - 1;
    rc = rt_solve_case_pol(&cs, &opts, &out);
    if (rc != RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE) {
        fprintf(stderr, "shape mismatch rc=%d expected=%d\n",
                rc, RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE);
        return 4;
    }

    opts.ext_bottom_n_mu = N_MU;
    opts.ext_bottom_per_m_U = NULL;
    rc = rt_solve_case_pol(&cs, &opts, &out);
    if (rc != RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE) {
        fprintf(stderr, "partial source rc=%d expected=%d\n",
                rc, RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE);
        return 5;
    }

    puts("PASS ext-bottom-mode-bounds zero-padding-and-shape-guard");
    return 0;
}
