#include "rt_atm.h"
#include "rt_value_phase.h"
#include "shared/mie_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static size_t kidx(int m, int j, int k, int n_mu)
{
    const int kw = 2 * n_mu + 1;
    return ((size_t)m * (size_t)(n_mu + 1) + (size_t)j) *
           (size_t)kw + (size_t)(k + n_mu);
}

int main(int argc, char **argv)
{
    if (argc != 7) {
        fprintf(stderr, "usage: %s MIE WL_NM N_MU M_COUNT NPHI OUT_CSV\n", argv[0]);
        return 2;
    }
    const char *mie_path = argv[1];
    const double wl_nm = atof(argv[2]);
    const int n_mu = atoi(argv[3]);
    const int m_count = atoi(argv[4]);
    const int nphi = atoi(argv[5]);
    const char *out_path = argv[6];
    if (!(wl_nm > 0.0) || n_mu < 1 || m_count < 1 || nphi < 16) return 3;

    mie_data_t mie;
    memset(&mie, 0, sizeof mie);
    if (read_mie_file(mie_path, &mie) != 0) return 4;
    double *p11 = (double *)malloc((size_t)mie.n_ang * sizeof(double));
    double *p12 = (double *)malloc((size_t)mie.n_ang * sizeof(double));
    double *p33 = (double *)malloc((size_t)mie.n_ang * sizeof(double));
    if (!p11 || !p12 || !p33) return 5;
    if (mie_phase_nodes_at_wavelength_linear(
            &mie, wl_nm / 1000.0, p11, p12, p33) != 0) return 6;

    rt_value_phase_interp_t phase;
    memset(&phase, 0, sizeof phase);
    if (rt_value_phase_interp_build(
            &phase, mie.angles, p11, p12, p33, mie.n_ang, 0) != 0) return 7;
    free(p11); free(p12); free(p33);
    mie_data_free(&mie);

    rt_atm_t atm;
    memset(&atm, 0, sizeof atm);
    if (rt_atm_alloc(&atm, 1, n_mu) != 0) return 8;
    if (rt_atm_build_rayleigh(&atm, 0.1, 0.0279,
                              cos(30.0 * M_PI / 180.0),
                              RT_RAYLEIGH_MODEL_BODHAINE_1999) != 0) return 9;

    const int kw = 2 * n_mu + 1;
    const size_t plane = (size_t)(n_mu + 1) * (size_t)kw;
    const size_t all = (size_t)m_count * plane;
    double *buf = (double *)calloc(6u * all, sizeof(double));
    if (!buf) return 10;
    if (rt_aerosol_value_phase_fourier_pol_allm(
            &atm, m_count, &phase, nphi,
            buf, buf + all, buf + 2u * all, buf + 3u * all,
            buf + 4u * all, buf + 5u * all) != 0) return 11;

    FILE *fp = fopen(out_path, "w");
    if (!fp) return 12;
    fprintf(fp, "record,m,j,k,rm,pfm,gr,gt,arr,art,att,n_mu,m_count,nphi,phase_n,is_fr631,norm_before,g,bb_b\n");
    for (int k = -n_mu; k <= n_mu; ++k) {
        fprintf(fp, "rm,-1,-1,%d,%.17g,0,0,0,0,0,0,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g\n",
                k, atm.rm[k], n_mu, m_count, nphi, phase.n, phase.is_fr631,
                phase.norm_before, phase.g_asym, phase.bb_b_ratio);
    }
    for (int m = 0; m < m_count; ++m) {
        for (int j = 0; j <= n_mu; ++j) {
            for (int k = -n_mu; k <= n_mu; ++k) {
                const size_t z = kidx(m, j, k, n_mu);
                fprintf(fp,
                        "kernel,%d,%d,%d,0,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g\n",
                        m, j, k,
                        buf[z], buf[all + z], buf[2u * all + z],
                        buf[3u * all + z], buf[4u * all + z], buf[5u * all + z],
                        n_mu, m_count, nphi, phase.n, phase.is_fr631,
                        phase.norm_before, phase.g_asym, phase.bb_b_ratio);
            }
        }
    }
    fclose(fp);
    free(buf);
    rt_atm_free(&atm);
    rt_value_phase_interp_free(&phase);
    return 0;
}
