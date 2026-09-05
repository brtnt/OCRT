#include "rt_value_phase.h"
#include "shared/mie_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int parse_bool(const char *s)
{
    return strcmp(s, "1") == 0 || strcmp(s, "true") == 0 || strcmp(s, "on") == 0;
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr, "usage: %s MIE_PATH WL_NM TRUNCATE_0_OR_1 OUT_CSV\n", argv[0]);
        return 2;
    }
    const char *mie_path = argv[1];
    const double wl_nm = strtod(argv[2], NULL);
    const int truncate = parse_bool(argv[3]);
    const char *out_path = argv[4];
    if (!(wl_nm > 0.0)) return 3;

    mie_data_t mie;
    memset(&mie, 0, sizeof mie);
    if (read_mie_file(mie_path, &mie) != 0) {
        fprintf(stderr, "failed to read Mie file: %s\n", mie_path);
        return 4;
    }
    double *p11 = (double *)malloc((size_t)mie.n_ang * sizeof(*p11));
    double *p12 = (double *)malloc((size_t)mie.n_ang * sizeof(*p12));
    double *p33 = (double *)malloc((size_t)mie.n_ang * sizeof(*p33));
    if (!p11 || !p12 || !p33) return 5;
    if (mie_phase_nodes_at_wavelength_linear(
            &mie, wl_nm / 1000.0, p11, p12, p33) != 0) {
        fprintf(stderr, "phase wavelength interpolation failed\n");
        return 6;
    }

    rt_value_phase_interp_t raw = {0}, out = {0};
    if (rt_value_phase_interp_build(
            &raw, mie.angles, p11, p12, p33, mie.n_ang, 0) != 0) {
        fprintf(stderr, "value-phase build failed\n");
        return 7;
    }
    free(p11); free(p12); free(p33);
    mie_data_free(&mie);

    double A = 0.0;
    const rt_value_phase_interp_t *selected = &raw;
    if (truncate) {
        if (rt_value_phase_interp_loglinear_truncate(
                &raw, 0.85, 0.92, 0.1, &out, &A) != 0) {
            fprintf(stderr, "truncation failed\n");
            return 8;
        }
        selected = &out;
    }

    FILE *fp = fopen(out_path, "w");
    if (!fp) return 9;
    fprintf(fp,
            "theta_deg,p11,p12,p33,A,norm_before,g,bb_b,is_fr631,n,truncate\n");
    for (int i = 0; i < selected->n; ++i) {
        fprintf(fp,
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d\n",
                selected->theta_deg[i], selected->p11[i], selected->p12[i],
                selected->p33[i], A, selected->norm_before,
                selected->g_asym, selected->bb_b_ratio,
                selected->is_fr631, selected->n, truncate);
    }
    fclose(fp);
    rt_value_phase_interp_free(&raw);
    rt_value_phase_interp_free(&out);
    return 0;
}
