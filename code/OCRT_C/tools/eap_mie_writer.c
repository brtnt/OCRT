#define _USE_MATH_DEFINES
/* eap_mie_writer.c -- generates legacy-format .mie files for the OCRT
 * water_iop path.  This is a build-time tool, not part of the public API:
 * rt_eap_mie_phase_compute() deliberately produces no files and no a/b/bb.
 *
 * usage: eap_mie_writer <species_id> <Deff_um> <out.mie>
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "../src/internal/rt_eap_species_catalog.h"
#include "../src/internal/rt_eap_coated_mie.h"

#define NMEDIA   1.334
#define VEFF     0.6
#define STEP_UM  0.05
#define MAXMULT  5.5
#define CI_MG_M3 2.0e6   /* 세포내 엽록소 밀도 2 kg/m3 */
#define NSPEC    101      /* 0.350..0.850 um, 5 nm */
#define NANGLE   361      /* 180..0 deg, 0.5 deg */
#define NPHASE   12
static const double PHASE_WL[NPHASE] = {350,400,443,470,490,510,555,620,660,680,745,850};

int main(int argc, char **argv)
{
    if (argc != 4) { fprintf(stderr, "usage: %s <species_id> <Deff_um> <out.mie>\n", argv[0]); return 2; }
    unsigned id = (unsigned)atoi(argv[1]);
    double deff = atof(argv[2]);
    const rt_eap_species_entry_t *sp = rt_eap_catalog_entry(id);
    if (!sp || !(deff > 0.0)) { fprintf(stderr, "bad species or Deff\n"); return 2; }
    const double vs = *sp->vs, fr = pow(1.0 - vs, 1.0 / 3.0);
    const size_t nd = (size_t)(MAXMULT * deff / STEP_UM + 0.5);

    double mu[NANGLE]; int col[NANGLE]; double th[NANGLE];
    for (int j = 0; j < NANGLE; ++j) {
        th[j] = 180.0 - 0.5 * j;                       /* descending, legacy order */
        double t = th[j] > 90.0 ? 180.0 - th[j] : th[j];
        col[j] = th[j] > 90.0 ? 1 : 0;
        mu[j] = cos(t * M_PI / 180.0);
    }
    double *rad_w = calloc(nd + 1, sizeof *rad_w);
    double *rad_area = calloc(nd + 1, sizeof *rad_area);
    double *rad_qext = calloc(nd + 1, sizeof *rad_qext);
    double *rad_qsca = calloc(nd + 1, sizeof *rad_qsca);
    double *rad_gqsc = calloc(nd + 1, sizeof *rad_gqsc);
    double *rad_civol = calloc(nd + 1, sizeof *rad_civol);
    int *rad_rc = calloc(nd + 1, sizeof *rad_rc);
    double *rad11 = malloc((nd + 1) * NANGLE * sizeof *rad11);
    double *rad12 = malloc((nd + 1) * NANGLE * sizeof *rad12);
    double *rad33 = malloc((nd + 1) * NANGLE * sizeof *rad33);
    if (!rad_w || !rad_area || !rad_qext || !rad_qsca || !rad_gqsc ||
        !rad_civol || !rad_rc || !rad11 || !rad12 || !rad33) return 6;
    static double P11[NPHASE][NANGLE], P12[NPHASE][NANGLE], P33[NPHASE][NANGLE];
    static double sQe[NSPEC], sQb[NSPEC], sG[NSPEC], sC[NSPEC], sS[NSPEC];

    /* Spectral block and phase block share the same size loop per wavelength.
     * Each radius is independent.  We calculate radius contributions in
     * parallel, then sum them in the original increasing-radius order so the
     * generated file remains bit-identical to the serial reference. */
    for (int is = -NPHASE; is < NSPEC; ++is) {
        int phase_idx = is < 0 ? is + NPHASE : -1;
        double wl = phase_idx >= 0 ? PHASE_WL[phase_idx] : 350.0 + 5.0 * is;
        double wl_ri = wl < rt_eap_catalog_wl_lo() ? rt_eap_catalog_wl_lo() : wl;
        double nsh, ksh, nco, kco;
        if (rt_eap_catalog_ri(id, wl_ri, &nsh, &ksh, &nco, &kco) != 0) return 3;
        double wvno = 2.0 * M_PI * NMEDIA / (wl * 1.0e-3);

        memset(rad_w, 0, (nd + 1) * sizeof *rad_w);
        memset(rad_area, 0, (nd + 1) * sizeof *rad_area);
        memset(rad_qext, 0, (nd + 1) * sizeof *rad_qext);
        memset(rad_qsca, 0, (nd + 1) * sizeof *rad_qsca);
        memset(rad_gqsc, 0, (nd + 1) * sizeof *rad_gqsc);
        memset(rad_civol, 0, (nd + 1) * sizeof *rad_civol);
        memset(rad_rc, 0, (nd + 1) * sizeof *rad_rc);

        int allocation_failed = 0;
        #pragma omp parallel
        {
            double *m1 = phase_idx >= 0 ? malloc(2 * NANGLE * sizeof *m1) : NULL;
            double *m2 = phase_idx >= 0 ? malloc(2 * NANGLE * sizeof *m2) : NULL;
            double *s21 = phase_idx >= 0 ? malloc(2 * NANGLE * sizeof *s21) : NULL;
            double *d21 = phase_idx >= 0 ? malloc(2 * NANGLE * sizeof *d21) : NULL;
            if (phase_idx >= 0 && (!m1 || !m2 || !s21 || !d21)) {
                #pragma omp atomic write
                allocation_failed = 1;
            }
            #pragma omp for schedule(static)
            for (long kk = 1; kk <= (long)nd; ++kk) {
                size_t k = (size_t)kk;
                if (allocation_failed) { rad_rc[k] = RT_EAP_MIE_ENOMEM; continue; }
                double d = (double)k * STEP_UM, r = 0.5 * d;
                double w = pow(r, (1.0 - 3.0 * VEFF) / VEFF) * exp(-r / (0.5 * deff * VEFF));
                if (!(w > 0.0) || !isfinite(w)) continue;
                rt_eap_mie_eff_t e;
                int mie_rc;
                if (phase_idx >= 0) {
                    mie_rc = rt_eap_coated_mie(r * fr, r, wvno,
                                               nsh - ksh * _Complex_I,
                                               nco - kco * _Complex_I,
                                               mu, NANGLE, m1, m2, s21, d21, &e);
                } else {
                    mie_rc = rt_eap_coated_mie(r * fr, r, wvno,
                                               nsh - ksh * _Complex_I,
                                               nco - kco * _Complex_I,
                                               NULL, 0, NULL, NULL, NULL, NULL, &e);
                }
                rad_rc[k] = mie_rc;
                if (mie_rc != RT_EAP_MIE_OK) continue;
                double rm = r * 1.0e-6;
                double area = w * rm * rm;
                rad_w[k] = w;
                rad_area[k] = area;
                rad_qext[k] = e.qext * area;
                rad_qsca[k] = e.qsca * area;
                rad_gqsc[k] = e.gqsc * area;
                rad_civol[k] = w * (4.0 / 3.0) * M_PI * rm * rm * rm;
                if (phase_idx >= 0) {
                    double *o11 = rad11 + k * NANGLE;
                    double *o12 = rad12 + k * NANGLE;
                    double *o33 = rad33 + k * NANGLE;
                    for (int j = 0; j < NANGLE; ++j) {
                        size_t q = (size_t)j * 2 + (size_t)col[j];
                        o11[j] = w * 0.5 * (m1[q] + m2[q]);
                        o12[j] = w * 0.5 * (m2[q] - m1[q]);
                        o33[j] = w * s21[q];
                    }
                }
            }
            free(m1); free(m2); free(s21); free(d21);
        }
        if (allocation_failed) return 6;

        double sumw = 0, sQext = 0, sQsca = 0, sGq = 0, civol = 0;
        double acc11[NANGLE] = {0}, acc12[NANGLE] = {0}, acc33[NANGLE] = {0};
        for (size_t k = 1; k <= nd; ++k) {
            if (rad_rc[k] != RT_EAP_MIE_OK) return 4;
            if (!(rad_w[k] > 0.0)) continue;
            sumw += rad_area[k];
            sQext += rad_qext[k];
            sQsca += rad_qsca[k];
            sGq += rad_gqsc[k];
            civol += rad_civol[k];
            if (phase_idx >= 0) {
                const double *i11 = rad11 + k * NANGLE;
                const double *i12 = rad12 + k * NANGLE;
                const double *i33 = rad33 + k * NANGLE;
                for (int j = 0; j < NANGLE; ++j) {
                    acc11[j] += i11[j];
                    acc12[j] += i12[j];
                    acc33[j] += i33[j];
                }
            }
        }
        if (phase_idx >= 0) {
            double s = 0;
            for (int j = NANGLE - 1; j > 0; --j) {
                double t0 = th[j] * M_PI / 180.0, t1 = th[j - 1] * M_PI / 180.0;
                s += 0.5 * (t1 - t0) * (acc11[j] * sin(t0) + acc11[j - 1] * sin(t1));
            }
            s *= 0.5;
            for (int j = 0; j < NANGLE; ++j) {
                P11[phase_idx][j] = acc11[j] / s;
                P12[phase_idx][j] = acc12[j] / s;
                P33[phase_idx][j] = acc33[j] / s;
            }
        } else {
            double sc = 1.0 / (civol * CI_MG_M3);
            sQe[is] = sQext / sumw; sQb[is] = sQsca / sumw;
            sG[is] = sGq / sQsca;
            sC[is] = M_PI * sQext * sc;
            sS[is] = M_PI * sQsca * sc;
        }
    }
    free(rad_w); free(rad_area); free(rad_qext); free(rad_qsca);
    free(rad_gqsc); free(rad_civol); free(rad_rc);
    free(rad11); free(rad12); free(rad33);
    FILE *f = fopen(argv[3], "w");
    if (!f) return 5;
    fprintf(f, " %d\n", NANGLE);
    fprintf(f, "   Wlgth  Nor_Ext_Co  Nor_Sca_Co  Sg_Sca_Alb  Asymm_Para  Extinct_Co  Scatter_Co  (%s)\n", sp->name);
    double c550 = 0;
    for (int i = 0; i < NSPEC; ++i) if (fabs(350.0 + 5.0 * i - 550.0) < 1e-9) c550 = sC[i];
    for (int i = 0; i < NSPEC; ++i)
        fprintf(f, "    %.4f     %.4f        %.4f        %.4f        %.4f     %.4E  %.4E\n",
                (350.0 + 5.0 * i) * 1e-3, sC[i] / c550, sS[i] / c550,
                sQb[i] / sQe[i], sG[i], sC[i], sS[i]);
    const char *hdr[3] = {"Phase Function (P11)", "Phase Function (P12 / Q-polarization)",
                          "Phase Function (P33 / U-polarization)"};
    for (int b = 0; b < 3; ++b) {
        fprintf(f, "%s\n   TETA ", hdr[b]);
        for (int k = 0; k < NPHASE; ++k) fprintf(f, "    %.4f", PHASE_WL[k] * 1e-3);
        fprintf(f, "\n");
        for (int j = 0; j < NANGLE; ++j) {
            fprintf(f, "  %6.2f", th[j]);
            for (int k = 0; k < NPHASE; ++k) {
                double v = b == 0 ? P11[k][j] : (b == 1 ? P12[k][j] : P33[k][j]);
                fprintf(f, "  %+.4E", v);
            }
            fprintf(f, "\n");
        }
        if (b < 2) fprintf(f, "\n");
    }
    fprintf(f, "\n");
    fclose(f);
    return 0;
}
