#include <stdio.h>
#include "rt_windows_compat.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "rt_io.h"
#include "rt_solver.h"
#include "rt_aerosol_runtime.h"
#include "shared/mie_io.h"

/* Internal per-row record. */
typedef struct {
    int    case_id;
    double sza_deg, vza_deg, raa_deg, wavelength_nm;

    int    have_aod;
    double aod_ref;
    double aod_ref_nm;

    int    have_tau_ref;
    int    have_rho_ref;
    int    have_rho_Q_ref;
    int    have_rho_U_ref;
    double tau_R_ref;
    double rho_I_ref;
    double rho_Q_ref;
    double rho_U_ref;

    int    rc;
    double rho_I_v2;
    double rho_Q_v2;
    double rho_U_v2;
    double tau_R_v2;
    int    n_orders_per_m[3];
    double walltime_ms;
} batch_row_t;

static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1.0e6;
}

static void rstrip(char *s) {
    if (!s) return;
    size_t n = strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r' ||
                     s[n - 1] == ' '  || s[n - 1] == '\t')) {
        s[--n] = '\0';
    }
}

static int count_fields(const char *line) {
    int n = 1;
    for (const char *p = line; *p; ++p) if (*p == ',') n++;
    return n;
}

int rt_io_run_batch(const char *in_csv, const char *out_csv,
                    const rt_options_t *opts,
                    const rt_case_t *cs_template,
                    const rt_batch_aerosol_config_t *aer_cfg) {
    if (!in_csv || !out_csv) return -1;

    rt_options_t opts_default = rt_options_default();
    if (!opts) opts = &opts_default;

    FILE *fp_in = fopen(in_csv, "r");
    if (!fp_in) {
        fprintf(stderr, "rt_io_run_batch: cannot open input '%s'\n", in_csv);
        return -1;
    }

    char header[2048];
    if (!fgets(header, sizeof header, fp_in)) {
        fprintf(stderr, "rt_io_run_batch: empty input file\n");
        fclose(fp_in);
        return -1;
    }
    rstrip(header);
    const int n_cols = count_fields(header);
    const int aerosol_layout = (n_cols == 9) &&
        ((strstr(header, "aod") != NULL) || (strstr(header, "os_I") != NULL) ||
         (strstr(header, "osoaa") != NULL) || (strstr(header, "decoup_I") != NULL));

    if (n_cols < 5 || n_cols > 9 || n_cols == 8) {
        fprintf(stderr, "rt_io_run_batch: header must have 5-7 or 9 columns (got %d in '%s')\n",
                n_cols, header);
        fclose(fp_in);
        return -1;
    }
    if (aerosol_layout && (!aer_cfg || !aer_cfg->mie_path)) {
        fprintf(stderr, "rt_io_run_batch: aerosol CSV requires --mie and aerosol config\n");
        fclose(fp_in);
        return -1;
    }

    const int six_is_rho = (n_cols == 6) && (strstr(header, "rho_I") != NULL);

    int capacity = 256;
    int n_rows = 0;
    batch_row_t *rows = (batch_row_t*)calloc((size_t)capacity, sizeof(batch_row_t));
    if (!rows) { fclose(fp_in); return -1; }

    char line[2048];
    while (fgets(line, sizeof line, fp_in)) {
        rstrip(line);
        if (line[0] == '\0' || line[0] == '#') continue;
        if (n_rows >= capacity) {
            capacity *= 2;
            batch_row_t *tmp = (batch_row_t*)realloc(rows, (size_t)capacity * sizeof(batch_row_t));
            if (!tmp) { free(rows); fclose(fp_in); return -1; }
            rows = tmp;
        }
        batch_row_t *r = &rows[n_rows];
        memset(r, 0, sizeof(*r));
        int parsed = 0;

        if (aerosol_layout) {
            /* case_id,sza_deg,vza_deg,raa_deg,wavelength_nm,aod_ref,rho_I_ref,rho_Q_ref,rho_U_ref */
            parsed = sscanf(line, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                            &r->case_id, &r->sza_deg, &r->vza_deg,
                            &r->raa_deg, &r->wavelength_nm, &r->aod_ref,
                            &r->rho_I_ref, &r->rho_Q_ref, &r->rho_U_ref);
            if (parsed == 9) {
                r->have_aod = 1;
                r->aod_ref_nm = (strstr(header, "aod_865") != NULL) ? 865.0 : 555.0;
                r->have_rho_ref = 1;
                r->have_rho_Q_ref = 1;
                r->have_rho_U_ref = 1;
            }
        } else if (n_cols == 9) {
            parsed = sscanf(line, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                            &r->case_id, &r->sza_deg, &r->vza_deg,
                            &r->raa_deg, &r->wavelength_nm,
                            &r->tau_R_ref, &r->rho_I_ref,
                            &r->rho_Q_ref, &r->rho_U_ref);
            if (parsed == 9) {
                r->have_tau_ref = 1;
                r->have_rho_ref = 1;
                r->have_rho_Q_ref = 1;
                r->have_rho_U_ref = 1;
            }
        } else if (n_cols == 7) {
            parsed = sscanf(line, "%d,%lf,%lf,%lf,%lf,%lf,%lf",
                            &r->case_id, &r->sza_deg, &r->vza_deg,
                            &r->raa_deg, &r->wavelength_nm,
                            &r->tau_R_ref, &r->rho_I_ref);
            if (parsed == 7) { r->have_tau_ref = 1; r->have_rho_ref = 1; }
        } else if (n_cols == 6 && six_is_rho) {
            parsed = sscanf(line, "%d,%lf,%lf,%lf,%lf,%lf",
                            &r->case_id, &r->sza_deg, &r->vza_deg,
                            &r->raa_deg, &r->wavelength_nm, &r->rho_I_ref);
            if (parsed == 6) r->have_rho_ref = 1;
        } else if (n_cols == 6) {
            parsed = sscanf(line, "%d,%lf,%lf,%lf,%lf,%lf",
                            &r->case_id, &r->sza_deg, &r->vza_deg,
                            &r->raa_deg, &r->wavelength_nm, &r->tau_R_ref);
            if (parsed == 6) r->have_tau_ref = 1;
        } else {
            parsed = sscanf(line, "%d,%lf,%lf,%lf,%lf",
                            &r->case_id, &r->sza_deg, &r->vza_deg,
                            &r->raa_deg, &r->wavelength_nm);
        }

        if (parsed < n_cols) {
            fprintf(stderr, "rt_io_run_batch: skip malformed row %d: '%s'\n",
                    n_rows + 1, line);
            continue;
        }
        n_rows++;
    }
    fclose(fp_in);

    if (n_rows == 0) {
        fprintf(stderr, "rt_io_run_batch: no data rows\n");
        free(rows);
        return -1;
    }

    mie_data_t batch_mie = {0};
    int use_batch_aerosol = 0;
    if (aerosol_layout) {
        use_batch_aerosol = 1;
        if (read_mie_file(aer_cfg->mie_path, &batch_mie) != 0) {
            fprintf(stderr, "rt_io_run_batch: cannot read MIE file '%s'\n", aer_cfg->mie_path);
            free(rows);
            return -1;
        }
    }

    int any_minus2 = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 8) reduction(|: any_minus2)
#endif
    for (int i = 0; i < n_rows; ++i) {
        batch_row_t *r = &rows[i];
        rt_case_t cs = {
            .sza_deg       = r->sza_deg,
            .vza_deg       = r->vza_deg,
            .raa_deg       = r->raa_deg,
            .wavelength_nm = r->wavelength_nm,
            .surface       = cs_template ? cs_template->surface       : RT_SURFACE_BLACK,
            .rayleigh_on   = cs_template ? cs_template->rayleigh_on   : 1,
            .aerosol_on    = cs_template ? cs_template->aerosol_on    : 0,
            .n_water           = cs_template ? cs_template->n_water           : 1.34,
            .wind_speed        = cs_template ? cs_template->wind_speed        : 0.0,
            .sigma_type        = cs_template ? cs_template->sigma_type        : 1,
            .q_convention      = cs_template ? cs_template->q_convention      : 1,
            .decouple_sunglint = cs_template ? cs_template->decouple_sunglint : 1,
            .tau_R_input       = r->have_tau_ref ? r->tau_R_ref : 0.0,
        };
        if (use_batch_aerosol) {
            cs.aerosol_on = 1;
        }
        rt_result_t res = {0};
        double t0 = now_ms();
        int rc = 0;

        /* v1.01 LUT mode: pre-allocate grid and use rt_solve_case_pol_lut(). */
        rt_lut_grid_out_t lut = {0};
        double *vza_g = NULL, *raa_g = NULL;
        if (opts->lut_enable && opts->vector_mode) {
            const double vmax  = (opts->lut_vza_max  > 0) ? opts->lut_vza_max  : 85.0;
            const double vstep = (opts->lut_vza_step > 0) ? opts->lut_vza_step : 2.5;
            const double rstep = (opts->lut_raa_step > 0) ? opts->lut_raa_step : 5.0;
            const int n_vza = (int)(vmax / vstep) + 1;
            const int n_raa = (int)(360.0 / rstep);
            vza_g = (double*)malloc(sizeof(double) * n_vza);
            raa_g = (double*)malloc(sizeof(double) * n_raa);
            for (int k = 0; k < n_vza; ++k) vza_g[k] = (double)k * vstep;
            for (int k = 0; k < n_raa; ++k) raa_g[k] = (double)k * rstep;
            lut.n_vza = n_vza; lut.vza_deg = vza_g;
            lut.n_raa = n_raa; lut.raa_deg = raa_g;
            const size_t N = (size_t)n_vza * n_raa;
            lut.rho_I = (double*)calloc(N, sizeof(double));
            lut.rho_Q = (double*)calloc(N, sizeof(double));
            lut.rho_U = (double*)calloc(N, sizeof(double));
            lut.T_diff_dn_dir  = (double*)calloc(N, sizeof(double));
            lut.T_sg_up_dir    = (double*)calloc(N, sizeof(double));
            lut.T_total_up_dir = (double*)calloc(N, sizeof(double));
        }

        if (use_batch_aerosol) {
            rt_aerosol_runtime_options_t ropts = {
                .aer_L_max = aer_cfg->aer_L_max,
                .theta_cut_deg = aer_cfg->theta_cut_deg,
                .delta_m_N = aer_cfg->delta_m_N,
                .apply_nt_tau = aer_cfg->apply_nt_tau,
                .use_loglin_trunc = aer_cfg->use_loglin_trunc,
                .loglin_mu1 = 0.8,
                .loglin_mu2 = 0.94,
                .loglin_threshold = 0.1
            };
            rt_aerosol_input_t aer = {0};
            if (!r->have_aod || rt_aerosol_runtime_prepare(&batch_mie,
                    r->wavelength_nm, r->aod_ref, r->aod_ref_nm, &ropts, &aer, NULL) != 0) {
                rc = -3;
            } else {
                if (opts->lut_enable && opts->vector_mode) {
                    rc = rt_solve_case_pol_lut(&cs, opts, &aer, &res, &lut);
                } else {
                    rc = rt_solve_case_pol_aerosol(&cs, opts, &aer, &res);
                }
                rt_aerosol_runtime_free(&aer);
            }
        } else {
            if (opts->lut_enable && opts->vector_mode) {
                rc = rt_solve_case_pol_lut(&cs, opts, NULL, &res, &lut);
            } else {
                /* 2026-07-14 M1 (Jae 지시): 스칼라 경로 완전 삭제. 항상 벡터. */
                rc = rt_solve_case_pol(&cs, opts, &res);
            }
        }

        /* Write per-case LUT CSV. */
        if (rc == 0 && opts->lut_enable && opts->vector_mode && opts->lut_output_dir) {
            char fpath[1024];
            snprintf(fpath, sizeof(fpath), "%s/case_%d.csv",
                     opts->lut_output_dir, r->case_id);
            FILE *fl = fopen(fpath, "w");
            if (fl) {
                /* v1.01 column naming (at 0+ side = atm side of surface):
                 *   direct_transmittance       = T_dir_dn  (= exp(-τ/μ_sun))
                 *   irradiance_transmittance   = T_dir_dn + T_diff_dn_hemi
                 *                                (= F_BOA(0+) / F_TOA, total flux)
                 *   diffuse_dn_hemi            = T_diff_dn_hemi (diffuse-only flux)
                 *   diffuse_dn_dir(μ_v,φ)      = T_diff_dn_dir   (directional radiance T)
                 *   sunglint_up_flux_ratio_diag = rho_glint * mu_sun (NOT atm T)
                 *   toa_up_flux_ratio_diag      = rho_TOA   * mu_sun (NOT atm T)
                 *
                 * "0+ side" = atm-side at the surface (just above water);
                 * 0- side (just below water) is not produced (black ocean). */
                fprintf(fl, "sza_deg,wavelength_nm,vza_deg,raa_deg,"
                            "rho_I,rho_Q,rho_U,"
                            "direct_transmittance,irradiance_transmittance,"
                            "diffuse_dn_hemi,diffuse_dn_dir,"
                            "sunglint_up_flux_ratio_diag,toa_up_flux_ratio_diag\n");
                const double T_irrad_case = lut.T_dir_dn + lut.T_diff_dn_hemi;
                for (int iv = 0; iv < lut.n_vza; ++iv) {
                    for (int ir = 0; ir < lut.n_raa; ++ir) {
                        const int idx = iv * lut.n_raa + ir;
                        fprintf(fl, "%.6f,%.3f,%.6f,%.6f,"
                                    "%.10e,%.10e,%.10e,"
                                    "%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
                                cs.sza_deg, cs.wavelength_nm,
                                lut.vza_deg[iv], lut.raa_deg[ir],
                                lut.rho_I[idx], lut.rho_Q[idx], lut.rho_U[idx],
                                lut.T_dir_dn, T_irrad_case,
                                lut.T_diff_dn_hemi,
                                lut.T_diff_dn_dir[idx], lut.T_sg_up_dir[idx],
                                lut.T_total_up_dir[idx]);
                    }
                }
                fclose(fl);
            }
        }

        if (opts->lut_enable && opts->vector_mode) {
            free(vza_g); free(raa_g);
            free(lut.rho_I); free(lut.rho_Q); free(lut.rho_U);
            free(lut.T_diff_dn_dir); free(lut.T_sg_up_dir); free(lut.T_total_up_dir);
        }
        double t1 = now_ms();
        r->rc = rc;
        r->rho_I_v2 = res.rho_I;
        r->rho_Q_v2 = res.rho_Q;
        r->rho_U_v2 = res.rho_U;
        r->tau_R_v2 = res.tau_R_total;
        r->walltime_ms = t1 - t0;
        for (int m = 0; m < 3; ++m) r->n_orders_per_m[m] = res.sos_orders_per_m[m];
        if (rc == -2 || rc == -3) any_minus2 = 1;
    }

    FILE *fp_out = fopen(out_csv, "w");
    if (!fp_out) {
        fprintf(stderr, "rt_io_run_batch: cannot open output '%s'\n", out_csv);
        if (use_batch_aerosol) mie_data_free(&batch_mie);
        free(rows);
        return -1;
    }

    if (aerosol_layout) {
        fprintf(fp_out, "case_id,sza_deg,vza_deg,raa_deg,wavelength_nm,aod_ref,");
    } else {
        fprintf(fp_out, "case_id,sza_deg,vza_deg,raa_deg,wavelength_nm,");
    }

    /* v1.01: output_mode_debug controls column set.
     *   SIMPLE (default): rho_I_v2, rho_Q_v2, rho_U_v2 only
     *                     (transmittance columns will be appended here in Step 3)
     *   DEBUG: full v0.9 column set with ref, delta, n_orders, walltime */
    const int debug_mode = opts->output_mode_debug;
    const int have_QU_ref = (n_rows > 0) && rows[0].have_rho_Q_ref && rows[0].have_rho_U_ref;

    if (debug_mode) {
        fprintf(fp_out,
            "rho_I_v2,rho_I_ref,delta_abs,delta_rel_pct,delta_rel_abs_pct,"
            "tau_R_v2,tau_R_ref,n_orders_m0,n_orders_m1,n_orders_m2,walltime_ms");
        if (opts->vector_mode && have_QU_ref) {
            fprintf(fp_out,
                ",rho_Q_v2,rho_Q_ref,delta_abs_Q,delta_rel_pct_Q,delta_rel_abs_pct_Q"
                ",rho_U_v2,rho_U_ref,delta_abs_U,delta_rel_pct_U,delta_rel_abs_pct_U");
        } else if (opts->vector_mode) {
            fprintf(fp_out, ",rho_Q_v2,rho_U_v2");
        }
    } else {
        /* SIMPLE mode: IQU reflectances only. */
        fprintf(fp_out, "rho_I_v2");
        if (opts->vector_mode) fprintf(fp_out, ",rho_Q_v2,rho_U_v2");
    }
    fprintf(fp_out, "\n");

    for (int i = 0; i < n_rows; ++i) {
        const batch_row_t *r = &rows[i];
        if (aerosol_layout) {
            fprintf(fp_out, "%d,%.6f,%.6f,%.6f,%.3f,%.10e,",
                    r->case_id, r->sza_deg, r->vza_deg, r->raa_deg,
                    r->wavelength_nm, r->aod_ref);
        } else {
            fprintf(fp_out, "%d,%.6f,%.6f,%.6f,%.3f,",
                    r->case_id, r->sza_deg, r->vza_deg, r->raa_deg,
                    r->wavelength_nm);
        }

        if (!debug_mode) {
            /* SIMPLE mode: just IQU reflectances. */
            fprintf(fp_out, "%.10e", r->rho_I_v2);
            if (opts->vector_mode) fprintf(fp_out, ",%.10e,%.10e", r->rho_Q_v2, r->rho_U_v2);
            fprintf(fp_out, "\n");
            continue;
        }

        /* DEBUG mode: full v0.9 column set. */
        fprintf(fp_out, "%.10e,", r->rho_I_v2);
        if (r->have_rho_ref) {
            const double da = r->rho_I_v2 - r->rho_I_ref;
            const double drs = (r->rho_I_ref != 0.0) ? 100.0 * da / r->rho_I_ref : 0.0;
            fprintf(fp_out, "%.10e,%.6e,%.6e,%.6e,",
                    r->rho_I_ref, da, drs, fabs(drs));
        } else {
            fprintf(fp_out, ",,,,");
        }
        fprintf(fp_out, "%.10e,", r->tau_R_v2);
        if (r->have_tau_ref) fprintf(fp_out, "%.10e,", r->tau_R_ref);
        else                 fprintf(fp_out, ",");
        fprintf(fp_out, "%d,%d,%d,%.3f",
                r->n_orders_per_m[0], r->n_orders_per_m[1],
                r->n_orders_per_m[2], r->walltime_ms);

        if (opts->vector_mode && have_QU_ref) {
            const double Iref_abs = fabs(r->rho_I_ref);
            const double daQ = r->rho_Q_v2 - r->rho_Q_ref;
            const double drsQ = (Iref_abs > 0.0) ? 100.0 * daQ / Iref_abs : 0.0;
            const double daU = r->rho_U_v2 - r->rho_U_ref;
            const double drsU = (Iref_abs > 0.0) ? 100.0 * daU / Iref_abs : 0.0;
            fprintf(fp_out,
                    ",%.10e,%.10e,%.6e,%.6e,%.6e"
                    ",%.10e,%.10e,%.6e,%.6e,%.6e",
                    r->rho_Q_v2, r->rho_Q_ref, daQ, drsQ, fabs(drsQ),
                    r->rho_U_v2, r->rho_U_ref, daU, drsU, fabs(drsU));
        } else if (opts->vector_mode) {
            fprintf(fp_out, ",%.10e,%.10e", r->rho_Q_v2, r->rho_U_v2);
        }
        fprintf(fp_out, "\n");
    }
    fclose(fp_out);

    int n_compared = 0;
    double max_drabs = 0.0;
    for (int i = 0; i < n_rows; ++i) {
        if (!rows[i].have_rho_ref || rows[i].rho_I_ref == 0.0) continue;
        n_compared++;
        const double dra = fabs(100.0 * (rows[i].rho_I_v2 - rows[i].rho_I_ref) / rows[i].rho_I_ref);
        if (dra > max_drabs) max_drabs = dra;
    }
    printf("rt_io_run_batch: n=%d, compared=%d, max|drel_I|=%.4f%%, output=%s\n",
           n_rows, n_compared, max_drabs, out_csv);

    if (use_batch_aerosol) mie_data_free(&batch_mie);
    free(rows);
    return any_minus2 ? -2 : 0;
}

void rt_io_print_single(const rt_case_t *cs, const rt_result_t *res) {
    const double pi = 3.14159265358979323846;
    if (cs->surface == RT_SURFACE_OCEAN) {
        const double ed0p = res->Ed_0plus_air;
        const double ed0m = res->Ed_0minus_water;
        const double RrsQ = res->R_rs_0plus_Q;
        const double RrsU = res->R_rs_0plus_U;
        const double rrsQ = res->r_rs_0minus_Q;
        const double rrsU = res->r_rs_0minus_U;
        printf("%.10e %.10e %.10e  "
               "TOA_rho_I=%.10e TOA_rho_Q=%.10e TOA_rho_U=%.10e "
               "TOA_rho_atm_I=%.10e TOA_rho_atm_Q=%.10e TOA_rho_atm_U=%.10e "
               "TOA_rho_water_direct_I=%.10e TOA_rho_water_direct_Q=%.10e TOA_rho_water_direct_U=%.10e "
               "TOA_rho_water_sky_I=%.10e TOA_rho_water_sky_Q=%.10e TOA_rho_water_sky_U=%.10e "
               "TOA_rho_water_total_I=%.10e TOA_rho_water_total_Q=%.10e TOA_rho_water_total_U=%.10e "
               "TOA_rho_glint_direct_I=%.10e TOA_rho_glint_direct_Q=%.10e TOA_rho_glint_direct_U=%.10e "
               "TOA_rho_glint_on_I=%.10e TOA_rho_glint_on_Q=%.10e TOA_rho_glint_on_U=%.10e "
               "T_dir_dn=%.10e T_diff_dn_hemi=%.10e T_total_dn_hemi=%.10e T_diff_dn_dir=%.10e "
               "T_dir_up_view=%.10e T_diff_up_view=%.10e T_total_up_view=%.10e TOA_water_signal_I=%.10e T_up_rt_valid=%d "
               "AOD_ref=%.10e AOD_ref_nm=%.3f AOD_band=%.10e AOD_ext_ratio=%.10e "
               "Lu0plus=%.6e Qu0plus=%.6e Uu0plus=%.6e Ed0plus=%.6e "
               "Ed0plus_direct=%.6e Ed0plus_diffuse=%.6e "
               "Rrs0plus_I=%.6e Rrs0plus_Q=%.6e Rrs0plus_U=%.6e "
               "BRF0plus_I=%.6e BRF0plus_Q=%.6e BRF0plus_U=%.6e "
               "Lu0minus=%.6e Qu0minus=%.6e Uu0minus=%.6e Ed0minus=%.6e Eu0minus=%.6e "
               "Ed0minus_direct=%.6e Ed0minus_diffuse=%.6e "
               "Eu0minus_direct=%.6e Eu0minus_diffuse=%.6e "
               "rrs0minus_I=%.6e rrs0minus_Q=%.6e rrs0minus_U=%.6e "
               "BRF0minus_I=%.6e BRF0minus_Q=%.6e BRF0minus_U=%.6e "
               "Kd0minus=%.6f Ku0minus=%.6f "
               "a_w=%.8e b_w=%.8e bb_w=%.8e a_dom=%.8e a_chl=%.8e "
               "a_pig=%.8e b_pig=%.8e bb_pig=%.8e a_min=%.8e b_min=%.8e bb_min=%.8e "
               "a_total=%.8e b_total=%.8e bb_total=%.8e omega_total=%.6f "
               "sza=%.3f vza=%.3f raa=%.3f wl=%.3f orders=%d conv=%d\n",
               res->rho_I, res->rho_Q, res->rho_U,
               res->rho_I, res->rho_Q, res->rho_U,
               res->rho_atm_path_I, res->rho_atm_path_Q, res->rho_atm_path_U,
               res->rho_water_direct_I, res->rho_water_direct_Q, res->rho_water_direct_U,
               res->rho_water_sky_I, res->rho_water_sky_Q, res->rho_water_sky_U,
               res->rho_water_total_I, res->rho_water_total_Q, res->rho_water_total_U,
               res->rho_glint_direct_I, res->rho_glint_direct_Q, res->rho_glint_direct_U,
               res->rho_I_glint, res->rho_Q_glint, res->rho_U_glint,
               res->T_dir_dn, res->T_diff_dn_hemi, res->T_total_dn_hemi, res->T_diff_dn_dir,
               res->T_dir_up_view, res->T_diff_up_view, res->T_total_up_view, res->I_TOA_water_signal, res->T_up_rt_valid,
               res->aerosol_aod_ref, res->aerosol_aod_ref_nm, res->aerosol_aod_band, res->aerosol_extinction_ratio,
               res->Lu_0plus_view, res->Qu_0plus_view, res->Uu_0plus_view, res->Ed_0plus_air,
               res->Ed_0plus_direct, res->Ed_0plus_diffuse,
               res->R_rs_0plus, RrsQ, RrsU, pi*res->R_rs_0plus, pi*RrsQ, pi*RrsU,
               res->Lu_0minus_view, res->Qu_0minus_view, res->Uu_0minus_view,
               res->Ed_0minus_water, res->Eu_0minus_water,
               res->Ed_0minus_direct, res->Ed_0minus_diffuse,
               res->Eu_0minus_direct, res->Eu_0minus_diffuse,
               res->r_rs_0minus, rrsQ, rrsU, pi*res->r_rs_0minus, pi*rrsQ, pi*rrsU,
               res->Kd_0minus, res->Ku_0minus,
               res->a_w_used, res->b_w_used, res->bb_w_used, res->a_cdom_used,
               res->a_chl_used,
               res->a_pig_used, res->b_pig_used, res->bb_pig_used,
               res->a_min_used, res->b_min_used, res->bb_min_used,
               res->a_total_used, res->b_total_used, res->bb_total_used, res->omega_water,
               cs->sza_deg, cs->vza_deg, cs->raa_deg, cs->wavelength_nm,
               res->n_orders_used, res->converged);
        return;
    }
    printf("%.10e %.10e %.10e  "
           "TOA_rho_I=%.10e TOA_rho_Q=%.10e TOA_rho_U=%.10e "
           "TOA_rho_glint_direct_I=%.10e TOA_rho_glint_direct_Q=%.10e TOA_rho_glint_direct_U=%.10e "
           "T_dir_dn=%.10e T_diff_dn_hemi=%.10e T_total_dn_hemi=%.10e T_diff_dn_dir=%.10e "
           "T_dir_up_view=%.10e T_diff_up_view=%.10e T_total_up_view=%.10e TOA_water_signal_I=%.10e T_up_rt_valid=%d "
           "diag_sunglint_up_flux_ratio=%.10e diag_TOA_up_flux_ratio=%.10e "
           "AOD_ref=%.10e AOD_ref_nm=%.3f AOD_band=%.10e AOD_ext_ratio=%.10e "
           "tau_R=%.10e sza=%.3f vza=%.3f raa=%.3f wl=%.3f orders=%d conv=%d\n",
           res->rho_I, res->rho_Q, res->rho_U,
           res->rho_I, res->rho_Q, res->rho_U,
           res->rho_glint_direct_I, res->rho_glint_direct_Q, res->rho_glint_direct_U,
           res->T_dir_dn, res->T_diff_dn_hemi, res->T_total_dn_hemi, res->T_diff_dn_dir,
           res->T_dir_up_view, res->T_diff_up_view, res->T_total_up_view, res->I_TOA_water_signal, res->T_up_rt_valid,
           res->T_sg_up_dir, res->T_total_up_dir,
           res->aerosol_aod_ref, res->aerosol_aod_ref_nm, res->aerosol_aod_band, res->aerosol_extinction_ratio,
           res->tau_R_total,
           cs->sza_deg, cs->vza_deg, cs->raa_deg, cs->wavelength_nm,
           res->n_orders_used, res->converged);
}
