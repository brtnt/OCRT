/* =============================================================================
 * rt_absorption.c — Trace-gas absorption optical depth (v1.01+, xsec-LUT based)
 * =============================================================================
 *
 * PHYSICS
 *   Vertical gas absorption optical depth by Beer-Lambert over a layered
 *   plane-parallel atmosphere:
 *
 *     tau_abs(lambda) = sum_gas sum_layers sigma_gas(lambda, z_layer) * N_layer
 *
 *   sigma [cm^2/molecule] comes from pre-tabulated cross-section LUTs
 *   (HITRAN LBL degraded to the working grid; for O3/NO2 the tables are the
 *   CONTINUUM-MERGED versions — Serdyuchenko-2014 / Bogumil-2000 continua
 *   merged over the LBL — the LBL-only tables have near-zero visible O3/NO2,
 *   see header check in the verification guideline).  N [molecules/cm^2]
 *   comes from AFGL standard-atmosphere mixing-ratio profiles.
 *
 * GRIDS / CONVENTIONS
 *   AFGL z_km ascending (0 = ground .. TOA); number density
 *   n_gas(z) = mr_ppmv(z) * 1e-6 * n_total(z) [cm^-3]; all column integrals
 *   are trapezoids on the level grid; all interpolation linear (project
 *   rule: never nearest-neighbor).  Unit conversions: dz km -> cm via 1e5;
 *   Dobson: 1 DU = 2.6867e16 molecules/cm^2 (O3 usstd76 column 9.321e18
 *   = 346.9 DU is a Tier-1 golden).
 *
 * COLUMN OVERRIDE
 *   A user-supplied total column (e.g. --o3-du) rescales every layer of
 *   that gas PROPORTIONALLY (profile shape preserved, total matched).
 * ========================================================================== */

#include "rt_absorption.h"
#include "rt_spectral_contract.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const char *RT_GAS_NAMES[RT_N_GAS] = {"h2o", "o3", "o2", "co2", "no2", "ch4"};

static const char *AFGL_NAMES[7] = {
    "usstd76", "tropical", "mlsumm", "mlwint", "sasumm", "sawint", "userdef"
};

/* === AFGL atm loader === */

int rt_afgl_load(const char *afgl_dir, rt_afgl_profile_t profile,
                 rt_afgl_atm_t *atm) {
    if (!atm || profile < 0 || profile > 6) return -1;
    memset(atm, 0, sizeof(*atm));

    char path[1024];
    snprintf(path, sizeof(path), "%s/afgl_%s.dat", afgl_dir, AFGL_NAMES[profile]);

    FILE *f = fopen(path, "r");
    if (!f) {
        if (profile == RT_AFGL_USERDEF) {
            fprintf(stderr,
                "rt_afgl_load: --atm-profile userdef requires %s\n"
                "  Template (빈 껍데기) 은 다음 turn에 제공될 예정.\n"
                "  현재는 inputs/afgl_atm/afgl_userdef.dat 파일을 수동으로\n"
                "  작성해야 합니다 (column 형식: usstd76 파일과 동일).\n",
                path);
        } else {
            fprintf(stderr, "rt_afgl_load: cannot open %s\n", path);
        }
        return -2;
    }

    /* Two-pass: count lines first */
    char line[1024];
    int n_data = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
        n_data++;
    }
    rewind(f);

    atm->n_levels = n_data;
    atm->z_km        = (double*)calloc(n_data, sizeof(double));
    atm->P_mbar      = (double*)calloc(n_data, sizeof(double));
    atm->T_K          = (double*)calloc(n_data, sizeof(double));
    atm->n_total_cm3  = (double*)calloc(n_data, sizeof(double));
    for (int g = 0; g < RT_N_GAS; ++g) {
        atm->mr[g] = (double*)calloc(n_data, sizeof(double));
    }

    /* Read data: z, P, T, n_total, mr_h2o, mr_o3, mr_o2, mr_co2, mr_no2, mr_ch4 */
    int idx = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
        double z, P, T, nt, h2o, o3, o2, co2, no2, ch4;
        int nread = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                           &z, &P, &T, &nt, &h2o, &o3, &o2, &co2, &no2, &ch4);
        if (nread != 10) continue;
        atm->z_km[idx] = z;
        atm->P_mbar[idx] = P;
        atm->T_K[idx] = T;
        atm->n_total_cm3[idx] = nt;
        atm->mr[RT_GAS_H2O][idx] = h2o;
        atm->mr[RT_GAS_O3][idx]  = o3;
        atm->mr[RT_GAS_O2][idx]  = o2;
        atm->mr[RT_GAS_CO2][idx] = co2;
        atm->mr[RT_GAS_NO2][idx] = no2;
        atm->mr[RT_GAS_CH4][idx] = ch4;
        idx++;
    }
    fclose(f);

    /* [v1.04+qaa_debug, 2026-05-24] 사용자 정의 profile sanity check:
     * n_total[k] vs ideal gas law P/(k_B T) 일관성 검증.
     *   n [cm⁻³] = P [mbar] × 7.2429e18 / T [K]
     *   (k_B = 1.380649e-23 J/K; 1 mbar = 100 Pa = 100 N/m²)
     * 1% 이상 차이 시 stderr 경고 (사용자 책임으로 둠, 결정 3).
     * AFGL standard profile 은 보통 < 0.1 % 일치.
     * T <= 0 또는 P <= 0 같은 비물리적 값도 함께 경고. */
    {
        const double K_B_SCALE = 7.2429e+18;   /* mbar·cm³·K⁻¹ */
        const double THRESH_REL = 0.01;        /* 1% relative tolerance */
        int n_warn = 0;
        const int MAX_WARN = 5;                /* 폭주 방지 */
        for (int k = 0; k < atm->n_levels; ++k) {
            double T = atm->T_K[k];
            double P = atm->P_mbar[k];
            double nt = atm->n_total_cm3[k];
            if (T <= 0.0 || P <= 0.0) {
                if (n_warn < MAX_WARN) {
                    fprintf(stderr,
                        "warning: %s layer %d (z=%.2f km): non-physical "
                        "T=%.3g K or P=%.3g mbar — OCRT 동작 보장 안 됨\n",
                        AFGL_NAMES[profile], k, atm->z_km[k], T, P);
                }
                n_warn++;
                continue;
            }
            double n_check = P * K_B_SCALE / T;
            if (n_check <= 0.0) continue;
            double rel = (nt - n_check) / n_check;
            if (fabs(rel) > THRESH_REL) {
                if (n_warn < MAX_WARN) {
                    fprintf(stderr,
                        "warning: %s layer %d (z=%.2f km): n_cm-3=%.3e "
                        "inconsistent with ideal-gas P/(k_B T)=%.3e "
                        "(rel diff %+.2f%%, threshold ±%.1f%%)\n",
                        AFGL_NAMES[profile], k, atm->z_km[k], nt, n_check,
                        rel * 100.0, THRESH_REL * 100.0);
                }
                n_warn++;
            }
        }
        if (n_warn > MAX_WARN) {
            fprintf(stderr,
                "warning: %s — %d additional sanity-check warnings suppressed\n",
                AFGL_NAMES[profile], n_warn - MAX_WARN);
        }
    }

    /* Compute default column abundances by trapezoid integration:
     *   N_gas_column [cm^-2] = ∫ n_gas(z) dz = ∫ mr(z)·1e-6 · n_total(z) dz
     *   dz in cm (z in km × 1e5)
     */
    for (int g = 0; g < RT_N_GAS; ++g) {
        double col = 0.0;
        for (int i = 0; i < atm->n_levels - 1; ++i) {
            double n1 = atm->mr[g][i]   * 1e-6 * atm->n_total_cm3[i];
            double n2 = atm->mr[g][i+1] * 1e-6 * atm->n_total_cm3[i+1];
            double dz_cm = (atm->z_km[i+1] - atm->z_km[i]) * 1.0e5;
            col += 0.5 * (n1 + n2) * dz_cm;
        }
        atm->column_default[g] = col;
    }

    return 0;
}

void rt_afgl_free(rt_afgl_atm_t *atm) {
    if (!atm) return;
    free(atm->z_km);
    free(atm->P_mbar);
    free(atm->T_K);
    free(atm->n_total_cm3);
    for (int g = 0; g < RT_N_GAS; ++g) free(atm->mr[g]);
    memset(atm, 0, sizeof(*atm));
}

/* === HITRAN cross section loader === */

static int xsec_grid_valid_for_ocrt(const rt_xsec_t *xsec)
{
    if (!xsec || xsec->n_wl < 2 || !xsec->wl_nm) return 0;
    for (int i = 0; i < xsec->n_wl; ++i) {
        if (!isfinite(xsec->wl_nm[i]) ||
            (i > 0 && !(xsec->wl_nm[i] > xsec->wl_nm[i - 1]))) return 0;
    }
    return rt_spectral_grid_covers(xsec->wl_nm[0],
                                    xsec->wl_nm[xsec->n_wl - 1]);
}

int rt_xsec_load(const char *xsec_dir, int gas_idx, rt_xsec_t *xsec) {
    if (!xsec || gas_idx < 0 || gas_idx >= RT_N_GAS) return -1;
    memset(xsec, 0, sizeof(*xsec));

    char path[1024];
    snprintf(path, sizeof(path), "%s/xsec_%s.dat", xsec_dir, RT_GAS_NAMES[gas_idx]);

    FILE *f = fopen(path, "r");
    if (!f) {
        /* Silent fallback: zero absorption for this gas */
        xsec->available = 0;
        return 1;
    }

    /* First pass: parse header for n_layers and n_wavelengths, count data rows. */
    char line[65536];
    int n_data = 0;
    int n_layers_decl = 1;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            const char *p;
            if ((p = strstr(line, "n_layers"))) {
                sscanf(p, "n_layers = %d", &n_layers_decl);
                if (n_layers_decl < 1) n_layers_decl = 1;
            }
            if ((p = strstr(line, "T_ref"))) sscanf(p, "T_ref = %lf", &xsec->T_ref_K);
            if ((p = strstr(line, "P_ref"))) sscanf(p, "P_ref = %lf", &xsec->P_ref_atm);
            continue;
        }
        if (line[0] == '\n' || line[0] == '\0') continue;
        n_data++;
    }
    rewind(f);

    if (n_layers_decl == 1) {
        /* Legacy format: each row = (wl, sigma). */
        xsec->n_layers   = 1;
        xsec->n_wl       = n_data;
        xsec->wl_nm      = (double*)calloc(n_data, sizeof(double));
        xsec->sigma_cm2  = (double*)calloc(n_data, sizeof(double));
        int idx = 0;
        while (fgets(line, sizeof(line), f)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
            double wl, sigma;
            if (sscanf(line, "%lf %lf", &wl, &sigma) == 2) {
                xsec->wl_nm[idx]     = wl;
                xsec->sigma_cm2[idx] = sigma;
                idx++;
            }
        }
        fclose(f);
        if (!xsec_grid_valid_for_ocrt(xsec)) {
            fprintf(stderr,
                    "rt_xsec_load: %s does not cover %.0f-%.0f nm\n",
                    path, RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM);
            rt_xsec_free(xsec);
            return -3;
        }
        xsec->available = 1;
        return 0;
    }

    /* Layered format. */
    xsec->n_layers = n_layers_decl;
    /* First data row tells us n_wl; need to parse it. */
    int got_header_row = 0;
    int n_wl = 0;
    /* Pre-scan first data row to count fields.
     * Outer loop condition MUST include `*p != '\n'` — otherwise the
     * inner skips never advance past a stray '\n' inside the line and
     * the loop spins forever. */
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
        char *p = line;
        while (*p && *p != '\n') {
            while (*p && (*p == ' ' || *p == '\t')) p++;
            if (*p && *p != '\n') {
                n_wl++;
                while (*p && *p != ' ' && *p != '\t' && *p != '\n') p++;
            }
        }
        got_header_row = 1;
        break;
    }
    if (!got_header_row || n_wl < 2) {
        fclose(f);
        xsec->available = 0;
        return -2;
    }
    rewind(f);

    xsec->n_wl = n_wl;
    xsec->wl_nm = (double*)calloc(n_wl, sizeof(double));
    xsec->sigma_cm2 = (double*)calloc((size_t)n_layers_decl * n_wl, sizeof(double));

    /* Re-parse: first non-comment data row = wavelength header; subsequent K rows = sigma per layer. */
    int row_state = 0;  /* 0 = need wl header, 1+ = layer index +1 */
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
        char *p = line;
        if (row_state == 0) {
            /* Parse wavelengths */
            for (int i = 0; i < n_wl; i++) {
                char *endp;
                double v = strtod(p, &endp);
                if (p == endp) break;
                xsec->wl_nm[i] = v;
                p = endp;
            }
            row_state = 1;
        } else {
            int layer = row_state - 1;
            if (layer < n_layers_decl) {
                for (int i = 0; i < n_wl; i++) {
                    char *endp;
                    double v = strtod(p, &endp);
                    if (p == endp) break;
                    xsec->sigma_cm2[(size_t)layer * n_wl + i] = v;
                    p = endp;
                }
            }
            row_state++;
        }
    }
    fclose(f);
    if (!xsec_grid_valid_for_ocrt(xsec)) {
        fprintf(stderr,
                "rt_xsec_load: %s does not cover %.0f-%.0f nm\n",
                path, RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM);
        rt_xsec_free(xsec);
        return -3;
    }
    xsec->available = 1;
    return 0;
}

void rt_xsec_free(rt_xsec_t *xsec) {
    if (!xsec) return;
    free(xsec->wl_nm);
    free(xsec->sigma_cm2);
    memset(xsec, 0, sizeof(*xsec));
}

double rt_xsec_interp(const rt_xsec_t *xsec, double wavelength_nm) {
    /* Legacy single-layer or fallback to layer 0 of layered. */
    return rt_xsec_interp_layer(xsec, 0, wavelength_nm);
}

double rt_xsec_interp_layer(const rt_xsec_t *xsec, int layer_idx,
                             double wavelength_nm) {
    if (!xsec || !xsec->available || xsec->n_wl < 2) return 0.0;
    int n_wl = xsec->n_wl;
    int n_layers = xsec->n_layers;
    if (layer_idx < 0) layer_idx = 0;
    if (layer_idx >= n_layers) layer_idx = n_layers - 1;

    const double *sig_row = xsec->sigma_cm2 + (size_t)layer_idx * n_wl;

    if (wavelength_nm < xsec->wl_nm[0] ||
        wavelength_nm > xsec->wl_nm[n_wl - 1]) return NAN;
    if (wavelength_nm == xsec->wl_nm[0]) return sig_row[0];
    if (wavelength_nm == xsec->wl_nm[n_wl - 1]) return sig_row[n_wl - 1];

    int lo = 0, hi = n_wl - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (xsec->wl_nm[mid] <= wavelength_nm) lo = mid;
        else hi = mid;
    }
    double t = (wavelength_nm - xsec->wl_nm[lo]) /
               (xsec->wl_nm[hi] - xsec->wl_nm[lo]);
    return (1.0 - t) * sig_row[lo] + t * sig_row[hi];
}

/* Bi-linear interp (altitude, wavelength). Linear in z, linear in wavelength.
 * Never uses nearest-neighbor — full interpolation between bracketing AFGL levels.
 *
 * z_grid_km must be MONOTONIC ASCENDING (z_grid_km[0] = ground, [N-1] = TOA).
 */
double rt_xsec_interp_layer_z(const rt_xsec_t *xsec,
                               const double *z_grid_km,
                               double z_km, double wavelength_nm) {
    if (!xsec || !xsec->available || !z_grid_km) return 0.0;
    int n_layers = xsec->n_layers;
    if (n_layers < 1) return 0.0;
    if (n_layers == 1) {
        return rt_xsec_interp_layer(xsec, 0, wavelength_nm);
    }

    /* Bracket z_km in z_grid_km. */
    if (z_km <= z_grid_km[0]) {
        return rt_xsec_interp_layer(xsec, 0, wavelength_nm);
    }
    if (z_km >= z_grid_km[n_layers - 1]) {
        return rt_xsec_interp_layer(xsec, n_layers - 1, wavelength_nm);
    }

    int lo = 0, hi = n_layers - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (z_grid_km[mid] <= z_km) lo = mid;
        else hi = mid;
    }
    /* Bracket: z_grid_km[lo] <= z_km <= z_grid_km[hi], hi = lo+1 */
    double dz = z_grid_km[hi] - z_grid_km[lo];
    if (dz <= 0.0) {
        return rt_xsec_interp_layer(xsec, lo, wavelength_nm);
    }
    double t = (z_km - z_grid_km[lo]) / dz;
    double sig_lo = rt_xsec_interp_layer(xsec, lo, wavelength_nm);
    double sig_hi = rt_xsec_interp_layer(xsec, hi, wavelength_nm);
    return (1.0 - t) * sig_lo + t * sig_hi;
}

/* === AFGL utility functions === */

/* Cumulative gas column ABOVE altitude z  [molecules/cm^2].
 *
 *   N_above(z) = int_z^TOA n_gas(z') dz'
 *
 * Implemented as a trapezoid over the AFGL levels, with the partially-cut
 * bottom layer handled by linear interpolation of n_gas at z (n varies
 * roughly exponentially, so linear-in-layer is a small approximation on the
 * 1-km AFGL grid; consistent with column_default which uses the same rule).
 * AFGL z_km is ASCENDING (index 0 = ground). */
double rt_afgl_cumulative_column(const rt_afgl_atm_t *atm, int gas_idx,
                                  double z_km) {
    if (!atm || gas_idx < 0 || gas_idx >= RT_N_GAS) return 0.0;
    if (!atm->z_km || atm->n_levels < 2) return 0.0;

    /* AFGL z_km is in ascending order (0 km, 1, 2, ..., 100). */
    double cum = 0.0;
    /* Integrate from z_km up to the top, layer by layer. */
    for (int i = 0; i < atm->n_levels - 1; ++i) {
        double z1 = atm->z_km[i];
        double z2 = atm->z_km[i+1];
        if (z2 <= z_km) continue;  /* layer entirely below z_km */
        double z_lo_eff = (z1 < z_km) ? z_km : z1;
        double z_hi_eff = z2;
        if (z_hi_eff <= z_lo_eff) continue;

        /* Linear interpolation of n_gas(z) at endpoints */
        double f1 = (z_lo_eff - z1) / (z2 - z1);
        double f2 = (z_hi_eff - z1) / (z2 - z1);
        double n1 = atm->mr[gas_idx][i] * 1e-6 * atm->n_total_cm3[i];
        double n2 = atm->mr[gas_idx][i+1] * 1e-6 * atm->n_total_cm3[i+1];
        double n_lo = (1.0 - f1) * n1 + f1 * n2;
        double n_hi = (1.0 - f2) * n1 + f2 * n2;
        double dz_cm = (z_hi_eff - z_lo_eff) * 1.0e5;
        cum += 0.5 * (n_lo + n_hi) * dz_cm;
    }
    return cum;
}

double rt_afgl_layer_column(const rt_afgl_atm_t *atm, int gas_idx,
                             double z_lo_km, double z_hi_km) {
    if (z_hi_km <= z_lo_km) return 0.0;
    double cum_lo = rt_afgl_cumulative_column(atm, gas_idx, z_lo_km);
    double cum_hi = rt_afgl_cumulative_column(atm, gas_idx, z_hi_km);
    return cum_lo - cum_hi;  /* N_above(z_lo) - N_above(z_hi) = column in [z_lo, z_hi] */
}

static double afgl_interp_field(const rt_afgl_atm_t *atm,
                                 const double *field, double z_km) {
    if (!atm || !field || atm->n_levels < 2) return 0.0;
    if (z_km <= atm->z_km[0]) return field[0];
    if (z_km >= atm->z_km[atm->n_levels - 1]) return field[atm->n_levels - 1];
    int lo = 0, hi = atm->n_levels - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (atm->z_km[mid] <= z_km) lo = mid;
        else hi = mid;
    }
    double t = (z_km - atm->z_km[lo]) / (atm->z_km[hi] - atm->z_km[lo]);
    return (1.0 - t) * field[lo] + t * field[hi];
}

double rt_afgl_T_at(const rt_afgl_atm_t *atm, double z_km) {
    return afgl_interp_field(atm, atm ? atm->T_K : NULL, z_km);
}

double rt_afgl_P_at(const rt_afgl_atm_t *atm, double z_km) {
    return afgl_interp_field(atm, atm ? atm->P_mbar : NULL, z_km);
}

/* === Combined init === */

int rt_absorption_init(const char *afgl_dir, const char *xsec_dir,
                       rt_afgl_profile_t profile,
                       const double gas_column_override[RT_N_GAS],
                       rt_absorption_t *st) {
    if (!st) return -1;
    memset(st, 0, sizeof(*st));

    int rc = rt_afgl_load(afgl_dir, profile, &st->atm);
    if (rc != 0) return rc;

    /* Load each xsec file (silent fallback if missing) */
    int n_missing = 0;
    for (int g = 0; g < RT_N_GAS; ++g) {
        int r = rt_xsec_load(xsec_dir, g, &st->xsec[g]);
        if (r > 0) n_missing++;  /* missing file */
    }
    if (n_missing == RT_N_GAS) {
        fprintf(stderr,
            "rt_absorption_init: WARNING all %d xsec files missing from '%s'.\n"
            "  Absorption will be 0 (no correction applied).\n"
            "  Run scripts/generate_xsec_tables.py on a host with HITRAN online access.\n",
            RT_N_GAS, xsec_dir);
    } else if (n_missing > 0) {
        fprintf(stderr,
            "rt_absorption_init: %d/%d xsec files missing (some gases will contribute zero).\n",
            n_missing, RT_N_GAS);
    }

    /* Effective column = default OR user override */
    for (int g = 0; g < RT_N_GAS; ++g) {
        if (gas_column_override && gas_column_override[g] >= 0.0) {
            st->column_eff[g] = gas_column_override[g];
            st->column_overridden[g] = 1;
        } else {
            st->column_eff[g] = st->atm.column_default[g];
            st->column_overridden[g] = 0;
        }
    }

    st->initialized = 1;
    return 0;
}

void rt_absorption_free(rt_absorption_t *st) {
    if (!st) return;
    rt_afgl_free(&st->atm);
    for (int g = 0; g < RT_N_GAS; ++g) rt_xsec_free(&st->xsec[g]);
    memset(st, 0, sizeof(*st));
}

/* === τ_abs computations === */

/* Per-gas vertical absorption optical depth  [dimensionless].
 *
 *   tau_gas(lambda) = sum_{k=0}^{n_lev-2} sigma(z_mid_k, lambda) * N_k
 *
 *   sigma(z_mid, lambda) : bilinear (z, lambda) interpolation of the xsec
 *       LUT at the layer midpoint — captures the T/p dependence of the
 *       cross section through the LUT's altitude layering (e.g. O3
 *       Huggins-band T-dependence).
 *   N_k : layer column from the AFGL profile (trapezoid, see
 *       rt_afgl_layer_column), optionally rescaled by
 *       column_eff/column_default when the user overrides the total.
 *
 * v1.03 anomaly fix (HISTORY — keep for regression archaeology): earlier
 * versions used surface-sigma x N_total here while the SOS internal gas
 * weighting used layered integration, making direct-beam attenuation
 * (exp(-tau/mu_sun)) inconsistent with diffuse Ed by 1-2 %pp.  Both paths
 * now share this layered tau. */
double rt_absorption_tau_per_gas(const rt_absorption_t *st,
                                  int gas_idx, double wavelength_nm) {
    if (!st || !st->initialized || gas_idx < 0 || gas_idx >= RT_N_GAS) return 0.0;
    const rt_xsec_t *xs = &st->xsec[gas_idx];
    if (!xs->available) return 0.0;

    int n_lev = st->atm.n_levels;
    if (n_lev < 2 || !st->atm.z_km) return 0.0;

    double tau = 0.0;
    double N_def_total = st->atm.column_default[gas_idx];
    int override_on = st->column_overridden[gas_idx] && N_def_total > 0.0;
    double scale = override_on ? (st->column_eff[gas_idx] / N_def_total) : 1.0;

    for (int k = 0; k < n_lev - 1; k++) {
        /* AFGL z grid is ASCENDING (st->atm.z_km[0]=ground, last=TOA) */
        double z_lo = st->atm.z_km[k];
        double z_hi = st->atm.z_km[k+1];
        double z_mid = 0.5 * (z_lo + z_hi);

        /* layered σ at z_mid */
        double sigma = rt_xsec_interp_layer_z(xs, st->atm.z_km,
                                                z_mid, wavelength_nm);

        /* layer column from AFGL profile */
        double N_layer = rt_afgl_layer_column(&st->atm, gas_idx, z_lo, z_hi);

        /* user-overridden column → proportional scaling */
        if (override_on) N_layer *= scale;

        tau += sigma * N_layer;
    }
    return tau;
}

double rt_absorption_tau_total(const rt_absorption_t *st,
                                double wavelength_nm) {
    if (!st || !st->initialized) return 0.0;
    double tau = 0.0;
    for (int g = 0; g < RT_N_GAS; ++g) {
        tau += rt_absorption_tau_per_gas(st, g, wavelength_nm);
    }
    return tau;
}

/* Two-path gas transmittance correction  [0..1].
 *
 *   T2(tau; mu_s, mu_v) = exp( -tau * (1/mu_s + 1/mu_v) )
 *
 * Plane-parallel airmass approximation: the DIRECT sun->surface->sensor
 * ray pair, each slant path scaling the vertical tau by 1/mu.  This is the
 * standard decoupled gas-correction of ocean-color AC (absorption applied
 * multiplicatively outside the scattering solve); exact only if the gas
 * is optically decoupled from scattering (weak absorption or gas above
 * the scattering layer, e.g. stratospheric O3).  Horizon guard: mu <= 1e-6
 * returns 0 (fully attenuated) rather than exp(+inf) garbage. */
double rt_absorption_correction_two_pass(double tau_abs,
                                          double mu_sun, double mu_v) {
    if (tau_abs <= 0.0) return 1.0;
    if (mu_sun <= 1e-6 || mu_v <= 1e-6) return 0.0;  /* horizon edge */
    return exp(-tau_abs * (1.0 / mu_sun + 1.0 / mu_v));
}

/* One-path (sun->surface only) variant: T1 = exp(-tau/mu_s).  Used where
 * the upward path is handled elsewhere (e.g. in-water quantities). */
double rt_absorption_correction_one_pass(double tau_abs, double mu_sun) {
    if (tau_abs <= 0.0) return 1.0;
    if (mu_sun <= 1e-6) return 0.0;
    return exp(-tau_abs / mu_sun);
}
