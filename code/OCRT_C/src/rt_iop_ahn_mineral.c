#define _USE_MATH_DEFINES
#include "rt_iop_ahn_mineral.h"
#include "rt_spectral_contract.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Supplied AHN files carry 757 bulk spectral rows. */
#define AHN_MAX_STAR 1024
#define AHN_MAX_PHW  1024
#define AHN_PATH_MAX 1024

typedef struct {
    int n_star;
    double wl_star[AHN_MAX_STAR];
    double astar[AHN_MAX_STAR];
    double bstar[AHN_MAX_STAR];

    int n_phw;
    double wl_phw[AHN_MAX_PHW];
    double bbratio[AHN_MAX_PHW];

    char mie_path[AHN_PATH_MAX];
    int loaded;
} ahn_species_data_t;

static ahn_species_data_t g_species[AHN_SPECIES_COUNT];
static int g_ready = 0;

static const char *const k_star_stem[AHN_SPECIES_COUNT] = {
    "redclay", "brownearth", "yellowclay", "calcareoussand"
};
static const char *const k_mie_stem[AHN_SPECIES_COUNT] = {
    "Red_clay", "Brown_earth", "Yellow_clay", "Calcareous_sand"
};
static const char *const k_species_name[AHN_SPECIES_COUNT] = {
    "red_clay", "brown_earth", "yellow_clay", "calcareous_sand"
};

static int species_valid(ahn_species_t species)
{
    return species >= AHN_RED_CLAY && species < AHN_SPECIES_COUNT;
}

static double interpolate(const double *x, const double *y, int n, double xq)
{
    if (!x || !y || n <= 0 || !isfinite(xq)) return 0.0;
    if (xq <= x[0]) return y[0];
    if (xq >= x[n - 1]) return y[n - 1];

    int lo = 0;
    int hi = n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (x[mid] <= xq) lo = mid;
        else hi = mid;
    }
    const double dx = x[hi] - x[lo];
    if (!(dx > 0.0)) return y[lo];
    const double t = (xq - x[lo]) / dx;
    return y[lo] + t * (y[hi] - y[lo]);
}

static int table_valid(const double *wl, const double *v, int n,
                       int strictly_positive)
{
    if (!wl || !v || n < 2) return 0;
    for (int i = 0; i < n; ++i) {
        if (!isfinite(wl[i]) || !isfinite(v[i]) || !(wl[i] > 0.0)) return 0;
        if (i > 0 && !(wl[i] > wl[i - 1])) return 0;
        if (strictly_positive) {
            if (!(v[i] > 0.0)) return 0;
        } else if (v[i] < 0.0) {
            return 0;
        }
    }
    return 1;
}

/* Optional two-column table. -2 means absent and permits .mie fallback. */
static int load_star_file(const char *path, double *wl, double *value, int cap)
{
    FILE *fp = fopen(path, "r");
    if (!fp) return (errno == ENOENT) ? -2 : -1;

    char line[1024];
    int n = 0;
    while (fgets(line, sizeof line, fp)) {
        double w = 0.0;
        double v = 0.0;
        if (sscanf(line, "%lf %lf", &w, &v) != 2) continue;
        if (w < 0.0) break;
        if (n >= cap) {
            fclose(fp);
            fprintf(stderr, "[ahn] too many rows in %s (capacity %d)\n", path, cap);
            return -3;
        }
        wl[n] = w;
        value[n] = v;
        ++n;
    }
    fclose(fp);
    return table_valid(wl, value, n, 0) ? n : -4;
}

/* Read the .mie bulk spectral block. Extinct_Co and Scatter_Co in the
 * delivered files are the mass-specific coefficients [m^2 g^-1]. */
static int load_mie_star(const char *path,
                         double *wl_nm, double *astar, double *bstar, int cap)
{
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[ahn] cannot open %s\n", path);
        return -1;
    }

    char line[8192];
    if (!fgets(line, sizeof line, fp) || !fgets(line, sizeof line, fp)) {
        fclose(fp);
        return -2;
    }

    int n = 0;
    while (fgets(line, sizeof line, fp)) {
        if (strstr(line, "Phase Function") != NULL) break;
        double wl_um, nor_ext, nor_sca, ssa, asym, ext, sca;
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
                   &wl_um, &nor_ext, &nor_sca, &ssa, &asym, &ext, &sca) != 7)
            continue;
        (void)nor_ext;
        (void)nor_sca;
        (void)ssa;
        (void)asym;
        if (n >= cap) {
            fclose(fp);
            return -3;
        }
        if (!(wl_um > 0.0) || !(ext >= 0.0) || !(sca >= 0.0) || sca > ext + 1e-10) {
            fclose(fp);
            return -4;
        }
        wl_nm[n] = 1000.0 * wl_um;
        astar[n] = fmax(0.0, ext - sca);
        bstar[n] = sca;
        ++n;
    }
    fclose(fp);

    if (!table_valid(wl_nm, astar, n, 0) || !table_valid(wl_nm, bstar, n, 1))
        return -5;
    return n;
}

/* Exact integral of a theta-linear phase segment times sin(theta).
 * This is the same integration contract used by rt_value_phase.c, kept local
 * so the standalone AHN smoke test does not acquire a solver dependency. */
static double ahn_segment_int_p_sin(double t0, double t1,
                                    double p0, double p1)
{
    const double h = t1 - t0;
    if (!(h > 0.0)) return 0.0;
    const double slope = (p1 - p0) / h;
    const double f0 = -p0 * cos(t0) + slope * sin(t0);
    const double f1 = -p1 * cos(t1) + slope * sin(t1);
    return f1 - f0;
}

/* Integrate the matching P11 over the backward and full hemispheres. */
static int load_mie_bbratio(const char *path,
                            double *wl_phw, double *bbratio, int cap)
{
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "[ahn] cannot open %s\n", path);
        return -1;
    }

    char line[16384];
    int found = 0;
    while (fgets(line, sizeof line, fp)) {
        if (strstr(line, "Phase Function (P11)") != NULL) {
            found = 1;
            break;
        }
    }
    if (!found || !fgets(line, sizeof line, fp)) {
        fclose(fp);
        return -2;
    }

    double wl_um[AHN_MAX_PHW];
    int nw = 0;
    char *token = strtok(line, " \t\r\n"); /* TETA */
    if (!token) {
        fclose(fp);
        return -3;
    }
    while ((token = strtok(NULL, " \t\r\n")) != NULL) {
        if (nw >= cap || nw >= AHN_MAX_PHW) {
            fclose(fp);
            return -4;
        }
        wl_um[nw++] = atof(token);
    }
    if (nw < 2) {
        fclose(fp);
        return -5;
    }

    double numerator[AHN_MAX_PHW] = {0};
    double denominator[AHN_MAX_PHW] = {0};
    double previous_theta = 0.0;
    double previous_p[AHN_MAX_PHW] = {0};
    int have_previous = 0;
    int n_angles = 0;

    while (fgets(line, sizeof line, fp)) {
        if (strstr(line, "Phase Function") != NULL) break;

        double row[1 + AHN_MAX_PHW];
        int count = 0;
        token = strtok(line, " \t\r\n");
        while (token && count < 1 + nw) {
            row[count++] = atof(token);
            token = strtok(NULL, " \t\r\n");
        }
        if (count != 1 + nw) continue;

        const double theta = row[0] * M_PI / 180.0;
        if (have_previous) {
            const int ascending = theta > previous_theta;
            const double t0 = ascending ? previous_theta : theta;
            const double t1 = ascending ? theta : previous_theta;
            const double half_pi = 0.5 * M_PI;
            for (int k = 0; k < nw; ++k) {
                const double p0 = ascending ? previous_p[k] : row[1 + k];
                const double p1 = ascending ? row[1 + k] : previous_p[k];
                const double contribution =
                    ahn_segment_int_p_sin(t0, t1, p0, p1);
                denominator[k] += contribution;
                if (t0 >= half_pi) {
                    numerator[k] += contribution;
                } else if (t1 > half_pi) {
                    const double weight = (half_pi - t0) / (t1 - t0);
                    const double pm = p0 + weight * (p1 - p0);
                    numerator[k] +=
                        ahn_segment_int_p_sin(half_pi, t1, pm, p1);
                }
            }
        }
        previous_theta = theta;
        for (int k = 0; k < nw; ++k) previous_p[k] = row[1 + k];
        have_previous = 1;
        ++n_angles;
    }
    fclose(fp);
    if (n_angles < 3) return -6;

    for (int k = 0; k < nw; ++k) {
        wl_phw[k] = 1000.0 * wl_um[k];
        const double den = fabs(denominator[k]);
        const double ratio = den > 0.0 ? fabs(numerator[k]) / den : -1.0;
        if (!(ratio > 0.0 && ratio < 0.5) || !isfinite(ratio)) return -7;
        bbratio[k] = ratio;
    }

    for (int i = 1; i < nw; ++i) {
        int j = i;
        while (j > 0 && wl_phw[j] < wl_phw[j - 1]) {
            double tmp = wl_phw[j];
            wl_phw[j] = wl_phw[j - 1];
            wl_phw[j - 1] = tmp;
            tmp = bbratio[j];
            bbratio[j] = bbratio[j - 1];
            bbratio[j - 1] = tmp;
            --j;
        }
    }
    return table_valid(wl_phw, bbratio, nw, 1) ? nw : -8;
}

void rt_iop_ahn_mineral_free(void)
{
    memset(g_species, 0, sizeof g_species);
    g_ready = 0;
}

int rt_iop_ahn_mineral_init(const char *data_dir)
{
    if (!data_dir || !data_dir[0]) data_dir = ".";
    rt_iop_ahn_mineral_free();

    char path_a[AHN_PATH_MAX];
    char path_b[AHN_PATH_MAX];
    for (int s = 0; s < AHN_SPECIES_COUNT; ++s) {
        ahn_species_data_t *data = &g_species[s];
        int n = snprintf(data->mie_path, sizeof data->mie_path,
                         "%s/%s_AHN.mie", data_dir, k_mie_stem[s]);
        if (n < 0 || (size_t)n >= sizeof data->mie_path) {
            rt_iop_ahn_mineral_free();
            return -1;
        }
        n = snprintf(path_a, sizeof path_a, "%s/astarmin_%s.txt",
                     data_dir, k_star_stem[s]);
        if (n < 0 || (size_t)n >= sizeof path_a) {
            rt_iop_ahn_mineral_free();
            return -1;
        }
        n = snprintf(path_b, sizeof path_b, "%s/bstarmin_%s.txt",
                     data_dir, k_star_stem[s]);
        if (n < 0 || (size_t)n >= sizeof path_b) {
            rt_iop_ahn_mineral_free();
            return -1;
        }

        double wa[AHN_MAX_STAR], av[AHN_MAX_STAR];
        double wb[AHN_MAX_STAR], bv[AHN_MAX_STAR];
        const int na = load_star_file(path_a, wa, av, AHN_MAX_STAR);
        const int nb = load_star_file(path_b, wb, bv, AHN_MAX_STAR);

        if (na > 0 && nb > 0) {
            data->n_star = na;
            for (int i = 0; i < na; ++i) {
                data->wl_star[i] = wa[i];
                data->astar[i] = av[i];
                data->bstar[i] = interpolate(wb, bv, nb, wa[i]);
            }
        } else if ((na == -2 || na > 0) && (nb == -2 || nb > 0)) {
            const int nm = load_mie_star(data->mie_path, data->wl_star,
                                         data->astar, data->bstar, AHN_MAX_STAR);
            if (nm <= 0) {
                rt_iop_ahn_mineral_free();
                return -2;
            }
            data->n_star = nm;
        } else {
            fprintf(stderr, "[ahn] malformed a*/b* data for %s\n", k_species_name[s]);
            rt_iop_ahn_mineral_free();
            return -2;
        }

        const int nphase = load_mie_bbratio(data->mie_path, data->wl_phw,
                                            data->bbratio, AHN_MAX_PHW);
        if (nphase <= 0) {
            rt_iop_ahn_mineral_free();
            return -3;
        }
        data->n_phw = nphase;
        if (!rt_spectral_grid_covers(data->wl_star[0],
                                     data->wl_star[data->n_star - 1]) ||
            !rt_spectral_grid_covers(data->wl_phw[0],
                                     data->wl_phw[data->n_phw - 1])) {
            fprintf(stderr,
                    "[ahn] %s data do not cover %.0f-%.0f nm "
                    "(scalar %.10g-%.10g; phase %.10g-%.10g)\n",
                    k_species_name[s], RT_SPECTRAL_MIN_NM, RT_SPECTRAL_MAX_NM,
                    data->wl_star[0], data->wl_star[data->n_star - 1],
                    data->wl_phw[0], data->wl_phw[data->n_phw - 1]);
            rt_iop_ahn_mineral_free();
            return -4;
        }
        data->loaded = 1;
    }

    g_ready = 1;
    return 0;
}

int rt_iop_ahn_mineral_ready(void)
{
    return g_ready;
}

ahn_species_t rt_iop_ahn_species_from_name(const char *name)
{
    if (!name || !name[0]) return AHN_RED_CLAY;

    char normalized[96];
    int n = 0;
    for (const unsigned char *p = (const unsigned char *)name;
         *p && n < (int)sizeof normalized - 1; ++p) {
        if (isalnum(*p)) normalized[n++] = (char)tolower(*p);
    }
    normalized[n] = '\0';

    if (strstr(normalized, "brownearth") || !strcmp(normalized, "brown"))
        return AHN_BROWN_EARTH;
    if (strstr(normalized, "yellowclay") || !strcmp(normalized, "yellow"))
        return AHN_YELLOW_CLAY;
    if (strstr(normalized, "calcareoussand") || strstr(normalized, "calcareous") ||
        !strcmp(normalized, "sand"))
        return AHN_CALCAREOUS_SAND;
    return AHN_RED_CLAY;
}

const char *rt_iop_ahn_species_name(ahn_species_t species)
{
    if (!species_valid(species)) species = AHN_RED_CLAY;
    return k_species_name[species];
}

const char *rt_iop_ahn_mineral_phase_path(ahn_species_t species)
{
    if (!g_ready) return NULL;
    if (!species_valid(species)) species = AHN_RED_CLAY;
    return g_species[species].loaded ? g_species[species].mie_path : NULL;
}

int rt_iop_ahn_mineral_eval(double lambda_nm, double tsm_g_m3,
                            ahn_species_t species, rt_iop_t *out,
                            const char **mie_phase_path)
{
    if (!out || !rt_spectral_wavelength_supported(lambda_nm) ||
        !isfinite(tsm_g_m3))
        return -1;
    if (!g_ready) return -2;
    if (!species_valid(species)) species = AHN_RED_CLAY;

    const ahn_species_data_t *data = &g_species[species];
    if (!data->loaded) return -3;
    if (mie_phase_path) *mie_phase_path = data->mie_path;

    out->a = 0.0;
    out->b = 0.0;
    out->bb = 0.0;
    if (tsm_g_m3 <= 0.0) return 0;

    const double astar = interpolate(data->wl_star, data->astar,
                                     data->n_star, lambda_nm);
    const double bstar = interpolate(data->wl_star, data->bstar,
                                     data->n_star, lambda_nm);
    const double bb_over_b = interpolate(data->wl_phw, data->bbratio,
                                         data->n_phw, lambda_nm);
    if (!(astar >= 0.0) || !(bstar > 0.0) ||
        !(bb_over_b > 0.0 && bb_over_b < 0.5)) return -4;

    out->a = astar * tsm_g_m3;
    out->b = bstar * tsm_g_m3;
    out->bb = out->b * bb_over_b;
    return 0;
}
