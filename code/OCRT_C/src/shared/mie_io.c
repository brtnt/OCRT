/* shared/mie_io.c
 *
 * Implementation of mie_data_t, read_mie_file, and build_aerosol_interpolators.
 * Combined from vrt_solver.c and ocean_solver.c:
 *   - read_mie_file: vrt version (handles OPAC layout: " 83\n column header\n
 *     spectral rows + Phase Function blocks")
 *   - build_aerosol_interpolators: vrt's angle-sort logic + ocean's
 *     diagnostic fields (ssa, g_asym, bb_b_ratio)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "mie_io.h"
#include "numerics.h"
#include "io_utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef LINE_BUFFER_SIZE
#define LINE_BUFFER_SIZE 131072
#endif

/* ─── mie_data_t ──────────────────────────────────────── */

void mie_data_init(mie_data_t* m, int n_ang, int n_wl, int n_phase_wl) {
    m->n_ang = n_ang;
    m->n_wl  = n_wl;
    m->n_phase_wl = n_phase_wl;
    m->wavelengths       = (double*)calloc(n_wl, sizeof(double));
    m->phase_wavelengths = (double*)calloc(n_phase_wl, sizeof(double));
    m->angles      = (double*)calloc(n_ang, sizeof(double));
    m->spectral    = (double*)calloc((size_t)n_wl * 6, sizeof(double));
    m->P11 = (double*)calloc((size_t)n_ang * n_phase_wl, sizeof(double));
    m->P12 = (double*)calloc((size_t)n_ang * n_phase_wl, sizeof(double));
    m->P33 = (double*)calloc((size_t)n_ang * n_phase_wl, sizeof(double));
}

void mie_data_free(mie_data_t* m) {
    free(m->wavelengths); free(m->phase_wavelengths); free(m->angles); free(m->spectral);
    free(m->P11); free(m->P12); free(m->P33);
    memset(m, 0, sizeof(*m));
}

/* ─── worker-private complete-model cache ─────────────── */

#define MIE_MODEL_CACHE_CAP 8

typedef struct {
    int valid;
    char path[1024];
    long long file_size;
    long long file_mtime;
    unsigned long long stamp;
    size_t resident_bytes;
    mie_data_t mie;
} mie_model_cache_entry_t;

static _Thread_local mie_model_cache_entry_t
    g_mie_model_cache[MIE_MODEL_CACHE_CAP];
static _Thread_local unsigned long long g_mie_model_tick;
static _Thread_local unsigned long long g_mie_model_loads;
static _Thread_local unsigned long long g_mie_model_hits;
static _Thread_local size_t g_mie_model_resident;

static size_t mie_data_resident_bytes_impl(const mie_data_t* m)
{
    if (!m) return 0;
    return (size_t)m->n_wl * sizeof(double) +
           (size_t)m->n_phase_wl * sizeof(double) +
           (size_t)m->n_ang * sizeof(double) +
           (size_t)m->n_wl * 6u * sizeof(double) +
           (size_t)m->n_ang * (size_t)m->n_phase_wl * 3u * sizeof(double);
}

static void mie_model_cache_entry_clear(mie_model_cache_entry_t* e)
{
    if (!e) return;
    if (e->valid) {
        if (g_mie_model_resident >= e->resident_bytes)
            g_mie_model_resident -= e->resident_bytes;
        else
            g_mie_model_resident = 0;
        mie_data_free(&e->mie);
    }
    memset(e, 0, sizeof(*e));
}

void mie_model_cache_clear_worker(void)
{
    for (int i = 0; i < MIE_MODEL_CACHE_CAP; ++i)
        mie_model_cache_entry_clear(&g_mie_model_cache[i]);
    g_mie_model_tick = 0;
    g_mie_model_loads = 0;
    g_mie_model_hits = 0;
    g_mie_model_resident = 0;
}

size_t mie_model_cache_resident_bytes(void) { return g_mie_model_resident; }
unsigned long long mie_model_cache_load_count(void) { return g_mie_model_loads; }
unsigned long long mie_model_cache_hit_count(void) { return g_mie_model_hits; }

int mie_model_cache_get(const char* path, const mie_data_t** out, int* cache_hit)
{
    if (out) *out = NULL;
    if (cache_hit) *cache_hit = 0;
    if (!path || !path[0] || !out) return -1;

    const char* disable = getenv("OCRT_DISABLE_MIE_MODEL_CACHE");
    if (disable && disable[0] && strcmp(disable, "0") != 0) return 1;

    struct stat st;
    if (stat(path, &st) != 0) return -2;
    const long long file_size = (long long)st.st_size;
    const long long file_mtime = (long long)st.st_mtime;

    for (int i = 0; i < MIE_MODEL_CACHE_CAP; ++i) {
        mie_model_cache_entry_t* e = &g_mie_model_cache[i];
        if (!e->valid || e->file_size != file_size ||
            e->file_mtime != file_mtime || strcmp(e->path, path) != 0) continue;
        e->stamp = ++g_mie_model_tick;
        ++g_mie_model_hits;
        *out = &e->mie;
        if (cache_hit) *cache_hit = 1;
        const char* trace = getenv("OCRT_MIE_MODEL_CACHE_TRACE");
        if (trace && trace[0] && strcmp(trace, "0") != 0)
            fprintf(stderr,
                    "[MIE-MODEL-CACHE] hit path=%s nang=%d nwl=%d nphase=%d "
                    "resident=%zu loads=%llu hits=%llu\n",
                    path, e->mie.n_ang, e->mie.n_wl, e->mie.n_phase_wl,
                    g_mie_model_resident, g_mie_model_loads, g_mie_model_hits);
        return 0;
    }

    int victim = 0;
    for (int i = 0; i < MIE_MODEL_CACHE_CAP; ++i) {
        if (!g_mie_model_cache[i].valid) { victim = i; break; }
        if (g_mie_model_cache[i].stamp < g_mie_model_cache[victim].stamp)
            victim = i;
    }
    mie_model_cache_entry_t* e = &g_mie_model_cache[victim];
    mie_model_cache_entry_clear(e);
    if (snprintf(e->path, sizeof(e->path), "%s", path) >= (int)sizeof(e->path))
        return -3;
    if (read_mie_file(path, &e->mie) != 0) {
        memset(e, 0, sizeof(*e));
        return -4;
    }
    e->valid = 1;
    e->file_size = file_size;
    e->file_mtime = file_mtime;
    e->stamp = ++g_mie_model_tick;
    e->resident_bytes = mie_data_resident_bytes_impl(&e->mie);
    g_mie_model_resident += e->resident_bytes;
    ++g_mie_model_loads;
    *out = &e->mie;
    const char* trace = getenv("OCRT_MIE_MODEL_CACHE_TRACE");
    if (trace && trace[0] && strcmp(trace, "0") != 0)
        fprintf(stderr,
                "[MIE-MODEL-CACHE] load path=%s nang=%d nwl=%d nphase=%d "
                "entry_bytes=%zu resident=%zu loads=%llu hits=%llu\n",
                path, e->mie.n_ang, e->mie.n_wl, e->mie.n_phase_wl,
                e->resident_bytes, g_mie_model_resident,
                g_mie_model_loads, g_mie_model_hits);
    return 0;
}

/* ─── read_mie_file (vrt format, OPAC layout) ─────────── */

int read_mie_file(const char* path, mie_data_t* out) {
    FILE* f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Cannot open .mie file: %s\n", path);
        return -1;
    }
    char buf[LINE_BUFFER_SIZE];

    /* Line 1: n_ang */
    int n_ang = 0;
    while (read_line(f, buf, sizeof(buf))) {
        if (strlen(buf) > 0 && line_starts_numeric(buf)) {
            n_ang = atoi(buf);
            break;
        }
    }
    if (n_ang <= 0) { fclose(f); return -1; }

    /* Read all remaining lines */
    typedef struct { char* text; } line_t;
    line_t* lines = NULL; int n_lines = 0, cap = 0;
    while (read_line(f, buf, sizeof(buf))) {
        if (n_lines >= cap) {
            cap = cap ? cap * 2 : 512;
            lines = (line_t*)realloc(lines, cap * sizeof(line_t));
        }
        lines[n_lines].text = my_strdup(buf);
        n_lines++;
    }
    fclose(f);

    /* Find spectral header (first non-numeric line) */
    int i = 0;
    while (i < n_lines && !lines[i].text[0]) i++;
    while (i < n_lines && line_starts_numeric(lines[i].text)) i++;  /* skip if starts numeric */
    while (i < n_lines && !line_starts_numeric(lines[i].text)) i++;  /* skip header lines */

    /* Now lines[i] is first spectral data row */
    int spec_start = i;
    while (i < n_lines && line_starts_numeric(lines[i].text)) i++;
    int n_wl = i - spec_start;
    if (n_wl <= 0) {
        for (int j = 0; j < n_lines; ++j) free(lines[j].text);
        free(lines);
        return -1;
    }

    /* Detect the phase-matrix wavelength grid.  OPAC-style files place a
     * "TETA <wl0> <wl1> ..." header before the phase blocks whose wavelength
     * grid is independent of (and usually coarser than) the spectral/IOP grid.
     * Legacy single-grid files have no TETA header and phase rows carry n_wl
     * values; for those we fall back to n_phase_wl = n_wl and reuse the
     * spectral grid for phase wavelength interpolation. */
    int    n_phase_wl = n_wl;
    int    have_teta  = 0;
    double phase_wl_tmp[1024];
    for (int j = i; j < n_lines; ++j) {
        const char* s = lines[j].text;
        while (*s == ' ' || *s == '\t') s++;
        if (strncmp(s, "TETA", 4) == 0) {
            int nt = parse_floats(s + 4, phase_wl_tmp, 1024);
            if (nt > 0) { n_phase_wl = nt; have_teta = 1; }
            break;
        }
        if (line_starts_numeric(s)) break;   /* first phase data row → legacy (no TETA) */
    }

    mie_data_init(out, n_ang, n_wl, n_phase_wl);
    for (int l = 0; l < n_wl; ++l) {
        double row[8];
        int nparsed = parse_floats(lines[spec_start + l].text, row, 7);
        if (nparsed < 7) {
            fprintf(stderr, "Spectral row %d has %d tokens\n", l, nparsed);
            for (int j = 0; j < n_lines; ++j) free(lines[j].text);
            free(lines);
            mie_data_free(out);
            return -1;
        }
        out->wavelengths[l] = row[0];
        for (int c = 0; c < 6; ++c) SPECTRAL_AT(out, l, c) = row[c + 1];
    }

    /* Phase wavelength grid: from the TETA header (separate grid) or a copy of
     * the spectral grid (legacy single-grid).  When single-grid, this leaves
     * every downstream phase λ-interpolation numerically identical to before. */
    for (int l = 0; l < n_phase_wl; ++l)
        out->phase_wavelengths[l] = have_teta ? phase_wl_tmp[l] : out->wavelengths[l];

    /* Parse phase blocks: each block has n_ang rows of (angle + n_phase_wl
     * values).  Use a dynamic row so files with arbitrary band counts (and a
     * phase grid distinct from the spectral grid) are accepted safely. */
    double* blocks_out[3] = { out->P11, out->P12, out->P33 };
    double* phase_row = (double*)malloc((size_t)(n_phase_wl + 1) * sizeof(double));
    if (!phase_row) {
        for (int j = 0; j < n_lines; ++j) free(lines[j].text);
        free(lines);
        mie_data_free(out);
        return -1;
    }
    int blocks_found = 0;
    while (i < n_lines && blocks_found < 3) {
        /* Skip until next (n_phase_wl+1)-column numeric row (skips the
         * "Phase Function (...)" and "TETA ..." header lines). */
        while (i < n_lines) {
            int nt = parse_floats(lines[i].text, phase_row, n_phase_wl + 1);
            if (nt == n_phase_wl + 1) break;
            i++;
        }
        if (i >= n_lines) break;
        int rows_in_block = 0;
        for (int k = 0; k < n_ang; ++k) {
            if (i + k >= n_lines) { break; }
            int nt = parse_floats(lines[i + k].text, phase_row, n_phase_wl + 1);
            if (nt < n_phase_wl + 1) break;
            if (blocks_found == 0) out->angles[k] = phase_row[0];
            for (int l = 0; l < n_phase_wl; ++l) {
                BLOCK_AT(blocks_out[blocks_found], n_phase_wl, k, l) = phase_row[l + 1];
            }
            rows_in_block++;
        }
        if (rows_in_block != n_ang) break;
        i += n_ang;
        blocks_found++;
    }

    free(phase_row);
    for (int j = 0; j < n_lines; ++j) free(lines[j].text);
    free(lines);
    if (blocks_found < 3) {
        fprintf(stderr, "Expected 3 phase blocks (P11, P12, P33), found %d\n", blocks_found);
        mie_data_free(out);
        return -1;
    }
    return 0;
}


int mie_phase_nodes_at_wavelength_pchip(const mie_data_t* mie,
                                        double wavelength_um,
                                        double* P11,
                                        double* P12,
                                        double* P33)
{
    if (!mie || !P11 || !P12 || !P33 || mie->n_ang < 2 ||
        mie->n_phase_wl < 1 || !mie_data_covers_phase_wavelength(mie, wavelength_um))
        return -1;
    const int nw = mie->n_phase_wl;
    double *row = (double*)malloc((size_t)nw * sizeof(double));
    if (!row) return -2;
    const double *blocks[3] = { mie->P11, mie->P12, mie->P33 };
    double *outs[3] = { P11, P12, P33 };
    for (int b = 0; b < 3; ++b) {
        for (int ia = 0; ia < mie->n_ang; ++ia) {
            for (int iw = 0; iw < nw; ++iw)
                row[iw] = BLOCK_AT(blocks[b], nw, ia, iw);
            outs[b][ia] = pchip_interp(mie->phase_wavelengths, row, nw,
                                       wavelength_um);
            if (!isfinite(outs[b][ia])) { free(row); return -3; }
        }
    }
    free(row);
    return 0;
}

int mie_phase_nodes_at_wavelength_linear(const mie_data_t* mie,
                                         double wavelength_um,
                                         double* P11,
                                         double* P12,
                                         double* P33)
{
    if (!mie || !P11 || !P12 || !P33 || mie->n_ang < 2 ||
        mie->n_phase_wl < 1 || !mie_data_covers_phase_wavelength(mie, wavelength_um))
        return -1;

    const int nw = mie->n_phase_wl;
    int lo = 0;
    double w = 0.0;
    if (nw > 1) {
        if (wavelength_um <= mie->phase_wavelengths[0]) {
            lo = 0; w = 0.0;
        } else if (wavelength_um >= mie->phase_wavelengths[nw - 1]) {
            lo = nw - 2; w = 1.0;
        } else {
            int hi = nw - 1;
            while (hi - lo > 1) {
                const int mid = (lo + hi) >> 1;
                if (mie->phase_wavelengths[mid] <= wavelength_um) lo = mid;
                else hi = mid;
            }
            const double h = mie->phase_wavelengths[lo + 1] -
                             mie->phase_wavelengths[lo];
            if (!(h > 0.0)) return -2;
            w = (wavelength_um - mie->phase_wavelengths[lo]) / h;
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;
        }
    }

    const double *blocks[3] = { mie->P11, mie->P12, mie->P33 };
    double *outs[3] = { P11, P12, P33 };
    for (int ia = 0; ia < mie->n_ang; ++ia) {
        for (int b = 0; b < 3; ++b) {
            const double a = BLOCK_AT(blocks[b], nw, ia, lo);
            const double v = (nw == 1) ? a :
                a + w * (BLOCK_AT(blocks[b], nw, ia, lo + 1) - a);
            if (!isfinite(v)) return -3;
            outs[b][ia] = v;
        }
    }
    return 0;
}

int mie_bulk_diagnostics_at_wavelength(const mie_data_t* mie,
                                       double wavelength_um,
                                       double* ssa,
                                       double* g_asym)
{
    if (!mie || mie->n_wl < 1 ||
        !mie_data_covers_bulk_wavelength(mie, wavelength_um)) return -1;
    double *ssa_arr = NULL, *g_arr = NULL;
    if (ssa) ssa_arr = (double*)malloc((size_t)mie->n_wl * sizeof(double));
    if (g_asym) g_arr = (double*)malloc((size_t)mie->n_wl * sizeof(double));
    if ((ssa && !ssa_arr) || (g_asym && !g_arr)) {
        free(ssa_arr); free(g_arr); return -2;
    }
    for (int iw = 0; iw < mie->n_wl; ++iw) {
        if (ssa_arr) ssa_arr[iw] = SPECTRAL_AT(mie, iw, 2);
        if (g_arr) g_arr[iw] = SPECTRAL_AT(mie, iw, 3);
    }
    if (ssa) *ssa = linterp(mie->wavelengths, ssa_arr, mie->n_wl, wavelength_um);
    if (g_asym) *g_asym = linterp(mie->wavelengths, g_arr, mie->n_wl, wavelength_um);
    free(ssa_arr); free(g_arr);
    return 0;
}

/* ─── build_aerosol_interpolators (combined vrt + ocean) ─── */

int mie_data_covers_bulk_wavelength(const mie_data_t* m, double wavelength_um)
{
    if (!m || m->n_wl < 2 || !m->wavelengths || !isfinite(wavelength_um)) return 0;
    return wavelength_um >= m->wavelengths[0] - 1.0e-12 &&
           wavelength_um <= m->wavelengths[m->n_wl - 1] + 1.0e-12;
}

int mie_data_covers_phase_wavelength(const mie_data_t* m, double wavelength_um)
{
    if (!m || m->n_phase_wl < 2 || !m->phase_wavelengths ||
        !isfinite(wavelength_um)) return 0;
    double lo = m->phase_wavelengths[0];
    double hi = m->phase_wavelengths[m->n_phase_wl - 1];
    if (lo > hi) { double tmp = lo; lo = hi; hi = tmp; }
    return wavelength_um >= lo - 1.0e-12 && wavelength_um <= hi + 1.0e-12;
}

int mie_data_covers_wavelength(const mie_data_t* m, double wavelength_um)
{
    return mie_data_covers_bulk_wavelength(m, wavelength_um) &&
           mie_data_covers_phase_wavelength(m, wavelength_um);
}

int mie_data_covers_range(const mie_data_t* m, double min_um, double max_um)
{
    return isfinite(min_um) && isfinite(max_um) && min_um <= max_um &&
           mie_data_covers_wavelength(m, min_um) &&
           mie_data_covers_wavelength(m, max_um);
}

int build_aerosol_interpolators(const mie_data_t* mie, double wavelength_um,
                                  aerosol_phase_interp_t* out)
{
    const int N = mie->n_ang;
    const int n_wl = mie->n_wl;
    const int n_phase_wl = mie->n_phase_wl;
    if (N < 2 || n_wl < 1 || n_phase_wl < 1) return -1;
    if (!mie_data_covers_wavelength(mie, wavelength_um)) {
        fprintf(stderr,
                "Mie wavelength %.9g um is outside the explicit bulk/phase range "
                "of the requested file\n",
                wavelength_um);
        return -2;
    }

    /* Sort angles ascending (OPAC files are 180→0 deg). */
    double* ang_asc = (double*)malloc(N * sizeof(double));
    double* P_at_wl = (double*)malloc(N * sizeof(double));
    int* order = (int*)malloc(N * sizeof(int));
    if (!ang_asc || !P_at_wl || !order) {
        free(ang_asc); free(P_at_wl); free(order); return -1;
    }
    for (int k = 0; k < N; ++k) { order[k] = k; ang_asc[k] = mie->angles[k]; }
    /* Insertion sort on (ang_asc, order) */
    for (int i = 1; i < N; ++i) {
        double t = ang_asc[i]; int vi = order[i]; int j = i - 1;
        while (j >= 0 && ang_asc[j] > t) {
            ang_asc[j + 1] = ang_asc[j]; order[j + 1] = order[j]; j--;
        }
        ang_asc[j + 1] = t; order[j + 1] = vi;
    }

    /* Build each of (P11, P12, P33): shape-preserving PCHIP over the PHASE
     * wavelength grid per angle, then PCHIP in angle.  This matches the OSOAA
     * conv_local spectral interpolation.  Spectral ssa/g below intentionally
     * retain their existing n_wl linear interpolation. */
    double* row = (double*)malloc((size_t)n_phase_wl * sizeof(double));
    if (!row) { free(ang_asc); free(P_at_wl); free(order); return -1; }
    const double* blocks[3] = { mie->P11, mie->P12, mie->P33 };
    pchip_t* pchips[3] = { &out->P11_pchip, &out->P12_pchip, &out->P33_pchip };
    for (int b = 0; b < 3; ++b) {
        for (int k = 0; k < N; ++k) {
            int k_orig = order[k];
            for (int l = 0; l < n_phase_wl; ++l) {
                row[l] = BLOCK_AT(blocks[b], n_phase_wl, k_orig, l);
            }
            P_at_wl[k] = pchip_interp(mie->phase_wavelengths, row, n_phase_wl, wavelength_um);
        }
        if (pchip_build(pchips[b], ang_asc, P_at_wl, N) != 0) {
            for (int bb = 0; bb < b; ++bb) pchip_free(pchips[bb]);
            free(row); free(ang_asc); free(P_at_wl); free(order);
            return -1;
        }
    }

    /* Diagnostic fields: ssa, g_asym (linear in wavelength) */
    double* ssa_arr = (double*)malloc(n_wl * sizeof(double));
    double* g_arr   = (double*)malloc(n_wl * sizeof(double));
    if (ssa_arr && g_arr) {
        for (int iw = 0; iw < n_wl; ++iw) {
            ssa_arr[iw] = SPECTRAL_AT(mie, iw, 2);
            g_arr[iw]   = SPECTRAL_AT(mie, iw, 3);
        }
        out->ssa     = linterp(mie->wavelengths, ssa_arr, n_wl, wavelength_um);
        out->g_asym  = linterp(mie->wavelengths, g_arr,   n_wl, wavelength_um);
    } else {
        out->ssa = 0.0;
        out->g_asym = 0.0;
    }
    free(ssa_arr); free(g_arr);

    /* bb/b ratio: numerical integration of P11 backward hemisphere on PCHIP */
    {
        int Nq = 200;
        double bb_int = 0.0, b_int = 0.0;
        for (int k = 0; k < Nq; ++k) {
            double theta = M_PI * (k + 0.5) / Nq;
            double sint  = sin(theta);
            double dtheta = M_PI / Nq;
            double p11 = pchip_eval(&out->P11_pchip, theta * 180.0 / M_PI);
            b_int  += p11 * sint * dtheta;
            if (theta > M_PI / 2.0) bb_int += p11 * sint * dtheta;
        }
        out->bb_b_ratio = (b_int > 0.0) ? bb_int / b_int : 0.0;
    }

    free(row); free(ang_asc); free(P_at_wl); free(order);
    return 0;
}

void eval_aerosol_phase(const aerosol_phase_interp_t* a, double theta_deg,
                          double* P11, double* P12, double* P33) {
    *P11 = pchip_eval(&a->P11_pchip, theta_deg);
    *P12 = pchip_eval(&a->P12_pchip, theta_deg);
    *P33 = pchip_eval(&a->P33_pchip, theta_deg);
}

void aerosol_phase_interp_free(aerosol_phase_interp_t* a) {
    pchip_free(&a->P11_pchip);
    pchip_free(&a->P12_pchip);
    pchip_free(&a->P33_pchip);
    memset(a, 0, sizeof(*a));
}
