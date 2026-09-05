/* shared/mie_io.h
 *
 * Mie data container, OPAC .mie file I/O, and PCHIP interpolators
 * for aerosol/particle phase functions. Shared by vrt_solver and
 * ocean_solver.
 *
 * Design:
 *   - mie_data_t holds wavelengths, angles, spectral coefficients, and
 *     phase matrix blocks (P11, P12, P33). Identical layout to OPAC files.
 *   - aerosol_phase_interp_t holds PCHIP interpolators for one wavelength,
 *     plus diagnostic quantities (ssa, g_asym, bb_b_ratio) used by ocean
 *     IOP code. vrt does not use the diagnostic fields but they cost little.
 *
 * File format: see vrt_solver.c::write_mie_file. read_mie_file is robust
 * against whitespace/comment variations.
 */
#ifndef OCRT_SHARED_MIE_IO_H
#define OCRT_SHARED_MIE_IO_H

#include <stddef.h>

#include "numerics.h"   /* pchip_t */

/* ─── mie_data_t ──────────────────────────────────────── */

typedef struct {
    int n_ang;
    int n_wl;                  /* # wavelengths in the IOP/spectral grid */
    int n_phase_wl;            /* # wavelengths in the phase-matrix grid; == n_wl
                                  for single-grid (legacy) files, may be coarser
                                  for OPAC-style files with a separate TETA grid */
    double* wavelengths;       /* [n_wl] in μm — IOP/spectral grid */
    double* phase_wavelengths; /* [n_phase_wl] in μm — phase-matrix grid
                                  (a copy of wavelengths when single-grid) */
    double* angles;            /* [n_ang] in degrees */
    double* spectral;          /* [n_wl * 6] row-major (col 0=Cext_norm,
                                   1=Csca_norm, 2=ssa, 3=g, 4=Cext, 5=Csca) */
    double* P11;               /* [n_ang * n_phase_wl] row-major */
    double* P12;               /* [n_ang * n_phase_wl] */
    double* P33;               /* [n_ang * n_phase_wl] */
} mie_data_t;

#define SPECTRAL_AT(m, iwl, col) ((m)->spectral[(iwl) * 6 + (col)])
#define BLOCK_AT(arr, n_wl, ia, iw) ((arr)[(ia) * (n_wl) + (iw)])

void mie_data_init(mie_data_t* m, int n_ang, int n_wl, int n_phase_wl);
void mie_data_free(mie_data_t* m);

/* Read .mie file produced by write_mie_file (OPAC layout).
 * Returns 0 on success, -1 on failure. */
int read_mie_file(const char* path, mie_data_t* out);

/* Worker-private parsed-model cache.  The returned pointer is borrowed and
 * remains valid until a later cache insertion in the same worker evicts it or
 * mie_model_cache_clear_worker() is called.  The complete .mie model (all
 * bulk and phase wavelengths) is retained in memory.
 *
 * Return values:
 *   0  success (possibly a cache hit)
 *   1  cache disabled by OCRT_DISABLE_MIE_MODEL_CACHE; caller should parse a
 *      temporary model with read_mie_file()
 *  <0  error
 */
int mie_model_cache_get(const char* path, const mie_data_t** out, int* cache_hit);
void mie_model_cache_clear_worker(void);
size_t mie_model_cache_resident_bytes(void);
unsigned long long mie_model_cache_load_count(void);
unsigned long long mie_model_cache_hit_count(void);


/* Exact spectral-coverage checks.  Runtime interpolation is allowed only when
 * the requested wavelength lies inside both the bulk and phase grids. */
int mie_data_covers_bulk_wavelength(const mie_data_t* m, double wavelength_um);
int mie_data_covers_phase_wavelength(const mie_data_t* m, double wavelength_um);
int mie_data_covers_wavelength(const mie_data_t* m, double wavelength_um);
int mie_data_covers_range(const mie_data_t* m, double min_um, double max_um);

/* ─── aerosol_phase_interp_t ─────────────────────────── */

typedef struct {
    pchip_t P11_pchip;
    pchip_t P12_pchip;
    pchip_t P33_pchip;
    /* Diagnostic fields (used by ocean IOP, ignored by vrt) */
    double  ssa;
    double  g_asym;
    double  bb_b_ratio;
} aerosol_phase_interp_t;

/* Build PCHIP interpolators (P11/P12/P33 vs theta) at given wavelength_um.
 * Also computes ssa, g_asym (linear interp from spectral table) and
 * bb_b_ratio (numerical integration of P11 backward hemisphere).
 * Caller must free via aerosol_phase_interp_free.
 * Returns 0 on success. */
int build_aerosol_interpolators(const mie_data_t* mie, double wavelength_um,
                                  aerosol_phase_interp_t* out);

/* Evaluate only the phase-matrix wavelength interpolation at every native
 * angle node.  No angle interpolation is performed.  This is the canonical
 * input to the Stage-3B theta-linear direct-LUT kernel.  Output arrays must
 * have mie->n_ang elements and retain the file angle order. */
int mie_phase_nodes_at_wavelength_pchip(const mie_data_t* mie,
                                        double wavelength_um,
                                        double* P11,
                                        double* P12,
                                        double* P33);

/* Evaluate P11/P12/P33 at every native angle node using one common linear
 * wavelength bracket/weight.  This is the Stage-3B direct-LUT spectral
 * contract: it preserves P11>=0 and |P12|,|P33|<=P11 between wavelength
 * nodes when the two endpoint matrices satisfy those bounds. */
int mie_phase_nodes_at_wavelength_linear(const mie_data_t* mie,
                                         double wavelength_um,
                                         double* P11,
                                         double* P12,
                                         double* P33);

/* Linear bulk-grid diagnostics used when a full angle interpolator is not
 * constructed.  Any output pointer may be NULL. */
int mie_bulk_diagnostics_at_wavelength(const mie_data_t* mie,
                                       double wavelength_um,
                                       double* ssa,
                                       double* g_asym);

/* Evaluate phase matrix elements at single scattering angle Theta (degrees). */
void eval_aerosol_phase(const aerosol_phase_interp_t* a, double theta_deg,
                          double* P11, double* P12, double* P33);

void aerosol_phase_interp_free(aerosol_phase_interp_t* a);

#endif /* OCRT_SHARED_MIE_IO_H */
