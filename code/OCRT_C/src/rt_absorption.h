#ifndef OCRT_RT_ABSORPTION_H
#define OCRT_RT_ABSORPTION_H

/* OCRT v1.01 — Atmospheric absorption module (post-RT Beer-Lambert correction).
 *
 * Status: MINIMAL VIABLE implementation. Applies a single column absorption
 * optical depth (τ_abs total) to RT output via Beer-Lambert two-pass:
 *   ρ_obs(μ_v, μ_sun, φ) = ρ_RT(μ_v, μ_sun, φ) · exp(-τ_abs · (1/μ_sun + 1/μ_v))
 *   T_dir_dn_obs = T_dir_dn_RT · exp(-τ_abs / μ_sun)
 *
 * This is an approximation. Full multiple-scattering coupling with absorption
 * requires per-layer τ_abs_layer in the SOS solver (ToDo: VLIDORT-style).
 *
 * Inputs:
 *   - AFGL atmospheric profile (one of 6 standard atmospheres)
 *   - HITRAN absorption cross section tables (xsec_<gas>.dat)
 *   - Optional per-gas column abundance override
 *
 * The 6 gases handled: H2O, O3, O2, CO2, NO2, CH4.
 */

#include <stddef.h>

typedef enum {
    RT_AFGL_USSTD76 = 0,
    RT_AFGL_TROPICAL = 1,
    RT_AFGL_MLSUMM = 2,
    RT_AFGL_MLWINT = 3,
    RT_AFGL_SASUMM = 4,
    RT_AFGL_SAWINT = 5,
    RT_AFGL_USERDEF = 6
} rt_afgl_profile_t;

#define RT_N_GAS 6
enum { RT_GAS_H2O = 0, RT_GAS_O3, RT_GAS_O2, RT_GAS_CO2, RT_GAS_NO2, RT_GAS_CH4 };

extern const char *RT_GAS_NAMES[RT_N_GAS];   /* {"h2o","o3","o2","co2","no2","ch4"} */

/* AFGL profile (loaded from inputs/afgl_atm/afgl_<name>.dat) */
typedef struct {
    int n_levels;
    double *z_km;            /* altitude grid (km), size n_levels */
    double *P_mbar;          /* pressure (mbar) */
    double *T_K;              /* temperature (K) */
    double *n_total_cm3;      /* total number density (cm^-3) */
    double *mr[RT_N_GAS];     /* mixing ratio (ppmv) for each gas, size n_levels */
    /* Default column abundances computed from profile (cm^-2 = molecules/cm²) */
    double column_default[RT_N_GAS];
} rt_afgl_atm_t;

/* HITRAN cross section LUT (one per gas).
 *
 * v1.02: supports BOTH legacy single-layer format AND layered format.
 *   legacy:  σ(λ) only, single T_ref/P_ref → n_layers=1
 *   layered: σ(layer_idx, λ) computed at AFGL layer-specific T_k, P_k
 *   The reader auto-detects format via "# n_layers" header tag.
 */
typedef struct {
    int n_wl;
    int n_layers;             /* 1 for legacy single-T/P, >1 for AFGL-layered */
    double *wl_nm;            /* wavelength grid (nm), size n_wl */
    double *sigma_cm2;        /* size n_layers × n_wl, layer-major
                               * sigma_cm2[k*n_wl + i] = σ at layer k, wavelength i */
    double T_ref_K;           /* (if legacy) reference T */
    double P_ref_atm;         /* (if legacy) reference P */
    int    available;         /* 1 if file was read successfully */
} rt_xsec_t;

/* Absorption runtime state */
typedef struct {
    rt_afgl_atm_t    atm;
    rt_xsec_t        xsec[RT_N_GAS];
    /* Effective column abundance (molecules/cm²) — default or user-override */
    double           column_eff[RT_N_GAS];
    /* Flags: 1 if user overrode the column for this gas */
    int              column_overridden[RT_N_GAS];
    int              initialized;
} rt_absorption_t;

/* === API === */

/* Load AFGL atmosphere from inputs/afgl_atm/afgl_<name>.dat
 * Returns 0 on success, non-zero on file error.
 * On success, computes column_default[] from integrated profile.
 */
int rt_afgl_load(const char *afgl_dir, rt_afgl_profile_t profile,
                 rt_afgl_atm_t *atm);

void rt_afgl_free(rt_afgl_atm_t *atm);

/* Load HITRAN cross section table from inputs/xsec/xsec_<gas>.dat
 * Sets xsec->available = 0 if file missing (silent, zero absorption used).
 * Returns 0 on success, non-zero on file error (with available=0).
 */
int rt_xsec_load(const char *xsec_dir, int gas_idx, rt_xsec_t *xsec);

void rt_xsec_free(rt_xsec_t *xsec);

/* Interpolate σ(λ) at given wavelength. Returns 0 if xsec not available. */
double rt_xsec_interp(const rt_xsec_t *xsec, double wavelength_nm);

/* Layered interpolation: σ at (layer, λ). Falls back to layer-0 if not layered. */
double rt_xsec_interp_layer(const rt_xsec_t *xsec, int layer_idx,
                             double wavelength_nm);

/* v1.02: BI-LINEAR interpolation in (layer altitude, wavelength) space.
 *
 * If xsec file stores σ at AFGL z-grid points, this function returns
 * σ at arbitrary z_km, wavelength_nm by *linear interpolation between the
 * two bracketing AFGL z-levels* AND wavelength interpolation.
 *
 * NEAREST-NEIGHBOR IS NEVER USED. Linear bracketing only.
 *
 * z_grid_km: altitude (km) of each xsec layer (size xsec->n_layers).
 *            Provided externally because xsec file format only stores the
 *            σ matrix; layer altitudes are in the AFGL profile.
 */
double rt_xsec_interp_layer_z(const rt_xsec_t *xsec,
                               const double *z_grid_km,
                               double z_km, double wavelength_nm);

/* AFGL altitude → cumulative gas column number density (mol/cm²).
 * cumulative_column(z) = ∫_z^∞ n_gas(z') dz' (mol/cm² above altitude z)
 * Linear interpolation between AFGL z grid points. */
double rt_afgl_cumulative_column(const rt_afgl_atm_t *atm, int gas_idx,
                                  double z_km);

/* Per-layer gas column within altitude band [z_lo, z_hi] (mol/cm²). */
double rt_afgl_layer_column(const rt_afgl_atm_t *atm, int gas_idx,
                             double z_lo_km, double z_hi_km);

/* Interpolate AFGL T (K) and P (mbar) at given altitude. */
double rt_afgl_T_at(const rt_afgl_atm_t *atm, double z_km);
double rt_afgl_P_at(const rt_afgl_atm_t *atm, double z_km);

/* Initialize full absorption state.
 *   afgl_dir, xsec_dir = paths to data directories
 *   profile = which AFGL atmosphere
 *   gas_column_override[i] = -1.0 (use default) or non-negative (override, mol/cm²)
 */
int rt_absorption_init(const char *afgl_dir, const char *xsec_dir,
                       rt_afgl_profile_t profile,
                       const double gas_column_override[RT_N_GAS],
                       rt_absorption_t *abs_state);

void rt_absorption_free(rt_absorption_t *abs_state);

/* Compute total atmospheric absorption optical depth at given wavelength.
 *   τ_abs(λ) = Σ_gas σ_gas(λ) · N_gas_column
 */
double rt_absorption_tau_total(const rt_absorption_t *abs_state,
                                double wavelength_nm);

/* Per-gas τ_abs at wavelength — useful for sensitivity analysis. */
double rt_absorption_tau_per_gas(const rt_absorption_t *abs_state,
                                  int gas_idx, double wavelength_nm);

/* Beer-Lambert correction factor for upward radiance at TOA (two-pass).
 *   factor = exp(-τ_abs · (1/μ_sun + 1/μ_v))
 * This is a post-RT approximation. Full layer-by-layer SOS coupling = ToDo.
 */
double rt_absorption_correction_two_pass(double tau_abs,
                                          double mu_sun, double mu_v);

/* One-pass downward correction (for T_dir_dn etc.):
 *   factor = exp(-τ_abs / μ_sun)
 */
double rt_absorption_correction_one_pass(double tau_abs, double mu_sun);

#endif /* OCRT_RT_ABSORPTION_H */
