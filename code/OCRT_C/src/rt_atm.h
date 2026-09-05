#ifndef OCRT_V2_RT_ATM_H
#define OCRT_V2_RT_ATM_H

#include <stdio.h>
#include "rt_types.h"
#include "rt_rayleigh.h"

/* Allocate atm arrays for given (n_layers, n_mu).
 * Returns 0 on success, -1 on alloc failure or bad args.
 * Caller must rt_atm_free() afterwards.
 *
 * Allocates: h, ch, xdel, ydel, rm_storage, gb_storage, xpl_storage.
 * Sets offset pointers rm/gb/xpl into their _storage backing arrays.
 */
int  rt_atm_alloc(rt_atm_t *atm, int n_layers, int n_mu);
void rt_atm_free(rt_atm_t *atm);

/* Build a Rayleigh-only atmosphere on a uniform layer subdivision.
 * atm must already be rt_atm_alloc()'d with the matching (n_layers, n_mu).
 *
 *   tau_total      total column Rayleigh optical depth (caller-computed
 *                  using the chosen model from rt_rayleigh.h)
 *   depol          raw depolarization factor δ (e.g. 0.0279)
 *   mu_sun         cos(SZA), in (0, 1]
 *   rayleigh_model metadata label of the model that produced tau_total.
 *                  OCRT production accepts BODHAINE_1999 only; the field is
 *                  retained for option-structure compatibility and dumps.
 *
 * Fills:
 *   h[k]    = k * tau_total / n_layers
 *   ch[k]   = exp(-h[k]/mu_sun) / 2
 *   xdel[k] = 0,  ydel[k] = 1
 *   rm[0]   = -mu_sun;  rm[+j], rm[-j] from Gauss-Legendre nodes
 *   gb[0]   = 0;        gb[±j]   from Gauss-Legendre weights
 *   xpl/xrl/xtl initially contain the m=0, l=2 scalar/spin-2 basis.
 *             The Fourier solver refreshes these work arrays for each m.
 *   beta0   = 1
 *   beta2   = 0.5 × ron, ron = 2(1-δ)/(2+δ)
 *
 * Returns 0 on success, -1 on bad args, -2 on quadrature failure.
 */
int rt_atm_build_rayleigh(rt_atm_t *atm, double tau_total,
                          double depol, double mu_sun,
                          rt_rayleigh_model_t rayleigh_model);

/* Build a zero-scattering U.S.-Standard-1962 altitude grid for direct-beam
 * extinction calculations.  This is used by the ocean-coupled PSSA helper
 * when Rayleigh and aerosol are both disabled but AFGL gas absorption is
 * present.  It establishes z_km_level[] and a zero optical-depth h[];
 * rt_atm_apply_gas_absorption() then fills the actual extinction profile. */
int rt_atm_build_extinction_grid_us62(rt_atm_t *atm, double mu_sun);

/* Build a mixed Aerosol+Rayleigh atmosphere on a uniform-τ layer subdivision.
 *
 * Caller pre-computes the aerosol expansion coefficients (e.g. via
 * rt_aerosol_compute_vector_legendre) and passes them in. This function
 * stores them in atm->{betal,gammal,alphal,zetal}_aer (allocated/copied
 * here; freed by rt_atm_free()).
 *
 * Layer mixing follows OCRT physical profiles:
 *   - Molecular/Rayleigh: US Standard Atmosphere 1962 pressure column.
 *   - Aerosol: exp(-z/H), production H=2 km.
 *
 * Inputs:
 *   atm          pre-allocated rt_atm_t
 *   tau_R        column Rayleigh optical depth
 *   depol        Rayleigh depolarization factor
 *   tau_a        column aerosol optical depth at λ
 *   ssa_a        aerosol single scattering albedo at λ
 *   L_max        max Legendre order in arrays (≥ 2)
 *   betal_aer / gammal_aer / alphal_aer / zetal_aer  size L_max+1 each
 *   mu_sun       cos(SZA), in (0, 1]
 *   rayleigh_model
 *   aer_h_km     positive aerosol scale height in km; production 2.0.
 */
int rt_atm_build_aerosol_rayleigh(
    rt_atm_t *atm,
    double tau_R, double depol,
    double tau_a, double ssa_a,
    int L_max,
    const double *betal_aer, const double *gammal_aer,
    const double *alphal_aer, const double *zetal_aer,
    double mu_sun,
    rt_rayleigh_model_t rayleigh_model,
    double aer_h_km
);

/* Human-readable dump of atm contents to fp (stdout typical). */
void rt_atm_dump(const rt_atm_t *atm, FILE *fp);

/* v1.02: Apply gas absorption to a pre-built atmosphere (SOS-integrated).
 *
 * Per-layer mechanics:
 *   1. From atm->z_km_level[k], compute layer altitude band [z_lo, z_hi]
 *   2. From AFGL profile (in abs_state), compute N_gas_layer (mol/cm²) per band
 *   3. From layered xsec table, get σ_gas(λ, T_k, P_k)
 *   4. dt_abs_layer[k] = Σ_g σ_g · N_g_layer
 *   5. h_new[k] = h_old[k] + Σ_{j<k} dt_abs_j
 *   6. ch[k] = exp(-h_new[k] / mu_sun) / 2
 *   7. xdel[k], ydel[k] re-normalized: dt_total_new = dt_ray + dt_aer + dt_abs
 *      keeping aerosol/Rayleigh scattering ODs absolute, only changing denom.
 *
 * Forward declaration here; concrete type rt_absorption_t in rt_absorption.h.
 * Caller must include both headers and pass a fully-initialized rt_absorption_t.
 *
 * Returns 0 on success, -1 on bad input.
 */
struct rt_absorption_state;
int rt_atm_apply_gas_absorption(rt_atm_t *atm,
                                 const struct rt_absorption_state *abs_state,
                                 double wavelength_nm);

#endif /* OCRT_V2_RT_ATM_H */
