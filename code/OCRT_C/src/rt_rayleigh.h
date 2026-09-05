#ifndef OCRT_V2_RT_RAYLEIGH_H
#define OCRT_V2_RT_RAYLEIGH_H

#include "rt_types.h"

/* OCRT Rayleigh optical depth.
 *
 * The production implementation follows the equations and constants of
 * Bodhaine et al. (1999): refractive-index dispersion, wavelength-dependent
 * King correction, hydrostatic molecular column, and local gravity.  The
 * model-selecting wrappers remain temporarily for source compatibility, but
 * only RT_RAYLEIGH_MODEL_BODHAINE_1999 is accepted.
 */

/* WGS84 surface gravity [m/s²] at given latitude/altitude. */
double rt_rayleigh_gravity_wgs84(double latitude_deg, double altitude_m);

/* King correction factor F(λ) for air, CO2-dependent (Bodhaine 1999). */
double rt_rayleigh_king_factor_air(double wavelength_nm, double co2_ppm);

/* Rayleigh depolarization ratio ρ(λ) = 6(F-1)/(7F+3) from King factor. */
double rt_rayleigh_depolarization_ratio(double wavelength_nm, double co2_ppm);

/* Bodhaine 1999 FIRST-PRINCIPLES tau (v1.09+; production default):
 *   Peck&Reeder-1972 n(λ) with CO2 correction, exact Lorentz-Lorenz
 *   cross-section, WAVELENGTH-DEPENDENT King factor (Bates 1984 species
 *   dispersions), hydrostatic column P*A/(m_a*g) with g at the [B99]
 *   mass-weighted column altitude.  See rt_rayleigh.c.
 *   The 4-arg form fixes CO2 = 360 ppm ([B99] Table 3 reference). */
double rt_rayleigh_tau_bodhaine1999_full(double wavelength_nm,
                                         double pressure_hpa,
                                         double latitude_deg,
                                         double altitude_m);
double rt_rayleigh_tau_bodhaine1999_full_co2(double wavelength_nm,
                                             double pressure_hpa,
                                             double latitude_deg,
                                             double altitude_m,
                                             double co2_ppm);

/* Dispatchers (sea-level standard: 1013.25 hPa, 45°N, 0 m). */
double rt_rayleigh_tau(double wavelength_nm);   /* default = Bodhaine */
double rt_rayleigh_tau_model(double wavelength_nm,
                             rt_rayleigh_model_t model);
double rt_rayleigh_tau_model_full(double wavelength_nm,
                                  double pressure_hpa,
                                  double latitude_deg,
                                  double altitude_m,
                                  rt_rayleigh_model_t model);

/* Backward-compatible alias for the Bodhaine full form (kept so existing
 * callers in main.c / dumps / tests don't all need updating in one go).
 * New code should use rt_rayleigh_tau_bodhaine1999_full() explicitly.
 */
double rt_rayleigh_tau_full(double wavelength_nm, double pressure_hpa,
                            double latitude_deg, double altitude_m);

#endif /* OCRT_V2_RT_RAYLEIGH_H */
