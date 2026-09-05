#ifndef OCRT_RT_MOLECULAR_PROFILE_H
#define OCRT_RT_MOLECULAR_PROFILE_H

/* OCRT molecular-column vertical profile.
 *
 * Production profile: U.S. Standard Atmosphere 1962 pressure-altitude
 * anchors documented by the OCRT methods specification.  The normalized
 * molecular column above altitude z is evaluated hydrostatically as
 *
 *   F(z) = [P_US62(z) / g_45(z)] / [P0 / g_45(0)]
 *
 * so F(0)=1 and F(z_top)=P_top/g_top divided by the sea-level column.
 * The finite residual above the highest tabulated level is assigned to the
 * top boundary; callers that discretize a complete column should use
 * rt_molecular_us62_layer_fraction(), which closes the column exactly.
 *
 * No external reference-code source is used by this module.  Its inputs are
 * the public U.S. Standard Atmosphere 1962 anchor table reproduced in the
 * OCRT methods document; interpolation is linear in log(P) versus altitude.
 */

#ifdef __cplusplus
extern "C" {
#endif

double rt_molecular_us62_top_km(void);
double rt_molecular_us62_pressure_pa(double altitude_km);
double rt_molecular_us62_fraction_above(double altitude_km);

/* Solver-grid cumulative fraction, renormalized so the finite US62 model
 * interval is exactly F_grid(0)=1 and F_grid(z_top)=0. */
double rt_molecular_us62_grid_fraction_above(double altitude_km);
double rt_molecular_us62_altitude_from_grid_fraction(double fraction_above);

/* Fraction of the complete molecular column contained in [z_lo,z_hi], with
 * z_hi >= z_lo.  The residual column above the US62 top anchor is included
 * when z_hi reaches the model top, which guarantees exact column closure. */
double rt_molecular_us62_layer_fraction(double z_lo_km, double z_hi_km);

#ifdef __cplusplus
}
#endif

#endif /* OCRT_RT_MOLECULAR_PROFILE_H */
