#ifndef OCRT_V2_RT_RAA_CONVENTION_H
#define OCRT_V2_RT_RAA_CONVENTION_H

#include <math.h>
#include "rt_types.h"

/*
 * Canonical OCRT public relative-azimuth convention (CLI / LUT metadata).
 *
 * Use the following operational definition; it avoids the ambiguity between
 * a solar *look* azimuth and the downward solar photon-propagation azimuth:
 *
 *   RAA_OCRT = 180 deg : direct rough-surface glint/specular branch when
 *                        VZA == SZA.
 *   RAA_OCRT =   0 deg : the opposite principal-plane branch.
 *   RAA_OCRT =  90 deg : convention-invariant quadrature branch used by the
 *                        OCRT--OSOAA comparison harness.
 *
 * For comparison with OSOAA, one physical mapping may be written in either
 * of the following mathematically equivalent forms:
 *
 *   A. Phi_OSOAA = (180 deg - RAA_OCRT) mod 360 deg, then reverse Stokes U;
 *   B. Phi_OSOAA = (RAA_OCRT + 180 deg) mod 360 deg, with no U reversal.
 *
 * Use exactly one form, never both.  The equivalence follows from the
 * axisymmetric identities I,Q(phi)=I,Q(360-phi) and U(phi)=-U(360-phi).
 * OCRT water outputs Rrs(0+) and rrs(0-) now use the same public RAA
 * reconstruction argument as the atmospheric output.  The legacy identity
 * mapping Phi_OSOAA,water=RAA_OCRT compensated a mirrored OCRT water output
 * and must not be used with this corrected source.
 *
 * Individual OCRT kernels use different internal azimuth coordinates.  All
 * public RAA conversions must go through these helpers.  Do not reintroduce
 * raw `raa_deg * pi/180` or `(180 - raa_deg)` expressions in solver code unless
 * the expression is explicitly a local analytic formula and covered by the
 * RAA audit test.
 */

static inline double rt_deg_to_rad(double deg) {
    return deg * RT_F_SOLAR_PI / 180.0;
}

static inline double rt_norm_raa_0_360(double raa_deg) {
    double x = fmod(raa_deg, 360.0);
    if (x < 0.0) x += 360.0;
    if (fabs(x - 360.0) < 1.0e-12) x = 0.0;
    return x;
}

/* Argument passed to the atmospheric Fourier reconstruction API.
 * The helper returns the public OCRT RAA in radians.  The stored atmospheric
 * modal field is in propagation azimuth and rt_fourier_reconstruct_* applies
 * the required +pi coordinate shift internally.  Therefore this return value
 * is not, by itself, the physical propagation azimuth.
 */
static inline double rt_raa_to_atm_fourier_phi(double raa_deg) {
    return rt_deg_to_rad(raa_deg);
}

/* Atmospheric solar-down -> view-up scattering-angle formula in the
 * public OCRT convention:
 *   RAA=0   (opposite principal-plane branch): cosTheta = -mu_s mu_v + ss sv
 *   RAA=180 (direct-glint branch)             : cosTheta = -mu_s mu_v - ss sv
 * For SZA=VZA, RAA=180 is exact backscatter (Theta=180 deg), consistent with
 * the specular rough-surface geometry.
 */
static inline double rt_atm_scatter_cos_from_public_raa(double mu_s, double mu_v,
                                                        double raa_deg) {
    const double ss = sqrt(fmax(0.0, 1.0 - mu_s * mu_s));
    const double sv = sqrt(fmax(0.0, 1.0 - mu_v * mu_v));
    double c = -mu_s * mu_v + ss * sv * cos(rt_deg_to_rad(raa_deg));
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    return c;
}

/* Water-SOS direct-scatter correction.  OCRT's water modal field uses the
 * internal view azimuth phi_water = pi - RAA_OCRT.  Substituting this into
 * cosTheta = -mu_s mu_v + ss sv cos(phi_water) gives the minus sign below.
 * This is an OCRT internal coordinate; it must not be confused with the OSOAA
 * water-output Phi used by the external comparison harness.
 */
static inline double rt_water_scatter_cos_from_public_raa(double mu_s, double mu_v,
                                                          double raa_deg) {
    const double ss = sqrt(fmax(0.0, 1.0 - mu_s * mu_s));
    const double sv = sqrt(fmax(0.0, 1.0 - mu_v * mu_v));
    double c = -mu_s * mu_v - ss * sv * cos(rt_deg_to_rad(raa_deg));
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    return c;
}

/* Local water single-scatter propagation-vector azimuth.
 *
 * This pi-RAA coordinate is used only when an explicit incoming/outgoing
 * propagation-vector pair is constructed for the exact single-scatter
 * correction.  It is NOT the reconstruction argument for reported water
 * Fourier outputs.  Reported Rrs(0+) and rrs(0-) must use
 * rt_raa_to_atm_fourier_phi(), matching the public OCRT RAA.
 */
static inline double rt_raa_to_water_scatter_phi(double raa_deg) {
    return rt_deg_to_rad(180.0 - raa_deg);
}

/* Direct rough-Fresnel surface helper.  Its local surface-view azimuth has
 * phi_surface=0 at RAA_OCRT=180, so VZA=SZA lies on the direct-glint peak.
 * It equals -phi_water modulo 2pi and remains a distinct helper to keep the
 * surface and water-SOS coordinates explicit at call sites.
 */
static inline double rt_raa_to_surface_view_phi(double raa_deg) {
    return rt_deg_to_rad(raa_deg - 180.0);
}

#endif /* OCRT_V2_RT_RAA_CONVENTION_H */
