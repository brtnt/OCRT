/* =============================================================================
 * rt_rayleigh.c — OCRT molecular (Rayleigh) column optical depth
 * =============================================================================
 *
 * Implements the first-principles method of Bodhaine et al. (1999):
 * refractive-index dispersion, the wavelength-dependent dry-air King factor,
 * an exact Lorentz-Lorenz molecular scattering cross section, hydrostatic
 * conversion from pressure to molecular column, and local-gravity scaling.
 *
 * The implementation is expressed directly from the published equations and
 * physical constants.  Reference-code comparison formulas that were not part
 * of the OCRT production model have been removed from this core module.
 * ========================================================================== */

#include <math.h>
#include "rt_rayleigh.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- shared physical helpers ---------- */

/* Local gravitational acceleration g(phi, z)  [m s^-2].
 *
 * Equation (normal gravity on the reference ellipsoid + free-air reduction):
 *
 *   g0(phi) = 9.780327 * (1 + 0.0053024*sin^2(phi) - 0.0000058*sin^2(2*phi))
 *   g(phi,z) = g0(phi) - 3.086e-6 * z
 *
 * Symbols / units:
 *   phi  = geographic latitude [rad]   (input latitude_deg is degrees)
 *   z    = altitude above sea level [m]  (input altitude_m)
 *   g0   = sea-level normal gravity [m s^-2]; 9.780327 is the equatorial
 *          value g_e of the 1967/1980 international gravity formula family
 *   3.086e-6 [s^-2] = free-air vertical gravity gradient dg/dz
 *
 * Reference: the same formula family used by [B99] to convert surface
 * pressure into molecular column ([B99] uses List 1968; differences vs this
 * form are < 1e-5 relative — see the tau routine notes).
 * Validity: |z| small compared to Earth radius (linear free-air term).
 *
 * Role in tau: tau scales as 1/g because the number of molecules above a
 * barometer reading P is N_col = P / (m_air * g); see the *_full() models.
 */
double rt_rayleigh_gravity_wgs84(double latitude_deg, double altitude_m) {
    double lat = latitude_deg * M_PI / 180.0;
    double s   = sin(lat) * sin(lat);            /* sin^2(phi)   */
    double s2  = sin(2.0 * lat) * sin(2.0 * lat);/* sin^2(2 phi) */
    double g0  = 9.780327 * (1.0 + 0.0053024 * s - 0.0000058 * s2);
    return g0 - 3.086e-6 * altitude_m;
}

/* King correction factor F(air)(lambda) of moist-free air  [dimensionless].
 *
 * The King factor corrects the ideal-dipole Rayleigh cross section for the
 * anisotropy of real molecules.  Composition-weighted mixing rule ([B99]
 * Eq. 23, species factors from [B84] as quoted in [B99] Eqs. 5-6):
 *
 *   F(N2)(lambda)  = 1.034 + 3.17e-4 / lambda_um^2
 *   F(O2)(lambda)  = 1.096 + 1.385e-3 / lambda_um^2 + 1.448e-4 / lambda_um^4
 *   F(Ar)          = 1.00                  (isotropic monatomic gas)
 *   F(CO2)         = 1.15                  (constant approximation, [B99])
 *
 *   F(air) = ( 78.084*F(N2) + 20.946*F(O2) + 0.934*F(Ar) + c_CO2*F(CO2) )
 *            / ( 78.084 + 20.946 + 0.934 + c_CO2 )
 *
 * Symbols / units:
 *   lambda_um = wavelength in micrometres (input wavelength_nm / 1000)
 *   c_CO2     = CO2 molar concentration in PERCENT by volume
 *             = co2_ppm / 1e4    (e.g. 400 ppm -> 0.04 %)
 *   78.084 / 20.946 / 0.934 = N2 / O2 / Ar molar percentages of dry air
 *
 * Validity: UV-A through NIR (the [B84] dispersion fits); dry air.
 * Note the denominator is NOT forced to 100%: adding CO2 slightly renormalizes
 * the mixture exactly as in [B99].
 */
double rt_rayleigh_king_factor_air(double wavelength_nm, double co2_ppm) {
    double wl_um = wavelength_nm / 1000.0;
    double wl2   = wl_um * wl_um;
    double wl4   = wl2 * wl2;
    double f_n2  = 1.034 + 3.17e-4 / wl2;
    double f_o2  = 1.096 + 1.385e-3 / wl2 + 1.448e-4 / wl4;
    double f_ar  = 1.0;
    double f_co2 = 1.15;
    double x_co2_pct = co2_ppm / 1.0e4;          /* ppm -> percent */
    double denom = 78.084 + 20.946 + 0.934 + x_co2_pct;
    return (78.084 * f_n2 + 20.946 * f_o2 + 0.934 * f_ar
            + x_co2_pct * f_co2) / denom;
}

/* Depolarization ratio rho(lambda) of air  [dimensionless].
 *
 * Exact algebraic inversion of the King factor definition
 *   F = (6 + 3*rho) / (6 - 7*rho)        ([B99] Eq. 3)
 * giving
 *   rho = 6*(F - 1) / (7*F + 3)
 *
 * Used by the atmospheric RT side to build the depolarized Rayleigh phase
 * matrix; NOT used inside the tau models below (Model 1 uses F(air)
 * directly).
 */
double rt_rayleigh_depolarization_ratio(double wavelength_nm, double co2_ppm) {
    double F = rt_rayleigh_king_factor_air(wavelength_nm, co2_ppm);
    return 6.0 * (F - 1.0) / (7.0 * F + 3.0);
}

/* ---------- Model 1: Bodhaine 1999 first-principles (production default) --- */

/* Refractive index of dry air, Bodhaine et al. (1999) prescription  [n - 1].
 *
 * Step 1 — Peck & Reeder (1972) dispersion of STANDARD air
 * (15 degC, 101325 Pa, 300 ppm CO2), [B99] Eq. (18):
 *
 *   (n_300 - 1) * 1e8 = 8060.51 + 2480990 / (132.274  - sigma^2)
 *                               + 17455.7 / (39.32957 - sigma^2)
 *   sigma = 1 / lambda_um   [um^-1]
 *
 * Step 2 — CO2 concentration correction, [B99] Eq. (19):
 *
 *   (n_CO2 - 1) = (n_300 - 1) * (1 + 0.54 * (C_CO2 - 0.0003))
 *   C_CO2 = CO2 concentration in parts per volume (ppv; 360 ppm -> 3.6e-4)
 *
 * Validity: Peck & Reeder fit range ~0.23-1.69 um.
 * CONSISTENCY: this n is defined at the same standard conditions as N_s
 * below (288.15 K, 1013.25 hPa); the pair must not be mixed with other
 * (T, P) references (sigma_R depends on ((n^2-1)/N_s)^2 to leading order).
 */
static double bodhaine1999_n_minus_1(double wavelength_nm, double co2_ppm) {
    double wl_um  = wavelength_nm / 1000.0;
    double sigma2 = 1.0 / (wl_um * wl_um);
    double n300_m1 = (8060.51
                      + 2480990.0 / (132.274  - sigma2)
                      + 17455.7   / (39.32957 - sigma2)) * 1.0e-8;
    double c_ppv = co2_ppm * 1.0e-6;
    return n300_m1 * (1.0 + 0.54 * (c_ppv - 0.0003));
}

/* Rayleigh optical depth, Bodhaine et al. (1999) first-principles route
 * (the method the paper RECOMMENDS over curve fits)  [dimensionless].
 *
 * Step 1 — scattering cross section per molecule, [B99] Eq. (2)
 * (exact Lorentz-Lorenz form, NOT the small-(n-1) limit):
 *
 *   sigma_R(lambda) = 24 * pi^3 * (n^2 - 1)^2
 *                     / ( lambda_cm^4 * N_s^2 * (n^2 + 2)^2 )
 *                     * F_air(lambda, C_CO2)                       [cm^2]
 *
 *   F_air = WAVELENGTH-DEPENDENT King factor from the per-species Bates
 *   (1984) dispersions, rt_rayleigh_king_factor_air() above ([B99] Eqs.
 *   (5), (6), (23)).  This is the paper's key correction over legacy
 *   constant-depolarization treatments (Penndorf 1.0608 / Young 1.0480):
 *   F_air runs from ~1.081 at 200 nm to ~1.047 in the NIR.
 *
 *   N_s = molecular number density of standard air, [B99]:
 *     N_s = (A / 22.4141 L*mol^-1) * (273.15 / 288.15) * 1e-3
 *         = 2.546899e19 cm^-3       (288.15 K, 1013.25 hPa)
 *   A   = 6.0221367e23 mol^-1 (Avogadro constant as used by [B99];
 *         the CODATA update changes tau by ~4e-6 relative — negligible,
 *         paper value kept for Table-3 reproducibility).
 *
 * Step 2 — molecular column above the site, [B99] Eqs. (4)-(5):
 *
 *   tau = sigma_R * P * A / (m_a * g)
 *
 *   P   = site pressure in dyn cm^-2  (hPa * 1000)
 *   m_a = mean molecular weight of dry air [g mol^-1], [B99]:
 *         m_a = 15.0556 * C_CO2 + 28.9595     (C_CO2 in ppv)
 *   g   = gravity REPRESENTATIVE OF THE AIR COLUMN, not the surface:
 *         [B99] evaluates g at the mass-weighted column altitude
 *
 *           z_c = 0.73737 * z_site + 5517.56        [m]   ([B99] appendix)
 *
 *         (sea-level site -> z_c ~= 5518 m; this raises tau by ~+0.174%
 *          relative to naive surface-g evaluation, since tau ~ 1/g and
 *          dg/dz = -3.086e-6 m s^-2 per m).  VALIDATION: with this term the
 *          routine reproduces the NOAA NEUBrew tabulation of [B99] at the
 *          five Brewer wavelengths 306.3-320.1 nm to <0.02% (see v1.09
 *          changelog); without it there is a flat -0.17% offset.
 *         Units: g in cm s^-2 (gal) = 100 * rt_rayleigh_gravity_wgs84().
 *         NOTE: [B99] uses the List (1968) gravity formula including
 *         quadratic/cubic altitude terms; our helper is the international-
 *         gravity-formula family with the linear free-air term only.
 *         Differences: <1e-5 relative in g0 at 45 deg; the quadratic
 *         altitude term at z_c ~ 5.5 km is ~2e-6 relative — negligible.
 *
 * Unit check: [cm^2] * [dyn cm^-2][mol^-1] / ([g mol^-1][cm s^-2])
 *           = [cm^2] * [g cm s^-2 cm^-2] / [g cm s^-2] = dimensionless.
 *
 * Cross-checks performed (v1.09, 2026-07-04):
 *   NOAA NEUBrew RayleighInBrewer.pdf Table 1 ([B99] at Brewer UV bands,
 *   sea level / 1013.25 mb / 45 deg): agreement < 0.02% at all 5 bands.
 *   (The colour-science v0.3.9 reference value 0.1004070 at 555 nm is NOT
 *    a valid anchor: that implementation passes CO2 in ppm where [B99]'s
 *    F_air mixing rule expects percent-by-volume, inflating F_air to an
 *    unphysical 1.1247 and tau by ~+7%.)
 */
double rt_rayleigh_tau_bodhaine1999_full_co2(double wavelength_nm,
                                             double pressure_hpa,
                                             double latitude_deg,
                                             double altitude_m,
                                             double co2_ppm) {
    const double A_AVOGADRO = 6.0221367e23;                /* /mol, [B99] */
    const double N_s = (A_AVOGADRO / 22.4141)
                       * (273.15 / 288.15) * 1.0e-3;       /* /cm^3 */

    double n_m1 = bodhaine1999_n_minus_1(wavelength_nm, co2_ppm);
    double n    = 1.0 + n_m1;
    double n2   = n * n;
    double F    = rt_rayleigh_king_factor_air(wavelength_nm, co2_ppm);

    double wl_cm  = wavelength_nm * 1.0e-7;
    double wl4_cm = wl_cm * wl_cm * wl_cm * wl_cm;
    double num    = (n2 - 1.0) * (n2 - 1.0);
    double den    = wl4_cm * N_s * N_s * (n2 + 2.0) * (n2 + 2.0);
    double sigma_R = 24.0 * M_PI * M_PI * M_PI * num / den * F;   /* cm^2 */

    double m_a   = 15.0556 * (co2_ppm * 1.0e-6) + 28.9595;  /* g/mol */
    /* [B99] appendix: gravity at the mass-weighted column altitude, not
     * at the surface.  z_c(z=0) = 5517.56 m.  See doc block above. */
    double z_col = 0.73737 * altitude_m + 5517.56;          /* m */
    double g_cgs = 100.0 * rt_rayleigh_gravity_wgs84(latitude_deg,
                                                     z_col);  /* gal */
    double P_dyn = pressure_hpa * 1000.0;                   /* dyn/cm^2 */

    return sigma_R * P_dyn * A_AVOGADRO / (m_a * g_cgs);
}

/* Production entry point (4-arg legacy signature).  CO2 fixed at 360 ppm =
 * the reference concentration of [B99] Table 3.  tau sensitivity to CO2 is
 * ~ +2e-5 relative per +100 ppm — irrelevant at our accuracy targets, so
 * the signature is kept CO2-free; use the _co2 variant when it matters. */
double rt_rayleigh_tau_bodhaine1999_full(double wavelength_nm,
                                         double pressure_hpa,
                                         double latitude_deg,
                                         double altitude_m) {
    return rt_rayleigh_tau_bodhaine1999_full_co2(wavelength_nm, pressure_hpa,
                                                 latitude_deg, altitude_m,
                                                 360.0);
}

/* ---------- Legacy: Hansen & Travis 1974 fit (regression A/B only) -------- */

/* The closed-form fit PREVIOUSLY (v1..v1.09) mislabeled "Bodhaine 1999
 * simplified"; the coefficients are Hansen & Travis (1974) Eq. (2.30):
 *
 *   tau_std(lambda) = 0.008569 * lambda_um^-4
 *                     * (1 + 0.0113*lambda_um^-2 + 0.00013*lambda_um^-4)
 *   tau = tau_std * (P / 1013.25 hPa) * (g_45 / g(phi,z))
 *
 * It embeds a FIXED effective King factor (no wavelength-dependent
 * depolarization).  Kept only so historical runs can be reproduced and the
 * H&T -> B99 delta can be measured; NOT selectable as production default.
 * Measured delta vs the first-principles route at 1013.25 hPa / 45 deg /
 * 360 ppm: 412 nm -0.004%, 555 nm +0.217%, 865 nm +0.335% (H&T high in
 * red/NIR; see v1.09 changelog entry). */
/* rt_rayleigh_tau_hansentravis1974_full: 2026-07-15 제거(Jae 지시).
 * 레거시 H&T-1974 근사식.  B99 대비 델타(1013.25 hPa, 45deg, 360 ppm):
 * 412 -0.004%, 555 +0.217%, 865 +0.335% -- v1.09 changelog 에 기록 보존. */

/* Backward-compat alias (v1 API name). */
double rt_rayleigh_tau_full(double wavelength_nm, double pressure_hpa,
                            double latitude_deg, double altitude_m) {
    return rt_rayleigh_tau_bodhaine1999_full(wavelength_nm, pressure_hpa,
                                             latitude_deg, altitude_m);
}

/* ---------- Public dispatchers ---------- */

/* Convenience wrapper: production model at standard reference conditions
 * (P = 1013.25 hPa, latitude 45 deg, sea level). */
double rt_rayleigh_tau(double wavelength_nm) {
    return rt_rayleigh_tau_bodhaine1999_full(wavelength_nm, 1013.25, 45.0, 0.0);
}

/* Model selector at standard reference conditions. */
double rt_rayleigh_tau_model(double wavelength_nm,
                             rt_rayleigh_model_t model) {
    return rt_rayleigh_tau_model_full(wavelength_nm, 1013.25, 45.0, 0.0, model);
}

/* Model selector, full site parameters.  The one-value enum is retained
 * only for API compatibility; any unknown value is treated as a programming
 * error and reported as NaN. */
double rt_rayleigh_tau_model_full(double wavelength_nm,
                                  double pressure_hpa,
                                  double latitude_deg,
                                  double altitude_m,
                                  rt_rayleigh_model_t model) {
    if (model != RT_RAYLEIGH_MODEL_BODHAINE_1999) return 0.0 / 0.0;
    return rt_rayleigh_tau_bodhaine1999_full(wavelength_nm, pressure_hpa,
                                             latitude_deg, altitude_m);
}
