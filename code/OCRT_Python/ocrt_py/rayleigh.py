"""Rayleigh optical depth, Bodhaine et al. (1999) first-principles route.

C parity port of rt_rayleigh.c (rt_rayleigh_tau_bodhaine1999_full and its
subroutines).  Float ops follow the C order for bit-identical results.

Production entry: tau_rayleigh_bodhaine(wl_nm, pressure_hpa, lat=45, alt=0,
co2=360) -- matches rt_rayleigh_tau_model_full(..., BODHAINE_1999).
"""
import math

M_PI = math.pi


def bodhaine_n_minus_1(wl_nm, co2_ppm=360.0):
    """n - 1 of air (Bodhaine 1999 / Peck-Reeder base at 300 ppm CO2)."""
    wl_um = wl_nm / 1000.0
    sigma2 = 1.0 / (wl_um * wl_um)
    n300_m1 = (8060.51
               + 2480990.0 / (132.274 - sigma2)
               + 17455.7 / (39.32957 - sigma2)) * 1.0e-8
    c_ppv = co2_ppm * 1.0e-6
    return n300_m1 * (1.0 + 0.54 * (c_ppv - 0.0003))


def king_factor_air(wl_nm, co2_ppm=360.0):
    """King correction factor F(air)(lambda)."""
    wl_um = wl_nm / 1000.0
    wl2 = wl_um * wl_um
    wl4 = wl2 * wl2
    f_n2 = 1.034 + 3.17e-4 / wl2
    f_o2 = 1.096 + 1.385e-3 / wl2 + 1.448e-4 / wl4
    f_ar = 1.0
    f_co2 = 1.15
    x_co2_pct = co2_ppm / 1.0e4
    denom = 78.084 + 20.946 + 0.934 + x_co2_pct
    return (78.084 * f_n2 + 20.946 * f_o2 + 0.934 * f_ar
            + x_co2_pct * f_co2) / denom


def depolarization_ratio(wl_nm, co2_ppm=360.0):
    """rho = 6(F-1)/(7F+3).  For the depolarized Rayleigh phase matrix.
    (The production RT uses the fixed 6SV delta = 0.0279 instead.)"""
    F = king_factor_air(wl_nm, co2_ppm)
    return 6.0 * (F - 1.0) / (7.0 * F + 3.0)


def gravity_wgs84(lat_deg, alt_m):
    """WGS84 gravity [m/s^2] at latitude and altitude."""
    lat = lat_deg * M_PI / 180.0
    s = math.sin(lat) * math.sin(lat)
    s2 = math.sin(2.0 * lat) * math.sin(2.0 * lat)
    g0 = 9.780327 * (1.0 + 0.0053024 * s - 0.0000058 * s2)
    return g0 - 3.086e-6 * alt_m


def tau_rayleigh_bodhaine(wl_nm, pressure_hpa=1013.25, lat_deg=45.0,
                          alt_m=0.0, co2_ppm=360.0):
    """Rayleigh optical depth (Bodhaine 1999 first-principles)."""
    A_AVOGADRO = 6.0221367e23
    N_s = (A_AVOGADRO / 22.4141) * (273.15 / 288.15) * 1.0e-3

    n_m1 = bodhaine_n_minus_1(wl_nm, co2_ppm)
    n = 1.0 + n_m1
    n2 = n * n
    F = king_factor_air(wl_nm, co2_ppm)

    wl_cm = wl_nm * 1.0e-7
    wl4_cm = wl_cm * wl_cm * wl_cm * wl_cm
    num = (n2 - 1.0) * (n2 - 1.0)
    den = wl4_cm * N_s * N_s * (n2 + 2.0) * (n2 + 2.0)
    sigma_R = 24.0 * M_PI * M_PI * M_PI * num / den * F

    m_a = 15.0556 * (co2_ppm * 1.0e-6) + 28.9595
    z_col = 0.73737 * alt_m + 5517.56
    g_cgs = 100.0 * gravity_wgs84(lat_deg, z_col)
    P_dyn = pressure_hpa * 1000.0

    return sigma_R * P_dyn * A_AVOGADRO / (m_a * g_cgs)
