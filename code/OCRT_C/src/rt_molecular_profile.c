#include "rt_molecular_profile.h"

#include <math.h>
#include <stddef.h>

#include "rt_rayleigh.h"

/* U.S. Standard Atmosphere 1962 anchors used by the OCRT methodology.
 * Altitude is geometric km; pressure is Pa. */
static const double k_us62_z_km[] = {
    0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852
};
static const double k_us62_p_pa[] = {
    101325.0, 22632.1, 5474.89, 868.019,
    110.906, 66.9389, 3.95642, 0.3734
};

enum { US62_N = (int)(sizeof(k_us62_z_km) / sizeof(k_us62_z_km[0])) };

static double clamp01(double x)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    return x;
}

double rt_molecular_us62_top_km(void)
{
    return k_us62_z_km[US62_N - 1];
}

double rt_molecular_us62_pressure_pa(double altitude_km)
{
    if (!isfinite(altitude_km)) return NAN;
    if (altitude_km <= k_us62_z_km[0]) return k_us62_p_pa[0];
    if (altitude_km >= k_us62_z_km[US62_N - 1]) return k_us62_p_pa[US62_N - 1];

    size_t lo = 0;
    size_t hi = US62_N - 1;
    while (hi - lo > 1) {
        const size_t mid = lo + (hi - lo) / 2;
        if (altitude_km < k_us62_z_km[mid]) hi = mid;
        else                                  lo = mid;
    }

    const double z0 = k_us62_z_km[lo];
    const double z1 = k_us62_z_km[hi];
    const double t = (altitude_km - z0) / (z1 - z0);
    const double lp0 = log(k_us62_p_pa[lo]);
    const double lp1 = log(k_us62_p_pa[hi]);
    return exp(lp0 + t * (lp1 - lp0));
}

/* Hydrostatic molecular column is proportional to P/g.  A fixed 45-degree
 * standard latitude is used for the normalized shape; the total Rayleigh
 * optical depth is separately scaled for the actual site pressure/gravity. */
double rt_molecular_us62_fraction_above(double altitude_km)
{
    if (!isfinite(altitude_km)) return NAN;
    if (altitude_km <= 0.0) return 1.0;

    const double z_top = rt_molecular_us62_top_km();
    if (altitude_km >= z_top) {
        const double g0 = rt_rayleigh_gravity_wgs84(45.0, 0.0);
        const double gt = rt_rayleigh_gravity_wgs84(45.0, z_top * 1000.0);
        return clamp01((k_us62_p_pa[US62_N - 1] / gt) /
                       (k_us62_p_pa[0] / g0));
    }

    const double p = rt_molecular_us62_pressure_pa(altitude_km);
    const double g0 = rt_rayleigh_gravity_wgs84(45.0, 0.0);
    const double gz = rt_rayleigh_gravity_wgs84(45.0, altitude_km * 1000.0);
    return clamp01((p / gz) / (k_us62_p_pa[0] / g0));
}

double rt_molecular_us62_grid_fraction_above(double altitude_km)
{
    if (!isfinite(altitude_km)) return NAN;
    const double z_top = rt_molecular_us62_top_km();
    if (altitude_km <= 0.0) return 1.0;
    if (altitude_km >= z_top) return 0.0;

    const double f_top = rt_molecular_us62_fraction_above(z_top);
    const double f = rt_molecular_us62_fraction_above(altitude_km);
    return clamp01((f - f_top) / (1.0 - f_top));
}

double rt_molecular_us62_altitude_from_grid_fraction(double fraction_above)
{
    if (!isfinite(fraction_above)) return NAN;
    const double z_top = rt_molecular_us62_top_km();
    if (fraction_above >= 1.0) return 0.0;
    if (fraction_above <= 0.0) return z_top;

    double lo = 0.0;
    double hi = z_top;
    /* Monotone bisection.  64 iterations exceed double precision for this
     * interval and are setup-only, not part of the SOS hot loop. */
    for (int it = 0; it < 64; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double f_mid = rt_molecular_us62_grid_fraction_above(mid);
        if (f_mid > fraction_above) lo = mid;
        else                         hi = mid;
    }
    return 0.5 * (lo + hi);
}

double rt_molecular_us62_layer_fraction(double z_lo_km, double z_hi_km)
{
    if (!isfinite(z_lo_km) || !isfinite(z_hi_km)) return NAN;
    if (z_lo_km > z_hi_km) {
        const double t = z_lo_km;
        z_lo_km = z_hi_km;
        z_hi_km = t;
    }

    const double z_top = rt_molecular_us62_top_km();
    if (z_hi_km <= 0.0 || z_lo_km >= z_top) return 0.0;
    if (z_lo_km < 0.0) z_lo_km = 0.0;
    if (z_hi_km > z_top) z_hi_km = z_top;

    const double f_lo = rt_molecular_us62_grid_fraction_above(z_lo_km);
    const double f_hi = rt_molecular_us62_grid_fraction_above(z_hi_km);
    const double layer = f_lo - f_hi;
    return (layer > 0.0) ? layer : 0.0;
}
