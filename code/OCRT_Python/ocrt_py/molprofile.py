"""rt_molecular_profile.c transliteration: U.S. Standard Atmosphere 1962
cumulative molecular column profile for Rayleigh layering (OCRT v1.2).

Anchors are geometric km / Pa.  fraction_above(z) = (P(z)/g(z)) / (P0/g0)
with the WGS84 international-gravity-formula family at 45 deg latitude
(rt_rayleigh_gravity_wgs84).  grid_fraction_above renormalizes so the model
top (84.852 km) maps to exactly zero; altitude_from_grid_fraction inverts it
by 64-step monotone bisection (setup-only, not in the SOS hot loop).
"""
from math import log, exp, isfinite

from .rayleigh import gravity_wgs84
from .backend import xp


def _clamp01_x(x):
    """clamp to [0,1], array-safe."""
    return xp.minimum(1.0, xp.maximum(0.0, x))

_US62_Z_KM = (0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852)
_US62_P_PA = (101325.0, 22632.1, 5474.89, 868.019,
              110.906, 66.9389, 3.95642, 0.3734)
_US62_N = len(_US62_Z_KM)


def _clamp01(x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    return x


def us62_top_km():
    return _US62_Z_KM[_US62_N - 1]


def us62_pressure_pa(altitude_km):
    """Log-linear pressure interpolation between US62 anchors.  Array-safe:
    scalar in -> float out, array in -> array out."""
    z = xp.asarray(altitude_km, dtype=float)
    scalar = (z.ndim == 0)
    zk = xp.asarray(_US62_Z_KM, dtype=float)
    lpk = xp.log(xp.asarray(_US62_P_PA, dtype=float))
    idx = xp.clip(xp.searchsorted(zk, z, side='right') - 1, 0, _US62_N - 2)
    z0 = zk[idx]; z1 = zk[idx + 1]; lp0 = lpk[idx]; lp1 = lpk[idx + 1]
    t = (z - z0) / (z1 - z0)
    p = xp.exp(lp0 + t * (lp1 - lp0))
    p = xp.where(z <= zk[0], _US62_P_PA[0], p)
    p = xp.where(z >= zk[-1], _US62_P_PA[-1], p)
    return float(p) if scalar else p


def us62_fraction_above(altitude_km):
    """Hydrostatic molecular column fraction above z: (P/g)/(P0/g0), 45 deg.
    Array-safe."""
    z = xp.asarray(altitude_km, dtype=float)
    scalar = (z.ndim == 0)
    z_top = us62_top_km()
    g0 = gravity_wgs84(45.0, 0.0)
    P0 = _US62_P_PA[0]
    p = us62_pressure_pa(z)
    gz = gravity_wgs84(45.0, z * 1000.0)
    frac = _clamp01_x((p / gz) / (P0 / g0))
    gt = gravity_wgs84(45.0, z_top * 1000.0)
    frac_top = _clamp01_x((_US62_P_PA[_US62_N - 1] / gt) / (P0 / g0))
    frac = xp.where(z >= z_top, frac_top, frac)
    frac = xp.where(z <= 0.0, 1.0, frac)
    return float(frac) if scalar else frac


def us62_grid_fraction_above(altitude_km):
    """fraction_above renormalized so the US62 top maps to exactly 0.
    Array-safe."""
    z = xp.asarray(altitude_km, dtype=float)
    scalar = (z.ndim == 0)
    z_top = us62_top_km()
    f_top = us62_fraction_above(z_top)
    f = us62_fraction_above(z)
    frac = _clamp01_x((f - f_top) / (1.0 - f_top))
    frac = xp.where(z <= 0.0, 1.0, frac)
    frac = xp.where(z >= z_top, 0.0, frac)
    return float(frac) if scalar else frac


def us62_altitude_from_grid_fraction(fraction_above):
    """Invert us62_grid_fraction_above by 64-step monotone bisection.
    Array-safe (all elements share the 64 iterations, branch-free)."""
    fa = xp.asarray(fraction_above, dtype=float)
    scalar = (fa.ndim == 0)
    z_top = us62_top_km()
    lo = xp.zeros_like(fa)
    hi = xp.full_like(fa, z_top)
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        f_mid = us62_grid_fraction_above(mid)
        gt = f_mid > fa
        lo = xp.where(gt, mid, lo)
        hi = xp.where(gt, hi, mid)
    z = 0.5 * (lo + hi)
    z = xp.where(fa >= 1.0, 0.0, z)
    z = xp.where(fa <= 0.0, z_top, z)
    return float(z) if scalar else z


def us62_layer_fraction(z_lo_km, z_hi_km):
    """Molecular column fraction inside [z_lo, z_hi] (grid-normalized)."""
    if not isfinite(z_lo_km) or not isfinite(z_hi_km):
        return float('nan')
    if z_lo_km > z_hi_km:
        z_lo_km, z_hi_km = z_hi_km, z_lo_km
    z_top = us62_top_km()
    if z_hi_km <= 0.0 or z_lo_km >= z_top:
        return 0.0
    if z_lo_km < 0.0:
        z_lo_km = 0.0
    if z_hi_km > z_top:
        z_hi_km = z_top
    f_lo = us62_grid_fraction_above(z_lo_km)
    f_hi = us62_grid_fraction_above(z_hi_km)
    layer = f_lo - f_hi
    return layer if layer > 0.0 else 0.0
