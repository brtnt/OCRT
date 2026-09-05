"""OCRT elastic-runtime spectral support contract.

Active runtime tables must explicitly cover the closed interval 330--1100 nm.
Interpolation is permitted only inside a table's native range.  Out-of-range
endpoint clamping and component-specific wavelength fallbacks are prohibited.
"""
from __future__ import annotations

import math
import numpy as np

SPECTRAL_MIN_NM = 330.0
SPECTRAL_MAX_NM = 1100.0
SPECTRAL_TOL_NM = 1.0e-9


def require_wavelength(wavelength_nm: float, *, context: str = "OCRT") -> float:
    wl = float(wavelength_nm)
    if (not math.isfinite(wl) or
            wl < SPECTRAL_MIN_NM - SPECTRAL_TOL_NM or
            wl > SPECTRAL_MAX_NM + SPECTRAL_TOL_NM):
        raise ValueError(
            f"{context}: wavelength {wl!r} nm is outside the supported "
            f"{SPECTRAL_MIN_NM:.0f}-{SPECTRAL_MAX_NM:.0f} nm interval")
    return min(max(wl, SPECTRAL_MIN_NM), SPECTRAL_MAX_NM)


def require_table_coverage(wavelengths, *, context: str) -> np.ndarray:
    wl = np.asarray(wavelengths, dtype=float)
    if wl.ndim != 1 or wl.size < 2 or not np.all(np.isfinite(wl)):
        raise ValueError(f"{context}: invalid wavelength axis")
    if not np.all(np.diff(wl) > 0.0):
        raise ValueError(f"{context}: wavelength axis is not strictly increasing")
    if (wl[0] > SPECTRAL_MIN_NM + SPECTRAL_TOL_NM or
            wl[-1] < SPECTRAL_MAX_NM - SPECTRAL_TOL_NM):
        raise ValueError(
            f"{context}: table range {wl[0]:.12g}-{wl[-1]:.12g} nm does not "
            f"cover {SPECTRAL_MIN_NM:.0f}-{SPECTRAL_MAX_NM:.0f} nm")
    return wl


def require_query_in_table(wavelength, wavelengths, *, context: str, unit: str = "nm") -> float:
    """Require a finite query to lie inside a table axis in the same units."""
    wl = float(wavelength)
    if not math.isfinite(wl):
        raise ValueError(f"{context}: non-finite wavelength query")
    axis = np.asarray(wavelengths, dtype=float)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError(f"{context}: invalid wavelength axis")
    if wl < axis[0] - SPECTRAL_TOL_NM or wl > axis[-1] + SPECTRAL_TOL_NM:
        raise ValueError(
            f"{context}: wavelength {wl:.12g} {unit} lies outside table range "
            f"{axis[0]:.12g}-{axis[-1]:.12g} {unit}")
    return wl
