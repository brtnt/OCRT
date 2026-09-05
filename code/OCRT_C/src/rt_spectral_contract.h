/* rt_spectral_contract.h
 *
 * OCRT elastic-runtime spectral support contract.
 *
 * All active input tables must explicitly cover this closed interval.  The
 * runtime may interpolate only inside the interval; it must not manufacture
 * missing component values with endpoint clamping or wavelength-specific
 * fallback branches.
 */
#ifndef OCRT_RT_SPECTRAL_CONTRACT_H
#define OCRT_RT_SPECTRAL_CONTRACT_H

#include <math.h>

#define RT_SPECTRAL_MIN_NM 330.0
#define RT_SPECTRAL_MAX_NM 1100.0
#define RT_SPECTRAL_TOL_NM 1.0e-9

static inline int rt_spectral_wavelength_supported(double wavelength_nm)
{
    return isfinite(wavelength_nm) &&
           wavelength_nm >= RT_SPECTRAL_MIN_NM - RT_SPECTRAL_TOL_NM &&
           wavelength_nm <= RT_SPECTRAL_MAX_NM + RT_SPECTRAL_TOL_NM;
}

static inline int rt_spectral_grid_covers(double first_nm, double last_nm)
{
    return isfinite(first_nm) && isfinite(last_nm) &&
           first_nm <= RT_SPECTRAL_MIN_NM + RT_SPECTRAL_TOL_NM &&
           last_nm >= RT_SPECTRAL_MAX_NM - RT_SPECTRAL_TOL_NM;
}

#endif /* OCRT_RT_SPECTRAL_CONTRACT_H */
