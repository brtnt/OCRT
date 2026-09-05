/* rt_eap_mie_phase.h
 *
 * Public OCRT adapter for generating one EAP phytoplankton species'
 * normalized 3-Stokes scattering-plane phase matrix.
 *
 * This interface deliberately exposes only the quantities required by the
 * current OCRT water-phase pipeline: P11, P12 and P33.  Species microphysics
 * (size distribution, core/shell refractive-index spectra, canonical water
 * refractive index and numerical Mie settings) remain implementation details.
 */
#ifndef OCRT_RT_EAP_MIE_PHASE_H
#define OCRT_RT_EAP_MIE_PHASE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RT_EAP_SPECIES_COUNT 17u

typedef uint32_t rt_eap_species_id_t;  /* Valid range: 0 .. 16. */

enum {
    RT_EAP_PHASE_OK        =  0,
    RT_EAP_PHASE_EINVAL    = -1,  /* Null pointer, invalid count or grid. */
    RT_EAP_PHASE_ESPECIES  = -2,  /* species_id is outside the 17-species catalog. */
    RT_EAP_PHASE_ERANGE    = -3,  /* Unsupported wavelength or angle range. */
    RT_EAP_PHASE_ENUMERIC  = -4   /* Mie/PSD integration or normalization failed. */
};

/* Generate normalized EAP phase-matrix elements for one species.
 *
 * Inputs
 * ------
 * species_id
 *     Stable index in the canonical 17-species EAP catalog [0, 16].
 *
 * wavelength_nm[n_wavelength]
 *     Vacuum wavelengths in nm.  Output preserves the supplied order.
 *     The OCRT EAP production range is 350–850 nm.
 *
 * theta_deg[n_theta]
 *     Caller-selected scattering-angle grid in degrees.  It must be strictly
 *     increasing, start at 0 deg and end at 180 deg.  n_theta must be >= 3.
 *
 * Outputs
 * -------
 * p11, p12, p33
 *     Caller-allocated arrays, each containing
 *         n_wavelength * n_theta
 *     doubles in wavelength-major row-major layout:
 *         out[iw * n_theta + itheta].
 *
 *     The returned scattering-plane elements follow the OCRT convention
 *
 *         P11 = (|S1|^2 + |S2|^2) / 2
 *         P12 = (|S2|^2 - |S1|^2) / 2
 *         P33 = Re(S1 * conj(S2))
 *
 *     and are normalized with one common factor so that, for every wavelength,
 *
 *         0.5 * integral_0^pi P11(theta) sin(theta) dtheta = 1.
 *
 * Contract
 * --------
 * - No memory is allocated for the caller.
 * - No file is written.
 * - No Chl, absorption, scattering coefficient, bb/b or RT moments are
 *   produced; OCRT derives those in their existing modules from this phase.
 * - The implementation shall be re-entrant and safe for concurrent calls.
 *
 * Return
 * ------
 * RT_EAP_PHASE_OK on success; a negative RT_EAP_PHASE_* code on failure.
 */
int rt_eap_mie_phase_compute(
    rt_eap_species_id_t species_id,
    const double *wavelength_nm,
    size_t n_wavelength,
    const double *theta_deg,
    size_t n_theta,
    double *p11,
    double *p12,
    double *p33
);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* OCRT_RT_EAP_MIE_PHASE_H */
