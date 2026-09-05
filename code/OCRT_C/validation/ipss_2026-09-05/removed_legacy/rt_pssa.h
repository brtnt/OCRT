#ifndef OCRT_V2_RT_PSSA_H
#define OCRT_V2_RT_PSSA_H

#include "rt_types.h"

/* =============================================================================
 * rt_pssa.h — Pseudo-Spherical Approximation (PSSA) for the DIRECT solar beam
 * =============================================================================
 *
 * Reference: He, Stamnes, Bai, Li, Wang (2018), "Effects of Earth curvature
 * on atmospheric correction for ocean color remote sensing", RSE 209, 118-133
 * (the PCOART-SA model).  Also Dahlback & Stamnes (1991) average-secant
 * pseudo-spherical technique.
 *
 * SCOPE (OCRT atmospheric-only subset of He et al. 2018, Section 3):
 *   - ONLY the attenuation of the DIRECT solar beam is computed along the
 *     curved spherical-shell path:
 *       (a) the DOWNWARD beam TOA -> level k          (xi_dn[k]),
 *       (b) the flat-sea-surface REFLECTED beam
 *           TOA -> surface (local incidence alpha_k) -> level k
 *                                                       (xi_refl[k]).
 *   - ALL multiple scattering stays locally plane-parallel (the SOS
 *     iteration, layer integrals and surface BC are untouched).
 *   - View-direction transmittances exp(-tau/mu_v) stay plane-parallel
 *     (outside the pseudo-spherical scope; PSSA is accurate for high SUN
 *     zenith angles, not extreme VIEW angles — Caudill et al. 1997).
 *   - IMPORTANT MODEL-SCOPE DECISION: the paper also notes that the same
 *     geometric construction can be applied to the refracted downward solar
 *     beam in the ocean. OCRT intentionally does NOT implement that underwater
 *     correction. The --pssa option ends at the air side of the interface
 *     (0+); after Snell/Fresnel transmission the water-column direct beam is
 *     propagated by the existing plane-parallel underwater solver. Over an
 *     OCRT open-ocean column (<100 m), the omitted curvature scale is depth/R
 *     (~1.6e-5 in geometry and much smaller in the resulting irradiance
 *     correction), so this is retained as an explicit documented limitation,
 *     not a hidden or pending automatic fallback.
 *   - VALIDATION STATUS (updated 2026-07-19): atmospheric downward-beam
 *     geometry has been checked against an independent continuous-US62
 *     integration.  Reflected-beam Eq. (8) geometry, layer-local beta phase
 *     source, endpoint attenuation and first-order up/down coupling have been
 *     corrected and covered by differential and invariant tests.  Gas-only
 *     ocean-coupled direct transmittance now uses the gas-layer PSSA path.
 *     This is still NOT a claim of complete PCOART-SA paper equivalence:
 *     CDISORT/AccuRT numerical benchmark data were not reproduced and
 *     underwater PSSA is intentionally absent.  See
 *     docs/PSSA_IMPLEMENTATION_NOTES_2026-07-18.md,
 *     docs/PSSA_NUMERIC_FIXES_2026-07-19.md and
 *     docs/PSSA_SCOPE_AND_RADIOMETRY_SPLIT_2026-07-18.md.
 *
 * GEOMETRY (paper Fig. 1, Eqs. 6-8):
 *   theta0 = solar zenith angle AT THE TARGET SURFACE POINT P (= the code's
 *   atm->mu_sun).  For target level k at altitude z_k:
 *     impact parameter D_k = (R + z_k) sin(theta0)
 *     slant length of layer i:  S_i = F(z_i) - F(z_{i+1}),
 *       F(z) = sqrt((R+z)^2 - D_k^2)
 *     layer secant sec_i(k) = S_i / (z_i - z_{i+1}),
 *     xi_dn[k] = sum_{i<k} dtau_i * sec_i(k)         (paper Eq. 7)
 *   Cancellation-free evaluation used in rt_pssa.c:
 *     F(z)^2 = ((R+z) mu0)^2 + (z - z_k)(2R + z + z_k) sin^2(theta0)
 *     sec_i(k) = (2R + z_i + z_{i+1}) / (F(z_i) + F(z_{i+1}))
 *   (exact identities; no large-number subtraction at grazing theta0).
 *
 *   Reflected beam (paper Eq. 8; the "2" is dropped by the journal's PDF
 *   text extraction — re-derived from Fig. 1b spherical triangle):
 *     sin(alpha)/sin(beta) = (R + H_k)/R,   2*alpha = beta + theta0
 *   alpha_k is solved by bisection of
 *     f(a) = (R+z_k) sin(2a - theta0) - R sin(a) = 0   on (theta0/2, theta0];
 *   f is monotone increasing there, f(theta0) = z_k sin(theta0) >= 0.
 *   xi_refl[k] = xi_dn(surface; incidence alpha_k)
 *              + xi_up(surface -> level k; impact parameter R sin(alpha_k)).
 *
 * WHAT rt_pssa_apply() DOES to atm:
 *   - fills atm->pssa_xi_dn / pssa_xi_refl / pssa_alpha / pssa_beta (lazy-allocated,
 *     freed by rt_atm_free()),
 *   - OVERWRITES atm->ch[k] = 0.5 * exp(-xi_dn[k])  (site A: the primary-
 *     scattering source attenuation; replaces 0.5*exp(-h[k]/mu_sun)),
 *   - sets atm->pssa_active = 1.
 *   The remaining consumers (sites B..H) branch on pssa_active inside
 *   rt_solver.c; with pssa_active == 0 every expression is bit-identical
 *   to the plane-parallel baseline.
 *
 * PRECONDITIONS: atm fully built — h[] must be FINAL (Rayleigh + aerosol +
 * gas absorption, i.e. call AFTER rt_atm_apply_gas_absorption) and
 * z_km_level[] populated by the builders (both rt_atm_build_rayleigh and
 * rt_atm_build_aerosol_rayleigh fill it).
 *
 * COST: O(nt^2) sqrt-level work + O(nt * 100) bisection steps, once per
 * (case, wavelength).  nt = 40 default -> well under 0.1 ms, i.e. invisible
 * next to the seconds-scale SOS solve (speed requirement).
 *
 * HIGH-SZA LAYER-RESOLUTION ADVISORY:
 *   PSSA is sensitive to the altitude span represented by each layer near the
 *   horizon.  rt_pssa_apply() emits a warning (no numerical auto-change) for:
 *     75 <= SZA < 80 deg and n_layers < 100,
 *     80 <= SZA < 84 deg and n_layers < 200,
 *     SZA >= 84 deg       and n_layers < 400.
 *   Near 85 deg, the current 412-nm US62 audit indicates >=200 layers for
 *   conservative total-TOA convergence, and 400-800 layers when the curvature
 *   correction magnitude itself must be sub-percent.  The warning is emitted
 *   outside the SOS hot loop and does not change numerical results.
 *
 * VALIDATION SUMMARY:
 *   - black and rough-surface unaffected controls remain byte-exact.
 *   - atmospheric downward direct-beam shell geometry passed independent
 *     US62 checks over the documented SZA 70-85 deg smoke range.
 *   - reflected-beam Eq. (8), surface and attenuation invariants pass at
 *     machine precision; local-beta basis construction is setup-only.
 *   - flat-surface reflected-beam and gas-only ocean fixes are quantified in
 *     validation/pssa_numeric_fixes_2026-07-19/.
 *   - no underwater spherical correction is implemented.
 *   - full paper-equivalence remains unclaimed because the published
 *     CDISORT/AccuRT numerical fields are not bundled as executable fixtures.
 */

/* Earth radius (km).  He et al. (2018) and the Adams & Kattawar (1978)
 * benchmarks both use 6371 km.  Fixed constant by design (no CLI option). */
#define RT_PSSA_EARTH_RADIUS_KM 6371.0

/* Apply the pseudo-spherical correction to a fully built atmosphere.
 * Returns 0 on success, -1 on bad args / allocation failure.
 * On failure atm is left with pssa_active == 0 and ch[] untouched. */
int rt_pssa_apply(rt_atm_t *atm);

#endif /* OCRT_V2_RT_PSSA_H */
