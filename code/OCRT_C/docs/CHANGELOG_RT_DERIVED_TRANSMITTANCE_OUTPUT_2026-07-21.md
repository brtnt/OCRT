# OCRT v1.2 RT-derived transmittance output fix

- Issue: `OCRT-OUTPUT-TRANS-RT-001`
- Date: 2026-07-21
- Target: `OCRT-v1.2-2026-07-20-KST-water-mmax-mode-bound-fix` plus Mie phase-wavelength PCHIP
- Scope: output diagnostics only

## Fixed

- Removed the exported upward diffuse-transmittance approximation
  `T_diff_up=T_diff_dn`.
- Removed the exported sum
  `T_total_up=exp(-tau/mu_view)+T_diff_dn`.
- Added exact effective upward atmospheric transmission for the actual
  ocean water-leaving BRDF:

  ```text
  TOA_water_signal_I = I_TOA_total - I_atm_path
  T_total_up_view = TOA_water_signal_I / Lu0plus
  T_diff_up_view = T_total_up_view - T_dir_up_view
  ```

- Added `T_up_rt_valid`.  When the rigorous atmospheric bottom-source SOS is
  unavailable, upward transmission outputs are `NaN` and the flag is zero.
- Added `T_total_dn_hemi`, the exact sum of direct and SOS-integrated diffuse
  downward irradiance transmission.
- Renamed `rho_glint*mu_sun` and `rho_TOA*mu_sun` public columns as explicit
  flux-ratio diagnostics; they are not atmospheric transmittance.
- Added exact transmittance fields to ocean full-grid CSV output.
- Removed the production AC LUT script's reciprocity reconstruction
  `Tup(vza)=Tdn(vza)`.  Black-surface LUT rows now report upward transmission as
  undefined (`NaN`, valid flag zero).

## Unchanged physics outputs

The patch does not modify the SOS field, atmosphere-ocean coupling, water RT,
TOA assembly, `Rrs`, `rrs`, IOPs or convergence controls.  The legacy fallback
used internally when a rigorous water bottom-source pass fails remains unchanged
so that existing TOA/Rrs behavior is not altered; the fallback is simply no
longer labelled as a valid transmittance.

## Validation

- no-atmosphere, Rayleigh and Rayleigh+aerosol ocean cases: common
  non-transmittance numeric outputs exactly unchanged versus the pre-fix binary;
- exact upward and downward identities: PASS;
- black-surface undefined upward output: PASS;
- ocean full-grid exact fields and identities: PASS;
- AC LUT reciprocity removal smoke: 324 rows PASS;
- ASan/UBSan low-resolution no-atmosphere and Rayleigh: no findings;
- PCHIP regression: PASS;
- radiometry flux split: PASS;
- water branch: 14/14 PASS.

See `validation/rt_derived_transmittance_output_2026-07-21/` and
`docs/OUTPUT_RT_DERIVATION_AUDIT_2026-07-21.md`.
