# Apply to the current validated OCRT hotfix source

## Target

Apply this output-only correction after the validated water-`m_max` hotfix and
the Mie phase-wavelength PCHIP correction.

Recommended order:

1. `OCRT-v1.2-2026-07-20-KST-water-mmax-mode-bound-fix`
2. Mie P11/P12/P33 wavelength PCHIP
3. `OCRT-OUTPUT-TRANS-RT-001` from this package

## Patch

```bash
cd <current-ocrt-source>/ocrt
git apply --check \
  patches/OCRT_v1.2_RT_DERIVED_TRANSMITTANCE_OUTPUT_CORE_20260721.patch
git apply \
  patches/OCRT_v1.2_RT_DERIVED_TRANSMITTANCE_OUTPUT_CORE_20260721.patch
```

The patch was generated from the PCHIP source.  Its production hunks are in
output/result definitions, final ocean-output assembly, CSV writers and scripts.
The water-`m_max` hotfix changes the external bottom-source shape contract; the
final integration session must still run `git apply --check` and the full hotfix
package regression suite.

## Required validation

```bash
./scripts/build_release_v1.2.sh build/ocrt
./scripts/regression_rt_derived_transmittance_output.sh build/ocrt
./scripts/test_phase_wavelength_pchip.sh
./scripts/smoke_radiometry_flux_split.sh build/ocrt
./scripts/smoke_cli_water_branch_v118.sh build/ocrt
./scripts/regression_ext_bottom_mode_bounds.sh build/ocrt
```

The last test is mandatory on the current hotfix source.  The validation binary
in this package is built from the 2026-07-19 source plus PCHIP and this output
fix; it is not the production hotfix binary.
