# Validation reference policy after the 2026-07-28 pure-water IOP replacement

The runtime pure-water table changed intentionally. Therefore, validation artifacts generated before
2026-07-28 are historical evidence, not active numeric anchors, whenever their result depends on
pure-water absorption or scattering.

Active current-data gates regenerated or verified in this package include:

- `scripts/test_pure_water_z09_table.sh`
- `scripts/regression_water_intrefl_msign.sh`
- `scripts/regression_iop_kd_fullgrid_csv.sh`
- `validation/stage2_coupling_clamp_fix_2026-07-25/reference/fullgrid_n{48,64,96}.csv`
- native coupled-LUT vs cell-replay equality tests

The previous coupling references are retained under
`validation/stage2_coupling_clamp_fix_2026-07-25/reference_pre_pure_water_z09_20260728/`.
Other dated validation folders remain for provenance and must not be interpreted as current pure-water
goldens unless explicitly regenerated after 2026-07-28.
