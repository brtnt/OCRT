# OCRT Rrs/rrs I/Q/U output completion and validation

Date: 2026-07-25

## 1. Scope and initial finding

The report distinguishes the two conventional reflectances:

- `Rrs(0+)_[I,Q,U] = [Lu,Qu,Uu](0+) / Ed(0+)` in air.
- `rrs(0-)_[I,Q,U] = [Lu,Qu,Uu](0-) / Ed(0-)` just below the sea surface.

The original situation was not simply “rrs was scalar everywhere.”

- C single-case text output already derived and printed `rrs0minus_I/Q/U`.
- Python single/coupled driver already returned `rrs0minus_I/Q/U`.
- The C public result structures contained only the scalar `rrs_I` field, and C ocean full-grid CSV emitted only `rrs_I`.
- The Python 1,000-case × 12-band production CSV retained the legacy scalar `Rrs` but did not emit vector `Rrs_I/Q/U` or `rrs_I/Q/U`.

Therefore the fix completes the public result schema and production/full-grid output. It does not add another radiative-transfer calculation.

## 2. C changes

Version: `OCRT-v1.2-2026-07-25-KST-rrs-iqu-output`

- Added explicit fields to `rt_result_t` and `rt_water_result_t`:
  - `r_rs_0minus_Q`, `r_rs_0minus_U`
  - `R_rs_0plus_Q`, `R_rs_0plus_U`
- Filled all six vector reflectances directly from the Stokes boundary radiances already computed by the coupled water/atmosphere solver.
- Added `rrs_Q` and `rrs_U` to the ocean full-grid CSV.
- Kept the existing single-case names and numerical values.
- Added `scripts/regression_rrs_iqu_output.sh`.

No SOS order, Fourier mode, quadrature node, atmosphere pass, water pass, Mie calculation, or coupling pass was added.

## 3. Python changes

Version: `pyOCRT-v1.2-2026-07-25-rrs-iqu-output`

- The existing batch water solve now exposes the already retained per-mode Stokes boundary fields as:
  - `Lu/Qu/Uu_0minus`
  - `Lu/Qu/Uu_0plus`
- `batch_driver.solve_r1_grid()` returns:
  - `Rrs_I`, `Rrs_Q`, `Rrs_U`
  - `rrs_I`, `rrs_Q`, `rrs_U`
- `produce_grid.py` emits these six columns in the 12,000-row production CSV.
- Legacy `Rrs` is retained as an exact backward-compatible alias of `Rrs_I`.
- Resume mode rejects an old-schema CSV instead of appending incompatible rows.
- The specialized atmosphere-free Tier-0 scalar LUT keeps a schema-complete `rrs_Q/U=0`; this route is scalar by design and is not the full-vector production path.
- Added an automatic call-count regression proving that the solve graph remains two atmosphere passes plus one water solve.

## 4. Numerical verification

Representative C full-grid cell (`490 nm`, `SZA=30`, `VZA=30`, `RAA=90`, `C50`, `AOD865=0.2`, mixed water):

- `rrs_I = 7.128574300484e-02`
- `rrs_Q = 5.890505446490e-05`
- `rrs_U = 1.794969309081e-03`

For every tested grid row:

- `rrs_I == Lu0minus_I / Ed0minus_water`
- `rrs_Q == Lu0minus_Q / Ed0minus_water`
- `rrs_U == Lu0minus_U / Ed0minus_water`

The old and new C CSVs have 58 common columns. Only three last-digit differences occurred after recompilation:

- maximum absolute difference: `1.00031e-13`
- maximum relative difference: `4.37646e-13`
- new columns only: `rrs_Q`, `rrs_U`

Representative Python full-vector batch output:

- legacy `Rrs = 1.1690135220725314e-02`
- `Rrs_I = 1.1690135220725314e-02` (exact equality)
- `Rrs_Q = 1.1033208286738593e-04`
- `Rrs_U = 1.1839328756155765e-03`
- `rrs_I = 2.4414108682258948e-02`
- `rrs_Q = -1.3809651892006778e-05`
- `rrs_U = 2.4726787313581683e-03`

The Python regression also verifies that Stokes reflectance ratios reproduce the already-solved boundary-radiance ratios without another solver call.

## 5. Runtime verification

### C full-grid, 36 cells, single thread, five repetitions

- baseline median: `2.48 s`
- modified median: `2.48 s`
- measured change: `0.00%`

The C addition consists of existing-field assignments, four scalar divisions per cell, and two extra CSV values.

### Python full-vector R1

Realistic one-case warm medians:

- baseline: `2.3086 s`
- modified: `2.2345 s`

The difference is timing noise; there is no measurable slowdown. A deliberately tiny low-resolution test showed `+0.0091 s` because the fixed output reconstruction is a larger fraction of an artificially short solve, but the RT call graph is unchanged and the realistic calculation showed no increase.

## 6. Regression results

C:

- Rrs/rrs I/Q/U output algebra: PASS
- IOP/Kd full-grid CSV: PASS (`60` columns)
- water internal-reflection msign regression: PASS
- RAA convention/direct-glint regression: PASS

Python:

- compileall: PASS
- pytest: `7 passed`
- no-extra-RT call-count test: PASS (`atmosphere=2`, `water=1`)
- legacy `Rrs == Rrs_I`: PASS
- production resume schema guard: PASS
- Tier-0 output schema/format check: PASS

## 7. Conclusion

The full-vector OCRT calculation already contained Q/U. The completed implementation exposes those existing values consistently in C result structures, C full-grid CSV, Python batch results, and the Python 12,000-row production CSV. No additional RT solve was introduced, and no meaningful runtime increase was measured.
