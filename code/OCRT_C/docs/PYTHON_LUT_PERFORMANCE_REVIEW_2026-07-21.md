# Python LUT performance review — 2026-07-21

## Scope reviewed

The production path requested by the user is `produce_grid.py`: 1,000 cases ×
12 bands = 12,000 outputs.  This is distinct from the diagnostic
`ocrt_fullgrid.py` / `ocrt_py.lut.run_ocean_aerosol_full_grid` angular-grid
helper.

## Production path findings

No C-style fix is required in the Python production path.

- `_mie_cache()` loads every unique aerosol `.mie` file once per production
  process.
- `_AER_PHASE_CACHE` is keyed by the Mie file/stat fingerprint, wavelength and
  numerical phase options.  PCHIP wavelength interpolation, log-linear
  truncation and Greek coefficients are reused independently of AOD.
- `produce_grid.py` processes one band at a time and groups pending rows into
  chunks.
- `solve_r1_grid`, `solve_r2_grid` and `solve_r3_grid` use grouped NumPy/CuPy
  batch solvers rather than launching one complete process per case.
- In R1, constituent IOPs are vectorized by band and water component moments
  are computed once per component/band and reused.
- There is no discarded single-case solve before the production batch.
- Resume/skip state is keyed by `(case_id, band_nm)`.

The existing aerosol-object regression suite passes: 3 tests in 0.41 s.

## Diagnostic full-grid helper

`ocrt_py.lut.run_ocean_aerosol_full_grid()` prepares one immutable
`AerosolRuntime`, but still calls the complete coupled solver for every
VZA/RAA cell.  It is a correctness/diagnostic helper and is not used by the
12,000-row production path.  Rewriting it to expose per-Fourier fields and
perform C-style RAA reconstruction would be a separate API-level optimization;
it was not applied because it would add risk without accelerating the requested
production workflow.

## Decision

Python production code modification: **not required**.

Future optional work: replace the diagnostic angular-grid helper with a batched
backend call if large VZA–RAA map generation becomes a production requirement.
