# Ahn four-species TSM optical data

This directory is the canonical OCRT mineral/TSM dataset introduced in v1.15.

Species and vector-phase files:

- `red_clay` → `Red_clay_AHN.mie` (default)
- `brown_earth` → `Brown_earth_AHN.mie`
- `yellow_clay` → `Yellow_clay_AHN.mie`
- `calcareous_sand` → `Calcareous_sand_AHN.mie`

The four `.mie` files are the files supplied for this integration. The
`astarmin_*` and `bstarmin_*` audit tables are derived without resampling from
their bulk spectral blocks:

- `a* = Extinct_Co - Scatter_Co`
- `b* = Scatter_Co`

OCRT uses these mass-specific coefficients with dry-weight TSM concentration
in g m^-3. The backscatter ratio is integrated from the matching P11 phase.
Set `OCRT_TSM_DIR` to override this directory.
