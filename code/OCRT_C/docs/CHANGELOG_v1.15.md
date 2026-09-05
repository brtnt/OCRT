# OCRT v1.15 changelog

- Added Ahn four-species dry-weight TSM optical-property adapter.
- Added `--tsm-species` with `red_clay` default.
- Reinterpreted the existing mineral concentration path as TSM `g m^-3` using
  species-specific mass-normalized absorption and scattering.
- Automatically connected each TSM species to its matching vector `.mie` phase.
- Removed the previous four mineral `.mie` files and installed the four supplied
  `*_AHN.mie` replacements under `inputs/tsm_ahn`.
- Added strict input conflict checks, data discovery, smoke tests and numerical
  audit tables.
- Preserved the OCRT radiative-transfer solver and legacy non-TSM results.
