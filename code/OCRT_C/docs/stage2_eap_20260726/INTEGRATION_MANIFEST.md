# Integration manifest

## C package

- 17 EAP representative `.mie` files under `inputs/water_iop/eap/`.
- Legacy EAP `.mie` files absent.
- Frozen default Chl absorption:
  `inputs/water_iop/phyto_absorption_default.csv`.
- Production Chl scattering disabled; detritus scattering retained.
- Recent C fixes and output/LUT features retained.

## Python package

- Same 17 EAP representative files under `data/water_iop/eap/`.
- Same frozen default Chl absorption CSV and SHA-256.
- Single/full-grid/production CLI Stage-2 option gate.
- Batch phase preparation skips zero-scattering phytoplankton.

## OSOAA package

- Public core source unchanged.
- Added `stage2_harness/`, reference CSV copies, source hashes and validated
  Linux executable.
- No claim of exact modified-OSOAA raw reproduction.
