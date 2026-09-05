# OCRT 330-1100 nm Acceptance Test Plan

## Gate A - Package/data

- installation map row count and SHA-256
- paired C/Python byte identity
- counts: aerosol 176, AHN phase 4, AHN scalar 8, gas 6, PLOPS 1, psi_T 1; detritus 1 conditional
- scalar grids and endpoint rows
- Mie node and angle counts

## Gate B - Loader/interpolation

- exact endpoints 330 and 1100 are in-range
- 330.5 and internal nodes interpolate normally
- 329.999 and 1100.001 fail loudly
- no NaN/Inf; no silent endpoint clamp accepted as extension

## Gate C - Physicality

- P11 positive
- |P12| <= P11; |P33| <= P11
- phase normalization validator pass
- 0 <= SSA <= 1
- gas xsec nonnegative
- a_w(T) positive for 0-30 C test grid

## Gate D - C/Python parity

Matched inputs at 330, 340, 349, 350, 443, 555, 750, 865, 940, 1100 nm.

Outputs: TOA I/Q/U; Rrs and rrs I/Q/U; Ed/Eu; Kd; component IOPs; phase metrics.

## Gate E - Component RT

- pure water T sensitivity
- Chl-only 0.1/1/3
- four AHN minerals
- detritus candidate on/off if adopted
- three aerosol families
- O3, NO2, H2O and all gas

## Gate F - Coupled RT

Use pass-region geometry and edge geometry. Confirm convergence, finite results, and no unexpected runtime penalty. OSOAA comparison is restricted to common supported wavelengths and matched option contract.

## Required deliverables

- test commands and input decks
- CSV results
- I/Q/U scatter plots
- spectral plots including 330-350 boundary
- timing table
- pass/fail summary with intentional model changes separated from regressions
