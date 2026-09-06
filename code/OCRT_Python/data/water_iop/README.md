# pyOCRT water IOP data

The 2026-07-26 EAP update replaces the obsolete three phytoplankton `.mie`
files with 17 regenerated representative files under `eap/`.

Compatibility aliases:

- `pico` -> `EAP_15_Synechococcus_D1p2.mie`
- `nano` -> `EAP_12_Hapto_Prymnesiaceae_D4.mie`
- `micro` -> `EAP_02_Diatoms_centric_D6.mie`

All 17 canonical `eap_*` names are exposed by the Python CLI. Original source
data and the reproducible layered-sphere generator are included in the package.
