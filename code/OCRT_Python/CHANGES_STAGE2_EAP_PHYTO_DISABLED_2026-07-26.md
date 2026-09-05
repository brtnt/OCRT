# Python Stage-2 EAP generator + phytoplankton-scattering-disabled integration

Version: `pyOCRT-v1.2-2026-07-26-stage2-eap-phyto-scattering-disabled`

- Uses the same frozen `phyto_absorption_default.csv` as the C implementation.
- Constituent phytoplankton b and bb are zero; detritus/mineral phase preparation remains.
- Skips phytoplankton phase-moment preparation in single and batch production paths.
- Explicit EAP species selection with Chl>0 fails loudly.
- Keeps the 17-species catalog and regenerated phase files for generator/future validation use.
