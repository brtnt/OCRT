# Stage-2 EAP generator + phytoplankton-scattering-disabled integration

Version: `OCRT-v1.2-2026-07-26-KST-stage2-eap-phyto-scattering-disabled-integrated`

- Retains the 17-species coated-sphere P11/P12/P33 generator and catalog.
- Production constituent Chl uses `phyto_absorption_default.csv` only.
- Production constituent phytoplankton b and bb are exactly zero; detritus scattering remains.
- Explicit species selection with Chl>0 fails loudly.
- Unsupported constituent truncation flags fail loudly; fixed-bulk IOP phase truncation remains.
- Preserves current Rrs/rrs I/Q/U output, black-Fresnel naming, water m-sign fix, RAA conventions and LUT cache optimizations.
