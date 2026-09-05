# Third-party material in this folder — RTSOS (CC BY-NC 4.0)

`main_pp.f90`, `main_ipss_gnd.f90`, `main_ipss_toa.f90`, `main_pssonly.f90` are **modified copies of
`src/main_program_monochromatic/main.f90` from RTSOS** — Radiative Transfer model based on Successive
Orders of Scattering, Copyright © 2025 Pengwang Zhai, released under the Creative Commons
Attribution–NonCommercial 4.0 International License (https://creativecommons.org/licenses/by-nc/4.0/).
`case*/ray.pmtx` and `case*/auxiliary_directory` are copied from RTSOS `validation/benchmark/Coulson_thick/`.

Changes made for the OCRT IPSS cross-check (2026-09-05), relative to the RTSOS original:

1. `fSnowBRDF_INPUT = .FALSE.` and `fRossLiBRDF_INPUT = .FALSE.` are initialized after `NMBIMINPUT = 0.0D0`
   (the original leaves them uninitialized; with `-O2` they read as .TRUE. and the run dies at the snow-BRDF
   file read).
2. Per variant, one of the following namelist-style flags is enabled:
   `SPHERICAL_SHELL_SINGLESCATTERING_CORRECTION = .true.` (ipss_gnd, ipss_toa),
   `IPSS_VIEWANGLE_GROUND = 1` (ipss_toa only), `PSEUDO_SPHERICAL_SHELL = .true.` (pssonly).
3. Nothing else. The RTSOS core (`rtsos_rao_dg.f90` etc.) is not redistributed here; obtain it from the author.

These files are used for non-commercial research verification only. Please cite:
Zhai, P. (2025). RTSOS: Radiative Transfer model based on Successive Orders of Scattering (Version 1.0).
