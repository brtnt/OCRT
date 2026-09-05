# pyOCRT rrs I/Q/U output completion — 2026-07-22

Version: `pyOCRT-v1.2-2026-07-25-rrs-iqu-output`

- Single-case and coupled full-grid paths already exposed Rrs/rrs I/Q/U.
- The 12,000-row production batch previously exported only Rrs I.
- Added vectorized reconstruction/export of Rrs(0+) and rrs(0-) I/Q/U from
  per-mode fields already retained by the water SOS.
- No additional RT solve, scattering order, kernel build, or quadrature pass.
