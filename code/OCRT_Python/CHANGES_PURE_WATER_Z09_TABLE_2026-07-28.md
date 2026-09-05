# Pure-water Z09 table update — 2026-07-28

Version: `pyOCRT-v1.2-2026-07-28-pure-water-z09-table-update`

- Replaced `data/water_iop/water_coef_z09_1nm.txt` with the same
  200–2449 nm, 1 nm band-average table used by OCRT C and OSOAA.
- The third column is total scattering `b_w`; Python returns
  `bb_w = 0.5*b_w`.
- No solver or native coupled-LUT routine was changed.
