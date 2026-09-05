# Pure-water Z09 table update — 2026-07-28

Version: `OCRT-v1.2-2026-07-28-KST-pure-water-z09-table-update`

- Replaced `inputs/water_iop/water_coef_z09_1nm.txt` with the supplied
  200–2449 nm, 1 nm band-average table derived from `zo0_iop.txt`.
- The third column remains total scattering `b_w`; OCRT computes
  `bb_w = 0.5*b_w` internally.
- No RT algorithm or LUT orchestration was changed.
- Updated stale warning text: pure-water IOPs now extend through 2449 nm;
  the independent temperature-correction LUT still spans 400–900 nm.
