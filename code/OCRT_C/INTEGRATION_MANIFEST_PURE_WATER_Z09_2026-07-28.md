# Integration manifest — pure-water Z09 table update (2026-07-28)

## Canonical package set

- OCRT C: `OCRT-v1.2-2026-07-28-KST-pure-water-z09-table-update`
- OCRT Python: `pyOCRT-v1.2-2026-07-28-pure-water-z09-table-update`
- OSOAA: `OSOAA_OCRT_reference_tool_2026-07-28-pure-water-z09`

## Runtime data replacement

| Package | Runtime file |
|---|---|
| OCRT C | `inputs/water_iop/water_coef_z09_1nm.txt` |
| OCRT Python | `data/water_iop/water_coef_z09_1nm.txt` |
| OSOAA | `fic/OSOAA_SEA_MOL_COEFFS_JUNE_2013.txt` |

All three use the identical 2,250-row numeric array over 200–2449 nm.

## Non-data edits

- Version strings and README descriptions
- stale 900 nm pure-water warning text
- data inventory range/row count/hash
- active regression anchors that necessarily changed with the data
- OSOAA active smoke reference and package manifest
- reproducibility tests and source-handoff copy

No RT solver equation, phase kernel, interface operator, or native coupled-LUT algorithm was changed.
