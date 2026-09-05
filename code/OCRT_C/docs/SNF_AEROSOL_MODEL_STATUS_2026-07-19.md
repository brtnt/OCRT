# SnF aerosol model status — 2026-07-19

Canonical runtime models:

```text
T50 T80 T90 T95
C50 C70 C80 C90 C95
M50C M70C M80C M90C M95C M98C
O99
```

All use the 20-wavelength OPAC spectral grid and a 361-point, 0.5° angular grid.
Source inputs are in `inputs/aerosol_snf_inp/`; the reproducible generator is in
`tools/snf_mie_generator/`.

Qualification:

- T50/C50/M80C/O99: exact regeneration regression passed.
- All 16: structural, phase-normalization, Mueller-physicality and OCRT runtime
  smoke checks passed.
- External OSOAA/6SV qualification remains the established T50/C50/M80C subset.

M50C/M95C/M98C legacy 83-angle files are audit-only under
`inputs/deprecated/legacy83/`. Do not use them as canonical runtime inputs.

Full report: `../../docs/SNF_AEROSOL_MODEL_INTEGRATION_2026-07-19.md` at package root.
