# pyOCRT v1.2 — DTPSIGN coupled-basis fix and phase diagnostic

Date: 2026-08-02

## Physical correction

The atmosphere-BOA Fourier field is stored in the legacy `phi+pi` basis,
while the underwater diffuse-top source uses `cos/sin(m*phi)`.  At the copied
hand-off arrays, odd Fourier modes are therefore multiplied by `(-1)^m`.

Applied to:

- `ocrt_py.driver.solve_case_ocean_coupled`
- `ocrt_py.driver.solve_case_ocean_coupled_lut`
- the batched R1 atmosphere-to-water hand-off in `ocrt_py.batch_driver`

The conversion is intentionally not placed inside the diffuse-top source
injector, preserving source/operator closure.  It is an in-place sign change
performed once per case and adds no radiative-transfer iteration.

## Diagnostic improvement

With both `OCRT_DEBUG=1` and `OCRT_DUMP_IOP=1`, a single active constituent
(e.g. TSM-only) now emits structured `OCRT_PHASE_COMPONENT` and
`OCRT_PHASE_MIX` lines.  Normal production execution performs no additional
phase evaluation or logging.

## Deferred items

- The reported 0.3–0.4 percentage-point dedicated-pure-water discrepancy was
  not reproduced in this canonical tree.  The measured reflectance difference
  was about 0.001–0.002%, so the default pure-water path was not changed.
- Constituent-model phase truncation remains fail-loud.  The validated
  fixed-bulk truncation route remains the OSOAA-parity path.
