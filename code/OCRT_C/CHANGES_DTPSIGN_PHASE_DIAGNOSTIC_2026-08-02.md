# OCRT v1.2 — DTPSIGN coupled-basis fix and single-component phase diagnostic

Date: 2026-08-02

## Physical bug fix

At the atmosphere-BOA → underwater diffuse-top hand-off, odd Fourier modes are
multiplied by `(-1)^m`.  The atmosphere stores coefficients in the legacy
`phi+pi` basis while the underwater source uses `cos/sin(m*phi)`.  The
conversion is applied to the copied hand-off arrays, not inside the diffuse-top
injector, preserving the source/operator closure identity.

Affected: coupled atmosphere–ocean runs, odd azimuthal modes only.
Unchanged: atmosphere-only, no-atmosphere ocean, and all even modes.

## Diagnostic improvement

`OCRT_DEBUG=1 OCRT_DUMP_IOP=1` now emits `OCRT_PHASE_COMPONENT` and
`OCRT_PHASE_MIX` for a single native constituent such as TSM-only.  This is a
debug-only output change and adds no production-path work.

## Deferred items

* The reported dedicated-pure-water vs zero-constituent mismatch was not
  reproduced in this canonical tree (reflectance differences were about
  0.001–0.002%, not 0.3–0.4 percentage points), so the default path is unchanged.
* Constituent-model phase truncation remains fail-loud.  The existing validated
  fixed-bulk `--iop-mie-truncation` path remains the OSOAA parity route; the
  constituent extension is optional and was not adopted without a fully
  validated transport/source bookkeeping implementation.
