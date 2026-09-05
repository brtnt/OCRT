# OCRT C LUT cache/performance fix — 2026-07-21

Version: `OCRT-v1.2-2026-07-21-KST-lut-cache-performance-fix`

## Scope

This change improves the existing coupled-ocean full-grid/LUT path.  It does
not add a second LUT solver and does not alter the RT equations.

## Fixes

1. **Removed discarded pre-grid solve.**
   `main.c` previously completed one full coupled single-view solve and then
   discarded it before entering `--output-full-grid` or `--batch-full-grid`.
   Grid dispatch now occurs immediately after absorption setup.

2. **Exact near-nadir RAA folding (S17a).**
   S7 correctly refuses Gauss-node extrapolation near nadir.  The fallback
   previously repeated atmosphere pass-1 for every RAA.  The first RAA now
   exports exact view-as-node Fourier samples and the BOA field; subsequent
   RAA values use the same Fourier reconstruction as the solver.  No VZA
   interpolation/extrapolation is introduced.

3. **Pass-2 S17 cache-key correctness.**
   The legacy S17 key omitted the external water-leaving bottom-source field.
   A preceding solve could therefore seed a stale pass-2 result.  The complete
   shaped I/Q/U source, mode bound and stride are now hashed into the key.

## Accuracy gates

- Single-view old/new stdout and stderr: byte-identical.
- 4-cell full-grid old/new: byte-identical.
- 36-cell optimized vs exact `OCRT_S17A_OFF=1`: byte-identical for all CSV
  fields.
- Pure-water 555 nm and mixed-water 490 nm optimized/reference full grids:
  byte-identical.
- Mie P11/P12/P33 PCHIP, external-bottom bounds, IOP/Kd full-grid regressions:
  PASS.

## Timing, one CPU thread

Representative mixed-water case at 490 nm, C50 AOD865=0.2, n_mu_water=16:

| Grid | Previous | Updated | Speedup |
|---|---:|---:|---:|
| 4 cells | 7.10 s | 3.79 s | 1.87x |
| 36 cells | 10.02 s | 4.00 s | 2.51x |
| 2,520 cells (2.5° VZA × 5° RAA) | >240 s; only 39 rows written | 18.95 s | >12.7x lower bound |

At n_mu_water=24, 4 cells improved from 10.82 s to 6.31 s (1.71x).
