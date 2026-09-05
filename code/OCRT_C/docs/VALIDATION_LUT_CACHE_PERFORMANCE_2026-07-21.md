# OCRT C LUT cache/performance validation — 2026-07-21

## Conclusion

The C code already had a LUT/full-grid implementation.  No second LUT solver
was added.  The existing path was corrected so that its intended caches are
used without redundant full solves and without near-nadir extrapolation.

## Root causes

1. `main.c` ran and discarded one complete coupled ocean case before entering
   `--output-full-grid` or `--batch-full-grid`.
2. S7 correctly bypassed near-nadir views to avoid out-of-range Gauss-node
   extrapolation, but the fallback repeated atmosphere pass-1 for every RAA.
3. Pass-2 S17 did not include the external water-leaving source in its cache
   key, allowing a preceding solve to seed a stale result.

## Implemented changes

- Early ocean full-grid/batch dispatch before ordinary single-case execution.
- S17a exact near-nadir cache: first RAA uses standard view-as-node solve;
  subsequent RAA values reconstruct the exact exported Fourier samples.
- Deep BOA field reuse and exact directional diffuse-transmittance
  reconstruction.
- Full FNV-1a fingerprint of the pass-2 external I/Q/U bottom source in S17.
- No interpolation or altered SOS tolerance/iteration/layer settings.

## Accuracy and regression results

| Test | Result |
|---|---|
| Existing single-view old vs new stdout | byte-identical |
| Existing single-view old vs new stderr | byte-identical |
| 4-cell full-grid old vs new | all CSV fields byte-identical |
| 36-cell optimized vs exact `OCRT_S17A_OFF=1` | all CSV fields byte-identical |
| Pure-water 555 nm cache on/off | byte-identical |
| Mixed-water 490 nm cache on/off | byte-identical |
| Mie P11/P12/P33 wavelength PCHIP | PASS |
| External bottom mode bounds/shape guard | PASS |
| IOP/Kd full-grid CSV regression | PASS |
| Python aerosol-object regression | 3/3 PASS |
| 2-row OpenMP batch-full-grid smoke | PASS, 0 failed |

The previous default S17 output for a 36-cell test differed in the first
near-nadir row because its cache key omitted the external source.  The corrected
output is byte-identical to the old executable with `OCRT_S17_OFF=1`, i.e. the
uncached exact reference.

## Timing results

Single CPU thread; mixed OCRT water, 490 nm, C50, AOD865=0.2, SZA=30°, wind=3.

| Configuration | Previous | Updated | Speedup |
|---|---:|---:|---:|
| n_mu_water=16, 4 cells | 7.10 s | 3.79 s | 1.87× |
| n_mu_water=16, 36 cells | 10.02 s | 4.00 s | 2.51× |
| n_mu_water=24, 4 cells | 10.82 s | 6.31 s | 1.71× |
| n_mu_water=16, 2,520 cells | >240 s; 39 rows written | 18.95 s | >12.7× lower bound |

For the 36-cell optimized run, disabling S17a increased time from 4.00 s to
7.24 s while producing the same CSV bytes.

## Remaining dominant cost

After the fix, the first cell still performs the necessary physical work:

- atmosphere pass-1: approximately 0.31 s
- water coupled SOS: approximately 2.78 s in the representative n_mu=16 case
- atmosphere pass-2: approximately 0.28 s

Later RAA cells are typically 1–4 ms.  Further acceleration now requires
optimizing the first water SOS itself or caching reusable surface/operators
across distinct cases; it is no longer a geometry-loop duplication problem.

## Python review

The 12,000-row Python production path already loads Mie models once, caches
phase/Greek data by file fingerprint and wavelength, batches rows by band/chunk,
and vectorizes IOP and RT work through NumPy/CuPy.  No corresponding production
code patch was required.  The diagnostic angular full-grid helper remains
cell-loop based but is not used by `produce_grid.py`.
