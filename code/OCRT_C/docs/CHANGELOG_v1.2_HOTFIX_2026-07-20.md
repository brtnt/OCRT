# OCRT v1.2 hotfix — 2026-07-20

## Fixed

- Fixed heap out-of-bounds reads in the coupled-ocean C3d external atmospheric
  bottom-source pass when `water_m_max < atmospheric_m_max`.
- Allocated water-leaving source storage to the full atmospheric Fourier bound
  and zero-padded modes not produced by the water solver.
- Added explicit external bottom-source shape metadata and error code
  `RT_SOLVER_ERR_BOTTOM_SOURCE_SHAPE`.
- Added mode bounds before external source slicing and cleared stale per-mode
  pointers.
- Added S7b Gauss-ring source projection, shape-aware fingerprinting, and safe
  exact-row fallback.

## Regression

- Stable `water_m_max >= 16` paths remain byte-identical.
- Low caps `2/4/8` now return finite converged results.
- ASan/UBSan and dedicated zero-padding tests pass.
