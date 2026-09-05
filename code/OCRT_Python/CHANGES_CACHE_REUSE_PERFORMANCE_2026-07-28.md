# pyOCRT cache-reuse performance update

Version: `pyOCRT-v1.2-2026-07-28-cache-reuse-performance`

- Shares an exact rough-Fresnel Fourier matrix cache between atmosphere pass 1
  and pass 2 of the native coupled angular LUT.
- The cache key includes the complete output/input cosine arrays and every
  surface option that affects the matrix.
- Stored arrays and returned hits are copied, so solver code does not own or
  mutate the cache.
- `OCRT_PY_AIR_SURFACE_CACHE_OFF=1` disables the cache for regression tests.
- Python layer transmissions were already vectorized/precomputed; no duplicate
  `exp(-delta_tau/|mu|)` loop corresponding to the C optimization existed.
