# OCRT immutable-cache performance update (2026-07-28)

This update removes two classes of repeated immutable work from the C vector
SOS implementation.  It does not alter radiative-transfer equations,
quadrature, convergence criteria, layer selection, the native coupled angular
LUT orchestration, or any public output definition.

## 1. Air-side rough-Fresnel Fourier cache

`surface_coxmunk_fourier_kernel()` previously repeated the complete azimuthal
integration of the air-side rough-Fresnel Mueller operator whenever atmosphere
pass 1/pass 2 requested the same angular slab.  The result depends only on:

- output and input cosine arrays,
- Fourier mode and azimuth quadrature size,
- wind speed and slope-variance selector,
- water refractive index and Q convention.

The implementation now keeps a bounded 128-slot exact cache per worker thread.
Every scalar and both cosine arrays are checked byte-for-byte before a hit is
served.  The cached slab is copied to the caller, so later caller-side mutation
cannot alter the stored value.  The cache can be disabled for authoritative
regression with:

```text
OCRT_AIR_SURFACE_CACHE_OFF=1
```

## 2. Layer transmission stencil

For a finalized medium, the Beer-Lambert term

```text
T(direction, layer) = exp(-delta_tau[layer] / abs(mu[direction]))
```

is invariant across Stokes I/Q/U and all scattering orders.  It was previously
recalculated in every formal-solution call.  `rt_atm_t` now owns the complete
`direction x layer` transmission stencil plus exact snapshots of `h[]` and
`rm[]`.  A cache hit is accepted only when both snapshots remain byte-identical;
view-node or absorption-grid changes therefore rebuild safely.  Allocation
failure falls back to the historical per-call `exp()` path.

Disable for regression with:

```text
OCRT_TRANSPORT_CACHE_OFF=1
```

At the 2026-07-28 high-turbidity production grid (`n_mu_water≈101` including
output nodes, 1,488 layers), the stencil requires roughly 2.4 MB for `T`, plus
about 14 kB for the exact keys.

## 3. Deliberately rejected changes

The following prototypes were measured but are not part of the canonical
update:

- ping-pong pointer swaps replacing historical field copies: small-grid tests
  passed, but long high-layer runs showed a severe post-order-10 performance
  regression;
- skipping full NaN/finite checks in the release kernel: short tests passed but
  long-order behavior was not robust;
- homogeneous-water phase-mix hoisting: about 6% faster in the reduced case,
  but changed a few 12-digit CSV values by up to `2.3e-13`.  The historical
  arithmetic is retained to preserve byte-identical validation data.

## 4. LUT contract

The native coupled atmosphere-ocean LUT remains solve-once:

- atmosphere pass 1: once per case/band,
- in-water cold SOS: once per case/band,
- water target extraction: once per VZA,
- atmosphere pass 2: once per case/band,
- RAA: Fourier reconstruction only.

The regression script `scripts/regression_immutable_cache_reuse.sh` compares
cache-on, surface-cache-off, transport-cache-off, and all-cache-off outputs and
requires byte identity.

## 5. Threading and lifetime

The surface cache is thread-local and bounded at 128 entries. The larger bound prevents pass-1/pass-2 thrashing when aerosol runs retain grid-grid and grid-solar slabs through Fourier mode 16.  The transport stencil is owned
by the corresponding `rt_atm_t` and released by `rt_atm_free()`.  No cache is
shared mutably between OpenMP workers, and no locks are introduced.
