# OCRT v1.16 — TSM vector-phase preparation and SOS cache

## 1. Scope

This release is a numerical-performance refactor only. It does **not** change:

- the four Ahn TSM species or their `.mie` files;
- TSM `a`, `b`, or `bb` conversion;
- P11/P12/P33 interpolation, normalization, or truncation;
- generalized vector moments (`beta`, `gamma`, `alpha`, `zeta`);
- SOS source equations, accumulation order, transport, convergence criteria,
  or output formatting;
- any public CLI option or default.

The physical result is therefore required to remain bit-identical when the
same compiler, optimization flags, input data and thread count are used.

## 2. Previous hot-path behavior

The vector phase functions and their generalized moments were already computed
outside the SOS scattering-order loop. However, the fast source builder still
performed two order-dependent operational tasks:

1. allocate and release a large transpose/accumulator arena once per SOS source
   application;
2. repack the same 18 particulate + 18 Rayleigh direction-pair kernel scalars
   for every outgoing quadrature direction at every scattering order.

These tasks did not change physics, but their cost scaled with the number of
SOS orders.

## 3. v1.16 execution model

For each Fourier mode `m`, after all phase tables are final:

```text
P11/P12/P33
    -> beta/gamma/alpha/zeta
    -> Fourier phase tables
    -> rt_solver_prepare_sos_source_operator()       [once per m]
         - pack all (k, jp) vector-kernel coefficients
         - allocate/reuse transpose + accumulator arena
    -> SOS order loop
         - transpose current I/Q/U field into reusable arena
         - apply the immutable packed operator
         - perform vertical transport and convergence test
```

The prepared operator is stored in `rt_legendre_workspace_t`. The cache is
invalidated whenever any generalized Legendre or phase-Fourier table is
recomputed. The existing fallback path retains the prior one-call allocation
behavior only for external low-level callers that bypass the preparation API.
The OCRT atmosphere and water production paths always call the preparation API.

## 4. Cache identity and lifetime

The effective phase state is determined by the containing solve and Fourier
mode, including:

```text
medium / selected phase source
wavelength
quadrature and Legendre settings
phase truncation / moment construction
Fourier mode m
```

Within one finalized state, TSM concentration changes scattering magnitude but
not the normalized species phase shape. The packed operator is read-only during
all SOS orders for that mode. It is released with the Legendre workspace.

## 5. Dynamic verification

With the representative high-scattering case below:

```text
species       red_clay
TSM           2 g m^-3
wavelength    443 nm
geometry      SZA/VZA/RAA = 30/30/90 deg
water m       0..8, early Fourier exit disabled
SOS cap       2 versus 100
OMP threads   1
```

`OCRT_DUMP_SOS_PHASE_CACHE=1` produced exactly nine build records in both runs:
`m=0` through `m=8`. The two trace files are byte-identical and have SHA-256:

```text
89146542445cd110649942f6b09ce559f1a95225ed0fbc9ca9b5e13d1994b833
```

A test-only allocation interposer measured:

| Build | SOS cap | malloc | calloc | realloc |
|---|---:|---:|---:|---:|
| v1.15 baseline | 2 | 1959 | 251 | 3 |
| v1.15 baseline | 100 | 2405 | 251 | 3 |
| v1.16 final | 2 | 1941 | 251 | 5 |
| v1.16 final | 100 | 1941 | 251 | 5 |

Thus the final production path adds **zero allocations** when the SOS cap rises
from 2 to 100. The two `realloc` calls relative to the old layout occur during
one-time workspace preparation, not during the order loop.

## 6. Result-invariance audit

### 6.1 Direct operator equivalence

`tests/test_sos_phase_cache_equivalence.c` applies the same synthetic I/Q/U
fields through both source-builder routes in one process:

1. the legacy per-call fallback packing path;
2. the prepared operator path;
3. a second prepared-path call using the same cached operator and arena.

The test covers 7 Fourier modes and 4 independent field seeds.  Every double
is compared by its raw 64-bit representation:

```text
modes                     7
field seeds               4
I/Q/U array comparisons   168
values per array          238
bit mismatches            0
```

This isolates the changed source-operator layer from the rest of the RT
solver and proves that both first use and cache reuse preserve the old source
field exactly.

### 6.2 Packaged v1.16 executable versus v1.15 baseline

The final packaged v1.16 executable was run directly against the v1.15
baseline under identical `OMP_NUM_THREADS=1` conditions over 124 cases.  The
matrix includes:

- all four Ahn species at low concentration: 5 wavelengths × 4 geometries
  (80 cases total);
- medium-concentration TSM cases;
- TSM + Rayleigh and TSM + Rayleigh + aerosol paths;
- pure-water, direct-IOP and atmosphere-only paths;
- explicit water SOS caps 2/30/100 and atmospheric SOS caps 2/20/40.

For every case, exit code, stdout and stderr are byte-identical:

```text
cases                     124
stdout exact              124
stderr exact              124
exit-code exact           124
max absolute difference   0
max relative difference   0
max ULP difference        0
comparison CSV SHA-256    f24b43648c7096a8b91d198f6f249f36f6a8bde395924d8c98267f0f642b3382
evidence SHA-256          96e0f7206dde20fbea221fde9336005a29def3b1db84d7e1eea00ef1e6699c4e
```

Evidence is stored in:

```text
validation/tsm_phase_cache_v116/final_v116_vs_v115_124.csv
validation/tsm_phase_cache_v116/final_v116_vs_v115_124_summary.txt
```

### 6.3 Broader linked numerical audit

Before changing the visible version identifier, the same numerical source was
run against the v1.15 baseline over 166 cases covering:

- all four TSM species;
- 412, 443, 555, 670 and 700 nm;
- nadir, off-plane, oblique and near-specular geometries;
- calm and rough surfaces;
- low, medium and high TSM concentrations;
- TSM-only, TSM + Rayleigh, and TSM + Rayleigh + aerosol paths;
- pure-water and direct-IOP branches;
- atmosphere-only Rayleigh and aerosol paths;
- explicit water SOS caps 2, 30 and 100;
- explicit atmospheric SOS caps 2, 20 and 40.

For every case, exit code, stdout and stderr were byte-identical:

```text
cases                     166
stdout exact              166
stderr exact              166
exit-code exact           166
max absolute difference   0
max relative difference   0
max ULP difference        0
artifact-matrix SHA-256   041078fd72b0299d8e3b8dfd2711fb1517b2bb632a40173849fbaa50f760d9d1
```

A separate non-fast strict-build audit covered six representative water and
atmosphere branches and was also byte-identical. An OpenMP four-thread audit
covered eight representative TSM, pure-water, direct-IOP, aerosol and high-order
water/atmosphere paths; exit code, stdout and stderr were byte-identical in all
eight cases. AddressSanitizer and UBSan reported no issue on TSM
low/high-scattering and atmospheric cases.

The 166-case executable retained the prior v1.15 version text.  As a linkage
check, the final source was rebuilt with only `V2_VERSION` changed back to that
text; the resulting executable reproduced the audited binary byte-for-byte
(SHA-256
`6bd15eb429635dddbbb3e1894b90c8108045855e07596941af10ef73c54f4a66`).
The numerical code in that broader matrix and the packaged v1.16 release is
therefore identical; only the embedded version string differs.  The expensive
default-resolution high-TSM cases are covered by this linked 166-case matrix,
while the packaged v1.16 executable itself was directly rerun for the 124-case
matrix above.

## 7. Build diagnostics and release reproducibility

The documented production FAST release build completes without errors and
without warning lines. The source list is sorted in the C locale before linking;
a clean final-source rebuild is therefore byte-identical to both packaged
executables. All three have SHA-256:

```text
26e46b00fdcd9b6d89c368efd4d333b4672c5cf77415f9885ef34da274369b8f
```

Separate `-Wall -Wextra` comparison builds retain the same pre-existing
diagnostics on both sides of the refactor:

```text
fast audit   baseline 51, candidate 51
strict audit baseline 54, candidate 54
```

No warning was introduced by v1.16; the audit warnings predate this
performance-only change.
