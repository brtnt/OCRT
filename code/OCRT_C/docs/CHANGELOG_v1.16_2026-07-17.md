# OCRT v1.16 changelog — TSM phase-operator cache

- Preserved the v1.15 Ahn four-species TSM optical model and all public CLI
  semantics.
- Added an explicit per-Fourier-mode preparation step for the vector SOS source
  operator after P11/P12/P33 have been converted to the OCRT
  beta/gamma/alpha/zeta representation.
- Packed the immutable aerosol/mineral and Rayleigh direction-pair kernels once
  per Fourier mode instead of rebuilding them at every SOS scattering order.
- Moved the vector-field transpose/accumulator work arena out of the SOS order
  loop and made it workspace-owned and reusable.
- Removed production-path `malloc/free` from the fast vector SOS source builder;
  the low-level compatibility fallback remains available for callers that do
  not invoke the preparation API.
- Added cache invalidation whenever generalized Legendre functions or vector
  Fourier phase tables are rebuilt.
- Added `OCRT_DUMP_SOS_PHASE_CACHE=1` diagnostic output and a regression script
  that proves phase-operator build count is independent of SOS order count.
- Corrected the TSM smoke-test conflict-message matcher; runtime behavior was
  already correct.
- Added a direct source-operator equivalence test: 7 Fourier modes × 4
  synthetic input fields × I/Q/U × fallback/prepared/reuse paths. All 168
  compared arrays are bit-identical.
- Compared the packaged final v1.16 executable directly with the v1.15
  baseline over 124 release cases. Exit status, stdout and stderr are
  byte-identical in every case; maximum absolute, relative and ULP differences
  are all zero.
- Retained the broader 166-case audit of the same numerical source built with
  the prior version text. It includes all four TSM species, high-TSM default
  resolution cases and water SOS caps up to 100; all 166 cases are exact.

- Proved release linkage: rebuilding the final source with only the version
  identifier reverted reproduced the 166-case-audited executable byte-for-byte.
- Rechecked the installed four Ahn `.mie` files against the user-supplied ZIP;
  all four SHA-256 values match and all four retired legacy mineral files remain
  absent.
- Verified eight representative production paths under `OMP_NUM_THREADS=4`;
  baseline and v1.16 exit code, stdout and stderr are byte-identical.
