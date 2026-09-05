# Python full-grid aerosol object fix — validation report

## Scope

C reference: `OCRT-v1.2-2026-07-21-KST-fullgrid-aerosol-object-fix`  
Python result: `pyOCRT-v1.2-2026-07-21-fullgrid-aerosol-object-fix`

The Python atmosphere-ocean full-grid path now prepares one immutable aerosol runtime before the geometry loop and passes that same object explicitly to every `(VZA, RAA)` cell and to both atmospheric passes in the coupled solver.

## Code changes

- Added frozen `AerosolRuntime` in `ocrt_py/aerosol.py`.
- Added source file fingerprint (`absolute path`, `mtime_ns`, `size`) to `MieData`.
- Centralized AOD spectral scaling, SSA interpolation, P11/P12/P33 wavelength PCHIP, log-linear truncation, tau/SSA transformation, and vector moments in `prepare_aerosol_runtime()`.
- Added an explicit `aerosol_runtime` argument to atmospheric and coupled-ocean solver paths.
- Added an immediate error when `AOD > 0` but neither Mie data nor a prepared runtime is supplied.
- Added `run_ocean_aerosol_full_grid()` and `ocrt_fullgrid.py`.
- Changed the production batch phase interpolation from linear to PCHIP and strengthened the phase-cache key with Mie file identity and numerical options.
- Added Rrs/rrs Q and U fields to the coupled result path.

## Validation results

| Test | Result | Numerical result |
|---|---|---|
| Python single vs full-grid, same cell | PASS | max absolute difference = `0.0`; all compared values bit-identical |
| AOD=0 original Python vs patched Python | PASS | max absolute difference = `0.0` |
| Execution-order independence A–D | PASS | max absolute difference = `0.0` |
| 100-cell object lifecycle | PASS | Mie reads `1`; runtime preparations `1`; cell solves `100`; one shared object |
| Aerosol model sensitivity | PASS | C50, M50C, r80f20v01 give distinct TOA I/Q/U |
| AOD scaling | PASS | AOD865 = 0, 0.05, 0.2, 0.5 give distinct outputs |
| C regression executables | PASS | PCHIP and external-bottom bounds tests pass |
| C vs Python aerosol runtime | PASS | max relative difference `1.034e-15` |

For C50 at 490 nm and AOD865=0.2:

| Runtime quantity | C | Python | absolute difference |
|---|---:|---:|---:|
| AOD at 490 nm | 0.27185337499960666 | 0.27185337499960666 | 0 |
| effective aerosol tau | 0.24387332404306053 | 0.24387332404306056 | 2.78e-17 |
| effective aerosol SSA | 0.95370972550096 | 0.95370972550096 | 0 |
| forward truncation coefficient | 0.21476494624842979 | 0.21476494624842957 | 2.22e-16 |

## C–Python whole-RT result

The aerosol-object change itself is aligned, but the complete C and Python coupled RT outputs are **not yet absolutely aligned**. This is not hidden in the result.

At the representative VZA=30°, RAA=90° diagnostic grid:

| AOD865 | Quantity | C | Python | relative difference |
|---:|---|---:|---:|---:|
| 0 | TOA I | 0.1107740 | 0.1052999 | 4.94% |
| 0 | Rrs I | 0.01674055 | 0.01484026 | 11.35% |
| 0 | rrs I | 0.03023010 | 0.02637354 | 12.76% |
| 0.2 | TOA I | 0.1332829 | 0.1285626 | 3.54% |
| 0.2 | Rrs I | 0.01564766 | 0.01383154 | 11.61% |
| 0.2 | rrs I | 0.02843798 | 0.02476513 | 12.92% |

Because the discrepancy already exists at AOD=0, it cannot be attributed to the aerosol-object forwarding fix. It is a pre-existing difference in the Rayleigh/surface/water coupling implementation or numerical grid path.

The aerosol response direction is nevertheless consistent:

| Quantity | C change, AOD865 0→0.2 | Python change | difference |
|---|---:|---:|---:|
| TOA I | +20.32% | +22.09% | +1.77 percentage points |
| Rrs I | -6.53% | -6.80% | -0.27 percentage points |
| rrs I | -5.93% | -6.10% | -0.17 percentage points |

These whole-RT comparisons used a deliberately reduced, non-converged grid (`water_orders=8`, `converged=0`) for execution-path regression. A production-resolution Python coupled solve exceeded the CPU execution limit, so no converged production-resolution whole-RT agreement is claimed here.

## Conclusion

The Python full-grid aerosol object defect class is closed:

- aerosol runtime is explicit;
- geometry-independent preparation occurs once;
- every cell uses the same immutable object;
- AOD>0 cannot silently fall back to `None`;
- output is independent of prior execution order;
- Python single and full-grid are bit-identical;
- C and Python aerosol preprocessing is effectively machine-identical.

Full C–Python coupled RT parity remains a separate follow-up defect because the AOD=0 baseline already differs.
