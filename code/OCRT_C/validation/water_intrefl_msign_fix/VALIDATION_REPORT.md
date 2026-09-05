# OCRT water internal-reflection Fourier-sign bugfix validation

Date: 2026-07-21  
C version: `OCRT-v1.2-2026-07-21-KST-water-intrefl-msign-fix`  
Python version: `pyOCRT-v1.2-2026-07-21-water-intrefl-msign-fix`

## 1. Root cause and scope

The water-side internal-reflection top boundary used `(-1)^m`.  In the OCRT
water coordinate convention the upward and downward hemispheres share the same
physical azimuth.  Specular reflection flips the vertical propagation component
but preserves azimuth, so every Fourier order must use `msign=+1`.

C changes:

- `src/rt_solver.c::add_flat_intrefl_downfield_order_complete`
- `src/rt_solver.c::rt_solver_sos_pol_intrefl`
- `src/rt_solver.c::rt_solver_sos_pol_intrefl_rough`

The previously unclassified use near v1.2 line 1996 was inspected.  It is the
same flat water-side internal-reflection boundary and was corrected as well.
`rt_fourier_pi_shift_sign()` itself was not changed.

Python changes:

- `ocrt_py/sos.py::sos_pol_intrefl_rough`
- `ocrt_py/atmos_batch.py::sos_water_intrefl_batch`
- `ocrt_py/lutbatch.py::solve_batch` water internal-reflection loop

## 2. Canonical validation configuration

- Red-clay inorganic particles, TSM = 5 g m-3
- Chl = 0, aDOM440 = 0
- wavelength = 555 nm
- SZA = 30 deg, RAA = 90 deg
- atmosphere isolated out: pressure = 0, AOD = 0
- water positive Gauss nodes = 48
- water Fourier maximum = 4
- wind = 3 m s-1 unless stated otherwise
- one CPU thread

OSOAA values were taken from the supplied work instruction; OSOAA was not rerun.

## 3. C vza sweep

| VZA air | Before DoLP | After DoLP | OSOAA DoLP | After-OSOAA | I unchanged | Q unchanged |
|---:|---:|---:|---:|---:|:---:|:---:|
| 20 | 0.029975 | 0.033513 | 0.033250 | +0.79% | yes | yes |
| 30 | 0.041310 | 0.047165 | 0.046650 | +1.10% | yes | yes |
| 45 | 0.059434 | 0.068813 | 0.067920 | +1.32% | yes | yes |
| 60 | 0.075737 | 0.088632 | 0.087310 | +1.51% | yes | yes |
| 70 | 0.083948 | 0.098910 | 0.097360 | +1.59% | yes | yes |
| 75 | 0.086960 | 0.102757 | 0.101110 | +1.63% | yes | yes |

All six angles are inside the supplied +/-2% acceptance band.  Maximum absolute
DoLP difference from OSOAA was 1.63%.

At VZA=60 deg, wind=3:

| Quantity | Before | After | OSOAA reference |
|---|---:|---:|---:|
| I = rrs0minus_I | 0.09145276 | 0.09145276 | -- |
| Q = rrs0minus_Q | -0.002221415 | -0.002221415 | -- |
| U = rrs0minus_U | 0.006560469 | 0.007795343 | -- |
| Q/I | -0.0242903 | -0.0242903 | -0.02370 |
| U/I | 0.0717362 | 0.0852390 | 0.08403 |
| DoLP | 0.0757370 | 0.0886324 | 0.08731 |

## 4. Required invariants

- C `rrs0minus_I`: exactly unchanged at all tested VZA values.
- C `rrs0minus_Q`: exactly unchanged at all tested VZA values.
- C `Rrs0plus_I/Q`: exactly unchanged at all tested VZA values.
- C nadir `rrs0minus_U`: exactly unchanged and effectively zero.
- Only off-nadir U and quantities derived from U changed.
- SOS order count and convergence flags were unchanged (198 orders, converged).

Nadir after-fix values:

- `rrs0minus_I = 8.402686e-02`
- `rrs0minus_Q = 1.420711e-03`
- `rrs0minus_U = -5.219607e-19`

The small difference from the document's last printed digits reflects the newer
LUT-cache-performance baseline, not this sign change.

## 5. Flat and rough water paths

VZA=60 deg:

| Surface path | Before DoLP | After DoLP | I/Q unchanged |
|---|---:|---:|:---:|
| rough, wind=3 | 0.075737 | 0.088632 | yes |
| flat, wind=0 | 0.074879 | 0.088720 | yes |

This confirms that both C internal-reflection implementations were corrected.

## 6. Above-water polarization

At VZA=60 deg, wind=3, above-water `Rrs(0+)` DoLP increased from 0.080346 to
0.092594.  I and Q stayed unchanged; only U changed, as expected from the
already-correct water-to-air transmission operator.

## 7. Python validation

Native Python water solver, VZA=60 deg, wind=3:

| Quantity | Before | After | Exact invariant? |
|---|---:|---:|:---:|
| I | 0.09150796848 | 0.09150796848 | yes |
| Q | -0.002229772815 | -0.002229772815 | yes |
| U | 0.006574623020 | 0.007812498543 | intentional change |
| DoLP | 0.07586711830 | 0.08878429493 | intentional change |
| Rrs(0+) I | 0.05992300289 | 0.05992300289 | yes |
| Rrs(0+) Q | 0.002166849971 | 0.002166849971 | yes |

Python nadir I/Q/U and Rrs I/Q/U were physically identical before and after.
Elapsed-time metadata differs and is not a physical output.

After the fix, C-Python below-water agreement at VZA=60 deg was:

- I relative difference: 0.060%
- Q relative difference: 0.376%
- U relative difference: 0.220%
- DoLP relative difference: 0.171%

## 8. Regression and build results

C:

- production release build: PASS
- new water-internal-reflection regression: PASS
- Mie P11/P12/P33 wavelength PCHIP: PASS
- external bottom mode bounds/shape guard: PASS
- IOP/Kd full-grid CSV regression: PASS
- LUT exact-cache byte comparison: PASS

Python:

- `compileall`: PASS
- unit tests: 6/6 PASS
- native single odd-m boundary sign: PASS
- NumPy/CuPy batch odd-m boundary sign: PASS
- LUT-batch source audit: PASS
- previous aerosol-object regression: 3/3 PASS

A repository-wide `-Wall -Wextra -Werror` build still fails on pre-existing
warnings in `rt_solver.c` and `rt_water_rt.c` (misleading indentation, unused
legacy helpers, and maybe-uninitialized diagnostics).  No new warning is caused
by this patch; the normal release build is clean and passes.  These unrelated
warnings were not modified as part of this narrowly scoped physics fix.

## 9. Runtime

The scalar multiply already existed; only its constant changed.  The canonical
C VZA=60 rough case measured 6.28 s both before and after.  No measurable runtime
cost was introduced.

## 10. Conclusion

The bugfix is accepted.  It changes only odd-m water internal-reflection
contributions, restores off-nadir Stokes U/DoLP, preserves I/Q/nadir outputs,
and brings the C and Python water polarization solutions into close agreement
with the supplied OSOAA references.
