# OCRT Native Coupled Angular LUT Recovery — Validation Report

Date: 2026-07-27 (KST)

## 1. Scope

This release restores the intended all-angle LUT structure for the fully coupled atmosphere–ocean path in both OCRT C and OCRT Python.

For one case and one wavelength, the native path now performs:

1. atmospheric pass-1 SOS once;
2. air-to-water coupling once;
3. water SOS once for all requested VZAs;
4. authoritative water target projection once per VZA;
5. water-to-air coupling and atmospheric pass-2 SOS once;
6. RAA dependence by Fourier reconstruction only.

The old per-cell path remains available only for regression comparison:

- C: `OCRT_NATIVE_COUPLED_LUT_OFF=1`
- Python: `OCRT_PY_NATIVE_COUPLED_LUT_OFF=1`

## 2. C implementation contract

Representative 36-cell diagnostic:

```text
[native-coupled-lut] calls=3 water_cold=1 water_views=3 exact_near_nadir=1 cells=36
```

Interpretation:

- coupled field-building phases: 3;
- cold water SOS: 1;
- water target projections: equal to VZA count, not VZA×RAA count;
- RAA cells: Fourier reconstruction and output assembly only.

## 3. C native-versus-authoritative replay

| Grid | Rows | Maximum difference, all CSV fields | Maximum difference, physical radiance/reflectance fields |
|---|---:|---:|---:|
| 3 VZA × 12 RAA | 36 | 0 | 0 |
| 7 VZA × 40 RAA | 280 | 9.9920e-14 (`T_total_up_view`) | 0 |
| 35 VZA × 72 RAA | 2520 | 1.0003e-13 (`T_total_up_view`) | 0 |

The following physical fields were unchanged in the comparisons:

- TOA rho I/Q/U;
- Rrs(0+) I/Q/U;
- rrs(0-) I/Q/U;
- Lu/Ed/Eu boundary fields;
- IOP and Kd values;
- convergence metadata.

The ~1e-13 residual at large grids is limited to a derived diagnostic transmittance and is floating-point reconstruction roundoff.

## 4. C RAA scaling and performance

### 4.1 RAA scaling

Seven VZAs were held fixed while the number of RAAs was increased by 20×.

| Grid | Median/representative runtime |
|---|---:|
| 7 VZA × 2 RAA | 2.36 s |
| 7 VZA × 40 RAA | 2.33 s |

The difference is measurement noise. Increasing RAA count does not trigger additional RT solves.

### 4.2 Native versus legacy replay

| Grid | Native | Legacy cell replay | Speedup |
|---|---:|---:|---:|
| 36 cells | 4.81 s | 6.46 s | 1.34× |
| 280 cells | 2.39 s | 2.50 s | 1.05× |
| 2520 cells | 7.74 s | 8.37 s | 1.08× |

Absolute timings depend strongly on warm caches, numerical settings, and output I/O. The structural acceptance criterion is that RAA growth is effectively free and that expensive SOS fields are solved once.

## 5. Python implementation and validation

Python now provides a native coupled-grid path with the same contract:

```text
atmos_pass1_solves = 1
water_solves       = 1
atmos_pass2_solves = 1
cell_solver_calls  = 0
```

### 5.1 Exact native-versus-legacy comparison

Representative 2 VZA × 2 RAA case:

- all output rows exactly equal;
- native median: 1.0950 s;
- legacy median: 4.0318 s;
- speedup: 3.68×.

### 5.2 RAA scaling

| Grid | Median runtime |
|---|---:|
| 7 VZA × 2 RAA | 1.3880 s |
| 7 VZA × 40 RAA | 1.3726 s |

Again, increasing RAA count by 20× did not increase runtime.

### 5.3 Python tests

```text
18 passed
```

The suite includes native coupled-LUT equality against the authoritative legacy cell replay, aerosol object lifetime, odd-m internal reflection, and existing Stage-2 production tests.

## 6. Frozen OCRT validation preservation

The recovered C source passed:

- public water RAA regression;
- direct-glint RAA branch regression;
- Stage-2 air–water FIX1/FIX2/FIX3;
- diffuse-top U-column closure;
- exact-pole surface limit;
- water internal-reflection m-sign regression;
- external bottom-source bounds and shape guards;
- Mie P11/P12/P33 wavelength PCHIP;
- Rrs/rrs IQU and IOP/Kd full-grid CSV regression;
- EAP Stage-2 production contract;
- coupling high-resolution checks at n_mu=64 and 96;
- 18 non-Chl frozen cases with stdout and stderr byte-identical.

Representative IOP/Kd regression:

```text
IOP_KD_FULLGRID_CSV PASS single_Kd=0.360659 grid_rows=4 cols=60
```

## 7. OSOAA reference validation

The latest canonical OSOAA package was not modified by this optimization.

Validation:

```text
SMOKE PASS — 11/11 reference points reproduced (rel < 1e-9)
```

The full OSOAA manifest also passed SHA-256 verification.

## 8. Performance-safety assessment

The optimization does not alter phase functions, IOPs, interface Mueller matrices, scattering orders, convergence thresholds, or RAA/U conventions. It changes orchestration and reuse of already validated Fourier fields.

No additional RT calculations are introduced. RAA reconstruction consists only of Fourier sums and CSV/output assembly.

## 9. Limitations

A complete full-grid ASan/UBSan run was not completed within the execution limit. Release builds, unit tests, frozen-value regressions, high-resolution coupling tests, and native-versus-legacy numerical comparisons all passed. The sanitizer limitation is recorded rather than treated as a pass.
