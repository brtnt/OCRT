# OCRT v1.2 IOP/Kd full-grid CSV and Kd(0−) correction integration

## Release

```text
OCRT-v1.2-2026-07-21-KST-iop-kd-fullgrid-csv-fix
```

Date: 2026-07-21

Immediate base:

```text
OCRT-v1.2-2026-07-21-KST-water-mmax-pchip-rt-trans-output-fix
```

## Work-order audit

The two uploaded work-order files were byte-identical:

```text
SHA-256: e48e798690f8cec2b4bc77daef9a805ad35dd1ecc456753f2d9e0c548f7f73e1
Size: 16,228 bytes each
```

One canonical copy is retained as:

```text
docs/C_WORKORDER_IOP_KD_FULLGRID_CSV_2026-07-21.md
```

## Confirmed defect

The former `Kd_0minus` diagnostic used the water-SOS field at the surface and first water level but omitted the unscattered transmitted atmospheric skylight term. The surface `Ed_0minus_water` output later received that term in the atmosphere–ocean orchestrator, while the first-level Ed used by Kd did not. This inconsistent pair could produce negative near-surface Kd even for physically absorbing water.

The corrected computation is:

```text
Ed_sky(k) = 2π (F_sun_water/π) Σ I_sky(μc) exp(-τk/μc) μc wc

Ed(0−) = Ed_water_field(0) + Ed_sky(0)
Ed(z1) = Ed_water_field(z1) + Ed_sky(z1)
Kd(0−) = -ln[Ed(z1)/Ed(0−)] / z1
```

## Source changes

- `src/rt_water_rt.h`
  - adds `a_chl_used`;
  - exposes `Ed_level1_water`, `tau_level1_used`, and `z_level1_used`.
- `src/rt_water_rt.c`
  - stores phytoplankton-only absorption `a_phyto` as `a_chl_used`;
  - stores the first-level Ed and its optical/geometric depth.
- `src/rt_types.h`
  - adds final-result field `a_chl_used`.
- `src/rt_solver.c`
  - reconstructs unscattered transmitted skylight on the water quadrature;
  - recomputes Kd from skylight-consistent surface and first-level Ed;
  - transfers `a_chl_used` to the public result.
- `src/rt_io.c`
  - adds `a_chl` to single-geometry stdout.
- `src/main.c`
  - preserves the original 46 ocean full-grid columns in their original order;
  - appends 12 columns:
    `a_w,b_w,bb_w,a_chl,a_phyto_detritus,b_phyto_detritus,bb_phyto_detritus,a_dom,a_min,b_min,bb_min,Kd0minus`.
- `scripts/regression_iop_kd_fullgrid_csv.sh`
  - pins the reference Kd, IOP values, Rrs/TOA/Ed/Lu invariance, and 58-column CSV schema.

## Reference-case correction

For Chl=1.2 mg m⁻³, TSM=3.5 g m⁻³, aDOM(440)=0.04 m⁻¹, 490 nm, C50 aerosol, SZA/VZA/RAA=30/30/90°, wind=3 m s⁻¹:

```text
Kd0minus before = -0.854219 m^-1
Kd0minus after  = +0.342726 m^-1
a_total + bb_total ≈ 0.325322 m^-1
```

The corrected value is positive and close to the expected `a+bb` scale.

## Rrs sensitivity test

Four representative cases were run with the immediate pre-fix binary and the corrected binary using the same optimized release build settings.

| Case | Rrs(0+) I before | Rrs(0+) I after | change | rrs(0−) I before | rrs(0−) I after | change | Kd before | Kd after |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `case_turbid_490` | 0.04313167 | 0.04313167 | 0.0e+00 | 0.07127285 | 0.07127285 | 0.0e+00 | -0.854219 | +0.342726 |
| `case_clear_490` | 0.02980273 | 0.02980273 | 0.0e+00 | 0.05091275 | 0.05091275 | 0.0e+00 | -0.046480 | +0.027803 |
| `case_pure_443` | 0.01677798 | 0.01677798 | 0.0e+00 | 0.03042000 | 0.03042000 | 0.0e+00 | +0.009509 | +0.010919 |
| `case_moderate_555` | 0.03418803 | 0.03418803 | 0.0e+00 | 0.05835532 | 0.05835532 | 0.0e+00 | -0.094352 | +0.247207 |

Result:

```text
max |Δ Rrs(0+) I| = 0
max |Δ rrs(0−) I| = 0
max |Δ TOA rho_I| = 0
max |Δ Ed(0−)| = 0
max |Δ Lu(0−)| = 0
```

All I/Q/U Rrs and rrs values, TOA reflectance, Ed and Lu were unchanged at the printed precision. This is expected: the correction changes the diagnostic derivation of Kd and adds output columns; it does not feed Kd back into the water RT solution.

## Full-grid CSV test

```text
columns before = 46
columns after  = 58
rows tested    = 4
```

The original 46-column prefix and all values in those columns were unchanged. The 12 new columns were appended, avoiding positional breakage for readers that still use the legacy prefix.

## Regression status

- new IOP/Kd/full-grid regression: PASS (`single_Kd=0.342726`, 4 grid rows, 58 columns);
- water external-bottom-source mode-bound: PASS;
- RT-derived transmittance output: PASS;
- Mie phase-wavelength PCHIP: PASS;
- radiometry flux split: PASS;
- water input branch: 14/14 PASS;
- CCRR Chl: 11/11 PASS;
- organic Chl: 10/10 PASS;
- Ahn TSM: 6/6 PASS;
- TSM phase cache: PASS;
- PSSA numeric fixes: 4/4 PASS;
- PSSA layer warning: PASS;
- US62 atmospheric profile: PASS;
- ASan/UBSan representative coupled-water case: PASS.

## Validation files

```text
validation/iop_kd_fullgrid_csv_fix_2026-07-21/
  rrs_before_after_O3.csv
  rrs_before_after_O3.json
  fullgrid_before.csv
  fullgrid_after.csv
  validation_summary.json
  regression and sanitizer logs
```


## Release binary

```text
SHA-256: a0fb765639e674c6f48d38f743db2216654ba725fc56ff254d951ee2eaa06528
Version: OCRT-v1.2-2026-07-21-KST-iop-kd-fullgrid-csv-fix
```

Fast package verification passed. The long reference IOP/Kd regression result is shipped under the validation directory and may be rerun with `OCRT_RUN_LONG_KD_REGRESSION=1 ./verify_package.sh`.

## Distribution

Canonical repository: `https://github.com/brtnt/OCRT`

The package remains subject to the root Academic and Non-Commercial Use License.
