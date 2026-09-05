# Stage-2 water RAA output fix validation package

## Scope

This directory validates the water-output RAA correction applied to `OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix`.

## Acceptance status

- 72/72 OCRT water execution units converged.
- 600/600 OCRT–OSOAA comparison cells were generated.
- Corrected water output equals the mirrored baseline exactly:
  `I_new(r)=I_old(180-r)`, `Q_new(r)=Q_old(180-r)`,
  `U_new(r)=-U_old(180-r)`.
- Same-RAA TOA I/Q/U and irradiances are unchanged.
- Reference-independent water/TOA ratio spread decreased from a factor of
  approximately 2.16 to 3.378% relative spread.
- Wind-zero and wind-positive paths passed.
- Full-grid and single-geometry outputs passed.
- Runtime is statistically unchanged.

## Official OSOAA comparison rule

Use one of these equivalent mappings:

```text
A: Phi_OSOAA=(180-RAA_OCRT) mod 360 and reverse OSOAA U
B: Phi_OSOAA=(RAA_OCRT+180) mod 360 and do not reverse U
```

Official `Rrs(0+)` must use the level-27 underwater Fourier field followed by
a separate TWA transfer calculation.  The attached level-26 water-minus-black
subtraction is retained only as a diagnostic method audit.

## Files

- `MATCHED_OCRT_OSOAA_600CELLS_CORRECTED_RAA.csv`: final corrected comparison.
- `COMPONENT_METRICS.csv`: scalar summaries; plots remain the primary review format.
- `RUN_MANIFEST_72CASES.csv`: exact 72-run input contract.
- `MIRROR_IDENTITY_AND_TOA_INVARIANCE_600CELLS.csv`: bit-level invariance audit.
- `INTERNAL_RAA_CONSISTENCY_B443_SZA40_VZA60.csv`: reference-independent check.
- `ATTACHED_OSOAA_METHOD_VS_APPROVED_TWA_40CELLS.csv`: method discrepancy audit.
- `ATTACHED_REPORT_METRICS_REPRODUCED.csv`: reproduction of the attached report's
  own max-reference-normalized statistics.
- `RUNTIME_BENCHMARK_RAW.csv`, `RUNTIME_BENCHMARK_SUMMARY.csv`: paired timings.
- `figures/`: independent I/Q/U scatterplots for Rrs and rrs.

## Scatterplot sets

- `all600_*`: all 600 cells.
- `purewater150_*`: pure water, all SZA.
- `purewater84_sza_le40_*`: pure water, SZA ≤ 40°.
- `internal_raa_consistency_before_after.png`: internal before/after diagnostic.

I intensity plots show MAPE. Q and U plots show
`RMS[100*(OCRT-OSOAA)/OSOAA_I]`, avoiding MAPE instability near zero.
