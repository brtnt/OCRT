# Stage-2 water RAA output fix — validation package

## Scope

This directory validates the water-output RAA correction applied in
`OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix`.

## Acceptance status

- 72/72 OCRT water execution units converged.
- 600/600 OCRT–reference (OSOAA) comparison cells were generated.
- The corrected water output equals the mirrored baseline exactly:
  `I_new(r) = I_old(180−r)`, `Q_new(r) = Q_old(180−r)`, `U_new(r) = −U_old(180−r)`.
- Same-RAA TOA I/Q/U and irradiances are unchanged.
- The reference-independent water/TOA ratio spread decreased from a factor of about 2.16 to a
  relative spread of 3.378 %.
- Wind-zero and wind-positive paths passed.
- Full-grid and single-geometry outputs passed.
- Runtime is statistically unchanged.

## Official comparison rule with the reference model

Use one of these equivalent mappings:

```text
A: Phi_ref = (180 − RAA_OCRT) mod 360 and reverse the reference U
B: Phi_ref = (RAA_OCRT + 180) mod 360 and do not reverse U
```

The official `Rrs(0+)` must use the level-27 underwater Fourier field followed by a separate TWA
transfer calculation. The attached level-26 water-minus-black subtraction is kept only as an audit
of the method discrepancy.

## Files

- `MATCHED_OCRT_OSOAA_600CELLS_CORRECTED_RAA.csv`: final corrected comparison.
- `COMPONENT_METRICS.csv`: scalar summaries; the plots remain the primary review format.
- `RUN_MANIFEST_72CASES.csv`: exact 72-run input contract.
- `MIRROR_IDENTITY_AND_TOA_INVARIANCE_600CELLS.csv`: bit-level invariance audit.
- `INTERNAL_RAA_CONSISTENCY_B443_SZA40_VZA60.csv`: reference-independent check.
- `ATTACHED_OSOAA_METHOD_VS_APPROVED_TWA_40CELLS.csv`: method discrepancy audit.
- `ATTACHED_REPORT_METRICS_REPRODUCED.csv`: reproduction of the attached report's own
  max-reference-normalized statistics.
- `RUNTIME_BENCHMARK_RAW.csv`, `RUNTIME_BENCHMARK_SUMMARY.csv`: paired timings.
- `figures/`: independent I/Q/U scatter plots of Rrs and rrs.

## Scatter-plot sets

- `all600_*`: all 600 cells.
- `purewater150_*`: pure water, all SZA.
- `purewater84_sza_le40_*`: pure water, SZA ≤ 40°.
- `internal_raa_consistency_before_after.png`: internal before/after diagnostic.

Intensity (I) plots show the MAPE. Q and U plots show `RMS[100 · (OCRT − ref) / ref_I]`, which
avoids the instability of the MAPE near zero.

---

# Stage-2 수중 RAA 출력 수정 — 검증 패키지

## 범위

이 디렉터리는 `OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix` 에 적용한 수중 출력 RAA 보정을
검증한다.

## 승인 상태

- OCRT 수중 실행 단위 72/72 수렴.
- OCRT–참조(OSOAA) 비교 셀 600/600 생성.
- 보정된 수중 출력은 거울상 기준값과 정확히 같다:
  `I_new(r) = I_old(180−r)`, `Q_new(r) = Q_old(180−r)`, `U_new(r) = −U_old(180−r)`.
- 같은 RAA 의 TOA I/Q/U 와 조도는 변하지 않는다.
- 참조값과 무관한 수중/TOA 비의 퍼짐이 약 2.16배에서 상대 퍼짐 3.378 % 로 줄었다.
- 풍속 0 과 풍속 양수 경로 통과.
- 전체격자와 단일기하 출력 통과.
- 실행 시간은 통계적으로 변화 없음.

## 참조 모델과의 공식 비교 규칙

다음 두 등가 매핑 중 하나를 쓴다.

```text
A: Phi_ref = (180 − RAA_OCRT) mod 360, 참조 U 부호 반전
B: Phi_ref = (RAA_OCRT + 180) mod 360, U 부호 반전 없음
```

공식 `Rrs(0+)` 는 준위 27 수중 푸리에 장에 별도의 TWA 전달 계산을 붙여 구해야 한다. 첨부된 준위 26
물−흑색 차감은 방법 차이 감사용으로만 보존한다.

## 파일

- `MATCHED_OCRT_OSOAA_600CELLS_CORRECTED_RAA.csv`: 최종 보정 비교.
- `COMPONENT_METRICS.csv`: 스칼라 요약. 그림이 주 검토 형식이다.
- `RUN_MANIFEST_72CASES.csv`: 정확한 72회 실행 입력 계약.
- `MIRROR_IDENTITY_AND_TOA_INVARIANCE_600CELLS.csv`: 비트 수준 불변성 감사.
- `INTERNAL_RAA_CONSISTENCY_B443_SZA40_VZA60.csv`: 참조값과 무관한 점검.
- `ATTACHED_OSOAA_METHOD_VS_APPROVED_TWA_40CELLS.csv`: 방법 차이 감사.
- `ATTACHED_REPORT_METRICS_REPRODUCED.csv`: 첨부 보고서의 최대 참조값 정규화 통계 재현.
- `RUNTIME_BENCHMARK_RAW.csv`, `RUNTIME_BENCHMARK_SUMMARY.csv`: 쌍 실행 시간.
- `figures/`: Rrs 와 rrs 의 독립 I/Q/U 산포도.

## 산포도 세트

- `all600_*`: 600 셀 전체.
- `purewater150_*`: 순수해수, 전 SZA.
- `purewater84_sza_le40_*`: 순수해수, SZA ≤ 40°.
- `internal_raa_consistency_before_after.png`: 내부 수정 전후 진단.

세기(I) 그림은 MAPE 를 보인다. Q 와 U 그림은 0 근처에서 MAPE 가 불안정한 문제를 피하기 위해
`RMS[100 · (OCRT − ref) / ref_I]` 를 보인다.
