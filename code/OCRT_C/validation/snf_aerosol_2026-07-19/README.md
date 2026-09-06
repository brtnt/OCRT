# SnF aerosol validation archive — 2026-07-19

This directory records the controlled conversion and integration of the 16 SnF atmospheric
aerosol inputs.

## Pass summary

| check | result |
|---|---:|
| baseline regeneration: T50/C50/M80C/O99 | PASS |
| baseline RT identity | 24/24 exact |
| deterministic regeneration | 16/16 byte-identical |
| structural / finite / SSA checks | 16/16 pass |
| Mueller physicality | 16/16 pass |
| OCRT aerosol-only runtime smoke | 48/48 pass |

## Files

- `baseline_reproduction.csv`: direct numeric comparison of the regenerated and the pre-existing
  high-resolution baseline tables.
- `baseline_rt_identity_24cases.csv`: I/Q/U and AOD equality of the baseline tables in two
  pressure regimes and three wavelengths.
- `generated_quality.csv`: structure, normalization, g-consistency and physicality metrics of all
  16 generated models.
- `legacy_replacement.csv`: phase and spectral differences between the displaced 83-angle
  M50C/M95C/M98C and the canonical regenerated lineage.
- `legacy83_vs_canonical_rt_18cases.csv`: representative I/Q/U impact of the three lineage
  replacements.
- `runtime_smoke_48cases.csv`: 16 models × 412/550/860 nm, aerosol only over a black surface,
  AOD555 = 0.15.
- `model_integration_matrix.csv`: final runtime action and qualification status.
- `generator_reproducibility.txt`, `generated_file_sha256.txt`: deterministic generation record.
- `validation.json`: detailed machine-readable phase-file checks.
- `non_pssa_regressions.log`: rerun of the CCRR Chl, OCRT organic Chl, Ahn TSM and TSM
  phase-cache regressions after the SnF integration.

## Interpretation boundary

These files establish generator lineage, table integrity, physical sanity and OCRT runtime
compatibility. They do not promote the newly added models into the external benchmark subset
compared against reference codes; that subset remains T50/C50/M80C until a separate full
reference-code campaign is run.

---

# SnF 에어로졸 검증 아카이브 — 2026-07-19

이 디렉터리는 SnF 대기 에어로졸 입력 16종의 통제된 변환과 통합을 기록한다.

## 통과 요약

| 검사 | 결과 |
|---|---:|
| 기준 재생성: T50/C50/M80C/O99 | PASS |
| 기준 RT 동일성 | 24/24 정확 |
| 결정적 재생성 | 16/16 바이트 동일 |
| 구조/유한값/SSA 검사 | 16/16 통과 |
| Mueller 물리성 | 16/16 통과 |
| OCRT 에어로졸 단독 런타임 스모크 | 48/48 통과 |

## 파일

- `baseline_reproduction.csv`: 재생성 표와 기존 고해상도 기준 표의 직접 수치 비교.
- `baseline_rt_identity_24cases.csv`: 두 기압 조건·세 파장에서 기준 표의 I/Q/U 와 AOD 동일성.
- `generated_quality.csv`: 생성 모델 16종 전체의 구조·정규화·g 일관성·물리성 지표.
- `legacy_replacement.csv`: 대체된 83각 M50C/M95C/M98C 와 표준 재생성 계보의 위상·분광 차이.
- `legacy83_vs_canonical_rt_18cases.csv`: 계보 교체 3건의 대표 I/Q/U 영향.
- `runtime_smoke_48cases.csv`: 16 모델 × 412/550/860 nm, 흑색 표면 위 에어로졸 단독, AOD555 = 0.15.
- `model_integration_matrix.csv`: 최종 런타임 조치와 자격 상태.
- `generator_reproducibility.txt`, `generated_file_sha256.txt`: 결정적 생성 기록.
- `validation.json`: 기계 판독용 상세 위상 파일 검사.
- `non_pssa_regressions.log`: SnF 통합 후 CCRR Chl, OCRT 유기물 Chl, Ahn TSM, TSM 위상 캐시 회귀 재실행.

## 해석 범위

이 파일들은 생성기 계보, 표 무결성, 물리적 타당성, OCRT 런타임 호환성을 확립한다. 새로 추가된 모델을 참조
코드와 비교하는 외부 벤치마크 부분집합에 올리지는 않는다. 그 부분집합은 별도의 전체 참조 코드 캠페인이
있을 때까지 T50/C50/M80C 로 유지한다.
