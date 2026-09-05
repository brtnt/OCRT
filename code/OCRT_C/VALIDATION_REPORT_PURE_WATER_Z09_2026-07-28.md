# OCRT C/Python 및 OSOAA 순수해수 IOP 표 교체 검증 보고서

작성일: 2026-07-28  
입력 handoff: `files (45).zip` / `FIX_pure_water_table_z09_2026-07-27.md`

## 1. 적용 결론

첨부 기준자료를 절대 기준으로 사용하여 세 패키지의 순수해수 광특성표를 교체했다.

- OCRT C: `inputs/water_iop/water_coef_z09_1nm.txt`
- OCRT Python: `data/water_iop/water_coef_z09_1nm.txt`
- OSOAA: `fic/OSOAA_SEA_MOL_COEFFS_JUNE_2013.txt`

세 코드가 사용하는 수치 배열은 완전히 같다.

- 범위: 200–2449 nm
- 간격: 1 nm
- 행수: 2,250
- 열: wavelength, `a_w`, total `b_w`
- `bb_w`: OCRT 내부에서 정확히 `0.5*b_w`
- 생성 규약: 원 0.1 nm 자료의 `[lambda-0.5, lambda+0.5]` 사다리꼴 띠평균
- 산란 상태값: T=20 °C, S=38.4 g kg⁻¹

RT 알고리즘, Fourier 재구성, native coupled LUT orchestration은 변경하지 않았다.

## 2. 기존 표의 문제와 변경 규모

기존 OCRT 표는 400–900 nm의 501행으로 제한됐고, 일부 앵커 사이 흡수계수가 잘못 재구성돼 있었다. 새 표와의 대표 차이는 다음과 같다.

| 파장 | a_w 변화 | b_w 변화 | 순수수 Rrs_I 변화* |
|---:|---:|---:|---:|
| 412 nm | +0.067% | -11.128% | -8.635% |
| 443 nm | -0.027% | -11.091% | -9.552% |
| 510 nm | +45.145% | -10.567% | -37.578% |
| 555 nm | +0.058% | -10.032% | -10.001% |
| 620 nm | +30.689% | -9.027% | -30.355% |
| 680 nm | -15.159% | -8.242% | +8.149% |
| 745 nm | +113.101% | -6.617% | -56.176% |
| 865 nm | +0.002% | -3.873% | -3.876% |
| 970 nm | +732.310% | -29.480% | -91.527% |
| 1100 nm | +230.921% | -57.889% | -87.274% |

\* 동일한 축소 순수수 RT 진단 조건에서 수정 전/후 C 결과 비교. 이는 의도된 물리 입력 변경이며 회귀 오류가 아니다.

380 nm도 이제 표 범위 안에 포함된다. 970–1100 nm에서 기존 900 nm 끝값 유지가 제거됐다.

## 3. 데이터 무결성

| 검증 | 결과 |
|---|---:|
| OCRT C 표 행수/범위 | 2,250행, 200–2449 nm PASS |
| OCRT Python 표 행수/범위 | 2,250행, 200–2449 nm PASS |
| OSOAA 표 행수/범위 | 2,250행, 200–2449 nm PASS |
| C ↔ Python 수치 배열 최대 절대차 | 0 |
| C ↔ OSOAA 수치 배열 최대 절대차 | 0 |
| C ↔ Python 대표 22파장 a_w/b_w/bb_w | 최대 절대차 0 |
| `bb_w=0.5*b_w` | PASS |

C runtime 로그는 `loaded 2250 rows ... (lambda 200.00..2449.00 nm)`를 확인했다.

## 4. LUT 동시계산 비회귀

### C

새 표가 실제로 사용되는 mixed-water 490 nm, 3 VZA × 4 RAA 조건에서:

- native coupled LUT: `water_cold=1`, `water_views=3`, `cells=12`
- native와 authoritative cell replay CSV: 모든 열 byte-identical
- pure-water 555 nm native/cell replay 결과도 수치상 완전 동일
- fixed-IOP 값 커널 native LUT 회귀시험: PASS

주의: pure-water-only full-grid가 `native coupled LUT unavailable rc=-5`로 legacy replay에 내려가는 기존 동작은 수정 전 기준본에도 동일했다. 이번 데이터 패치로 새로 발생한 회귀가 아니다. Mixed/aerosol coupled 경로의 native solve-once 구조는 정상 유지됐다.

### Python

- native coupled LUT vs authoritative cell replay: exact equality
- 진단 카운터: atmosphere pass-1=1, water=1, atmosphere pass-2=1, cell solver=0
- 값 커널 native LUT 보존시험 포함 전체 pytest PASS

## 5. 성능

동일한 fixed-bulk IOP를 사용해 순수 데이터 변경의 코드 오버헤드를 분리했다.

| C full-grid, 3 VZA × 40 RAA | 중앙값 |
|---|---:|
| 수정 전 | 0.415321 s |
| 수정 후 | 0.415175 s |
| 변화 | -0.035% |

출력 CSV는 byte-identical이다. 즉 더 큰 표를 배포해도 측정 가능한 실행시간 저하는 없다.

## 6. OSOAA 검증

### 6.1 빌드 및 smoke

- 수정 OSOAA source rebuild: PASS
- 새 순수수 표 unit test: PASS
- 새 표 기반 active smoke: 11/11 PASS, relative error < 1e-9
- 기존 P6 smoke 기준은 물리 입력 변경 때문에 평균 약 -0.602% 이동했다.
- 과거 `reference_values_P6.csv`는 역사자료로 보존하고, active 기준은 `reference_values_Z09_20260728.csv`로 교체했다.

### 6.2 OCRT–OSOAA 교차검증

Red-clay TSM=5 g m⁻³, 555 nm, SZA=30°, RAA=90°, water-side VZA corresponding to air VZA=60°:

| 항목 | OCRT | OSOAA | 차이 |
|---|---:|---:|---:|
| DoLP | 0.08863865 | 0.08731997 | +1.51% |
| Q/I | -0.024290 | -0.023709 | — |
| |U|/I | 0.085239 | 0.084040 | — |

기존 m-sign 허용기준 ±2% 안이다.

CDOM aDOM440=0.1, SZA=40°, VZA={0,30,60}의 TOA I 교차검증:

| 파장 | TOA I MAPE | rrs I MAPE |
|---:|---:|---:|
| 443 nm | 0.082% | 0.868% |
| 555 nm | 0.204% | 0.446% |

세부 옵션 차이를 포함해도 모두 1% 이내이며, OCRT와 OSOAA가 같은 순수수 표를 읽는 상태를 확인했다.

## 7. 갱신한 golden/reference

의도된 데이터 변경 때문에 다음 active 기준을 새 표로 재생성했다.

- C water internal-reflection m-sign canonical anchors
- C mixed-water IOP/Kd full-grid anchors
- C coupling clamp `n_mu_water={48,64,96}` full-grid references
- OSOAA active smoke reference

이전 coupling reference는 `reference_pre_pure_water_z09_20260728/`에 보존했다. 그 밖의 과거 날짜 validation 자료는 provenance 용도이며, 순수수 IOP에 의존하면 current golden으로 사용하면 안 된다. 각 패키지의 `VALIDATION_REFERENCE_POLICY_2026-07-28.md`에 이 정책을 기록했다.

## 8. 회귀시험

### C

통과한 주요 항목:

- release build and version
- pure-water Z09 table
- value-kernel native LUT
- water internal-reflection m-sign
- RAA/direct-glint convention
- Mie P11/P12/P33 wavelength PCHIP
- external bottom-source bounds/shape guard
- Stage-2 interface FIX1–FIX3
- diffuse-top U closure
- exact surface pole limit
- Stage-2 water RAA
- EAP scattering-disabled production contract
- organic NIR gate
- IOP/Kd full-grid CSV
- coupling n_mu_water 48/64/96 and single/full parity

### Python

- compileall: PASS
- full pytest: 26/26 PASS
- targeted pure-water/native-LUT/value-kernel/m-sign: 7/7 PASS

## 9. 남은 제한

순수수 `a_w/b_w` 본표는 2449 nm까지 확장됐지만, 별도 온도보정 `psi_T` 표는 400–900 nm에 머문다. 따라서 20 °C가 아닌 수온으로 400 nm 미만 또는 900 nm 초과를 계산하려면 향후 `psi_T` 확장이 필요하다. 20 °C 검증에서는 보정항이 0이므로 영향이 없다.
