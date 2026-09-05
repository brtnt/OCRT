# OCRT–OSOAA 2단계 누적 코드 업데이트 작업지시서

문서 기준일: 2026-07-23  
적용 시작점: 세션 시작 시점 OCRT 소스 `OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22`  
적용 종료점: `OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix`  
본 문서는 마지막 수정분만을 기술하지 않는다. **세션 시작 소스에서 최종 수정본까지의 전체 누적 절차**를 기술한다.

## 1. 변경관리 원칙

1. 기준 소스 ZIP의 SHA-256은 다음과 같아야 한다.

```text
bb6ecca0a4bab5463832dae92416b5e7cd3d28eb03da150398c32625a89ce377
```

2. 기준선 재현 단계에서는 소스를 수정하지 않는다. 기존 72개 실행 단위와 600개 비교 셀을 먼저 재현하고, 그 결과를 변경 전 동결 기준으로 보존한다.
3. OCRT와 OSOAA의 입력 데이터 및 옵션은 완전히 일치시킨다. 불명확한 입력은 임의 선택하지 않는다.
4. RAA 보정과 U 부호 보정을 중복 적용하지 않는다.
5. 공식 OSOAA `Rrs(0+)`는 level-27 수중 Fourier 장과 별도 TWA 전달 계산으로 산출한다. level-26 water-minus-black 차감값을 공식 기준으로 사용하지 않는다.
6. 코드 변경 전후 실행시간을 동일 조건에서 교차 측정한다. 평균 또는 중앙값이 5% 이상 증가하면 프로파일링을 수행한 뒤 승인한다.
7. 정합성 결과는 Rrs I/Q/U와 rrs I/Q/U를 각각 독립 산포도로 작성한다. Q/U는 MAPE를 주 지표로 사용하지 않고 기준 I 정규화 차이를 함께 제시한다.

## 2. 단계 0 — 세션 시작 기준선 동결

### 2.1 기준 소스 준비

```bash
unzip "OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22 (1)(1).zip"
cd OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22
```

원본은 읽기 전용 보관본으로 유지하고, 수정 작업은 복사본에서 수행한다.

```bash
cp -a OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22       OCRT_v1.2_STAGE2_WATER_RAA_OUTPUT_FIX_WORK
cd OCRT_v1.2_STAGE2_WATER_RAA_OUTPUT_FIX_WORK
```

### 2.2 기준 바이너리 빌드

```bash
./scripts/build_release_v1.2.sh build/ocrt_stage2_base
./build/ocrt_stage2_base --version
sha256sum build/ocrt_stage2_base
```

기준 바이너리의 확인값은 다음과 같다.

```text
1ff23e3fd53dd38faf6d7069efef2b2e0cbe02acafe8979010c450360dfe8b28
```

### 2.3 기준선 재현 조건

다음 입력계약을 고정한다.

```text
Atmosphere:
  Rayleigh on
  pressure = 1013.25 hPa
  aerosol = 0
  gas columns H2O/O3/NO2/O2/CO2/CH4 = 0
  n_mu = 48
  n_layers = 40
  m_max = 2
  sos_max_orders = 100
  water max orders = 100 via OCRT_DEBUG=1 and --debug-water-max-orders 100

Ocean:
  pure water and CDOM a440 = 0.01, 0.1, 1.0 m^-1
  CDOM slope = 0.014 nm^-1
  Chl = 0
  TSM = 0
  n_water = 1.34
  temperature = 20 degC
  salinity = 38.4 PSU
  seawater depolarization = 0.039
  black bottom
  water shared grid
  wind speed = 3 m/s
  direct sunglint decoupled

Spectral/geometric grid:
  wavelength = 412, 443, 490, 555, 660, 865 nm
  SZA = 0, 40, 80 deg
  VZA = 0, 30, 60 deg
  RAA = 0, 45, 90, 135, 180 deg
```

72개 실행 명령은 최종 패키지의 다음 파일에 수록한다.

```text
validation/stage2_water_raa_output_fix_2026-07-23/RUN_MANIFEST_72CASES.csv
```

기준선의 승인 조건은 72/72 실행 수렴과 600/600 비교 셀 생성이다.

## 3. 단계 1 — 수중 공개 RAA 출력 규약 수정

### 3.1 `src/rt_water_rt.c`

#### 3.1.1 주 수중 Fourier 재구성

다음 변경을 적용한다.

```diff
- double dphi_rad = rt_raa_to_water_view_phi(raa_recon_deg);
+ double dphi_rad = rt_raa_to_atm_fourier_phi(raa_recon_deg);

- double U_diffuse_view = -rt_solver_reconstruct_phi_sin(...);
+ double U_diffuse_view =  rt_solver_reconstruct_phi_sin(...);
```

목적은 `Rrs(0+)`와 `rrs(0-)`의 보고 방위각을 대기 출력과 동일한 공개 OCRT RAA로 통일하는 것이다.

#### 3.1.2 정확 단일산란 보정

정확 단일산란 전파벡터 기하는 출력 재구성과 구분한다.

```diff
- double phv = rt_raa_to_water_view_phi(raa_deg);
+ double phv = rt_raa_to_water_scatter_phi(raa_deg);
```

이 지점의 `pi-RAA`는 물리적 입사·출사 벡터 구성에 필요한 국소 좌표이므로 제거하지 않는다.

정확 단일산란 U의 내부 Fourier 재구성에서는 추가 음수를 제거한다.

```diff
- double Uss_int = -rt_solver_reconstruct_phi_sin(...);
+ double Uss_int =  rt_solver_reconstruct_phi_sin(...);
```

#### 3.1.3 TWA 노드별 수면 위 U 재구성

동일 패턴의 세 실행문을 모두 수정한다.

```diff
- double Lw_U = -rt_solver_reconstruct_phi_sin(...) * f_scale;
+ double Lw_U =  rt_solver_reconstruct_phi_sin(...) * f_scale;
```

세 지점 중 하나라도 남으면 branch 또는 보간 경로에 따라 U 부호가 다시 갈릴 수 있으므로 일괄 적용한다.

#### 3.1.4 풍속 양수 rough-interface 출력

```diff
- const double phi_v = rt_raa_to_water_view_phi(raa_deg);
+ const double phi_v = rt_raa_to_atm_fourier_phi(raa_deg);
```

이 변경은 풍속 양수 Cox–Munk/TWA 분기의 `Rrs(0+)` 공개 방위각을 바로잡는다.

#### 3.1.5 결합 Fourier U 재구성

```diff
- LU += -f * Ua[mm] * s;
+ LU +=  f * Ua[mm] * s;
```

### 3.2 `src/rt_solver.c`

D3 외부 하향장 보정 경로에서 다음을 적용한다.

```diff
- double dphi = rt_raa_to_water_view_phi(cs->raa_deg);
+ double dphi = rt_raa_to_atm_fourier_phi(cs->raa_deg);

- double dvU = -rt_solver_reconstruct_phi_sin(...);
+ double dvU =  rt_solver_reconstruct_phi_sin(...);
```

### 3.3 변경 지점 계수

실제 실행문 변경은 총 10개이다.

- RAA 변환: 3개
- U 부호: 7개
- 논리 패턴으로 묶을 경우: 8종

검토 문서에 “아홉 곳”이라고 적지 않는다.

## 4. 단계 2 — 규약 helper와 주석의 강제 분리

`src/rt_raa_convention.h`에서 다음을 적용한다.

```diff
- rt_raa_to_water_view_phi
+ rt_raa_to_water_scatter_phi
```

새 helper의 주석에는 다음 제한을 명시한다.

```text
이 helper는 정확 단일산란에서 실제 전파벡터를 구성할 때만 사용한다.
보고되는 Rrs(0+) 및 rrs(0-) Fourier 출력 재구성에는 사용하지 않는다.
```

OSOAA 대응 주석은 다음 두 동등한 방식으로 갱신한다.

```text
A. Phi_OSOAA = (180 - RAA_OCRT) mod 360, then U_OSOAA *= -1
B. Phi_OSOAA = (RAA_OCRT + 180) mod 360, no U reversal
```

정확히 하나만 사용한다.

소스와 시험 디렉터리에서 legacy helper의 실행 코드 잔존을 금지한다.

```bash
grep -RIn --include='*.c'   'rt_raa_to_water_view_phi' src tests
```

결과가 없어야 한다.

## 5. 단계 3 — 회귀시험 추가

### 5.1 단위시험 확장

`tests/test_raa_convention.c`에 다음을 추가한다.

- `rt_fourier.h` 포함
- `rt_raa_to_water_scatter_phi()`의 0°/180° 국소 기하 검사
- 합성 cosine 모드로 공개 RAA의 I/Q branch 검사
- 합성 sine 모드로 U의 양의 sine 재구성 검사

### 5.2 소스 감사 스크립트

다음 파일을 추가한다.

```text
scripts/test_stage2_water_raa_fix.sh
```

실행:

```bash
bash scripts/test_stage2_water_raa_fix.sh
```

기대 출력:

```text
PASS: OCRT public RAA convention helpers
PASS: stage-2 water RAA source audit
```

### 5.3 기준코드 비의존 실행 회귀시험

다음 파일을 추가한다.

```text
tests/test_water_raa_output_regression.py
```

검사 항목은 다음과 같다.

1. 순수해수 443 nm, SZA 40°, VZA 60°, RAA 0/45/90/135/180°에서 `rrs_I / TOA_water_increment_I` 상대 범위가 8% 미만
2. RAA 0°와 180° 주평면에서 `Rrs_U/I`, `rrs_U/I`가 \(10^{-10}\) 미만
3. full-grid와 단일 기하 출력이 허용오차 내 일치
4. 풍속 0 분기 실행 및 수렴
5. 풍속 0 주평면 U 검사

실행:

```bash
python3 tests/test_water_raa_output_regression.py build/ocrt_stage2_final
```

기대 결과:

```text
PASS: water RAA output regression
```

## 6. 단계 4 — 버전과 문서

`src/main.c`의 버전을 다음과 같이 변경한다.

```c
#define V2_VERSION "OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix"
```

`README.md`의 Current validated release를 동일하게 변경하고, 다음 세 항목을 추가한다.

- water-output RAA convention correction
- public output와 local single-scatter geometry 분리
- RAA/full-grid/wind-zero 회귀시험 추가

다음 문서를 포함한다.

```text
docs/STAGE2_WATER_RAA_OUTPUT_FIX_REVIEW_2026-07-23_KO.md
docs/OCRT_STAGE2_CUMULATIVE_CODE_UPDATE_INSTRUCTION_2026-07-23_KO.md
```

## 7. 단계 5 — 최종 빌드

```bash
./scripts/build_release_v1.2.sh build/ocrt_stage2_final
./build/ocrt_stage2_final --version
sha256sum build/ocrt_stage2_final
```

기대 버전:

```text
ocrt OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix
```

본 작업에서 생성한 최종 바이너리 SHA-256:

```text
41391518435c788c4d99f798c82fd39a016c8e4227ce220cedc2096f001838df
```

## 8. 단계 6 — 수치 검증

### 8.1 72개 실행 재실행

`RUN_MANIFEST_72CASES.csv`의 각 명령에서 다음 placeholder를 치환한다.

```text
<OCRT_BINARY> = build/ocrt_stage2_final
<OUTPUT_CSV>  = 실행별 full-grid CSV 경로
```

환경은 다음과 같이 고정한다.

```bash
export OCRT_ADVANCED=1
export OCRT_DEBUG=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C
```

승인 조건:

```text
72/72 runs complete
72/72 water_converged = 1
600/600 comparison cells present
```

### 8.2 거울변환 항등식

모든 비교 셀에 대해 다음을 검사한다.

```text
new I(r) == old I(180-r)
new Q(r) == old Q(180-r)
new U(r) == -old U(180-r)
```

검사 대상은 Rrs와 rrs의 I/Q/U 전부이다. 최대 절대차 기대값은 0이다.

### 8.3 비영향성

동일 RAA에서 다음 값은 기준선과 동일하여야 한다.

```text
TOA_rho_I
TOA_rho_Q
TOA_rho_U
Ed0plus_air
Ed0minus_water
```

최대 절대차 기대값은 0이다.

### 8.4 OSOAA 공식 비교

공식 OSOAA reference 생성은 다음 순서로 고정한다.

```text
1. level-27 underwater Fourier field extraction
2. Snell-mapped underwater VZA reconstruction
3. separate TWA water-to-air transfer
4. Ed(0+) normalization for Rrs(0+)
5. Ed(0-) normalization for rrs(0-)
6. one and only one RAA/U mapping from section 4
```

첨부 `compare_entry.py`의 level-26 water-minus-black `Rrs(0+)`는 진단용 보조값으로만 보존한다.

## 9. 단계 7 — 산포도 작성 규칙

각 비교 범위에 대해 다음 6개 파일을 각각 생성한다.

```text
Rrs0plus_I_scatter.png
Rrs0plus_Q_scatter.png
Rrs0plus_U_scatter.png
rrs0minus_I_scatter.png
rrs0minus_Q_scatter.png
rrs0minus_U_scatter.png
```

본 패키지는 세 범위에 대해 총 18개 그림을 제공한다.

```text
all600_*
purewater150_*
purewater84_sza_le40_*
```

그림 작성 규칙:

- x축은 OSOAA, y축은 OCRT
- 1:1 선 포함
- SZA를 서로 다른 marker로 표시
- I에는 MAPE 병기
- Q/U에는 `RMS[100*(OCRT-OSOAA)/OSOAA_I]` 병기
- Q/U 축은 0을 중심으로 대칭
- 표만으로 정합성을 판정하지 않음

## 10. 단계 8 — 계산시간 검증

다음 대표 조건을 기준과 수정본에서 각각 1회 예열 후 5회 교차 실행한다.

```text
pure water
443 nm
SZA 40°
VZA grid 0/30/60°
RAA step 45°
n_mu 48
m_max 2
wind 3 m/s
gas off
```

본 검증 결과:

```text
baseline mean  = 4.848507 s
updated mean   = 4.848084 s
mean change    = -0.0087%

baseline median = 4.827955 s
updated median  = 4.828572 s
median change   = +0.0128%
```

승인 규칙:

- 평균·중앙값 모두 5% 미만 증가: 성능 PASS
- 5% 이상 증가: `perf`, `gprof` 또는 sampling profiler로 hotspot 확인
- 새 반복문·할당·적분이 추가된 경우에는 증가율과 무관하게 hotspot 검토

본 수정은 단순 방위각 argument와 U 부호의 치환이며 실행시간 증가는 관측되지 않았다.

## 11. 단계 9 — 승인 및 보류 항목

### 11.1 승인

- 수중 공개 RAA와 대기 공개 RAA의 통일
- 추가 U 부호반전 제거
- 국소 단일산란 `pi-RAA` helper 분리
- 풍속 0 및 풍속 양수 분기 회귀시험
- full-grid/단일 기하 일치
- TOA 비영향성
- 실행시간 비증가

### 11.2 보류

- SZA 80° 편광 잔차의 물리 원인
- Chl/TSM 위상함수 포함 정합성
- aerosol 포함 coupled-field 정합성
- 누락된 Tier-0 phase LUT 자산을 사용하는 정규 Tier-0 실행

## 12. 롤백

수정 패치를 제거할 때는 세션 시작 기준 소스로 되돌린다.

```bash
rm -rf OCRT_v1.2_STAGE2_WATER_RAA_OUTPUT_FIX_WORK
cp -a OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22       OCRT_v1.2_STAGE2_WATER_RAA_OUTPUT_FIX_WORK
```

부분 롤백은 금지한다. RAA argument만 되돌리거나 U 부호만 되돌리면 동일 결함이 다른 형태로 재발한다.

## 13. 납품물 확인

최종 납품에는 다음을 포함한다.

```text
updated full source tree
final binary
cumulative patch
cumulative work instruction
review and validation report
600-cell corrected matched CSV
72-case run manifest
runtime raw/summary CSV
attached-method vs approved-TWA audit CSV
mirror identity/TOA invariance CSV
Rrs/rrs IQU scatterplots
SHA-256 manifest
```
