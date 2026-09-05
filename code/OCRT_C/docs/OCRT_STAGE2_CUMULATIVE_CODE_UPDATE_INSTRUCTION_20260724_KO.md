# OCRT–OSOAA 2단계 누적 코드 업데이트 작업지시서

문서 기준일: 2026-07-24  
적용 시작점: 세션 시작 시점 OCRT 소스 `OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22`  
적용 종료점: `OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix` + 2026-07-24 TAW 진단도구/검증계약  
본 문서는 마지막 수정분만을 기술하지 않는다. **세션 시작 소스에서 최종 수정본까지의 전체 누적 절차**를 기술한다.

## 1. 변경관리 원칙

1. 기준 소스 ZIP의 SHA-256은 다음과 같아야 한다.

```text
bb6ecca0a4bab5463832dae92416b5e7cd3d28eb03da150398c32625a89ce377
```

2. 기준선 재현 단계에서는 소스를 수정하지 않는다. 기존 72개 실행 단위와 600개 비교 셀을 먼저 재현하고, 그 결과를 변경 전 동결 기준으로 보존한다.
3. OCRT와 OSOAA의 입력 데이터 및 옵션은 완전히 일치시킨다. 불명확한 입력은 임의 선택하지 않는다.
4. RAA 보정과 U 부호 보정을 중복 적용하지 않는다.
5. OSOAA `Rrs(0+)`는 최종 전달규약 확정 전까지 두 branch를 병행 보존한다. level-27 수중 Fourier 장 + 별도 TWA를 잠정 기준으로 사용하고, level-26 water-minus-black 차감값은 독립 진단 branch로 유지한다. 어느 branch도 TWA/TAW 물리 불변량 검증 전에 유일한 정본으로 선언하지 않는다.
6. 코드 변경 전후 실행시간을 동일 조건에서 교차 측정한다. 평균 또는 중앙값이 5% 이상 증가하면 프로파일링을 수행한 뒤 승인한다.
7. 정합성 결과는 Rrs I/Q/U와 rrs I/Q/U를 각각 독립 산포도로 작성한다. Q/U는 MAPE를 주 지표로 사용하지 않고 기준 I 정규화 차이를 함께 제시한다.
8. 비교코드에 대한 상대 오차 감소만으로 물리코드를 수정하지 않는다. 단일 연산자 대조, contraction closure, reciprocity/energy/étendue 또는 독립 해석해 중 하나 이상이 수정 방향을 지지하여야 한다.

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

잠정 기준 OSOAA reference branch는 다음 순서로 생성한다.

```text
1. level-27 underwater Fourier field extraction
2. Snell-mapped underwater VZA reconstruction
3. separate TWA water-to-air transfer
4. Ed(0+) normalization for Rrs(0+)
5. Ed(0-) normalization for rrs(0-)
6. one and only one RAA/U mapping from section 4
```

첨부 `compare_entry.py`의 level-26 water-minus-black `Rrs(0+)`도 독립 진단 branch로 보존하며, 두 branch의 차이를 매 gate에서 기록한다.

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


## 14. 단계 10 — 첨부 계면 분석 재계산 및 범위 감사

### 14.1 원자료 단일화

첨부 패키지의 본문·그림·CSV가 서로 다른 지표 또는 분석 스냅샷을 혼합할 수 있으므로, 판정은 다음 두 CSV에서 직접 재계산한다.

```text
gate5_2026-07-22.csv
interface_testbed_2026-07-22.csv
```

재계산 결과와 그림은 다음 디렉터리에 보존한다.

```text
validation/stage2_taw_operator_audit_2026-07-24/recalculated/
```

순수해수 gate는 Rrs/rrs I/Q/U를 6개 독립 산포도로 작성한다. 계면 시험대는 I, Q/I, U/I를 각각 독립 산포도로 작성한다.

### 14.2 첨부 gate의 적용 범위

첨부 gate를 다음의 **국소 회귀시험**으로 정의한다.

```text
wavelength = 412, 443, 490, 555, 660 nm
SZA = 40 deg
VZA = 60 deg
RAA = 0,45,90,135 deg
wind = 3 m/s
pure water
aerosol = 0
```

이 gate의 결과를 SZA 80°, 다른 VZA, 865 nm, 풍속 0 또는 입자·에어로솔 조건으로 일반화하지 않는다.

### 14.3 지표 고정

각 결과 파일에는 지표 정의를 문자열로 기록한다.

```text
I: MAPE 또는 NMAE/max|I_ref| 중 사용한 식을 명시
Q/U: RMS[(OCRT-OSOAA)/OSOAA_I]
보조: NMAE/max|Q_ref|, NMAE/max|U_ref|
```

signed bias, MAE, RMSE를 혼용하지 않는다.

## 15. 단계 11 — OSOAA TAW Fourier 연산자 직접 감사

### 15.1 추가 도구

다음 도구를 추가한다.

```text
tools/stage2_taw_operator_audit.py
scripts/run_stage2_taw_operator_audit.sh
```

도구는 생산 solver 호출경로에 연결하지 않는다. `src/` 물리코드는 이 단계에서 수정하지 않는다.

### 15.2 입력

```text
OSOAA TAW:
  TAW-1.340-03.0-RadMU48-NB200-SZA40.000-TSZA28.665
angles:
  pure_water_b443_s40/Advanced_outputs/RAD_UsedAngles.txt
OCRT:
  surface_T_aw_coxmunk_fourier_kernel
wind = 3.0
n_water = 1.34
sigma_type = 1
q_convention = 1
n_phi = 1024
m = 0..4
active nodes = weight > 0, 48 nodes
```

### 15.3 OSOAA→OCRT 정규화

OSOAA Step 7과 OCRT contraction을 다음과 같이 맞춘다.

```text
K_OSOAA^m = 2 * TAW^m / (C_m * mu_water * mu_air)
C_0 = 2*pi
C_m = pi  for m>0
```

공개 OCRT RAA와 OSOAA azimuth의 180° 차이는 mode별 `(-1)^m`으로 한 번만 적용한다.

### 15.4 baseline signature

현재 수정 전 TAW 편광 잔차의 진단 signature를 다음과 같이 기록한다.

```text
M11: best-fit scale 0.999988..1.000003, max/maxref 0.011..0.091%
M12: OCRT/OSOAA ~1.010..1.011
M21: OCRT/OSOAA ~1.017..1.018
M23: sign/shape nearly identical
M32: nearly equal magnitude, opposite sign
M13/M31: not explained by one global scale/sign
```

이는 승인된 최종 상태가 아니라, 수정 전 국소화 baseline이다.

### 15.5 배제시험

다음을 함께 수행한다.

- n_phi 128/256/512/1024 convergence
- q_convention 0/1
- input/output rotation sine sign diagnostic variants
- zero-weight auxiliary node 포함/제외
- 기본 TAW와 SZA-specific TAW의 active submatrix identity

현재 결과에서 n_phi 128→1024 최대 차이는 0.000768% 이하이고, q_convention=0은 M12/M21을 약 12% 수준으로 악화시킨다.

## 16. 단계 12 — 해수 깊이 입력계약 보정

### 16.1 기존 상태 기록

기준 OSOAA raw run은 `SEA.Depth=1000 m`를 사용한다. OCRT 기존 순수해수 기준은 adaptive depth로 443 nm에서 약 `tau=20`, `z=1674.8 m`를 사용하였다. 두 조건은 사실상 반무한 수심으로 결과가 수렴하지만 입력값은 문자 그대로 동일하지 않다.

### 16.2 공식 parity gate

향후 공식 OCRT–OSOAA gate에서는 다음을 양쪽에 고정한다.

```text
physical water depth = 1000 m
black bottom
```

OCRT 현 CLI에서는 검증용으로 다음을 명시한다.

```bash
export OCRT_DEBUG=1
--debug-water-max-depth 1000
```

이 옵션은 비교 입력계약을 고정하는 검증용이며 운영 기본값으로 승격하지 않는다.

### 16.3 수렴 근거

pure water 443 nm, SZA40, VZA60, RAA90 sweep 결과:

```text
200 m  (tau 2.388): deep reference 대비 I 약 -0.7%
500 m  (tau 5.971): 약 -0.0005%
1000 m (tau 11.941): 약 -0.00002%
adaptive no-cap (tau 20): reference
```

따라서 200 m cap의 위험은 인정하되, 고정 `minimum tau=30`을 일반 규칙으로 도입하지 않는다. 매 IOP별 depth convergence로 결정한다.

## 17. 단계 13 — contraction closure

물리코드 수정 전 다음 세 입력 Stokes basis를 mode/node별로 주입한다.

```text
I basis = [1,0,0]^T
Q basis = [0,1,0]^T
U basis = [0,0,1]^T
```

OSOAA Step 7과 OCRT coupling contraction이 각 basis에서 생성하는 in-water Stokes를 직접 대조한다. 이후 실제 level-26 Rayleigh BOA Fourier field를 입력한다.

승인 조건:

1. 행렬 차이가 level-27 Q/U 차이를 수치적으로 재현한다.
2. residual이 다른 solver block에서 추가로 생성되지 않음을 확인한다.
3. node ordering, `mu*w`, `2/mu_water`, Fourier normalization을 각각 독립 검증한다.

## 18. 단계 14 — 물리 불변량 판정

다음 gate 중 최소 두 개가 동일 수정 방향을 지지할 때만 생산 코드를 변경한다.

- flat-interface analytical Mueller limit
- wind→0 continuity
- TAW/TWA reciprocity
- flux energy conservation
- Snell–Bouguer/étendue `L/n^2` invariance
- principal-plane U=0
- independent Jones-vector basis calculation

OSOAA를 무조건 정답으로 두지 않는다. OCRT 또는 OSOAA 어느 쪽도 수정 후보가 될 수 있다.

## 19. 단계 15 — 후속 생산 패치 승인 절차

### 19.1 최소 패치

불일치 성분과 물리 원인이 확정된 뒤 해당 함수·성분만 수정한다. 전역 U 부호반전 또는 회전행렬 전체 반전은 금지한다.

### 19.2 필수 재검증

```text
20-cell SZA40 interface gate
150-cell pure-water stress gate
66-cell SZA80 subset
600-cell full comparison
Rrs/rrs IQU six scatterplots per scope
mirror identity and TOA invariance
flat/wind-zero/reciprocity/energy gates
```

### 19.3 계산시간

생산 패치 전후 대표 조건을 각각 1회 예열 후 5회 교차 측정한다. 평균 또는 중앙값이 5% 이상 증가하면 profiling을 수행한다.

2026-07-24 진단 도구 추가는 생산 binary를 변경하지 않는다.

```text
production binary SHA-256:
41391518435c788c4d99f798c82fd39a016c8e4227ce220cedc2096f001838df
runtime impact on production solver: exactly zero (byte-identical binary)
standalone TAW audit elapsed: about 5.83 s for 48x48, m=0..4, n_phi=1024
```

## 20. 2026-07-24 승인 및 보류

### 20.1 승인

- 첨부 계면 국소화 가설을 후속 우선 분석축으로 채택
- OSOAA TAW decoder와 OCRT kernel 직접 감사 도구
- q convention/sigma/n_phi/depth의 manifest 의무화
- water depth 1000 m parity gate
- Rrs 두 extraction branch 병행 보존
- production physics code 무변경

### 20.2 보류

- TAW sine/reference-plane 불일치의 물리적 정답
- OCRT 또는 OSOAA 중 수정 대상 결정
- SZA80 잔차가 TAW 수정만으로 해소되는지 여부
- level-27+TWA와 level-26 subtraction 중 최종 Rrs 정본
- 입자 phase-function truncation 정합성
- aerosol 포함 coupled-field 정합성
