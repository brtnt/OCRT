# OCRT–OSOAA 2단계 누적 코드 업데이트 작업지시서 — 개정 6

문서 기준일: 2026-07-26  
적용 시작점: 세션 시작 시점 OCRT 소스 `OCRT_v1.2_RAA_CONVENTION_COMMENT_AUDIT_2026-07-22`  
적용 종료점: `OCRT-v1.2-2026-07-26-KST-stage2-eap-generator-phyto-scattering-disabled`  
본 문서는 마지막 수정분만을 기술하지 않는다. **세션 시작 소스에서 최종 수정본까지의 전체 누적 절차**를 기술한다.

2026-07-24의 FIX1-only TAW 회전 후보와 FIX123 시점의 diffuse-top U 보류 결론은 역사적 checkpoint로 보존한다. 본 개정본의 단계 21 이후가 FIX1-only checkpoint를 대체하고, 단계 25 이후가 FIX123 시점의 보류 결론을 대체한다. 단계 31 이후는 정확한 수면 극한 회전 및 재현 빌드 규약을 추가한다. 단계 39 이후는 water-to-air coupling 보간의 64절점 silent clamp를 제거한다. 최종 생산 범위는 water RAA 출력 수정 + FIX1+FIX2+FIX3 + diffuse-top incoming-U column closure + rough-interface exact-pole limit + coupling full-grid interpolation 보존이다 + EAP 17종 phase generator 통합 + constituent phytoplankton scattering 비활성화까지 포함한다.

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
9. 릴리스 빌드는 `-march=native`를 사용하지 않는다. 기본 ISA는 `cascadelake`로 고정하고, 다른 ISA를 사용할 때에는 `OCRT_MARCH` 값과 새 바이너리 해시를 별도 기록한다. 동일 해시 주장은 동일 compiler/linker/toolchain에 한정한다.
10. 비트 동일 CSV 검사는 원시 바이트 또는 Python 내장 `float()` 기반 binary64 비교를 사용한다. `pandas`로 이미 파싱된 값을 비트 동일 판정의 근거로 사용하지 않는다.

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

## 20. 2026-07-24 1차 진단 시점 승인 및 보류 — 역사적 checkpoint

### 20.1 승인

- 첨부 계면 국소화 가설을 후속 우선 분석축으로 채택
- OSOAA TAW decoder와 OCRT kernel 직접 감사 도구
- q convention/sigma/n_phi/depth의 manifest 의무화
- water depth 1000 m parity gate
- Rrs 두 extraction branch 병행 보존
- 해당 checkpoint의 production physics code 무변경(FIX1+FIX2+FIX3 판정 전)

### 20.2 보류

- TAW sine/reference-plane 불일치의 물리적 정답
- OCRT 또는 OSOAA 중 수정 대상 결정
- SZA80 잔차가 TAW 수정만으로 해소되는지 여부
- level-27+TWA와 level-26 subtraction 중 최종 Rrs 정본
- 입자 phase-function truncation 정합성
- aerosol 포함 coupled-field 정합성

## 21. 단계 16 — 추가 결함 패키지 검토와 원자적 변경 범위 확정

### 21.1 첨부 주장 분해

추가 분석 패키지의 세 수정은 다음과 같이 분리한다.

```text
FIX1: surface_T_aw_coxmunk_trig에서 cphi_r=-cphi, sphi_r=-sphi
FIX2: 동일 회전 블록에서 mu_i_s=-mu_i, mu_o_s=-mu_o
FIX3: surface_T_aw_coxmunk_fourier_kernel 정규화 후 m>0의 M13/M23 부호 반전
```

세 수정은 독립적으로 선택하지 않는다. 다음 이유로 하나의 원자적 변경으로 적용한다.

1. FIX1은 실제 광자 진행방향에 대응하는 상대방위각의 π 이동이다.
2. FIX2는 공기→물 투과광의 입사·출사 방향여현을 signed propagation cosine으로 변환한다.
3. FIX3은 각도공간 Mueller 행렬이 아니라 sine-Fourier 저장·축약 convention을 공통화한다.
4. FIX2 없이 FIX3을 적용하면 이미 잘못된 각도공간 회전에 저장 부호만 덧씌우게 된다.
5. FIX2를 적용한 뒤 FIX3을 생략하면 M31/M32는 교정되지만 m>0 입사-U column contraction이 OSOAA Step 7과 달라진다.

### 21.2 기존 FIX1-only 후보의 지위

다음 버전은 최종 생산본이 아니라 중간 진단본으로 분류한다.

```text
OCRT-v1.2-2026-07-24-KST-stage2-taw-rotation-fix
```

FIX1-only 결과와 바이너리는 회귀 비교를 위해 보존하되 신규 생산 실행에는 사용하지 않는다.

## 22. 단계 17 — FIX1+FIX2 물리 회전식 적용

### 22.1 수정 파일

```text
src/shared/surface.c
```

### 22.2 `surface_T_aw_coxmunk_trig()` 수정

기존 FIX1-only 블록을 다음과 같이 확장한다.

```c
const double cphi_r = -cphi;
const double sphi_r = -sphi;
const double mu_i_s = -mu_i;
const double mu_o_s = -mu_o;

double cosPsi = mu_i_s * mu_o_s + si * so * cphi_r;
cosPsi = fmax(-1.0, fmin(1.0, cosPsi));

double ci1 = (mu_o_s - mu_i_s * cosPsi) / (dnm * s1_safe);
double si1 = so * sphi_r / dnm;

double ci2 = (mu_i_s - mu_o_s * cosPsi) / (dnm * s2_safe);
double si2 = si * sphi_r / dnm;
```

미시면 법선/Fresnel 입사각 계산에 사용하는 away-from-interface 벡터는 변경하지 않는다. signed cosine은 Stokes 회전 블록에만 적용한다.

### 22.3 독립 물리 근거

공기→물 진행벡터는 다음과 같다.

\[
\mathbf{k}_i=(-\sin\theta_i,0,-\mu_i),
\qquad
\mathbf{k}_o=(\sin\theta_o\cos\phi,\sin\theta_o\sin\phi,-\mu_o).
\]

따라서

\[
\mathbf{k}_i\cdot\mathbf{k}_o
=\mu_i\mu_o-\sin\theta_i\sin\theta_o\cos\phi.
\]

코드의 `mu_i_s`, `mu_o_s`, `cphi_r` 조합은 이 내적을 직접 재현해야 한다.

## 23. 단계 18 — FIX3 Fourier 입사-U 열 저장

### 23.1 적용 위치

`surface_T_aw_coxmunk_fourier_kernel()`에서 azimuth 적분과 mode 정규화를 완료한 직후 적용한다.

```c
for (int kl = 0; kl < 9; kl++) T_pair[kl] *= norm * dphi;

if (m > 0) {
    T_pair[0*3 + 2] = -T_pair[0*3 + 2];
    T_pair[1*3 + 2] = -T_pair[1*3 + 2];
}
```

### 23.2 금지사항

다음을 적용하지 않는다.

```text
M31/M32/M33 U-output row 전체 반전
M33 반전
m=0 반전
downstream contraction에서 동일 부호를 다시 적용
```

FIX3의 의미는 `m>0` sine-Fourier field에서 **incoming U column**에 필요한 음수를 저장 operator에 한 번 접어 넣는 것이다. 이후 contraction은 일반 3×3 행렬-벡터 곱으로 유지한다.

### 23.3 OSOAA Step 7 계약

OSOAA TAW contraction은 다음 구조로 해석한다.

```text
Iout = T11*Iin + T12*Qin + T13*Uin
Qout = T21*Iin + T22*Qin + T23*Uin
Uout = T31*Iin + T32*Qin + T33*Uin
```

따라서 저장 규약 보정은 입력 column 기준으로 판정한다.

## 24. 단계 19 — FIX123 독립 회귀시험

### 24.1 신규 단위시험

다음 파일을 추가한다.

```text
tests/test_surface_aw_interface_conventions.c
```

시험 1: outgoing-pole 해석 극한

```text
air incidence = 40 deg
water outgoing zenith = 1e-4 deg
physical azimuth offset = 30 deg
wind = 3 m/s
n_water = 1.34
```

다음 정규화 성분을 독립 Fresnel 해석식과 비교한다.

```text
M21/M11, M22/M11, M23/M11,
M31/M11, M32/M11, M33/M11
```

허용오차:

```text
2e-6 absolute
```

시험 2: Fourier incoming-U column

1. `surface_T_aw_coxmunk_trig_test()`를 `n_phi=4096`, `m=1`로 직접 적분한다.
2. `surface_T_aw_coxmunk_fourier_kernel()`의 저장값과 비교한다.
3. M13/M23만 direct integral의 음수이고 나머지 7성분은 동일해야 한다.

### 24.2 실행 스크립트

다음 파일을 추가한다.

```text
scripts/test_stage2_interface_fix123.sh
```

필수 검사:

```bash
grep -F 'const double cphi_r = -cphi;' src/shared/surface.c
grep -F 'const double sphi_r = -sphi;' src/shared/surface.c
grep -F 'const double mu_i_s = -mu_i;' src/shared/surface.c
grep -F 'const double mu_o_s = -mu_o;' src/shared/surface.c
grep -F 'T_pair[0*3 + 2] = -T_pair[0*3 + 2];' src/shared/surface.c
grep -F 'T_pair[1*3 + 2] = -T_pair[1*3 + 2];' src/shared/surface.c
```

또한 단계 26에서 기각한 확산상단 단순 부호반전이 소스에 존재하면 실패시킨다.

기대 출력:

```text
PASS: air-to-water signed-rotation and Fourier U-column regression
PASS: stage-2 interface FIX1+FIX2+FIX3 source audit
```

## 25. 단계 20 — OSOAA TAW 전 성분 감사

### 25.1 고정 입력

```text
TAW file = TAW-1.340-03.0-RadMU48-NB200-SZA40.000-TSZA28.665
active nodes = 48 x 48
m = 0..4
matrix = 3 x 3 full
wind = 3 m/s
n_water = 1.34
sigma_type = 1
q_convention = 1
n_phi = 1024
```

### 25.2 승인 기준

각 성분과 mode에 대해 다음을 저장한다.

```text
mean(abs(delta))/max(abs(reference))
RMS(delta)/max(abs(reference))
max(abs(delta))/max(abs(reference))
correlation
best-fit scale
```

FIX123 승인 결과는 다음 범위여야 한다.

```text
worst mean/max_ref <= 0.0023272%
worst RMS/max_ref  <= 0.020464%
minimum correlation >= 0.999995
```

극단적 grazing node의 단일점 max는 약 0.653%까지 허용하되, 해당 위치와 절대 커널 크기를 보고서에 기록한다.

### 25.3 과거 해석 정정

FIX1-only 감사에서 M31/M32의 반대부호를 OSOAA의 “U 행 저장부호”로 해석한 결론은 폐기한다. FIX123 이후의 올바른 공통계약은 다음과 같다.

```text
physical angle-space rotation: FIX1+FIX2
m>0 Fourier storage: incoming-U column M13/M23 sign folded once
```

## 26. 단계 21 — 확산상단 U 가설의 역사적 보류와 해제 조건

### 26.1 FIX123 시점의 보류 근거

443 nm, SZA40, VZA60, RAA90, wind3, `aDOM440=20 m^-1`, water `n_mu=64`, `m_max=4`의 한 종단 시험에서는 molecular incoming-U 세 항의 부호를 바꿀 때 Q가 개선되고 U가 악화되었다. 이 시점에는 외부 U field의 저장규약과 주 SOS 연산자의 홀수 signed-index 채널을 완전히 추적하지 못했으므로 생산 변경을 보류하였다.

이 보류는 안전한 잠정조치였지만, 최종 물리 판정은 아니었다. 종단 Q/U 비는 여러 잔여 오차의 상쇄를 포함하므로 해당 값만으로 국소 source operator의 부호를 판정해서는 안 된다.

### 26.2 보류 해제에 요구된 시험

각 Fourier mode, incident node 및 I/Q/U basis에 대해 다음 두 경로를 같은 source normalization으로 직접 비교한다.

1. `rt_sos_operator_apply_vector()`
2. `ocrt_add_diffuse_top_primary()`

비교 manifest에는 signed input/output index, `xpl/xrl/xtl`, particle/molecular weight, Gauss weight, attenuation, U-output wrapper sign 및 external U storage convention을 기록한다.

### 26.3 등가빔 제한

`mu_sun_water_override>0`에서 `atm.beam_q=0`을 사용하는 equivalent-beam path는 확산 하늘광의 full Q/U Fourier field를 보존하지 않는다. 이 경로는 intensity fallback 또는 진단용으로만 유지하며 편광 정본으로 사용하지 않는다.

## 27. 단계 22 — FIX123 72개 실행 및 600셀 재검증

### 27.1 실행계약

단계 0의 72개 실행 manifest와 OSOAA 기준표를 그대로 사용한다. FIX1-only와 FIX123 사이에서 입력 데이터, depth, quadrature, Fourier order, SOS cap 또는 Rrs extraction branch를 변경하지 않는다.

### 27.2 성공조건

```text
72/72 runs complete
24 full-grid rows per run
water_converged = 1 for every run
600/600 matched cells
```

### 27.3 필수 산포도

다음 각 scope마다 Rrs/rrs I/Q/U 6개 산포도를 독립 생성한다.

```text
all600
purewater150
purewater84_sza_le40
purewater66_sza80
gate_sza40_vza60
```

Q/U의 주 통계는 `RMS(delta component / OSOAA I)`로 고정한다. I는 MAPE를 병기한다.

### 27.4 승인된 결과 해석

FIX123은 FIX1-only보다 최종 aggregate 오차를 소폭 악화시킨다.

```text
all600 Rrs Q: 9.6507 -> 9.8699%
all600 Rrs U: 7.2323 -> 7.3925%
all600 rrs Q: 9.6010 -> 9.8236%
all600 rrs U: 7.1738 -> 7.3361%

purewater SZA80 Rrs Q: 13.0636 -> 13.3605%
purewater SZA80 Rrs U:  9.7614 ->  9.9786%
```

이는 rollback 조건이 아니다. 극한해와 TAW operator가 개선되므로 기존 FIX1-only 결함이 다른 미해결 오차와 상쇄되고 있었던 것으로 기록한다.

## 28. 단계 23 — 계산시간 및 병목 판정

### 28.1 benchmark 조건

```text
wavelength = 443 nm
pure water
SZA/VZA/RAA = 40/60/90 deg
wind = 3 m/s
water n_mu = 48
water m_max = 2
1 warm-up + 5 interleaved measured runs per binary
```

### 28.2 결과

```text
FIX1-only mean   = 2.477739 s
FIX123 mean      = 2.476390 s
mean change      = -0.0544%

FIX1-only median = 2.455856 s
FIX123 median    = 2.479298 s
median change    = +0.9546%
```

5% 증가 기준에 미달하므로 별도 병목 프로파일링은 수행하지 않는다. FIX123은 상수 부호와 두 방향여현의 사용만 바꾸며 신규 반복문, 적분, 할당을 추가하지 않는다.

## 29. 단계 24 — FIX123 checkpoint 버전, 빌드 및 동결

### 29.1 버전 문자열

`src/main.c`의 버전을 다음으로 갱신한다.

```text
OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123
```

### 29.2 빌드

```bash
./scripts/build_release_v1.2.sh build/ocrt_interface_fix123
./build/ocrt_interface_fix123 --version
sha256sum build/ocrt_interface_fix123
```

승인 바이너리 SHA-256:

```text
4a2269a13f81f4112dfb351ec3cac3dc5548915a8913670efa6efe2aa79180c8
```

### 29.3 필수 회귀 실행

```bash
bash scripts/test_stage2_interface_fix123.sh
bash scripts/test_stage2_water_raa_fix.sh
python tests/test_water_raa_output_regression.py \
  --binary build/ocrt_interface_fix123
```

기존 build/release smoke와 패키지 ZIP integrity 검사도 수행한다.

## 30. FIX123 checkpoint 승인 및 보류 상태

### 30.1 생산 승인

- 수중 공개 RAA 및 U positive-sine 출력 수정
- 공기→물 rough-interface FIX1+FIX2
- m>0 TAW incoming-U column FIX3
- 극한해 및 Fourier storage 회귀시험
- OSOAA TAW 9성분 감사도구
- 72/72, 600/600 재검증
- 실행시간 증가 없음

### 30.2 분석 승인

- 확산상단 polarized boundary-source path를 다음 최우선 분석축으로 설정
- equivalent-beam `beam_q=0` 편광 제한 명문화
- SZA40과 SZA80를 동시 gate로 유지
- Rrs 두 extraction branch 병행 유지

### 30.3 보류 또는 금지

- `ocrt_add_diffuse_top_primary` molecular U coefficients의 단순 전역 부호반전
- M31/M32/M33 U-row 전체 반전
- FIX2 없는 FIX3 단독 적용
- 최종 aggregate MAPE 감소만을 기준으로 FIX123 rollback
- level-27+TWA 또는 level-26 subtraction을 유일한 Rrs 정본으로 조기 확정

## 31. FIX123 checkpoint 납품물 요구사항

최종 배포묶음에는 다음을 포함한다.

```text
full source + approved binary
this cumulative instruction
additional defect review report
session-start-to-final cumulative patch
FIX1-only-to-FIX123 delta patch
TAW/polar/DTP validation package
72-run status + 600-cell matched table
Rrs/rrs IQU scatterplots
runtime raw and summary CSV
SHA-256 manifest
```

## 32. 단계 25 — 다른 세션 인수인계 묶음의 독립 감사

### 32.1 검토 자료

다음 첨부를 읽기 전용으로 해제한다.

```text
ocrt_osoaa_session_2026-07-22.tar.gz
```

핵심 자료는 `HANDOFF_FINAL_2026-07-22.md`, `PATCH_water_FIX4_candidate.diff`, `compare_entry.py`, Gate B 자료 및 계면 분석문서이다. 문서의 결론을 그대로 적용하지 않고 실제 patch와 production source를 대조한다.

### 32.2 주장 분류

- 결함 A, B, C: 기존 독립 검증과 일치하므로 유지한다.
- 결함 D 코드 대조: 유력하며 직접 closure가 필요하다.
- 결함 D “적용 금지”: 한 종단 비율 시험에 의존하므로 재판정한다.
- “계면 0.05% 종결”: SZA40 강흡수 계면 gate의 국소 결론으로만 인정한다.
- “최소 수중 광학두께 30”: 보수 guard일 뿐 보편 물리 상수로 사용하지 않는다.

### 32.3 첨부 하네스 입력계약 감사

`compare_entry.py`는 OSOAA에 `depth_m=48/(a+b)`를 전달하지만 OCRT 명령에는 동일 물리 깊이를 전달하지 않는다. 또한 412 nm 순수해수는 수중 SOS 20차에서 종료한다. 따라서 해당 하네스는 방위각·계면 원인진단에는 사용할 수 있으나, 입력 완전일치가 요구되는 최종 1% parity gate로 승인하지 않는다.

후속 하네스는 `water_depth_m`, `water_max_orders`, `n_mu`, `m_max`, surface matrix order, phase function, depolarization 및 Rrs extraction branch를 행 단위 manifest에 명시한다.

## 33. 단계 26 — external U Fourier 규약의 소스 추적

다음 경로를 줄 단위로 확인한다.

```text
atmospheric total_u at BOA
 -> boa_export->U_per_m
 -> rt_air_water_couple_atm_to_water ordinary 3x3 contraction
 -> coupled.U_inwater_per_m
 -> memcpy to ocrt_dt_U
 -> w_opts.ext_top_U
 -> ocrt_add_diffuse_top_primary
```

어느 단계에서도 추가 U 부호변환이 없어야 한다. FIX3은 m>0 incoming-U storage sign을 TAW `M13/M23`에 이미 한 번 접어 넣으므로, downstream field는 water SOS와 동일한 sine-Fourier U convention이다.

검증 결과 이 조건이 충족되었다. 따라서 `ext_top_U`를 반대부호 field로 해석하여 FIX4를 기각해서는 안 된다.

## 34. 단계 27 — diffuse-top source/operator mode–node–basis closure

### 34.1 신규 시험

다음 파일을 추가한다.

```text
tests/test_diffuse_top_primary_operator_closure.c
scripts/test_stage2_diffuse_top_uclosure.sh
```

시험은 pure Rayleigh와 synthetic Rayleigh–aerosol 혼합에서 m=0…4, 모든 incident node 및 I/Q/U basis를 순회한다.

Discrete incident field는 injector에서 다음과 같이 표현된다.

```text
AI/AQ/AU = 2 * Gauss_weight * incident_coefficient
ch_c     = 0.5 * exp(-tau/mu)
```

따라서 두 인자의 곱은 production SOS operator의 `weight * field * attenuation`과 정확히 동일하다.

### 34.2 수정 전 결과

FIX123은 I/Q basis에서 machine precision으로 닫히지만 incoming-U basis에서 다음 불일치를 보인다.

```text
m=1 max abs = 4.512674445328e-2
m=2 max abs = 2.239091465088e-2
```

external U를 인위적으로 반전하면 닫힌다는 사실은 injector의 incoming-U column sign이 반대임을 뜻한다.

### 34.3 FIX4 적용

`src/rt_water_rt.c::ocrt_add_diffuse_top_primary()`의 molecular incoming-U 세 항을 다음과 같이 변경한다.

```diff
- const double ruI = ray_on ? +gamma2 * atm->xpl[j] * xtlc_u : 0.0;
- const double ruQ = ray_on ? +alpha2 * atm->xrl[j] * xtlc_u : 0.0;
- const double ruU = ray_on ? +alpha2 * atm->xtl[j] * xtlc_u : 0.0;
+ const double ruI = ray_on ? -gamma2 * atm->xpl[j] * xtlc_u : 0.0;
+ const double ruQ = ray_on ? -alpha2 * atm->xrl[j] * xtlc_u : 0.0;
+ const double ruU = ray_on ? -alpha2 * atm->xtl[j] * xtlc_u : 0.0;
```

`src_u`에 존재하는 외부 전역 음수는 유지한다. particle/table kernel 또는 `raU/rqU`는 변경하지 않는다.

### 34.4 수정 후 승인 기준과 결과

```text
pure Rayleigh max abs closure = 2.776e-17
mixed medium max abs closure  = 5.551e-17
```

허용 기준 `1e-13`을 충분히 만족한다. 이 시험은 OSOAA를 사용하지 않으므로 FIX4를 내부 연산자 불변량으로 승인한다.

## 35. 단계 28 — FIX4 72개 실행 및 600셀 재검증

### 35.1 실행계약

단계 0과 FIX123 checkpoint의 72개 manifest를 그대로 사용한다.

```text
atmosphere n_mu=48, n_layers=40, m_max=2
water shared grid
sos_max_orders=100
OCRT_DEBUG=1 and --debug-water-max-orders 100
wind=3 m/s
n_water=1.34
T=20 C, S=38.4 PSU, dpol_water=0.039
gas columns=0, aerosol=0
```

### 35.2 성공조건

```text
72/72 physical runs
24 full-grid rows per run
water_converged=1
600/600 matched cells
pure-water 412-nm water orders = 29/29/28 for SZA 0/40/80
```

모든 조건을 만족하였다.

### 35.3 결과 해석

FIX4는 전체 및 SZA80 편광 잔차를 크게 감소시킨다.

```text
all600 Rrs Q: 9.8699 -> 3.6677%
all600 Rrs U: 7.3925 -> 3.2794%
all600 rrs Q: 9.8236 -> 3.4496%
all600 rrs U: 7.3361 -> 3.1436%

pure-water SZA80 Rrs Q: 13.3605 -> 4.5349%
pure-water SZA80 Rrs U:  9.9786 -> 3.9543%
pure-water SZA80 rrs Q: 13.3118 -> 4.3228%
pure-water SZA80 rrs U:  9.9037 -> 3.7766%
```

반면 pure-water SZA≤40 및 I MAPE는 일부 악화된다. 이는 rollback 조건이 아니다. 내부 source/operator closure가 machine precision으로 개선되므로 기존 결함이 저 SZA의 다른 잔여 오차를 상쇄하고 있었던 것으로 기록한다.

### 35.4 필수 산포도

최종 보고에는 다음 6개를 독립 그림으로 포함한다.

```text
purewater150_Rrs0plus_I_scatter_dtp_uclosure.png
purewater150_Rrs0plus_Q_scatter_dtp_uclosure.png
purewater150_Rrs0plus_U_scatter_dtp_uclosure.png
purewater150_rrs0minus_I_scatter_dtp_uclosure.png
purewater150_rrs0minus_Q_scatter_dtp_uclosure.png
purewater150_rrs0minus_U_scatter_dtp_uclosure.png
```

Q/U는 `RMS(delta component / OSOAA I)`, I는 MAPE를 사용하며 각 그림에 FIX123과 FIX4 값을 병기한다.

## 36. 단계 29 — FIX4 계산시간 및 병목 판정

동일 443 nm 순수해수 조건에서 1회 예열 후 3회씩 교차 측정한다.

```text
FIX123 mean       = 2.649036 s
FIX123+FIX4 mean  = 2.640879 s
mean change       = -0.3079%

FIX123 median      = 2.653287 s
FIX123+FIX4 median = 2.640050 s
median change      = -0.4989%
```

계산시간 증가는 없다. FIX4는 상수 부호 세 개만 변경하며 신규 반복문, 적분, 메모리 할당 또는 phase-function 평가를 추가하지 않는다. 5% 병목 프로파일링 기준에 해당하지 않는다.

## 37. 단계 30 — 최종 버전, 빌드 및 회귀시험

### 37.1 버전

```text
OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure
```

### 37.2 빌드

```bash
bash scripts/build_release_v1.2.sh build/ocrt_stage2_dtp_uclosure
./build/ocrt_stage2_dtp_uclosure --version
sha256sum build/ocrt_stage2_dtp_uclosure
```

승인 바이너리 SHA-256:

```text
7e5b2938e3529cc068a70cd11514dfa22cbbc3563aa5ec8a910ce4a6aa258fd7
```

### 37.3 필수 회귀

```bash
bash scripts/test_stage2_interface_fix123.sh
bash scripts/test_stage2_diffuse_top_uclosure.sh
bash scripts/test_stage2_water_raa_fix.sh
python tests/test_water_raa_output_regression.py \
  --binary build/ocrt_stage2_dtp_uclosure
```

최종 ZIP을 별도 디렉터리에 다시 풀어 동일 회귀와 `--version`을 반복한다.

## 38. 최종 승인·보류 상태와 다음 분석축

### 38.1 생산 승인

- water public RAA/U output correction
- air-to-water FIX1+FIX2
- m>0 TAW incoming-U FIX3
- diffuse-top molecular incoming-U FIX4
- source/operator closure regression
- 72/72 및 600/600 재검증
- 계산시간 증가 없음

### 38.2 계속 보류

- level-27+TWA와 level-26 subtraction 중 하나를 유일한 Rrs 정본으로 조기 확정
- 해수 최소 광학두께 30을 보편 상수로 강제
- 입력 물리 깊이가 다른 비교 하네스의 최종 1% gate 사용
- SZA40 한 조건의 종단 Q/U 비만을 근거로 부호 수정 승인 또는 기각
- hydrosol truncation 입력을 맞추지 않은 Chl/TSM 비교

### 38.3 후속 최우선 분석

1. SZA≤40에서 FIX4 이후 증가한 I/Q/U 잔차를 Ed 정규화, depth, Rrs branch 및 source I/Q 열로 분리한다.
2. SZA80 잔여 3–4%를 direct source, BOA diffuse field, TAW contraction 및 water scattering order별로 분해한다.
3. 두 코드에 동일한 physical depth와 SOS cap을 강제하는 `compare_entry` 개정판을 만든다.
4. Rrs 두 extraction branch를 같은 결과행에 병기하여 branch 차이를 degree of freedom에서 제거한다.

## 39. 개정 3 최종 납품물 요구사항

```text
full source + approved binary
this revision-3 cumulative instruction
other-session review report
session-start-to-final cumulative patch
FIX123-to-DTP-U-closure delta patch
source/operator closure raw and regression log
72-run status + 600-cell matched table
Rrs/rrs IQU scatterplots
runtime raw and summary CSV
SHA-256 manifest
```



## 40. 단계 31 — 첨부 극한 패치의 무결성 및 적용 전 판정

### 40.1 첨부 무결성

`files (41).zip`에서 실제 제공된 다음 네 파일의 SHA-256은 manifest와 일치하여야 한다.

```text
NOTE_double_pole_and_build_2026-07-24.md
PATCH_surface_POLE_2026-07-24.diff
surface.c
surface_double_pole_perturbed.c
```

manifest에 기재되었으나 ZIP에 포함되지 않은 다음 다섯 파일은 외부 미제공 자료로 기록한다.

```text
ADDENDUM2_atmval_pole_2026-07-24.md
atmval_v4_2026-07-24.csv
atmval_v4_2026-07-24.png
polecont.c
probe_nadir_pole_constraint.py
```

따라서 첨부의 종단 통계나 224조건 전수검사 주장을 그대로 인용하지 않고, 포함된 코드 패치의 물리식과 별도 회귀시험을 독립 검증한다.

### 40.2 적용 승인 기준

수정은 다음 세 조건을 모두 만족할 때 승인한다.

1. 정확한 극과 일방 극한의 Mueller 행렬 연속성이 복원될 것.
2. 연직 출력의 Q/U가 불필요한 m=0 누출 없이 spin-2 공변성을 만족할 것.
3. 기존 비극 식, FIX1–FIX4 및 대기 하네스가 변경되지 않을 것.

종단 Rrs/rrs 또는 TOA 상대오차 감소 여부는 승인 기준으로 사용하지 않는다.

## 41. 단계 32 — `src/shared/surface.c`의 정확한 극 한계 적용

### 41.1 공통 문턱

파일 상단에 다음 상수를 추가한다.

```c
#define OCRT_SURF_POLE_EPS 1.0e-8
```

주석에는 일반 회전식이 `sin(theta_self)->0`에서 0/0이 되며, 해석 극한이 다음임을 명시한다.

```text
cos(sigma) -> -sign(mu_self) cos(phi)
sin(sigma) -> sin(phi)
```

### 41.2 수정 대상 네 함수

다음 네 함수의 입사측·출사측 회전각 계산에 동일 원칙을 적용한다.

```text
surface_R_coxmunk_trig
surface_T_coxmunk_trig
surface_R_ww_coxmunk_trig
surface_T_aw_coxmunk_trig
```

각 방향의 `sin(theta_self)`가 문턱 이하면 극한식을 사용하고, 그보다 크면 기존 일반식을 글자 단위로 유지한다.

대표 패턴은 다음과 같다.

```c
if (s_self <= OCRT_SURF_POLE_EPS) {
    ci = -copysign(1.0, mu_self) * cos_phi;
    si = sin_phi;
} else {
    /* 기존 일반식 및 정규화 */
}
```

공기→물 FIX1+FIX2 경로에서는 이미 정의된 signed propagation cosine과 `pi` 이동 상대방위각을 그대로 사용한다. 극한 패치가 FIX1+FIX2를 다시 해석하거나 부호를 중복 적용해서는 안 된다.

### 41.3 금지사항

- 비극 구간의 회전식을 변경하지 않는다.
- 정확한 극에서 임의로 identity rotation을 강제하지 않는다.
- 양쪽 극 축퇴를 종단 오차 감소 목적으로 조정하지 않는다.
- FIX3의 m>0 incoming-U 열 및 FIX4의 diffuse-top U 열을 변경하지 않는다.
- 대기 입력과 1단계 하네스 설정을 변경하지 않는다.

## 42. 단계 33 — 극한 및 양쪽 극 회귀시험

### 42.1 새 단위시험

다음 파일을 추가한다.

```text
tests/test_surface_pole_limit.c
scripts/test_stage2_surface_pole_limit.sh
```

시험은 네 kernel 각각에서 입사 또는 출사 방향 하나만 극인 8개 경우를 포함한다.

수정 전 대표 결과:

```text
exact-pole / one-sided-limit relative discontinuity ≈ 1.49995–1.50005
```

수정 후 승인 기준과 실측 결과:

```text
R_air max relative discontinuity  = 9.43e-6
T_wa  max relative discontinuity  = 4.92e-5
R_ww  max relative discontinuity  = 5.31e-6
T_aw  max relative discontinuity  = 4.92e-5
acceptance threshold              = 1e-4
```

연직 spin-2 시험은 다음을 요구한다.

```text
m=0 leakage / spin-2 amplitude < 1e-12
paired Q/U m=2 amplitude mismatch < 1e-12
```

수정 후 네 kernel의 m=0 누출은 약 `6e-17–8e-17`, Q/U pair mismatch는 `0–1.6e-15`이다.

### 42.2 양쪽 극 축퇴

양쪽 방향이 모두 극이면 상대방위각 기준면은 접근경로에 의존한다. 생산 격자에서 해당 쌍은 삽입 special node이며 경계 적분에서 zero quadrature weight와 곱해진다.

첨부의 perturbation source를 사용하되, 본 세션에서는 최소 다음 독립 범위를 재실행한다.

```text
pure water
SZA = 0 deg
bands = 412, 443, 490, 555, 660, 865 nm
full grid = 24 rows per band
```

정본과 double-pole identity-perturbation 판의 6개 CSV가 모두 byte-identical이어야 한다. 첨부가 주장한 224조건 전체를 독립 반복하지 않았으면 그 사실을 보고서에 명시한다.

## 43. 단계 34 — 호스트 독립 빌드 규약

### 43.1 릴리스 스크립트

다음 파일에서 `-march=native`를 제거한다.

```text
scripts/build_release_v1.1.sh
scripts/build_release_v1.2.sh
scripts/regression_ext_bottom_mode_bounds.sh
```

공통 구현:

```bash
OCRT_MARCH=${OCRT_MARCH:-cascadelake}
...
-march="$OCRT_MARCH"
```

기본 릴리스는 `cascadelake`로 빌드한다. 다른 CPU 기준이 필요한 경우 명시적인 `OCRT_MARCH`를 사용하고 기존 승인 해시를 재사용하지 않는다.

### 43.2 재현성 시험

다음 스크립트를 추가한다.

```text
scripts/test_release_build_reproducibility.sh
```

동일 소스와 동일 toolchain으로 독립 디렉터리에 두 번 빌드하여 binary SHA-256과 byte identity를 비교한다.

본 개정 승인값:

```text
version = OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure-pole-limit
march   = cascadelake
sha256  = 6c7f4b61d311323374233e0be9643541570d3b15289457ca9336d14d8f12ffde
```

두 독립 빌드는 동일 SHA-256을 생성하여야 한다. 이 주장은 동일 compiler/linker 버전과 환경에 한정한다.

### 43.3 CSV bit-exact 도구

다음 도구를 추가한다.

```text
scripts/compare_csv_bitexact.py
```

- `--raw`: 파일 바이트 직접 비교
- 기본: Python 내장 `float()`로 numeric cell을 binary64로 변환하여 비트 비교
- 문자열 cell은 정확한 텍스트 비교

`pandas` 파싱 결과를 bit-exact assertion에 사용하지 않는다.

## 44. 단계 35 — 기존 회귀의 비회귀 확인

다음 시험을 모두 실행한다.

```bash
bash scripts/test_stage2_surface_pole_limit.sh
bash scripts/test_stage2_interface_fix123.sh
bash scripts/test_stage2_diffuse_top_uclosure.sh
bash scripts/test_stage2_water_raa_fix.sh
bash scripts/regression_ext_bottom_mode_bounds.sh
bash scripts/regression_raa_convention.sh
bash scripts/test_release_build_reproducibility.sh
```

`regression_raa_convention.sh`는 `tests/test_raa_convention.c` 링크 시 `src/rt_fourier.c`를 포함하도록 수정한다. 이는 물리 변경이 아니라 기존 회귀스크립트의 누락 링크 수정이다.

필수 결과:

```text
FIX123 interface regression: PASS
FIX4 pure-Rayleigh closure: 2.776e-17
FIX4 mixed closure:         5.551e-17
water public RAA regression: PASS
ext-bottom bounds regression: PASS
RAA glint anchors: unchanged
```

## 45. 단계 36 — 72개 실행 및 종단 영향 범위

### 45.1 입력계약

단계 0의 72개 manifest를 그대로 사용한다. 대기 설정, Rayleigh optical depth, depolarization, layer/Fourier/SOS 옵션을 변경하지 않는다.

승인 결과:

```text
72/72 runs
24 rows/run
water_converged = 1 for all runs
max water order = 29
600/600 regression-reference cells
```

### 45.2 동일 아키텍처 수정 전·후 비교

수정 전 FIX123+FIX4와 수정 후 pole-limit 바이너리를 모두 `-march=cascadelake`로 빌드한다.

SZA 40°, VZA 0°, 순수해수 6개 파장·8개 RAA에서:

```text
TOA rho_I: unchanged
TOA rho_Q/U: exact-pole reference-plane component corrected
Rrs I/Q/U: unchanged
rrs I/Q/U: unchanged
```

443 nm의 최대 절대 변화는 대략 다음이다.

```text
TOA rho_Q max |delta| = 4.48e-4
TOA rho_U max |delta| = 3.78e-4
```

따라서 이 패치는 수면 연산자의 exact-nadir 편광 기준면을 교정하지만, 기존 Rrs/rrs 해양 결과를 우연한 오차상쇄 목적으로 변경하지 않는다.

### 45.3 Rrs/rrs 산포도

동일 OSOAA 회귀 기준표를 사용하여 다음 6개 그림을 각각 생성한다.

```text
purewater150_Rrs0plus_I_scatter.png
purewater150_Rrs0plus_Q_scatter.png
purewater150_Rrs0plus_U_scatter.png
purewater150_rrs0minus_I_scatter.png
purewater150_rrs0minus_Q_scatter.png
purewater150_rrs0minus_U_scatter.png
```

현재 기준표는 최신 OSOAA 수심 정정 이전 자료이므로 물리 정확도 최종 승인값이 아니라 코드 비회귀 대조용으로만 사용한다. pole-limit 패치 전후 지표는 동일하다.

## 46. 단계 37 — 계산시간 및 병목 판정

동일 compiler, `-march=cascadelake`, `OMP_NUM_THREADS=1`, 443 nm 순수해수 SZA 40° full-grid 조건에서 각 variant 1회 예열 후 5회 교차 측정한다.

```text
before mean   = 4.831627 s
pole mean     = 4.856947 s
mean change   = +0.5241%

before median = 4.819624 s
pole median   = 4.873949 s
median change = +1.1272%
```

5% 미만이므로 병목 프로파일링을 수행하지 않는다. 새 반복문·적분·동적 할당은 추가되지 않았다.

## 47. 단계 38 — 최종 버전 및 패키지 검증

### 47.1 버전과 바이너리

```text
OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure-pole-limit
```

승인 바이너리 SHA-256:

```text
6c7f4b61d311323374233e0be9643541570d3b15289457ca9336d14d8f12ffde
```

### 47.2 최종 ZIP 검증

최종 ZIP을 빈 디렉터리에 풀고 다음을 수행한다.

```bash
./scripts/build_release_v1.2.sh build/ocrt_stage2_pole_limit_clean
cmp build/ocrt_stage2_pole_limit build/ocrt_stage2_pole_limit_clean
bash scripts/test_stage2_surface_pole_limit.sh
bash scripts/test_stage2_interface_fix123.sh
bash scripts/test_stage2_diffuse_top_uclosure.sh
bash scripts/test_stage2_water_raa_fix.sh
```

포함 바이너리와 clean rebuild 바이너리가 byte-identical이어야 한다.

## 48. 개정 4 최종 승인·보류 상태

### 48.1 동결 모듈

- 1단계 대기 모듈과 대기 하네스
- public water RAA/U output reconstruction
- air-to-water FIX1+FIX2+FIX3
- diffuse-top incoming-U FIX4
- 네 rough-interface kernel의 exact-pole meridian limit
- fixed-ISA release build contract

### 48.2 계속 미해결

- corrected-depth OSOAA 기준의 SZA 80° `rrs(0-)`
- corrected-depth OSOAA `Rrs(0+)` 두 extraction branch
- rough direct-solar boundary field와 in-water solver의 완전 분리 폐합
- water-side `R_ww` 전 성분·전 mode 독립 감사
- 동일 물리수심·층경계·n_mu·m_max·SOS tolerance를 강제한 최종 1% gate

극한 패치가 이 미해결 항목을 해결했다고 주장하지 않는다.

## 49. 개정 4 최종 납품물 요구사항

```text
full source + approved fixed-ISA binary
revision-4 cumulative instruction
surface pole-limit review report
session-start-to-final cumulative patch
DTP-U-closure-to-pole-limit incremental patch
exact-pole and spin-2 regression logs/figures
same-architecture before/after audit
72-run status + 600-cell regression table
Rrs/rrs IQU scatterplots
runtime raw/summary
clean-zip rebuild verification
SHA-256 manifest
```


## 50. 단계 39 — water-to-air coupling 64절점 silent clamp 제거

### 50.1 수정대상

```text
src/rt_air_water_coupling.c
aw_interp_on_unsorted()
```

기존의 `ord[64]`, `ms[64]`, `tmp[64]` 및 `if (n > 64) n = 64`를 제거한다. 이 함수의 `n`은 requested Gauss node 수가 아니라 관측·태양·천저 zero-weight node가 추가된 filled-grid 크기이므로 절단하면 안 된다.

### 50.2 구현계약

```text
AW_INTERP_STACK_N = 320
n <= 320 : stack workspace
n > 320  : heap fallback
allocation failure : abort; silent clamp 금지
```

`n <= 64`에서는 기존 산술순서를 유지한다.

### 50.3 변경금지 범위

다음 파일과 입력계약을 수정하지 않는다.

```text
atmosphere solver and atmosphere harness
FIX1+FIX2+FIX3 surface kernel
FIX4 diffuse-top source
exact-pole equations
public water RAA/U reconstruction
in-water SOS
```

## 51. 단계 40 — 전용 회귀시험 등록

다음 파일을 추가한다.

```text
scripts/test_stage2_coupling_clamp_fix.py
scripts/test_stage2_coupling_clamp_fix.sh
scripts/test_stage2_coupling_clamp_single_parity.sh
scripts/check_stage2_coupling_clamp_single_parity.py
tests/test_aw_interp_workspace.c
scripts/test_aw_interp_workspace_sanitized.sh
```

필수 게이트:

1. source audit: 64절점 clamp 완전 제거
2. `n_mu_water=48`: 수정 전·후 byte-identical
3. `n_mu_water=64/96`: 천저 Rrs/TOA Q의 forbidden `m=0` 성분 0
4. 모든 node 수에서 `rrs(0−)` I/Q/U 불변
5. stack/heap workspace ASAN/UBSAN 통과
6. single-target/full-grid Rrs/rrs tolerance 일치

## 52. 단계 41 — 기존 동결 모듈 비회귀

다음 기존 시험을 그대로 재실행한다.

```bash
bash scripts/test_stage2_water_raa_fix.sh
bash scripts/test_stage2_interface_fix123.sh
bash scripts/test_stage2_diffuse_top_uclosure.sh
bash scripts/test_stage2_surface_pole_limit.sh
bash scripts/regression_raa_convention.sh
bash scripts/regression_ext_bottom_mode_bounds.sh
```

FIX4 closure는 `2.776e-17` 및 `5.551e-17`을 유지하여야 한다. 대기 하네스는 변경하지 않는다.

## 53. 단계 42 — 72실행·600셀 비회귀

세션의 동결 manifest를 그대로 사용한다.

```text
72/72 successful
72/72 water_converged
72/72 CSV byte-identical to pole-limit baseline
maximum water SOS order = 29
600-cell numeric columns max absolute difference = 0
```

이 게이트는 `n_mu_water=48`이므로 clamp 미발동 영역의 완전 비회귀를 보장한다. OSOAA 표는 corrected-depth 이전 자료이므로 정확도 최종 판정에는 사용하지 않는다.

## 54. 단계 43 — 생산격자 영향검사

443 nm, SZA 40°, VZA 0/30/60°, 8개 RAA, 24셀 full-grid에서:

```text
n_mu_water=48 : byte-identical
n_mu_water=64 : rrs I/Q/U unchanged; forbidden nadir m0 removed
n_mu_water=96 : rrs I/Q/U unchanged; forbidden nadir m0 removed
```

96절점에서 수정 전 clamp 영향은 Rrs 성분최대 대비 I 9.568%, Q 30.066%, U 35.718%까지 나타날 수 있다. 이는 reference-code 오차가 아니라 잘못된 기존 OCRT 출력의 변경량이다.

## 55. 단계 44 — 성능 및 병목

동일 fixed-ISA/단일스레드 조건에서 variant별 3회 교차 측정한다.

```text
n_mu=48 mean -1.765%, median -1.974%
n_mu=64 mean -7.751%, median -1.122%  (before outlier 포함)
n_mu=96 mean +0.638%, median +0.254%
```

5% 성능저하가 없으므로 병목 프로파일링은 수행하지 않는다.

## 56. 단계 45 — C/Python 동시개발 판정

현재 패키지에는 이 coupling interpolation의 독립 Python 물리 구현이 없다. Python 하네스·분석도구에는 동일 수정 지점이 없으므로 C만 변경한다. 향후 Python core가 도입되면 전체 filled-grid 보존 및 silent-clamp 금지 계약을 동일하게 구현한다.

## 57. 단계 46 — 최종 버전 및 패키지

```text
OCRT-v1.2-2026-07-25-KST-stage2-coupling-clamp-fix
```

승인 바이너리 SHA-256:

```text
328f28b8f50b6d98207eabe57309dde3538e4658a1ceec8f53d970b9dda3587a
```

최종 ZIP을 빈 디렉터리에 풀어 fixed-ISA clean rebuild를 수행하고 포함 바이너리와 byte-identical인지 확인한다. 세션 시작 소스→최종 cumulative patch와 pole-limit→clamp-fix incremental patch를 각각 fresh tree에 적용하여 최종 소스와 비교한다.

## 58. 개정 5 최종 상태

### 58.1 동결 모듈

- 개정 4의 모든 동결 모듈
- water-to-air coupling interpolation full-grid preservation
- n_mu_water 48/64/96 및 forbidden-nadir-m0 regression

### 58.2 계속 미해결

- corrected-depth OSOAA SZA 80° `rrs(0−)`
- corrected-depth OSOAA `Rrs(0+)` two-branch extraction
- exact matched-discretization final 1% gate
- direct boundary/in-water separation closure
- water-side `R_ww` full operator audit


# 개정 6 추가 누적 단계 — EAP 17종 generator 및 species scattering 비활성화

## 59. 단계 47 — 첨부 작업지시서와 자료 무결성 확인

1. 작업지시서의 최종 코드 기준선은 과거 v3/v5 계통이므로 파일 전체를 교체하지 않는다.
2. 현재 개정 5 기준선에 EAP 관련 변경만 선택적으로 이식한다.
3. 첨부 11,154행 검증 CSV와 50행 summary를 독립 재계산한다.
4. summary 모든 수치열의 최대 절대차가 0인지 확인한다.
5. 제공 Rrs/rrs I/Q/U 산포도는 원자료에서 18개 독립 그림으로 다시 생성한다.
6. OSOAA 실행파일·누적 patch·surface matrix 해시가 불완전하면 원시 OSOAA 재실행을 주장하지 않는다.

확인값:

```text
validation rows     11154
validation columns  37
summary rows         50
summary columns      26
summary max abs diff 0
```

## 60. 단계 48 — 17종 EAP coated-sphere generator 통합

다음 파일을 추가한다.

```text
src/rt_eap_mie_phase.{c,h}
src/internal/rt_eap_coated_mie.{c,h}
src/internal/rt_eap_species_catalog.{c,h}
src/generated/rt_eap_species_data.inc
tools/eap_mie_writer.c
tools/eap_generate_all.sh
tools/generate_eap_species_catalog.py
tools/eap_catalog_freeze.py
```

공개 API는 하나만 둔다.

```c
int rt_eap_mie_phase_compute(
    rt_eap_species_id_t species_id,
    const double *wavelength_nm,
    size_t n_wavelength,
    const double *theta_deg,
    size_t n_theta,
    double *p11,
    double *p12,
    double *p33);
```

계약:

- species ID 0–16
- 350–850 nm
- theta strictly increasing, 0°–180°
- wavelength-major P11/P12/P33
- `0.5*integral(P11 sin(theta) dtheta)=1`
- caller-owned memory
- no Chl/a/b/bb/moment/file-I/O in public API
- re-entrant and concurrent-call safe

## 61. 단계 49 — EAP generator 자료 및 provenance

17종 대표 파일을 생성한다.

```bash
gcc -std=c11 -O3 -march=cascadelake -Isrc \
    tools/eap_mie_writer.c \
    src/internal/rt_eap_species_catalog.c \
    src/internal/rt_eap_coated_mie.c \
    -lm -o build/eap_mie_writer

sh tools/eap_generate_all.sh representative inputs/water_iop
```

다음 원자료는 catalog 주석에 SHA-256만 있고 이번 패키지에는 파일 자체가 없다.

```text
EAP_invivo_means.csv
92051b77cec90744294536dbfe611b65350532b98d7d0f8c32545bb5e4f273a8

501nm_extended_e1701000.mat
50638c22d4596c38da5f34e10d846188faa225d6fb262409f858fce6e9725c02
```

따라서 raw source→generated catalog 단계는 미재현으로 표시한다. catalog→P11/P12/P33
및 catalog→`.mie` 단계는 현재 패키지에서 재현한다.

## 62. 단계 50 — EAP API 단위시험

다음 조건을 검사한다.

```text
species 11
wavelengths 412, 443 nm
theta 0..180 by 0.5 deg
```

필수 판정:

```text
normalization PASS
P11 >= 0 PASS
|P12| <= P11 PASS
|P33| <= P11 PASS
repeat-call bit identity PASS
invalid species/range error PASS
```

species 11/Deff 0.5 µm를 독립 2회 생성하고 배포 `.mie`와 byte-identical인지 확인한다.

## 63. 단계 51 — canonical 17종 이름 판정

`rt_iop_organic_group_parse()`는 canonical exact match를 legacy substring alias보다 먼저
수행해야 한다. 그렇지 않으면 `eap_synechococcus`가 pico로, `eap_diatoms_pennate`가
micro로 잘못 흡수된다.

기존 pico/nano/micro enum 번호는 유지하고 canonical 17종은 뒤에 추가한다.

## 64. 단계 52 — constituent phytoplankton scattering 비활성화

생산 기본값:

```c
organic_phyto_scattering = 0;
```

Chl>0에서 사용자가 `--ocrt-phyto-group`을 명시하면 종료코드 2로 실패한다. 임의로
기본 phase로 대체해서는 안 된다.

species를 생략한 Chl 실행은 다음과 같이 처리한다.

```text
phytoplankton absorption: active
b_phyto: 0
bb_phyto: 0
particle phase: detritus only
```

비활성화 근거는 현재 L=200 phase moment representation이 강한 forward peak를 가진
대부분 species의 중간각 P11을 음수로 reconstruct할 수 있다는 모듈 단위 진단이다.
종단 reference 오차 감소는 gate 해제 근거로 사용하지 않는다.

## 65. 단계 53 — absorption-only fast path

비활성화 상태에서 discarded phytoplankton phase moment와 bb/b를 계산하지 않는다.
다음 함수를 추가한다.

```c
int rt_iop_eap_phyto_absorption_eval(
    double lambda_nm,
    double chl_mg_m3,
    organic_phyto_group_t group,
    double *a_out_m_inv);
```

또한 생산 단일 Chl 실행은 `rt_iop_organic_init_selected()`로 선택된 absorption table과
detritus만 로드한다. catalog 전체 로딩은 catalog-wide diagnostic에 한정한다.

## 66. 단계 54 — constituent truncation fail-loud

다음 두 플래그는 constituent model에서 구현되어 있지 않으므로 조용히 무시하지 않는다.

```text
--ocrt-mie-truncation
--ocrt-mie-ss-mode
```

오류문은 검증된 fixed-bulk 경로를 명시한다.

```text
--water-model iop --iop-mie-phase FILE --iop-mie-truncation
```

component-level truncation을 복원하려면 removed mass, b-rescaling, truncated bb/b,
truncation-aware cache key를 함께 구현한다. pure-water Rayleigh는 particle truncation에
포함하지 않는다.

## 67. 단계 55 — detritus-only component count

기존 `n_component>=2` 전제는 phyto+detritus 항상 활성 상태에만 유효하다. species
scattering 비활성화 후에는 detritus 한 성분이 정상이다.

```text
minimum active particle phase components: 1
```

0개일 때만 hard error로 처리한다.

## 68. 단계 56 — EAP 비활성화 회귀시험

다음을 자동화한다.

1. canonical species+Chl>0 → rc=2
2. legacy pico/nano/micro+Chl>0 → rc=2
3. species 생략 Chl=0.3 → a_phyto>0, b_phyto=bb_phyto=0
4. detritus b/bb>0
5. phase component=detritus only
6. constituent truncation flags fail-loud
7. fixed-bulk truncation path succeeds

## 69. 단계 57 — 비-Chl byte regression

CDOM 12조건과 TSM 6조건, 총 18조건을 구/신 바이너리에서 실행한다.

```text
stdout byte identity: 18/18
stderr byte identity: 18/18
return code 0: 18/18
```

이 게이트는 이미 검증된 비-Chl 해양경로가 변경되지 않았음을 보장한다.

## 70. 단계 58 — 동결 모듈 집중 회귀

전체 72/600을 반복하여 세션 자원을 소모하지 않고, 변경 영향경계에 맞춘 다음 시험을
재실행한다.

```text
water RAA PASS
FIX1+FIX2+FIX3 PASS
FIX4 pure-Rayleigh closure 2.776e-17
FIX4 Rayleigh+aerosol closure 5.551e-17
exact-pole PASS
coupling clamp n_mu=48/64/96 PASS
single/full-grid parity PASS
direct-glint RAA anchors PASS
```

대기 소스 및 대기 하네스는 변경하지 않는다.

## 71. 단계 59 — 제공 수중자료 정합성 판정

제공 원자료에서 독립 재계산한 수면 아래 I 평균오차 범위:

```text
CDOM 0.0346–0.3409%
TSM  0.0446–0.4943%
Chl  0.1321–0.2483%
```

Chl 자료는 production constituent path가 아니라 다음 fixed-bulk validation bypass다.

```text
constituent scalar IOP extraction
-> detritus-only truncation
-> pure-water Rayleigh re-mixing
-> b = b_w + b_det(1-A/2)
-> fixed-bulk RT
```

따라서 Chl 수치를 production constituent end-to-end 정합성으로 표기하지 않는다.

## 72. 단계 60 — 성능

동일 물리 무-Chl CDOM 단일기하를 구/신 바이너리 각각 5회 교차 측정한다.

```text
old mean   0.228257 s
new mean   0.226025 s
mean change   -0.978%
old median 0.227057 s
new median 0.226742 s
median change -0.139%
```

5% 병목조사 기준 미만이다. Chl absorption-only 단일기하의 절대 실행시간은 별도로
기록하되, 구버전과 물리가 다르므로 성능회귀 비율로 직접 비교하지 않는다.

## 73. 단계 61 — 최종 버전

```text
OCRT-v1.2-2026-07-26-KST-stage2-eap-generator-phyto-scattering-disabled
```

최종 clean build, ZIP 재추출 build, 증분 patch 및 세션 시작점 cumulative patch의
적용성을 검증한다.

## 74. 개정 6 최종 동결 및 미결

### 74.1 추가 동결

- 17종 EAP P11/P12/P33 generator public API
- representative `.mie` generation contract
- explicit species+Chl fail-loud gate
- absorption-only phytoplankton constituent contract
- detritus-only Chl-linked particle phase
- constituent truncation fail-loud contract
- fixed-bulk truncation diagnostic path

### 74.2 계속 미결

- component-level phytoplankton forward-peak truncation
- species-specific EAP scattering 복원
- `FIXED_BULK_BETAL_LMAX` 옵션화 또는 higher-order phase representation
- raw EAP source files provenance completion
- modified OSOAA executable/source provenance completion
- corrected-depth SZA 80° 및 Rrs two-branch final 1% gate
