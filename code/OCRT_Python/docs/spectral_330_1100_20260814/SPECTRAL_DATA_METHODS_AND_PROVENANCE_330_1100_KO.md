# OCRT 330–1100 nm 분광 입력자료 구축·통합 방법론

## 초록

본 문서는 OCRT의 elastic atmosphere–ocean coupled radiative-transfer 계산에 사용되는 분광 입력자료를 330–1100 nm로 확장한 과정, 자료출처, 변환규칙, 보간계약, 물리성 검사, C/Python paired integration 및 제한사항을 기록한다. 자료범위 확장은 runtime에서 누락값을 임의 생성하는 방식이 아니라, 활성 자료가 명시적 파장축을 보유하고 solver가 그 범위 안에서만 보간하도록 구현하였다. 다만 organic detritus vector phase의 330–1100 nm 후보자료는 기본 물리성 검사만 통과했으나 현재 OCRT의 L=200 coefficient representation과 SOS solver에서 부적합하였으므로 production 채택에서 제외하였다.

## 1. 범위와 원칙

### 1.1 목표범위

```text
elastic OCRT: 330–1100 nm
C/Python paired implementation
atmosphere + surface + water coupled RT
```

### 1.2 제외범위

- Raman inelastic coupling
- EAP 17-species production scattering
- legacy CCRR branch의 전면 확장
- OSOAA 자체의 330 nm 확장
- O2-O2/O4 CIA

### 1.3 runtime 데이터 원칙

1. component-specific runtime `if`로 과학값을 생성하지 않는다.
2. 활성 표는 명시적 wavelength axis와 실제 endpoint를 가진다.
3. 보간은 자료범위 내부에서만 수행한다.
4. 범위 밖 요청은 endpoint hold가 아니라 오류로 처리한다.
5. C와 Python에 동일 파일을 byte-identical하게 설치한다.
6. 자료변경과 solver 변경을 분리해 검증한다.
7. 후보자료는 파일 형식만 맞는다는 이유로 production에 채택하지 않는다.

## 2. Pure-water IOP

### 2.1 물리량

순수해수 자료는 흡수계수 `a_w(lambda)`와 총 산란계수 `b_w(lambda)`를 제공한다. OCRT는 분자산란 위상의 대칭성을 전제로

```text
bb_w(lambda) = 0.5 b_w(lambda)
```

를 사용한다.

### 2.2 자료

- 파일: `water_coef_z09_1nm.txt`
- 범위: 200–2449 nm
- 해상도: 1 nm
- 행수: 2,250
- 원 기술: Pope & Fry visible absorption, Kou NIR absorption, Zhang et al. total scattering lineage
- 변환기록: 0.1 nm 기준자료를 `[lambda-0.5,lambda+0.5]` 구간에서 사다리꼴 band-average

첨부된 순수해수 자료와 기존 C/Python canonical 파일은 SHA-256이 동일했다. 따라서 이번 통합에서는 수치를 바꾸지 않았다. 이 항목은 330–1100 nm 전체에 대해 이미 충분한 범위를 제공한다.

## 3. Pure-water absorption temperature coefficient `psi_T`

### 3.1 계산식

```text
a_w(lambda,T) = a_w(lambda,20 C) + psi_T(lambda) [T-20 C]
```

`psi_T`의 단위는 `m^-1 degC^-1`이다.

### 3.2 최종자료

- 범위: 330–1100 nm
- 해상도: 1 nm
- 330–440 nm: exact zero
- 441–1100 nm: Röttgers/WOPP-family 자료를 handoff 규칙대로 PCHIP 계열로 준비

이전 코드의 400–900 nm endpoint clamp를 제거하였다. 현재는 330 또는 1100 nm가 정확한 endpoint로 처리되고 329.999/1100.001 nm는 오류다.

### 3.3 해석상의 주의

330–440 nm zero는 “측정상 모든 수온효과가 물리적으로 정확히 0”이라는 독립적 보편명제가 아니라 handoff에서 채택한 project closure이다. 논문에는 source-derived 구간과 project closure 구간을 구분해 표시해야 한다.

## 4. Default phytoplankton absorption

### 4.1 물리량과 명칭

파일의 값은 chlorophyll-a 분자용액 흡광이 아니라 chlorophyll-specific phytoplankton absorption coefficient `a_ph*(lambda)`이다.

```text
a_ph(lambda) = a_ph*(lambda) TChl-a
```

### 4.2 PLOPS 자료처리

handoff 문서에 기록된 대표스펙트럼 구성은 다음과 같다.

1. 배양 phytoplankton `aph_spec`와 HPLC `TChl a`가 연결되는 자료 선별
2. phase/treatment/concentration/nonnegative QC
3. UV MAA-like outlier 분리
4. 34개 자료, 6분류군 사용
5. 분류군별 파장 중앙값 계산
6. 6개 class median의 중앙값으로 class-balanced spectrum 구성
7. 1 nm 재격자화
8. 5-point centered moving average, reflected padding

따라서 최종 기본값은 single diatom spectrum이 아니다. 합성결과와 가장 가까운 실제시료가 특정 diatom이었다는 사실과 모델 자체가 single-species라는 주장을 구분한다.

### 4.3 장파장 closure

```text
330–750 nm : source-derived class-balanced spectrum
750–800 nm : 750-nm 값에서 800-nm exact zero까지 선형 taper
800–1100 nm: exact zero rows
```

이전의 runtime cutoff branch는 제거하였다. 800–1100 nm의 0은 파일 자체에 기록되므로 C와 Python은 동일 보간기만 사용한다.

### 4.4 production scattering 계약

Stage-2 constituent model에서 phytoplankton 자체의 `b`와 `bb`는 0이다. EAP catalog는 생성·진단용으로 보존하지만 production scattering으로 활성화하지 않는다.

## 5. AHN mineral scalar IOP

### 5.1 성분

- brown earth
- yellow clay
- red clay
- calcareous sand

각 성분은 mass-specific absorption `a*`와 scattering `b*`를 사용한다.

```text
a_min = a*_min TSM
b_min = b*_min TSM
```

### 5.2 330–399 nm 확장

400–420 nm의 canonical 21개 값에 ordinary least squares를 적용하여 330–399 nm를 생성하였다.

```text
y(lambda) = slope lambda + intercept
```

400 nm 이상은 기존값을 보존하였다. PCHIP endpoint extrapolation은 단파장에서 기울기 증폭과 음수 가능성이 있어 채택하지 않았다.

### 5.3 C/Python 계약

8개 scalar 파일은 paired exact 파일이며 C/Python IOP 계산은 테스트 파장집합에서 최대 상대차 `5.49e-15`였다.

## 6. AHN mineral vector Mie

### 6.1 source-informed reconstruction

handoff는 Ahn thesis의 PSD, 상대 복소굴절률, 효율·backscatter 도식을 디지타이징한 source-informed seawater Mie 복원을 제공한다.

핵심 계산계약:

```text
host medium      seawater
n_water          1.34
particle RI      seawater-relative complex RI
size parameter   pi D n_water / lambda_vac
PSD integration  logarithmic
internal angle   0.1 deg
output angle     0.5 deg, 361 points
matrix           P11/P12/P33
normalization    0.5 integral P11 dmu = 1
```

최종 bulk grid는 330–1100 nm이고 phase node는 16개이다.

### 6.2 과학적 caveat

원 Coulter bin table과 원 계산 front-end가 없기 때문에 exact historical reproduction은 주장하지 않는다. Ahn 도식의 Qa/Qb/Qbb/bb/b와 잔차가 남을 수 있으며, 이를 줄이기 위한 inverse tuning은 수행하지 않았다.

### 6.3 검증

4개 파일 모두 explicit range, finite values, `P11>=0`, `|P12|<=P11`, `|P33|<=P11`을 통과하였다.

## 7. Organic detritus candidate와 비채택 결정

### 7.1 handoff candidate

레시피:

```text
homogeneous sphere
dN/dD proportional to D^-4
D = 0.05–500 um
n_real = 1.04 relative
k(lambda) = 0.010658 exp(-0.007186 lambda_nm) * 0.75
n_water = 1.334
```

candidate는 330–1100 nm bulk와 22 phase nodes를 갖고 기본 파일-level 물리성 검사를 통과했다.

### 7.2 handoff의 조건부 상태

handoff 문서는 candidate가 legacy 550 nm phase와 `g`, `bb/b`, DoLP90 잔차를 가지므로 자동설치에서 제외하고 matched-input I/Q/U regression 이후 명시적으로 채택하도록 규정했다.

### 7.3 본 통합의 acceptance 결과

candidate를 current constituent path에 직접 넣어 443 nm Chl-only 조건을 시험하였다.

- L=200 복원위상 최소값: 약 `-8.070e-02`
- solver warning: negative reconstructed phase
- 축소 설정에서도 `rt_solver_sos_pol m=0` 실패
- 정상 RT 결과 생성 실패

따라서 candidate는 production에 부적합하다고 판정하였다.

### 7.4 최종정책

- validated legacy detritus 350–850 nm를 active로 유지
- candidate는 candidates directory에 provenance와 함께 보존
- 350–850 nm 밖 Chl>0 constituent request는 fail-loud
- Chl absorption-only 표는 330–1100 nm이지만, 전체 Chl constituent 지원범위를 absorption 자료만으로 주장하지 않음

이 결정은 향후 다음 중 하나가 완료되면 재검토할 수 있다.

- higher-order/value-kernel-compatible phase treatment
- validated truncation policy
- source recipe 재조정 후 matched-input I/Q/U 통과
- 별도의 detritus vector phase reference 확보

## 8. Atmospheric aerosol Mie

### 8.1 모델군

- SnF/OPAC 16
- Ahmad2010 paper family 80
- Ahmad2010 AccuRT family 80

총 176개이다.

### 8.2 endpoint 확장

기존 RI nodes를 유지하면서:

```text
330 nm  : 350–400 nm 두 점 선형 외삽
1100 nm : 860–1240 nm 선형 내삽
```

을 추가하였다. n과 k는 독립 처리된다.

### 8.3 검증

176개가 모두 330/1100 bulk 및 phase endpoint를 갖고 C/Python 파일이 동일하다. 대표 C50 aerosol coupled smoke는 330/1100 nm에서 finite I/Q/U를 산출하였다.

## 9. Six-gas cross sections

### 9.1 계약

```text
H2O, O2, CO2, CH4, O3, NO2
330–1100 nm
1 nm, 771 points
40 AFGL levels
cm2 molecule-1
```

350–1100 canonical 부분은 보존하고 330–349 nm를 추가하였다.

### 9.2 성분별 처리

- H2O: HITRAN/HAPI Voigt 계열
- O2, CO2, CH4: 해당 UV 범위에 유효 전이선이 없음을 확인한 explicit zero
- O3: Serdyuchenko-Gorshelev 온도자료, 음수 측정잡음 0 제한, 1 nm boxcar 및 층온도 보간
- NO2: Vandaele 220/294 K 상대분광형을 각 층의 canonical 350 nm 값에 정합

### 9.3 제한

O2-O2/O4 CIA는 포함하지 않았다. 따라서 “six-gas contract complete”와 “모든 대기흡수과정 완전”은 같은 표현이 아니다.

## 10. Loader와 보간 계약

### 10.1 C

- central spectral header defines closed interval
- each loader checks monotonic finite axis and coverage
- target/ref wavelengths must be inside both bulk and phase axes
- exact endpoints are not extrapolated
- range errors stop the run with component/path information

### 10.2 Python

- `spectral_contract.py` provides global and table-query guards
- NumPy `interp` is called only after an in-range check
- Mie wavelength axes in um and environmental axes in nm are distinguished in diagnostics

### 10.3 Component-gated semantics

A global envelope does not override narrower active components. A requested simulation is valid only when every selected component covers that wavelength. This is why Chl=0 TSM/CDOM/aerosol runs support 330–1100 while Chl>0 currently fails outside 350–850.

## 11. C/Python paired reproducibility

- handoff paired payload equality: 197/197
- accepted runtime installation: 196/196 exact in C/Python
- detritus candidate archive: exact in C/Python and handoff
- scalar IOP matched-input max relative difference: `5.49e-15`
- native LUT parity at 330 and 1100 verified independently in C and Python

## 12. 논문작성 권고

논문에서는 각 스펙트럼 구간을 다음 네 종류로 구분해 표기해야 한다.

1. direct source-derived
2. source-node interpolation/regridding
3. project closure/extrapolation
4. conditional or rejected candidate

특히 PLOPS 750–800 taper, PLOPS >=800 zero, psi_T shortwave zero, AHN extrapolated RI/scalar, aerosol endpoint RI, NO2 scaled UV shape를 source measurement와 같은 범주로 기술하면 안 된다.
