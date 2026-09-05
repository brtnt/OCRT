# OCRT v1.2 PSSA 수치 오류 수정 및 영향 분석

작성일: 2026-07-19  
실행 식별자: `OCRT-v1.2-2026-07-19-KST-pssa-numeric-fixes`

## 1. 수정 범위

He et al. (2018)의 PCOART-SA는 대기에서 하향 태양 직달빔과 평면 해수면 반사 직달빔의 감쇠 및 단일산란 source를 구면기하로 계산하고, 그 이후 고차산란은 locally plane-parallel로 계산한다. OCRT는 다음 네 수치 문제를 순차적으로 수정했다.

1. 반사 직달빔의 lower endpoint 평면평행 감쇠 부호
2. Rayleigh=0, aerosol=0, gas>0인 해양 결합 PSSA helper의 조기 평면평행 반환
3. 반사빔 기원 downward first-order field의 고차 SOS 결합 누락
4. 반사빔의 층별 국소 방향 `beta`를 phase source에 적용하지 않던 문제

수중 굴절 직달빔 PSSA는 사용자 결정에 따라 구현하지 않는다. `--pssa` 적용 범위는 TOA에서 0+까지이며 0− 이하 수중 RT는 평면평행이다.

## 2. 수정 내용

### 2.1 Lower endpoint 감쇠

TOA에서 아래 방향으로 측정한 optical depth를 `h`, 전체 대기 optical depth를 `tau_s`라 하면 평면 해수면에서 반사된 직달빔이 level `h`까지 통과한 optical path는

\[
2\tau_s-h
\]

이다. 따라서 평면평행 reflected-beam amplitude는

\[
B(h)=\frac12\exp\left[-\frac{2\tau_s-h}{\mu_0}\right]
\]

이다. 기존 lower endpoint의 `exp(-h/mu0)`를 `exp(+h/mu0)`로 수정했다. PSSA가 활성화된 경우에는 수치적으로 더 안정하게 baseline×correction 비를 곱하지 않고

\[
B_{SA}(h_k)=\frac12\exp[-\xi_{refl}(k)]
\]

를 직접 평가한다. 이 방식은 강한 흡수대·고 SZA에서 발생할 수 있는 `0*inf`를 제거한다.

### 2.2 Gas-only ocean helper

기존 helper는 gas absorption을 layer에 넣기 전에 `tau_R+tau_a==0`이면 평면평행 값을 반환했다. 이제 산란이 0이어도 US Standard Atmosphere 1962 고도 grid를 만들고 AFGL gas optical depth를 layer별로 적분한 뒤 `rt_pssa_apply()`를 호출한다.

산란 source가 0인 것과 direct-beam extinction이 0인 것은 별개로 처리한다.

### 2.3 Reflected downward order-1 결합

반사 직달빔의 단일산란으로 생성된 upward 및 downward field를 모두 `order1_*`와 누적 field에 넣는다. 이후 SOS source operator가 두 방향을 동일하게 소비하므로 별도 반복 loop나 추가 solve는 없다.

수정 전에는

```text
reflected beam -> upward order 1 -> TOA
```

만 존재했고,

```text
reflected beam -> downward order 1 -> order 2+ -> TOA/surface
```

경로가 빠져 있었다.

### 2.4 Layer-local beta phase source

He et al. Eq. (8)의 관계

\[
\frac{\sin\alpha_k}{\sin\beta_k}=\frac{R+H_k}{R},
\qquad 2\alpha_k=\beta_k+\theta_0
\]

에서

\[
\beta_k=2\alpha_k-\theta_0
\]

를 계산하여 `rt_atm_t::pssa_beta[k]`에 저장한다. 반사빔 source는 더 이상 고정 태양방향 basis를 사용하지 않고 각 layer endpoint의 `beta_k`에서 scalar 및 spin-2 generalized angular basis를 평가한다.

모든 beta basis와 P11/P12/P33 moment contraction은 Fourier mode별로 SOS order loop 진입 전에 계산한다. 순수 Rayleigh에서는 l=2 항만 계산하여 setup overhead를 줄였다.

## 3. 기하·수치 불변조건

`tests/test_pssa_geometry.c`는 400-layer, SZA 85°에서 다음을 검사한다.

- Eq. (8) 비율 최대 상대잔차: `2.2204396324461602e-16`
- `2 alpha - beta - theta0` 최대 절대잔차: `0`
- `ch[k] - 0.5 exp(-xi_dn[k])` 최대 절대잔차: `0`
- surface에서 `alpha=beta=theta0`
- surface에서 `xi_refl=xi_dn`

 arbitrary-direction basis evaluator는 기존 quadrature-table builder와 raw-bit identical임을 별도 unit test로 확인했다.

## 4. 결과 차이

### 4.1 Rayleigh + flat Fresnel + PSSA

조건 행렬:

```text
72 cases
SZA       = 70, 75, 80, 85 deg
wavelength= 412, 443 nm
VZA       = 0, 30, 60 deg
RAA       = 0, 90, 180 deg
layers    = 100/100/200/400
absorbing gas = off
```

아래 수치는 각 수정만 되돌린 ablation과 최종 코드를 비교한 최대 변화다. Q/U가 0을 지나는 경우 단순 상대오차가 발산하므로, 편광은 절대차와 `|delta Q,U|/|I|`도 함께 제시한다.

| 수정 | max `|delta rho_I|` | max `|delta rho_I/rho_I|` | max `|delta Q|` | max `|delta Q|/|I|` | max `|delta U|` | max `|delta U|/|I|` |
|---|---:|---:|---:|---:|---:|---:|
| endpoint 부호 | 6.61698e-3 | 1.11957% | 2.37969e-3 | 0.49779% | 4.01455e-4 | 0.09840% |
| downward order-1 결합 | 9.73842e-3 | 1.72091% | 3.33648e-3 | 0.76855% | 6.13789e-4 | 0.15423% |
| local beta | 2.44829e-4 | 0.03992% | 2.47254e-4 | 0.03953% | 4.36136e-4 | 0.07954% |
| 전체 수정 대 이전 패키지 | 9.82533e-3 | 1.73843% | 3.47660e-3 | 0.75964% | 8.42747e-4 | 0.20751% |

`T_diff_dn_hemi`의 최대 변화는 endpoint 6.54%, downward order-1 8.88%, 전체 8.88%였다. 이는 빠져 있던 아래방향 1차장이 확산 하향장에 직접 기여하기 때문이다.

### 4.2 Aerosol + Rayleigh + flat Fresnel + PSSA

조건 행렬:

```text
12 cases
M80C, AOD555=0.15, wavelength=443 nm
SZA=75,80,85; VZA=0,60; RAA=0,90
layers=100/200/400; gas off
```

| 수정 | max `|delta rho_I/rho_I|` | max `|delta Q|/|I|` | max `|delta U|/|I|` | max `|delta T_diff_dn_hemi/T|` |
|---|---:|---:|---:|---:|
| endpoint 부호 | 0.67913% | 0.27932% | 0.02845% | 3.23601% |
| downward order-1 결합 | 1.00846% | 0.42862% | 0.05418% | 3.81477% |
| local beta | 0.00877% | 0.00756% | 0.01754% | 0.00586% |
| 전체 수정 대 이전 패키지 | 1.01217% | 0.42663% | 0.06407% | 3.81458% |

### 4.3 Plane-parallel flat Fresnel

endpoint 부호와 downward order-1 결합은 PSSA 전용 문제가 아니라 평면평행 flat reflected-beam source에도 존재했다. 27개 PP 조건에서 전체 수정의 `rho_I` 최대 변화는 1.48113%, `T_diff_dn_hemi` 최대 변화는 6.13639%였다.

black surface와 wind>0 Cox–Munk rough-surface control은 이전 패키지와 stdout/stderr byte-exact였다. 이 네 수정은 해당 경로를 건드리지 않는다.

### 4.4 Gas-only ocean PSSA

24개 조건(SZA 70/75/80/85°, 412/555/670/687/760/820 nm)에서 `T_dir_dn` 변화는 파장과 흡수 포화도에 따라 달랐다.

| 조건 | PP/fallback `T_dir_dn` | 수정 PSSA `T_dir_dn` | 변화 |
|---|---:|---:|---:|
| 412 nm, SZA 85° | 0.979557587 | 0.985673704 | +0.6244% |
| 555 nm, SZA 85° | 0.688187138 | 0.761423966 | +10.6420% |
| 670 nm, SZA 85° | 0.837510219 | 0.878596008 | +4.9057% |
| 687 nm, SZA 85° | 7.69664e-6 | 2.79026e-5 | +262.53% |
| 820 nm, SZA 85° | 0.0230063 | 0.0268969 | +16.9111% |

687/760 nm 같은 포화대는 분모가 매우 작아 상대변화가 크게 보인다. 760 nm, SZA 85°에서는 `7.31e-52 -> 2.17e-46`이지만 두 값 모두 사실상 0이다. 절대값과 포화 여부를 함께 해석해야 한다.

`Rrs0plus_I`는 gas-only helper 수정 전후 동일했다. 직달광량의 분자·분모에 같은 인자가 들어가 상쇄되는 현재 순수해수 조건의 결과다.

### 4.5 PSSA correction magnitude 갱신

flat Fresnel, gas off, 400 layers, nadir의 기존/수정 PSSA correction은 다음과 같다.

| wavelength | SZA | 이전 | 수정 |
|---:|---:|---:|---:|
| 412 | 70° | 0.347651% | 0.358203% |
| 412 | 75° | 0.854337% | 0.892984% |
| 412 | 80° | 2.556747% | 2.674818% |
| 412 | 85° | 11.069548% | 11.293166% |
| 443 | 70° | 0.257797% | 0.264268% |
| 443 | 75° | 0.685831% | 0.717856% |
| 443 | 80° | 2.260429% | 2.382414% |
| 443 | 85° | 10.869931% | 11.258206% |

이 표는 각 버전의 자체 PP 결과를 기준으로 계산했다. endpoint/downward-order 수정이 PP와 PSSA 양쪽에 영향을 주므로 서로 다른 버전의 PP 분모를 혼용하면 안 된다.

## 5. 산란차수 확인

412 nm, SZA 85°, VZA 30°, RAA 90°, flat PSSA, 400 layers에서:

| SOS cap | `rho_I` |
|---:|---:|
| 1 | 0.22076059458 |
| 2 | 0.29703812927 |
| 20 | 0.34989871055 |

누락된 downward order-1 field는 order 1 TOA에는 직접 영향을 주지 않고 order 2부터 차이를 만든다. 이 특성이 단계별 A/B에서 확인됐다.

## 6. 성능

beta phase source는 setup 단계에서만 구성되며 SOS order loop에서는 기존 operator 적용만 수행한다. 순수 Rayleigh는 l=2 전용 경로를 사용한다.

`OMP_NUM_THREADS=1`, 이전 stage3와 교대 실행:

| 조건 | 이전 median | 최종 median | 비율 |
|---|---:|---:|---:|
| Rayleigh PSSA, 400 layers | 0.13593 s | 0.13446 s | 0.9892 |
| M80C PSSA, 200 layers | 1.57456 s | 1.57627 s | 1.0011 |

측정 가능한 성능 저하는 없다. aerosol 사례의 0.11% 차이는 실행 잡음 범위다. 추가 SOS solve, order, phase file I/O 또는 order별 동적할당은 없다.

## 7. 검증 자산

```text
tests/test_pssa_geometry.c
scripts/smoke_pssa_numeric_fixes.sh
validation/pssa_numeric_fixes_2026-07-19/
```

주요 CSV:

```text
PSSA_NUMERIC_FIXES_RAYLEIGH_MATRIX_2026-07-19.csv
PSSA_NUMERIC_FIXES_AEROSOL_MATRIX_2026-07-19.csv
PSSA_NUMERIC_FIXES_PP_FLAT_MATRIX_2026-07-19.csv
PSSA_NUMERIC_FIXES_GAS_ONLY_OCEAN_MATRIX_2026-07-19.csv
PSSA_CORRECTION_BEFORE_AFTER_2026-07-19.csv
```

## 8. 남은 범위

- 수중 PSSA: 의도적으로 미구현
- SZA 85° 부근 층수 민감도: 경고 유지; 자동 layer 변경 없음
- CDISORT/AccuRT 논문 수치 field 직접 재현: 미완료
- SZA 85° 초과: 미검증 외삽

따라서 현재 구현은 **He et al. (2018)의 대기 하향·평면해수면 반사 direct-beam PSSA를 OCRT SOS 구조에 통합한 구현**으로 기술할 수 있다. 수중 보정과 논문 benchmark 재현이 없으므로 전체 PCOART-SA paper-equivalence는 주장하지 않는다.
