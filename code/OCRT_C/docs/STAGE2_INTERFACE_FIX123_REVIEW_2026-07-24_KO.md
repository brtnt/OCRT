# OCRT–OSOAA 2단계 추가 코드 결함 검토 보고서

문서 기준일: 2026-07-24  
검토 대상: `files (34).zip`의 `ANALYSIS3_interface_closed_2026-07-22.md`, `gateB_2026-07-22.csv`, 그림 및 비교 진입점  
비교 기준 코드: `OCRT-v1.2-2026-07-24-KST-stage2-taw-rotation-fix`와 세션 시작 기준 소스  
검토 후 생산 후보: `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123`

## 1. 결론

첨부 분석은 두 부분으로 분리하여 판정하여야 한다.

첫째, 공기→물 rough-interface 투과 연산자에 제안된 **FIX1, FIX2, FIX3의 결합은 타당하며 생산 코드에 반영해야 한다.** 기존 세션에서 적용한 FIX1만으로는 회전각 산란각은 교정되지만, 실제 광자 진행방향에 대한 자오면 회전과 Fourier U 저장 규약이 완결되지 않는다. 독립 극한해, OSOAA TAW 원시행렬의 48×48 전 각도망, Fourier mode `m=0…4`, 3×3 전 성분 대조가 세 수정을 하나의 원자적 변경으로 지지한다.

둘째, 첨부가 후속 결함으로 지목한 `ocrt_add_diffuse_top_primary()`의 확산 하늘광 U 소스는 **유력한 잔여 원인 후보이지만, 제안된 단순 부호수정은 아직 생산 패치로 승인할 수 없다.** 해당 함수와 주 SOS 연산자 사이에 입사 U 열의 부호 구조가 다르게 보이는 지점은 확인하였다. 그러나 분자 U 입력계수 세 항을 단순 반전한 시험은 강흡수 단일산란 조건에서 U를 약 5.4배 과보정하고 Q 응답을 훼손하였다. 따라서 이 부분은 mode/node별 boundary-source contraction closure를 먼저 수행한 뒤 수정하여야 한다.

최종적으로 본 검토에서 생산 코드에 반영한 범위는 다음과 같다.

- FIX1: 실제 광자 진행방향에 대응하는 상대방위각의 π 이동
- FIX2: 공기측 입사와 수중측 출사에 대한 부호 있는 방향여현 적용
- FIX3: `m>0` TAW Fourier 연산자의 **입사 U 열** `M13`, `M23` 부호 저장
- 독립 극한 및 Fourier 저장 규약 회귀시험
- 기존 수중 공개 RAA 출력 수정 유지
- 확산상단 U 단순 부호반전은 명시적으로 제외

## 2. 첨부 주장별 판정

### 2.1 “FIX1, FIX2, FIX3이 함께 필요하다”

**확인됨.**

세 수정은 동일한 부호를 반복 보정하는 것이 아니라 서로 다른 층위를 교정한다.

- FIX1은 microfacet 기하에서 사용하는 away-from-interface 벡터와 실제 광자 진행방향 사이의 상대방위각 차이를 교정한다.
- FIX2는 Stokes 자오면 회전에 사용하는 입사·출사 방향여현의 부호를 교정한다.
- FIX3은 각도공간 Mueller 행렬 자체가 아니라, OCRT와 OSOAA의 sine-Fourier 계수 저장 및 일반 3×3 contraction 규약을 일치시킨다.

FIX2를 적용하면 물리적 각도공간 행렬의 U 출력행이 바로잡히지만, `m>0` sine 계수의 입사 U 열에는 convolution convention에 따른 추가 음수가 필요하다. 따라서 FIX3은 FIX2와 함께 적용되어야 하며, FIX3 단독 적용은 올바른 수정이 아니다.

### 2.2 “계면이 0.05% 이내로 닫혔다”

**핵심 주장은 강하게 지지되나, 첨부의 정확한 전 각도 0.05% 대표값은 원시 CSV 부재로 완전 독립 재계산하지 못하였다.**

첨부 ZIP에는 수정 후 전 각도 계면 시험대 원시 CSV가 포함되어 있지 않다. 다만 다음 독립 검증이 동일 결론을 지지한다.

- 단일 미시면 극한에서 6개 편광 성분이 해석값과 약 `3×10^-7…6×10^-7` 이내로 일치
- OSOAA TAW 48×48, `m=0…4`, 3×3 전 성분의 최악 mode RMS가 성분 최대 기준 0.020463% 이하
- 최악 mode 평균절대차가 0.002327% 이하
- 모든 성분의 최소 상관계수가 0.999995 이상
- 첨부의 대표 강흡수 단일기하 OCRT 출력이 수치적으로 재현됨

따라서 “계면 연산자 결함이 FIX1+2+3으로 종결된다”는 물리적 결론은 승인할 수 있다.

### 2.3 “남은 최종 오차는 확산 하늘광 U 하나이다”

**국소화 가설은 유력하나 ‘하나로 확정’되었다는 표현은 현재 증거를 초과한다.**

확산상단 경로의 U 응답이 의심스럽다는 관측은 재현된다. 또한 등가빔 fallback이 `beam_q=0`으로 설정되어 확산광 Q/U를 온전히 운반하지 못한다는 문서화된 제한도 확인된다.

그러나 다음 이유로 단일 원인 확정은 보류한다.

1. 첨부 자체가 Q 해석해를 닫지 못하였다.
2. 단순 분자 U 열 부호반전은 종단 결과를 개선하지 않고 과보정한다.
3. SZA 80°의 큰 Q/U 잔차는 첨부 검증범위에 포함되지 않는다.
4. 경계 주입 소스와 주 SOS operator 사이에는 Fourier mode, 상·하향 parity, 전역 U wrapper 부호가 동시에 존재하므로 항별 부호만 보고 수정할 수 없다.

따라서 올바른 후속 표현은 다음과 같다.

> 확산 하늘광의 편광 경계주입 경로는 다음 우선 분석대상이다. 다만 생산 수정은 실제 하향 Fourier field를 이용한 mode/node별 contraction closure 후 결정한다.

### 2.4 “최종 순수해수 gate는 거의 개선되지 않는다”

**확인됨.**

FIX1+2+3은 계면 연산자 자체를 물리적으로 교정하지만, 최종 Rrs/rrs에서는 기존 오차상쇄가 제거되어 일부 집합의 통계가 소폭 악화된다. 이는 수정이 틀렸다는 뜻이 아니다. 물리적으로 틀린 계면 연산자가 다른 미해결 오차를 일부 상쇄하고 있었음을 의미한다.

## 3. FIX1과 FIX2의 물리적 유도

`surface_aw_microfacet_geometry_()`의 입사·출사 벡터는 계면에서 각 매질 바깥쪽을 향하도록 구성된다. 이 표현은 미시면 법선과 Fresnel 입사각을 구할 때 유효하지만, Stokes 기준면 회전은 광자의 실제 진행방향을 사용하여야 한다.

공기에서 물로 진행하는 입사광과 수중 출사광의 실제 진행방향을 다음과 같이 두면,

\[
\mathbf{k}_i=(-\sin\theta_i,0,-\mu_i),
\]

\[
\mathbf{k}_o=(\sin\theta_o\cos\phi,
               \sin\theta_o\sin\phi,
              -\mu_o),
\]

두 방향의 내적은

\[
\mathbf{k}_i\cdot\mathbf{k}_o
=
\mu_i\mu_o-
\sin\theta_i\sin\theta_o\cos\phi
\]

가 된다. 기존 FIX1 이전 식의 `+ sin_i sin_o cos(phi)`는 이 진행벡터 내적과 일치하지 않는다.

코드에서는 다음과 같이 구현한다.

```c
const double cphi_r = -cphi;
const double sphi_r = -sphi;
const double mu_i_s = -mu_i;
const double mu_o_s = -mu_o;

double cosPsi = mu_i_s * mu_o_s + si * so * cphi_r;

double ci1 = (mu_o_s - mu_i_s * cosPsi) / (dnm * s1_safe);
double si1 = so * sphi_r / dnm;

double ci2 = (mu_i_s - mu_o_s * cosPsi) / (dnm * s2_safe);
double si2 = si * sphi_r / dnm;
```

이 구조는 OCRT 방법론 문서의 signed-μ 진행방향과 자오면 Stokes 회전 정의에 부합한다.

## 4. 독립 극한검사

공기측 입사각 40°, 수중 출사 천정각 `10^-4°`, 실제 진행방향 상대방위각 30°, 풍속 3 m/s, 굴절률 1.34에서 극한해를 계산하였다.

FIX1만 적용한 경우 `M23`, `M31`, `M32`의 부호가 해석값과 반대이다.

- `M23`: 해석 −0.836153, FIX1 +0.836153
- `M31`: 해석 +0.225494, FIX1 −0.225494
- `M32`: 해석 +0.866025, FIX1 −0.866025

FIX1+FIX2 적용 후 해당 세 성분의 절대오차는 각각 약 `3.54×10^-7`, `2.71×10^-7`, `2.71×10^-7`로 감소한다. `M21`, `M22`, `M33`도 `6.2×10^-7` 이하로 일치한다.

이 검사는 OSOAA 수치에 의존하지 않으므로 FIX2의 물리적 필요성을 독립적으로 판정한다.

원자료: `POLAR_LIMIT_FIX1_VS_FIX123.csv`

## 5. FIX3의 의미: U 행이 아니라 입사 U 열

OSOAA Step 7은 TAW를 일반 행렬-벡터 곱으로 축약한다.

\[
I_o=T_{11}I_i+T_{12}Q_i+T_{13}U_i,
\]

\[
Q_o=T_{21}I_i+T_{22}Q_i+T_{23}U_i,
\]

\[
U_o=T_{31}I_i+T_{32}Q_i+T_{33}U_i.
\]

따라서 Fourier storage convention을 공통화할 때 보정 대상은 `U_o` 행 전체가 아니라, `m>0`에서 sine 계수로 저장되는 **입사 U 열**이다. OCRT에서는 정규화 직후 다음 두 성분만 반전한다.

```c
if (m > 0) {
    T_pair[0*3 + 2] = -T_pair[0*3 + 2]; /* M13 */
    T_pair[1*3 + 2] = -T_pair[1*3 + 2]; /* M23 */
}
```

`M33`은 반전하지 않는다. `M31`, `M32`는 FIX2의 물리적 회전 교정으로 바로잡힌다.

이 결과는 이전 세션의 “OSOAA M31/M32 저장 부호를 공통물리 규약으로 행 단위 변환한다”는 해석을 수정한다. 운영 contraction 기준에서는 FIX3을 입사 U 열에 적용하는 것이 옳다.

## 6. OSOAA TAW 9성분 감사

검증 조건은 다음과 같다.

```text
OSOAA TAW: TAW-1.340-03.0-RadMU48-NB200-SZA40.000-TSZA28.665
active grid: 48 × 48
Fourier mode: m = 0, 1, 2, 3, 4
Mueller elements: 3 × 3 전체
n_phi: 1024
wind speed: 3 m/s
n_water: 1.34
sigma_type: 1
q_convention: 1
```

FIX1만 적용한 경우 M31과 M32는 거의 정확한 반대부호였다.

- M31: 최소 상관계수 −0.999996, 최악 최대오차 200.445%
- M32: 최소 상관계수 −1.000000, 최악 최대오차 199.9997%

FIX1+FIX2+FIX3 적용 후에는 전 9성분이 동시에 닫힌다.

- 전 성분 최악 평균절대차/성분최대: 0.002327% 이하
- 전 성분 최악 RMS/성분최대: 0.020463% 이하
- 전 성분 최악 단일셀 최대차/성분최대: 0.652976% 이하
- 전 성분 최소 상관계수: 0.999995 이상

최대 단일셀 차이는 공기측 `mu≈0.0163`의 극단적 grazing node에 위치하며, 전체 연산자 평균 및 RMS는 0.02% 이하이다.

## 7. 첨부 Gate B 재계산

첨부 `gateB_2026-07-22.csv`는 준위별 20셀이다.

```text
bands: 412, 443, 490, 555, 660 nm
SZA: 40°
VZA: 60°
RAA: 0°, 45°, 90°, 135°
wind: 3 m/s
water n_mu: 64
water m_max: 4
```

첨부와 동일한 `mean(|ΔX|)/max(|X_OSOAA|)` 지표를 재계산하면 다음과 같다.

- Rrs(0+) I: 0.5726%, Q: 2.1086%, U: 2.1607%
- rrs(0−) I: 0.3933%, Q: 1.5230%, U: 2.0310%

각 셀의 OSOAA I로 정규화한 RMS는 다음과 같다.

- Rrs(0+) Q: 1.5160%, U: 1.9270%
- rrs(0−) Q: 1.0478%, U: 1.8725%

이 결과는 계면 연산자 수정이 최종 순수해수 Rrs/rrs 잔차 전체를 제거하지 않는다는 첨부 결론과 일치한다.

산포도 파일:

```text
gateB_Rrs0plus_I_scatter.png
gateB_Rrs0plus_Q_scatter.png
gateB_Rrs0plus_U_scatter.png
gateB_rrs0minus_I_scatter.png
gateB_rrs0minus_Q_scatter.png
gateB_rrs0minus_U_scatter.png
```

## 8. 확산상단 U 가설의 검토

### 8.1 현재 경로

`ocrt_add_diffuse_top_primary()`는 외부 확산 하향 Fourier field의 I, Q, U를 수중 1차 산란 소스로 주입한다. 이 함수에는 다음 부호 계층이 동시에 존재한다.

- 입사 방향 `cc=-c`
- 상향/하향 출사 parity
- value-kernel table의 I/Q/U 성분 배치
- U 출력행에 대한 전역 `ru += -ch_c*(...)` wrapper
- molecular Rayleigh 열의 `xpl`, `xrl`, `xtl` 조합

따라서 함수 내부의 특정 U 계수만 주 SOS 연산자와 육안 대조하여 반전하면 중복 부호가 발생할 수 있다.

### 8.2 재현된 기본 결과

443 nm, SZA 40°, VZA 60°, RAA 90°, 풍속 3 m/s, 강흡수 `aDOM440=20 m^-1`, water `n_mu=64`, `m_max=4`에서 FIX123 기본 결과는 첨부와 일치한다.

```text
rrs I = 1.345482e-05
rrs Q = -1.125497e-06
rrs U = -4.169513e-06
```

확산상단 경로를 차단하면 다음과 같다.

```text
rrs I = 1.137670e-05
rrs Q = -8.611914e-07
rrs U = -4.111171e-06
```

따라서 확산광 기여는 다음과 같다.

```text
ΔI = +2.07812e-06
ΔQ = -2.64306e-07
ΔU = -5.83420e-08
ΔQ/ΔI = -0.127185
ΔU/ΔI = -0.028074
```

첨부의 해석 목표 `Q/I=-0.17270`, `U/I=+0.03945`와 비교하면 U 부호가 반대라는 관측은 재현된다.

### 8.3 단순 부호반전 시험

분자 입사 U 열의 세 계수를 단순 반전한 진단 빌드는 다음을 산출하였다.

```text
ΔQ/ΔI = +0.004356
ΔU/ΔI = +0.213157
```

이는 목표 U 응답 `+0.03945`보다 약 5.4배 크며, Q 응답도 `−0.127`에서 거의 0으로 붕괴한다. 종단 rrs U의 OSOAA 정합성도 개선되지 않는다.

따라서 다음 패치는 **금지**한다.

```text
ocrt_add_diffuse_top_primary의 molecular U-input 계수를 근거 없이 일괄 부호반전
```

### 8.4 등가빔 경로

등가빔 fallback은 `mu_sun_water_override>0`일 때 `atm.beam_q=0`으로 설정한다. 이 경로는 하나의 등가빔으로 확산 하늘광의 전체 angular Q/U Fourier field를 표현할 수 없으므로, 편광 정합성의 생산 기준으로 사용할 수 없다.

## 9. 72개 실행 및 600셀 재검증

FIX123 생산 후보로 기존 입력계약의 72개 물리 실행을 전부 재계산하였다.

```text
실행 성공: 72/72
각 실행 full-grid 행수: 24
water convergence: 전 실행 1
비교 셀: 600/600
```

### 9.1 전체 600셀

FIX1-only 대비 FIX123에서 다음 변화가 발생하였다.

- Rrs I MAPE: 2.1903% → 2.1978%
- Rrs Q RMS(ΔQ/I): 9.6507% → 9.8699%
- Rrs U RMS(ΔU/I): 7.2323% → 7.3925%
- rrs I MAPE: 2.0960% → 2.1073%
- rrs Q RMS(ΔQ/I): 9.6010% → 9.8236%
- rrs U RMS(ΔU/I): 7.1738% → 7.3361%

### 9.2 순수해수 150셀

- Rrs I MAPE: 1.8870% → 1.8942%
- Rrs Q RMS(ΔQ/I): 8.7358% → 8.9329%
- Rrs U RMS(ΔU/I): 6.5356% → 6.6795%
- rrs I MAPE: 1.8250% → 1.8341%
- rrs Q RMS(ΔQ/I): 8.6881% → 8.8884%
- rrs U RMS(ΔU/I): 6.4778% → 6.6238%

### 9.3 순수해수 SZA ≤ 40° 84셀

I는 약 0.002 percentage point 개선되지만 Q/U는 약 0.01–0.02 percentage point 악화되는 미소한 변화이다.

- Rrs I: 1.0063% → 1.0046%
- Rrs Q: 1.4795% → 1.4974%
- Rrs U: 1.1869% → 1.1981%
- rrs I: 1.1250% → 1.1231%
- rrs Q: 1.3399% → 1.3591%
- rrs U: 1.1209% → 1.1327%

### 9.4 순수해수 SZA = 80° 66셀

고 SZA에서는 기존 오차상쇄 제거로 Q/U가 약 0.22–0.30 percentage point 악화된다.

- Rrs I: 3.0078% → 3.0264%
- Rrs Q: 13.0636% → 13.3605%
- Rrs U: 9.7614% → 9.9786%
- rrs I: 2.7160% → 2.7390%
- rrs Q: 13.0102% → 13.3118%
- rrs U: 9.6834% → 9.9037%

이 결과는 FIX123의 롤백 근거가 아니다. 연산자 및 물리 극한이 명백히 개선되므로, 남은 고 SZA 결함을 독립적으로 찾아야 한다.

각 범위의 Rrs/rrs I/Q/U 산포도 6개는 `OCRT_STAGE2_INTERFACE_FIX123_72RUNS_20260724/figures`에 수록하였다.

## 10. 생산 코드 변경

### 10.1 `src/shared/surface.c`

`surface_T_aw_coxmunk_trig()`에 FIX1과 FIX2를 적용하였다.

`surface_T_aw_coxmunk_fourier_kernel()`에서 `m>0`의 M13/M23에 FIX3을 적용하였다.

### 10.2 생산 코드에서 제외한 변경

`src/rt_water_rt.c`의 `ocrt_add_diffuse_top_primary()`는 변경하지 않았다. 진단용 단순 U 열 부호반전은 최종 트리에서 제거하였고, 회귀 스크립트가 해당 식의 재유입을 차단한다.

### 10.3 버전

```text
OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123
```

생산 바이너리 SHA-256:

```text
4a2269a13f81f4112dfb351ec3cac3dc5548915a8913670efa6efe2aa79180c8
```

## 11. 회귀시험

다음 시험을 통과하였다.

```text
PASS: air-to-water signed-rotation and Fourier U-column regression
PASS: stage-2 interface FIX1+FIX2+FIX3 source audit
PASS: OCRT public RAA convention helpers
PASS: stage-2 water RAA source audit
PASS: water RAA output regression
72/72 physical runs complete
600/600 matched cells generated
```

신규 독립시험 `tests/test_surface_aw_interface_conventions.c`는 다음을 검사한다.

1. 극한 해석행렬 6성분
2. `m=1` 각도공간 적분과 저장 Fourier 커널의 비교
3. FIX3이 M13/M23에만 적용되는지 여부
4. 확산상단 U의 기각된 단순 부호반전이 소스에 재유입되지 않았는지 여부

## 12. 계산시간

443 nm, 순수해수, SZA 40°, VZA 60°, RAA 90°, water `n_mu=48`, `m_max=2`, 확산상단 주입 제외 조건에서 1회 예열 후 FIX1-only와 FIX123을 각각 5회 교차 실행하였다.

```text
FIX1-only mean   = 2.477739 s
FIX123 mean      = 2.476390 s
mean change      = -0.0544%

FIX1-only median = 2.455856 s
FIX123 median    = 2.479298 s
median change    = +0.9546%
```

변화는 측정 변동 범위이다. 신규 적분, 반복문, 동적 메모리 할당 또는 phase-function 평가가 추가되지 않았으므로 병목 프로파일링 조건인 5% 증가에 해당하지 않는다.

## 13. 후속 분석방향

### 13.1 확산상단 boundary-source contraction closure

각 `m`, 입사 node, 출사 node에 다음 basis를 주입한다.

```text
[1,0,0]^T
[0,1,0]^T
[0,0,1]^T
```

다음 세 출력을 직접 대조한다.

1. 주 SOS operator의 단일입사 contraction
2. `ocrt_add_diffuse_top_primary()`의 경계주입 source
3. 동일 하향 Fourier field를 사용하는 독립 Python/NumPy contraction

비교 시 다음 부호를 각 항별로 기록한다.

```text
incident direction parity
outgoing direction parity
cosine/sine Fourier storage
U-output wrapper minus
m=0 vs m>0
particle vs molecular block
```

### 13.2 실제 OSOAA 하향장 적용

OSOAA level-26 또는 수면 위 하향 Fourier field의 I/Q/U를 동일 quadrature에 매핑한 뒤, TAW 및 수중 단일산란 source를 순서대로 축약한다. 단일 equivalent beam으로 대체하지 않는다.

### 13.3 Q 해석 폐쇄

첨부 해석해는 직달광 Q를 재현하지 못하였다. 수중 depolarized Rayleigh의 F11/F12/F22/F33, Stokes 회전 순서 및 `q_convention`을 소스와 동일하게 구현한 뒤 Q를 다시 판정한다.

### 13.4 SZA 40°와 80° 병행

모든 후속 패치는 다음 두 gate를 동시에 통과하여야 한다.

- SZA 40° 강흡수 계면/단일산란 gate
- SZA 80° 순수해수 66셀 stress gate

SZA 40°만 개선하고 SZA 80°를 악화시키는 경험적 부호수정은 승인하지 않는다.

### 13.5 Rrs 산출 branch

level-27 + 별도 TWA와 level-26 water-minus-black branch를 계속 병행 보존한다. TWA/TWA reciprocity, flux/étendue 및 Snell–Bouguer convention이 완결되기 전까지 어느 한 branch를 유일한 정본으로 선언하지 않는다.

## 14. 최종 판정

첨부 분석은 공기→물 rough-interface의 추가 결함을 정확히 찾아냈으며, FIX2와 FIX3은 기존 FIX1-only 후보에 반드시 추가되어야 한다. 이 세 수정은 독립 물리 극한과 OSOAA TAW 전 성분 대조로 승인되었고 생산 후보 코드에 반영하였다.

반면 확산상단 U 문제는 아직 분석단계이다. 관측 및 국소화는 유효하지만 단순 분자 U 열 부호반전은 종단 검증에서 실패하였다. 따라서 해당 변경은 생산 코드에 넣지 않았으며, 다음 단계는 경계주입 source와 주 SOS operator의 mode/node별 contraction closure이다.
