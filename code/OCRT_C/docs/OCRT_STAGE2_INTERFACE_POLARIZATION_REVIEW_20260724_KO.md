# OCRT–OSOAA 2단계 공기→물 편광 투과 분석 검토 보고서

문서 기준일: 2026-07-24  
검토 대상: `files (32).zip`의 `ANALYSIS_interface_polarization_2026-07-22.md`, `gate5_2026-07-22.csv`, `interface_testbed_2026-07-22.csv` 및 부속 그림  
OCRT 기준 소스: `OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix`

## 1. 종합 판정

첨부 분석의 **핵심 분석방향은 타당하며, 본 세션의 후속 우선순위로 채택한다.** 특히 다음 가설은 제공 CSV의 재계산과 OSOAA `TAW` 이진행렬–OCRT 투과 커널 직접 대조에 의해 지지된다.

> SZA 40° 순수해수 gate에서 방위각 출력 결함을 제거한 뒤 남은 편광 잔차는 공기→물 rough-interface 투과 연산자의 Fourier/자오면 기준 조립 경로에서 최초로 나타난다.

그러나 첨부 문서의 다음 강한 표현은 현재 증거범위를 초과한다.

- “남은 오차가 전부 계면 하나에 있다.”
- “대기는 완벽하고 수중도 완벽하다.”
- “OSOAA 쪽에 버그가 없다.”
- level-26 water-minus-black `Rrs(0+)`가 유일한 정본이라는 판정

첨부 시험은 SZA 40°, VZA 60°, 풍속 3 m/s, 5개 파장, 에어로솔 없음에 국한된다. 본 세션의 가장 큰 잔차가 집중된 SZA 80°는 포함되지 않았다. 따라서 현재의 정확한 판정은 **“국소화 가설이 강하게 지지되며, 전역 원인 판정은 아직 미완료”**이다.

## 2. 첨부 CSV의 독립 재계산

### 2.1 순수해수 gate

`gate5_2026-07-22.csv`의 40행은 준위별 20셀로 구성된다.

- 파장: 412, 443, 490, 555, 660 nm
- SZA: 40°
- VZA: 60°
- OCRT RAA: 0°, 45°, 90°, 135°
- 준위: `0plus`, `0minus`

첨부 보고서가 사용한 `mean(|ΔX|)/max(|X_OSOAA|)` 및 `max(|ΔX|)/max(|X_OSOAA|)`를 재계산한 결과는 보고서의 순수해수 gate 수치와 일치하였다. 따라서 gate CSV와 해당 표의 기본 계산은 재현된다.

별도로 각 셀의 OSOAA I로 정규화한 `RMS(ΔQ/I)`와 `RMS(ΔU/I)`도 산출하였다. 이는 Q/U가 0에 가까운 셀에서 MAPE가 발산하는 문제를 피하면서, 본 세션의 150셀 결과와 동일한 정의로 비교하기 위한 것이다.

재생성 그림:

- `attached_gate_0plus_I_scatter.png`
- `attached_gate_0plus_Q_scatter.png`
- `attached_gate_0plus_U_scatter.png`
- `attached_gate_0minus_I_scatter.png`
- `attached_gate_0minus_Q_scatter.png`
- `attached_gate_0minus_U_scatter.png`

### 2.2 파장 의존성

`rrs(0-)`의 첨부 정규화 평균 오차는 412→660 nm 순으로 급감하였다.

- 412 nm: I 1.118%, Q 4.516%, U 6.008%
- 443 nm: I 0.589%, Q 2.354%, U 3.130%
- 490 nm: I 0.221%, Q 0.684%, U 0.935%
- 555 nm: I 0.049%, Q 0.081%, U 0.122%
- 660 nm: I 0.003%, Q 0.004%, U 0.006%

이 패턴은 해당 좁은 gate에서 잔차가 Rayleigh 하늘광 기여와 함께 증가한다는 해석을 지지한다. 다만 파장과 순수해수 흡수·산란계수도 동시에 변하므로, “Rayleigh optical depth에 정확히 비례한다”는 표현은 회귀계수와 독립축 실험 없이 확정하지 않는다.

### 2.3 강흡수 계면 시험대

`interface_testbed_2026-07-22.csv` 52셀을 재계산하였다.

- I의 최적 배율 OCRT/OSOAA: 0.999999990
- 배율 제거 후 I 형상 MAE: 0.00676%
- 배율 제거 후 I 형상 RMSE: 0.00952%
- Q/I RMSE, `max|OSOAA Q/I|` 정규화: 1.178%
- U/I RMSE, `max|OSOAA U/I|` 정규화: 1.707%

따라서 **세기 형상은 약 0.01% 수준으로 일치하지만 편광 비는 약 1–2% 수준으로 다르다**는 핵심 관찰은 유효하다.

첨부 문서 및 그림에 적힌 U/I 편향 −2.30%, RMSE 2.34%는 제공 CSV 전체 52셀을 동일 정의로 계산한 값과 일치하지 않는다. U/I의 음·양 zenith branch는 반대 부호의 편향을 가지므로 전체 signed mean은 거의 0이며, 전체 MAE와 RMSE는 각각 1.588%, 1.707%이다. 따라서 첨부 패키지는 일부 그림·본문·CSV가 서로 다른 분석 스냅샷 또는 지표를 혼합한 것으로 판단된다.

재생성 그림:

- `interface_testbed_I_scatter.png`
- `interface_testbed_Q_over_I_scatter.png`
- `interface_testbed_U_over_I_scatter.png`

## 3. OSOAA TAW–OCRT 커널 직접 대조

### 3.1 목적

첨부 문서가 제안한 다음 작업을 본 세션에서 실제로 수행하였다.

- OSOAA `TAW-1.340-03.0-RadMU48-NB200-SZA40.000-TSZA28.665` Fortran unformatted binary를 Fourier mode별로 해독
- OCRT `surface_T_aw_coxmunk_fourier_kernel()`을 동일 각도망에서 직접 호출
- m=0…4, 3×3 전 성분, 48×48 active Gauss node를 비교

### 3.2 입력 및 규약

대조 조건은 다음과 같이 고정하였다.

```text
wind_speed = 3.0 m/s
n_water = 1.34
sigma_type = 1  (sigma^2 = 0.003 + 0.00512*wind)
q_convention = 1
n_phi = 1024
m = 0,1,2,3,4
active angular nodes = 48/51
```

OSOAA Step 7과 OCRT contraction 식의 정규화를 직접 맞추었다.

```text
K_OSOAA^m = 2 * TAW^m / (C_m * mu_water * mu_air)
C_0 = 2*pi
C_m = pi  (m>0)
```

공개 방위각의 180° 이동에 해당하는 `(-1)^m`을 한 번 적용하였다. SZA40 전용 TAW 파일과 기본 TAW 파일의 active 48×48 부분행렬이 동일함도 확인하여, 추가 zero-weight solar/view node가 비교를 오염시키지 않음을 확인하였다.

### 3.3 스칼라 성분

M11의 OCRT/OSOAA 최적 배율은 m=0…4에서 0.999988–1.000003이다. 최대 절대차를 각 mode의 최대 기준값으로 정규화하면 0.011–0.091%이다.

이는 다음 항목이 상호 일치한다는 강한 근거이다.

- TAW binary decoder와 축 순서
- Gaussian slope density와 wind variance
- rough-interface 기하 및 공통 Jacobian
- 주된 scalar Fresnel transmission factor
- Fourier normalization과 공개 RAA 이동

따라서 첨부 분석의 “세기 공통인자는 정상이고 편광 조립을 우선 점검해야 한다”는 방향은 타당하다.

### 3.4 cosine 편광 블록

M12와 M21은 형상 상관이 0.9988 이상으로 높지만 약한 체계적 배율 차이를 보였다.

- M12: OCRT/OSOAA 최적 배율 약 1.010–1.011
- M21: OCRT/OSOAA 최적 배율 약 1.017–1.018
- M22, M33: 약 0.995–0.998

이는 계면 시험대의 Q/I 약 1% 차이와 정성적으로 일치한다. 단순 scalar Fresnel 오류라면 M11도 같은 비율로 이동해야 하므로, 현재 패턴은 편광 회전 또는 편광 성분별 조립 경로를 더 강하게 지시한다.

### 3.5 sine/U 결합 블록

m>0에서 다음 비대칭이 확인되었다.

- M23: 형상·부호·크기가 거의 일치
- M32: 크기는 거의 같으나 부호가 반대
- M13, M31: 단일 배율 또는 단일 전역 부호로 설명되지 않는 형상 차이

따라서 문제를 “U 관련 네 성분을 모두 뒤집으면 된다”로 단순화할 수 없다. `phi -> -phi`, 공통 U basis 반전, 두 회전각의 동시 반전은 모두 네 sine 성분을 일괄 변환하므로 M23의 이미 좋은 일치를 파괴한다.

실제로 OCRT 회전행렬의 입력측·출력측 sine 부호를 각각 또는 동시에 반전한 진단 빌드를 비교하였다. 일부 성분은 개선되었으나 다른 성분 및 M22/M33가 악화되어 어느 변형도 3×3 전체를 동시에 복원하지 못하였다. 첨부 세션에서 임시 부호 패치를 철회한 결정은 타당하다.

### 3.6 배제된 후보

- `n_phi=128`과 `1024`의 최대 차이는 전체 mode/성분에서 0.000768% 이하이다. 방위각 적분 해상도는 현재 잔차의 원인이 아니다.
- `q_convention=0`은 M12/M21의 평균 정규화 RMSE를 약 12% 수준으로 악화시킨다. 현재 비교에는 `q_convention=1`이 적합하다.
- 단순 입력측/출력측 회전 sine 부호 반전은 전 성분 동시 해를 제공하지 않는다.

## 4. 첨부 주장별 판정

### 4.1 승인하는 내용

1. 강흡수 수중매질로 수중 산란·내부반사 오염을 억제한 계면 전용 시험대는 적절하다.
2. level-26과 level-27을 직접 비교하여 오차가 최초로 나타나는 경계를 찾는 방식은 적절하다.
3. OSOAA TAW binary와 OCRT Fourier kernel의 mode-by-mode 3×3 직접 비교가 다음 우선순위라는 제안은 정확하다.
4. 임의의 전역 부호 패치를 먼저 적용하지 말아야 한다는 결론은 타당하다.
5. 풍속 0, 입자, 에어로솔, SZA 80°가 미검증이라는 자기 제한은 맞다.

### 4.2 표현을 제한해야 하는 내용

1. “모든 남은 오차”라는 표현은 SZA 40° gate 범위를 넘는다.
2. 대기장 level-26 일치가 해당 조건의 대기 solver를 지지하지만, SZA 80° 및 전체 geometry에서 대기가 완전하다는 증거는 아니다.
3. 수중 산란을 강흡수로 억제한 시험은 계면 국소화에는 유효하지만, 실제 순수해수 다중산란 solver 전부가 무오류임을 증명하지 않는다.
4. OSOAA의 local Fresnel 계수가 해석식과 맞는다는 사실은 generalized-spherical/Fourier 조립 전체가 무오류임을 증명하지 않는다.
5. 어느 코드가 물리적으로 옳은지는 reciprocity, energy/étendue, flat-interface limit 및 독립 Jones-basis 계산으로 판정해야 한다.

### 4.3 Rrs 산출 경로

첨부 문서는 level-26 full-water minus black-water 차감을 정본으로 사용한다. 본 세션의 기존 기준선은 level-27 수중 Fourier field를 별도 TWA로 전달한다. 두 방법은 동일 20셀에서 Rrs I/Q/U에 유의한 차이를 보이며, 특히 Q는 두 큰 항의 차감으로 얻어져 조건수가 나쁘다.

따라서 현 단계에서는 다음과 같이 관리한다.

- level-27 + 별도 TWA: **잠정 기준 branch**
- level-26 water-minus-black: **독립 진단 branch**
- 두 branch의 최종 우선순위: TWA/TAW 전달 연산자 및 Snell–Bouguer convention 검증 후 확정

방법론 문서 역시 수중장과 계면 전달을 분리하여 기술하고, Snell–Bouguer factor의 구현 해석을 open cross-check로 남겨 둔다. 따라서 첨부 방법을 즉시 정본으로 교체하지 않는다.

## 5. 입력 완전 일치 관점에서 추가로 확인한 사항

### 5.1 해수 깊이

첨부 문서의 “200 m cap은 맑은 물에서 충분히 깊지 않을 수 있다”는 경고는 타당하다. 현재 OCRT에서 pure water 443 nm, SZA40, VZA60, RAA90을 깊이 cap별로 실행한 결과:

- 200 m, tau=2.388: 깊은 기준 대비 Rrs I −0.714%, rrs I −0.699%
- 500 m, tau=5.971: 깊은 기준 대비 약 −0.0005%
- 1000 m, tau=11.941: 깊은 기준 대비 약 −0.00002%
- adaptive no-cap: tau=20, z=1674.8 m

따라서 **200 m가 부족할 수 있다는 진단은 맞지만 “최소 tau=30”은 제공 근거보다 강하다.** 이 사례에서는 tau≈6부터 사실상 수렴하였다.

OSOAA 기준 실행은 `SEA.Depth=1000 m`이고 OCRT 기존 비교는 adaptive tau=20을 사용하였다. 출력 차이는 무시 가능한 수준이나 사용자의 입력 완전 일치 원칙상 향후 공식 gate에서는 양쪽 모두 1000 m로 명시한다.

### 5.2 누락된 메타데이터

첨부 `gate5`에는 다음 값이 명시적으로 기록되지 않았다.

- OCRT `q_convention`
- OCRT `sigma_type`
- OCRT 투과 kernel `n_phi`
- 사용한 `OCRT_DUMP_SKY_FULL` 진단 패치의 정확한 소스

현재 배포된 OCRT 수정본에는 `OCRT_DUMP_SKY_FULL`이 없고 `OCRT_DUMP_SKY`만 존재한다. 따라서 첨부 계면장 dump는 제공된 소스만으로 완전 재현되지 않는다. 후속 gate는 소스 SHA, binary SHA, 진단 patch SHA와 모든 옵션을 manifest에 포함해야 한다.

## 6. 본 세션의 적용 결정

### 6.1 생산 물리코드

**이번 검토만으로 생산 물리코드를 수정하지 않는다.** 이유는 다음과 같다.

- 불일치가 TAW 편광 Fourier block에 존재함은 확인되었다.
- 그러나 M23은 이미 맞고 M32만 반대이며 M13/M31은 단순 부호 문제가 아니다.
- 어느 구현이 물리적으로 옳은지 아직 독립 불변량으로 판정되지 않았다.
- 임의 패치는 특정 gate를 개선하면서 다른 성분·기하를 악화시킬 위험이 높다.

따라서 현재 생산 바이너리와 버전은 유지한다.

```text
version: OCRT-v1.2-2026-07-23-KST-stage2-water-raa-output-fix
binary SHA-256: 41391518435c788c4d99f798c82fd39a016c8e4227ce220cedc2096f001838df
```

### 6.2 진단 코드 및 작업지시서

다음 항목은 즉시 본 세션에 적용한다.

1. OSOAA TAW decoder와 OCRT kernel 직접 대조 도구 추가
2. 입력 manifest에 `q_convention`, `sigma_type`, `n_phi`, water depth 추가
3. SZA40 interface gate와 SZA80 stress gate 분리
4. Rrs 두 산출 branch를 동시에 보존
5. 후속 패치 전 contraction-closure 및 독립 물리 불변량 gate 의무화
6. 세션 시작 소스부터 이어지는 누적 작업지시서 갱신

진단 도구는 생산 solver 호출경로에 연결되지 않으므로 OCRT 계산시간에는 영향이 없다. 생산 바이너리는 byte-identical이다. 별도 TAW audit 자체의 48×48, m=0…4, n_phi=1024 전체 실행시간은 약 5.83초이며, 이는 solver runtime이 아니라 독립 검증 도구의 비용이다.

## 7. 후속 분석 순서

### 7.1 P0 — 재현계약 고정

- 첨부 세션의 정확한 source/binary/diagnostic patch 확보
- water depth 1000 m 양쪽 고정
- q convention, sigma law, n_mu, m_max, n_phi, dpol, IOP 명시
- 모든 지표식을 manifest에 기록

### 7.2 P1 — contraction closure

각 Fourier mode와 단일 Gauss node에 대해 입력 Stokes basis vector를 주입한다.

```text
[1,0,0]^T
[0,1,0]^T
[0,0,1]^T
```

OSOAA Step 7과 OCRT coupling contraction 결과를 직접 비교하여, TAW 행렬 차이가 level-27 Q/U 차이를 정량적으로 재현하는지 확인한다. 그 다음 실제 level-26 Rayleigh BOA field를 입력하여 전체 계면 시험대 잔차를 닫는다.

### 7.3 P2 — 물리 불변량 판정

- flat-interface limit
- wind→0 연속성
- TAW/TWA reciprocity
- unpolarized flux energy conservation
- étendue/radiance n² law
- principal-plane U=0
- independent Jones-vector basis 계산

이 단계에서 OCRT 또는 OSOAA 중 어느 쪽의 sine/reference-plane 조립을 수정할지 결정한다.

### 7.4 P3 — SZA 80° 전이

SZA40에서 확정한 변환 또는 코드 수정이 SZA80 66셀의 Q/U 잔차를 실제로 감소시키는지 확인한다. SZA80 잔차가 남으면 대기 직달 경로, pseudo-spherical 여부, 고 SZA angular resolution을 별도 축으로 분리한다.

### 7.5 P4 — 패치 승인

물리 판정이 끝난 뒤 최소 변경만 적용한다. 변경 후 다음을 모두 재실행한다.

- 첨부 20셀 interface gate
- 순수해수 150셀 stress gate
- 전체 600셀
- Rrs/rrs I/Q/U 독립 산포도 6개씩
- mirror identity, TOA invariance
- 실행시간 1회 예열 + 5회 교차 측정

평균 또는 중앙 실행시간이 5% 이상 증가하면 hotspot profiling을 수행한다.

## 8. 결론

첨부 분석은 **후속 연구 방향으로 채택할 가치가 충분하다.** 본 세션의 독립 행렬 감사는 scalar TAW가 사실상 일치하는 반면 polarized Fourier block, 특히 sine/U 결합 성분에서 구조적 차이가 존재함을 확인하였다. 따라서 다음 작업을 수중 SOS나 Fresnel 진폭 전체에 넓게 분산시키지 않고, 공기→물 투과의 Fourier/reference-plane assembly와 contraction closure에 집중하는 것이 적절하다.

다만 현 시점에서 즉시 적용 가능한 것은 **진단·검증 체계의 업데이트**이며, 생산 물리코드의 부호 또는 회전식 수정은 아직 승인하지 않는다.
