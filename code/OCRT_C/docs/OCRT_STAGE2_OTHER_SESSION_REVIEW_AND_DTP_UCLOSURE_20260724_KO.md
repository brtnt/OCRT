# 다른 세션 OCRT–OSOAA 분석묶음 검토 및 FIX4 독립 검증 보고서

문서 기준일: 2026-07-24  
검토 대상: `ocrt_osoaa_session_2026-07-22.tar.gz`  
검토 기준 소스: `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123`  
검토 후 생산 후보: `OCRT-v1.2-2026-07-24-KST-stage2-interface-fix123-dtp-uclosure`

## 1. 결론

첨부 분석은 중요한 원인 규명을 다수 포함하고 있으며, A–C 결함에 관한 결론은 기존 세션의 독립 검증과 일치한다. 가장 중요한 추가 판정은 후보 D, 즉 `ocrt_add_diffuse_top_primary()`의 분자 산란 입사-U 열이다.

첨부 문서는 코드 대조를 통해 `ruI`, `ruQ`, `ruU`의 부호가 주 SOS 연산자와 다르다고 진단했으나, 하나의 종단 비율 시험에서 U가 악화되었다는 이유로 “적용 금지”를 선언하였다. 본 검토에서는 첨부가 다음 단계로 요구한 mode–node–basis contraction closure를 실제로 수행하였다. 그 결과 기존 FIX123은 I 및 Q 입사 basis에서는 닫히지만, U 입사 basis의 m=1 및 m=2에서 각각 최대 0.0451267과 0.0223909의 절대 불일치를 보였다. 세 계수의 부호를 함께 교정한 뒤에는 순수 Rayleigh와 합성 Rayleigh–aerosol 혼합 모두 최대 `5.551×10^-17`로 닫혔다.

따라서 첨부의 **결함 진단은 맞고, “적용 금지”라는 최종 판정은 폐기하여야 한다.** FIX4는 OSOAA 수치에 맞추는 경험적 조정이 아니라 OCRT 내부의 동일 산란 연산자가 서로 다른 경로에서 동일 결과를 내도록 하는 필수 폐합 수정이다.

다만 첨부의 다른 주장 일부는 검증범위를 초과하거나 입력 완전일치 원칙을 충족하지 못한다. 특히 해수 깊이, 412 nm 수중 SOS 수렴, `Rrs(0+)` 추출 branch, SZA=80° 검증범위에 관한 정정이 필요하다.

## 2. 검토 방법과 증거 우선순위

판정은 다음 순서로 수행하였다.

1. 실제 패치와 생산 소스를 줄 단위로 대조하였다.
2. OSOAA와 무관한 독립 물리·수치 불변량을 우선 적용하였다.
3. 같은 mode, 입사 node, 출력 node 및 I/Q/U basis를 사용하여 boundary-source injector와 production SOS operator를 직접 축약하였다.
4. 동일 입력계약으로 72개 물리 실행과 600개 비교 셀을 다시 계산하였다.
5. 수정 전후 계산시간을 교차 측정하였다.

이 순서는 종단 오차가 다른 결함과 상쇄될 수 있다는 점을 고려한 것이다. 종단 평균이 나빠졌다는 사실만으로 독립 폐합을 만족시키는 물리 수정 자체를 기각하지 않았다.

## 3. 첨부 주장별 판정

### 3.1 결함 A — 수중 공개 RAA 및 U 재구성

**판정: 승인, 기존 수정 유지.**

수중 `Rrs(0+)`와 `rrs(0−)`의 공개 출력이 대기 출력과 다른 거울 RAA를 사용하던 결함 및 추가 U 부호반전은 이미 독립적인 거울변환 항등식과 water/TOA 일관성 검사로 확정하였다. 이 수정은 출력 기하 라벨을 바로잡으며, 실제 단일산란 전파벡터에서 필요한 `pi−RAA` 기하와 분리되어야 한다.

### 3.2 결함 B·C — 공기→물 rough-interface 회전 및 Fourier 저장

**판정: 승인, 기존 FIX1+FIX2+FIX3 유지.**

실제 광자 진행방향에 따른 방위각 π 이동, 부호 있는 입사·출사 방향여현, m>0 입사-U 열 저장부호는 단일 미시면 해석극한과 OSOAA TAW 전 9성분 감사로 독립 확인되었다. 세 변경은 원자적으로 적용되어야 하며, FIX3 단독 적용 또는 U 출력행 전체 반전은 금지한다.

첨부가 말하는 “계면 0.05% 종결”은 **SZA=40°, 강흡수 계면 시험대에서의 국소 판정**으로는 지지된다. 그러나 첨부 gate는 SZA=80°를 포함하지 않으므로 전역 Rrs/rrs 오차 종결을 의미하지 않는다.

### 3.3 결함 D — 확산상단 primary source의 입사-U 열

**판정: 코드 진단 승인, 첨부의 적용 금지 판정 철회, 생산 후보에 적용.**

수정 위치는 `src/rt_water_rt.c::ocrt_add_diffuse_top_primary()`의 분자 성분 세 항이다.

```c
const double ruI = ray_on ? -gamma2 * atm->xpl[j] * xtlc_u : 0.0;
const double ruQ = ray_on ? -alpha2 * atm->xrl[j] * xtlc_u : 0.0;
const double ruU = ray_on ? -alpha2 * atm->xtl[j] * xtlc_u : 0.0;
```

기존 `src_u` 조립에 이미 전역 음수가 있으므로, 동일한 downward incident U field를 주 SOS 연산자와 같은 규약으로 처리하려면 위 세 항이 음수이어야 한다.

## 4. U Fourier 규약의 소스 추적

`ext_top_U`가 별도의 반대부호 규약일 가능성을 소스 전체 경로에서 확인하였다.

- 대기 solver는 U를 sine Fourier 계수로 저장하며, BOA downward `total_u`를 `boa_export->U_per_m`에 그대로 복사한다.
- `rt_air_water_couple_atm_to_water()`는 이 값을 TAW 3×3 행렬과 통상적인 행렬–벡터 곱으로 축약하고, 결과 `Uw`를 `U_inwater_per_m`에 그대로 기록한다.
- `rt_solver.c`는 `U_inwater_per_m`을 `ocrt_dt_U`로 그대로 복사하고 `w_opts.ext_top_U`에 직접 연결한다.
- 이 사이에 별도의 U 부호변환은 없다.

따라서 `ext_top_U`는 수중 SOS field와 동일한 sine-Fourier Stokes-U 규약이다. 첨부가 남긴 두 번째 미확인 사항도 본 검토에서 닫혔다.

상세 추적은 `DTP_U_CONVENTION_SOURCE_TRACE.txt`에 수록하였다.

## 5. 기준코드 비의존 source/operator 폐합

시험은 각 Fourier mode와 입사 node에서 다음 basis를 각각 주입하였다.

\[
[1,0,0]^T,\qquad [0,1,0]^T,\qquad [0,0,1]^T.
\]

비교 대상은 다음과 같다.

1. `rt_sos_operator_apply_vector()`
2. `ocrt_add_diffuse_top_primary()`

두 경로의 Gauss 가중치와 깊이 감쇠는 동일하게 구성하였다. injector의 `AI=2wX`와 `ch_c=0.5 exp(-tau/mu)`의 곱이 SOS 연산자의 `wX exp(-tau/mu)`와 정확히 같아진다.

수정 전에는 U basis m=1에서 최대 0.0451267, m=2에서 0.0223909의 불일치가 있었다. 수정 후에는 다음과 같이 닫혔다.

- 순수 Rayleigh: 최대 `2.776×10^-17`
- Rayleigh–aerosol 혼합: 최대 `5.551×10^-17`
- 시험 mode: m=0…4
- 모든 입사 node 및 I/Q/U basis 포함

![DTP U source/operator closure](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/dtp_u_basis_source_operator_closure.png)

이 결과가 FIX4 승인에 대한 가장 강한 근거이다.

## 6. 첨부 하네스에서 정정해야 할 사항

### 6.1 해수 깊이는 두 코드에 동일하게 전달되지 않는다

첨부 `compare_entry.py`는 OCRT를 먼저 실행하여 IOP를 읽은 다음 `depth_m = 48/(a+b)`를 계산하고, 이 깊이를 OSOAA의 `-SEA.Depth`에만 전달한다. OCRT 명령에는 대응하는 물리 깊이 인수가 없다. 따라서 하네스가 “모든 입력을 명시했다”고 주장하더라도 해수 물리 깊이는 문자 그대로 동일하지 않다.

후속 공식 gate에서는 동일한 `depth_m`을 두 코드에 명시적으로 전달하고, 그 깊이에서 두 코드가 모두 수렴함을 별도로 확인해야 한다.

### 6.2 “최소 광학두께 30”은 보편 불변량이 아니다

반무한 수심을 보장하기 위한 보수적 guard로는 사용할 수 있으나, 물리적 필수값으로 선언할 근거는 부족하다. 본 세션의 443 nm 깊이 sweep에서는 약 `tau≈6`부터 결과가 사실상 수렴하였다. 정확한 원칙은 임의의 고정 tau가 아니라 다음 두 조건이다.

- 두 코드에 동일한 물리 깊이를 입력할 것
- 관심 출력의 depth convergence를 수치적으로 확인할 것

### 6.3 412 nm 수중 SOS가 충분히 수렴하지 않았다

첨부 gate는 412 nm 순수해수에서 기본 수중 SOS 상한 20으로 종료한 기록을 포함한다. 본 기준선에서는 같은 계열이 28–29차에서 수렴한다. 따라서 첨부의 412 nm 값은 방위각·계면 진단에는 사용할 수 있으나, 1% 정합성 최종 gate로 사용하려면 `--debug-water-max-orders 100` 등 수중 상한을 명시하여 재계산해야 한다.

### 6.4 Rrs(0+) 차감 branch는 진단용으로 유지한다

첨부는 level-26 full-water에서 black-water를 차감하여 `Rrs(0+)`를 만든다. 이 방식은 Q에서 두 큰 항의 차를 취하므로 조건수가 나쁘고, level-27 수중 Fourier 장 + 별도 TWA branch와 유의한 차이를 보인다.

따라서 두 branch를 병행 보존한다.

- level-27 + 별도 TWA: 잠정 기준 branch
- level-26 water-minus-black: 독립 진단 branch

최종 정본은 reciprocity, energy/étendue 및 Snell–Bouguer convention 검증 후 확정한다.

### 6.5 극 제약의 표현

첨부의 `Q₂=U₂`는 극에서의 spin-2 field에 대해 정의된 **짝지어진 m=2 Fourier 진폭**의 관계이다. 공간상 한 점에서 Stokes Q와 U가 항상 같다는 뜻은 아니다. 회귀시험에 유지하되 계수 정의와 정규화를 명시해야 한다.

## 7. 72개 실행 및 600셀 재검증

최종 후보 소스로 기존 기준선 입력계약을 그대로 사용하였다.

- 물리 실행: 72/72 성공
- 각 실행 full-grid 행수: 24
- 수중 수렴: 전 실행 `water_converged=1`
- 비교 셀: 600/600
- 412 nm 순수해수 수중 SOS: SZA 0°/40°/80°에서 29/29/28차

FIX4는 SZA=80° 편광 잔차를 크게 감소시켰으나, SZA≤40°와 I 성분의 일부 잔차를 악화시켰다.

전체 600셀에서 Q/U는 다음과 같이 감소하였다.

- Rrs Q: 9.870% → 3.668%
- Rrs U: 7.393% → 3.279%
- rrs Q: 9.824% → 3.450%
- rrs U: 7.336% → 3.144%

순수해수 SZA=80°에서는 다음과 같다.

- Rrs Q: 13.360% → 4.535%
- Rrs U: 9.979% → 3.954%
- rrs Q: 13.312% → 4.323%
- rrs U: 9.904% → 3.777%

반면 순수해수 SZA≤40°에서는 Q/U가 약 0.4–0.7 percentage point 악화되었고, I MAPE도 증가하였다. 이는 FIX4가 틀렸다는 뜻이 아니라 기존 U-column 결함이 다른 하류 오차를 부분 상쇄하고 있었음을 의미한다. 다음 분석에서는 저 SZA/I 잔차와 고 SZA 잔차를 독립적으로 분리해야 한다.

## 8. 순수해수 150셀 Rrs·rrs I/Q/U 산포도

### 8.1 Rrs(0+) I

![Rrs I](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/purewater150_Rrs0plus_I_scatter_dtp_uclosure.png)

I MAPE는 1.894%에서 2.277%로 증가하였다. 편광 폐합 수정과 별도로 intensity source 또는 전달 정규화 잔차가 남아 있음을 뜻한다.

### 8.2 Rrs(0+) Q

![Rrs Q](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/purewater150_Rrs0plus_Q_scatter_dtp_uclosure.png)

RMS(ΔQ/I)는 8.933%에서 3.340%로 감소하였다.

### 8.3 Rrs(0+) U

![Rrs U](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/purewater150_Rrs0plus_U_scatter_dtp_uclosure.png)

RMS(ΔU/I)는 6.680%에서 2.959%로 감소하였다.

### 8.4 rrs(0−) I

![rrs I](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/purewater150_rrs0minus_I_scatter_dtp_uclosure.png)

I MAPE는 1.834%에서 2.137%로 증가하였다.

### 8.5 rrs(0−) Q

![rrs Q](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/purewater150_rrs0minus_Q_scatter_dtp_uclosure.png)

RMS(ΔQ/I)는 8.888%에서 3.150%로 감소하였다.

### 8.6 rrs(0−) U

![rrs U](../OCRT_STAGE2_DTP_UCLOSURE_72RUNS_20260724/figures/purewater150_rrs0minus_U_scatter_dtp_uclosure.png)

RMS(ΔU/I)는 6.624%에서 2.824%로 감소하였다.

## 9. 계산시간과 병목

동일 443 nm 순수해수 조건에서 1회 예열 후 3회씩 교차 측정하였다.

- FIX123 평균: 2.649036 s
- FIX123+FIX4 평균: 2.640879 s
- 평균 변화: -0.3079%
- FIX123 중앙값: 2.653287 s
- FIX123+FIX4 중앙값: 2.640050 s
- 중앙값 변화: -0.4989%

신규 반복문, 적분 또는 동적 할당은 추가되지 않았고 계산시간 증가는 없다. 5% 병목 조사 기준에 해당하지 않는다.

## 10. 코드 반영 내용

생산 후보에 다음을 반영하였다.

1. 기존 water RAA/U output correction 유지
2. 기존 air→water FIX1+FIX2+FIX3 유지
3. diffuse-top molecular incoming-U `ruI/ruQ/ruU` 부호 교정
4. pure-Rayleigh 및 mixed Rayleigh–aerosol mode/node/basis closure test 추가
5. 기존 회귀 스크립트의 “FIX4 금지” 조건을 폐기하고 source contract 검사로 교체
6. 버전 문자열 및 README 갱신
7. 72/72 실행, 600/600 비교, 독립 산포도 및 runtime 자료 추가

생산 후보 바이너리 SHA-256:

```text
7e5b2938e3529cc068a70cd11514dfa22cbbc3563aa5ec8a910ce4a6aa258fd7
```

## 11. 다음 분석방향

FIX4 이후의 우선순위는 다음과 같다.

1. **저 SZA intensity/Q/U 잔차의 분리**: source term의 I/Q 열, Ed 정규화, water depth 및 Rrs branch를 독립축으로 시험한다.
2. **SZA 80° 잔여 3–4%의 order/mode 분해**: direct atmospheric source, diffuse BOA field, TAW contraction, water scattering order별 기여를 각각 비교한다.
3. **입력 완전일치 하네스 개정**: 물리 수심, 수중 SOS 상한, `n_mu`, `m_max`, surface matrix order, phase function 및 depolarization을 양쪽 코드에 명시한다.
4. **Rrs 두 branch 동시 출력**: level-27+TWA와 level-26 subtraction을 같은 행에 기록하여 차이를 별도 오차축으로 관리한다.
5. **입자 케이스 전 선행 gate**: OSOAA의 hydrosol truncation과 OCRT phase-function 처리를 완전히 일치시킨 뒤 Chl/TSM으로 확장한다.

## 12. 최종 판정

첨부 세션은 수중 공개 RAA, air→water 회전 및 diffuse-top U-column 결함을 정확히 찾아냈다. A–C는 이미 반영된 수정과 일치한다. D에 대해서는 첨부의 코드 진단이 맞았고, 첨부 스스로 제안한 직접 폐합시험을 수행한 결과 적용 필요성이 기계정밀도로 확정되었다.

따라서 본 세션에서는 FIX4를 누적 코드에 반영한다. 다만 첨부 하네스의 수심·수렴·Rrs 추출 범위는 공식 1% 정합성 판정에 그대로 사용할 수 없으며, 위 정정사항을 반영한 뒤 다음 단계로 진행하여야 한다.
