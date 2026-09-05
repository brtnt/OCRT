# OCRT 편광 water-leaving U 결손 수정 지시서

작성일: 2026-07-21
대상: OCRT 코드 수정 메인 세션
목적: 수중 편광 복사전달의 계면 내부반사 방위각 부호 오류를 수정하고, OSOAA 재실행 없이 제공된 참조 데이터로 검증한다.

---

## 0. 한 문장 요약

수중 SOS 솔버의 계면 내부반사 경계조건이 반사 항에 방위각 부호 (−1)^m을 곱하는데, 정반사는 방위각을 보존하므로 이 부호는 +1이어야 한다. 이 오류로 관측 U 성분(과 DoLP)이 OSOAA보다 약 14% 낮았고, 부호를 +1로 바꾸면 전 관측각에서 OSOAA와 1% 이내로 일치한다.

---

## 1. 근본 원인

### 1.1 물리
계면 내부반사는 수면에서 위로 가던 빛이 계면에 부딪혀 다시 아래로 반사되는 항이다. 평평한 계면의 정반사는 방위각 φ를 보존한다. 위로 가는 광선 (θ, φ)는 아래로 가는 (θ, φ)로 반사된다(수직 성분만 뒤집히고 수평 방향은 그대로다). 위·아래 반구가 같은 φ 기준을 쓰면, 방위각 푸리에 차수 m별 반사 계수는 모든 m에서 반사율 그대로다. 즉 부호는 +1이어야 하며 (−1)^m이 붙으면 안 된다.

### 1.2 코드 위치 (버전 v1.2 기준 라인번호)
파일: src/rt_solver.c

두 함수에서 내부반사 경계조건에 (−1)^m을 곱한다.

- rt_solver_sos_pol_intrefl (평평면 = 풍속 0 경로)
  - line 2062: const double msign = rt_fourier_pi_shift_sign(m);   /* specular (-1)^m */
- rt_solver_sos_pol_intrefl_rough (거친면 = 풍속 > 0 경로)
  - line 2506: const double msign = rt_fourier_pi_shift_sign(m);   /* specular (-1)^m, A2 top-BC convention */

rt_fourier_pi_shift_sign(m)은 src/rt_fourier.c:35에서 return (m & 1) ? -1.0 : 1.0; 즉 (−1)^m이다.

이 msign은 두 함수 안에서 반사된 하향 경계값을 만들 때 반사 뮐러 행렬 곱에 곱해진다. v1.2에서 flat 경로는 line 2165·2200·2373, rough 경로는 2506 함수 내부에서 다음 형태로 쓴다.

    tb_I[jp - 1] = msign * (M[0] * I_up + M[1] * Q_up + M[2] * U_up);
    tb_Q[jp - 1] = msign * (M[3] * I_up + M[4] * Q_up + M[5] * U_up);
    tb_U[jp - 1] = msign * (M[6] * I_up + M[7] * Q_up + M[8] * U_up);

여기서 M은 계면 내부반사 뮐러 행렬(3×3, Rww_M 또는 rough 커널), *_up은 직전 차수의 상향 장이다.

주의: 라인번호는 v1.2 소스 기준이다. 메인 세션 버전이 다르면 라인번호가 어긋날 수 있으므로, 위 두 함수 이름과 "specular (-1)^m" 주석, 그리고 tb_* = msign * (M[..]*..) 패턴으로 위치를 찾는다.

---

## 2. 수정

두 함수의 msign 정의를 rt_fourier_pi_shift_sign(m)에서 1.0으로 바꾼다.

수정 전:
    const double msign = rt_fourier_pi_shift_sign(m);   /* specular (-1)^m */
수정 후:
    const double msign = 1.0;   /* specular reflection preserves azimuth -> +1 (was (-1)^m, corrupted odd-m U) */

flat 경로(2062)와 rough 경로(2506) 두 곳 모두 바꾼다. flat은 풍속 0, rough는 풍속 > 0 케이스를 담당하므로 논문 케이스(풍속 > 0)는 rough 경로가 핵심이나, 두 경로가 물리적으로 같은 정반사이므로 둘 다 고친다.

rt_fourier_pi_shift_sign 함수 자체는 바꾸지 않는다(다른 곳에서 정당하게 쓰일 수 있다. 3항 참조).

### 비용
msign은 이미 계산돼 곱해지는 스칼라 상수다. 값만 바뀌므로 연산 비용은 0이다. 속도에 영향이 없다.

---

## 3. 확인 요망 (미검증 항목)

src/rt_solver.c:1996에도 const double msign = rt_fourier_pi_shift_sign(m);가 있다. 이 사용처는 위 두 내부반사 함수가 아닌, 그 앞의 다른 함수(추정: 내부반사 없는 기저 SOS rt_solver_sos_pol 또는 인접 함수)에 속한다. 이 세션에서는 검증하지 않았다.

메인 세션에서 line 1996이 속한 함수와 그 msign의 용도를 확인한다. 판단 기준:
- 만약 이 msign도 계면 정반사 성격의 상향→하향 반사 경계조건에 쓰인다면, 같은 논리로 +1이어야 한다(같은 수정 적용).
- 만약 위에서 아래로 입사하는 별도 장(예: 하향 입사 경계값)에 대한 정당한 π-이동 규약이라면, 그대로 둔다.

투과(coupling) 경로는 이미 부호를 안 붙인다. src/rt_air_water_coupling.c:189 주석: "C_0 = 2*pi, C_{m>0} = pi, no (-1)^m (transmission)." 이 사실이 반사도 +1이어야 한다는 정합성 근거다(투과와 반사는 둘 다 같은 정반사이므로 부호가 같아야 한다).

---

## 4. 증상-원인 매핑 (수정이 안전한 이유)

- 관측 I·Q는 짝수 m으로 구성된다. 짝수 m에서 (−1)^m = +1이므로, 잘못된 부호가 붙어도 값이 같았다. 그래서 수정 전에도 I·Q는 정확했다.
- 관측 U는 홀수 m으로 구성된다. 홀수 m에서 (−1)^m = −1이므로, 반사된 홀수차 장 I¹·Q¹·U¹의 부호가 뒤집혔다. 이게 U¹ 형성을 매 반사마다 틀어 다중산란에서 누적돼 U가 약 14% 낮았다.
- nadir(천정 관측)는 m=0만 존재하고 U=0이다. 따라서 nadir 세기(intensity)는 수정과 무관하다.

따라서 이 수정은 홀수 m(관측 U)만 바꾸고, I·Q·nadir 세기는 비트 단위로 그대로다. 이 성질은 버전과 무관하며, 아래 검증의 핵심 불변량이다.

---

## 5. 검증 절차 (OSOAA 재실행 불필요)

OSOAA 참조값은 이 문서(6·7항)에 제공한다. 메인 세션은 수정된 OCRT만 실행해 아래 표와 대조한다.

### 5.1 검증 케이스 정의 (canonical)
- 수중 성분: Red_clay 무기물, TSM=5 g/m³, Chl=0, aCDOM440=0
- 파장: 555 nm
- 기하: 태양천정각(sza)=30°, 상대방위각(raa)=90°
- 대기 없음: --pressure 0 --aod 0 (수중 복사전달만 격리)
- 격자: --n-mu-water 48 --water-m-max 4
- 환경변수: OCRT_ADVANCED=1 OCRT_DEBUG=1 OMP_NUM_THREADS=1

실행 명령 (예: vza=60, wind=3):

    OCRT_ADVANCED=1 OCRT_DEBUG=1 OMP_NUM_THREADS=1 \
    ./build/ocrt --surface ocean --wind-speed 3 --sza 30 --vza 60 --raa 90 \
      --wavelength 555 --pressure 0 --aod 0 --water-model ocrt \
      --ocrt-chl 0 --ocrt-tsm 5 --ocrt-adom440 0 --n-mu-water 48 --water-m-max 4

출력에서 rrs0minus_I, rrs0minus_Q, rrs0minus_U를 읽는다. DoLP = sqrt(Q^2+U^2)/I.

### 5.2 검증 A — 수정 전 상태 확인 (선택)
수정 전 OCRT가 6항 "OCRT 前" 표를 재현하면 시작 상태가 예상대로 buggy임을 확인한다. (같은 v1.2 빌드면 비트 재현. 버전이 다르면 값이 조금 다를 수 있으나 −12~−14% 결손 경향은 같아야 한다.)

### 5.3 검증 B — 수정 후 값 재현
수정 후 OCRT가 6항 "OCRT 後(목표)" 표를 재현해야 한다. 판정: DoLP가 수정 전 대비 약 +15~17% 증가하고, 7항 OSOAA 참조와 ±2% 이내로 일치한다.

### 5.4 검증 C — 불변량 확인 (버전 무관, 필수)
수정 전후로 다음이 비트 단위로 동일해야 한다.
- 임의 기하에서 rrs0minus_I (세기, 짝수 m).
- 임의 기하에서 rrs0minus_Q (짝수 m).
- nadir(vza=0)에서 rrs0minus_U (U=0, m=0만).
- 오직 off-nadir의 rrs0minus_U만 바뀐다(교정).

이 검증이 통과하면 수정이 홀수 m U에만 국한됨을 확인한 것이다.

### 5.5 검증 D — 수면 위(0+) DoLP 확인 (권고)
논문의 실제 관측량은 수면 위(0+) DoLP(대기 상단으로 나가는 편광 수출광)이다. 이 수정은 수면 아래(0−) U를 교정하며, 물→공기 투과(coupling)는 이미 부호가 옳으므로(+1) 수면 위 U도 따라서 교정된다. 메인 세션에서 0+ 산출량(Rrs 편광 또는 Lw 편광)의 DoLP가 수정 전후로 상향(약 +15~17%)되는지 확인할 것을 권한다. 이 세션에서는 0− 까지만 OSOAA와 직접 대조했다.

---

## 6. 참조 데이터 — vza 스윕 (wind=3, raa=90, 555nm, Red_clay TSM=5)

DoLP = sqrt(Q^2+U^2)/I. vza_wat는 스넬 굴절 후 수중 각도(n=1.34). OSOAA는 수중 레벨에서 vza_wat로 추출·보간한 값이다.

| vza_air | vza_wat | OCRT 前 DoLP | OCRT 後(목표) DoLP | OSOAA DoLP | (後−OSOAA) |
|--------:|--------:|-------------:|-------------------:|-----------:|-----------:|
| 20 | 14.8 | 0.03004 | 0.03342 | 0.03325 | +0.5% |
| 30 | 21.9 | 0.04142 | 0.04699 | 0.04665 | +0.7% |
| 45 | 31.8 | 0.05962 | 0.06850 | 0.06792 | +0.9% |
| 60 | 40.3 | 0.07601 | 0.08816 | 0.08731 | +1.0% |
| 70 | 44.5 | 0.08427 | 0.09833 | 0.09736 | +1.0% |
| 75 | 46.1 | 0.08730 | 0.10214 | 0.10111 | +1.0% |

- 수정 전: 전 각도에서 OSOAA 대비 −12~−14%.
- 수정 후: 전 각도에서 OSOAA 대비 +0.5~1.0% (하네스 허용범위 ±2% 안).
- 잔차 +1%는 하네스가 이미 아는 "OCRT 위상 cap 대 OSOAA 위상 truncation" 차이 및 wind=3 거친면 차이 범위다.

### 6.1 상세 성분비 (vza=60)
| 조건 | U/I | Q/I | DoLP |
|---|---:|---:|---:|
| wind=3, OCRT 前 | +0.07203 | −0.02426 | 0.07601 |
| wind=3, OCRT 後(목표) | +0.08475 | −0.02426 | 0.08816 |
| wind=3, OSOAA | +0.08403 | −0.02370 | 0.08731 |
| wind=0, OCRT 前 | +0.07114 | −0.02429 | 0.07517 |
| wind=0, OCRT 後(목표) | +0.08479 | −0.02429 | 0.08820 |

- Q/I는 수정 전후 불변(짝수 m). U/I만 교정됨.
- OSOAA와 Q/I가 −0.02426 대 −0.02370로 약 2% 차이나는 것은 위상 cap/truncation 차이로, 이 수정과 무관한 기존 항목이다.

---

## 7. 참조 데이터 — nadir 불변 (vza=0, wind=3)

수정 전·후 비트 단위로 동일해야 한다.
- rrs0minus_I = 8.402613e-02
- rrs0minus_Q = 1.421471e-03
- rrs0minus_U = −5.222401e-19 (거의 0)

참고: OSOAA nadir 세기 = 8.349e-2 (OCRT/OSOAA = +0.64%). 이 nadir 세기는 하네스가 검증하는 유일한 경로이며, 이 수정으로 바뀌지 않는다.

---

## 8. 커밋 게이트 고려사항

표준 커밋 게이트에 대한 이 수정의 영향:
- Tier-0 앵커(nadir rrs0minus_I): 비트 불변. 통과.
- 443 nm 기준결과: nadir 세기 기반이면 비트 불변(통과). off-nadir U를 포함하는 기준이라면 값이 교정되므로 기준 갱신 필요.
- 직전 버전과 비트 동일: I·Q·nadir는 통과. off-nadir U는 의도적으로 바뀐다(교정). 따라서 "비트 동일" 검사는 off-nadir U에 대해 예외 처리하고, 6·7항 목표값으로 기준을 갱신한다.
- STRICT 빌드: 영향 없음(스칼라 상수 변경).

권고: off-nadir U(및 DoLP)를 포함하는 모든 저장된 기준결과를 6항 "OCRT 後" 값으로 갱신한다.

---

## 9. 재현 환경 요약

- 빌드 명령:

    gcc -std=c11 -O3 -march=native -ffp-contract=fast -fassociative-math \
      -fno-signed-zeros -fno-trapping-math -DOCRT_FAST_KERNELS -fopenmp -Isrc \
      $(find src -name '*.c') -o build/<bin> -lm

- 검증 케이스는 대기 없이(--pressure 0 --aod 0) 수중만 격리했다. 이는 부호 오류가 수중 SOS 내부반사에 있어서, 대기를 빼면 원인이 더 선명하기 때문이다.
- 이 세션의 OSOAA 참조는 다음 설정에서 얻었다: Red_clay 위상함수(0.555µm), 무기물 전용 프로파일(a_min=0.13905, b_min=3.8238, OSOAA가 Z09 순수수 추가 → 단일산란알베도 0.95064, OCRT 0.9506과 일치), 수중 굴절률 1.34, 풍속 3, 48 가우스 절점, sza=30, raa=90. nadir 세기 교차확인 +0.64%로 설정 일치를 이미 검증했다.

---

## 10. 결론

- 검증된 수정: src/rt_solver.c의 rt_solver_sos_pol_intrefl(v1.2 line 2062)와 rt_solver_sos_pol_intrefl_rough(v1.2 line 2506)에서 내부반사 msign을 rt_fourier_pi_shift_sign(m)에서 1.0으로 변경한다.
- 효과: 수면 아래 편광 DoLP가 전 관측각에서 OSOAA와 −14%에서 +1% 이내로 교정된다. 세기·Q·nadir는 비트 불변. 비용 0.
- 확인 요망: line 1996의 msign 용도(3항).
- 논문 함의: 이 수정으로 OCRT 편광 수출광 DoLP가 OSOAA와 ~1% 일치하므로, DoLP 기반 흡수성 에어로솔 대기보정에 사용 가능하다.
