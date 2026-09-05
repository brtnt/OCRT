# 330–1100 nm Detritus Mie 후보자료 채택 판정

## 1. 질문

handoff가 제공한 `Detritus_Stramski2001.mie` 330–1100 nm candidate를 현행 OCRT constituent production path의 새 canonical phase로 채택할 수 있는가?

## 2. handoff의 사전상태

handoff 문서는 이 파일을 조건부로 분류했다.

- default installer에서 제외
- legacy 550 nm 대비 `g`, `bb/b`, DoLP90 잔차 기록
- source-recipe regeneration이지 exact legacy reproduction은 아님
- 현재 본 개발 세션에서 명시적으로 채택하고 matched-input I/Q/U regression을 통과시킬 때만 설치

## 3. 파일단위 검사

candidate:

- SHA-256: `3c88824db36a098cf8b58e3d36f160ed7a44b61570ec46e99c9277470ac21a60`
- bulk range: 330–1100 nm
- phase range: 330–1100 nm
- finite: PASS
- P11 nonnegative on stored grid: PASS
- `|P12|<=P11`: PASS
- `|P33|<=P11`: PASS

이 단계만으로는 OCRT production solver 적합성을 보장하지 않는다.

## 4. coefficient representation 검사

OCRT constituent phase는 현재 주로 L=200 generalized-moment representation을 사용한다. Candidate를 water+organic mixture에 넣어 재구성한 결과 sampled back-angle 구간의 최소 P11이 약:

```text
-8.070e-02
```

이었다. 물리적으로 P11은 음수가 될 수 없다.

## 5. SOS 검사

443 nm, Chl=1 mg m-3, atmosphere off, rough ocean 조건에서 candidate를 active filename으로 임시 교체하여 실행했다.

축소 설정에서도:

```text
S-009: L=200 reconstruction negative
rt_solver_sos_pol FAIL m=0 rc=-2
rt_water_rt_sos_pure FAIL rc=-3
```

으로 종료했다. 정상적인 TOA/Rrs/rrs I/Q/U를 생성하지 못했다.

더 큰 설정에서는 300초 제한 전에 결과를 반환하지 못했다.

## 6. 최종 판정

```text
REJECT FOR STAGE-2 PRODUCTION
```

이유:

1. handoff가 이미 conditional로 지정
2. current L=200 solver representation에서 phase negativity
3. SOS failure/nonconvergence
4. matched-input I/Q/U acceptance 미통과

## 7. 최종 파일배치

Active:

```text
water_iop/Detritus_Stramski2001.mie
SHA-256 0202c29d7a3f050d88715766a37ef7080cbf709ef5a6e69c1ab645171e61f516
native phase range 350–850 nm
```

Archived candidate:

```text
water_iop/candidates/Detritus_Stramski2001_330_1100_CANDIDATE_REJECTED_L200.mie
```

C/Python에 동일하게 배치했다.

## 8. runtime 정책

- Chl<=0: detritus phase가 필요 없으므로 330–1100 nm 다른 component 계산 가능
- Chl>0, 350–850 nm: active validated detritus 사용
- Chl>0, 330–349 또는 851–1100 nm: fail-loud
- generic endpoint hold 사용 금지

## 9. 재검토 조건

다음 중 하나가 충족되면 candidate 또는 후속자료를 다시 평가할 수 있다.

- value-kernel 경로의 production 안정성 검증
- constituent multi-component truncation의 정식 구현·검증
- L order 확장과 수렴·성능 검증
- reference phase 또는 원 microphysics 확보
- C/Python/OSOAA matched I/Q/U 회귀 통과
