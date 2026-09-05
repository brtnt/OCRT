# OCRT v1.2 PSSA 범위와 0+/0− 직달·확산 조도 분해

## 1. PSSA 적용 범위

OCRT의 `--pssa`는 **대기 중 태양 직달빔**에만 적용한다.

적용 범위는 다음과 같다.

```text
TOA → 대기층 → 해수면 직상 0+
```

대기에서 보정하는 항은 다음 두 가지다.

1. 태양 하향 직달빔의 구면 껍질 경사 광학경로
2. 평면 해수면에서 반사된 태양 직달빔의 대기 왕복 경사 광학경로

고차 다중산란은 locally plane-parallel SOS를 유지한다.

### 1.1 수중 PSSA는 의도적으로 구현하지 않음

He et al. (2018)은 대기에서 굴절된 태양 직달빔의 수중 경로에도 같은 구면기하를 적용할 수 있다고 기술한다. OCRT v1.2에서는 이를 **의도적으로 구현하지 않는다**.

```text
0+ → Snell/Fresnel → 0− → 수중 RT
                         └─ plane-parallel 유지
```

따라서 `--pssa`의 적용은 0+에서 종료되고, 굴절 후 수중 태양 직달빔은 기존 평면평행 수중 RT로 전달된다. 이는 숨은 fallback이나 미완성 자동 분기가 아니라 명시적인 모델 범위 결정이다.

현재 OCRT가 대상으로 하는 수십–100 m 규모 수심은 지구반지름에 비해 매우 작으므로, 수중 곡률효과는 대기 PSSA 효과보다 훨씬 작을 것으로 판단한다. 다만 수중 PSSA를 실제로 구현·검증한 것은 아니므로, 본 패키지는 `PCOART-SA paper-equivalent full implementation`을 주장하지 않는다.

## 2. 새 조도 출력의 정의

이번 릴리스는 다음 여섯 변수를 출력한다.

```text
Ed0plus_direct
Ed0plus_diffuse
Ed0minus_direct
Ed0minus_diffuse
Eu0minus_direct
Eu0minus_diffuse
```

이 분해는 **표준 각도공간 radiometric direct/diffuse 분해**다.

- `direct`: 해당 경계면에 도달한 비산란 collimated solar beam
- `diffuse`: 그 외 모든 방향분포 복사장의 반구 적분

이는 광자의 최초 기원을 추적하는 `sun-origin / sky-origin` 분해가 아니다. 예를 들어 직달 태양광이 수중에서 한 번 이상 산란된 뒤 하향으로 진행하면 `diffuse`에 포함된다.

### 2.1 0+ 공기측

태양천정각을 \(\theta_0\), 공기측 태양 코사인을 \(\mu_0=\cos\theta_0\), BOA 직달 태양복사 플럭스를 \(F_{\mathrm{sun,BOA}}\)라 하면:

\[
E_d^{\mathrm{direct}}(0+)
=F_{\mathrm{sun,BOA}}\mu_0
\]

\[
E_d^{\mathrm{diffuse}}(0+)
=E_d(0+)-E_d^{\mathrm{direct}}(0+)
\]

출력 대응:

```text
Ed0plus          = 총 하향조도 E_d(0+)
Ed0plus_direct   = 비산란 태양 직달 하향조도
Ed0plus_diffuse  = 대기 확산 하향조도
```

PSSA가 켜지면 `F_sun_BOA` 계산에 대기 구면 직달경로가 반영된다.

### 2.2 0− 수중측 하향조도

공기–해수면 Fresnel 투과와 Snell 굴절을 적용한 수중 직달빔의 플럭스 스케일을 \(F_{\mathrm{sun,water}}\), 수중 태양 코사인을 \(\mu_{0w}\)라 하면:

\[
E_d^{\mathrm{direct}}(0-)
=F_{\mathrm{sun,water}}\mu_{0w}
\]

\[
E_d^{\mathrm{diffuse}}(0-)
=E_d(0-)-E_d^{\mathrm{direct}}(0-)
\]

출력 대응:

```text
Ed0minus          = 총 하향조도 E_d(0−)
Ed0minus_direct   = 해수면을 통과한 비산란 태양 직달 하향조도
Ed0minus_diffuse  = 수중 하향 확산조도
```

`Ed0minus_diffuse`에는 다음이 포함될 수 있다.

- 대기 skylight가 해수면을 통과한 성분
- 수중 산란으로 생성된 하향 성분
- 수면 내부반사 후 다시 하향하는 확산 성분

### 2.3 0− 수중측 상향조도

현재 OCRT는 깊은 수체와 black bottom을 사용하며, 수면 직하에서 상향하는 collimated source를 두지 않는다. 따라서:

\[
E_u^{\mathrm{direct}}(0-)=0
\]

\[
E_u^{\mathrm{diffuse}}(0-)=E_u(0-)
\]

출력 대응:

```text
Eu0minus          = 총 상향조도 E_u(0−)
Eu0minus_direct   = 0
Eu0minus_diffuse  = 총 상향 확산조도
```

향후 specular bottom 또는 별도 상향 collimated source를 추가하면 이 정의를 확장해야 한다.

## 3. 보존식

내부 double 값에서 다음 항등식이 성립한다.

\[
E_d(0+)
=E_d^{\mathrm{direct}}(0+)+E_d^{\mathrm{diffuse}}(0+)
\]

\[
E_d(0-)
=E_d^{\mathrm{direct}}(0-)+E_d^{\mathrm{diffuse}}(0-)
\]

\[
E_u(0-)
=E_u^{\mathrm{direct}}(0-)+E_u^{\mathrm{diffuse}}(0-)
\]

단일 실행 stdout은 일부 항목을 `%.6e`로 출력하므로, 텍스트를 다시 더하면 마지막 출력 자리 수준의 반올림 오차가 보일 수 있다. Full-grid CSV는 `%.12e`로 출력한다.

## 4. 성능 구현 원칙

여섯 변수는 이미 계산된 다음 값으로 조립한다.

```text
F_sun_BOA
F_sun_water
mu_sun_air
mu_sun_water
Ed0plus total
Ed0minus total
Eu0minus total
```

추가하지 않은 작업:

```text
추가 SOS solve
추가 산란차수
추가 각도 적분
추가 phase/moment 계산
추가 파일 I/O
동적 메모리 할당
```

따라서 계산 복잡도 증가는 상수 개수의 뺄셈·곱셈·구조체 저장뿐이다. 출력 문자열이 여섯 항목 늘어나는 비용 외에 RT solver hot path의 연산량은 변하지 않는다.

## 5. 검증 범위

- 기존 공개 출력과 수치 결과: 대표 5개 경로에서 byte-exact
- 단일 실행 3개 조건: 변수 존재 및 합 보존 PASS
- Full-grid CSV: 4개 셀에서 합 보존 PASS
- ASan/UBSan: 수중 순수해수 및 대기 PSSA 경로 finding 0
- 성능: 교대 실행 20회에서 평균 실행시간 차이 약 0.3%, 중앙값 차이 약 1.2%로 측정 잡음 범위

## 6. 별도 기능으로 남는 분해

다음 분해는 이번 변수와 다르며, 필요할 경우 별도 설계가 필요하다.

```text
직달 태양 기원 광자 vs 대기 skylight 기원 광자
각 기원의 모든 다중산란 후 기여
```

이를 엄밀하게 계산하려면 복사원별 독립 RHS 또는 dual solve를 유지해야 하므로, 현재의 zero-cost radiometric 분해와 달리 계산시간 및 메모리 증가가 발생할 수 있다.

## 7. 고태양천정각 PSSA 층수 경고

2026-07-19부터 `--pssa` 실행 시 다음 저해상도 조합에 대해 stderr 경고를 출력한다.

```text
75 <= SZA < 80 and n_layers < 100
80 <= SZA < 84 and n_layers < 200
84 <= SZA      and n_layers < 400
```

85도 부근에서는 412 nm US62 수렴시험을 근거로, 총 TOA radiance에는 200층 이상을 보수적으로 권장하고 PSSA correction magnitude 자체의 sub-percent 수렴에는 400~800층을 권장한다. 경고는 계산을 자동 변경하지 않으며 SOS 반복 내부에 추가 연산을 넣지 않는다. 상세 근거와 잔여 PSSA 이슈는 `PSSA_KNOWN_ISSUES_AND_LAYER_WARNING_2026-07-19.md`를 참조한다.
