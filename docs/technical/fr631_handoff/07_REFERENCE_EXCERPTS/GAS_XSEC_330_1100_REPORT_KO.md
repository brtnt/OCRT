# OCRT 기체 흡수단면적 330–1100 nm 생성·검증 보고서

## 1. 입력 원자료 검증

사용자가 반환한 `OCRT_HITRAN_RAW_330_1100_20260813T144408Z.zip`은 내부 검증보고서 기준 PASS였다.

- H2O: 226,887 HITRAN 전이선
- O2: 1,290 전이선
- CO2: 10,731 전이선
- CH4: 26,332 전이선
- O3: Serdyuchenko–Gorshelev 213–1100 nm, 193–293 K 11온도
- NO2: HITRAN/Vandaele XSC 220 K 및 294 K
- API key 저장: 없음

## 2. 최종 데이터 계약

- 파장: 330–1100 nm, 1 nm, 771점
- 대기격자: AFGL US Standard 1976, 40 level
- 단위: cm² molecule⁻¹
- C/Python 설치파일: byte-identical
- 기존 검증된 350–1100 nm 값: 6종 모두 정확히 보존
- 330–349 nm: 파일에 명시적 데이터행으로 추가; runtime endpoint hold 불필요

## 3. 330–349 nm 생성정책

### H2O

clean global-isotopologue HITRAN line list에서 UV 구간 9,456개 전이선을 추출하고 HAPI 1.3.0.0 Voigt, 0.05 cm⁻¹ 조건으로 계산했다. 350 nm 독립 overlap의 평균 상대차는 1.772e-04, 최대 상대차는 8.861e-04이다.

### O2, CO2, CH4

330–350 nm 구간과 25 cm⁻¹ 안전여유에 해당하는 전이선이 0개임을 clean line list에서 확인하였다. 따라서 330–349 nm를 명시적 0으로 저장했다.

### O3

Serdyuchenko–Gorshelev 원자료의 유한 음수 측정잡음을 0으로 제한한 뒤 ±0.5 nm boxcar와 층온도 선형내삽을 적용했다. 기존 canonical O3와 350–700 nm 평균 상대차는 5.832e-05이다. 최종파일은 신규 330–349 nm만 사용하고 350–1100 nm canonical 값을 그대로 유지한다.

### NO2

현재 canonical 350–1100 nm는 Bogumil/SCIAMACHY 5온도 계열이다. 반환 HITRAN raw package에는 Vandaele 220/294 K XSC만 포함되어 있으므로, 330–349 nm에는 Vandaele의 상대 분광형을 사용하고 각 층의 350 nm canonical 값에 단일 scale로 연결했다. scale 범위는 1.0099–1.0875이다. 이 정책은 값의 연속성을 보장하면서 350 nm 이상 검증결과를 변경하지 않는다. 파장별 임의 blending이나 코드분기는 사용하지 않았다.

## 4. 검증

- 구조: 6/6 PASS (40 × 771)
- finite/nonnegative: PASS
- canonical 350–1100 exact preservation: 6/6 PASS
- OCRT C reader: PASS
- OCRT Python reader: PASS
- C/Python 반환값 일치: PASS
- ZIP 무결성: PASS

## 5. 제한

- O2–O2(O4) collision-induced absorption은 현재 OCRT 6-gas 계약에 포함하지 않았다.
- NO2 330–349 nm는 canonical Bogumil 원 수치표 자체가 아니라 Vandaele 상대형 기반의 경계정합 확장이다. 향후 Bogumil 203/223/243/273/293 K 원파일이 확보되면 해당 20개 파장만 source-identical 방식으로 교체할 수 있다.
- 본 패키지는 설치 준비 후보이며 active OCRT tree는 아직 변경하지 않았다.
