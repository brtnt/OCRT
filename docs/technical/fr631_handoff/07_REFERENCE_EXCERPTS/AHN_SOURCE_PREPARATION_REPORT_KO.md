# AHN 광물입자 해수 배경 Mie 파일 생성 보고서

## 1. 상태

- 대상: Brown earth, Yellow clay, Red clay, Calcareous sand
- 분광범위: **330–1100 nm**
- 상태: **검토용 생성 완료, OCRT active tree에는 아직 설치하지 않음**
- 계산방식: Ahn 논문에서 디지타이징한 PSD 및 상대 복소굴절률을 이용한 homogeneous-sphere forward Mie
- legacy phase와의 blending: **사용하지 않음**
- inverse Mie: **사용하지 않음**

## 2. 해수 배경 처리

사용자가 지적한 “배경이 공기가 아니라 해수” 조건을 다음과 같이 구현했다.

```text
m_rel(lambda) = n(lambda) - i n_imag(lambda)
x = pi D n_water / lambda_vacuum
n_water = 1.34
```

`n`, `n_imag`는 **해수에 대한 상대굴절률**이다. 입자 RI를 공기 기준 절대굴절률로 잘못 취급하지 않았고, `n_water=1.34`를 size parameter에 명시적으로 포함했다.

Mie phase 계산에는 “배경 질량밀도”가 직접 입력되지 않는다. 필요한 배경 물성은 광학적 매질 굴절률이다. 입자 또는 해수의 질량밀도는 number cross section을 질량특이계수로 새로 변환할 때 필요하지만, 이번 `.mie` bulk Ext/Sca에는 이미 승인된 OCRT Step-4의 mass-specific `a*`, `b*` 표를 사용했으므로 별도의 밀도 가정이 들어가지 않았다.

## 3. 입력자료 및 처리

### 3.1 굴절률

- 400–750 nm: 새 디지타이징 자료
- 330–399 nm: 승인된 단파장 선형외삽
- 751–1100 nm: 승인된 장파장 선형외삽
- Calcareous sand 허수부: 1100 nm에서 다른 세 종 평균으로 수렴
- 실수부와 허수부 모두 약 50 nm low-pass 적용
- 최종 2.5 nm 표를 PCHIP로 1 nm 격자화; PCHIP는 2.5→1 nm 내부 interpolation에만 사용

### 3.2 입경분포

- Figure 3-2를 상대 number distribution `dN/dD`로 해석
- equivalent-sphere diameter 사용
- 광물별 범위: Table 3.1의 `Dmin=0.3 µm`, `Dmax`
- Figure 3-2는 Junge-like 단조감소 분포이므로 디지타이징의 작은 역증가를 monotone-decreasing least-squares projection으로 제거
- 마지막 판독 가능점 이후 중앙정책: **0**
- 별도 민감도: 마지막 판독값에서 Table 3.1 `m`으로 `Dmax`까지 Junge continuation

monotone projection의 phase 영향은 매우 작았다. 최대 변화는 `|Δg|=1.442e-06`, `|Δbb/b|=0.0032%`, `|ΔDoLP90|=7.847e-05`이다.

## 4. Mie 계산규격

| 항목 | 적용값 |
|---|---|
| 입자형상 | homogeneous equivalent sphere |
| host medium | seawater |
| host RI | 1.34 |
| 파장 정의 | vacuum wavelength |
| particle RI | seawater-relative complex RI |
| PSD 중앙안 | Figure 3-2, zero tail |
| 입경 적분 | logarithmic diameter grid, Δlog10D≈0.0025 |
| 내부 각도격자 | 0.1° |
| 출력 각도격자 | 361개, 180°→0°, 0.5° |
| phase node | 330, 350, 400, 412, 443, 470, 488, 515, 550, 590, 633, 670, 694, 760, 860, 1100 nm |
| bulk grid | 330–1100 nm, 1 nm |
| 출력요소 | P11, P12, P33 |
| 정규화 | 0.5 × integral(P11 dµ) = 1 |
| bulk Ext/Sca | 승인된 OCRT Step-4 a*, b* |
| Asymm_Para | 신규 seawater Mie의 g |

## 5. 생성 파일

| 광물 | 파일 | bulk 행 | phase node | 각도 수 | 최대 정규화 오차 | 최소 P11 | SHA-256 |
|---|---|---|---|---|---|---|---|
| Brown earth | Brown_earth_AHN.mie | 771 | 16 | 361 | 4.470e-10 | 0.00464664 | 18439d8b8062b788b7d7359706c119063c4237f1aab7105d1606d458c737c7da |
| Yellow clay | Yellow_clay_AHN.mie | 771 | 16 | 361 | 2.195e-10 | 0.00405305 | b4a6aefa7cf9d80d9995e81d03ffb2489e40bfae0014f7dc030228aac69c8901 |
| Red clay | Red_clay_AHN.mie | 771 | 16 | 361 | 2.514e-10 | 0.0093284 | 48ac4bb17c1874174a125f3306ae493663bb19b9be221e936c5f271ca715000c |
| Calcareous sand | Calcareous_sand_AHN.mie | 771 | 16 | 361 | 4.314e-10 | 0.00581678 | db9d91e150c6794663f2a74c18ff74017fb15900f3b21d9f3357b39dc2948a3f |

C와 Python에는 동일한 파일을 복사해야 하며, 별도 변환본을 만들지 않는다.

## 6. 대표 phase 지표

| 광물 | λ (nm) | g | bb/b | -P12/P11 @90° | P33/P11 @90° |
|---|---|---|---|---|---|
| Brown earth | 330 | 0.931614 | 0.01259126 | 0.193572 | 0.771099 |
| Brown earth | 550 | 0.957706 | 0.00443844 | 0.639482 | 0.604844 |
| Brown earth | 1100 | 0.904527 | 0.01109923 | 0.789605 | 0.428237 |
| Yellow clay | 330 | 0.959352 | 0.00497112 | 0.469332 | 0.715228 |
| Yellow clay | 550 | 0.952451 | 0.00452178 | 0.766099 | 0.484216 |
| Yellow clay | 1100 | 0.879010 | 0.01484088 | 0.851479 | 0.340049 |
| Red clay | 330 | 0.936111 | 0.00864229 | 0.381311 | 0.720475 |
| Red clay | 550 | 0.911825 | 0.00989599 | 0.700518 | 0.521310 |
| Red clay | 1100 | 0.779487 | 0.03376146 | 0.858957 | 0.316555 |
| Calcareous sand | 330 | 0.954865 | 0.00541814 | 0.527109 | 0.682165 |
| Calcareous sand | 550 | 0.927604 | 0.00833168 | 0.651188 | 0.570436 |
| Calcareous sand | 1100 | 0.825989 | 0.02494453 | 0.824867 | 0.372144 |

## 7. 검증

### 7.1 구조 및 물리성

- 4개 파일 모두 OCRT 표준 validator 통과
- 771개 bulk row, 16개 phase node, 361개 각도 확인
- 모든 P11 양수
- `|P12|≤P11`, `|P33|≤P11` 충족
- 직렬화 후 최대 phase 정규화 오차: `4.470e-10`
- SSA 전 구간 0–1

### 7.2 실제 OCRT reader 및 실행

- OCRT C `read_mie_file`/interpolator smoke: PASS
- OCRT Python `read_mie` 및 PCHIP phase evaluation: PASS
- OCRT C CLI: 4종 × 330/750/1100 nm = **12/12 PASS**
- C/Python single-sphere + PSD integration parity 최대 절대차: `5.684e-14`

### 7.3 수치수렴

- 입경격자 fine 대비 기준:
  - max `|Δg| = 2.357e-04`
  - max relative `Δ(bb/b) = 1.251e-03`
  - max `|ΔDoLP90| = 2.510e-04`
  - max `|ΔP33r90| = 2.573e-04`
- 내부 각도격자 0.25° 대비 0.1°:
  - max `|Δg| ≈ 1.670e-05`
  - max relative `Δ(bb/b) ≈ 2.403e-04`
- 1 nm phase PCHIP:
  - finite/positive/Mueller bounds PASS
  - node 사이 P11 half-integral drift 최대 `0.169%`
- deterministic replay: **28/28 byte-identical**

## 8. 주요 불확실성

### 8.1 미판독 PSD tail

최종파일은 사용자 결정에 따라 zero-tail을 사용한다. 보수적인 Junge-tail sensitivity는 다음과 같다.

| 광물 | max |Δg| | max |Δbb/b| (%) | max |ΔDoLP90| | max |ΔP33r90| |
|---|---|---|---|---|
| Brown earth | 0.01476 | 14.3 | 0.04793 | 0.04254 |
| Calcareous sand | 0.03632 | 21.6 | 0.05196 | 0.04037 |
| Red clay | 0.04633 | 22.6 | 0.03879 | 0.03154 |
| Yellow clay | 0.02145 | 18.7 | 0.03685 | 0.03511 |

따라서 tail은 number fraction으로 작지만, 장파장 phase의 `g`와 `bb/b`에는 무시할 수 없는 영향을 줄 수 있다. 이 민감도표는 최종파일의 uncertainty 범위로 유지한다.

### 8.2 논문 도식과의 비교

새 forward-Mie 결과와 Figure 3-3~3-6 디지타이징의 차이는 아래와 같다.

| 광물 | 지표 | MAPE (%) | 평균 bias (%) |
|---|---|---|---|
| Brown earth | Qa | 35.2 | -35.2 |
| Brown earth | Qb | 30.0 | +29.4 |
| Brown earth | Qbb | 35.9 | -25.9 |
| Brown earth | bb_over_b | 43.6 | -43.6 |
| Calcareous sand | Qa | 24.3 | -24.3 |
| Calcareous sand | Qb | 18.0 | -17.6 |
| Calcareous sand | Qbb | 41.6 | -41.6 |
| Calcareous sand | bb_over_b | 27.4 | -26.2 |
| Red clay | Qa | 24.5 | -24.5 |
| Red clay | Qb | 24.6 | -20.1 |
| Red clay | Qbb | 16.7 | +16.7 |
| Red clay | bb_over_b | 50.7 | +50.7 |
| Yellow clay | Qa | 24.6 | -24.6 |
| Yellow clay | Qb | 15.7 | -9.5 |
| Yellow clay | Qbb | 60.0 | -60.0 |
| Yellow clay | bb_over_b | 54.9 | -54.9 |

일부 지표는 16–60% 수준의 차이가 남는다. 따라서 이번 결과를 **exact historical reproduction**으로 부르지 않는다. 정확한 명칭은 “Ahn source-informed seawater forward-Mie reconstruction candidate”이다. 이 차이를 줄이기 위해 RI 또는 PSD를 legacy 출력에 맞춰 임의 조정하지 않았다.

## 9. 최종 판정

- 해수 배경 Mie 계산계약: PASS
- 파일 구조·물리성·C/Python parity·수치수렴: PASS
- 330–1100 nm 전구간 phase 생성: 완료
- OCRT active 설치: 아직 미실시
- 최종 채택 전 검토 포인트: PSD tail sensitivity와 논문 도식 대비 잔차
