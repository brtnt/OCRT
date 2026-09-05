# OCRT 유기 성분 IOP 모델 명세

GOCI-III 편광 대기보정 논문을 위한 수중 유기 입자(식물플랑크톤 + 유기 detritus)의
고유광학특성(IOP) 산출 모델을 정리한다. 무기물(TSM)은 별도 성분이며 본 문서 범위 밖이다.
CDOM은 OCRT 독립 인터페이스가 담당하므로 본 문서에서 다루지 않는다.

본 문서는 OCRT 통합의 명세이며, 각 성분의 흡광 a·산란 b·위상함수 P를 어떻게
산출하고 병합하는지, 그 근거 레퍼런스가 무엇인지 기록한다.

---

## 1. 성분 구성

수중 입자는 세 성분으로 나뉜다. 중복 계산을 막기 위해 각 성분의 경계를 명확히 한다.

| 성분 | 모델 | 상대 굴절률 | 크기분포 | 구동 변수 | 상태 |
|---|---|---|---|---|---|
| 식물플랑크톤 | EAP 이중층 구(coated sphere) | shell 1.10 / core 1.02 | Bernard 2007 감마형 | Chl | 완료 |
| 유기 detritus | 균질 Mie | 1.04 | Junge slope 4, 0.05-500µm | Chl covary | 완료 |
| 무기물 TSM | Ahn 균질 Mie | 1.10-1.135 | Junge 0.3-few µm | TSM 농도(독립) | 완료(별도) |

식물플랑크톤과 유기 detritus는 Chl과 covary한다. 무기물은 Chl과 독립이다.
이 독립성이 Case-1(무기물 없는 외해)과 Case-2(무기물 있는 연안)를 나누는 정의이다.

---

## 2. 후방산란 bb 모델 — 방향 B

### 2.1 문제

순수 식물플랑크톤(EAP) 후방산란 bb는 세포 색소 흡광이 각인되어 파장에 따라
울퉁불퉁한 분광을 보인다(변곡 있음). 실제 해수에서는 Chl과 covary하는 유기 detritus가
후방산란을 지배하며, detritus bb는 넓은 크기분포에서 유래한 부드러운 멱함수(λ^-0.9)이다.
따라서 두 성분을 합하면 합산 bb 분광이 부드러워진다(관측 사실).

### 2.2 방향 A vs 방향 B — 방향 B 채택

방향 A("detritus = Huot 총 − EAP phyto 절대값")는 실패한다. EAP 이중층 구 bb가
실측 식물플랑크톤 비중(총 후방산란의 2-3%, Stramski 2004)보다 과대평가되기 때문이다.
EAP coated sphere는 균질 구보다 후방산란이 크게 나오는 알려진 특성이 있다.

방향 B(채택)는 절대 크기와 분광 형태를 분리한다.

```
총 절대 크기:   bbp(λ,Chl) = Huot 2008
분배 비중:      phyto 비율 f_ph(Chl), 나머지 detritus
분광 형태:      phyto → EAP 형태(각인 있음)
                detritus → Stramski λ^-0.9 형태
```

식물플랑크톤의 가치는 후방산란 절대량이 아니라 편광 위상함수(IQU)의 각도 구조이다.
따라서 절대 비중은 실측(Stramski 2004)으로, 각도·편광 형태는 EAP로 담당한다.

### 2.3 총 bbp — Huot 2008

Huot, Y. et al. (2008), Biogeosciences 5:495-507. 광범위 현장 관측(BIOSOPE 등)
기반 Case-1 총 particulate 후방산란의 Chl 함수. 고Chl까지 적용 검증됨.

```
bbp(λ,Chl) = α₁(λ) · Chl^β₁(λ)
  α₁(λ) = 2.267e-3 − 5.058e-6·(λ−550)   [m²/mg Chl]
  β₁(λ) = 0.565 + 0.000486·(λ−550)
```

이 모델은 bbp를 상수 배경 성분(서브마이크론 = detritus)과 phyto Chl-specific으로
분해한다. 배경 성분이 우리 detritus에 해당하므로 물리적으로 정합된다.

### 2.4 phyto 비율 f_ph — 연속 Chl 함수

Stramski et al. (2004), Prog. Oceanogr. 61:27-56. 개방양에서 식물플랑크톤 포함
미생물이 후방산란의 2-3%. 이를 f_ph의 기준으로 삼는다.

Chl 구간별 계단 함수는 배제한다. 이유는 고Chl에서 파장방향 거칠기가 2배 악화되기
때문이다(경계 불연속 12%는 f_ph 계단이 아니라 Huot의 정상 Chl^0.565 의존성으로 확인됨).

채택: 연속 함수.
```
f_ph(Chl) = 0.035 + 0.015·tanh(log10(Chl))
  저Chl(0.1) → ~0.020,  Chl=1 → 0.035,  고Chl(10) → ~0.050
```

거칠기가 전 Chl에서 phyto 단독 대비 1/25~1/40로 완화된다. 고Chl 잔여 각인은
무기물(TSM) 독립 시나리오를 얹으면 무기물의 부드러운 λ^-0.9(후방산란 30배 강)가 덮는다.

### 2.5 TSM covary 없음

무기물과 Chl은 covary하지 않는다(Case 정의). TSM은 Chl 함수가 아니라 독립 시나리오
파라미터(예: 0/1/5 g/m³)로 얹는다.

---

## 3. σbb 정규화 — wiscombe(Bohren-Huffman) 규약

단일입자 후방산란 단면적 σbb는 절대 진폭 규약이 필요하다. 위상함수 정규화 규약
(norm='one' 등)을 쓰면 λ² 인자가 어긋나 분광 기울기 부호가 뒤집힌다.

올바른 규약: miepython norm='wiscombe'가 정확한 BH 절대 진폭이다.
검증: ∫_0^π (|S1|²+|S2|²) sinθ dθ = qsca·x² (비 1.0000).

```
dCsca/dΩ = (|S1|² + |S2|²) / (2 k²),   k = 2π n_w / λ_vac  [µm⁻¹]
σbb = ∫_{90°..180°} dCsca/dΩ dΩ = (π/k²) ∫_back (|S1|²+|S2|²) sinθ dθ
```

이 정규화로 detritus 멱함수 γ=0.93(목표 0.90), 무기물(n=1.18) γ=0.93(목표 0.91),
무기물/detritus 배율 28.3배(목표 30배)를 재현한다. 위상행렬 요소:
```
P11 = (|S1|²+|S2|²)/2,  P12 = (|S2|²-|S1|²)/2,  P33 = Re(S1·S2*)
```

---

## 4. 식물플랑크톤 (EAP)

### 4.1 모델

Lain, Kravitz, Matthews, Bernard (2023), Sci Data 10:412. 이중층 구(core=세포질,
shell=엽록체). d'Milay(Toon-Ackerman 1981) Fortran 코드로 전 Mueller 행렬 산출.
공개 CSV는 스칼라 a/b/bb만 제공하나, d'Milay가 m1/m2/s21/d21을 반환하여 IQU 유도 가능.

EAP a_phyto는 순수 식물플랑크톤 흡광(in vivo 색소, detritus 배제)이다. 따라서
detritus 흡광을 별도로 더해도 중복되지 않는다.

### 4.2 굴절률 recipe

```
kshell_base: EAP_invivo_means.csv 그룹 컬럼 (np.interp)
kcore: .mat 파일 col0
Vs=0.2, Vc=0.8, FR=(1-Vs)^(1/3), ci=2e6 mg/m³ (2 kg/m³)
kshell_norm = (6.75e-7/nmedia)·(0.027·ci/Vs)/(4π)   [Johnsen 0.027 = 675nm 최대 Chl 흡광]
kshell = kshell_base · (kshell_norm / kshell_base[675nm])
nshell = 1.10 + Hilbert_imag(kshell),  ncore = 1.02 + Hilbert_imag(kcore)
nmedia = 1.334, wavelength = λ/nmedia, wvno = 2π/wavelength
```

### 4.3 크기분포 (Bernard 2007 감마형)

```
psd2 = 1e20·(psd/2)^((1-3·V_eff)/V_eff)·exp((-psd/2)/((Deff/2)·V_eff)),  V_eff=0.6
civol = π/6·Σ(psd2·psdm1³·Δd),  psdm2 = psd2/(civol·ci)   [Chl 1 mg/m³ 정규화]
```

### 4.4 Chl → 절대 IOP 변환

크기분포를 Chl 1 mg/m³ 기준으로 정규화하므로 산출 a*, b*는 Chl-specific(m²/mg)이다.
따라서 절대 변환은 단순 곱이다.
```
a_phyto(λ) = a*(λ) · Chl,   b_phyto(λ) = b*(λ) · Chl
```
Chl=1이면 a*, b* 값이 그대로 절대값이다. 별도 변환 수식 불필요.

### 4.5 다중 그룹 옵션 (Brewin 크기계급)

Brewin 2010 3-성분 모델로 Chl → pico/nano/micro 분율. 각 계급 대표 그룹:
- pico:  Synechococcus (Deff~1µm)  — a*(443)=0.093, bb/b(550)=0.0144, DoLP90=0.960
- nano:  Haptophytes (Deff~5µm)    — a*(443)=0.061, bb/b(550)=0.0043, DoLP90=0.849
- micro: Diatoms centric (Deff~20µm) — a*(443)=0.034, bb/b(550)=0.0035, DoLP90=0.844

소형 세포일수록 a* 크고, bb/b 높고, DoLP 높다(package effect 약함). 물리적으로 타당.

### 4.6 350nm 확장 (Hilbert edge 통제)

OCRT 시작 밴드는 350nm. EAP kshell 데이터는 380-900nm이므로 350-380 외삽 필요.
단순 PCHIP 외삽은 발산한다(nshell<1, 겹침구간 왜곡 26.9%). Hilbert가 전역 변환이라
끝단 오염이 안쪽으로 번지기 때문이다.

해결: 넓은 범위(280-1020nm) 지수감쇠 외삽 + 감쇠율 τ 최적화 + 신뢰구간 추출.
```
kshell 외삽: k(경계)·exp(-Δλ/τ),  τ = scipy minimize_scalar (bounds 15-400)
목적함수: 겹침구간(400-850) nshell 매끄러움 + 물리제약(nshell>1) 페널티
Hilbert를 280-1020 넓은 범위에서 계산 → 350-850 신뢰구간만 추출
```
결과: 겹침구간 왜곡 26.9% → 0.45%, nshell(350)=1.089(>1), 발산 없음.

---

## 5. 유기 detritus

### 5.1 산란 — Stramski 2001

Stramski, Bricaud & Morel (2001), Appl. Opt. 40:2929-2945.
파라미터 출처: Stramski et al. (2004), Prog. Oceanogr. 61:27-56, p.42.
```
크기분포: Junge slope 4, D 0.05-500µm, 로그격자
실수 굴절률(상대): 1.04, 전 파장 일정
허수 굴절률: n_imag(λ) = 0.010658·exp(-0.007186·λ)·0.75  [λ nm, 25%↓ 오타정정]
```
검증: σbb ∝ λ^-0.9 재현(γ=0.93). 무기물(n=1.18)은 σbb ∝ λ^-0.91, 후방산란 30배 강.

위상함수 IQU(P11/P12/P33)를 wiscombe 정규화로 크기분포 적분(350-850nm, 361각).
검증(550nm): g=0.967, bb/b=0.0058, DoLP(90°)=0.979, 4π 정규화=1.0000.
detritus DoLP(0.98)는 식물플랑크톤(0.15-0.84)보다 훨씬 높아 편광 구분 신호가 된다.

### 5.2 흡광 — Bricaud & Stramski 1990

Bricaud & Stramski (1990), Limnol. Oceanogr. 35(3):562-582.
DOI 10.4319/lo.1990.35.3.0562. 비조류 biogenous matter를 phytoplankton과 분리 측정.
"biogenous"이므로 무기물 배제, "particulate"이므로 CDOM 배제 → 우리 원칙에 부합.

지수식 부여(굴절률 Mie 대신). 이유: S_d가 크기분포 의존 가변량이므로 사용자 옵션으로
노출하려면 지수식이 자연스럽다.
```
a_d(λ) = a_d(440) · exp(-S_d·(λ-440))
기본 S_d = 0.0109 nm⁻¹  (Sargasso Sea 평균, SD 0.0019)
옵션 범위 = 0.0024-0.017 nm⁻¹  (Peru+Sargasso 전체 실측; 크기분포 의존)
```
Sargasso 채택 이유: 우리 detritus 굴절률이 Iturriaga & Siegel 1989(Sargasso detritus)
유래라 같은 해역 계열로 일관.
물리 검증: 우리 굴절률 Mie 계산값 S_d=0.0088이 실측 범위 정중앙 → 지수식 기본값 0.0109
타당함을 뒷받침.

### 5.3 detritus 흡광의 Chl 의존성 — B안

Bricaud & Stramski 1990 그림 18: a_d(440)/a_p(440)이 Chl 증가에 따라 감소
(Sargasso 12-82%, Peru 5-42%, 산포 큼). detritus 흡광은 Chl covary하되 지수가
phyto 흡광(~0.68)보다 작다(저Chl에서 상대적으로 큼, 고Chl에서 phyto에 밀림).

산점도 산포가 커 단일 정량 지수 추출은 과대정밀. 정성 방향만 채택하고, 절대 스케일은
지역 보정 대상으로 둔다(covary 구조는 산란 배경 성분과 동일).

주의: detritus는 산란과 흡광의 Chl 의존성이 다르다. 산란은 배경 지배(Chl^0.565 계열,
약한 의존), 흡광은 covary(Chl 지수 phyto보다 작으나 산란보다 강함). 같은 성분이라도
a와 b의 Chl 의존성이 다르므로 독립 스케일해야 한다(6.2 결정 2).

---

## 6. OCRT 병합 — 절대 스케일 연결 (방향 B → OCRT)

### 6.1 OCRT 병합 구조

OCRT는 각 성분의 a, b, 위상함수 P를 받아 병합한다.
```
a_total = Σ a_i                      (흡광 단순 합)
b_total = Σ b_i                      (산란 단순 합)
P_total = Σ (b_i · P_i) / b_total    (위상함수 산란 가중 평균)
```
방향 B는 bb 기준으로 절대 스케일을 정했으나, OCRT 병합은 b 가중이다. bb ≠ b이므로
bb에서 b로 역산해야 한다.

### 6.2 세 결정

**결정 1 (a안): phyto의 모든 절대 크기를 f_ph 비중에 종속.**
phyto bb를 f_ph로 억누른 뒤 EAP bb/b로 나눠 b를 구한다. b도 함께 억눌린다.
EAP는 형태(분광·각도·편광)만 제공하고, 절대 크기(a, b, bb)는 실측 비중에 종속.
```
bb_phyto = f_ph(Chl) · bbp_Huot(550)   [절대, 방향 B]
b_phyto  = bb_phyto / (bb/b)_EAP        [bb → b 역산]
```

**결정 2 (분리): 흡광과 산란의 Chl 의존성을 독립 스케일.**
a와 b는 다른 Chl 지수를 따른다. 분리하여 각각 독립 스케일한다.
```
산란·후방산란: 방향 B (Huot bbp + f_ph 분배)
흡광 phyto:    EAP a*·Chl (Chl-specific, 4.4)
흡광 detritus: a_d(440) Chl covary × 지수식 형태 (5.2, 5.3)
```

**b 모델은 독립 모델이 아니다.** b는 Huot bbp에서 성분별 bb/b로 역산된다.
```
b_i(λ,Chl) = bb_i(λ,Chl) / (bb/b)_i(λ)
  bb_i: 방향 B (Huot 총 + f_ph 분배)
  (bb/b)_i: 각 성분 .mie의 계산값 (EAP, Stramski Mie)
```

**결정 3: detritus 흡광 절대 스케일 a_d(440).**
Chl covary로 결정(지역 보정). 분광 형태만 Bricaud & Stramski 1990 지수식.

### 6.3 성분별 IOP 산출 요약

| 성분 | a(흡광) | b(산란) | P(위상함수) |
|---|---|---|---|
| 식물플랑크톤 | EAP a*·Chl | bb_phyto / (bb/b)_EAP | EAP .mie P11/P12/P33 |
| 유기 detritus | a_d(440)·exp(-S_d·(λ-440)), Chl covary | bb_det / (bb/b)_Stramski | Stramski .mie P11/P12/P33 |
| 무기물 TSM | Ahn a*·conc | Ahn b*·conc | Ahn .mie P11/P12/P33 |

최종: a_total, b_total 단순 합, P_total은 b 가중 병합.

---

## 7. 산출물 (.mie, 모두 350-850nm, 361각)

| 파일 | 성분 | 비고 |
|---|---|---|
| pico_Synechococcus_EAP.mie | 식물플랑크톤 pico | Brewin 소형 |
| nano_Haptophytes_EAP.mie | 식물플랑크톤 nano | Brewin 중형 |
| Diatoms_centric_EAP.mie | 식물플랑크톤 micro | Brewin 대형, 연안 우점 |
| Detritus_Stramski2001.mie | 유기 detritus | Chl covary |
| Red_clay_AHN.mie 등 4종 | 무기물 TSM | 독립 농도 |

모든 성분이 350nm에서 시작해 파장 정합. OCRT 시작 밴드(350nm)와 일치.

---

## 8. 레퍼런스

- Lain, Kravitz, Matthews, Bernard (2023). EAP. Sci Data 10:412. DOI 10.1038/s41597-023-02310-z
- Bernard et al. (2007). EAP 크기분포 감마형.
- Stramski, Bricaud & Morel (2001). detritus 산란. Appl. Opt. 40:2929-2945.
- Stramski et al. (2004). detritus 파라미터·phyto 비중. Prog. Oceanogr. 61:27-56.
- Iturriaga & Siegel (1989). detritus 굴절률(Sargasso). Limnol. Oceanogr. 34:1706-1726.
- Bricaud & Stramski (1990). detritus 흡광. Limnol. Oceanogr. 35(3):562-582. DOI 10.4319/lo.1990.35.3.0562
- Huot et al. (2008). 총 bbp Chl 함수. Biogeosciences 5:495-507.
- Brewin et al. (2010). 3-성분 크기계급.
- Toon & Ackerman (1981). d'Milay 이중층 구 코드.
