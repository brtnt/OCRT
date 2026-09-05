# OCRT 330–1100 nm 광특성 분광자료 및 FR631 Mie 재생성 기술문서

**문서 버전:** 2026-08-20  
**범위:** elastic coupled ocean-atmosphere OCRT, 330–1100 nm

## 1. 설계 목표

이 작업의 목표는 OCRT의 모든 production 광특성 입력을 330–1100 nm로 통일하고, 입자 위상함수의 전방 산란 표현을 고해상도 고정격자로 전환하는 것이다. 최종 data contract는 “코드에서 범위 밖 값을 임의 보정”하는 방식이 아니라 “입력자료가 지원범위를 명시적으로 포함”하는 방식이다.

## 2. 최종 자료 인벤토리

| 분류 | 수량 | 최종 범위 |
|---|---:|---|
| 대기 SNF/OPAC Mie | 16 | 330–1100 nm, 1 nm, FR631 |
| Ahmad2010 paper Mie | 80 | 330–1100 nm, 1 nm, FR631 |
| Ahmad2010 AccuRT Mie | 80 | 330–1100 nm, 1 nm, FR631 |
| AHN 광물 Mie | 4 | 330–1100 nm, 1 nm, FR631 |
| detritus Mie | 1 | 330–1100 nm, 1 nm, FR631 |
| EAP phytoplankton Mie | 17 | 330–1100 nm, 1 nm, FR631 |
| PLOPS Chl absorption | 1 | 330–1100 nm, 1 nm |
| AHN scalar a*, b* | 8 | 330–1100 nm 이상 |
| absorbing-gas LUT | 6 | 330–1100 nm, 1 nm, 40 levels |
| pure-water ψT | 1 | 330–1100 nm, 1 nm |

## 3. FR631 각도 정본

FR631은 631점의 내림차순 angle grid이다.

| 구간 | 간격 | interval 수 |
|---|---:|---:|
| 180→20° | 0.5° | 320 |
| 20→5° | 0.1° | 150 |
| 5→1° | 0.05° | 80 |
| 1→0.2° | 0.02° | 40 |
| 0.2→0° | 0.005° | 40 |

경계각은 한 번만 포함한다. 첨부 규격문서의 explicit 631-line list와 segment definitions를 정본으로 사용했다. 문서 본문의 요약 문장에 적힌 `[0°,1°)` 및 `[0°,10°)` 점수는 explicit list와 일치하지 않으므로, 납품파일은 explicit list를 문자 단위로 따랐다. 실제 점수는 각각 80, 210이다.

각 source wavelength에서 single-particle 또는 ensemble Mie phase를 FR631의 모든 631각에서 직접 계산했다. 0.5° 또는 0.1° 자료를 계산한 뒤 0.005° 전방구간을 보간으로 채우지 않았다.

## 4. 공통 Mie phase 물리·수치 계약

### 4.1 정규화

각 wavelength column은

```text
(1/2) ∫[-1,1] P11(μ) dμ = 1
```

로 정규화한다. 전달된 bulk asymmetry parameter는 동일한 FR631 P11에서

```text
g = (1/2) ∫[-1,1] μ P11(μ) dμ
```

로 다시 계산한다.

### 4.2 1 nm phase materialization

Mie solver는 계열별 승인된 source wavelength nodes에서 직접 계산하고, 그 사이 정수 파장 column은 다음 방법으로 materialize했다.

- P11: wavelength 방향 log-PCHIP
- P12/P11: wavelength 방향 PCHIP 후 [-1,1] 제한
- P33/P11: wavelength 방향 PCHIP 후 [-1,1] 제한
- P12, P33 재구성 후 각 1 nm column 재정규화
- g를 최종 FR631 P11 column에서 재계산

이는 **wavelength interpolation**이며, 금지된 것은 **angle interpolation에 의한 FR631 생성**이다. 런타임에서는 이미 1 nm column이 존재하므로 fractional wavelength만 선형 보간한다.

### 4.3 Mueller 물리성

모든 column에서 다음을 검사했다.

```text
P11 > 0
|P12| <= P11
|P33| <= P11
```

## 5. 대기 에어로졸 176종

### 5.1 계열

- SnF/OPAC 16종
- Ahmad et al. (2010) paper-centric 80종
- Ahmad2010 AccuRT-source 80종

Ahmad 80종은 8 RH × 10 fine-volume fraction 조합이다.

### 5.2 굴절률 파장 확장

기존 input card 20-node를 출발점으로 330 nm와 1100 nm를 추가했다.

- 330 nm: 350–400 nm 선형외삽
- 1100 nm: 860–1240 nm 선형내삽
- real/imaginary part 독립 처리

직접 phase source nodes:

```text
330, 350, 400, 412, 443, 470, 488, 515,
550, 590, 633, 670, 694, 760, 860, 1100 nm
```

bulk extinction은 positive/log interpolation, SSA는 direct interpolation 후 [0,1] 제한, scattering은 extinction×SSA로 재구성했다.

## 6. AHN 광물 4종

### 6.1 원자료

Ahn (1990) 박사학위논문의 Chapter 3을 직접 근거로 했다.

- Fig. 3-2: 상대 입경분포
- Table 3.1: m, dmin, dmax
- Fig. 3-7a: imaginary relative refractive index
- Fig. 3-7b: real relative refractive index
- Fig. 3-3–3-6: Qa, Qb, Qbb, bb/b 검증자료

논문은 size parameter를 진공파장과 해수 굴절률 `nw=1.34`로 정의하고, 입자 굴절률을 해수에 대한 상대값으로 다룬다.

### 6.2 PSD 처리

Fig. 3-2의 디지타이징 곡선을 상대 number distribution으로 사용했다. 원 그림에서 x축과 구분되지 않는 마지막 point 이후의 tail은 central 계산에서 0으로 두었다. Junge continuation은 sensitivity용으로만 유지했다.

### 6.3 굴절률

- 400–750 nm: 재디지타이징 자료
- 330–399 nm: 400–420 nm 경향 선형외삽
- 751–1100 nm: 장파장 선형외삽
- real/imaginary 모두 50 nm low-pass
- calcareous sand imaginary part는 1100 nm에서 나머지 3종 평균으로 수렴

### 6.4 Mie 조건

- homogeneous equivalent sphere
- host medium: seawater
- host refractive index: 1.34
- vacuum wavelength
- size integration: log diameter grid
- direct FR631 calculation

이 자료는 Ahn source-informed reconstruction이며 원 역사적 front-end의 exact byte reproduction을 주장하지 않는다.

## 7. Detritus

`Detritus_Stramski2001.mie`는 Stramski, Bricaud and Morel (2001)의 detailed-composition IOP framework에 기반한 보존 recipe를 사용했다.

- homogeneous spheres
- dN/dD ∝ D⁻⁴
- diameter 0.05–500 μm
- real relative index 1.04
- wavelength-dependent imaginary index from preserved OCRT recipe
- seawater host n=1.334
- direct phase source nodes 22개: 330, 340, 350, 380, 400, 412, 443, 490, 510, 550, 620, 660, 680, 709, 745, 765, 780, 850, 865, 940, 1030, 1100 nm

## 8. EAP phytoplankton 17종

EAP는 현재 production default에서 비활성이다. 다만 advanced selection과 향후 활성화를 위해 17종 자료를 같은 규격으로 재생성했다.

- 기존 350–850 nm complex RI exact 보존
- 330–349 nm 및 851–1100 nm extension
- real n: edge 50 nm linear least-squares extrapolation
- imaginary k: edge 50 nm log-linear least-squares extrapolation
- direct FR631 phase at 16 source nodes
- 1 nm phase materialization

EAP를 기본 Chl scattering에 자동 연결하지 않는다.

## 9. PLOPS phytoplankton absorption

Lomas et al. (2024) PLOPS dataset은 50개 이상의 배양주에 대해 hyperspectral UV–VIS particle absorption, backscatter, VSF, PSD, pigments를 제공한다. OCRT 기본 `a*ph`는 diatom 계열 representative spectrum을 기반으로 준비했다.

- 330–750 nm: source-derived representative spectrum
- 750–800 nm: 0으로 선형 taper
- 800–1100 nm: exact zero
- 단위: Chl-a 1 mg m⁻³에 대한 specific absorption

## 10. AHN scalar a*, b*

네 광물의 350–400 nm 기존 값이 이미 외삽자료였기 때문에, 원 측정영역의 400–420 nm를 기준으로 330–399 nm를 선형회귀 외삽했다. 400 nm 이상 기존 canonical 값은 보존했다.

## 11. 흡광성 기체

6개 LUT는 330–1100 nm, 1 nm, AFGL US Standard 40 levels, cm² molecule⁻¹이다.

- H2O/O2/CO2/CH4: clean HITRAN transition downloads and HAPI-based line treatment
- O3: Serdyuchenko–Gorshelev 213–1100 nm, 193–293 K
- NO2: Vandaele 220/294 K UV spectral shape, canonical 350 nm layer value에 정합
- 350–1100 nm 기존 canonical 수치 보존
- O2–O2/O4 CIA는 이번 scope에서 제외

## 12. 순수해수 흡광 온도계수

```text
aw(λ,T) = aw(λ,20°C) + ψT(λ)(T-20°C)
```

- 330–440 nm: ψT = 0 exact
- 441–1100 nm: Röttgers/WOPP source family
- 2 nm source nodes 보존, 홀수 파장 PCHIP
- 20°C baseline exact 보존

## 13. 최종 전수검증

| 항목 | 결과 |
|---|---:|
| canonical Mie | 198/198 PASS |
| total uncompressed size | 4,941,848,335 bytes |
| phase normalization 최대오차 | 1.0952×10⁻⁹ |
| phase-derived g vs bulk g 최대오차 | 5.8675×10⁻¹⁰ |
| P12 bound violation | 0 |
| P33 bound violation | 0 |
| SSA identity 최대오차 | 5.0233×10⁻⁹ |
| C full release build | PASS, warning 0 |
| C reader representative six families | PASS |
| Python reader representative six families | PASS |

## 14. 알려진 제한

- AHN은 exact historical reproduction이 아니라 source-informed forward reconstruction이다.
- EAP는 data-complete이나 default disabled이다.
- phase는 direct source wavelength nodes 사이에서 spectral materialization되었다.
- O4 CIA와 Raman은 scope 밖이다.
- 3개 legacy alias는 재생성하지 않았다.

## 15. 참고문헌

1. Ahn, Y.-H. (1990). *Propriétés optiques des particules biologiques et minérales présentes dans l'océan; application: inversion de la réflectance*. Thèse de doctorat, Université Paris VI.
2. Lomas, M. W., Neeley, A. R., Vandermeulen, R., Mannino, A., Thomas, C., Novak, M. G., & Freeman, S. A. (2024). Phytoplankton optical fingerprint libraries for development of phytoplankton ocean color satellite products. *Scientific Data*, 11, 168. doi:10.1038/s41597-024-03001-z.
3. Ahmad, Z., Franz, B. A., McClain, C. R., Kwiatkowska, E. J., Werdell, J., Shettle, E. P., & Holben, B. N. (2010). New aerosol models for the retrieval of aerosol optical thickness and normalized water-leaving radiances from the SeaWiFS and MODIS sensors over coastal regions and open oceans. *Applied Optics*, 49, 5545–5560. doi:10.1364/AO.49.005545.
4. Hess, M., Koepke, P., & Schult, I. (1998). Optical properties of aerosols and clouds: The software package OPAC. *Bulletin of the American Meteorological Society*, 79, 831–844. doi:10.1175/1520-0477(1998)079<0831:OPOAAC>2.0.CO;2.
5. Stramski, D., Bricaud, A., & Morel, A. (2001). Modeling the inherent optical properties of the ocean based on the detailed composition of the planktonic community. *Applied Optics*, 40, 2929–2945. doi:10.1364/AO.40.002929.
6. Röttgers, R., McKee, D., & Utschig, C. (2014). Temperature and salinity correction coefficients for light absorption by water in the visible to infrared spectral region. *Optics Express*, 22, 25093–25108. doi:10.1364/OE.22.025093.
7. Sullivan, J. M., Twardowski, M. S., Zaneveld, J. R. V., Moore, C. M., Barnard, A. H., Donaghay, P. L., & Rhoades, B. (2006). Hyperspectral temperature and salt dependencies of absorption by water and heavy water in the 400–750 nm spectral range. *Applied Optics*, 45, 5294–5309. doi:10.1364/AO.45.005294.
8. Kochanov, R. V., Gordon, I. E., Rothman, L. S., Wcisło, P., Hill, C., & Wilzewski, J. S. (2016). HITRAN Application Programming Interface (HAPI). *Journal of Quantitative Spectroscopy and Radiative Transfer*, 177, 15–30. doi:10.1016/j.jqsrt.2016.03.005.
9. Gorshelev, V., Serdyuchenko, A., Weber, M., Chehade, W., & Burrows, J. P. (2014). High spectral resolution ozone absorption cross-sections – Part 1. *Atmospheric Measurement Techniques*, 7, 609–624. doi:10.5194/amt-7-609-2014.
10. Serdyuchenko, A., Gorshelev, V., Weber, M., Chehade, W., & Burrows, J. P. (2014). High spectral resolution ozone absorption cross-sections – Part 2. *Atmospheric Measurement Techniques*, 7, 625–636. doi:10.5194/amt-7-625-2014.
11. Vandaele, A. C., et al. (1998). Measurements of the NO2 absorption cross-section from 42,000 to 10,000 cm⁻¹ (238–1000 nm) at 220 K and 294 K. *JQSRT*, 59, 171–184. doi:10.1016/S0022-4073(97)00168-4.
12. Bohren, C. F., & Huffman, D. R. (1983). *Absorption and Scattering of Light by Small Particles*. Wiley.
13. van de Hulst, H. C. (1957). *Light Scattering by Small Particles*. Wiley.
14. Fritsch, F. N., & Carlson, R. E. (1980). Monotone piecewise cubic interpolation. *SIAM Journal on Numerical Analysis*, 17, 238–246.

## 16. Raw FR631과 truncation 해석 (2026-08-21 확정)

FR631은 0°와 0.005°를 포함하는 명시적 고해상도 격자이다. `0–0.005°` local cap은 시험한 수중 Mie에서 raw 결과와 machine precision 수준으로 동일했다. 그러나 OSOAA-style broad truncation은 broad forward lobe에 영향을 주고 일부 조건에서 입자산란의 큰 비율을 transformed term으로 이동시키며 directional `Rrs/rrs`를 수 % 변경했다.

따라서 FR631 broad truncation OFF는 영향이 작기 때문이 아니라 raw phase가 원 복사전달방정식의 canonical representation이기 때문이다. Raw/raw와 matched broad/broad는 별도 validation mode이며 raw/broad 혼합 비교는 `REFERENCE_MISMATCH`이다.

Non-FR631을 일괄적으로 coarse라 부르지 않는다. Exact legacy 361×0.5°만 별도로 식별하여 FR631 교체를 우선 권고하고, 기타 custom grid는 forward-angle resolution과 RT convergence 확인 전까지 unvalidated이다.

## 17. Figure–CSV provenance

본 문서 및 후속 검증문서의 모든 정량 figure는 package 내 authoritative source CSV와 연결한다. 각 campaign은 `PLOT_DATA_MAP.csv`, reproduction script, campaign metadata 및 SHA-256 manifest를 포함한다.

입력자료 또는 자료처리가 갱신되면 source CSV, derived metrics 및 모든 관련 figure를 같은 revision으로 재생성하고 교체한다. CSV-only 또는 figure-only 갱신은 금지한다.

저장 위치와 감사 도구:

```text
05_VALIDATION/artifacts/<campaign_id>/
python 05_VALIDATION/tools/audit_validation_artifacts.py
```

상세 정책: `01_OCRT_C/docs/OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21.md`

