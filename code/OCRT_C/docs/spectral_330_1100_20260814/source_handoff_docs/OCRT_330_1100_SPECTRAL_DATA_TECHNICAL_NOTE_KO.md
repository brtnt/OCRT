# OCRT 330-1100 nm 광특성 분광자료 구축 기술문서

**대상:** OCRT C/Python paired implementation의 elastic atmosphere-ocean radiative-transfer production data

**작성 목적:** 교체자료의 물리량, 원자료, 처리식, 외삽·보간·스무딩 정책, Mie 계산계약, 검증결과, 한계를 재현 가능한 수준으로 기록한다.

## 1. 요약

이번 작업은 기존 OCRT의 분광자료를 330-1100 nm로 확장했다. 핵심 원칙은 코드 안의 파장별 예외분기 대신 자료 파일에 명시적 값을 기록하는 것이었다. 각 자료는 source-backed 구간과 OCRT project policy 구간을 구분하고, C와 Python에 byte-identical하게 배포한다.

최종 범위는 기본 phytoplankton absorption, AHN mineral scalar/vector optics, organic detritus phase candidate, atmospheric aerosol 176종, six-gas absorption, pure-water absorption temperature coefficient를 포함한다. EAP, Raman, CCRR legacy, O4/CIA는 범위 밖이다.

## 2. 공통 설계계약

### 2.1 명시적 자료계약

- runtime 목표범위 330-1100 nm를 자료가 직접 덮는다.
- source 범위 밖 처리는 데이터 생성단계에서 결정하고 헤더에 기록한다.
- constituent-specific `if (lambda...)`로 값을 대체하지 않는다.
- generic interpolation은 선언된 범위 내부에서만 사용한다.
- 범위 밖 요청은 fail-loud이며 endpoint hold는 과학적 외삽으로 인정하지 않는다.

### 2.2 C/Python paired maintenance

동일한 원파일을 C `inputs/`와 Python `data/`에 복사한다. 별도 직렬화나 반올림본을 만들지 않는다. 각 파일의 SHA-256을 설치맵에 고정한다.

### 2.3 Source와 project policy 구분

| 구분 | 예 |
|---|---|
| 직접 source-derived | PLOPS 330-750, Ahn 400-750 RI, Röttgers/WOPP psi_T, HITRAN line lists |
| 내부 interpolation | 2.5->1 nm PCHIP, source node 사이 linear/boxcar |
| project closure | PLOPS 750-800 taper, >=800 zero; psi_T 330-440 zero; selected RI extrapolations |
| 조건부 candidate | Detritus forward Mie; AHN source-informed phase caveat |


### 2.4 최종 교체자료 inventory

| 구성요소 | 파일 수 / 구현 | 기본 적용 | 핵심 범위·계약 |
|---|---:|---|---|
| PLOPS 기반 `a_ph*` | 1 | 예 | 330–750 source-derived; 750–800 zero taper; 800–1100 exact zero |
| AHN mineral scalar `a*`,`b*` | 8 | 예 | 330–399 OLS from canonical 400–420; 400 nm 이상 보존 |
| AHN mineral vector Mie | 4 | 예, caveat | 330–1100 bulk; 16 phase nodes; 361 angles; seawater host |
| Organic detritus Mie | 1 | 조건부 | 330–1100 bulk; 22 phase nodes; legacy residual 존재 |
| Atmospheric aerosol Mie | 176 | 예 | 기존 20 nodes + explicit 330/1100; 361 angles |
| Absorbing-gas xsec | 6 | 예 | 330–1100, 1 nm, 40 AFGL levels |
| Pure-water `psi_T` | 1 | 예 | 330–440 exact zero; 441–1100 WOPP/Röttgers family |

각 구현에는 총 197개 파일이 포함된다. 기본 installer는 조건부 detritus 1개를 제외한 196개를 설치하며, C와 Python payload는 파일 단위로 byte-identical하다. 개별 목적지와 SHA-256은 `05_VALIDATION_AND_MANIFESTS/INSTALLATION_MAP_C_PYTHON.csv`가 단일 기준이다.

## 3. 기본 phytoplankton absorption `a_ph*`

### 3.1 물리량

OCRT의 legacy 명칭 `a_chla*`는 순수 chlorophyll-a 분자 용액흡광이 아니라 chlorophyll-a-specific phytoplankton absorption coefficient이다.

```text
a_ph(lambda) [m-1]
  = a_ph_star(lambda) [m2 mg-1 Chl-a]
  * TChl-a [mg m-3]
```

따라서 `TChl-a=1 mg m-3`이면 두 값은 수치적으로 같다.

### 3.2 원자료

Lomas et al. (2024)의 PLOPS 자료는 50종 이상의 배양 phytoplankton에 대해 UV-VIS absorption, backscatter, size distribution, pigments 등을 제공한다. Filter-pad spectrophotometer는 290-850 nm를 측정했고 공개 workbook에 sample metadata와 `aph_spec`, HPLC `TChl a`가 포함된다.

`aph_spec`은 total particulate absorption에서 depigmented-particle absorption을 뺀 phytoplankton pigment absorption이다. Accessory pigments, package effect, 세포상태의 영향이 포함될 수 있다.

### 3.3 대표값 구성

1. `aph_spec`와 HPLC `TChl a`가 정합되는 90개 자료 확보
2. exponential phase, 비혼합, 일반 treatment, 농도 QC, nonnegative spectrum 등으로 선별
3. 강한 MAA-like UV outlier를 분리
4. 최종 34개, 6분류군 사용: Cyanobacteria, Diatom, Dinoflagellate, Haptophyte, Prasinophyte, Raphidophyte
5. 각 분류군 파장별 중앙값 계산
6. 6개 class median의 파장별 중앙값 계산
7. 1 nm 재격자화 후 5-point centered moving average, reflected padding

이는 single diatom 대표가 아니다. 합성 스펙트럼에 가장 가까운 실측자료가 `CCMP1316 Chaetoceros muelleri`였다.

### 3.4 장파장 정책

```text
330-750 nm  PLOPS-derived class-balanced a_ph*
750-800 nm  a_ph*(750)에서 800 nm exact zero로 선형 taper
800-1100 nm exact zero
```

750-800 및 zero 구간은 PLOPS 측정이 아니라 OCRT project closure이다.

![PLOPS default spectrum](../03_SOURCE_AND_PROVENANCE/phytoplankton_absorption/PLOPS_APH_DEFAULT_330_1100.png)

## 4. AHN mineral mass-specific absorption/scattering

### 4.1 원자료와 물리량

Ahn (1990)은 brown earth, yellow clay, red clay, calcareous sand의 해수현탁 광학특성, 상대 입경분포, 복소상대굴절률을 제시한다. 논문에서 `a*`, `b*`는 건조질량 농도에 정규화한 계수이다.

```text
a_min = a_star_min * TSM
b_min = b_star_min * TSM
```

단위는 `a*`, `b*`: m2 g-1; TSM: g m-3.

### 4.2 330-399 nm 확장

기존 350-399 값도 이미 외삽된 자료임을 확인했기 때문에 400-420 nm 21개 canonical 값에 ordinary least squares를 적용해 330-399를 모두 다시 계산했다.

```text
y(lambda) = slope * lambda + intercept
```

400 nm 이상 값은 그대로 보존했다. PCHIP 외삽은 endpoint derivative가 단파장으로 크게 증폭되고 calcareous sand에서 비물리적 음수를 만들 수 있어 폐기했다.

![AHN scalar b-star](../03_SOURCE_AND_PROVENANCE/ahn_tsm_scalar/AHN_TSM_B_STAR_OFFICIAL_300_1100.png)

## 5. AHN mineral vector Mie

### 5.1 디지타이징 자료

Ahn (1990) Chapter 3에서 다음을 재디지타이징하고 원 그림 overlay로 검수했다.

- Fig. 3-2: relative PSD `F_r(D)`
- Fig. 3-7a: relative imaginary RI `n'(lambda)`
- Fig. 3-7b: relative real RI `n(lambda)`
- Fig. 3-3~3-6: Qa, Qb, Qbb, bb/b 검증곡선
- Table 3.1: `m`, `Dmin`, `Dmax`

Ahn thesis는 size parameter를 `alpha=pi D n_w/lambda_vac`, 해수굴절률 `n_w=1.34`, 입자굴절률을 외부매질에 대한 relative complex RI로 정의한다.

### 5.2 굴절률 준비

```text
400-750 nm  high-precision digitization
330-399 nm  linear extrapolation from the approved shortwave trend
751-1100 nm linear extrapolation from the longwave trend
```

Calcareous sand의 허수부가 장파장에서 비정상 증가하므로 750 nm를 보존하고 1100 nm에서 다른 3종의 평균으로 수렴하도록 선형 연결했다. 디지타이징 noise를 줄이기 위해 실수부·허수부 모두 약 50 nm low-pass를 적용하고, 2.5 nm 표를 PCHIP로 1 nm 내부보간했다.

![AHN real RI](../03_SOURCE_AND_PROVENANCE/ahn_mineral_mie/full_extracted/plots/01_RI_REAL_INPUT_330_1100.png)

![AHN imaginary RI](../03_SOURCE_AND_PROVENANCE/ahn_mineral_mie/full_extracted/plots/02_RI_IMAG_INPUT_330_1100.png)

### 5.3 PSD와 tail

Fig. 3-2를 relative number distribution `dN/dD`로 해석했다. Figure의 마지막 판독점은 약 1-1.5 pixel 검출한계이며 값은 peak의 약 0.4-0.6%였다. 중앙계산은 이후 tail을 zero로 놓았다. 보수적 sensitivity에서는 마지막 값에서 Table 3.1 Junge exponent로 `Dmax`까지 감소시켰다.

### 5.4 Mie 계산계약

```text
particle        homogeneous equivalent sphere
host medium     seawater
n_water         1.34
particle RI     seawater-relative n - i n'
wavelength      vacuum
size parameter  x = pi D n_water / lambda_vac
PSD integration logarithmic, Delta log10 D ~ 0.0025
internal angle  0.1 degree
output angle    361, 180->0 degree, 0.5 degree
matrix          P11, P12, P33
normalization   0.5 integral P11 dmu = 1
```

Bulk 330-1100 nm 1 nm, phase node 16개: 330, 350, 400, 412, 443, 470, 488, 515, 550, 590, 633, 670, 694, 760, 860, 1100 nm.

![AHN asymmetry](../03_SOURCE_AND_PROVENANCE/ahn_mineral_mie/full_extracted/plots/04_MIE_G_330_1100.png)

### 5.5 판정과 한계

구조, normalization, Mueller bounds, C/Python loader, CLI smoke, 수치수렴은 통과했다. 그러나 Ahn 도식의 Qa/Qb/Qbb/bb/b와 일부 16-60% 차이가 남는다. 원 Coulter bin table과 원 front-end가 없으므로 exact historical reproduction을 주장하지 않는다. RI/PSD를 legacy 출력에 맞추는 inverse fitting은 수행하지 않았다.

## 6. Organic detritus phase candidate

### 6.1 레시피

```text
shape           homogeneous sphere
dN/dD           proportional to D^-4
D               0.05-500 um
n_real          1.04 relative
k(lambda)       0.010658 exp(-0.007186 lambda_nm) * 0.75
n_water         1.334
size parameter  2 pi n_water r / lambda_vac
```

Bulk 330-1100 nm 1 nm, vector phase 22 nodes, 361 angles, P11/P12/P33로 생성했다.

### 6.2 조건부 이유

Legacy 550 nm와 비교해 `g` 약 +0.00194, `bb/b` 약 -6.6%, DoLP90 약 -0.0278의 잔차가 남는다. 이는 source recipe를 재계산한 candidate이며 exact legacy reproduction이 아니다. 따라서 패키지에는 포함하지만 기본 자동설치에서 제외한다.

![Detritus phase comparison](../03_SOURCE_AND_PROVENANCE/detritus_mie/full_extracted/validation/04_P11_550_generated_vs_legacy.png)

## 7. Atmospheric aerosol 176-model Mie

### 7.1 모델군

- SnF/OPAC family 16
- Ahmad et al. (2010) paper-centric 80
- Ahmad AccuRT-source family 80

모든 `.inp`가 존재했고 1:1 대응을 확인했다.

### 7.2 RI node 확장

기존 20-node RI를 보존하고 다음 2개를 추가했다.

```text
330 nm  350-400 nm 두 점 선형외삽, n/k 독립
1100 nm 860-1240 nm 선형내삽, n/k 독립
```

Generator는 20->22 node로 확장하고 550 nm normalization index를 파장값 검색으로 변경했다. 기존 kernel, PSD, mixing, angle grid는 유지했다.

### 7.3 검증

176/176 생성·reader·physicality를 통과했다. SnF/OPAC 16과 AccuRT 80은 기존 20 node를 exact 재현했다. Paper-centric 80은 보존 generator/front-end 차이로 최대 spectral relative L2 약 1.33%, P11 NRMSE 약 1.59%였고 사용자가 수용했다.

![Aerosol g](../03_SOURCE_AND_PROVENANCE/aerosol_mie/full_extracted/plots/aerosol_ahmad2010_paper_mie_g_330_1240.png)

## 8. Absorbing-gas cross sections

### 8.1 데이터 계약

```text
six gases       H2O, O2, CO2, CH4, O3, NO2
wavelength      330-1100 nm, 1 nm, 771 points
vertical grid   AFGL US Standard 1976, 40 levels
unit            cm2 molecule-1
```

350-1100 canonical 값은 전 기체에서 보존하고 330-349만 새로 추가했다.

### 8.2 H2O/O2/CO2/CH4

H2O는 clean global-isotopologue HITRAN line list와 HAPI 1.3 Voigt 계산으로 생성했다. O2, CO2, CH4는 UV 구간과 safety margin에 전이선이 없음을 확인해 explicit zero를 기록했다.

### 8.3 O3

Serdyuchenko-Gorshelev 213-1100 nm, 193-293 K 11온도 자료를 사용했다. 유한 음수 측정잡음을 0으로 제한한 후 1 nm boxcar와 층온도 선형보간을 적용했다. 350 nm 이상은 canonical을 그대로 유지했다.

### 8.4 NO2

HITRAN XSC의 Vandaele 220/294 K 자료는 238-666.6 nm를 덮는다. 그러나 기존 canonical은 Bogumil/SCIAMACHY 5온도 family이다. 따라서 330-349 nm에는 Vandaele의 온도별 상대분광형을 사용하고 각 AFGL 층의 canonical 350 nm 값에 scale을 맞췄다. source-identical Bogumil UV 확장은 아니며 향후 원 5온도 파일 확보 시 대체할 수 있다.

### 8.5 미포함

O2-O2/O4 CIA는 현재 six-gas 계약에 포함하지 않았다.

![Gas cross sections](../03_SOURCE_AND_PROVENANCE/gas_absorption/full_extracted/plots/01_GROUND_CROSS_SECTIONS_330_1100.png)

![Gas vertical optical depth](../03_SOURCE_AND_PROVENANCE/gas_absorption/full_extracted/plots/05_VERTICAL_OPTICAL_DEPTH_330_1100.png)

## 9. Pure-water absorption temperature coefficient `psi_T`

### 9.1 식

```text
a_w(lambda,T) = a_w(lambda,20 C) + psi_T(lambda)*(T-20 C)
```

`psi_T` 단위는 m-1 degC-1이다.

### 9.2 source and policy

Röttgers, McKee & Utschig (2014)는 400~2700 nm에서 temperature/salinity coefficients를 측정했고 400-700 nm에는 PSICAM 정밀자료를 제공한다. WOPP `purewater_abs_coefficients_v3.dat`의 2 nm node를 사용했다.

논문은 300-440 nm에서 temperature coefficient가 현재 측정정밀도상 zero와 유의하게 다르지 않다고 보고한다. 따라서 330-440 nm를 exact zero로 자료에 기록했다. 441-1100은 440 zero anchor와 WOPP 2 nm nodes 사이 PCHIP 내부보간이며, 442-1100 even source nodes는 exact 보존했다. Negative `psi_T`는 물리적인 band-shift response이므로 유지했다.

![psi_T final](../03_SOURCE_AND_PROVENANCE/pure_water_temperature/full_extracted/plots/01_PSI_T_FINAL_330_1100.png)

![a_w temperature spectra](../03_SOURCE_AND_PROVENANCE/pure_water_temperature/full_extracted/plots/04_AW_TEMPERATURE_SPECTRA_330_1100.png)

## 10. 검증체계

### 10.1 Data-level

- range, monotonic wavelength, row/node count
- SHA-256 and deterministic replay
- exact source-node preservation where claimed
- boundary continuity

### 10.2 Mie physicality

- `P11>0`
- `|P12|<=P11`
- `|P33|<=P11`
- half-integral normalization
- finite `g`, `bb/b`, SSA
- size/angular quadrature convergence

### 10.3 Cross-language

C/Python은 같은 data file을 읽고 matched wavelength/component에서 IOP 및 phase outputs를 비교했다. 최종 본 OCRT 통합 후에는 RT I/Q/U까지 같은 입력으로 재검증해야 한다.

## 11. 불확실성과 향후 교체 우선순위

1. Detritus: legacy residual 때문에 explicit adoption gate 필요
2. AHN PSD tail: zero-tail과 Junge sensitivity 차이가 장파장 `g`, `bb/b`에 영향
3. AHN exact history: 원 Coulter channels/front-end 부재
4. NO2: Bogumil 5온도 source-identical 330-349 자료 확보 시 교체 가능
5. O4: UV 완결성을 높이려면 별도 CIA 수직적분 구현 필요
6. EAP/Raman: 현재 범위 밖

## 12. 참고문헌과 구성요소별 근거맵

### 12.1 구성요소별 직접 근거

| 구성요소 | 직접 원자료·주 레퍼런스 | 보조 물리·수치 레퍼런스 |
|---|---|---|
| PLOPS `a*ph` | Lomas et al. (2024); Dryad dataset; Neeley et al. (2022) | Stramski et al. (2015) filter-pad correction; Kishino et al. (1985) depigmentation |
| AHN mineral scalar/vector | Ahn (1990) thesis | Morel & Ahn (1991); van de Hulst (1957); Bohren & Huffman (1983); Fritsch & Carlson (1980) |
| Detritus vector Mie | OCRT documented Stramski-type recipe | Stramski et al. (2001, 2004); Bricaud & Stramski (1990); Huot et al. (2008) |
| Aerosol 176 models | Ahmad et al. (2010); OPAC/Hess et al. (1998) | Vermote et al. (1997) 6S/EXSCPHASE lineage |
| Gas cross sections | HITRAN2020; HAPI; Serdyuchenko/Gorshelev O3; Vandaele NO2; AFGL profiles | Bogumil et al. (2003) canonical NO2 lineage |
| Pure-water `psi_T` | Röttgers et al. (2014); WOPP 2016 table | Sullivan et al. (2006); Pegau et al. (1997) |
| Unchanged pure-water scattering | Zhang et al. (2009) | 기존 OCRT canonical 유지 |

### 12.2 전체 서지

아래 목록은 본 인계자료의 직접 원자료, 물리모델, 독립검증 및 수치방법을 구분해 수록한다. Zero 구간, taper endpoint, low-pass 폭, accepted overlap tolerance 등은 문헌 측정값이 아니라 OCRT project policy이며 각 data header와 기술문서에 별도로 표시한다.
1. Ahn, Y.-H. (1990). *Propriétés optiques des particules biologiques et minérales présentes dans l'océan; application: inversion de la réflectance*. Doctoral thesis, Université Paris VI. Direct source for the four mineral PSDs, complex relative refractive indices, and bulk optical quantities.
2. Morel, A., & Ahn, Y.-H. (1991). Optics of heterotrophic nanoflagellates and ciliates: A tentative assessment of their scattering role in oceanic waters compared to those of bacterial and algal cells. *Journal of Marine Research, 49*(1), 177-202. DOI: 10.1357/002224091784968639. Contextual particle-optics reference; the thesis remains the direct mineral digitization source.
3. van de Hulst, H. C. (1957). *Light Scattering by Small Particles*. Wiley.
4. Bohren, C. F., & Huffman, D. R. (1983). *Absorption and Scattering of Light by Small Particles*. Wiley.
5. Lomas, M. W., Neeley, A. R., Vandermeulen, R., Mannino, A., Thomas, C., Novak, M. G., & Freeman, S. A. (2024). Phytoplankton optical fingerprint libraries for development of phytoplankton ocean color satellite products. *Scientific Data, 11*, 168. DOI: 10.1038/s41597-024-03001-z.
6. Lomas, M. W., et al. (2023). *Phytoplankton optical fingerprint libraries for development of phytoplankton ocean color satellite products* [Dataset]. Dryad. DOI: 10.5061/dryad.rbnzs7hfg.
7. Neeley, A. R., Lomas, M. W., Mannino, A., Thomas, C. S., & Vandermeulen, R. A. (2022). Impact of growth phase, pigment adaptation, and climate change conditions on the cellular pigment and carbon content of fifty-one phytoplankton isolates. *Journal of Phycology, 58*, 669-690. DOI: 10.1111/jpy.13279.
8. Stramski, D., Reynolds, R. A., Kaczmarek, S., Uitz, J., & Zheng, G. (2015). Correction of pathlength amplification in the filter-pad technique for measurements of particulate absorption coefficient in the visible spectral region. *Applied Optics, 54*(22), 6763-6782. DOI: 10.1364/AO.54.006763.
9. Kishino, M., Takahashi, M., Okami, N., & Ichimura, S. (1985). Estimation of the spectral absorption coefficients of phytoplankton in the sea. *Bulletin of Marine Science, 37*(2), 634-642.
10. Stramski, D., Bricaud, A., & Morel, A. (2001). Modeling the inherent optical properties of the ocean based on the detailed composition of the planktonic community. *Applied Optics, 40*, 2929-2945. DOI: 10.1364/AO.40.002929.
11. Stramski, D., Boss, E., Bogucki, D., & Voss, K. J. (2004). The role of seawater constituents in light backscattering in the ocean. *Progress in Oceanography, 61*, 27-56. DOI: 10.1016/j.pocean.2004.07.001.
12. Bricaud, A., & Stramski, D. (1990). Spectral absorption coefficients of living phytoplankton and nonalgal biogenous matter: a comparison between the Peru upwelling area and the Sargasso Sea. *Limnology and Oceanography, 35*, 562-582. DOI: 10.4319/lo.1990.35.3.0562.
13. Huot, Y., Morel, A., Twardowski, M. S., Stramski, D., & Reynolds, R. A. (2008). Particle optical backscattering along a chlorophyll gradient in the upper layer of the eastern South Pacific Ocean. *Biogeosciences, 5*, 495-507. DOI: 10.5194/bg-5-495-2008.
14. Fritsch, F. N., & Carlson, R. E. (1980). Monotone piecewise cubic interpolation. *SIAM Journal on Numerical Analysis, 17*(2), 238-246. DOI: 10.1137/0717021.
15. Ahmad, Z., Franz, B. A., McClain, C. R., Kwiatkowska, E. J., Werdell, J., Shettle, E. P., & Holben, B. N. (2010). New aerosol models for the retrieval of aerosol optical thickness and normalized water-leaving radiances from the SeaWiFS and MODIS sensors over coastal regions and open oceans. *Applied Optics, 49*, 5545-5560. DOI: 10.1364/AO.49.005545.
16. Hess, M., Koepke, P., & Schult, I. (1998). Optical properties of aerosols and clouds: The software package OPAC. *Bulletin of the American Meteorological Society, 79*, 831-844. DOI: 10.1175/1520-0477(1998)079<0831:OPOAAC>2.0.CO;2.
17. Vermote, E. F., Tanré, D., Deuzé, J.-L., Herman, M., & Morcrette, J.-J. (1997). Second Simulation of the Satellite Signal in the Solar Spectrum, 6S: an overview. *IEEE Transactions on Geoscience and Remote Sensing, 35*, 675-686. DOI: 10.1109/36.581987.
18. Gordon, I. E., et al. (2022). The HITRAN2020 molecular spectroscopic database. *Journal of Quantitative Spectroscopy and Radiative Transfer, 277*, 107949. DOI: 10.1016/j.jqsrt.2021.107949.
19. Kochanov, R. V., Gordon, I. E., Rothman, L. S., Wcisło, P., Hill, C., & Wilzewski, J. S. (2016). HITRAN Application Programming Interface (HAPI): A comprehensive approach to working with spectroscopic data. *Journal of Quantitative Spectroscopy and Radiative Transfer, 177*, 15-30. DOI: 10.1016/j.jqsrt.2016.03.005.
20. Gorshelev, V., Serdyuchenko, A., Weber, M., Chehade, W., & Burrows, J. P. (2014). High spectral resolution ozone absorption cross-sections - Part 1: Measurements, data analysis and comparison with previous measurements around 293 K. *Atmospheric Measurement Techniques, 7*, 609-624. DOI: 10.5194/amt-7-609-2014.
21. Serdyuchenko, A., Gorshelev, V., Weber, M., Chehade, W., & Burrows, J. P. (2014). High spectral resolution ozone absorption cross-sections - Part 2: Temperature dependence. *Atmospheric Measurement Techniques, 7*, 625-636. DOI: 10.5194/amt-7-625-2014.
22. Serdyuchenko, A., & Gorshelev, V. *Ozone absorption cross-section dataset, 213-1100 nm, 193-293 K*. Zenodo. DOI: 10.5281/zenodo.5793207.
23. Vandaele, A. C., Hermans, C., Simon, P. C., Carleer, M., Colin, R., Fally, S., Mérienne, M.-F., Jenouvrier, A., & Coquart, B. (1998). Measurements of the NO2 absorption cross-section from 42000 cm-1 to 10000 cm-1 (238-1000 nm) at 220 K and 294 K. *Journal of Quantitative Spectroscopy and Radiative Transfer, 59*, 171-184. DOI: 10.1016/S0022-4073(97)00168-4.
24. Bogumil, K., et al. (2003). Measurements of molecular absorption spectra with the SCIAMACHY pre-flight model: instrument characterization and reference data for atmospheric remote-sensing in the 230-2380 nm region. *Journal of Photochemistry and Photobiology A: Chemistry, 157*, 167-184. DOI: 10.1016/S1010-6030(03)00062-5.
25. Anderson, G. P., Clough, S. A., Kneizys, F. X., Chetwynd, J. H., & Shettle, E. P. (1986). *AFGL Atmospheric Constituent Profiles (0-120 km)*. AFGL-TR-86-0110.
26. Röttgers, R., McKee, D., & Utschig, C. (2014). Temperature and salinity correction coefficients for light absorption by water in the visible to infrared spectral region. *Optics Express, 22*, 25093-25108. DOI: 10.1364/OE.22.025093.
27. Röttgers, R. (2016). *Combined data for Water Radiance pure-water absorption, salinity and temperature coefficients and related uncertainties*, `purewater_abs_coefficients_v3.dat`, HZG/WOPP data record, Nov. 2016.
28. Sullivan, J. M., Twardowski, M. S., Zaneveld, J. R. V., Moore, C. M., Barnard, A. H., Donaghay, P. L., & Rhoades, B. (2006). Hyperspectral temperature and salt dependencies of absorption by water and heavy water in the 400-750 nm spectral range. *Applied Optics, 45*, 5294-5309. DOI: 10.1364/AO.45.005294.
29. Pegau, W. S., Gray, D., & Zaneveld, J. R. V. (1997). Absorption and attenuation of visible and near-infrared light in water: dependence on temperature and salinity. *Applied Optics, 36*, 6035-6046. DOI: 10.1364/AO.36.006035.
30. Zhang, X., Hu, L., & He, M.-X. (2009). Scattering by pure seawater: Effect of salinity. *Optics Express, 17*, 5698-5710. DOI: 10.1364/OE.17.005698. Baseline pure-water scattering reference; unchanged by this patch.

패키지 내부 생성·검증 보고서, source hash, 원자료 파일은 `03_SOURCE_AND_PROVENANCE/`, `04_REPRODUCIBILITY_ARCHIVES/`, `07_REFERENCES/`에 보존했다.
