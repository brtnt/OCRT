# OCRT ↔ reference model underwater consistency validation data — 2026-07-26

Prepared for the paper session. This is the complete set of underwater radiative-transfer
consistency data produced in that session. The reference model is OSOAA.

---

## 1. Files

| file | contents |
|---|---|
| `OCRT_OSOAA_underwater_validation_2026-07-26.csv` | per-direction raw data, 11,154 rows × 37 columns |
| `OCRT_OSOAA_underwater_summary_2026-07-26.csv` | per-combination score summary, 50 rows |
| `figures/` | all scatter plots made in the session |
| `raw/` | per-series source CSV before merging |

---

## 2. Common settings

| item | value |
|---|---|
| viewing zenith angle | 0–60°, 5° step (13) |
| relative azimuth | 0–180°, 15° step (13) |
| directions | 169 per combination |
| solar zenith angle | 40° unless stated |
| wind speed | 3 m/s unless stated |
| atmosphere | pressure 1013.25 hPa, aerosol optical depth 0, gas absorption off |
| sea-water refractive index | 1.34 fixed |
| in-water nodes | 64 |
| sun glint | decoupled |

---

## 3. Columns of the raw data

**Classification**

| column | meaning |
|---|---|
| `campaign` | `cdom`, `tsm`, `chl`, `cdom_scan_wind`, `cdom_scan_sza` |
| `ocrt_build` | `v5`, `v11` (chlorophyll), `v3` (control) |
| `note` | description of the combination |

**Conditions**

`chl_mg_m3`, `tsm_g_m3`, `adom440_m_inv`, `band_nm`, `sza_deg`, `vza_deg`, `raa_deg`, `wind_ms`,
`n_mu_water`

**Medium** (may be empty depending on the series)

`a_total_m_inv`, `b_total_m_inv`, `omega`, `g_eff`, `tau_water`, `sea_depth_m`,
`tau_rayleigh_atm`, `ocrt_trunc_A`, `osoaa_trunc_A`

**Reflectances** — the eight key columns.

| column | meaning |
|---|---|
| `OCRT_rrs0minus_I/Q/U` | OCRT remote-sensing reflectance below the surface |
| `OSOAA_rrs0minus_I/Q/U` | the same quantity from the reference model |
| `OCRT_Rrs0plus_I/Q/U` | OCRT remote-sensing reflectance above the surface |
| `OSOAA_Rrs0plus_I/Q/U` | the same quantity from the reference model |

**Diagnostics**

`black_frac_I/Q/U` — size ratio of the black-ocean component subtracted to obtain the
above-surface value. Read Section 6.

---

## 4. The two levels are produced differently

| level | OCRT | reference model |
|---|---|---|
| below surface `rrs(0⁻)` | directly from one grid run | one run; level 27 read at the **refracted underwater angles** and divided by the downward irradiance of the same run |
| above surface `Rrs(0⁺)` | directly from one grid run | **two runs and a subtraction**: (level26_actual − level26_black) ÷ ed26 |

Level 26 above the surface is dominated by sky light reflected by the surface (96 % of the
intensity, 99.9 % of Q). A second run with a black ocean under the same atmosphere, surface and
azimuth must therefore be subtracted to isolate the water-leaving light.

Azimuth convention: `phi(reference) = raa(OCRT) + 180°`.

---

## 5. Error normalization

Each component is normalized by **the maximum magnitude of the reference over that combination**.

```
error[%] = |OCRT − ref| / max|ref over the 169 directions| × 100
```

- **Do not normalize across mixed combinations.** Signal magnitudes differ by up to a factor of
  200 between combinations, so the largest one would dominate the denominator. The summary table
  is computed per combination.
- The degree of linear polarization is `sqrt(Q² + U²) / I` and is compared as an absolute
  difference (percentage points).

---

## 6. Caution when using the above-surface scores

The reference above-surface value is the difference of two large values. Where the water-leaving
light is small the residual becomes small and the relative error is inflated.

| black-subtraction share (I, median) | combinations | above-surface intensity error | below-surface intensity error |
|---|---|---|---|
| 0.11–0.13 | CDOM 0.01, 412·443·490 nm | 0.12–0.24 % | 0.08–0.12 % |
| 0.30–0.52 | CDOM 0.1, 412·490·555 nm | 0.42–0.89 % | 0.05–0.12 % |
| 0.74–0.79 | all 660 nm levels | 3.37–4.21 % | 0.12–0.34 % |
| 0.88–0.92 | CDOM 1.0, 412·443 nm | 8.47–13.78 % | 0.13 % |

Below the surface, where no subtraction is needed, the error is flat regardless of the share.
**When quoting above-surface numbers in a paper, also give `black_frac_*`, or use only
combinations with a share below 0.3.**

---

## 7. Scores by series (below-surface intensity)

| series | combinations | directions | mean error range |
|---|---|---|---|
| CDOM (3 levels × 5 bands) | 15 | 2,535 | 0.04–0.34 % |
| TSM (2 levels × 3 bands) | 6 | 1,014 | 0.05–0.49 % |
| chlorophyll (3 levels × 3 bands) | 9 | 1,521 | 0.13–0.25 % |

The absolute difference of the degree of linear polarization is 0.04–0.16 percentage points on
average in all three series.

---

## 8. Notes per series

**CDOM** — the simplest case. There are no particles, so no phase truncation applies. The v3
control is included, showing the effect of the v4·v5 corrections (above-surface U maximum error
15.9 % → 0.78 %).

**TSM** — one mineral phase function was given to the reference model as is. The truncation
coefficient is below the threshold, so **neither code truncates**. The water optical depth follows
the transport-optical-depth convention.

**Chlorophyll** — three things to know.

1. **Phytoplankton scattering is off.** The species phase functions do not fit into the L = 200
   expansion, and the reconstructed phase becomes negative. Chlorophyll currently affects
   **absorption only**; the only particle phase function is detritus.
2. Hence **the species choice does not affect the result**. Care is needed when a paper mentions
   phytoplankton species.
3. The OCRT side used the fixed-bulk bypass because the constituent path does not support
   truncation. The reproduction limit of that bypass is **0.18 %**, so the 0.13–0.25 % above is of
   the same size as that floor; the real difference between the codes may be smaller.

**Wind and solar zenith angle scans** — CDOM 0.1 basis. At wind 10 m/s the below-surface maximum
is 0.76 %, and at SZA 60° it is 1.05 %. At SZA 0° the U of both codes is 0, so the relative error
is meaningless and **must be excluded from the scores**.

---

## 9. Figures

| file | contents |
|---|---|
| `cdom_all_IQU_v5.png`, `cdom_all_IQU_v3.png` | CDOM six-component 1:1 scatter plots |
| `cdom_all_DoLP.png` | CDOM degree of linear polarization |
| `cdom_all_v3_vs_v5.png` | v3·v5 error by viewing zenith angle |
| `cdom_subtraction_diag.png` | diagnostic of the black-subtraction cancellation |
| `tsm_IQU.png`, `tsm_DoLP.png` | TSM |
| `tsm_depth_rule.png` | water optical depth convention |
| `chl_final_IQU.png`, `chl_final_DoLP.png` | chlorophyll |
| `chl_final_summary.png` | effect of the water-truncation defect fix and scores of the nine combinations |
| `scan_wind.png`, `scan_sza.png`, `scan_sza0_symmetry.png` | wind and SZA scans |
| `eap_oscillation.png`, `eap_new_vs_old.png` | EAP phase-data diagnostics |

---

# OCRT ↔ 참조모델 수중 정합 검증 자료 — 2026-07-26

논문 세션 전달용. 그 세션에서 만든 수중 복사전달 정합 검증 자료 전부다. 참조 모델은 OSOAA 다.

---

## 1. 들어 있는 파일

| 파일 | 내용 |
|---|---|
| `OCRT_OSOAA_underwater_validation_2026-07-26.csv` | 방향별 원자료 11,154행 × 37열 |
| `OCRT_OSOAA_underwater_summary_2026-07-26.csv` | 조합별 성적 요약 50행 |
| `figures/` | 그 세션에서 만든 산포도 전부 |
| `raw/` | 계열별 원본 CSV(합치기 전) |

---

## 2. 공통 설정

| 항목 | 값 |
|---|---|
| 관측천정각 | 0~60도, 5도 간격(13개) |
| 상대방위각 | 0~180도, 15도 간격(13개) |
| 방향 수 | 조합당 169개 |
| 태양천정각 | 40도(별도 표기 없으면) |
| 풍속 | 3 m/s(별도 표기 없으면) |
| 대기 | 기압 1013.25 hPa, 에어로졸 광학두께 0, 기체 흡수 끔 |
| 해수 굴절률 | 1.34 고정 |
| 수중 절점 | 64 |
| 선글린트 | 분리(decouple) |

---

## 3. 원자료 열 설명

**분류**

| 열 | 뜻 |
|---|---|
| `campaign` | `cdom`, `tsm`, `chl`, `cdom_scan_wind`, `cdom_scan_sza` |
| `ocrt_build` | `v5`, `v11`(엽록소), `v3`(대조군) |
| `note` | 조합 설명 |

**조건**

`chl_mg_m3`, `tsm_g_m3`, `adom440_m_inv`, `band_nm`, `sza_deg`, `vza_deg`, `raa_deg`, `wind_ms`,
`n_mu_water`

**매질**(계열에 따라 비어 있을 수 있다)

`a_total_m_inv`, `b_total_m_inv`, `omega`, `g_eff`, `tau_water`, `sea_depth_m`, `tau_rayleigh_atm`,
`ocrt_trunc_A`, `osoaa_trunc_A`

**반사도** — 여덟 열이 핵심이다.

| 열 | 뜻 |
|---|---|
| `OCRT_rrs0minus_I/Q/U` | OCRT 의 수면 아래 원격반사도 |
| `OSOAA_rrs0minus_I/Q/U` | 참조 모델의 같은 양 |
| `OCRT_Rrs0plus_I/Q/U` | OCRT 의 수면 위 원격반사도 |
| `OSOAA_Rrs0plus_I/Q/U` | 참조 모델의 같은 양 |

**진단**

`black_frac_I/Q/U` — 수면 위 값을 얻을 때 뺀 흑색 해수 성분의 크기 비다. 6절을 반드시 읽을 것.

---

## 4. 두 준위의 산출 방식이 다르다

| 준위 | OCRT | 참조 모델 |
|---|---|---|
| 수면 아래 `rrs(0⁻)` | 격자 실행 1회에서 직접 | 실행 1회. 준위 27 을 **굴절된 수중 각도**에서 읽고 같은 실행의 하향 조도로 나눔 |
| 수면 위 `Rrs(0⁺)` | 격자 실행 1회에서 직접 | **실행 2회와 차감.** (준위26_실제 − 준위26_흑색) ÷ ed26 |

수면 위 준위 26 은 해면이 되비친 하늘광이 지배한다(세기의 96 %, Q 의 99.9 %). 그래서 같은 대기·해면·
방위각으로 흑색 해수를 한 번 더 돌려 빼야 수출광만 남는다.

방위각 규약은 `phi(참조) = raa(OCRT) + 180도` 다.

---

## 5. 오차 정규화 규약

성분별로 **그 조합의 참조 최대 크기**로 나눈다.

```
error[%] = |OCRT − ref| / max|ref over the 169 directions| × 100
```

- **조합을 섞어 한 번에 정규화하면 안 된다.** 조합마다 신호 크기가 200배까지 차이 나므로 큰 조합이 분모를
  독점한다. 요약표는 조합별로 계산해 두었다.
- 선형편광도는 `sqrt(Q² + U²) / I` 로 정의하고 절대차(퍼센트포인트)로 본다.

---

## 6. 수면 위 성적을 쓸 때의 주의

참조 모델의 수면 위 값은 두 큰 값의 차이다. 수출광이 작은 조합에서는 잔차가 작아져 상대오차가 부풀려진다.

| 흑색 차감 비중(I, 중앙값) | 해당 조합 | 수면 위 세기 오차 | 수면 아래 세기 오차 |
|---|---|---|---|
| 0.11~0.13 | CDOM 0.01, 412·443·490 nm | 0.12~0.24 % | 0.08~0.12 % |
| 0.30~0.52 | CDOM 0.1, 412·490·555 nm | 0.42~0.89 % | 0.05~0.12 % |
| 0.74~0.79 | 660 nm 전 단계 | 3.37~4.21 % | 0.12~0.34 % |
| 0.88~0.92 | CDOM 1.0, 412·443 nm | 8.47~13.78 % | 0.13 % |

차감이 필요 없는 수면 아래는 비중과 무관하게 평평하다. **수면 위 숫자를 논문에 실을 때는 `black_frac_*` 를
함께 적거나, 비중 0.3 미만 조합만 쓰기를 권한다.**

---

## 7. 계열별 성적(수면 아래 세기)

| 계열 | 조합 수 | 방향 수 | 평균 오차 범위 |
|---|---|---|---|
| CDOM(3단계 × 5밴드) | 15 | 2,535 | 0.04~0.34 % |
| 총부유물(2단계 × 3밴드) | 6 | 1,014 | 0.05~0.49 % |
| 엽록소(3단계 × 3밴드) | 9 | 1,521 | 0.13~0.25 % |

선형편광도 절대차는 세 계열 모두 평균 0.04~0.16 퍼센트포인트다.

---

## 8. 계열별로 알아 둘 사항

**CDOM** — 가장 단순하다. 입자가 없어 위상 절단이 걸리지 않는다. v3 대조군이 함께 들어 있어 v4·v5 수정
효과를 볼 수 있다(수면 위 U 최대 오차 15.9 % → 0.78 %).

**총부유물** — 무기물 위상 한 개를 참조 모델에 그대로 넣었다. 절단 계수가 문턱 미달이라 **두 코드 모두
절단하지 않는다.** 해수 광학두께는 수송 광학두께 규약으로 잡았다.

**엽록소** — 세 가지를 알아야 한다.

1. **식물플랑크톤 산란이 꺼져 있다.** 종별 위상이 L = 200 전개에 담기지 않아 되살린 위상이 음수가 되는
   문제가 확인됐다. 현재 엽록소는 **흡수에만** 영향을 준다. 입자 산란 위상은 쇄설물 하나다.
2. 그 결과 **종 선택이 결과에 관여하지 않는다.** 논문에서 식물플랑크톤 종을 언급할 때 주의가 필요하다.
3. OCRT 쪽은 고정 벌크 우회 경로로 돌렸다. 성분 경로가 절단을 지원하지 않기 때문이다. 이 우회의 재현 정확도
   한계가 **0.18 %** 이므로, 위 0.13~0.25 % 는 그 바닥값과 같은 크기다. 실제 코드 간 차이는 더 작을 수 있다.

**풍속·태양천정각 변화** — CDOM 0.1 기준이다. 풍속 10 m/s 에서도 수면 아래 최대 0.76 %, 태양천정각 60도에서
최대 1.05 % 다. 태양천정각 0도의 U 는 두 코드 모두 0 이라 상대오차가 무의미하므로 **성적에서 제외해야 한다.**

---

## 9. 그림 목록

| 파일 | 내용 |
|---|---|
| `cdom_all_IQU_v5.png`, `cdom_all_IQU_v3.png` | CDOM 6성분 1:1 산포도 |
| `cdom_all_DoLP.png` | CDOM 선형편광도 |
| `cdom_all_v3_vs_v5.png` | 관측천정각별 v3·v5 오차 |
| `cdom_subtraction_diag.png` | 흑색 차감 상쇄 진단 |
| `tsm_IQU.png`, `tsm_DoLP.png` | 총부유물 |
| `tsm_depth_rule.png` | 해수 광학두께 규약 |
| `chl_final_IQU.png`, `chl_final_DoLP.png` | 엽록소 |
| `chl_final_summary.png` | 물 절단 결함 수정 효과와 아홉 조합 성적 |
| `scan_wind.png`, `scan_sza.png`, `scan_sza0_symmetry.png` | 풍속·태양천정각 변화 |
| `eap_oscillation.png`, `eap_new_vs_old.png` | EAP 위상 자료 진단 |
