# BFO(흑체 프레넬 해양) Rayleigh IQU 비교 가이드 — OCRT·OSOAA 기준자료 및 RTSOS 삼자비교용 상세 명세 (2026-08-25)

용도: 첨부 `atm_rhoIQU_bfo.csv`(84기하, OCRT vs OSOAA의 TOA ρ_I·ρ_Q·ρ_U)를 제3의 RT 코드(RTSOS)와 비교할 때 필요한 모든 옵션·정의·규약을 이 문서 하나로 확정한다. 수치는 전부 저장된 검증 런에서 나온 것이며(2026-08-18 F9 대기검증 1단계 실행분), 본 문서 작성 시 신규 RT 실행은 없다 — 설정값은 소스 코드·실행 로그·출력 헤더에서 재확인했다(§2 근거 열).

## 1. 케이스 정의 (물리 구성)

Rayleigh 단독 대기 + 거친 프레넬 해면 + 흑체 해양(수중 복사 전달 없음)의 TOA 상향 반사도. 직달 선글린트(태양 직달빔의 경면 반사)는 **양 코드 모두 제외**하고, 대기 산란광의 해면 반사(확산 글린트)는 포함한다. 흡수 기체·에어로졸 없음.

## 2. 실행 구성 (전 항목, 근거 병기)

| 항목 | 값 | OCRT 근거 | OSOAA 근거 |
|---|---|---|---|
| 파장 | 555 nm | `--wavelength 555` | `-Wa 0.555` |
| SZA | 40° | `--sza 40` | `-Thetas 40` (출력 RAD_UsedAngles: 40.000) |
| τ_R (레일리 광학두께) | **0.0935485118** | 내부 Bodhaine 1999 full model(1013.25 hPa, lat 45, alt 0); 실행출력 `[TAUR] wl=555 tau_R=0.0935485118` | `-AP.MOT 0.0935485118`; PROFILE_ATM 표면 0.09355 |
| 공기 depol (δ_air) | **0.0279** (양 코드 동일) | `rt_solver.c:644 const double depol = 0.0279` + 실행출력 `depol=0.02790000` | `inc/OSOAA.h:266 #define CTE_MDF_AIR 0.0279` |
| 흡수 기체 | 0 | `--gas-column-{h2o,o3,no2,o2,co2,ch4} 0` | (없음 — Rayleigh만) |
| 에어로졸 | 0 | (미지정) | PROFILE_ATM AER_PC=0 전층 |
| 해면 | Cox–Munk 등방, wind 3 m/s | `--surface black_fresnel_ocean --wind-speed 3` | `-SEA.Wind 3` |
| σ² 법칙 | **σ² = 0.003 + 0.00512·max(0.01, W)** → W=3에서 0.01836 | `shared/surface.c:58` | `OSOAA_SURF_MATRICES.F:2738` (P5 parity patch, 동일식) |
| 해면 shadowing | **Sancer(1969) bistatic, 반사 경로 상시 적용** (양 코드 동일 구현) | `shared/surface.c` Sancer S_bi | P6 patch (`SURF_MATRICES.F:3377` "identical form to OCRT") |
| 직달 글린트 | **제외** | `--decouple-sunglint` | env `OSOAA_NO_DIRECT_GLINT=1` (parity patch) |
| 수중 굴절률 | n_w = 1.34 | `--n-water 1.34` | (동일; RAD_UsedAngles 수중 태양천정각 28.665° = asin(sin40°/1.34) 일치) |
| 흑체 해양 구현 | 수중 RT 자체 없음 | `black_fresnel_ocean` 모드 | 사용자 해수 프로파일 a=29 m⁻¹, b=1e-9, depth 1 m (`prof.txt`) — 불투명 수체 |
| 대기 이산화 | — | n_layers 40 (기본) | 26층, Δτ≈0.0036 균일, TOA 300 km |
| 각 구적 | — | n_mu 24/반구(기본) + view node 삽입(OSOAA IMUS 규칙: 정렬 위치 삽입·GL 노드와 dedup) | NbGauss 48 (`NB_GAUSS_ANGLES: 48`) |
| 푸리에 모드 | m_max 2 (레일리 단독 자동 디스패치 — 레일리 위상행렬은 m≤2가 정확) | `main.c` m_max dispatch | OS_NM 296 (내부) |
| 산란 차수 | max_orders 20, SOS tol 1e-7 | `rt_types.h` 기본값 | OS_NB 200 |
| 출력 격자 | VZA 0–60°×5°, RAA 0–345°×15° (fullgrid) | `--lut-vza-max 60 --lut-vza-step 5 --lut-raa-step 15` | View.Phi 5종(0/60/90/120/180) × 가우스 VZA(±) |

OCRT 재현 명령(전문):

```
./ocrt --surface black_fresnel_ocean --wind-speed 3 --decouple-sunglint \
  --n-water 1.34 --wavelength 555 --sza 40 \
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 \
  --lut-vza-max 60 --lut-vza-step 5 --lut-raa-step 15 --output-full-grid <csv>
```

OSOAA: 표준 인자 + `-AP.MOT 0.0935485118 -SEA.Wind 3 -SEA.Depth 1 -SEA.Ind 1.34 -Thetas 40 -Wa 0.555 -NbGauss 48` + 흑수 프로파일(ExtData) + env `OSOAA_NO_DIRECT_GLINT=1`, View.Phi ∈ {0,60,90,120,180}로 5회 실행.

## 3. 물리량 정의 (정규화)

- **ρ_X ≡ πL_X / (E_sun · cosθ_s)**, X ∈ {I,Q,U}. TOA 상향 복사휘도 L을 대기 외 태양 복사조도와 태양천정각 코사인으로 정규화한 반사도. 무차원.
- OCRT `fullgrid.csv`의 `rho_I/rho_Q/rho_U`가 바로 이 정의.
- OSOAA와의 등가: Standard 출력 REFL = πL/E_d(TOA)이고 TOA에서 E_d = E_sun·cosθ_s이므로 REFL ≡ ρ_I. Advanced 출력(RESLUM_Adv_UP, level 0) I,Q,U는 πL/E_sun 정규화이므로 **cosθ_s(=cos40°)로 나눠** ρ로 환산했다(본 CSV의 osoaa_* 열).
- Stokes 기준면: 양 코드 모두 관측 자오면(meridian plane) 기준. V 성분은 양 코드 모두 미사용(레일리+프레넬 체계에서 0).

## 4. 각도·부호 규약 (RTSOS 매핑 시 가장 주의할 부분)

- VZA: 천정각, TOA 상향 관측. CSV의 `vza_signed`는 OSOAA 출력 관례(한 실행이 φ와 φ+180 두 반평면을 ±VZA로 담음)의 부호를 유지한 것으로, 물리 VZA는 |vza_signed|.
- OCRT RAA 규약: **raa=180°가 태양과 같은 반평면(글린트 쪽), raa=0°가 반대 반평면** (`rt_raa_convention.h`: RAA=180 direct-glint branch; SZA=VZA·RAA=180에서 정확 후방산란 Θ=180°).
- OSOAA φ 규약(출력 헤더 명문): φ=180° = 위성·태양 같은 반평면, φ=0° = 반대 반평면.
- **실측 확정 매핑: raa_OCRT = (φ_OSOAA + 180°) mod 360°.** 본 병합에서 U 부호로 판정했다: 후보 (φ+180) 채택 시 U 잔차 max 0.091 %/max|I|, 거울 후보 (180−φ) 사용 시 86.5 %/max|I|로 파탄 — I·Q는 φ-우함수라 두 후보를 구분하지 못하고 U(φ-기함수)만이 판정한다. (기존 2026-08-18 F9 보고서의 "raa=180−φ" 기재는 I 기반 판정이었고 거울쌍 한쪽을 임의 선택한 것 — I·Q 비교에는 유효하나 **부호 있는 U 비교에는 (φ+180)이 옳다.** 수중장 WL 대조에서 확립된 φ=raa+180 규약과도 일치.)
- 주의: 두 코드의 문서상 규약(둘 다 "180=같은 반평면")대로라면 항등 매핑이어야 하나 실측은 +180 이동이다. 즉 **어느 한쪽의 문서 규약과 실제 출력이 반대**다(어느 쪽인지는 제3 코드 대조로 판정 가능 — RTSOS 비교의 부수 성과가 된다). RTSOS 매핑 시 문서를 믿지 말고 §6의 단일산란 폐형식 앵커로 직접 판정하라.
- U 부호: 위 매핑 적용 시 양 코드 U가 부호까지 일치(잔차 0.091 %/max|I|). 주평면(φ 0/180)에서 U=0 확인(수치상 ≤1e-17).

## 5. CSV 컬럼 사전 (`atm_rhoIQU_bfo.csv`, 84행)

| 열 | 의미 |
|---|---|
| plane_phi | OSOAA 실행의 View.Phi (0/60/90/120/180) |
| osoaa_phi_view | 해당 행의 실제 관측 방위(φ_view): vza_signed>0이면 plane_phi, <0이면 plane_phi+180 |
| vza_signed | OSOAA 부호 VZA (±5…±60, 5° 간격; 물리 VZA=|·|) |
| ocrt_raa_deg | 쌍을 이룬 OCRT 격자 RAA = (osoaa_phi_view+180) mod 360 |
| ocrt_rho_I/Q/U | OCRT TOA ρ (fullgrid 직독, 10유효자리) |
| osoaa_rho_I/Q/U | OSOAA TOA ρ (Adv_UP level 0을 cos40°로 나눔; 가우스각→부호 VZA축 선형보간) |
| dev_I_pct | 100·(ocrt/osoaa − 1) [%] |
| dev_Q_pctmaxI, dev_U_pctmaxI | 100·(ocrt−osoaa)/max\|osoaa_rho_I\| [%/max\|I\|], max\|I\|=0.093915 |

행 구성: 5개 평면의 VZA>0 분기 + 주평면(0/180)의 VZA<0 분기 = 7그룹 × 12 VZA = 84. 주평면 두 실행(p0, p180)은 같은 물리 기하를 서로 다른 실행에서 재표본한 교차확인 쌍이다.

보간 프로토콜: OSOAA는 가우스 각도(48노드/반구)에서 출력되므로, 부호 VZA축에서 성분별 선형보간(np.interp)으로 OCRT 정수각(5° 격자)에 맞췄다. OCRT는 view-as-node 방식이라 보간 없음. 5° 간격의 레일리 장에서 선형보간 오차는 잔차 통계(≤0.29%)에 포함된 상태다.

## 6. 정합 실측치 및 게이트 (이 CSV가 담고 있는 상태)

| 성분 | median | max\|·\| | 단위 | 게이트(제안) |
|---|---|---|---|---|
| ρ_I | −0.148 % | 0.286 % | % | mean\|dev\| ≤ 0.5 (기존 하니스) |
| ρ_Q | +0.001 | 0.048 | %/max\|I\| | ≤ 0.20 |
| ρ_U | −0.000 | 0.091 | %/max\|I\| | ≤ 0.15 |

- I의 −0.15% 계열 잔차는 기지 항목(후보: OSOAA 26층 이산화 vs OCRT 40층; F9 보고서 §결과). Q·U는 사실상 완전 정합.
- 검증 게이트(병합 스크립트 `merge_atm_iqu_bfo.py` 출력): G1 osoaa_I가 기존 I-only CSV(REFL 기반)와 2.6e-6 이내 재현, G2 ocrt_I가 기존 CSV와 인쇄 정밀도(5e-9) 이내 일치 + OCRT 거울 raa 노드 대칭 동시 검증, G3 U 분기 판정(위), G4 주평면 U=0. 전부 PASS.
- 1:1 산포도: `figS6_atm_rhoIQU_1to1.png` (x=OSOAA 참값, y=OCRT, 대각선=1:1, 하단 잔차 스트립).

## 7. RTSOS 삼자비교 매핑 지침 (RTSOS_RAYBENCH_STATUS_2026-08-23 미확정 항목 회답 포함)

1. **케이스 대응**: RTSOS monochromatic 코어 IPT=−101(수중 RT 없는 프레넬 하부경계) + RSR0P=0 = 본 BFO 케이스. 확인됨.
2. **depol**: **0.0279로 확정** (양 코드 소스·실행 로그 근거, §2). 0.0284가 아님. RTSOS .pmtx를 0.0279로 생성할 것.
3. **τ_R**: 555 nm에서 **0.0935485118** (Bodhaine 1999, 1013.25 hPa, lat 45). RTSOS 시연의 443 nm 값(0.23584)이 아니라 이 값으로 555 nm 런을 만들어야 본 CSV와 직접 비교된다.
4. **정규화**: RTSOS ESUN=π 설정 시 ρ=L/μ0 → 본 ρ 정의와 동일. 그대로 사용.
5. **바람·글린트**: 본 CSV는 **W=3 m/s + 직달 글린트 제외**다. RTSOS에서 (a) 직달 글린트 항을 제거할 수 있으면 W=3으로 직접 비교, (b) 제거 불가면 글린트 로브에서 먼 기하(반태양 반평면, 본 CSV의 ocrt_raa 0/60 근방 행)로 한정 비교하거나 W=0 flat run과의 별도 대조를 설계하라. W=0 델타 글린트의 이산 표현 규약(격자점 스파이크) 통일 문제는 본 CSV에는 해당 없음(직달 제외 상태).
6. **σ²·shadowing**: σ²=0.003+0.00512W는 삼자 동일 계수. 단 shadowing은 OCRT·OSOAA가 **Sancer(1969) bistatic 상시 적용**인 반면 RTSOS 기본은 Smith 1/(1+Λ₁+Λ₂)이므로 **모델이 다르다**. W=3 확산 글린트에서는 소차이지만, 비교 시 RTSOS shadowing off/on 양쪽을 돌려 민감도를 정량하고 명기할 것.
7. **방위·U 부호 규약**: §4의 실측 매핑 절차를 그대로 적용하라 — 문서 규약을 신뢰하지 말고, (i) I로 방위쌍 판정 → (ii) U 부호로 거울 분기 판정. 판정 앵커로는 레일리 단일산란 폐형식(ρ⁽¹⁾ = ω τ_R P(Θ)/(4(μ_s+μ_v)) × [1−exp(−τ_R(1/μ_s+1/μ_v))] × 회전행렬, 얇은 τ 극한)이 유용하다. RTSOS는 Siewert 기반이라 U 부호가 반대일 가능성이 있다 — 그 경우 U에 −1을 곱해 정렬하고 명기.
8. **비교 대상 열**: RTSOS 산출을 §3 정의의 ρ로 환산 후 osoaa_rho_* 열(참값측)과 대조하는 것을 권장(OSOAA가 본 체인의 기준 코드). dev 규약은 §5와 동일하게 I는 %, Q/U는 %/max|I|.

## 8. 파일 위치·재현

- 본 CSV·그림·병합 스크립트: `docs/validation/cross_check_fr631_20260823/` (atm_rhoIQU_bfo.csv, figS6_atm_rhoIQU_1to1.png, merge_atm_iqu_bfo.py). 기존 I-only `atm_rhoI_bfo.csv`는 그대로 보존(상위집합 관계).
- 원시 저장 런: 샌드박스 `05_VALIDATION/runs/atm_ocrt_bfo/fullgrid.csv`(OCRT), `runs/atm_osoaa_bfo_p{0,60,90,120,180}/`(OSOAA; Advanced_outputs/RESLUM_Adv_UP.txt level 0 사용). 로컬 미러: `code/validation_05/` 동일 경로.
- OCRT 바이너리는 2026-08-25 T_wa 수정판과 무관(BFO는 수중장이 없어 T_wa 미사용 — 수정 전후 fullgrid 비트 동일 확인됨).
