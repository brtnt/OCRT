# OCRT–reference cross-check residual scatter plots (2026-08-23)

Conditions: 555 nm, SZA 40°, wind 3 m/s, n_mu_water 48. OCRT = FR631 canonical data + 08-21
overlay runtime (the same run as the cross_check_fr631 harness). Reference = OSOAA with NT_SEA
640; the red-clay case was regenerated with the canonical FR631 phase in
`runs/nt640_osoaa_tsm555_fr631`.

| figure | contents | source CSV |
|---|---|---|
| fig1_rrs_underwater_residual.png | underwater rrs(0−) residual vs VZA (RAA 90°, 12 points × 3 series: pure / red_clay direct / red_clay moment) | rrs_underwater_raa90.csv |
| fig2_Rrs_I_residual_156geom.png | above-water Rrs(0+) I residual, 156 geometries (VZA 5–60 × RAA 0–180°), colour = RAA | Rrs_above_156geom.csv |
| fig3_Rrs_QU_residual.png | Rrs Q·U residual, %/max\|I\| convention, 156 geometries × 2 cases | Rrs_above_156geom.csv |
| fig4_atm_rhoI_residual.png | atmospheric (Rayleigh, black Fresnel ocean) TOA ρ_I residual, 84 points in 5 azimuth planes | atm_rhoI_bfo.csv |

Gates: rrs ±0.15 %, Rrs ±0.60 %, atmosphere ±0.5 % (unchanged harness values). The κ-residual
envelope ±0.41 % is shaded in fig1.

## 1:1 scatter plots (formal validation plots — x = reference, y = OCRT, diagonal = 1:1)

| figure | contents |
|---|---|
| figS1_rrs_1to1.png | underwater rrs(0−) 1:1 + residual strip (pure 12 points, red_clay 12 points, RAA 90°) |
| figS2_RrsI_1to1.png | Rrs(0+) I 1:1 + residual strip (156 geometries per case) |
| figS3_RrsQU_1to1.png | Rrs Q·U 1:1 + residual strip (%/max\|I\|) |
| figS4_atm_rhoI_1to1.png | atmospheric ρ_I 1:1 + residual strip (84 points; marker shape = azimuth plane) |

fig1–fig4 (residual-versus-angle structure) are kept as supporting material of the 1:1 plots.
Dotted lines around the diagonal = gate envelope (y = x × (1 ± g)).

## Validation after the fix (2026-08-25, T_wa U-column folding fix)

| figure | contents | source CSV |
|---|---|---|
| figS5_RrsIQU_1to1_postfix.png | Rrs I·Q·U 1:1 after the fix (156 geometries per case; grey = before the fix) | Rrs_above_156geom_postfix.csv |
| figQ3_interface_transfer.png | (diagnostic) interface Q transfer ratio — the defect shape before the fix | Q_diagnosis_156.csv |

Summary after the fix: Q max deviation 5.28 → 0.09 (red_clay) and 7.12 → 0.13 (pure) %/max|I|;
U unchanged (0.091); I keeps the κ residual (+0.24 / +0.65 % median; separate issue).

---

# OCRT–참조모델 교차검증 잔차 산포도 (2026-08-23)

조건: 555 nm, SZA 40°, 풍속 3 m/s, n_mu_water 48. OCRT = FR631 표준 자료 + 08-21 overlay 런타임
(cross_check_fr631 하니스와 같은 실행). 참조 = OSOAA, NT_SEA 640; red_clay 는 표준 FR631 위상으로
`runs/nt640_osoaa_tsm555_fr631` 에서 재생성.

| 그림 | 내용 | 원본 CSV |
|---|---|---|
| fig1_rrs_underwater_residual.png | 수중 rrs(0−) 잔차 vs VZA (RAA 90°, 12점 × 3계열: pure / red_clay direct / red_clay moment) | rrs_underwater_raa90.csv |
| fig2_Rrs_I_residual_156geom.png | 수면 위 Rrs(0+) I 잔차, 156 기하(VZA 5–60 × RAA 0–180°), 색 = RAA | Rrs_above_156geom.csv |
| fig3_Rrs_QU_residual.png | Rrs Q·U 잔차, %/max\|I\| 규약, 156 기하 × 2 케이스 | Rrs_above_156geom.csv |
| fig4_atm_rhoI_residual.png | 대기(Rayleigh, 흑색 Fresnel 해양) TOA ρ_I 잔차, 방위평면 5개 84점 | atm_rhoI_bfo.csv |

게이트: rrs ±0.15 %, Rrs ±0.60 %, 대기 ±0.5 % (기존 하니스 값 유지). κ 잔차 포락선 ±0.41 % 는 fig1 에
음영으로 표시.

## 1:1 산포도 (정식 검증 산포도 — x = 참조값, y = OCRT 값, 대각선 = 1:1)

| 그림 | 내용 |
|---|---|
| figS1_rrs_1to1.png | 수중 rrs(0−) 1:1 + 잔차 스트립 (pure 12점, red_clay 12점, RAA 90°) |
| figS2_RrsI_1to1.png | Rrs(0+) I 1:1 + 잔차 스트립 (케이스당 156 기하) |
| figS3_RrsQU_1to1.png | Rrs Q·U 1:1 + 잔차 스트립 (%/max\|I\|) |
| figS4_atm_rhoI_1to1.png | 대기 ρ_I 1:1 + 잔차 스트립 (84점, 방위평면 5개 = 마커 형태) |

fig1~fig4(잔차–각도 구조도)는 1:1 산포도의 보조 자료로 유지한다. 대각선 주변 점선 = 게이트 포락선
(y = x × (1 ± g)).

## 수정 후 검증 (2026-08-25, T_wa U 열 접기 수정)

| 그림 | 내용 | 원본 CSV |
|---|---|---|
| figS5_RrsIQU_1to1_postfix.png | 수정 후 Rrs I·Q·U 1:1 산포도(케이스당 156 기하, 회색 = 수정 전) | Rrs_above_156geom_postfix.csv |
| figQ3_interface_transfer.png | (진단) 계면 Q 전달비 — 수정 전 결함 형상 | Q_diagnosis_156.csv |

수정 후 요약: Q 최대 편차 5.28 → 0.09(red_clay)·7.12 → 0.13(pure) %/max|I|, U 무회귀(0.091), I 는 κ 잔차
유지(+0.24 / +0.65 % 중앙값, 별건).
