# OCRT-OSOAA 교차검증 잔차 산포도 (2026-08-23)

조건: 555 nm, SZA 40°, 풍속 3 m/s, n_mu_water 48. OCRT = FR631 canonical 자료 + 08-21 오버레이 런타임(cross_check_fr631 하니스와 동일 실행). OSOAA 기준 = NT_SEA 640, red_clay는 canonical FR631 위상으로 재생성한 `runs/nt640_osoaa_tsm555_fr631`.

| 그림 | 내용 | 소스 CSV |
|---|---|---|
| fig1_rrs_underwater_residual.png | 수중 rrs(0−) 잔차 vs VZA (RAA 90°, 12점×3계열: pure / red_clay direct / red_clay moment) | rrs_underwater_raa90.csv |
| fig2_Rrs_I_residual_156geom.png | 수면 위 Rrs(0+) I 잔차, 156기하(VZA 5–60×RAA 0–180°), 색=RAA | Rrs_above_156geom.csv |
| fig3_Rrs_QU_residual.png | Rrs Q·U 잔차, %/max\|I\| 규약, 156기하×2케이스 | Rrs_above_156geom.csv |
| fig4_atm_rhoI_residual.png | 대기(Rayleigh-BFO) TOA ρ_I 잔차, 5개 방위평면 84점 | atm_rhoI_bfo.csv |

게이트: rrs ±0.15%, Rrs ±0.60%, 대기 ±0.5% (기존 하니스 값 유지). κ-잔차 포락선 ±0.41%는 fig1에 음영으로 표시.

## 1:1 산포도 (정식 검증 산포도 — x=OSOAA 참값, y=OCRT 추정값, 대각선=1:1)

| 그림 | 내용 |
|---|---|
| figS1_rrs_1to1.png | 수중 rrs(0−) 1:1 + 잔차 스트립 (pure 12점, red_clay 12점, RAA 90°) |
| figS2_RrsI_1to1.png | Rrs(0+) I 1:1 + 잔차 스트립 (케이스당 156기하) |
| figS3_RrsQU_1to1.png | Rrs Q·U 1:1 + 잔차 스트립(%/max|I|) |
| figS4_atm_rhoI_1to1.png | 대기 ρ_I 1:1 + 잔차 스트립 (84점, 5개 방위평면 = 마커 형태) |

fig1~fig4(잔차-각도 구조도)는 1:1 산포도의 보조 자료로 유지한다. 대각선 주변 점선 = 게이트 포락선(y=x×(1±g)).

## 수정 후 검증 (2026-08-25, T_wa U-열 접기 수정)

| 그림 | 내용 | 소스 CSV |
|---|---|---|
| figS5_RrsIQU_1to1_postfix.png | 수정 후 Rrs I·Q·U 1:1 산포도(케이스당 156기하, 회색=수정 전) | Rrs_above_156geom_postfix.csv |
| figQ3_interface_transfer.png | (진단) 계면 Q 전달비 — 수정 전 상태의 결함 형상 | Q_diagnosis_156.csv |

수정 후 요약: Q dev max 5.28→0.09(red_clay)·7.12→0.13(pure) %/max|I|, U 무회귀(0.091), I는 κ-잔차 유지(+0.24/+0.65% 중앙, 별건).
