# 작업지시서: T_wa 푸리에 커널 U-열 부호접기 수정 (근천저 수출광 Q 결함) — 2026-08-25

대상: OCRT C/Python 페어 개발 트리 (MIGRATION_PKG_2026-08-19 계열). 본 지시서는 독립 세션이 이 문서만 읽고 이어갈 수 있도록 작성하며, **검증결과 파일을 동봉한 단일 압축 패키지**(`OCRT_WORKORDER_TWA_QFOLD_FIX_2026-08-25.zip`)로 배포한다 — §5의 동봉물과 §6의 절차만으로 수정 반영·검증이 완결된다. 상세 판별 이력은 `docs/validation/cross_check_fr631_20260823/Q_TRUNCATION_FINDINGS_KO.md` §8–9.

## 1. 결함 요약 (수정 완료됨 — 본 지시서는 반영·후속용)

- 증상: OSOAA 대비 수면 위 수출광 Rrs Q가 근천저에서 진폭 0.57×로 감쇠(VZA≤10°, wind 3). 수중장(0⁻)은 I/Q/U 완전 정합, U(0⁺)도 정합 — Q(0⁺)만 결손.
- 원인: `src/shared/surface.c`의 `surface_T_wa_coxmunk_fourier_kernel`(거친면 물→공기 투과 푸리에 연산자)에 **m>0 U-열 부호접기 누락**. 동일 접기가 `surface_T_aw_coxmunk_fourier_kernel`(공기→물)에는 존재하며 주석에 "OSOAA Step-7 storage/contraction 일치"로 명시되어 있음. 접기 없이는 m≥1 편광 모드의 스핀-2 폐쇄가 깨져 Q만 감쇠한다(콘볼루션 대수와 단위시험으로 확정).
- 수정: 접기를 헬퍼 `ocrt_twa_fold_ucol_m()`로 추가하고 **직접 경로와 FKC 캐시 경로 양쪽 공통**으로 적용(래퍼 말단 1곳 + FKC 분기 1곳). diff는 surface.c 내 "FIX 2026-08-25" 주석 블록으로 식별 가능. 패치 후 소스 SHA-256 `6a7da737aa906498…`.

## 2. 검증 배터리 (전부 통과, 2026-08-25)

| 게이트 | 결과 |
|---|---|
| 단위 게이트(신규 `tests/test_twa_spin2_gate.c`): 합성 m=2 장 → rough/flat 전달비 | (T_I, T_Q2, T_U2) = **(1.000, 1.000, 1.000)** @ VZA 2–40°, FKC ON/OFF 동일 |
| PSSA full-grid 스모크 (BFO, 555/SZA75) | CSV SHA **98e4b2433be300eb… 비트 불변** (T_wa 미사용 경로 무회귀) |
| wind=0 red_clay fullgrid | 수정 전과 **비트 동일** (평면 경로 불변, G1 성격) |
| red_clay w3 vs OSOAA(FR631 기준): Rrs Q | slope(VZA≤10) 0.573→**0.996**, dev max 5.28→**0.09**, mean 0.542→**0.031** %/max\|I\| |
| pure w3 vs OSOAA: Rrs Q | slope **0.994**, dev max 7.12→**0.13** %/max\|I\| |
| Rrs U (양 케이스) | max 0.091%/max\|I\| — **무회귀** |
| 공식 하니스(run_cross_check_fr631.sh, I-게이트) | pure PASS(+0.245%), rayleigh PASS(0.162%), red_clay I 잔차 −0.258/+0.659%는 기존 κ-이슈로 **별건 유지**(§4-c) |

## 3. 왜 기존 정합성 테스트가 이 결함을 통과시켰나 (수정 과정에서 무비용 판명)

이 성분을 보는 게이트가 하나도 없었다: (i) 기존 게이트는 전부 I 기반이고 I는 m=0 지배라 둔감(수정 후에도 I 변화 ≤0.02%p), (ii) U-행은 접기가 원래 불필요한 구조라 08-18 U 절대검증이 통과, (iii) G1 비트게이트·과거 비교의 상당수는 wind=0 평면 경로였고 이 커널은 wind>0 전용, (iv) Q-vs-OSOAA는 08-18에 "판정 부적합"으로 분류되어 게이트에서 제외돼 있었다. 즉 결함은 "거친면 × 편광 Q × 근천저" 교집합에만 살았고 그 교집합은 무게이트였다.

## 4. 후속 작업 (우선순위순)

- (a) **Python 페어 구현 점검**: `ocrt_py/coupling.py`·`surface.py`의 fourier_kernel 계열이 같은 접기 규약을 갖는지 확인, C와 동일 단위게이트(합성 m=2 장, rough/flat=(1,1,1)) 이식. C–Py 1e-15 정합 원칙 적용.
- (b) **기존 편광 캠페인 산출물 재평가**: 이 결함은 wind>0 수출광 Q(근천저)에 영향 — PSSA 캠페인·과거 8,045-run 중 Q를 소비한 산출물의 재실행 필요 범위 판정(스모크·I 산출물은 영향 없음).
- (c) **I-κ 잔차(별건)**: red_clay Rrs_I +0.66%/rrs −0.26%는 이 수정과 무관한 기존 κ-계열 항목. 게이트(0.15/0.60%) 재조정 또는 추가 환원 결정 대기(FINDINGS §2·§6).
- (d) **게이트 추가**: `tests/test_twa_spin2_gate.c`를 회귀 게이트에 편입, 하니스에 Q 게이트(예: near-nadir devQ ≤0.2%/max|I|) 신설.
- (e) G-b2 폐형식 검증에 근천저 투과 케이스(θ→0, T_Q/T_I=1) 추가 — 이번 결함류의 상시 감시.

## 5. 동봉물 (압축 패키지 구성)

```
OCRT_WORKORDER_TWA_QFOLD_FIX_2026-08-25.zip
├── OCRT_WORKORDER_TWA_QFOLD_FIX_2026-08-25.md   ← 본 문서
├── patch/
│   ├── surface.c                     패치 완료본 전체 (SHA-256 6a7da737aa906498…)
│   └── surface_twa_qfold.diff        수정 전→후 unified diff (47줄) — 다른 트리에 적용용
├── tests/
│   └── test_twa_spin2_gate.c         단위 게이트 소스 (합성 m=2 장, rough/flat=(1,1,1))
├── verify/
│   ├── EXPECTED_GATES.md             기대 게이트 값·허용오차 전체 명세
│   ├── verify_twa_qfold_fix.py       자동 검증기 (아래 §6)
│   ├── Rrs_above_156geom_postfix.csv 수정 직후 OCRT vs OSOAA, 156기하×2케이스×I/Q/U (기준 CSV)
│   ├── rrs_underwater_raa90_postfix.csv  수중 rrs(0⁻) 직접커널 2케이스 (무회귀 확인용)
│   └── atm_rhoI_bfo_postfix.csv      대기 ρ_I 84점 (무회귀 확인용)
└── figs/
    ├── figS5_RrsIQU_1to1_postfix.png 수정 후 Rrs I·Q·U 1:1 산포도 (회색=수정 전)
    └── figQ3_interface_transfer.png  (진단) 수정 전 결함 형상 — 계면 Q 전달비
```

## 6. 검증 절차 (수정 반영 여부 확인)

1. `patch/surface_twa_qfold.diff` 적용(또는 `patch/surface.c`로 교체, SHA-256 대조) 후 재빌드.
2. **단위 게이트**: `tests/test_twa_spin2_gate.c` 빌드·실행 → (T_I, T_Q2, T_U2)=(1,1,1)±0.005, FKC ON/OFF 동일. 결함 상태면 T_Q2≈0.57.
3. **무회귀 비트 게이트**: PSSA 스모크 CSV SHA `98e4b2433be300eb…` 비트 동일 + wind=0 red_clay fullgrid 비트 동일.
4. **OSOAA 대조**: §2 조건으로 pure/red_clay fullgrid 산출 후
   `python3 verify/verify_twa_qfold_fix.py --pure <pure.csv> --tsm <tsm.csv>` → `VERIFY: PASS` (게이트 상세는 `verify/EXPECTED_GATES.md`).
   OSOAA 재실행 없이 동봉 기준 CSV의 osoaa 열을 참값으로 쓰므로 OCRT 빌드만 있으면 된다. 인수 없이 실행하면 동봉 CSV 자체 점검(설치 무결성).

## 7. 반영 위치

- 샌드박스 런타임: `~/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C` — 패치·재빌드 완료(binary sha 70833d4d…). OSOAA 기준·비교 스크립트 포함 전체 환경 유효.
- 로컬(<LOCAL_ROOT>): `code/OCRT_C/src/shared/surface.c`(패치본), `code/OCRT_C/tests/test_twa_spin2_gate.c`, `code/OCRT_C/build/ocrt.exe`·`ocrt_x86-64-baseline.exe`(수정 트리로 재컴파일, -march=x86-64-v3/x86-64, static) 커밋됨.
- 검증 CSV 최신본(수정 후): `docs/validation/cross_check_fr631_20260823/`의 `Rrs_above_156geom.csv`·`rrs_underwater_raa90.csv`·`atm_rhoI_bfo.csv`는 2026-08-25 수정 후 버전으로 갱신됨(수정 전 Rrs는 `Rrs_above_156geom_prefix_20260824.csv`로 보존).
- 본 압축 패키지: `docs/technical/workorders_2026-07/OCRT_WORKORDER_TWA_QFOLD_FIX_2026-08-25.zip`.
