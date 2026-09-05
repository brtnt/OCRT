# TSM mineral a*/b* Mie-graft 개정 반영 기록 — 2026-08-29

전달 패키지 `OCRT_TSM_astar_bstar_MIEgraft_2026-08-29_KST` (4광물 × {a*,b*} = 8파일)를 검증하고 OCRT 트리에 반영했다. 판정: **수용·반영 완료, 기존 정합성/속도 기준 무효화 없음**.

## 1. 개정 내용 (패키지 요지)

- 330–400 nm: 구 OLS 선형외삽 → **Mie graft 교체** (구 방식이 400 nm 부근 b* 극소를 놓치고 db*/dλ 부호를 반대로 만들었음).
- 750–1100 nm: `*_AHN.mie` bulk 컬럼 유래 값 → **Mie graft 교체**.
- 400–450 nm, 700–750 nm: **선형 blend 창 신설**.
- **450–700 nm 및 >1100 nm sparse 노드: bit-identical 불변.** 파장 grid·포맷·종단자 불변(drop-in).

## 2. 독립 검증 (샌드박스, 전 항목 통과)

| 검정 | 결과 |
|---|---|
| 패키지 자체 gate (`src/verify.py`) | **GATE PASS — 8/8** (grid 보존, 450–700 불변 0.0e+00, >1100 불변, 400/750 접붙임 연속, 양수, 690–760 단조) |
| SHA256SUMS 전수 | **34/34 OK** |
| 패키지 `data_original` vs 현재 설치본 | 8/8 **완전 일치** — 개정의 기준선이 우리 트리와 동일함을 확인 |
| 450–700 nm 재검산(독립 스크립트) | 8파일 전부 최대 상대차 **0.0** |
| red_clay 555 nm 상수 | a*=0.027810, b*=0.764760 — **정합성 검증 상수와 동일(불변)** |
| 밴드별 변화율 재검산 | VALIDATION §4 표와 일치 (412 최대 −2.66%, 865 최대 ±2.48%, 490–709 0.00%) |

## 3. OCRT 실행 검증 (결정적)

신 a*/b* 설치 후 실제 실행:

- **555 nm 정례 red_clay 케이스(312기하): 기준 산출과 비트 동일** → `task\ref\reference_redclay_case.csv`, 교차검증 게이트, 속도 벤치마크 **전부 유효 유지**.
- 공식 하니스 재실행: pure PASS +0.245% / red_clay −0.258·+0.659%(기존 κ-별건) / rayleigh_bfo PASS 0.162% — **수정 전과 자릿수까지 동일**.
- 865 nm 대조군: b_total **+2.021%**, Rrs_I +1.95~2.02% 변화 — 신 파일이 실제로 소비되고 있으며 VALIDATION 표(red_clay b* 865 nm 2.02%)와 일치.

## 4. a*/b* 소비 경로 확인 (패키지 경고 #2의 실제 영향 범위)

`src/rt_iop_ahn_mineral.c` 확인 결과 **`astarmin_*.txt`/`bstarmin_*.txt`가 production 경로의 유일한 a*/b* 출처**이고, `.mie` bulk 컬럼(Extinct_Co/Scatter_Co)을 읽는 `load_mie_star()`는 **.txt 파일이 없을 때만 쓰이는 fallback**이다. `.mie` phase 블록(P11/P12/P33)은 위상함수 전용이며 a*/b*와 무관하다.

→ 따라서 패키지가 지적한 **`.mie` re-mirroring 미완은 production 수치에 영향이 없다**. 다만 fallback 경로·외부 소비자에게 구/신 불일치가 남으므로 후속 과제로 유지한다(§5-a).

## 5. 미해결·후속 (패키지 지시 포함)

- (a) **`.mie` bulk 컬럼 re-mirroring**: 미수행. mirror 범위 선정 필요(`.mie` phase grid는 350 nm 시작, 개정은 330 nm부터). production 영향 없음(§4).
- (b) **330–399 nm 채움 방식 최종 결정**: 구 OLS / `.mie` 자체값 / Mie graft 세 후보 중 미확정. 본 반영은 패키지 채택안(Mie graft)을 따랐고 재검토 가능.
- (c) **R_eff는 물리 입자크기가 아니다** — 분광기울기 모형의 유효 파라미터. **phase function·P11/P12/P33·bb/b 계산에 사용 금지**(METHOD §6).
- (d) **문서 정정 3건**(CHANGELOG): "H5.0 spline + linear-to-zero tail" 서술 및 `a* zero@2000 nm`/`b* zero@1926 nm` 류 수치 **폐기 대상**; miepython 재현 실패는 `.inp % density`를 number fraction으로 오독한 것으로 원인 규명됨; Junge ξ=4는 이 4종 mineral에 부적합. → 논문·기술문서에서 해당 서술 수정 필요.
- (e) bb/b 상수(Yellow 0.0092 / Brown 0.0099 / Calcareous 0.0122 / Red 0.0067) 미변경.
- (f) 700 nm 기울기 불연속은 **구파일 자체 특성**(개정이 만든 것 아님). 원본 생성 이력 확인은 별도 항목.

## 6. 반영 위치

- 로컬: `code\OCRT_C\inputs\tsm_ahn\` 8파일, `code\OCRT_Python\data\tsm_ahn\` 8파일 — 양 트리 SHA-256 8/8 일치 확인.
- 원본 패키지: `package\OCRT_TSM_astar_bstar_MIEgraft_2026-08-29_KST.tar.gz` (버전관리 보관).
- 샌드박스 런타임(C·Python 양쪽) 동일 반영, 위 §3 검증 수행.
- 구파일은 패키지 `data_original\`에 보존되어 있으므로 별도 백업 불필요.
