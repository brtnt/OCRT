# AHN mineral phase contract — 2026-08-16

The active four-mineral phase files keep the validated July-2026 AHN angular phase matrices and
adopt the fixed OCRT FR631 phase-angle grid.

- The historical generator / `.inp` recipe was not in the handoff; these files are therefore a
  **lossless theta-linear remap**, not a new Mie calculation.
- All original 0.5° phase nodes at the common 350–860 nm phase wavelengths are kept exactly.
- The phase-angle grid is FR631 (`0–0.2°: 0.005°`, `0.2–1°: 0.02°`, `1–5°: 0.05°`,
  `5–20°: 0.1°`, `20–180°: 0.5°`).
- The 330 nm phase holds the 350 nm end point, because no validated historical generator input
  below 350 nm is available.
- The 1100 nm phase is a common-weight linear interpolation between the historical 860 and
  1240 nm matrices.
- The 330–1100 nm scalar `a*` and `b*` tables remain the accepted spectral-extension tables.
- The superseded source-informed reconstruction is kept under
  `validation/mie_candidates/tsm_source_informed_330_1100_rejected_20260816/` and is not active.

Reproduction tool: `tools/remap_legacy_tsm_ahn_fr631.py`.

---

# AHN 광물 위상 계약 — 2026-08-16

활성 광물 4종 위상 파일은 검증된 2026년 7월 AHN 각도 위상행렬을 유지하면서 OCRT 고정 FR631 위상각
격자를 채택한다.

- 과거 생성기 / `.inp` 레시피가 인계에 없었으므로 이 파일들은 새 Mie 계산이 아니라 **손실 없는 theta-linear
  재배치**다.
- 공통 위상 파장 350–860 nm 의 원래 0.5° 위상 절점은 모두 정확히 유지한다.
- 위상각 격자는 FR631 이다(`0–0.2°: 0.005°`, `0.2–1°: 0.02°`, `1–5°: 0.05°`, `5–20°: 0.1°`,
  `20–180°: 0.5°`).
- 330 nm 위상은 350 nm 끝값을 유지한다. 350 nm 아래의 검증된 과거 생성기 입력이 없기 때문이다.
- 1100 nm 위상은 과거 860 nm 와 1240 nm 행렬 사이의 공통 가중 선형 보간이다.
- 330–1100 nm 스칼라 `a*`·`b*` 표는 승인된 분광 확장표로 유지한다.
- 대체된 원자료 기반 재구성은 `validation/mie_candidates/tsm_source_informed_330_1100_rejected_20260816/`
  에 보존하며 활성이 아니다.

재현 도구: `tools/remap_legacy_tsm_ahn_fr631.py`.
