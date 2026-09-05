# AHN TSM phase contract — 2026-08-16

The active four-mineral phase files preserve the validated July-2026 AHN angular phase matrices while adopting the fixed OCRT FR631 phase-angle grid.

- The historical generator/`.inp` recipe was not present in the handoff; these files are therefore a **lossless theta-linear remap**, not a new Mie regeneration.
- All original 0.5° phase nodes at the common 350–860 nm phase wavelengths are preserved exactly.
- The phase-angle grid is FR631 (`0–0.2°:0.005°`, `0.2–1°:0.02°`, `1–5°:0.05°`, `5–20°:0.1°`, `20–180°:0.5°`).
- The 330 nm phase uses a conservative 350 nm endpoint hold because no validated historical sub-350-nm generator input is available.
- The 1100 nm phase is common-weight linear interpolation between the historical 860 and 1240 nm matrices.
- The 330–1100 nm scalar `a*` and `b*` tables remain the accepted spectral-extension tables.
- The superseded source-informed reconstruction is retained under `validation/mie_candidates/tsm_source_informed_330_1100_rejected_20260816/` and is not active.

Reproduction tool: `tools/remap_legacy_tsm_ahn_fr631.py`.
