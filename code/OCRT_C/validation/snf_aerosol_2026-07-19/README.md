# SnF aerosol validation archive — 2026-07-19

This directory records the controlled conversion and integration of the 16 SnF
atmospheric aerosol inputs.

## Pass summary

| Check | Result |
|---|---:|
| Baseline regeneration: T50/C50/M80C/O99 | PASS |
| Baseline RT identity | 24/24 exact |
| Deterministic regeneration | 16/16 byte-identical |
| Structural/finite/SSA checks | 16/16 pass |
| Mueller physicality | 16/16 pass |
| OCRT aerosol-only runtime smoke | 48/48 pass |

## Files

- `baseline_reproduction.csv`: direct numeric comparison of regenerated and
  pre-existing high-resolution baseline caches.
- `baseline_rt_identity_24cases.csv`: I/Q/U and AOD equality for the baseline
  caches in two pressure regimes and three wavelengths.
- `generated_quality.csv`: structural, normalization, g-consistency and
  physicality metrics for all 16 generated models.
- `legacy_replacement.csv`: phase and spectral differences between displaced
  83-angle M50C/M95C/M98C and the canonical regenerated lineage.
- `legacy83_vs_canonical_rt_18cases.csv`: representative I/Q/U impact of the
  three lineage replacements.
- `runtime_smoke_48cases.csv`: 16 models × 412/550/860 nm, aerosol-only black
  surface, AOD555=0.15.
- `model_integration_matrix.csv`: final runtime action and qualification status.
- `generator_reproducibility.txt`, `generated_file_sha256.txt`: deterministic
  generation record.
- `validation.json`: detailed machine-readable phase-file checks.

## Interpretation boundary

These files establish generator lineage, cache integrity, physical sanity and
OCRT runtime compatibility. They do not promote the newly added models into the
external OSOAA/6SV benchmark subset. That subset remains T50/C50/M80C pending a
separate full reference-code campaign.
- `non_pssa_regressions.log`: rerun of CCRR Chl, OCRT organic Chl, Ahn TSM
  and TSM phase-cache regressions after the SnF integration.
