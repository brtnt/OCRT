# OCRT–OSOAA Stage-2 EAP update validation report

Date: 2026-07-26

## 1. Scope and provenance

The handoff defines an OCRT target state and supplies internal OCRT–OSOAA
underwater comparison data. It does **not** include the modified OSOAA
executable or a complete cumulative OSOAA patch that produced those CSVs.
Accordingly:

- OCRT C was updated to the handoff's Stage-2 contract while preserving later
  approved OCRT features already present in this development session.
- The paired Python implementation was updated to the same production
  constituent contract.
- The public OSOAA `src/*.F` core was not changed by inference. A portable
  build/smoke/reference-audit/performance harness was added around the
  unchanged public core.
- Exact raw re-generation of the supplied OSOAA columns is not claimed.

## 2. OCRT C changes

Version:

`OCRT-v1.2-2026-07-26-KST-stage2-eap-phyto-scattering-disabled-integrated`

### 2.1 EAP generator/catalog

- Retains the 17-species coated-sphere P11/P12/P33 API and generated data.
- Preserves phase normalization, Mueller physicality and deterministic output
  gates.
- Keeps the optimized writer path: spectral-only calculations skip angular
  matrix work; independent radius bins may be calculated with OpenMP and are
  accumulated in deterministic radius order.

### 2.2 Production constituent Chl contract

- `Chl > 0` without explicit species: phytoplankton absorption is retained.
- The validated default absorption is read from
  `inputs/water_iop/phyto_absorption_default.csv`.
- Phytoplankton `b=0` and `bb=0` exactly.
- Chl-linked detritus scattering remains active.
- `Chl > 0` with an explicit `--ocrt-phyto-group`: fail-loud, exit code 2.
- Unsupported constituent truncation controls fail loudly rather than being
  ignored.
- Fixed-bulk IOP phase/truncation remains available.

The absorption spectrum was deliberately separated from the EAP `.mie`
catalog. This prevents regenerated species phase files from silently changing
Stage-2 production Chl absorption.

### 2.3 Preserved later OCRT functionality

The integration preserves the existing:

- `Rrs(0+)` and `rrs(0-)` I/Q/U output,
- water internal-reflection odd-m sign correction,
- public water RAA/U convention,
- air–water FIX1–FIX4 and exact-pole treatment,
- water-to-air full-grid coupling clamp fix,
- `black_fresnel_ocean` surface naming,
- full-grid aerosol-object lifecycle fix,
- LUT exact-cache performance improvements,
- PCHIP P11/P12/P33 wavelength interpolation.

## 3. Python paired update

Version:

`pyOCRT-v1.2-2026-07-26-stage2-eap-phyto-scattering-disabled`

- Uses the byte-identical C/Python
  `phyto_absorption_default.csv` file.
- Returns phytoplankton `b=bb=0` in scalar and batch IOP assembly.
- Retains detritus/mineral scattering.
- Rejects explicit EAP species selection for `Chl > 0` in single, full-grid and
  12,000-case production entry points.
- Skips phytoplankton phase-moment preparation in the production batch path;
  no extra RT call was introduced.
- Retains the 17-species phase catalog for generator/future validation use.

## 4. OSOAA public package

The packaged OSOAA core is the handoff's public reference source, unchanged.
Added `stage2_harness/` includes:

- source-hash verification,
- public build and demo smoke scripts,
- supplied reference-data audit,
- warm runtime benchmark,
- Korean scope/provenance documentation.

The included Linux executable was built from that unchanged public core.

## 5. Validation results

### 5.1 C physics and interface gates

| Gate | Result |
|---|---|
| 17-species EAP API normalization/physicality/determinism | PASS |
| Stage-2 Chl absorption-only, phyto b/bb zero, detritus phase | PASS |
| Explicit species fail-loud | PASS |
| Unsupported constituent truncation fail-loud | PASS |
| Water RAA convention | PASS |
| FIX1+FIX2+FIX3 interface contraction | PASS |
| FIX4 diffuse-top U-column closure | PASS |
| Exact-pole surface limit | PASS |
| Water internal-reflection m-sign | PASS |
| Coupling full-grid n_mu_water 48/64/96 | PASS |
| Single/full-grid parity | PASS |
| Rrs/rrs I/Q/U and IOP/Kd CSV | PASS |
| Mie phase-wavelength PCHIP/source audit | PASS |
| Frozen-module source audit | PASS |

### 5.2 Baseline and expected-final parity

- Non-Chl baseline regression: **18/18** conditions are byte-identical in both
  stdout and stderr.
- Chl expected-final parity: **39/39** recorded values have absolute difference
  0.0. The compared fields cover TOA I/Q/U, Rrs I/Q/U, rrs I/Q/U, Kd,
  `a_chl`, particle `b` and `bb` at 443/555/660 nm.
- Supplied underwater validation/summary audit:
  - validation shape: **11,154 × 37**,
  - summary shape: **50 × 26**,
  - recomputed summary maximum absolute difference: **0.0**.

### 5.3 C–Python IOP parity

At Chl=0.3 mg m-3 and 443/490/555/660 nm:

- maximum C–Python absolute difference in `a_phyto`:
  **1.735e-18 m-1**,
- Python `b_phyto=bb_phyto=0` at all four wavelengths,
- C and Python absorption CSV SHA-256:
  `c257981a8835858e7c059f36cc478d995b9fcfd3be8be637de3400bd5bf9f02e`.

### 5.4 Python tests

- `pytest`: **17/17 PASS**.
- Includes default absorption value, scalar/batch zero scattering, explicit
  species rejection and existing RAA/m-sign/Rrs-output regressions.

### 5.5 OSOAA public-core checks

- All 12 audited `src/*.F` SHA-256 values: PASS after build and smoke.
- Public build: PASS.
- Public demo smoke: PASS.
- Reference audit: PASS, max summary difference 0.0.
- Three-run warm benchmark: median **1.03958 s**.

This is a public-core operational check, not a raw reproduction of the supplied
modified-OSOAA validation values.

## 6. Runtime regression

### 6.1 OCRT C

Nine alternating old/new single-thread runs:

- old mean: 0.232458 s,
- new mean: 0.233394 s,
- mean change: **+0.403%**,
- old median: 0.233399 s,
- new median: 0.234207 s,
- median change: **+0.346%**.

This is below the handoff's 5% profiling threshold.

### 6.2 Python phase preparation

Three alternating four-band measurements:

- prior median: 4.47261 s,
- Stage-2 median: 2.52363 s,
- change: **-43.576%**.

The gain comes from not generating phytoplankton phase moments when its
production scattering coefficient is exactly zero.

### 6.3 EAP writer

The optimized writer reproduced the supplied sample files byte-for-byte. The
large 24-µm representative species can use OpenMP radius-bin parallelism while
retaining deterministic serial accumulation.

## 7. Remaining limitation

The modified OSOAA binary/full cumulative patch used to create the supplied
Stage-2 comparison CSVs was not included in the handoff. A future exact OSOAA
code update requires that provenance. Until then, the supplied CSVs are treated
as immutable reference data and the packaged public OSOAA core remains
unchanged.
