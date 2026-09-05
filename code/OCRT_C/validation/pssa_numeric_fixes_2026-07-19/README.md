# PSSA numeric fixes validation assets

Date: 2026-07-19

## Scope

The matrices in this directory quantify four changes relative to the preceding
`pssa-aodref-radiometry-split-layerwarn` package:

1. reflected-beam lower-endpoint attenuation sign,
2. gas-only ocean PSSA direct-transmittance helper,
3. reflected downward first-order coupling into higher SOS orders,
4. layer-local reflected-beam `beta` phase source.

## Source sequence

```text
b92a476 baseline package
f07a0e7 endpoint sign
4a831e5 gas-only ocean helper
8a6a380 downward first-order coupling
ff46823 layer-local beta source
f284cb4 direct spherical reflected amplitude
27134f1 pure-Rayleigh beta setup optimization
```

Intermediate binaries were used to generate differential matrices and were
removed from the distribution. Their SHA-256 identifiers were:

```text
baseline        db3eba472d4fd63e94268d0a649c958181590ed734e36dbc1e48801ba6dc8a95
endpoint        7366ede52ce97e24930d3be62aa49b02f52f97e257aa2d383fb29dc084916ad2
gas-only        a15449178e15faf1e4ab2c5755384bc979860bc71c6971f8470304072d19a6eb
downward-order  2c3627a69d46bb4bcb44c5ec85df0d3291392ac09a1e289390855a281816695b
local-beta      6bf513bf28705302d44a65c32a1c4efe952788eecfc9808af597cae6a48b157a
stable-amplitude bc178a40488e819369922ae6721717a92ba2cb2e34a3700142c187b093564093
basis-opt       f8a332e3806bf27c5af9acbd0f11e47c05f71723ade51b57301580ec711f486e
rayleigh-fast   53d68444e1f62ccf36309be3669d2083ba0b698be5348de87510113c00bee933
```

## Matrices

```text
PSSA_NUMERIC_FIXES_RAYLEIGH_MATRIX_2026-07-19.csv       72 cases
PSSA_NUMERIC_FIXES_AEROSOL_MATRIX_2026-07-19.csv        12 cases
PSSA_NUMERIC_FIXES_PP_FLAT_MATRIX_2026-07-19.csv        27 cases
PSSA_NUMERIC_FIXES_GAS_ONLY_OCEAN_MATRIX_2026-07-19.csv 24 cases
PSSA_CORRECTION_BEFORE_AFTER_2026-07-19.csv               8 cases
```

JSON summaries retain maxima and representative-case metadata. For Q/U, use
absolute differences or the normalized `|delta Q,U|/|I|` values in
`PSSA_NUMERIC_FIXES_ROBUST_SUMMARY_2026-07-19.json`; raw relative percentages
are unstable when Q or U crosses zero.

## Reproducible current-code checks

```bash
./scripts/smoke_pssa_numeric_fixes.sh ./build/ocrt
./scripts/test_pssa_internal.sh
```

The historical differential matrices require the previous/intermediate
binaries identified above and are therefore audit records rather than a
standalone rerun harness.
