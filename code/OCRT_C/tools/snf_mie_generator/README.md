# SnF aerosol Mie generator

This directory contains the reproducible generator used for the 2026-07-19 SnF
atmospheric aerosol integration. It implements the OCRT V1/6SV-lineage
Bohren-Huffman Mie calculation, log-normal radius integration, component mixing
and OPAC/OCRT P11/P12/P33 writer.

## Canonical output

- wavelengths: 20 OPAC wavelengths, 0.350–3.750 µm;
- angles: 180° to 0° at 0.5° spacing, 361 points;
- radius step: `rlogpas=0.011`;
- output blocks: P11, P12 and P33;
- model set: T50/T80/T90/T95, C50/C70/C80/C90/C95,
  M50C/M70C/M80C/M90C/M95C/M98C, O99.

`O99` is the canonical model name supplied by the source archive; no `O99C`
input exists.

## Build and generate

Run from `ocrt/`:

```bash
make -C tools/snf_mie_generator
OMP_NUM_THREADS=32 ./scripts/generate_snf_aerosols.sh /tmp/snf_mie
```

The source `.inp` files are preserved under `inputs/aerosol_snf_inp/`.

## Runtime smoke

```bash
python3 tools/snf_mie_generator/run_snf_smoke.py \
  --exe ./build/ocrt \
  --mie-dir /tmp/snf_mie \
  --out /tmp/snf_smoke.csv \
  --workers 16
```

## Controlled regression result

The pre-existing high-resolution T50, C50, M80C and O99 caches were regenerated
at file precision. Spectral values, P11, P12 away from 0°/180° and P33 were
exact; the only differences were sub-2.3e-15 floating-point residuals at the
algebraic-zero P12 endpoints. The runtime package therefore keeps these four
pre-existing files byte-for-byte. The other twelve canonical caches are direct
generator outputs.

A complete validation archive is under
`validation/snf_aerosol_2026-07-19/`.
