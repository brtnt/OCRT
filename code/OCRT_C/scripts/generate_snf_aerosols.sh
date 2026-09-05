#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
OUT=${1:-"$ROOT/generated_snf_mie"}
GEN="$ROOT/tools/snf_mie_generator/snf_mie_gen"
INP="$ROOT/inputs/aerosol_snf_inp"
MODELS=(T50 T80 T90 T95 C50 C70 C80 C90 C95 M50C M70C M80C M90C M95C M98C O99)

mkdir -p "$OUT"
make -C "$ROOT/tools/snf_mie_generator" >/dev/null
for model in "${MODELS[@]}"; do
  "$GEN" "$INP/$model.inp" "$OUT/$model.mie" 0.5 0.011
done
count=$(find "$OUT" -maxdepth 1 -type f -name '*.mie' | wc -l)
printf 'generated %s SnF aerosol files in %s\n' "$count" "$OUT"
[[ "$count" -eq 16 ]]
