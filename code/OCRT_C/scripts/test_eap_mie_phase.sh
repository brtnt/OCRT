#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-"$ROOT/build/test_eap_mie_phase"}
mkdir -p "$(dirname "$OUT")"
gcc -std=c11 -O2 -march=cascadelake -I"$ROOT/src" \
  "$ROOT/tests/test_eap_mie_phase.c" \
  "$ROOT/src/rt_eap_mie_phase.c" \
  "$ROOT/src/internal/rt_eap_species_catalog.c" \
  "$ROOT/src/internal/rt_eap_coated_mie.c" -lm -o "$OUT"
"$OUT"
