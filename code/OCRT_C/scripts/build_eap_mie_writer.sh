#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
CC=${CC:-gcc}
mkdir -p "$ROOT/build"
"$CC" -std=c11 -O3 -march=native -fopenmp -I"$ROOT/src" \
  -o "$ROOT/build/eap_mie_writer" \
  "$ROOT/tools/eap_mie_writer.c" \
  "$ROOT/src/internal/rt_eap_species_catalog.c" \
  "$ROOT/src/internal/rt_eap_coated_mie.c" -lm
printf '%s\n' "$ROOT/build/eap_mie_writer"
