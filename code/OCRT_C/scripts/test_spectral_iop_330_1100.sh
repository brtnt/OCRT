#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror -I"$ROOT/src" \
  "$ROOT/tests/test_spectral_iop_330_1100.c" \
  "$ROOT/src/rt_water_iop.c" "$ROOT/src/rt_iop_organic.c" \
  "$ROOT/src/rt_iop_ahn_mineral.c" "$ROOT/src/shared/mie_io.c" \
  "$ROOT/src/shared/numerics.c" "$ROOT/src/shared/io_utils.c" -lm \
  -o "$TMP/test_spectral_iop"
"$TMP/test_spectral_iop" "$ROOT/inputs/water_iop" "$ROOT/inputs/tsm_ahn"
