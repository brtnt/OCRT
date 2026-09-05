#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
test -f "$ROOT/inputs/water_iop/candidates/Detritus_Stramski2001_330_1100_CANDIDATE_REJECTED_L200.mie"
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror -I"$ROOT/src" \
  "$ROOT/tests/test_detritus_candidate_gate.c" \
  "$ROOT/src/rt_iop_organic.c" "$ROOT/src/shared/mie_io.c" \
  "$ROOT/src/shared/numerics.c" "$ROOT/src/shared/io_utils.c" -lm \
  -o "$TMP/test_detritus_gate"
"$TMP/test_detritus_gate" "$ROOT/inputs/water_iop"
