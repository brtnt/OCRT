#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

"$CC_BIN" -std=c11 -O2 -Isrc \
  tests/test_angular_basis_single.c src/rt_angular_basis.c \
  -o "$TMP/test_angular_basis_single" -lm
"$TMP/test_angular_basis_single"

mapfile -t SRCS < <(find src -name '*.c' ! -name main.c -print | LC_ALL=C sort)
"$CC_BIN" -std=c11 -O2 -fopenmp -Isrc \
  tests/test_pssa_geometry.c "${SRCS[@]}" \
  -o "$TMP/test_pssa_geometry" -lm
OCRT_ADVANCED=1 OMP_NUM_THREADS=1 "$TMP/test_pssa_geometry"
