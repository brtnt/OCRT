#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
OUT=$(mktemp); trap 'rm -f "$OUT"' EXIT
mapfile -t SOURCES < <(find src -name '*.c' ! -name main.c -print | LC_ALL=C sort)
"$CC_BIN" -std=c11 -O2 -fopenmp -Isrc tests/test_sos_zero_col_skip.c "${SOURCES[@]}" -o "$OUT" -lm
"$OUT"
