#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
OUT="$ROOT/build/test_aw_interp_workspace_sanitized"
mapfile -t SRC < <(find "$ROOT/src" -name '*.c' ! -name 'main.c' -print | LC_ALL=C sort)
${CC:-gcc} -std=c11 -O1 -g -D_POSIX_C_SOURCE=200809L -DOCRT_TEST_AW_INTERP \
  -fopenmp -fsanitize=address,undefined -fno-omit-frame-pointer -I"$ROOT/src" \
  "$ROOT/tests/test_aw_interp_workspace.c" "${SRC[@]}" -o "$OUT" -lm
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
  OMP_NUM_THREADS=1 "$OUT"
