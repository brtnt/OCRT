#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
OUT=build/test_diffuse_top_primary_operator_closure
mkdir -p build
mapfile -t SOURCES < <(find src -name '*.c' ! -name 'main.c' -print | LC_ALL=C sort)

# Compile the production source with a test-only wrapper exposing the static
# diffuse-top injector. No production symbol is added in normal builds.
"$CC_BIN" -std=c11 -O2 -D_POSIX_C_SOURCE=200809L \
  -DOCRT_TEST_DTP_CLOSURE -fopenmp -Isrc \
  tests/test_diffuse_top_primary_operator_closure.c \
  "${SOURCES[@]}" -o "$OUT" -lm
"$OUT"
