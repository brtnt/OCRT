#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
OUT=build/test_surface_pole_limit
mkdir -p build
"$CC_BIN" -std=c11 -O2 -Isrc \
  tests/test_surface_pole_limit.c \
  src/shared/surface.c src/shared/surface_multibounce.c \
  -o "$OUT" -lm
"$OUT"

# Four physical rough-interface kernels must carry the explicit pole branch.
grep -Fq '#define OCRT_SURF_POLE_EPS 1.0e-8' src/shared/surface.c
count=$(grep -c '<= OCRT_SURF_POLE_EPS' src/shared/surface.c)
[[ "$count" -eq 8 ]]

echo 'PASS: stage-2 surface pole-limit source audit'
