#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

CC_BIN=${CC:-gcc}
OUT=build/test_stage2_water_raa_fix
mkdir -p build

"$CC_BIN" -std=c11 -O2 -Isrc   tests/test_raa_convention.c src/rt_fourier.c   -o "$OUT" -lm

"$OUT"

if grep -RIn --include='*.c'     'rt_raa_to_water_view_phi' src tests >/tmp/ocrt_water_raa_legacy_hits.txt; then
  cat /tmp/ocrt_water_raa_legacy_hits.txt >&2
  echo "FAIL: legacy water-view RAA helper remains in source" >&2
  exit 1
fi

echo "PASS: stage-2 water RAA source audit"
