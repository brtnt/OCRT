#!/usr/bin/env bash
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/.." && pwd)
PREBUILT="$ROOT/build/test_ext_bottom_mode_bounds"
if [[ -x "$PREBUILT" && "${OCRT_REBUILD_EXT_TEST:-0}" != 1 ]]; then
  "$PREBUILT"
  exit 0
fi
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
OCRT_MARCH=${OCRT_MARCH:-cascadelake}
mapfile -t TEST_SRC < <(find "$ROOT/src" -name '*.c' ! -name 'main.c' | LC_ALL=C sort)
${CC:-gcc} -std=c11 -O3 -march="$OCRT_MARCH" \
  -ffp-contract=fast -fassociative-math -fno-signed-zeros -fno-trapping-math \
  -DOCRT_FAST_KERNELS -fopenmp -I"$ROOT/src" \
  "$ROOT/tests/test_ext_bottom_mode_bounds.c" "${TEST_SRC[@]}" \
  -o "$TMP/test_ext_bottom_mode_bounds" -lm
"$TMP/test_ext_bottom_mode_bounds"
