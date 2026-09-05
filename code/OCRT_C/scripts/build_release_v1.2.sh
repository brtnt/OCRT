#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

CC_BIN=${CC:-gcc}
OCRT_MARCH=${OCRT_MARCH:-cascadelake}
OUT=${1:-build/ocrt_v1.2}
mkdir -p "$(dirname "$OUT")"

mapfile -t SOURCES < <(find src -name '*.c' -print | LC_ALL=C sort)

export LC_ALL=C
export TZ=UTC
export SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH:-1784332800}

"$CC_BIN" \
  -std=c11 -O3 -march="$OCRT_MARCH" \
  -ffp-contract=fast -fassociative-math \
  -fno-signed-zeros -fno-trapping-math \
  -DOCRT_FAST_KERNELS -fopenmp \
  -Isrc "${SOURCES[@]}" \
  -o "$OUT" -lm

chmod 0755 "$OUT"
for target in build/ocrt build/v2_solver_vk build/ocrt_v1.2; do
  if [[ "$OUT" != "$target" ]]; then cp -f "$OUT" "$target"; fi
done

printf 'built %s (march=%s)\n' "$OUT" "$OCRT_MARCH"
"$OUT" --version
sha256sum "$OUT"
