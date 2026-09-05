#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
OCRT_MARCH=${OCRT_MARCH:-cascadelake}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mapfile -t SOURCES < <(find src -name '*.c' -print | LC_ALL=C sort)
export LC_ALL=C TZ=UTC SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH:-1784332800}
for n in 1 2; do
  "$CC_BIN" -std=c11 -O3 -march="$OCRT_MARCH" \
    -ffp-contract=fast -fassociative-math \
    -fno-signed-zeros -fno-trapping-math \
    -DOCRT_FAST_KERNELS -fopenmp -Isrc "${SOURCES[@]}" \
    -o "$TMP/ocrt.$n" -lm
  sha256sum "$TMP/ocrt.$n"
done
h1=$(sha256sum "$TMP/ocrt.1" | awk '{print $1}')
h2=$(sha256sum "$TMP/ocrt.2" | awk '{print $1}')
[[ "$h1" == "$h2" ]]
"$TMP/ocrt.1" --version
printf 'PASS: reproducible release build march=%s sha256=%s\n' "$OCRT_MARCH" "$h1"
