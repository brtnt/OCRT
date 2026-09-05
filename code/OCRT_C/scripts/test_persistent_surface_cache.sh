#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN=${OCRT_BIN:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/cache" "$TMP/empty" "$TMP/out"
cd "$ROOT"
COMMON=(--surface ocean --wind-speed 3 --decouple-sunglint --n-water 1.34
 --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0.5 --ocrt-adom440 0
 --ocrt-tsm-species red_clay --sza 40 --wavelength 555 --pressure 0 --aod 0
 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
 --n-mu-water 16 --water-m-max 4
 --lut-vza-step 20 --lut-vza-max 60 --lut-raa-step 60)
ENV=(OCRT_ADVANCED=1 OMP_NUM_THREADS=1 OCRT_WATER_PARTICLE_KERNEL=direct OCRT_WATER_VALUE_NPHI=128)
env "${ENV[@]}" OCRT_SURFACE_PERSIST_CACHE_MODE=off \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/out/off.csv" >/dev/null 2>"$TMP/off.err"
env "${ENV[@]}" OCRT_SURFACE_PERSIST_CACHE_MODE=build \
  OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/cache" \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/out/build.csv" >/dev/null 2>"$TMP/build.err"
env "${ENV[@]}" OCRT_SURFACE_PERSIST_CACHE_MODE=required \
  OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/cache" \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/out/required.csv" >/dev/null 2>"$TMP/required.err"
cmp "$TMP/out/off.csv" "$TMP/out/build.csv"
cmp "$TMP/out/off.csv" "$TMP/out/required.csv"
[[ $(find "$TMP/cache" -name 'ocrt_surface_v*.bin' -type f | wc -l) -gt 0 ]]

set +e
env "${ENV[@]}" OCRT_SURFACE_PERSIST_CACHE_MODE=required \
  OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/empty" \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/out/missing.csv" >/dev/null 2>"$TMP/missing.err"
rc=$?
set -e
[[ $rc -ne 0 ]]
grep -F '[surface-pcache] miss:' "$TMP/missing.err" >/dev/null

cat > "$TMP/batch.csv" <<CSV
out,wavelength
$TMP/out/forbidden.csv,555
CSV
set +e
env "${ENV[@]}" OCRT_SURFACE_PERSIST_CACHE_MODE=build \
  OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/cache" \
  "$BIN" "${COMMON[@]}" --batch-full-grid "$TMP/batch.csv" >/dev/null 2>"$TMP/batch.err"
rc=$?
set -e
[[ $rc -ne 0 ]]
grep -F 'BUILD mode is forbidden inside --batch-full-grid' "$TMP/batch.err" >/dev/null

echo 'PASS: persistent surface cache off/build/required bit-identical; required miss fail-loud; batch build forbidden'
