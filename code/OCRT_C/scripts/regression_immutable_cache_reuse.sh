#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
BIN=${OCRT_BIN:-./build/ocrt}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
COMMON=(--surface ocean --wind-speed 3 --decouple-sunglint \
 --sza 40 --wavelength 555 --pressure 1013.25 --aod 0 \
 --water-model ocrt --ocrt-chl 0 --ocrt-tsm 50 --ocrt-adom440 0 \
 --ocrt-tsm-species calcareous_sand \
 --n-mu 48 --n-layers 26 --sos-max-orders 100 --sos-tolerance 1e-3 \
 --n-mu-water 48 --debug-water-tau-max-target 30 \
 --debug-water-n-layers 80 --debug-water-max-orders 700 \
 --lut-vza-step 30 --lut-vza-max 60 --lut-raa-step 90 "${GAS[@]}")
run_case() {
  local name=$1; shift
  env OCRT_ADVANCED=1 OCRT_DEBUG=1 OMP_NUM_THREADS=1 "$@" \
    "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/$name.csv" \
    >"$TMP/$name.out" 2>"$TMP/$name.err"
}
run_case on
run_case surface_off OCRT_AIR_SURFACE_CACHE_OFF=1
run_case transport_off OCRT_TRANSPORT_CACHE_OFF=1
run_case all_off OCRT_AIR_SURFACE_CACHE_OFF=1 OCRT_TRANSPORT_CACHE_OFF=1
cmp "$TMP/on.csv" "$TMP/surface_off.csv"
cmp "$TMP/on.csv" "$TMP/transport_off.csv"
cmp "$TMP/on.csv" "$TMP/all_off.csv"
grep -F '[native-coupled-lut] calls=3 water_cold=1 water_views=3' \
  "$TMP/on.err" >/dev/null
python3 - "$TMP/on.csv" <<'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1])))
assert len(rows)==12, len(rows)
assert max(int(r['water_orders']) for r in rows) > 1
assert all(int(r['water_converged']) == 1 for r in rows)
print('PASS immutable cache reuse: 12-cell CSV byte-identical; native water solve once')
PY
