#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
BIN=${OCRT_BIN:-./build/ocrt}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
COMMON=(--surface ocean --water-model ocrt --ocrt-chl 1 --ocrt-tsm 0 \
 --ocrt-adom440 0 --wind-speed 3 --sza 30 --vza 30 --raa 90 \
 --pressure 0 "${GAS[@]}")
OCRT_ADVANCED=1 OCRT_DUMP_IOP=1 "$BIN" "${COMMON[@]}" --wavelength 850 \
  >"$TMP/850.out" 2>"$TMP/850.err"
OCRT_ADVANCED=1 OCRT_DUMP_IOP=1 "$BIN" "${COMMON[@]}" --wavelength 865 \
  >"$TMP/865.out" 2>"$TMP/865.err"
OCRT_ADVANCED=1 OCRT_DUMP_IOP=1 "$BIN" "${COMMON[@]}" --wavelength 1100 \
  >"$TMP/1100.out" 2>"$TMP/1100.err"
python3 - "$TMP/850.out" "$TMP/865.out" "$TMP/1100.out" "$TMP/865.err" "$TMP/1100.err" <<'PY'
import re,sys

def v(path,k):
 s=open(path).read();m=re.search(r'\b'+k+r'=([-+0-9.eE]+)',s)
 if not m: raise SystemExit('missing '+k)
 return float(m.group(1))
assert v(sys.argv[1],'a_chl') > 0.0
assert v(sys.argv[2],'a_chl') == 0.0
assert v(sys.argv[3],'a_chl') == 0.0
for path in sys.argv[4:]:
 s=open(path).read()
 assert 'setting Chl absorption to 0' not in s
print('PASS: 850-nm anchor retained; 865/1100-nm Chl absorption uses explicit zero extension')
PY
