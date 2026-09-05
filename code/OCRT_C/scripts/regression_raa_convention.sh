#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
mkdir -p build validation/raa_convention

"$CC_BIN" -std=c11 -O2 -Isrc tests/test_raa_convention.c src/rt_fourier.c -o build/test_raa_convention -lm
./build/test_raa_convention

BIN=${OCRT_BIN:-./build/ocrt}
GAS_OFF=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
         --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
COMMON=(--surface black_fresnel_ocean --wind-speed 3 --sza 30 --vza 30 \
        --wavelength 865 --pressure 1013.25 --aod 0 --n-mu 24 --m-max 16 \
        "${GAS_OFF[@]}")

extract_glint() {
  local raa=$1
  OMP_NUM_THREADS=1 OCRT_ADVANCED=1 "$BIN" "${COMMON[@]}" --raa "$raa" \
    2>"validation/raa_convention/raa${raa}.stderr" \
    | tee "validation/raa_convention/raa${raa}.stdout" \
    | sed -n 's/.*TOA_rho_glint_direct_I=\([^ ]*\).*/\1/p'
}

g0=$(extract_glint 0)
g90=$(extract_glint 90)
g180=$(extract_glint 180)
python3 - "$g0" "$g90" "$g180" <<'PY'
import math, sys
g0,g90,g180=map(float,sys.argv[1:])
assert all(math.isfinite(x) and x >= 0 for x in (g0,g90,g180))
assert g180 > g90 > g0, (g0,g90,g180)
assert g180 > 1.0e6 * max(g0,1.0e-300), (g0,g180)
print(f"PASS: direct-glint branch RAA180; g0={g0:.12e} g90={g90:.12e} g180={g180:.12e}")
PY
