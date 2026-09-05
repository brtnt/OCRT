#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
CANDIDATE=${1:-build/ocrt}
BASELINE=${2:-build/ocrt_baseline}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

"$CC_BIN" -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror -Isrc \
  tests/test_molecular_profile_us62.c \
  src/rt_molecular_profile.c src/rt_rayleigh.c -lm -o "$TMP/test_us62"
"$TMP/test_us62"

[[ -x "$CANDIDATE" ]] || { echo "missing candidate: $CANDIDATE" >&2; exit 2; }
[[ -x "$BASELINE" ]] || { echo "missing baseline: $BASELINE" >&2; exit 2; }

gas_off=(
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
)

run_exact() {
  local name=$1; shift
  "$BASELINE" "$@" >"$TMP/$name.base.out" 2>"$TMP/$name.base.err"
  "$CANDIDATE" "$@" >"$TMP/$name.cand.out" 2>"$TMP/$name.cand.err"
  cmp "$TMP/$name.base.out" "$TMP/$name.cand.out"
  cmp "$TMP/$name.base.err" "$TMP/$name.cand.err"
  echo "PASS $name byte-exact"
}

run_exact rayleigh_black_gas_off \
  --surface black --wind-speed 0 --sza 40 --vza 30 --raa 90 \
  --wavelength 443 --pressure 1013.25 "${gas_off[@]}"
run_exact rayleigh_surface_gas_off \
  --surface coxmunk --wind-speed 3 --sza 60 --vza 50 --raa 120 \
  --wavelength 555 --pressure 1013.25 "${gas_off[@]}"

set +e
"$CANDIDATE" --surface black --wind-speed 0 --sza 40 --vza 30 --raa 90 \
  --wavelength 443 --pressure 1013.25 --mie inputs/M80C.mie --aod 0.3 \
  --aer-h-km 0 >"$TMP/aer0.out" 2>"$TMP/aer0.err"
rc=$?
set -e
[[ $rc -eq 2 ]]
grep -q -- '--aer-h-km must be positive' "$TMP/aer0.err"
echo 'PASS legacy non-positive aerosol scale height rejected'

echo 'ATMOSPHERE_PROFILE_US62 PASS'
