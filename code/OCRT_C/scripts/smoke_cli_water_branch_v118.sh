#!/usr/bin/env bash
# OCRT v1.18 water-input branch contract smoke test.
set -euo pipefail

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/.." && pwd)
BIN=${1:-$ROOT/build/v2_solver_vk_v1.18}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

[[ -x "$BIN" ]] || { echo "missing executable: $BIN" >&2; exit 2; }

pass=0
fail=0
ok(){ printf 'PASS %s\n' "$1"; pass=$((pass+1)); }
bad(){ printf 'FAIL %s\n' "$1"; fail=$((fail+1)); }

GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
GEOM=(--surface ocean --wind-speed 0 --sza 30 --vza 20 --raa 90
      --wavelength 443 --pressure 0 "${GAS[@]}")
FAST=(--debug-water-max-orders 4 --n-mu-water 12 --water-m-max 4)

run_fast(){
  OCRT_DEBUG=1 OCRT_ADVANCED=1 "$BIN" "$@" "${FAST[@]}"
}

expect_reject(){
  local name=$1 pattern=$2
  shift 2
  set +e
  "$BIN" "$@" >"$TMP/$name.out" 2>"$TMP/$name.err"
  local rc=$?
  set -e
  if [[ $rc -eq 2 ]] && grep -Eq -- "$pattern" "$TMP/$name.err"; then
    ok "$name"
  else
    bad "$name"
    printf 'rc=%s expected=%s\n' "$rc" "$pattern" >&2
    cat "$TMP/$name.err" >&2
  fi
}

expect_reject no_model 'requires exactly one --water-model' "${GEOM[@]}"
expect_reject duplicate_model 'must be specified exactly once' \
  --water-model ocrt --water-model ccrr \
  --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 "${GEOM[@]}"
expect_reject cross_prefix 'another water branch' \
  --water-model ocrt --ccrr-chl 0 --ccrr-tsm 0 --ccrr-adom440 0 "${GEOM[@]}"
expect_reject missing_required 'requires explicit --ocrt-chl, --ocrt-tsm and --ocrt-adom440' \
  --water-model ocrt --ocrt-chl 0 --ocrt-adom440 0 "${GEOM[@]}"
expect_reject legacy_unprefixed 'legacy/unprefixed water option' \
  --water-model ocrt --chl 0 "${GEOM[@]}"
expect_reject non_ocean_model 'require --surface ocean' \
  --surface black --water-model ocrt \
  --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 \
  --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0
expect_reject species_without_tsm 'requires --ocrt-tsm > 0' \
  --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 \
  --ocrt-tsm-species brown_earth "${GEOM[@]}"
expect_reject iop_incomplete 'requires explicit --iop-a, --iop-b and --iop-bb' \
  --water-model iop --iop-a 0.1 --iop-b 0.2 "${GEOM[@]}"
expect_reject iop_phase_conflict 'select at most one IOP phase source' \
  --water-model iop --iop-a 0.1 --iop-b 0.2 --iop-bb 0.01 \
  --iop-phase-lut dummy.csv --iop-mie-phase dummy.mie "${GEOM[@]}"
expect_reject orphan_iop_mie_control 'require --iop-mie-phase' \
  --water-model iop --iop-a 0.1 --iop-b 0.2 --iop-bb 0.01 \
  --iop-mie-truncation "${GEOM[@]}"

if run_fast --water-model ocrt \
     --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 "${GEOM[@]}" \
     >"$TMP/ocrt_zero.out" 2>"$TMP/ocrt_zero.err"; then
  ok ocrt_zero_pure
else bad ocrt_zero_pure; fi

if run_fast --water-model ccrr \
     --ccrr-chl 0 --ccrr-tsm 0 --ccrr-adom440 0 "${GEOM[@]}" \
     >"$TMP/ccrr_zero.out" 2>"$TMP/ccrr_zero.err"; then
  ok ccrr_zero_pure
else bad ccrr_zero_pure; fi

if cmp -s "$TMP/ocrt_zero.out" "$TMP/ccrr_zero.out"; then
  ok zero_branches_stdout_bit_identity
else bad zero_branches_stdout_bit_identity; fi

if run_fast --water-model iop \
     --iop-a 0.1762 --iop-b 2.618309 --iop-bb 0.038578 "${GEOM[@]}" \
     >"$TMP/iop.out" 2>"$TMP/iop.err"; then
  ok iop_branch
else bad iop_branch; fi

printf 'WATER_BRANCH_V118_SMOKE pass=%d fail=%d\n' "$pass" "$fail"
[[ $fail -eq 0 ]]
