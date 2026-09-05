#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN=${1:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
GAS_OFF=(
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
)
COMMON=(--surface black --vza 0 --raa 0 --wavelength 412 "${GAS_OFF[@]}")

cd "$ROOT"
"$BIN" "${COMMON[@]}" --sza 85 --pssa \
  >"$TMP/sza85_n40.out" 2>"$TMP/sza85_n40.err"
OCRT_ADVANCED=1 "$BIN" "${COMMON[@]}" --sza 85 --pssa --n-layers 400 \
  >"$TMP/sza85_n400.out" 2>"$TMP/sza85_n400.err"
"$BIN" "${COMMON[@]}" --sza 85 \
  >"$TMP/sza85_pp.out" 2>"$TMP/sza85_pp.err"
"$BIN" "${COMMON[@]}" --sza 80 --pssa \
  >"$TMP/sza80_n40.out" 2>"$TMP/sza80_n40.err"

# Numerical output must not be modified by the advisory.
# The n=40 and n=400 values intentionally differ; only warning presence is tested.
grep -q '^warning: PSSA vertical resolution may be insufficient' "$TMP/sza85_n40.err"
grep -q '400-800 layers' "$TMP/sza85_n40.err"
if grep -q '^warning: PSSA vertical resolution may be insufficient' "$TMP/sza85_n400.err"; then
  echo 'unexpected warning for SZA=85, n_layers=400' >&2
  exit 1
fi
if grep -q '^warning: PSSA vertical resolution may be insufficient' "$TMP/sza85_pp.err"; then
  echo 'unexpected PSSA warning when --pssa is off' >&2
  exit 1
fi
grep -q '^warning: PSSA vertical resolution may be insufficient' "$TMP/sza80_n40.err"
grep -q 'use --n-layers >=200' "$TMP/sza80_n40.err"

printf 'PSSA_LAYER_WARNING PASS\n'
