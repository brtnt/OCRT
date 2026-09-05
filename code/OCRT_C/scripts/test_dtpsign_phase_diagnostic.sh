#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export TERM=xterm

# Source-placement gate: conversion belongs at the coupled hand-off.
grep -F 'DTPSIGN (2026-08-01' "$ROOT/src/rt_solver.c" >/dev/null
grep -F 'ocrt_dt_I[z_] = -ocrt_dt_I[z_];' "$ROOT/src/rt_solver.c" >/dev/null

# Existing source/operator closure must remain exact.
"$ROOT/scripts/test_stage2_diffuse_top_uclosure.sh" "$BIN" >/dev/null

# Single TSM component must produce structured phase diagnostics.
OCRT_ADVANCED=1 OCRT_DEBUG=1 OCRT_DUMP_IOP=1 OMP_NUM_THREADS=1 "$BIN" \
  --surface ocean --water-model ocrt \
  --ocrt-chl 0 --ocrt-tsm 5 --ocrt-adom440 0 --ocrt-tsm-species red_clay \
  --wind-speed 3 --sza 30 --vza 30 --raa 90 --wavelength 555 \
  --pressure 0 --aod 0 \
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 \
  --n-mu-water 12 --water-m-max 4 \
  --debug-water-n-layers 60 --debug-water-max-orders 20 \
  --debug-water-tau-max-target 6 \
  >"$TMP/out" 2>"$TMP/err"
grep -F 'OCRT_PHASE_COMPONENT name=tsm weight_b=1 ' "$TMP/err" >/dev/null
grep -F 'OCRT_PHASE_MIX components=1 ' "$TMP/err" >/dev/null
grep -F 'Rrs0plus_I=' "$TMP/out" >/dev/null

echo 'PASS: DTPSIGN hand-off placement, closure, and single-component diagnostic'
