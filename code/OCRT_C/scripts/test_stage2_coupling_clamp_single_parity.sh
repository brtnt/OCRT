#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REFROOT="$ROOT/validation/stage2_coupling_clamp_fix_2026-08-02_dtpsign"
OUTDIR="$REFROOT/single"
mkdir -p "$OUTDIR"
export OCRT_ADVANCED=1 OCRT_DEBUG=1 OMP_NUM_THREADS=1
"$ROOT/build/ocrt" \
  --sza 40 --wavelength 443 --surface ocean --wind-speed 3 --pressure 1013.25 \
  --n-mu 24 --n-layers 40 --m-max 2 --sos-max-orders 100 --debug-water-max-orders 100 \
  --n-water 1.34 --water-temperature 20 --water-salinity 38.4 \
  --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 --ocrt-adom-slope 0.014 \
  --water-m-max 4 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 --decouple-sunglint \
  --n-mu-water 64 --vza 30 --raa 45 > "$OUTDIR/single_n64_vza30_raa45.out"
python "$ROOT/scripts/check_stage2_coupling_clamp_single_parity.py" \
  --single "$OUTDIR/single_n64_vza30_raa45.out" \
  --full-grid "$REFROOT/reference/fullgrid_n64.csv" \
  | tee "$OUTDIR/SINGLE_FULLGRID_CONSISTENCY_RESULTS.txt"
