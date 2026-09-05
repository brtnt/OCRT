#!/usr/bin/env bash
set -euo pipefail
PKG_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
RUNTIME_DIR=${1:-$PKG_ROOT/runtime}
OUT_DIR=${2:-$PKG_ROOT/results/pure_rayleigh_5band_pssa}
python3 "$PKG_ROOT/runner/run_campaign.py" \
  --runtime-root "$RUNTIME_DIR/MIGRATION_PKG_2026-08-19" \
  --matrix "$PKG_ROOT/matrices/PURE_RAYLEIGH_15RUN_MATRIX.csv" \
  --output-dir "$OUT_DIR" --workers 16
