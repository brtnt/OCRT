#!/usr/bin/env bash
set -euo pipefail
if [ "$#" -lt 2 ]; then
  echo "usage: $0 <matrix.csv> <output-dir> [runtime-dir]" >&2
  exit 2
fi
PKG_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
MATRIX=$1
OUT=$2
RUNTIME_DIR=${3:-$PKG_ROOT/runtime}
python3 "$PKG_ROOT/runner/run_campaign.py" \
  --runtime-root "$RUNTIME_DIR/MIGRATION_PKG_2026-08-19" \
  --matrix "$MATRIX" --output-dir "$OUT" --workers 16
