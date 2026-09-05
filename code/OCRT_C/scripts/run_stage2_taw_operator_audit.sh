#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 <OSOAA_TAW_binary> <RAD_UsedAngles.txt> [output_dir]" >&2
  exit 2
fi

ROOT=$(cd "$(dirname "$0")/.." && pwd)
TAW=$1
ANGLES=$2
OUT=${3:-"$ROOT/validation/stage2_taw_operator_audit_2026-07-24/run"}

python3 "$ROOT/tools/stage2_taw_operator_audit.py" \
  --ocrt-root "$ROOT" \
  --osoaa-taw "$TAW" \
  --angles "$ANGLES" \
  --outdir "$OUT" \
  --m-max 4 \
  --n-phi 1024 \
  --wind 3.0 \
  --n-water 1.34 \
  --sigma-type 1 \
  --q-convention 1
