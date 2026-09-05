#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN=${1:-"$ROOT/build/ocrt"}
OUT=${2:-"$(mktemp)"}
CLEAN_OUT=0
if [[ $# -lt 2 ]]; then CLEAN_OUT=1; fi
trap 'if [[ $CLEAN_OUT -eq 1 ]]; then rm -f "$OUT"; fi' EXIT
python3 "$ROOT/tools/snf_mie_generator/run_snf_smoke.py" \
  --exe "$BIN" \
  --mie-dir "$ROOT/inputs" \
  --out "$OUT" \
  --workers "${SNF_SMOKE_WORKERS:-16}" \
  --timeout "${SNF_SMOKE_TIMEOUT:-60}"
