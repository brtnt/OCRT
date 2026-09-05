#!/usr/bin/env bash
# GOCI-III 대기보정 자료생산 (CPU). 사용법: ./run_produce.sh <grid.csv> <out.csv> [옵션...]
set -e
[ -z "$2" ] && { echo "사용법: ./run_produce.sh <grid.csv> <out.csv> [옵션...]"; exit 1; }
GRID="$1"; OUT="$2"; shift 2
OCRT_PY_GPU=0 python3 produce_grid.py --grid "$GRID" --out "$OUT" --data data \
  --n-mu-water 24 --nt-atm 400 --m-max 16 --max-it-water 500 "$@"
