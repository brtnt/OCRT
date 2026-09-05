#!/usr/bin/env bash
# GOCI-III 대기보정 자료생산 (GPU, CuPy 필요). 사용법: ./run_produce_gpu.sh <grid.csv> <out.csv> [옵션...]
set -e
[ -z "$2" ] && { echo "사용법: ./run_produce_gpu.sh <grid.csv> <out.csv> [옵션...]"; exit 1; }
GRID="$1"; OUT="$2"; shift 2
OCRT_PY_GPU=1 python3 produce_grid.py --grid "$GRID" --out "$OUT" --data data \
  --n-mu-water 24 --nt-atm 400 --m-max 16 --max-it-water 500 "$@"
