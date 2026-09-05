#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
mkdir -p build
gcc -std=c11 -O2 -Isrc tests/test_pure_water_z09_table.c src/rt_water_iop.c -o build/test_pure_water_z09_table -lm
./build/test_pure_water_z09_table inputs/water_iop/water_coef_z09_1nm.txt
