#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
mkdir -p build
mapfile -t SRC < <(find src -name '*.c' ! -name main.c | LC_ALL=C sort)
${CC:-gcc} -std=c11 -O3 -DOCRT_FAST_KERNELS -fopenmp -Isrc \
  tests/test_stage3b_direct_kernel_dump.c "${SRC[@]}" \
  -o build/test_stage3b_direct_kernel_dump -lm
