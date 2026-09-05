#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

CC_BIN=${CC:-gcc}
OUT=${1:-build/test_phase_wavelength_pchip}
mkdir -p "$(dirname "$OUT")"

"$CC_BIN" -std=c11 -O2 -Wall -Wextra -pedantic \
  -Isrc \
  tests/test_phase_wavelength_pchip.c \
  src/shared/numerics.c src/shared/io_utils.c src/shared/mie_io.c \
  -o "$OUT" -lm

"$OUT"

if grep -RFn "linterp(mie->phase_wavelengths" \
    src/shared/mie_io.c src/rt_aerosol_runtime.c >/dev/null; then
  echo "FAIL: linear phase-wavelength interpolation call remains" >&2
  exit 1
fi

count=$(grep -Rh "pchip_interp(mie->phase_wavelengths" \
    src/shared/mie_io.c src/rt_aerosol_runtime.c | wc -l)
if [[ "$count" -ne 5 ]]; then
  echo "FAIL: expected five legacy/helper Mie phase-wavelength PCHIP call sites, found $count" >&2
  exit 1
fi

grep -F "mie_phase_nodes_at_wavelength_linear" src/rt_water_rt.c >/dev/null
echo "PASS source audit: 5 legacy/helper PCHIP call sites; Stage-3B direct path uses common-weight linear wavelength interpolation"
