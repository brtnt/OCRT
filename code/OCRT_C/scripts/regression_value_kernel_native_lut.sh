#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
BIN=${OCRT_BIN:-./build/ocrt}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
COMMON=(--surface ocean --wind-speed 3 --sza 30 --wavelength 555 \
 --pressure 1013.25 --aod 0 --water-model iop \
 --iop-a 0.1 --iop-b 0.5 --iop-bb 0.01 \
 --iop-mie-phase inputs/tsm_ahn/Red_clay_AHN.mie \
 --iop-mie-moment-mode gauss --iop-mie-moment-nmu 400 \
 --n-mu-water 8 --water-m-max 4 --lut-vza-step 30 --lut-vza-max 60 \
 --lut-raa-step 90 "${GAS[@]}")
ENV=(OMP_NUM_THREADS=1 OCRT_ADVANCED=1 OCRT_WATER_VALUE_KERNEL_POL=1 \
     OCRT_WATER_VALUE_NPHI=256)
env "${ENV[@]}" "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/native.csv" \
  >"$TMP/native.out" 2>"$TMP/native.err"
env "${ENV[@]}" OCRT_NATIVE_COUPLED_LUT_OFF=1 OCRT_ALLOW_LEGACY_CELL_REPLAY=1 \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/legacy.csv" \
  >"$TMP/legacy.out" 2>"$TMP/legacy.err"
python - "$TMP/native.csv" "$TMP/legacy.csv" <<'PY'
import sys
import numpy as np
import pandas as pd
a = pd.read_csv(sys.argv[1])
b = pd.read_csv(sys.argv[2])
if list(a.columns) != list(b.columns) or a.shape != b.shape:
    raise SystemExit('native/legacy CSV structure mismatch')
num = a.select_dtypes(include=[np.number]).columns
x = a[num].to_numpy(float); y = b[num].to_numpy(float)
if not np.allclose(x, y, rtol=5e-13, atol=5e-15, equal_nan=True):
    d = np.nanmax(np.abs(x-y))
    raise SystemExit(f'native/legacy numeric mismatch max_abs={d:.3e}')
print(f'max_abs={np.nanmax(np.abs(x-y)):.3e}')
PY
grep -F '[native-coupled-lut] calls=3 water_cold=1 water_views=3' "$TMP/native.err" >/dev/null
grep -F 'replay=0' "$TMP/native.err" >/dev/null
echo 'PASS: polarized value kernel preserves native coupled LUT within floating-point tolerance; cell replay=0'
