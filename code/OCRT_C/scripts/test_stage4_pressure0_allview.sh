#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
BIN=${OCRT_BIN:-./build/ocrt_v1.2}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
COMMON=(--surface ocean --wind-speed 3 --sza 30 --wavelength 443 \
 --pressure 0 --aod 0 --water-model iop \
 --iop-a 0.1 --iop-b 0.5 --iop-bb 0.01 \
 --iop-mie-phase inputs/tsm_ahn/Red_clay_AHN.mie \
 --iop-mie-moment-mode gauss --iop-mie-moment-nmu 400 \
 --n-mu-water 24 --water-m-max 8 \
 --lut-vza-step 30 --lut-vza-max 60 --lut-raa-step 90 \
 --decouple-sunglint "${GAS[@]}")
ENV=(OMP_NUM_THREADS=1 OCRT_ADVANCED=1 OCRT_WATER_PARTICLE_KERNEL=direct \
     OCRT_WATER_VALUE_NPHI=512)
env "${ENV[@]}" "$BIN" "${COMMON[@]}" \
  --output-full-grid "$TMP/native.csv" >"$TMP/native.out" 2>"$TMP/native.err"
env "${ENV[@]}" OCRT_NATIVE_COUPLED_LUT_OFF=1 \
  OCRT_ALLOW_LEGACY_CELL_REPLAY=1 \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/replay.csv" \
  >"$TMP/replay.out" 2>"$TMP/replay.err"
grep -E '\[native-coupled-lut\].*calls=1.*water_cold=1.*water_views=3.*cells=12.*replay=0' \
  "$TMP/native.err" >/dev/null
grep -F 'diagnostic legacy cell replay enabled' "$TMP/replay.err" >/dev/null
grep -E 'replay=12' "$TMP/replay.err" >/dev/null
python - "$TMP/native.csv" "$TMP/replay.csv" <<'PY'
from pathlib import Path
import csv, math, sys
a=Path(sys.argv[1]); b=Path(sys.argv[2])
if a.read_bytes() == b.read_bytes():
    print('pressure0 native/replay CSV byte-identical')
    raise SystemExit(0)
with a.open(newline='') as fa, b.open(newline='') as fb:
    ra=list(csv.DictReader(fa)); rb=list(csv.DictReader(fb))
if len(ra)!=len(rb) or (ra and ra[0].keys()!=rb[0].keys()):
    raise SystemExit('pressure0 native/replay structure mismatch')
max_abs=0.0
for xa,xb in zip(ra,rb):
    for key in xa:
        try: x=float(xa[key]); y=float(xb[key])
        except ValueError:
            if xa[key]!=xb[key]: raise SystemExit(f'text mismatch {key}')
            continue
        if math.isnan(x) and math.isnan(y): continue
        d=abs(x-y); max_abs=max(max_abs,d)
        if d!=0.0: raise SystemExit(f'pressure0 mismatch {key}: {x} vs {y}')
print(f'pressure0 native/replay max_abs={max_abs:.3e}')
PY
echo 'PASS: pressure=0 native all-view uses one water cold solve; no production cell replay'
