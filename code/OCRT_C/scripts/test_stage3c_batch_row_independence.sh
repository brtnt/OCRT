#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
BIN=${OCRT_BIN:-./build/ocrt_v1.2}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/serial" "$TMP/parallel"
cat > "$TMP/serial.csv" <<CSV
out,wavelength
$TMP/serial/row1_443.csv,443
$TMP/serial/row2_555.csv,555
$TMP/serial/row3_660.csv,660
$TMP/serial/row4_443.csv,443
CSV
cat > "$TMP/parallel.csv" <<CSV
out,wavelength
$TMP/parallel/row1_443.csv,443
$TMP/parallel/row2_555.csv,555
$TMP/parallel/row3_660.csv,660
$TMP/parallel/row4_443.csv,443
CSV
COMMON=(--surface ocean --wind-speed 3 --sza 40 --wavelength 443 \
 --pressure 1013.25 --aod 0 --n-water 1.34 \
 --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0.5 --ocrt-adom440 0 \
 --ocrt-adom-slope 0.014 --ocrt-tsm-species red_clay \
 --n-mu-water 24 --water-m-max 8 \
 --lut-vza-step 10 --lut-vza-max 60 --lut-raa-step 30 \
 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 \
 --decouple-sunglint)
ENV=(OCRT_ADVANCED=1 OCRT_WATER_PARTICLE_KERNEL=direct \
     OCRT_WATER_VALUE_NPHI=128 OCRT_BATCH_DIAG=1)
env OMP_NUM_THREADS=1 "${ENV[@]}" "$BIN" "${COMMON[@]}" \
  --batch-full-grid "$TMP/serial.csv" >"$TMP/serial.out" 2>"$TMP/serial.err"
env OMP_NUM_THREADS=4 "${ENV[@]}" "$BIN" "${COMMON[@]}" \
  --batch-full-grid "$TMP/parallel.csv" >"$TMP/parallel.out" 2>"$TMP/parallel.err"
grep -F 'BATCH_DIAG VERDICT: PARALLEL_OK threads overlap as expected' \
  "$TMP/parallel.err" >/dev/null
grep -F '# batch-full-grid done: 4 rows, 0 failed' "$TMP/parallel.err" >/dev/null
# Every independent all-view row must use the native path, never cell replay.
if [[ $(grep -c '\[native-coupled-lut\].*replay=0' "$TMP/serial.err") -ne 4 ]]; then
  echo 'ERROR: serial batch did not complete four native all-view rows' >&2
  cat "$TMP/serial.err" >&2
  exit 1
fi
if [[ $(grep -c '\[native-coupled-lut\].*replay=0' "$TMP/parallel.err") -ne 4 ]]; then
  echo 'ERROR: parallel batch did not complete four native all-view rows' >&2
  cat "$TMP/parallel.err" >&2
  exit 1
fi
python - "$TMP" <<'PY'
from pathlib import Path
import csv, math, sys
root=Path(sys.argv[1])
for i,wl in ((1,443),(2,555),(3,660),(4,443)):
    a=root/'serial'/f'row{i}_{wl}.csv'
    b=root/'parallel'/f'row{i}_{wl}.csv'
    if a.read_bytes()!=b.read_bytes():
        with a.open(newline='') as fa, b.open(newline='') as fb:
            ra=list(csv.DictReader(fa)); rb=list(csv.DictReader(fb))
        if len(ra)!=len(rb) or (ra and ra[0].keys()!=rb[0].keys()):
            raise SystemExit(f'row {i}: structure mismatch')
        max_abs=0.0
        for xa,xb in zip(ra,rb):
            for key in xa:
                try: x=float(xa[key]); y=float(xb[key])
                except ValueError:
                    if xa[key]!=xb[key]: raise SystemExit(f'row {i}: text mismatch {key}')
                    continue
                if math.isnan(x) and math.isnan(y): continue
                max_abs=max(max_abs,abs(x-y))
        raise SystemExit(f'row {i}: serial/parallel mismatch max_abs={max_abs:.3e}')
print('serial/parallel rows byte-identical')
PY
echo 'PASS: independent batch rows are serial/parallel byte-identical; all-view only; no cross-row cache state'
