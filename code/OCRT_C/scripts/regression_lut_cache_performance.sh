#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-$ROOT/build/ocrt}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
COMMON=(--surface ocean --water-model ocrt --ocrt-chl 1.2 --ocrt-tsm 3.5
 --ocrt-adom440 0.04 --ocrt-tsm-species red_clay
 --sza 30 --wind-speed 3 --wavelength 490 --pressure 1013.25
 --mie "$ROOT/inputs/C50.mie" --aod-865 0.2 --n-mu-water 16
 --n-layers 80 --m-max 3 --n-mu 16 --lut-vza-step 30
 --lut-vza-max 30 --lut-raa-step 180)
(
 cd "$ROOT"
 OCRT_ADVANCED=1 OMP_NUM_THREADS=1 "$BIN" "${COMMON[@]}" \
   --output-full-grid "$TMP/on.csv" >/dev/null 2>"$TMP/on.err"
 OCRT_ADVANCED=1 OMP_NUM_THREADS=1 OCRT_S17A_OFF=1 "$BIN" "${COMMON[@]}" \
   --output-full-grid "$TMP/off.csv" >/dev/null 2>"$TMP/off.err"
)
python3 - "$TMP/on.csv" "$TMP/off.csv" <<'PY'
import csv
import sys

on_path, off_path = sys.argv[1:]
# This diagnostic is populated only by the cache-assisted path; the forced
# cache-off reference intentionally leaves it at zero.  All physical outputs
# and convergence metadata must remain byte-identical.
ignored = {"T_diff_dn_dir"}
with open(on_path, newline="") as f:
    on_rows = list(csv.DictReader(f))
with open(off_path, newline="") as f:
    off_rows = list(csv.DictReader(f))
if len(on_rows) != len(off_rows):
    raise SystemExit(f"row-count mismatch: {len(on_rows)} != {len(off_rows)}")
if on_rows and off_rows and list(on_rows[0]) != list(off_rows[0]):
    raise SystemExit("CSV header mismatch")
for row_index, (on, off) in enumerate(zip(on_rows, off_rows)):
    for key in on:
        if key in ignored:
            continue
        if on[key] != off[key]:
            raise SystemExit(
                f"row {row_index} field {key} mismatch: {on[key]} != {off[key]}"
            )
print(
    "PASS LUT cache physical-output byte comparison "
    "(diagnostic T_diff_dn_dir intentionally excluded)"
)
PY
