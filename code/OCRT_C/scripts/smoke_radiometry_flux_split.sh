#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN=${1:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OCRT_ADVANCED=1

GAS_OFF=(
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
)
COMMON_WATER=(
  --surface ocean --water-model ocrt
  --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0
  --wind-speed 0 --vza 20 --raa 90
)

(
  cd "$ROOT"
  "$BIN" "${COMMON_WATER[@]}" --sza 30 --wavelength 555 --rayleigh off \
    "${GAS_OFF[@]}" >"$TMP/noatm.out" 2>"$TMP/noatm.err"

  "$BIN" "${COMMON_WATER[@]}" --sza 40 --wavelength 555 \
    --surface-pressure 1013.25 "${GAS_OFF[@]}" \
    >"$TMP/atm.out" 2>"$TMP/atm.err"

  "$BIN" "${COMMON_WATER[@]}" --sza 80 --wavelength 443 \
    --surface-pressure 1013.25 --pssa "${GAS_OFF[@]}" \
    >"$TMP/pssa.out" 2>"$TMP/pssa.err"

  "$BIN" --surface ocean --water-model ocrt \
    --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 --wind-speed 0 \
    --sza 30 --wavelength 555 --rayleigh off "${GAS_OFF[@]}" \
    --lut-vza-step 85 --lut-raa-step 180 \
    --output-full-grid "$TMP/grid.csv" \
    >"$TMP/grid.stdout" 2>"$TMP/grid.err"
)

python3 - "$TMP" <<'PY'
from __future__ import annotations
import csv
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
keys = {
    "Ed0plus_direct", "Ed0plus_diffuse",
    "Ed0minus_direct", "Ed0minus_diffuse",
    "Eu0minus_direct", "Eu0minus_diffuse",
    "T_dir_dn", "T_diff_dn_hemi", "T_total_dn_hemi",
    "T_dir_up_view", "T_diff_up_view", "T_total_up_view",
    "TOA_water_signal_I", "T_up_rt_valid", "Lu0plus",
}

def parse_line(path: pathlib.Path) -> dict[str, float]:
    line = path.read_text().strip().splitlines()[-1]
    out: dict[str, float] = {}
    for token in line.split():
        if "=" not in token:
            continue
        k, v = token.split("=", 1)
        try:
            out[k] = float(v)
        except ValueError:
            pass
    missing = sorted(keys - out.keys())
    if missing:
        raise SystemExit(f"{path.name}: missing fields: {missing}")
    return out

def check_printed_sum(d: dict[str, float], total: str, direct: str, diffuse: str) -> None:
    # Single-line diagnostics use %.6e. Permit only formatting-level residual.
    err = abs(d[total] - (d[direct] + d[diffuse]))
    tol = 5.0e-6 * max(1.0, abs(d[total]))
    if err > tol:
        raise SystemExit(f"{total} != {direct}+{diffuse}: err={err:.3e}, tol={tol:.3e}")

for name in ("noatm", "atm", "pssa"):
    d = parse_line(root / f"{name}.out")
    check_printed_sum(d, "Ed0plus", "Ed0plus_direct", "Ed0plus_diffuse")
    check_printed_sum(d, "Ed0minus", "Ed0minus_direct", "Ed0minus_diffuse")
    check_printed_sum(d, "Eu0minus", "Eu0minus_direct", "Eu0minus_diffuse")
    if d["Eu0minus_direct"] != 0.0:
        raise SystemExit(f"{name}: Eu0minus_direct must be zero in the current black-bottom model")
    if abs(d["T_total_dn_hemi"] - (d["T_dir_dn"] + d["T_diff_dn_hemi"])) > 2.0e-10:
        raise SystemExit(f"{name}: T_total_dn_hemi identity failed")
    if int(d["T_up_rt_valid"]) != 1:
        raise SystemExit(f"{name}: expected rigorous upward RT output")
    if abs(d["T_total_up_view"] * d["Lu0plus"] - d["TOA_water_signal_I"]) > 8.0e-10:
        raise SystemExit(f"{name}: T_total_up_view RT identity failed")

noatm = parse_line(root / "noatm.out")
atm = parse_line(root / "atm.out")
if abs(noatm["Ed0plus_diffuse"]) > 1.0e-12:
    raise SystemExit("no-atmosphere case must have Ed0plus_diffuse=0")
if not (atm["Ed0plus_diffuse"] > 0.0):
    raise SystemExit("atmosphere case must have positive Ed0plus_diffuse")

with (root / "grid.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
if not rows:
    raise SystemExit("full-grid output is empty")
required = {
    "Ed0plus_air", "Ed0plus_direct", "Ed0plus_diffuse",
    "Ed0minus_water", "Ed0minus_direct", "Ed0minus_diffuse",
    "Eu0minus_water", "Eu0minus_direct", "Eu0minus_diffuse",
}
missing = sorted(required - rows[0].keys())
if missing:
    raise SystemExit(f"full-grid missing columns: {missing}")
for row in rows:
    for total, direct, diffuse in (
        ("Ed0plus_air", "Ed0plus_direct", "Ed0plus_diffuse"),
        ("Ed0minus_water", "Ed0minus_direct", "Ed0minus_diffuse"),
        ("Eu0minus_water", "Eu0minus_direct", "Eu0minus_diffuse"),
    ):
        a, b, c = map(float, (row[total], row[direct], row[diffuse]))
        if abs(a - (b + c)) > 2.0e-11 * max(1.0, abs(a)):
            raise SystemExit(f"grid identity failed: {total}={a}, {direct}+{diffuse}={b+c}")

print(f"RADIOMETRY_FLUX_SPLIT PASS single_cases=3 grid_rows={len(rows)}")
PY
