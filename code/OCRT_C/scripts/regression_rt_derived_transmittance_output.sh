#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
FIXED=${1:-"$ROOT/build/ocrt"}
BASELINE=${2:-}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OCRT_ADVANCED=1
export OCRT_WATER_GRID_CACHE=0
export OCRT_ATM_GRID_CACHE=0

GAS_OFF=(
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
)
COMMON=(
  --surface ocean --water-model ocrt
  --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0
  --wind-speed 3 --sza 40 --vza 30 --raa 90 --wavelength 555
  --water-m-max 16
)

run_set() {
  local bin=$1 prefix=$2
  (
    cd "$ROOT"
    "$bin" "${COMMON[@]}" --rayleigh off "${GAS_OFF[@]}" \
      >"$TMP/${prefix}_noatm.out" 2>"$TMP/${prefix}_noatm.err"
    "$bin" "${COMMON[@]}" "${GAS_OFF[@]}" \
      >"$TMP/${prefix}_rayleigh.out" 2>"$TMP/${prefix}_rayleigh.err"
    "$bin" "${COMMON[@]}" \
      --mie inputs/M80C.mie --aod-865 0.1 --aer-l-max 80 --m-max 16 \
      --aer-h-km 2 --n-layers 40 "${GAS_OFF[@]}" \
      >"$TMP/${prefix}_rayaer.out" 2>"$TMP/${prefix}_rayaer.err"
    "$bin" --surface black --wind-speed 0 --sza 40 --vza 30 --raa 90 \
      --wavelength 555 "${GAS_OFF[@]}" \
      >"$TMP/${prefix}_black.out" 2>"$TMP/${prefix}_black.err"
  )
}

run_set "$FIXED" fixed
if [[ -n "$BASELINE" ]]; then run_set "$BASELINE" baseline; fi

(
  cd "$ROOT"
  OCRT_WATER_GRID_CACHE=1 OCRT_ATM_GRID_CACHE=1 \
  "$FIXED" --surface ocean --water-model ocrt \
    --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 --wind-speed 3 \
    --sza 40 --wavelength 555 --n-mu 8 --n-mu-water 8 --water-m-max 4 --rayleigh off \
    --lut-vza-step 60 --lut-vza-max 60 --lut-raa-step 180 \
    --output-full-grid "$TMP/ocean_grid.csv" "${GAS_OFF[@]}" \
    >"$TMP/ocean_grid.out" 2>"$TMP/ocean_grid.err"
)

python3 - "$TMP" "$BASELINE" "$ROOT/scripts/ocrt_ac_lut.py" <<'PY'
from __future__ import annotations
import csv
import math
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
have_baseline = bool(sys.argv[2])
lut_script = pathlib.Path(sys.argv[3])


def parse(path: pathlib.Path) -> dict[str, float | list[float]]:
    line = path.read_text().strip().splitlines()[-1]
    toks = line.split()
    out: dict[str, float | list[float]] = {"_lead": [float(x) for x in toks[:3]]}
    for tok in toks[3:]:
        if "=" not in tok:
            continue
        key, value = tok.split("=", 1)
        try:
            out[key] = float(value)
        except ValueError:
            pass
    return out


def close(a: float, b: float, rel: float = 2e-9, abs_: float = 5e-12) -> bool:
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


for name in ("noatm", "rayleigh", "rayaer"):
    d = parse(root / f"fixed_{name}.out")
    required = {
        "T_dir_dn", "T_diff_dn_hemi", "T_total_dn_hemi", "T_diff_dn_dir",
        "T_dir_up_view", "T_diff_up_view", "T_total_up_view",
        "TOA_water_signal_I", "T_up_rt_valid", "Lu0plus",
    }
    missing = sorted(required - d.keys())
    if missing:
        raise SystemExit(f"{name}: missing fields {missing}")
    if int(d["T_up_rt_valid"]) != 1:
        raise SystemExit(f"{name}: expected T_up_rt_valid=1")
    if not close(d["T_total_dn_hemi"], d["T_dir_dn"] + d["T_diff_dn_hemi"], 2e-10, 2e-11):
        raise SystemExit(f"{name}: downward RT identity failed")
    # stdout Lu has fewer printed digits than the transmission fields.
    if not close(d["T_total_up_view"] * d["Lu0plus"], d["TOA_water_signal_I"], 5e-6, 8e-10):
        raise SystemExit(f"{name}: upward RT identity failed")
    if name == "noatm":
        for key, expected in (("T_total_dn_hemi", 1.0), ("T_total_up_view", 1.0), ("T_diff_up_view", 0.0)):
            if not close(d[key], expected, 0.0, 2e-12):
                raise SystemExit(f"noatm: {key}={d[key]} expected={expected}")

black = parse(root / "fixed_black.out")
if int(black.get("T_up_rt_valid", -1)) != 0:
    raise SystemExit("black: upward T must be invalid without an upward BOA source")
for key in ("T_diff_up_view", "T_total_up_view", "TOA_water_signal_I"):
    if not math.isnan(float(black[key])):
        raise SystemExit(f"black: {key} must be NaN")

if have_baseline:
    ignored_prefixes = ("T_", "diag_", "TOA_water_signal_I")
    for name in ("noatm", "rayleigh", "rayaer", "black"):
        old = parse(root / f"baseline_{name}.out")
        new = parse(root / f"fixed_{name}.out")
        if old["_lead"] != new["_lead"]:
            raise SystemExit(f"{name}: leading TOA rho changed")
        for key in sorted(set(old) & set(new)):
            if key == "_lead" or key.startswith(ignored_prefixes):
                continue
            a, b = float(old[key]), float(new[key])
            if math.isnan(a) and math.isnan(b):
                continue
            if a != b:
                raise SystemExit(f"{name}: non-trans output changed: {key}: {a} -> {b}")

with (root / "ocean_grid.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
if not rows:
    raise SystemExit("ocean full-grid output is empty")
required_grid = {
    "T_dir_dn", "T_diff_dn_hemi", "T_total_dn_hemi", "T_diff_dn_dir",
    "T_dir_up_view", "T_diff_up_view", "T_total_up_view",
    "TOA_water_signal_I", "T_up_rt_valid", "Lu0plus_I",
}
missing = sorted(required_grid - rows[0].keys())
if missing:
    raise SystemExit(f"ocean grid missing exact T columns: {missing}")
for i, row in enumerate(rows):
    vals = {k: float(row[k]) for k in required_grid}
    if int(vals["T_up_rt_valid"]) != 1:
        raise SystemExit(f"grid row {i}: invalid upward RT")
    if not close(vals["T_total_dn_hemi"], vals["T_dir_dn"] + vals["T_diff_dn_hemi"], 5e-10, 5e-12):
        raise SystemExit(f"grid row {i}: downward identity failed")
    if not close(vals["T_total_up_view"] * vals["Lu0plus_I"], vals["TOA_water_signal_I"], 5e-10, 5e-12):
        raise SystemExit(f"grid row {i}: upward identity failed")

text = lut_script.read_text()
if "Tup=T[vza]" in text or "reciprocity THEOREM, rigorous" in text:
    raise SystemExit("production AC LUT script still contains reciprocity-derived Tup")
if 'Tup_rt_valid=0' not in text:
    raise SystemExit("production AC LUT script does not invalidate undefined upward T")

print(f"RT_DERIVED_TRANSMITTANCE_OUTPUT PASS cases=4 grid_rows={len(rows)} baseline_compare={int(have_baseline)}")
PY
