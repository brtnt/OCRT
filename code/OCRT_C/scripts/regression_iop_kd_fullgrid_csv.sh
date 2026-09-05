#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN=${1:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

export OCRT_ADVANCED=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

COMMON=(
  --surface ocean --water-model ocrt
  --ocrt-chl 1.2 --ocrt-tsm 3.5 --ocrt-adom440 0.04
  --ocrt-tsm-species red_clay
  --sza 30 --wind-speed 3 --wavelength 490 --pressure 1013.25
  --mie inputs/C50.mie --aod-865 0.2
  --n-mu-water 16 --n-layers 100 --m-max 3
)

(
  cd "$ROOT"
  "$BIN" "${COMMON[@]}" --vza 30 --raa 90 \
    >"$TMP/single.out" 2>"$TMP/single.err"
  "$BIN" "${COMMON[@]}" --n-mu 16 \
    --lut-vza-step 30 --lut-vza-max 30 --lut-raa-step 180 \
    --output-full-grid "$TMP/grid.csv" \
    >"$TMP/grid.out" 2>"$TMP/grid.err"
)

python3 - "$TMP" <<'PY'
from __future__ import annotations
import csv
import math
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
text = (root / "single.out").read_text()


def val(name: str) -> float:
    m = re.search(rf"(?<![A-Za-z0-9_]){re.escape(name)}=([+-]?[0-9.]+(?:[eE][+-]?[0-9]+)?)", text)
    if not m:
        raise SystemExit(f"single output missing {name}")
    return float(m.group(1))

# Stage-2 absorption-only phytoplankton + detritus scattering reference.
# Pure-water a_w/b_w anchors were regenerated for the corrected 2026-07-28
# 200--2449 nm Z09-derived table; non-water constituent anchors are unchanged.
kd = val("Kd0minus")
if not (kd > 0.0 and math.isclose(kd, 0.360547, rel_tol=0.0, abs_tol=8e-4)):
    raise SystemExit(f"Kd0minus={kd:.9g}, expected positive ~0.360547")
if not math.isclose(val("a_chl"), 3.0672e-2, rel_tol=2e-9, abs_tol=2e-11):
    raise SystemExit("a_chl does not match absorption-only phytoplankton")
if not math.isclose(val("a_pig"), 3.0672e-2, rel_tol=2e-9, abs_tol=2e-11):
    raise SystemExit("a_phyto_detritus/a_pig reference mismatch")
for name, expected, tol in (
    ("TOA_rho_I", 1.7907922413e-1, 2e-11),
    ("Rrs0plus_I", 3.837254e-2, 5e-9),
    ("rrs0minus_I", 6.428147e-2, 5e-9),
    ("Ed0minus", 8.272498e-1, 5e-7),
    ("Lu0minus", 5.317683e-2, 5e-8),
):
    got = val(name)
    if not math.isclose(got, expected, rel_tol=0.0, abs_tol=tol):
        raise SystemExit(f"{name}={got:.12g}, expected {expected:.12g}")

with (root / "grid.csv").open(newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
if len(rows) != 4:
    raise SystemExit(f"full-grid row count={len(rows)}, expected 4")

old_prefix = [
    "case_kind","sza_deg","wavelength_nm","wind_speed_ms","n_mu","mu_index",
    "mu_view","vza_deg","raa_deg","Ed0plus_air","Ed0plus_direct",
    "Ed0plus_diffuse","Lu0plus_I","Lu0plus_Q","Lu0plus_U","Rrs_I","Rrs_Q",
    "Rrs_U","Ed0minus_water","Ed0minus_direct","Ed0minus_diffuse",
    "Eu0minus_water","Eu0minus_direct","Eu0minus_diffuse","Lu0minus_I",
    "Lu0minus_Q","Lu0minus_U","rrs_I","rrs_Q","rrs_U","TOA_rho_I","TOA_rho_Q","TOA_rho_U",
    "T_dir_dn","T_diff_dn_hemi","T_total_dn_hemi","T_diff_dn_dir",
    "T_dir_up_view","T_diff_up_view","T_total_up_view","TOA_water_signal_I",
    "T_up_rt_valid","a_total","b_total","bb_total","omega_water",
    "water_orders","water_converged",
]
new_tail = [
    "a_w","b_w","bb_w","a_chl","a_phyto_detritus","b_phyto_detritus",
    "bb_phyto_detritus","a_dom","a_min","b_min","bb_min","Kd0minus",
]
if reader.fieldnames != old_prefix + new_tail:
    raise SystemExit("full-grid header/order mismatch")

expected = {
    "a_w": 1.502097e-2,
    "b_w": 2.824033e-3,
    "bb_w": 1.4120165e-3,
    "a_chl": 3.0672e-2,
    "a_phyto_detritus": 3.0672e-2,
    "a_dom": 1.986341215166e-2,
    "a_min": 1.904e-1,
    "b_min": 2.916655,
    "bb_min": 6.534415253044e-2,
}
for i, row in enumerate(rows):
    for name, ref in expected.items():
        got = float(row[name])
        if not math.isclose(got, ref, rel_tol=3e-9, abs_tol=3e-11):
            raise SystemExit(f"grid row {i} {name}={got} expected={ref}")
    if not float(row["Kd0minus"]) > 0.0:
        raise SystemExit(f"grid row {i} Kd0minus is not positive")

print(f"IOP_KD_FULLGRID_CSV PASS single_Kd={kd:.6f} grid_rows={len(rows)} cols={len(reader.fieldnames)}")
PY
