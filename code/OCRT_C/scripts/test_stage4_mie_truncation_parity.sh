#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PYROOT=${1:-${OCRT_PY_ROOT:-}}
if [[ -z "$PYROOT" || ! -d "$PYROOT/ocrt_py" ]]; then
  echo "usage: $0 /path/to/pyOCRT (or set OCRT_PY_ROOT)" >&2
  exit 2
fi
cd "$ROOT"
mkdir -p build
mapfile -t SRC < <(find src -name '*.c' ! -name main.c | LC_ALL=C sort)
${CC:-gcc} -std=c11 -O3 -DOCRT_FAST_KERNELS -fopenmp -Isrc \
  tests/dump_stage4_mie_truncation.c "${SRC[@]}" \
  -o build/dump_stage4_mie_truncation -lm
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

CASES=(
  "inputs/water_iop/Detritus_Stramski2001.mie|443|detritus"
  "inputs/tsm_ahn/Red_clay_AHN.mie|443|tsm"
  "inputs/water_iop/eap/EAP_06_Dinoflagellates_D24.mie|443|eap_d24"
)
for spec in "${CASES[@]}"; do
  IFS='|' read -r mie wl name <<<"$spec"
  ./build/dump_stage4_mie_truncation "$mie" "$wl" 1 "$TMP/$name.csv"
done
PYTHONPATH="$PYROOT" python3 - "$PYROOT" "$TMP" <<'PY'
from __future__ import annotations
import csv
from pathlib import Path
import sys
import numpy as np

pyroot = Path(sys.argv[1])
tmp = Path(sys.argv[2])
sys.path.insert(0, str(pyroot))
from ocrt_py.aerosol import read_mie
from ocrt_py.value_phase import value_phase_at_wavelength

cases = [
    ("detritus", pyroot / "data/water_iop/Detritus_Stramski2001.mie", 443.0),
    ("tsm", pyroot / "data/tsm_ahn/Red_clay_AHN.mie", 443.0),
    ("eap_d24", pyroot / "data/water_iop/eap/EAP_06_Dinoflagellates_D24.mie", 443.0),
]
max_all = 0.0
for name, path, wl in cases:
    with (tmp / f"{name}.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    theta = np.array([float(r["theta_deg"]) for r in rows])
    c_arrays = [np.array([float(r[k]) for r in rows]) for k in ("p11", "p12", "p33")]
    c_A = float(rows[0]["A"]); c_norm = float(rows[0]["norm_before"])
    c_g = float(rows[0]["g"]); c_bb = float(rows[0]["bb_b"])
    mie = read_mie(str(path))
    raw = value_phase_at_wavelength(mie, wl)
    phase, A = raw.loglinear_truncated()
    errs = [float(np.max(np.abs(theta - phase.theta_deg)))]
    errs.extend(float(np.max(np.abs(c - p))) for c, p in zip(c_arrays, (phase.p11, phase.p12, phase.p33)))
    errs.extend([abs(c_A - A), abs(c_norm - phase.norm_before),
                 abs(c_g - phase.g_asym), abs(c_bb - phase.bb_b_ratio)])
    err = max(errs); max_all = max(max_all, err)
    tol = 2.0e-11 if name == "detritus" else 5.0e-12
    if err > tol:
        raise SystemExit(f"FAIL {name}: max_abs={err:.17g}, errors={errs}")
    print(f"PASS {name}: A={A:.15g} g={phase.g_asym:.15g} bb_b={phase.bb_b_ratio:.15g} max_abs={err:.3e}")
print(f"PASS C/Python actual-file truncation parity max_abs={max_all:.3e}")
PY
