#!/usr/bin/env bash
set -euo pipefail

BIN=${1:-./build/ocrt}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"

export OCRT_ADVANCED=1
export OMP_NUM_THREADS=1
export LC_ALL=C

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

GAS_OFF=(
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
)

get_key() {
  local key=$1 file=$2
  grep -o "${key}=[^ ]*" "$file" | head -n1 | cut -d= -f2
}

close() {
  python3 - "$1" "$2" "$3" <<'PY'
import math,sys
x=float(sys.argv[1]); y=float(sys.argv[2]); tol=float(sys.argv[3])
if not math.isfinite(x) or abs(x-y) > tol*max(1.0,abs(y)):
    raise SystemExit(f"not close: got {x:.17g}, expected {y:.17g}, tol={tol}")
PY
}

# 1. Gas-only ocean helper must use spherical gas-layer path, not return the
# plane-parallel shortcut when Rayleigh and aerosol are absent.
COMMON_OCEAN=(
  --surface ocean --wind-speed 0
  --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0
  --sza 85 --vza 0 --raa 0 --wavelength 555
  --surface-pressure 0 --n-layers 400
)
"$BIN" "${COMMON_OCEAN[@]}" >"$TMP/gas_pp.out" 2>"$TMP/gas_pp.err"
"$BIN" "${COMMON_OCEAN[@]}" --pssa >"$TMP/gas_pssa.out" 2>"$TMP/gas_pssa.err"
T_PP=$(get_key T_dir_dn "$TMP/gas_pp.out")
T_SA=$(get_key T_dir_dn "$TMP/gas_pssa.out")
close "$T_PP" 6.8818713780e-01 5e-10
close "$T_SA" 7.6142396567e-01 5e-10
python3 - "$T_PP" "$T_SA" <<'PY'
import sys
pp,sa=map(float,sys.argv[1:])
if not sa > pp:
    raise SystemExit(f"gas-only PSSA did not shorten the grazing path: PP={pp}, PSSA={sa}")
PY

# 2. Flat reflected-beam source, including the downward first-order field and
# its higher-order coupling.  The first-order value is a useful boundary: the
# added downward field cannot affect TOA until order 2.
BASE=(
  --surface flat --wind-speed 0
  --sza 85 --vza 30 --raa 90 --wavelength 412
  --pssa --n-layers 400 "${GAS_OFF[@]}"
)
"$BIN" "${BASE[@]}" --sos-max-orders 1 >"$TMP/o1.out" 2>"$TMP/o1.err"
"$BIN" "${BASE[@]}" --sos-max-orders 2 >"$TMP/o2.out" 2>"$TMP/o2.err"
"$BIN" "${BASE[@]}" --sos-max-orders 20 >"$TMP/o20.out" 2>"$TMP/o20.err"
close "$(get_key TOA_rho_I "$TMP/o1.out")" 2.2076059458e-01 5e-10
close "$(get_key TOA_rho_I "$TMP/o2.out")" 2.9703812927e-01 5e-10
close "$(get_key TOA_rho_I "$TMP/o20.out")" 3.4989871055e-01 5e-10
close "$(get_key TOA_rho_Q "$TMP/o20.out")" 2.6477000904e-01 5e-10
close "$(get_key TOA_rho_U "$TMP/o20.out")" -2.2308875845e-02 5e-10

# 3. Direct spherical reflected-beam evaluation must remain finite in a
# saturated O2 band where a baseline*inverse-correction product can form 0*inf.
"$BIN" --surface flat --wind-speed 0 \
  --sza 85 --vza 60 --raa 90 --wavelength 760 \
  --surface-pressure 1013.25 --pssa --n-layers 400 \
  >"$TMP/saturated.out" 2>"$TMP/saturated.err"
python3 - "$TMP/saturated.out" <<'PYCHECK'
import math, pathlib, sys
line = pathlib.Path(sys.argv[1]).read_text().strip().splitlines()[-1]
vals = {}
for tok in line.split()[3:]:
    if "=" not in tok:
        continue
    k, v = tok.split("=", 1)
    try:
        vals[k] = float(v)
    except ValueError:
        pass
physical = (
    "TOA_rho_I", "TOA_rho_Q", "TOA_rho_U",
    "T_dir_dn", "T_diff_dn_hemi", "T_total_dn_hemi", "T_diff_dn_dir",
    "T_dir_up_view", "tau_R", "AOD_ref", "AOD_band",
)
for key in physical:
    if key not in vals or not math.isfinite(vals[key]):
        raise SystemExit(f"non-finite physical value in saturated-band PSSA output: {key}={vals.get(key)}")
# The RT-derived transmittance output fix deliberately invalidates upward
# transmission when no upward BOA source was solved. These NaNs are metadata,
# not a radiative instability.
if int(vals.get("T_up_rt_valid", -1)) != 0:
    raise SystemExit("standalone flat case must report T_up_rt_valid=0")
for key in ("T_diff_up_view", "T_total_up_view", "TOA_water_signal_I"):
    if key not in vals or not math.isnan(vals[key]):
        raise SystemExit(f"undefined upward diagnostic must be NaN: {key}={vals.get(key)}")
PYCHECK

# 4. Controls: no reflected-beam correction for a black lower boundary, and
# no beta correction outside PSSA.  Pin representative current results.
"$BIN" --surface black --sza 85 --vza 60 --raa 90 --wavelength 412 \
  --pssa --n-layers 400 "${GAS_OFF[@]}" >"$TMP/black.out" 2>"$TMP/black.err"
close "$(get_key TOA_rho_I "$TMP/black.out")" 5.6296257913e-01 5e-10
"$BIN" --surface flat --wind-speed 0 --sza 70 --vza 30 --raa 90 \
  --wavelength 443 --n-layers 100 "${GAS_OFF[@]}" \
  >"$TMP/pp.out" 2>"$TMP/pp.err"
close "$(get_key TOA_rho_I "$TMP/pp.out")" 1.6789952072e-01 5e-10

printf 'PSSA_NUMERIC_FIXES pass=4 fail=0\n'
