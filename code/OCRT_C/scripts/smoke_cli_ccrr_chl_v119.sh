#!/usr/bin/env bash
# OCRT v1.19 CCRR Chl table, interpolation and execution smoke test.
set -euo pipefail

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/.." && pwd)
BIN=${1:-$ROOT/build/v2_solver_vk_v1.19}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

[[ -x "$BIN" ]] || { echo "missing executable: $BIN" >&2; exit 2; }
BIN=$(realpath "$BIN")

pass=0
fail=0
ok(){ printf 'PASS %s\n' "$1"; pass=$((pass+1)); }
bad(){ printf 'FAIL %s\n' "$1"; fail=$((fail+1)); }

CANON=$ROOT/inputs/water_iop/aph_ccrr_morel1988_mm01.txt
LEGACY=$ROOT/inputs/water_iop/aph_bricaud_1998.txt
SOURCE=$ROOT/inputs/water_iop/source/apstarchl_morel1988_normalized_user_supplied.txt

if cmp -s "$CANON" "$LEGACY"; then ok canonical_legacy_byte_identity; else bad canonical_legacy_byte_identity; fi

python3 "$ROOT/scripts/build_ccrr_chl_table.py" "$SOURCE" "$TMP/generated.txt"
if cmp -s "$CANON" "$TMP/generated.txt"; then ok table_reproducible; else bad table_reproducible; fi

if python3 - "$SOURCE" "$CANON" <<'PY'
import math, sys
from pathlib import Path
src, out = map(Path, sys.argv[1:])
source=[]
for line in src.read_text().splitlines():
    f=line.split()
    if len(f)<2: continue
    try: wl=float(f[0]); a=float(f[1])
    except ValueError: continue
    if wl<0 and a<0: break
    source.append((wl,a))
table=[]
for line in out.read_text().splitlines():
    if not line.strip() or line.lstrip().startswith('!'): continue
    wl,ap,ep,A,E=map(float,line.split())
    table.append((wl,ap,ep,A,E))
assert len(source)==len(table)==141
for (wl,norm),(wl2,ap,ep,A,E) in zip(source,table):
    assert wl==wl2 and ap==0 and ep==0
    assert math.isclose(A,0.06*norm,rel_tol=0,abs_tol=5e-13)
    assert E==0.65
assert next(A for wl,_,_,A,_ in table if wl==440)==0.06
PY
then ok table_values; else bad table_values; fi

GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
FAST=(--debug-water-max-orders 4 --n-mu-water 12 --water-m-max 4)

run_case(){
  local wl=$1 chl=$2 out=$3
  (cd "$ROOT" && OCRT_DEBUG=1 OCRT_ADVANCED=1 "$BIN" \
    --surface ocean --water-model ccrr --wind-speed 0 \
    --ccrr-chl "$chl" --ccrr-tsm 0 --ccrr-adom440 0 \
    --sza 30 --vza 20 --raa 90 --wavelength "$wl" --pressure 0 \
    "${GAS[@]}" "${FAST[@]}") >"$out" 2>"$out.err"
}

# Verify canonical-first loading and legacy-name fallback in isolated working trees.
mkdir -p "$TMP/canonical" "$TMP/legacy"
cp -a "$ROOT/inputs" "$TMP/canonical/inputs"
cp -a "$ROOT/inputs" "$TMP/legacy/inputs"
rm -f "$TMP/canonical/inputs/water_iop/aph_bricaud_1998.txt"
rm -f "$TMP/legacy/inputs/water_iop/aph_ccrr_morel1988_mm01.txt"
loader_case(){
  local wd=$1 out=$2
  (cd "$wd" && OCRT_DEBUG=1 OCRT_ADVANCED=1 "$BIN" \
    --surface ocean --water-model ccrr --wind-speed 0 \
    --ccrr-chl 0.3 --ccrr-tsm 0 --ccrr-adom440 0 \
    --sza 30 --vza 20 --raa 90 --wavelength 443 --pressure 0 \
    "${GAS[@]}" "${FAST[@]}") >"$out" 2>"$out.err"
}
if loader_case "$TMP/canonical" "$TMP/canonical.out"; then ok canonical_only_load; else bad canonical_only_load; fi
if loader_case "$TMP/legacy" "$TMP/legacy.out"; then ok legacy_fallback_load; else bad legacy_fallback_load; fi
if cmp -s "$TMP/canonical.out" "$TMP/legacy.out" && cmp -s "$TMP/canonical.out.err" "$TMP/legacy.out.err"; then
  ok canonical_legacy_runtime_identity
else
  bad canonical_legacy_runtime_identity
fi

for spec in '440 1.0' '443 0.3' '555 0.3' '670 0.3'; do
  set -- $spec
  if run_case "$1" "$2" "$TMP/case_${1}_${2}.out"; then ok "run_${1}_${2}"; else bad "run_${1}_${2}"; fi
done

if python3 - "$TMP" <<'PY'
import math,re,sys
from pathlib import Path
p=Path(sys.argv[1])
def value(path,key):
    txt=path.read_text()
    m=re.search(rf'(?:^|\s){re.escape(key)}=([+-]?[0-9.]+(?:e[+-]?[0-9]+)?)',txt,re.I)
    assert m, (path,key)
    return float(m.group(1))
cases=[
    ('case_440_1.0.out',0.060000000),
    ('case_443_0.3.out',0.06*(1+3/5*(0.8995-1))*0.3**0.65),
    ('case_555_0.3.out',0.06*0.2071*0.3**0.65),
    ('case_670_0.3.out',0.06*0.5930*0.3**0.65),
]
for name,expected in cases:
    got=value(p/name,'a_pig')
    assert math.isclose(got,expected,rel_tol=3e-9,abs_tol=5e-11),(name,got,expected)
PY
then ok runtime_absorption_values; else bad runtime_absorption_values; fi

printf 'CCRR_CHL_V119_SMOKE pass=%d fail=%d\n' "$pass" "$fail"
[[ $fail -eq 0 ]]
