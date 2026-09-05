#!/usr/bin/env bash
# OCRT v1.18 OCRT-prefixed organic Chl smoke and phase-cache invariants.
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/.." && pwd)
BIN=${1:-$ROOT/build/v2_solver_vk_v1.18}
DATA=$ROOT/inputs/water_iop
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1

pass=0
fail=0
ok(){ printf 'PASS %s\n' "$1"; pass=$((pass+1)); }
bad(){ printf 'FAIL %s\n' "$1"; fail=$((fail+1)); }

[[ -x "$BIN" ]] || { echo "missing executable: $BIN" >&2; exit 2; }

# 1. Standalone module: strict compile + model anchors.
cat >"$TMP/test_org.c" <<'C'
#include <stdio.h>
#include "rt_iop_organic.h"
int main(int argc, char **argv) {
    if (argc != 2 || rt_iop_organic_init(argv[1]) != 0) return 2;
    printf("huot %.17g fph %.17g\n", rt_iop_huot_bbp(550.0,1.0), rt_iop_fph_fraction(1.0));
    for (int g=0; g<3; ++g) {
        rt_iop_t x={0};
        if (rt_iop_eap_phyto_eval(443.0,1.0,(organic_phyto_group_t)g,&x) != 0) return 3;
        printf("%s %.17g %.17g %.17g %.17g\n",
               rt_iop_organic_group_name((organic_phyto_group_t)g),
               x.a,x.b,x.bb,x.bb/x.b);
    }
    rt_iop_t d={0};
    if (rt_iop_detritus_eval(443.0,1.0,0.01,0.0109,&d) != 0) return 4;
    printf("det %.17g %.17g %.17g %.17g\n",d.a,d.b,d.bb,d.bb/d.b);
    return 0;
}
C
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -Werror -I"$ROOT/src" \
  "$TMP/test_org.c" "$ROOT/src/rt_iop_organic.c" \
  "$ROOT/src/shared/mie_io.c" "$ROOT/src/shared/numerics.c" \
  "$ROOT/src/shared/io_utils.c" -lm -o "$TMP/test_org"
"$TMP/test_org" "$DATA" >"$TMP/unit.out"
if python3 - "$TMP/unit.out" <<'PY'
import math,sys
rows=[x.split() for x in open(sys.argv[1],encoding='utf-8')]
assert math.isclose(float(rows[0][1]),2.267e-3,rel_tol=0,abs_tol=1e-15)
assert math.isclose(float(rows[0][3]),0.035,rel_tol=0,abs_tol=1e-15)
d={r[0]:list(map(float,r[1:])) for r in rows[1:]}
assert math.isclose(d['micro'][0],0.034612,rel_tol=0,abs_tol=1e-14)
assert math.isclose(d['det'][0],0.01*math.exp(-0.0109*3),rel_tol=2e-15,abs_tol=2e-15)
assert d['pico'][0] > d['nano'][0] > d['micro'][0] > 0
assert d['pico'][3] > d['nano'][3] > d['micro'][3] > 0
assert 0 < d['det'][3] < 0.5
PY
then ok 'organic module strict build and model anchors'; else bad 'organic module strict build and model anchors'; fi

GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
GEOM=(--surface ocean --wind-speed 0 --sza 30 --vza 30 --raa 90 \
      --wavelength 443 --pressure 0 "${GAS[@]}")
FAST=(--debug-water-max-orders 4 --n-mu-water 12 --water-m-max 4)
run_fast(){
  OCRT_DEBUG=1 OCRT_ADVANCED=1 "$BIN" "$@" "${FAST[@]}"
}

# 2. Canonical contract: Chl, TSM and aDOM440 are explicit; zero is valid.
set +e
"$BIN" --water-model ocrt --ocrt-chl 0.1 --ocrt-adom440 0 "${GEOM[@]}" \
  >"$TMP/missing.out" 2>"$TMP/missing.err"
rc=$?
set -e
if [[ $rc -eq 2 ]] && grep -q 'requires explicit --ocrt-chl, --ocrt-tsm and --ocrt-adom440' "$TMP/missing.err"; then
  ok 'canonical required-field gate'
else bad 'canonical required-field gate'; fi

if run_fast --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 "${GEOM[@]}" \
     >"$TMP/zero.out" 2>"$TMP/zero.err"; then
  ok 'explicit-zero constituent input'
else bad 'explicit-zero constituent input'; fi

# 3. Three EAP groups are operational and distinct.
for g in pico nano micro; do
  OCRT_DUMP_IOP=1 run_fast --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 \
    --ocrt-phyto-group "$g" "${GEOM[@]}" >"$TMP/$g.out" 2>"$TMP/$g.err"
done
if python3 - "$TMP/pico.out" "$TMP/nano.out" "$TMP/micro.out" <<'PY'
import re,sys
vals=[]
for p in sys.argv[1:]:
    t=open(p,encoding='utf-8').read().splitlines()[-1]
    m=re.search(r'a_pig=([0-9.eE+\-]+)',t); assert m
    vals.append(float(m.group(1)))
assert vals[0] > vals[1] > vals[2] > 0, vals
PY
then ok 'pico/nano/micro group selection'; else bad 'pico/nano/micro group selection'; fi

# 4. Default detritus slope is exactly the documented 0.0109.
OCRT_DUMP_IOP=1 run_fast --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 \
  --ocrt-detritus-a440 0.01 "${GEOM[@]}" >"$TMP/sdef.out" 2>"$TMP/sdef.err"
OCRT_DUMP_IOP=1 run_fast --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 \
  --ocrt-detritus-a440 0.01 --ocrt-detritus-slope 0.0109 "${GEOM[@]}" \
  >"$TMP/sexp.out" 2>"$TMP/sexp.err"
if cmp -s "$TMP/sdef.out" "$TMP/sexp.out"; then ok 'detritus slope default bit identity';
else bad 'detritus slope default bit identity'; fi

# 5. Organic-data wavelength boundary.
set +e
"$BIN" --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 \
  --surface ocean --wind-speed 0 --sza 30 --vza 30 --raa 90 \
  --wavelength 900 --pressure 0 "${GAS[@]}" >"$TMP/range.out" 2>"$TMP/range.err"
rc=$?
set -e
if [[ $rc -eq 2 ]] && grep -q 'validated only over 350-850 nm' "$TMP/range.err"; then
  ok 'organic wavelength range gate'
else bad 'organic wavelength range gate'; fi

# 6. Chl + TSM + aDOM: three vector particulate components and exact IOP/phase bb/b.
OCRT_DUMP_IOP=1 run_fast --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0.01 --ocrt-adom440 0.02 \
  --ocrt-phyto-group micro --ocrt-tsm-species yellow_clay "${GEOM[@]}" \
  >"$TMP/mixed.out" 2>"$TMP/mixed.err"
if python3 - "$TMP/mixed.err" <<'PY'
import re,sys,math
text=open(sys.argv[1],encoding='utf-8').read()
assert 'OCRT_PHASE_MIX components=3' in text
rows=re.findall(r'OCRT_PHASE_COMPONENT name=(\w+).*?bb_b_iop=([0-9.eE+\-]+) bb_b_phase=([0-9.eE+\-]+)',text)
assert {r[0] for r in rows}=={'phyto','detritus','tsm'}, rows
for name,a,b in rows:
    assert float(a)==float(b),(name,a,b)
PY
then ok 'three-component vector phase and exact bb/b identity'; else bad 'three-component vector phase and exact bb/b identity'; fi

# 7. Phase/moment work is independent of SOS order count.
for n in 2 12; do
  OCRT_DUMP_WATER_COMPONENT_PHASE_CACHE=1 OCRT_DEBUG=1 OCRT_ADVANCED=1 \
    "$BIN" --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 "${GEOM[@]}" \
    --debug-water-max-orders "$n" --n-mu-water 12 --water-m-max 4 >"$TMP/order$n.out" 2>"$TMP/order$n.err"
  grep -E '^(WCPHASE build|WCRATIO build|WCMIX build)' "$TMP/order$n.err" >"$TMP/order$n.trace"
  [[ $(grep -c '^WCPHASE build' "$TMP/order$n.trace") -eq 2 ]]
  [[ $(grep -c '^WCRATIO build' "$TMP/order$n.trace") -eq 2 ]]
  [[ $(grep -c '^WCMIX build' "$TMP/order$n.trace") -eq 1 ]]
done
if cmp -s "$TMP/order2.trace" "$TMP/order12.trace"; then ok 'phase preparation count independent of SOS orders';
else bad 'phase preparation count independent of SOS orders'; fi

# 8. Cache only changes reuse, not numerical output.
OCRT_DEBUG=1 OCRT_ADVANCED=1 "$BIN" --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 \
  "${GEOM[@]}" "${FAST[@]}" >"$TMP/cache_on.out" 2>"$TMP/cache_on.err"
OCRT_DISABLE_WATER_COMPONENT_PHASE_CACHE=1 OCRT_DEBUG=1 OCRT_ADVANCED=1 \
  "$BIN" --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 \
  "${GEOM[@]}" "${FAST[@]}" >"$TMP/cache_off.out" 2>"$TMP/cache_off.err"
if cmp -s "$TMP/cache_on.out" "$TMP/cache_off.out"; then ok 'cache enabled/disabled stdout byte identity';
else bad 'cache enabled/disabled stdout byte identity'; fi

# 9. P12 and P33 are active in I/Q/U while P11-derived IOP remains unchanged.
mkdir -p "$TMP/org_orig" "$TMP/org_p12zero" "$TMP/org_p33half"
cp "$DATA"/{pico_Synechococcus_EAP.mie,nano_Haptophytes_EAP.mie,Diatoms_centric_EAP.mie,Detritus_Stramski2001.mie} "$TMP/org_orig/"
cp "$TMP/org_orig"/* "$TMP/org_p12zero/"
cp "$TMP/org_orig"/* "$TMP/org_p33half/"
python3 - "$TMP/org_p12zero" "$TMP/org_p33half" <<'PY'
from pathlib import Path
import sys

def transform(path, block, factor):
    lines=path.read_text().splitlines()
    active=False; header=False; out=[]
    marker=f'Phase Function ({block}'
    for line in lines:
        if line.startswith('Phase Function ('):
            active=line.startswith(marker); header=False
            out.append(line); continue
        if active and line.lstrip().startswith('TETA'):
            header=True; out.append(line); continue
        if active and header and line.strip():
            p=line.split()
            try: ang=float(p[0])
            except Exception:
                out.append(line); continue
            vals=[float(x)*factor for x in p[1:]]
            out.append(f'{ang:8.2f}' + ''.join(f'  {v:+.7E}' for v in vals))
        else:
            out.append(line)
    path.write_text('\n'.join(out)+'\n')

for name in ['Diatoms_centric_EAP.mie','Detritus_Stramski2001.mie']:
    transform(Path(sys.argv[1])/name,'P12',0.0)
    transform(Path(sys.argv[2])/name,'P33',0.5)
PY
for variant in orig p12zero p33half; do
  dir="$TMP/org_$variant"
  OCRT_ORGANIC_DIR="$dir" OCRT_DUMP_IOP=1 run_fast \
    --water-model ocrt --ocrt-chl 0.1 --ocrt-tsm 0 --ocrt-adom440 0 "${GEOM[@]}" \
    >"$TMP/$variant.out" 2>"$TMP/$variant.err"
  grep '^OCRT_IOP ' "$TMP/$variant.err" >"$TMP/$variant.iop"
done
if python3 - "$TMP/orig.out" "$TMP/p12zero.out" "$TMP/p33half.out" \
                    "$TMP/orig.iop" "$TMP/p12zero.iop" "$TMP/p33half.iop" <<'PY'
import sys

def stokes(p):
    row=open(p,encoding='utf-8').read().strip().splitlines()[-1].split()
    return tuple(map(float,row[:3]))
o,p12,p33=map(stokes,sys.argv[1:4])
iops=[open(p,'rb').read() for p in sys.argv[4:7]]
assert iops[0]==iops[1]==iops[2]
assert abs(o[1]-p12[1])+abs(o[2]-p12[2]) > 1e-8,(o,p12)
assert abs(o[1]-p33[1])+abs(o[2]-p33[2]) > 1e-8,(o,p33)
PY
then ok 'P12/P33 affect IQU with invariant P11-derived IOP'; else bad 'P12/P33 affect IQU with invariant P11-derived IOP'; fi

printf 'ORGANIC_CHL_SMOKE pass=%d fail=%d\n' "$pass" "$fail"
[[ $fail -eq 0 ]]
