#!/usr/bin/env bash
set -euo pipefail

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/.." && pwd)
BIN=${1:-$ROOT/build/v2_solver_vk_v1.18}
DATA=$ROOT/inputs/tsm_ahn
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

pass=0
fail=0
ok() { printf 'PASS %s\n' "$1"; pass=$((pass+1)); }
bad() { printf 'FAIL %s\n' "$1"; fail=$((fail+1)); }

if [ ! -x "$BIN" ]; then
  echo "missing executable: $BIN" >&2
  exit 2
fi

cat > "$TMP/test_ahn.c" <<'C'
#include <stdio.h>
#include "rt_iop_ahn_mineral.h"
int main(int argc, char **argv) {
    if (argc != 2 || rt_iop_ahn_mineral_init(argv[1]) != 0) return 2;
    for (int s = 0; s < AHN_SPECIES_COUNT; ++s) {
        rt_iop_t out = {0};
        if (rt_iop_ahn_mineral_eval(443.0, 2.0, (ahn_species_t)s, &out, NULL) != 0) return 3;
        printf("%s %.12g %.12g %.12g\n", rt_iop_ahn_species_name((ahn_species_t)s),
               out.a, out.b, out.bb/out.b);
    }
    return 0;
}
C

gcc -std=c11 -Wall -Wextra -Wpedantic -Werror -I"$ROOT/src" \
  "$TMP/test_ahn.c" "$ROOT/src/rt_iop_ahn_mineral.c" -lm -o "$TMP/test_ahn"
"$TMP/test_ahn" "$DATA" > "$TMP/unit.txt"
if python - "$TMP/unit.txt" <<'PY'
import math,sys
# FR631 theta-linear integration of the validated legacy AHN phase.
# These differ slightly from the old 361-node trapezoid values while the
# underlying values at every historical 0.5-degree node remain exact.
exp={
'red_clay':(0.15594,1.72062,0.02247242315369668),
'brown_earth':(0.20748,1.56368,0.02016791165937695),
'yellow_clay':(0.07834,1.77760,0.02012441125324285),
'calcareous_sand':(0.05120,2.07360,0.02288574528582990),
}
rows={}
for line in open(sys.argv[1]):
    p=line.split(); rows[p[0]]=tuple(map(float,p[1:]))
assert rows.keys()==exp.keys()
for k,want in exp.items():
    got=rows[k]
    for a,b in zip(got,want):
        assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-12),(k,got,want)
PY
then ok 'four-species 443-nm optical coefficients'; else bad 'four-species 443-nm optical coefficients'; fi

mkdir "$TMP/mie_only"
cp "$DATA"/*_AHN.mie "$TMP/mie_only/"
"$TMP/test_ahn" "$TMP/mie_only" > "$TMP/fallback.txt"
if diff -u "$TMP/unit.txt" "$TMP/fallback.txt" >/dev/null; then
  ok 'Mie-only a*/b* fallback'
else
  bad 'Mie-only a*/b* fallback'
fi

COMMON=(--surface ocean --wind-speed 0 --water-model ocrt \
        --ocrt-chl 0 --ocrt-tsm 0.01 --ocrt-adom440 0 \
        --sza 30 --vza 0 --raa 90 --wavelength 443 --pressure 0)
"$BIN" "${COMMON[@]}" --ocrt-tsm-species red_clay >"$TMP/red.out" 2>"$TMP/red.err"
if grep -q 'Red_clay_AHN.mie' "$TMP/red.err" && grep -q 'a_min=' "$TMP/red.out"; then
  ok 'CLI automatic red-clay vector-phase routing'
else
  bad 'CLI automatic red-clay vector-phase routing'
fi

"$BIN" "${COMMON[@]}" >"$TMP/default.out" 2>"$TMP/default.err"
if cmp -s "$TMP/red.out" "$TMP/default.out"; then
  ok 'default species equals explicit red_clay'
else
  bad 'default species equals explicit red_clay'
fi

set +e
"$BIN" --surface ocean --wind-speed 0 --water-model ocrt \
       --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 \
       --ocrt-tsm-species brown_earth \
       --sza 30 --vza 0 --raa 90 --wavelength 443 --pressure 0 \
       >"$TMP/no_min.out" 2>"$TMP/no_min.err"
rc=$?
set -e
if [ "$rc" -eq 2 ] && grep -q 'requires --ocrt-tsm > 0' "$TMP/no_min.err"; then
  ok 'species without positive TSM is rejected'
else
  bad 'species without positive TSM is rejected'
fi

set +e
"$BIN" --surface ocean --wind-speed 0 --water-model ocrt \
       --ocrt-chl 0 --ocrt-tsm 0.01 --ocrt-adom440 0 \
       --ocrt-tsm-species brown_earth \
       --iop-mie-phase "$DATA/Red_clay_AHN.mie" \
       --sza 30 --vza 0 --raa 90 --wavelength 443 --pressure 0 \
       >"$TMP/conflict.out" 2>"$TMP/conflict.err"
rc=$?
set -e
if [ "$rc" -eq 2 ] && grep -Eq 'cannot be combined|another water branch' "$TMP/conflict.err"; then
  ok 'explicit phase override conflict is rejected'
else
  bad 'explicit phase override conflict is rejected'
fi

printf 'TSM_AHN_SMOKE pass=%d fail=%d\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
