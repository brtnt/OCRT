#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-$ROOT/build/ocrt}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
run_case() {
  local vza=$1 wind=$2 out=$3 nmu=${4:-48} maxord=${5:-500}
  (cd "$ROOT" && OCRT_ADVANCED=1 OCRT_DEBUG=1 OMP_NUM_THREADS=1 "$BIN" \
    --surface ocean --wind-speed "$wind" --sza 30 --vza "$vza" --raa 90 \
    --wavelength 555 --pressure 0 --aod 0 --water-model ocrt \
    --ocrt-chl 0 --ocrt-tsm 5 --ocrt-adom440 0 \
    --n-mu-water "$nmu" --water-m-max 4 --max-orders "$maxord") > "$out" 2> "$out.err"
}
run_case 0 3 "$TMP/nadir.out"
run_case 60 3 "$TMP/rough.out"
run_case 60 0 "$TMP/flat.out" 12 80
python3 - "$TMP/nadir.out" "$TMP/rough.out" "$TMP/flat.out" <<'PY'
import math,re,sys

def parse(path):
    s=open(path).read()
    def g(k):
        m=re.search(r'\b'+re.escape(k)+r'=([-+0-9.eE]+)',s)
        if not m: raise SystemExit(f'missing {k}: {path}')
        return float(m.group(1))
    return {k:g(k) for k in ('rrs0minus_I','rrs0minus_Q','rrs0minus_U','Rrs0plus_I','Rrs0plus_Q','Rrs0plus_U')}
n,r,f=map(parse,sys.argv[1:])
# Nadir and even-m invariants.
assert abs(n['rrs0minus_U']) < 1e-15, n
assert abs(n['rrs0minus_I'] - 8.400221e-2) < 2e-7, n
assert abs(n['rrs0minus_Q'] - 1.420403e-3) < 2e-8, n
# Rough canonical with the corrected 2026-07-28 pure-water table.
# OSOAA was rerun with the same table: DoLP=0.0873199661.  Stage-2
# public water U uses the positive-sine output convention, so RAA=90 is
# negative in OCRT; DoLP is convention-invariant and must agree within 2%.
dolp=math.hypot(r['rrs0minus_Q'],r['rrs0minus_U'])/abs(r['rrs0minus_I'])
assert abs(dolp/0.0873199660635599-1.0) < 0.02, dolp
assert r['rrs0minus_U']/r['rrs0minus_I'] < -0.084, r
# Flat path is the same azimuth-preserving physics; low-resolution smoke guard.
assert f['rrs0minus_U']/f['rrs0minus_I'] < -0.08, f
print(f'PASS water internal-reflection msign: rough_DoLP={dolp:.9f}')
PY
# The generic helper is intentionally retained, but this water solver file must
# not apply it to internal-reflection boundary conditions.
if grep -q 'const double msign = rt_fourier_pi_shift_sign(m)' "$ROOT/src/rt_solver.c"; then
  echo 'FAIL: stale (-1)^m water internal-reflection sign remains' >&2
  exit 1
fi
