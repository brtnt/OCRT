#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export TERM=xterm
COMMON=(
  --surface ocean --water-model ocrt
  --ocrt-chl 0.3 --ocrt-tsm 0 --ocrt-adom440 0
  --wind-speed 3 --sza 40 --vza 30 --raa 90 --wavelength 443
  --pressure 0
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0
  --decouple-sunglint
)

expect_group_reject() {
  local name=$1
  set +e
  "$BIN" "${COMMON[@]}" --ocrt-phyto-group "$name" >"$TMP/$name.out" 2>"$TMP/$name.err"
  local rc=$?
  set -e
  [[ $rc -eq 2 ]] || { echo "expected rc=2 for $name, got $rc" >&2; return 1; }
  grep -F -- "--ocrt-phyto-group is not supported yet" "$TMP/$name.err" >/dev/null
  grep -F -- "requested '$name'" "$TMP/$name.err" >/dev/null
}
expect_group_reject eap_synechococcus
expect_group_reject eap_diatoms_pennate
expect_group_reject micro

OCRT_DUMP_IOP=1 "$BIN" "${COMMON[@]}" >"$TMP/chl.out" 2>"$TMP/chl.err"
python - "$TMP/chl.err" <<'PY'
import re, sys
s=open(sys.argv[1], encoding='utf-8', errors='replace').read()
m=re.search(r'OCRT_IOP .*?phyto\[a=([0-9eE+.-]+) b=([0-9eE+.-]+) bb=([0-9eE+.-]+)\] '
            r'detritus\[a=([0-9eE+.-]+) b=([0-9eE+.-]+) bb=([0-9eE+.-]+)\]', s)
if not m:
    raise SystemExit('OCRT_IOP line not found')
a,b,bb,ad,bd,bbd=map(float,m.groups())
assert a > 0.0, (a,b,bb)
assert b == 0.0 and bb == 0.0, (a,b,bb)
assert bd > 0.0 and bbd > 0.0, (ad,bd,bbd)
if 'OCRT_PHASE_COMPONENT name=phyto' in s:
    raise SystemExit('disabled phyto phase entered the RT mixture')
if 'OCRT_PHASE_COMPONENT name=detritus' not in s:
    raise SystemExit('detritus-only particle phase not observed')
print('PASS: Chl absorption retained; EAP b/bb zero; detritus-only phase')
PY

OCRT_DUMP_IOP=1 "$BIN" "${COMMON[@]}" --ocrt-mie-truncation \
  >"$TMP/trunc.out" 2>"$TMP/trunc.err"
grep -F 'Rrs0plus_I=' "$TMP/trunc.out" >/dev/null
python - "$TMP/trunc.err" <<'PYTRUNC'
import re, sys
s=open(sys.argv[1], encoding='utf-8', errors='replace').read()
mc=re.search(r'OCRT_TRUNC_COMPONENT name=detritus A=([0-9eE+.-]+) f=([0-9eE+.-]+) '
             r'b_raw=([0-9eE+.-]+) b_eff=([0-9eE+.-]+) '
             r'bb_b_residual=([0-9eE+.-]+)', s)
mp=re.search(r'OCRT_PHASE_COMPONENT name=detritus .*?b=([0-9eE+.-]+) '
             r'bb_b_iop=([0-9eE+.-]+) bb_b_phase=([0-9eE+.-]+) '
             r'A=([0-9eE+.-]+)', s)
mi=re.search(r'OCRT_IOP .*?phyto\[a=([0-9eE+.-]+) b=([0-9eE+.-]+) bb=([0-9eE+.-]+)\] '
             r'detritus\[a=([0-9eE+.-]+) b=([0-9eE+.-]+) bb=([0-9eE+.-]+)\]', s)
if not (mc and mp and mi):
    raise SystemExit('missing constituent truncation diagnostics')
A,f,braw,beff,bbres=map(float,mc.groups())
bphase,bb_iop_ratio,bb_phase_ratio,Aphase=map(float,mp.groups())
a,b,bb,ad,bd,bbd=map(float,mi.groups())
assert 0.1 < A < 2.0 and abs(f - A/2.0) < 1e-10, (A,f)
assert abs(beff - braw*(1.0-f)) <= 2e-10*max(1.0,abs(beff)), (braw,beff,f)
assert b == 0.0 and bb == 0.0, (a,b,bb)
assert bd > 0.0 and bbd > 0.0, (ad,bd,bbd)
assert abs(bd-beff) <= 2e-10*max(1.0,abs(beff)), (bd,beff)
assert abs(bb_iop_ratio-bb_phase_ratio) <= 2e-10*max(1.0,abs(bb_iop_ratio)), (bb_iop_ratio,bb_phase_ratio)
assert abs(bbres-bb_phase_ratio) <= 2e-10*max(1.0,abs(bbres)), (bbres,bb_phase_ratio)
assert abs(A-Aphase) <= 2e-10, (A,Aphase)
print('PASS: constituent detritus direct truncation is active and bb/b is consistent')
PYTRUNC

"$BIN" --surface ocean --water-model iop \
  --iop-a 0.02 --iop-b 0.2 --iop-bb 0.002 \
  --iop-mie-phase "$ROOT/inputs/water_iop/Detritus_Stramski2001.mie" \
  --iop-mie-truncation --wind-speed 3 --sza 40 --vza 30 --raa 90 \
  --wavelength 443 --pressure 0 \
  --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
  --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 \
  --decouple-sunglint >"$TMP/iop.out" 2>"$TMP/iop.err"
grep -F 'Rrs0plus_I=' "$TMP/iop.out" >/dev/null

echo 'PASS: EAP species scattering disabled and fixed-bulk truncation preserved'
