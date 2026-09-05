#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
BIN=${OCRT_BIN:-./build/ocrt_v1.2}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
GAS=(--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 \
     --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)
COMMON=(--surface ocean --wind-speed 3 --sza 30 --wavelength 555 \
 --pressure 1013.25 --aod 0 --water-model iop \
 --iop-a 0.1 --iop-b 0.5 --iop-bb 0.01 \
 --iop-mie-phase inputs/tsm_ahn/Red_clay_AHN.mie \
 --iop-mie-moment-mode gauss --iop-mie-moment-nmu 400 \
 --n-mu-water 8 --water-m-max 4 --lut-vza-step 30 --lut-vza-max 60 \
 --lut-raa-step 90 "${GAS[@]}")
ENV=(OMP_NUM_THREADS=1 OCRT_ADVANCED=1 OCRT_WATER_PARTICLE_KERNEL=direct \
     OCRT_WATER_VALUE_NPHI=128)
env "${ENV[@]}" OCRT_MIE_MODEL_CACHE_TRACE=1 OCRT_VALUE_KERNEL_CACHE_TRACE=1 \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/cache.csv" \
  >"$TMP/cache.out" 2>"$TMP/cache.err"
env "${ENV[@]}" OCRT_NO_PHASE_CACHE=1 OCRT_DISABLE_MIE_MODEL_CACHE=1 \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/no_cache.csv" \
  >"$TMP/no_cache.out" 2>"$TMP/no_cache.err"

grep -E '\[MIE-MODEL-CACHE\] load .*loads=1 hits=0' "$TMP/cache.err" >/dev/null
grep -E '\[VALUE-POL-CACHE\].*entry_full_builds=1.*entry_beam_hits=[1-9][0-9]*' "$TMP/cache.err" >/dev/null
# One all-view case with one active Mie phase must retain exactly one complete
# Fourier-kernel entry.  A larger count means the beam slot leaked into the
# beam-independent direction-grid hash or another false-miss key regressed.
if grep -E '\[VALUE-POL-CACHE\].*entries=([2-9]|[1-9][0-9]+)([^0-9]|$)' "$TMP/cache.err" >/dev/null; then
  echo 'ERROR: one all-view case created more than one direct kernel entry' >&2
  cat "$TMP/cache.err" >&2
  exit 1
fi
grep -E '\[VALUE-POL-CACHE\].*entries=1([^0-9]|$)' "$TMP/cache.err" >/dev/null
grep -E '\[native-coupled-lut\].*water_cold=1.*replay=0' "$TMP/cache.err" >/dev/null

python - "$TMP/cache.csv" "$TMP/no_cache.csv" <<'PY'
import csv, math, sys
with open(sys.argv[1], newline='') as f: a=list(csv.DictReader(f))
with open(sys.argv[2], newline='') as f: b=list(csv.DictReader(f))
if len(a)!=len(b) or (a and a[0].keys()!=b[0].keys()):
    raise SystemExit('cache/no-cache structure mismatch')
max_abs=0.0
for ra,rb in zip(a,b):
    for k in ra:
        try:
            x=float(ra[k]); y=float(rb[k])
        except ValueError:
            if ra[k]!=rb[k]: raise SystemExit(f'text mismatch {k}')
            continue
        if math.isnan(x) and math.isnan(y): continue
        d=abs(x-y); max_abs=max(max_abs,d)
        if d!=0.0: raise SystemExit(f'cache changed output {k}: {x} vs {y}')
print(f'cache/no-cache max_abs={max_abs:.3e}')
PY

set +e
env "${ENV[@]}" OCRT_NATIVE_COUPLED_LUT_OFF=1 \
  "$BIN" "${COMMON[@]}" --output-full-grid "$TMP/prohibited.csv" \
  >"$TMP/prohibited.out" 2>"$TMP/prohibited.err"
rc=$?
set -e
if [[ $rc -eq 0 ]]; then
  echo 'ERROR: per-cell replay was not prohibited' >&2
  exit 1
fi
grep -F 'Per-cell VZA/RAA solver replay is prohibited' "$TMP/prohibited.err" >/dev/null

echo 'PASS: worker-private Mie model cache + full-direction/beam kernel cache; output invariant; cell replay prohibited'
