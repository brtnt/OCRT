#!/usr/bin/env bash
# Verify that the prepared water phase operator is built once per Fourier mode,
# independent of the number of SOS scattering orders.
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/.." && pwd)
BIN=${1:-$ROOT/build/v2_solver_vk_v1.18}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

if [[ ! -x "$BIN" ]]; then
  echo "missing executable: $BIN" >&2
  exit 2
fi

# Apply identical synthetic I/Q/U fields through both execution modes inside one process: (a) per-call fallback packing
# and (b) the prepared per-mode immutable operator cache.  The C test compares every output double by raw bits and then
# repeats the cached call to verify arena/operator reuse.
mapfile -t TEST_SRC < <(find "$ROOT/src" -name '*.c' ! -name 'main.c' | sort)
gcc -std=c11 -O0 -g0 -fopenmp -DOCRT_FAST_KERNELS \
  -ffunction-sections -fdata-sections -I"$ROOT/src" \
  "$ROOT/tests/test_sos_phase_cache_equivalence.c" "${TEST_SRC[@]}" \
  -Wl,--gc-sections -o "$TMP/test_sos_phase_cache_equivalence" -lm \
  >"$TMP/direct-build.out" 2>"$TMP/direct-build.err"
"$TMP/test_sos_phase_cache_equivalence" >"$TMP/direct.out" 2>"$TMP/direct.err"
grep -Fq 'PASS direct-cache-equivalence' "$TMP/direct.out"
cat "$TMP/direct.out"

COMMON=(--surface ocean --wind-speed 0 --water-model ocrt
        --ocrt-chl 0 --ocrt-tsm 2 --ocrt-adom440 0
        --ocrt-tsm-species red_clay
        --sza 30 --vza 30 --raa 90 --wavelength 443 --pressure 0
        --water-m-max 8
        --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
        --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0)

for n in 2 100; do
  OMP_NUM_THREADS=1 OCRT_ADVANCED=1 OCRT_DEBUG=1 \
  OCRT_WATER_M_NO_EARLY_EXIT=1 OCRT_DUMP_SOS_PHASE_CACHE=1 \
    "$BIN" "${COMMON[@]}" --debug-water-max-orders "$n" \
    >"$TMP/$n.out" 2>"$TMP/$n.err"
  grep '^SOS_PHASE_CACHE_BUILD' "$TMP/$n.err" >"$TMP/$n.trace"
  [[ $(wc -l <"$TMP/$n.trace") -eq 9 ]]
  grep -q 'orders='"$n" "$TMP/$n.out"
done

cmp "$TMP/2.trace" "$TMP/100.trace"
printf 'TSM_PHASE_CACHE PASS builds=9 trace_sha256='
sha256sum "$TMP/2.trace" | awk '{print $1}'
