#!/usr/bin/env bash
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OCRT=$(CDPATH= cd -- "$HERE/.." && pwd)
ROOT=$(CDPATH= cd -- "$OCRT/.." && pwd)
BIN="$OCRT/build/ocrt"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
export OMP_NUM_THREADS=1
FULL_IDENTITY=0
for arg in "$@"; do
  case "$arg" in
    --full-identity) FULL_IDENTITY=1 ;;
    -h|--help)
      echo "Usage: verify_release_v1.1.sh [--full-identity]"
      exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

[[ -x "$BIN" ]] || { echo "missing executable: $BIN" >&2; exit 2; }
[[ "$($BIN --version)" == 'ocrt OCRT-v1.1-2026-07-18-KST-release' ]]
cmp "$OCRT/build/ocrt" "$OCRT/build/ocrt_v1.1"
cmp "$OCRT/build/ocrt" "$OCRT/build/v2_solver_vk"
echo 'PASS version and canonical binary copies'

"$ROOT/verify_data.sh"

cd "$OCRT"
"$HERE/smoke_cli_water_branch_v118.sh" build/ocrt >"$TMP/water.log" 2>"$TMP/water.err"
grep -q 'WATER_BRANCH_V118_SMOKE pass=14 fail=0' "$TMP/water.log"
cat "$TMP/water.log"

"$HERE/smoke_cli_ccrr_chl_v119.sh" build/ocrt >"$TMP/ccrr.log" 2>"$TMP/ccrr.err"
grep -q 'CCRR_CHL_V119_SMOKE pass=11 fail=0' "$TMP/ccrr.log"
cat "$TMP/ccrr.log"

"$HERE/smoke_cli_organic_chl.sh" build/ocrt >"$TMP/organic.log" 2>"$TMP/organic.err"
grep -q 'ORGANIC_CHL_SMOKE pass=10 fail=0' "$TMP/organic.log"
cat "$TMP/organic.log"

"$HERE/smoke_cli_tsm_ahn.sh" build/ocrt >"$TMP/tsm.log" 2>"$TMP/tsm.err"
grep -q 'TSM_AHN_SMOKE pass=6 fail=0' "$TMP/tsm.log"
cat "$TMP/tsm.log"

"$HERE/regression_tsm_phase_cache_v116.sh" build/ocrt >"$TMP/cache.log" 2>"$TMP/cache.err"
grep -q 'TSM_PHASE_CACHE PASS builds=9' "$TMP/cache.log"
cat "$TMP/cache.log"

IDENTITY_SUMMARY="$OCRT/validation/v1.1_release/numerical_identity_v1.1_vs_dev_v1.19.txt"
IDENTITY_CSV="$OCRT/validation/v1.1_release/numerical_identity_v1.1_vs_dev_v1.19.csv"
grep -q '^cases=47$' "$IDENTITY_SUMMARY"
grep -q '^exact=47$' "$IDENTITY_SUMMARY"
grep -q '^fail=0$' "$IDENTITY_SUMMARY"
expected=$(awk -F= '$1=="csv_sha256"{print $2}' "$IDENTITY_SUMMARY")
actual=$(sha256sum "$IDENTITY_CSV" | awk '{print $1}')
[[ -n "$expected" && "$expected" == "$actual" ]]
echo "PASS prequalified 47-case identity record sha256=$actual"

if [[ $FULL_IDENTITY -eq 1 ]]; then
  [[ -x "$ROOT/history/binaries/ocrt_dev_v1.19" ]] || { echo 'missing retained v1.19 binary' >&2; exit 2; }
  python3 tools/compare_v1.1_to_dev_v1.19.py \
    --baseline "$ROOT/history/binaries/ocrt_dev_v1.19" \
    --candidate "$BIN" \
    --csv "$TMP/identity.csv" --summary "$TMP/identity.txt" \
    >"$TMP/identity.log" 2>"$TMP/identity.err"
  grep -q '^fail=0$' "$TMP/identity.txt"
  cat "$TMP/identity.txt"
fi

echo 'OCRT_V1.1_RELEASE_VERIFY PASS'
