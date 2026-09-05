#!/bin/bash
# IPSS 게이트 회귀 실행기.
#   1) 기하/단일산란/스케일링/앵커 게이트 (G-I0 ~ G-I3)
#   2) 삭제된 --pssa-mode 재진입 차단
#   3) --pssa on/off 가 실제로 다른 값을 내는지 (보정이 살아있는지)
#   4) 비-PSSA 경로 비트 불변
#   5) 결합 해양 + --pssa fail-loud
# 사용: bash scripts/run_ipss_gates.sh   (OCRT_C 루트에서)
set -u
cd "$(dirname "$0")/.."
fail=0
GAS="--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0"

echo "== 1) 기하/단일산란/스케일링/앵커 게이트 =="
gcc -std=c11 -O2 -Isrc tests/test_ipss_gates.c src/rt_ipss.c -o build/test_ipss -lm || exit 1
./build/test_ipss || fail=1

echo "== 2) 삭제된 --pssa-mode 재진입 차단 =="
if env OMP_NUM_THREADS=1 ./build/ocrt --sza 40 --wavelength 555 --surface black \
      --vza 20 --raa 90 --pssa --pssa-mode legacy $GAS >/dev/null 2>&1; then
  echo "[FAIL] --pssa-mode legacy 가 통과됨 (v1.11 에서 삭제되어야 함)"; fail=1
else echo "[PASS] --pssa-mode 거부 (구 legacy 경로 도달 불가)"; fi
# 코드 심볼만 검사한다(제거 경위를 적어둔 주석은 남겨 둔다).
if grep -rqsE 'rt_pssa_apply *\(|RT_PSSA_MODE_|#include *"rt_pssa\.h"|->pssa_active|\.pssa_active|pssa_xi_dn|pssa_alpha|pssa_beta|pssa_mode' src/ tests/ tools/ \
   || [ -e src/rt_pssa.c ] || [ -e src/rt_pssa.h ]; then
  echo "[FAIL] 소스 트리에 legacy PSSA 심볼/파일이 남아 있음"; fail=1
else echo "[PASS] 소스 트리에 legacy PSSA 심볼/파일 없음"; fi

echo "== 3) --pssa on/off 가 서로 달라야 (보정 활성 확인) =="
SM=(--sza 75 --wavelength 555 --pressure 1013.25 --surface black_fresnel_ocean
    --wind-speed 5 --decouple-sunglint --lut-vza-step 5 --lut-vza-max 85
    --lut-raa-step 5)
env OMP_NUM_THREADS=1 ./build/ocrt "${SM[@]}" \
    --output-full-grid /tmp/ipss_gate_off.csv 2>/dev/null
env OMP_NUM_THREADS=1 ./build/ocrt "${SM[@]}" --pssa \
    --output-full-grid /tmp/ipss_gate_on.csv 2>/dev/null
if cmp -s /tmp/ipss_gate_off.csv /tmp/ipss_gate_on.csv; then
  echo "[FAIL] --pssa 가 결과를 바꾸지 않음"; fail=1
else echo "[PASS] --pssa on != off"; fi

echo "== 4) 비-PSSA 경로 불변 =="
# 기준 CSV 는 legacy 삭제 직전 바이너리(동일 머신)와 비트 동일함을 확인한 산출이다.
# 같은 -march 로 빌드하면 SHA 까지 일치하고, 다른 -march 에서는 -fassociative-math
# 재결합 순서가 달라져 ULP 수준 차이가 정상이므로 수치 허용오차로 판정한다.
REF=tests/fixtures/ipss_nonpssa_reference.csv
got=$(sha256sum /tmp/ipss_gate_off.csv | cut -c1-16)
want="4ced4f18d50e4ac9"
if [ "$got" = "$want" ]; then
  echo "[PASS] non-PSSA 비트 동일 (SHA $got)"
else
  worst=$(python3 - "$REF" /tmp/ipss_gate_off.csv <<'PY'
import csv, sys
a = list(csv.reader(open(sys.argv[1]))); b = list(csv.reader(open(sys.argv[2])))
if len(a) != len(b) or a[0] != b[0]:
    print("shape"); raise SystemExit
# 반사도 1e-12 미만은 수치적 0 이다(천저에서 Q,U 는 해석적으로 0, 계산값 ~1e-17).
# 절대차 1e-12 이내이거나 상대차 1e-12 이내면 통과로 본다.
ATOL = 1e-12
w = 0.0
for ra, rb in zip(a[1:], b[1:]):
    for x, y in zip(ra, rb):
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            continue
        ad = abs(fx - fy)
        if ad <= ATOL:
            continue
        d = ad / max(abs(fx), abs(fy), ATOL)
        if d > w: w = d
print(f"{w:.3e}")
PY
)
  if [ "$worst" = "shape" ]; then
    echo "[FAIL] non-PSSA 기준 CSV 와 형식이 다름"; fail=1
  elif python3 -c "import sys;sys.exit(0 if float('$worst')<1e-12 else 1)"; then
    echo "[PASS] non-PSSA 수치 일치 (최대 상대차 $worst; SHA $got != $want 는 -march 차이)"
  else
    echo "[FAIL] non-PSSA 최대 상대차 $worst (허용 1e-12)"; fail=1
  fi
fi

echo "== 5) 결합 해양 오진입 차단 =="
if env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 OCRT_DEBUG=1 ./build/ocrt --water-model ocrt \
      --ocrt-chl 0 --ocrt-tsm 0.5 --ocrt-tsm-species red_clay --ocrt-adom440 0 \
      --surface ocean --wind-speed 3 --n-water 1.34 --wavelength 555 --sza 40 \
      --vza 20 --raa 90 --pssa $GAS >/dev/null 2>&1; then
  echo "[FAIL] 결합 해양 + --pssa 가 통과됨 (fail-loud 이어야 함)"; fail=1
else echo "[PASS] 결합 해양 + --pssa fail-loud"; fi

echo
[ $fail -eq 0 ] && echo "IPSS REGRESSION: ALL PASS" || echo "IPSS REGRESSION: FAILURES PRESENT"
exit $fail
