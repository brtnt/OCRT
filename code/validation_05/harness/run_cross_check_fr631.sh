#!/bin/sh
# OCRT-OSOAA 자동 교차 하니스: 패키지 루트에서 실행. OCRT 3런 후 저장 OSOAA 참조와 자동 대조.
set -e
PKG="$(cd "$(dirname "$0")/../.." && pwd)"
B="$PKG/01_OCRT_C"
[ -x "$B/build/ocrt" ] || (cd "$B" && make all)
TAU=16.6358629872763
run(){ d="$1"; shift; mkdir -p "$d"; (cd "$B" && env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 OCRT_DEBUG=1 ./build/ocrt "$@" --output-full-grid "$d/fullgrid.csv" >/dev/null 2>&1); }
O="$PKG/05_VALIDATION/harness/out"; rm -rf "$O"; mkdir -p "$O"
run "$O/pure" --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0 --ocrt-adom440 0 --surface ocean --wind-speed 3 --n-water 1.34 --n-mu-water 48 --wavelength 555 --sza 40 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 --debug-water-tau-max-target $TAU --debug-water-max-tau $TAU --lut-vza-max 60 --lut-vza-step 5 --lut-raa-step 15
run "$O/tsm"  --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0.5 --ocrt-tsm-species red_clay --ocrt-adom440 0 --surface ocean --wind-speed 3 --n-water 1.34 --n-mu-water 48 --wavelength 555 --sza 40 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 --debug-water-tau-max-target $TAU --debug-water-max-tau $TAU --lut-vza-max 60 --lut-vza-step 5 --lut-raa-step 15
run "$O/bfo"  --surface black_fresnel_ocean --wind-speed 3 --decouple-sunglint --n-water 1.34 --wavelength 555 --sza 40 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 --lut-vza-max 60 --lut-vza-step 5 --lut-raa-step 15
python3 "$PKG/05_VALIDATION/harness/compare_osoaa_fr631.py" "$PKG" "$O/pure/fullgrid.csv" "$O/tsm/fullgrid.csv" "$O/bfo/fullgrid.csv"
