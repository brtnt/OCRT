#!/bin/bash
# 비트 동일성 배터리: legacy 제거/앵커 수정이 IPSS·비-PSSA 경로를 건드리지 않는지 확인.
BIN="$1"; OUT="$2"; mkdir -p "$OUT"
V=/root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C
cd "$V"
GAS="--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0"
r(){ tag="$1"; shift; env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 "$BIN" "$@" 2>/dev/null | tail -1 > "$OUT/$tag.txt"; }
# 1) PSSA 미사용 경로 (절대 불변이어야 함)
r a01 --wavelength 555 --sza 60 --vza 40 --raa 90 --surface black --n-layers 60 $GAS
r a02 --wavelength 412 --sza 30 --vza 55 --raa 90 --surface black_fresnel_ocean --wind-speed 5 --n-layers 60
r a03 --wavelength 865 --sza 70 --vza 20 --raa 120 --surface black_fresnel_ocean --wind-speed 7 --mie inputs/M80C.mie --aod-555 0.15 --n-layers 60
r a04 --wavelength 555 --sza 45 --vza 30 --raa 60 --surface ocean --wind-speed 5 --water-model ocrt --ocrt-chl 1 --ocrt-tsm 2 --ocrt-adom440 0.05 --n-layers 40
r a05 --wavelength 490 --sza 55 --vza 50 --raa 150 --surface ocean --wind-speed 3 --water-model ccrr --ccrr-chl 0.5 --ccrr-tsm 1 --ccrr-adom440 0.02 --n-layers 40
r a06 --wavelength 660 --sza 20 --vza 0 --raa 0 --surface flat --n-layers 60 $GAS
# 2) LUT 격자 (PSSA 미사용)
env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 "$BIN" --wavelength 555 --sza 60 --surface black --n-layers 60 $GAS \
  --lut-vza-max 70 --lut-vza-step 10 --lut-raa-step 30 --output-full-grid "$OUT/a07.csv" >/dev/null 2>&1
# 3) IPSS 경로
r b01 --wavelength 555 --sza 60 --vza 55 --raa 90 --surface black --n-layers 60 --pssa $GAS
r b02 --wavelength 412 --sza 80 --vza 55 --raa 90 --surface black_fresnel_ocean --wind-speed 5 --mie inputs/M80C.mie --aod-555 0.1 --n-layers 400 --pssa
env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 "$BIN" --wavelength 555 --sza 70 --surface black --n-layers 60 --pssa $GAS \
  --lut-vza-max 70 --lut-vza-step 10 --lut-raa-step 30 --output-full-grid "$OUT/b03.csv" >/dev/null 2>&1
