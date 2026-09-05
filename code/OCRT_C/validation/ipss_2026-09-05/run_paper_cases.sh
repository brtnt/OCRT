#!/bin/bash
# 논문(Zhai & Hu 2022 §3.1) 균질 보존 Rayleigh 검증 케이스를 OCRT 3모드로 실행.
#   대기: 0-100 km, 흑면(black), 기체흡수 OFF, tau_total = 0.25 / 1.0
#   기하: SZA 0 / 70.47 / 84.26,  VZA 0-70(5도, TOA 정의),  RAA 0/5/.../180
#   모드: pp(보정 없음) / legacy(평균할선 Chapman PSSA) / ipss(신규 기본)
set -e
V=/root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C
D=/root/ocrt/ipss_bench
cd "$V"
GAS="--gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0"
GRID="--lut-vza-max 70 --lut-vza-step 5 --lut-raa-step 5"
# 층수는 legacy PSSA 의 고SZA 경고 기준(84도 이상 >=400)을 만족시켜 두 모드를
# 같은 이산화에서 비교한다.
NL="--n-layers 400"

run () {                      # $1=tag $2=pressure $3=sza $4=mode-args
  local out="$D/case_$1.csv"
  [ -f "$out" ] && return 0
  env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 ./build/ocrt \
      --wavelength 555 --pressure "$2" --sza "$3" --surface black \
      $NL $GAS $GRID $4 --output-full-grid "$out" 2>/dev/null
}

for tau in 025 100; do
  case $tau in
    025) P=2707.819666 ;;
    100) P=10831.278663 ;;
  esac
  for sza in 0.0 70.47 84.26; do
    s=${sza/./p}
    run "t${tau}_s${s}_pp"     "$P" "$sza" ""
    run "t${tau}_s${s}_legacy" "$P" "$sza" "--pssa --pssa-mode legacy"
    run "t${tau}_s${s}_ipss"   "$P" "$sza" "--pssa"
  done
done
echo "RUNS_DONE"
ls -1 "$D"/case_*.csv | wc -l
