#!/bin/bash
# 사용자 지정 비교 구성 (2026-09-05, v1.11 지상앵커 확정판)
#   고정: VZA 55, RAA 90, AOD(555) 0.1, 에어로졸 M80C, 표면 black_fresnel_ocean,
#         풍속 5 m/s, US62, 기체흡수 ON, 층수 400
#   스윕: SZA 0..85 (5도)
#   모드: pp(--pssa 없음) / ipss(--pssa)
set -u
V=/root/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C
D=/root/ocrt/ipss_bench/sza_sweep
mkdir -p "$D"; cd "$V"
run () {   # $1=wl $2=sza $3=tag $4=modeargs
  local out="$D/wl$1_sza$2_$3.txt"
  [ -s "$out" ] && return 0
  env OMP_NUM_THREADS=1 OCRT_ADVANCED=1 ./build/ocrt \
      --wavelength "$1" --sza "$2" --vza 55 --raa 90 \
      --surface black_fresnel_ocean --wind-speed 5 \
      --mie inputs/M80C.mie --aod-555 0.1 \
      --n-layers 400 $4 2>/dev/null | tail -1 > "$out"
}
N=0
for wl in 412 555 865; do
  for sza in 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85; do
    run "$wl" "$sza" pp   ""        &
    run "$wl" "$sza" ipss "--pssa"  &
    N=$((N+2)); if [ $((N % 12)) -eq 0 ]; then wait; fi
  done
done
wait
echo "RUNS_DONE files=$(ls -1 $D/*.txt | wc -l)"
