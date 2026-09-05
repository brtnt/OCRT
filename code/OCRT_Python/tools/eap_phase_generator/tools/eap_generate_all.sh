#!/bin/sh
# usage: sh tools/eap_generate_all.sh {representative|all} <outdir>
set -e
MODE="${1:-representative}"; OUT="${2:-inputs/water_iop}"
BIN="${EAP_WRITER:-build/eap_mie_writer}"
mkdir -p "$OUT"
gen() { echo "  $3"; "$BIN" "$1" "$2" "$OUT/$3"; }
# id  Deff  filename                                (대표 조합 = 카탈로그 동결값)
REP="0 6 EAP_00_Diatoms_pennate_D6.mie
1 8 EAP_01_Chlorophytes_D8.mie
2 6 EAP_02_Diatoms_centric_D6.mie
3 6 EAP_03_Cryptophytes_D6.mie
4 6 EAP_04_Cyano_blue_D6.mie
5 6 EAP_05_Cyano_red_D6.mie
6 24 EAP_06_Dinoflagellates_D24.mie
7 6 EAP_07_Eustigmatophytes_D6.mie
8 6 EAP_08_Hapto_Pavlovaceae_D6.mie
9 3 EAP_09_Pelagophytes_D3.mie
10 3 EAP_10_Prasinophytes_D3.mie
11 0.5 EAP_11_Prochlorococcus_D0p5.mie
12 4 EAP_12_Hapto_Prymnesiaceae_D4.mie
13 24 EAP_13_Raphidophytes_D24.mie
14 6 EAP_14_Rhodophytes_D6.mie
15 1.2 EAP_15_Synechococcus_D1p2.mie
16 5 EAP_16_Microcystis_D5.mie"
ALL_EXTRA="0 12 |0 24 |0 48 |1 2 |1 4 |1 6 |2 12 |2 24 |2 48 |3 2 |3 12 |3 24 |3 48 |4 2 |4 12 |4 24 |5 2 |5 12 |5 24 |6 2 |6 6 |6 12 |7 2 |7 12 |7 24 |8 2 |8 12 |8 24 |9 1 |9 2 |9 4 |10 2 |10 4 |10 5 |11 0.4 |11 0.7 |11 0.9 |12 1 |12 2 |12 3 |13 12 |13 48 |13 60 |14 2 |14 12 |14 24 |14 48 |15 0.4 |15 0.8 |15 1.8 |16 3 |16 4 |16 6 "
echo "$REP" | while read id de fn; do gen "$id" "$de" "$fn"; done
if [ "$MODE" = "all" ]; then
  echo "$ALL_EXTRA" | tr '|' '\n' | while read id de; do
    [ -z "$id" ] && continue
    d=$(echo "$de" | tr '.' 'p')
    gen "$id" "$de" "EAP_$(printf '%02d' "$id")_D${d}.mie"
  done
fi
echo "완료: $OUT"
