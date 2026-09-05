#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$ROOT"
CC_BIN=${CC:-gcc}
OUT=build/test_surface_aw_interface_conventions
mkdir -p build
"$CC_BIN" -std=c11 -O2 -Isrc \
  tests/test_surface_aw_interface_conventions.c \
  src/shared/surface.c \
  src/shared/surface_multibounce.c \
  -o "$OUT" -lm
"$OUT"

# FIX1/FIX2 source contract.
grep -Fq 'const double cphi_r = -cphi;' src/shared/surface.c
grep -Fq 'const double sphi_r = -sphi;' src/shared/surface.c
grep -Fq 'const double mu_i_s = -mu_i;' src/shared/surface.c
grep -Fq 'const double mu_o_s = -mu_o;' src/shared/surface.c
# FIX3 source contract: incoming-U column only, m>0.
grep -Fq 'T_pair[0*3 + 2] = -T_pair[0*3 + 2];' src/shared/surface.c
grep -Fq 'T_pair[1*3 + 2] = -T_pair[1*3 + 2];' src/shared/surface.c
# FIX4 source contract: diffuse-top molecular incoming-U column must close
# against rt_sos_operator_apply_vector under the same sine-Fourier U field.
grep -Fq 'const double ruI = ray_on ? -gamma2' src/rt_water_rt.c
grep -Fq 'const double ruQ = ray_on ? -alpha2' src/rt_water_rt.c
grep -Fq 'const double ruU = ray_on ? -alpha2' src/rt_water_rt.c

echo 'PASS: stage-2 interface FIX1+FIX2+FIX3 + diffuse-top U-column source audit'
