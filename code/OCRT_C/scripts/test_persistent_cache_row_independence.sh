#!/usr/bin/env bash
set -euo pipefail
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BIN=${OCRT_BIN:-"$ROOT/build/ocrt"}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/cache" "$TMP/prebuild" "$TMP/serial" "$TMP/parallel" "$TMP/process"
cd "$ROOT"
COMMON=(--surface ocean --wind-speed 3 --sza 40 --pressure 1013.25 --aod 0 --n-water 1.34
 --water-model ocrt --ocrt-chl 0 --ocrt-tsm 0.5 --ocrt-adom440 0
 --ocrt-adom-slope 0.014 --ocrt-tsm-species red_clay
 --n-mu-water 24 --water-m-max 8 --lut-vza-step 10 --lut-vza-max 60 --lut-raa-step 30
 --gas-column-h2o 0 --gas-column-o3 0 --gas-column-no2 0
 --gas-column-o2 0 --gas-column-co2 0 --gas-column-ch4 0 --decouple-sunglint)
BASEENV=(OCRT_ADVANCED=1 OCRT_WATER_PARTICLE_KERNEL=direct OCRT_WATER_VALUE_NPHI=128)
for wl in 443 555 660; do
  env OMP_NUM_THREADS=1 "${BASEENV[@]}" OCRT_SURFACE_PERSIST_CACHE_MODE=build \
    OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/cache" "$BIN" "${COMMON[@]}" \
    --wavelength "$wl" --output-full-grid "$TMP/prebuild/$wl.csv" >/dev/null 2>"$TMP/prebuild/$wl.err"
done
cat > "$TMP/serial.csv" <<CSV
out,wavelength
$TMP/serial/row1_443.csv,443
$TMP/serial/row2_555.csv,555
$TMP/serial/row3_660.csv,660
$TMP/serial/row4_443.csv,443
CSV
cat > "$TMP/parallel.csv" <<CSV
out,wavelength
$TMP/parallel/row1_443.csv,443
$TMP/parallel/row2_555.csv,555
$TMP/parallel/row3_660.csv,660
$TMP/parallel/row4_443.csv,443
CSV
env OMP_NUM_THREADS=1 "${BASEENV[@]}" OCRT_BATCH_DIAG=1 \
  OCRT_SURFACE_PERSIST_CACHE_MODE=required OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/cache" \
  "$BIN" "${COMMON[@]}" --wavelength 443 --batch-full-grid "$TMP/serial.csv" \
  >/dev/null 2>"$TMP/serial.err"
env OMP_NUM_THREADS=4 "${BASEENV[@]}" OCRT_BATCH_DIAG=1 \
  OCRT_SURFACE_PERSIST_CACHE_MODE=required OCRT_SURFACE_PERSIST_CACHE_DIR="$TMP/cache" \
  "$BIN" "${COMMON[@]}" --wavelength 443 --batch-full-grid "$TMP/parallel.csv" \
  >/dev/null 2>"$TMP/parallel.err"
for spec in '1 443' '2 555' '3 660' '4 443'; do
  set -- $spec
  cmp "$TMP/serial/row${1}_${2}.csv" "$TMP/parallel/row${1}_${2}.csv"
done
grep -F '# batch-full-grid done: 4 rows, 0 failed' "$TMP/parallel.err" >/dev/null
! grep -F '[surface-pcache] miss:' "$TMP/parallel.err" >/dev/null
! grep -R '#pragma omp critical(p2b_rowcache)' "$ROOT/src" >/dev/null
python - "$TMP/parallel.err" <<'PYDIAG'
import re, sys
rows=[]
for line in open(sys.argv[1], encoding='utf-8', errors='replace'):
    m=re.search(r'row\s+\d+: tid (\d+)\s+([0-9.]+)\.\.([0-9.]+)', line)
    if m: rows.append((int(m.group(1)),float(m.group(2)),float(m.group(3))))
if len(rows)!=4: raise SystemExit(f'missing BATCH_DIAG timeline: {rows}')
if len({r[0] for r in rows}) < 2: raise SystemExit(f'only one worker used: {rows}')
if max(r[1] for r in rows) >= min(r[2] for r in rows):
    raise SystemExit(f'rows did not overlap: {rows}')
print('OpenMP rows overlap without p2b cache lock')
PYDIAG

python - "$ROOT" "$BIN" "$TMP" <<'PY'
from pathlib import Path
import os, subprocess, sys, time
root=Path(sys.argv[1]); binary=str(Path(sys.argv[2]).resolve()); tmp=Path(sys.argv[3])
common=['--surface','ocean','--wind-speed','3','--sza','40','--pressure','1013.25','--aod','0','--n-water','1.34',
 '--water-model','ocrt','--ocrt-chl','0','--ocrt-tsm','0.5','--ocrt-adom440','0','--ocrt-adom-slope','0.014','--ocrt-tsm-species','red_clay',
 '--n-mu-water','24','--water-m-max','8','--lut-vza-step','10','--lut-vza-max','60','--lut-raa-step','30',
 '--gas-column-h2o','0','--gas-column-o3','0','--gas-column-no2','0','--gas-column-o2','0','--gas-column-co2','0','--gas-column-ch4','0','--decouple-sunglint']
env=os.environ.copy(); env.update({'OMP_NUM_THREADS':'1','OCRT_ADVANCED':'1','OCRT_WATER_PARTICLE_KERNEL':'direct','OCRT_WATER_VALUE_NPHI':'128','OCRT_SURFACE_PERSIST_CACHE_MODE':'required','OCRT_SURFACE_PERSIST_CACHE_DIR':str(tmp/'cache')})
rows=[(1,443),(2,555),(3,660),(4,443)]; ps=[]; t0=time.monotonic()
for i,wl in rows:
    out=tmp/'process'/f'row{i}_{wl}.csv'
    f1=open(tmp/f'p{i}.out','wb'); f2=open(tmp/f'p{i}.err','wb')
    p=subprocess.Popen([binary,*common,'--wavelength',str(wl),'--output-full-grid',str(out)],cwd=root,env=env,stdout=f1,stderr=f2)
    ps.append((i,wl,p,time.monotonic()-t0,f1,f2))
ends=[]
for i,wl,p,start,f1,f2 in ps:
    rc=p.wait(); ends.append(time.monotonic()-t0); f1.close(); f2.close()
    if rc: raise SystemExit(f'process row {i} rc={rc}')
    if (tmp/'serial'/f'row{i}_{wl}.csv').read_bytes() != (tmp/'process'/f'row{i}_{wl}.csv').read_bytes():
        raise SystemExit(f'process row {i} mismatch')
if max(x[3] for x in ps) >= min(ends): raise SystemExit('process launch overlap failed')
print('independent processes overlap and are byte-identical')
PY

echo 'PASS: immutable persistent cache has no row build dependency; serial, 4-thread and 4-process outputs are byte-identical'
