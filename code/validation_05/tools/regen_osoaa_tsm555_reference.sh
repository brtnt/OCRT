#!/bin/bash
set -e
M=~/ocrt/runtime/MIGRATION_PKG_2026-08-19
O=$M/02_OSOAA/OSOAA_OCRT_REFERENCE
CONV=$M/01_OCRT_C/scripts/mie_to_osoaa_extdata.py
MIE=$1; RES=$2
mkdir -p $RES
python3 $CONV "$MIE" 555 $RES/extdata_555.txt
# mineral-only profile: a_min=0.5*astar(555), b_min=0.5*bstar(555)
python3 - "$RES" <<'PY'
import sys,os,numpy as np
C=os.path.expanduser('~/ocrt/runtime/MIGRATION_PKG_2026-08-19/01_OCRT_C')
def at555(p):
    d=[l.split() for l in open(p) if l and not l.lstrip().startswith('#')]
    d=[(float(r[0]),float(r[1])) for r in d if len(r)>=2]
    import numpy as np; a=np.array(d); i=np.argmin(abs(a[:,0]-555)); return a[i,1]
astar=at555(f'{C}/inputs/tsm_ahn/astarmin_redclay.txt'); bstar=at555(f'{C}/inputs/tsm_ahn/bstarmin_redclay.txt')
am=0.5*astar; bm=0.5*bstar
open(sys.argv[1]+'/prof.txt','w').write(f"# red_clay Csed=0.5 555nm mineral-only\n#\n#\n# depth a b\n#\n0.0 {am:.8f} {bm:.8f}\n36.343 {am:.8f} {bm:.8f}\n")
print(f"prof: a_min={am:.8f} b_min={bm:.8f}")
PY
cd $RES
env OSOAA_ROOT=$O HOME=/tmp OSOAA_NO_DIRECT_GLINT=1 OCRT_WL_FILE=$RES/wl.txt \
 $O/exe/OSOAA_MAIN.exe -OSOAA.ResRoot $RES -OSOAA.Log M.Log \
 -OSOAA.Wa 0.55500 -ANG.Thetas 40. -ANG.Rad.NbGauss 48 -ANG.Mie.NbGauss 100 \
 -AP.MOT 0.0935485 -AP.HR 8.0 -AER.AOTref 0.0 \
 -HYD.Model 3 -HYD.ExtData $RES/extdata_555.txt -HYD.UserProfile $RES/prof.txt \
 -PHYTO.Chl 0.0 -SED.Csed 0.0 -YS.Abs440 0.0 -DET.Abs440 0.0 \
 -SEA.Depth 36.343 -SEA.Ind 1.34 -SEA.Wind 3 -SEA.Dir $O/DATABASE/SURF_MATR \
 -SEA.SurfAlb 0.0 -SEA.BotType 1 -SEA.BotAlb 0.0 \
 -OSOAA.View.Phi 90. -OSOAA.View.Level 1 \
 -OSOAA.ResFile.vsVZA RESLUM_vsVZA.txt -OSOAA.ResFile.Adv.Up RESLUM_Adv_UP.txt > run.log 2>&1
echo "RUN_DONE $RES"
