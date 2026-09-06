#!/usr/bin/env python3
"""item#1 aerosol IQU: OSOAA(black ocean TOA) ↔ OCRT. raa=90. 절차문서 §2 준수.
- OSOAA: glint-free, REFL, matched(ExtData ssa/phase, AOTref=밴드AOD, AP.HA 2.0, Tronca1). 6밴드×3aot×3sza.
- 추출: ρ_I=REFL, ρ_Q=Q/μsun, ρ_U=U/μsun (Adv_UP level0). PCHIP→vza{0,30,60}. U: OCRT_U=−OSOAA_U.
"""
import subprocess, os, csv, math, sys
import numpy as np
from scipy.interpolate import PchipInterpolator as PCHIP
ROOT="/home/user/osoaa_build"; EXE=f"{ROOT}/exe/OSOAA_MAIN.exe"; SURF=f"{ROOT}/DATABASE/SURF_MATR"
BANDS=[412,443,490,555,660,865]; WUM={b:f"{b/1000:.5f}" for b in BANDS}
RATIO={412:1.1309,443:1.1152,490:1.0937,555:1.0691,660:1.0409,865:1.0000}
TAU={412:0.318540221,443:0.236054530,490:0.155974381,555:0.093751620,660:0.046362496,865:0.015540855}
AOT865=[0.05,0.3,1.0]; SZA=[0,40,80]; VZA_T=[0,30,60]; PHI=90
os.makedirs("/tmp/orayaer",exist_ok=True); os.makedirs("/tmp/mie_dir",exist_ok=True)

def run(band,sza,aot865):
    aot=aot865*RATIO[band]; res=f"/tmp/orayaer/{band}_{sza}_{aot865}"
    if os.path.exists(f"{res}/Standard_outputs/RESLUM_vsVZA.txt"): return res  # resume
    os.makedirs(res,exist_ok=True)
    cmd=[EXE,"-OSOAA.ResRoot",res,"-OSOAA.Log","M.Log","-OSOAA.Wa",WUM[band],
         "-ANG.Thetas",f"{sza}.","-ANG.Rad.NbGauss","48","-ANG.Mie.NbGauss","100",
         "-AP.MOT",f"{TAU[band]}","-AP.HR","8.0","-AP.HA","2.0","-SOS.IGmax","100",
         "-AER.Model","4","-AER.ExtData",f"/tmp/M80C_{band}.extdata","-AER.DirMie","/tmp/mie_dir",
         "-AER.Tronca","1","-AER.Waref",WUM[band],"-AER.AOTref",f"{aot:.5f}",
         "-PHYTO.Chl","0.0","-SED.Csed","0.0","-YS.Abs440","1000.0","-DET.Abs440","0.0",
         "-SEA.Depth","1000.0","-SEA.Ind","1.34","-SEA.Wind","3","-SEA.Dir",SURF,
         "-SEA.SurfAlb","0.0","-SEA.BotType","1","-SEA.BotAlb","0.0",
         "-OSOAA.View.Phi",f"{PHI}.","-OSOAA.View.Level","1",
         "-OSOAA.ResFile.vsVZA","RESLUM_vsVZA.txt","-OSOAA.ResFile.Adv.Up","RESLUM_Adv_UP.txt"]
    subprocess.run(cmd,capture_output=True,text=True,timeout=250,env={**os.environ,"OSOAA_NO_DIRECT_GLINT":"1"})
    return res if os.path.exists(f"{res}/Standard_outputs/RESLUM_vsVZA.txt") else None

def parse(res,sza):
    mu=math.cos(math.radians(sza)); vz,refl=[],[]
    for ln in open(f"{res}/Standard_outputs/RESLUM_vsVZA.txt"):
        s=ln.split()
        if len(s)==7:
            try: v=float(s[0]); R=float(s[3].replace('D','E'))
            except: continue
            if v>=0: vz.append(v); refl.append(R)
    vq,Q,U=[],[],[]
    for ln in open(f"{res}/Advanced_outputs/RESLUM_Adv_UP.txt"):
        s=ln.split()
        if len(s)>=7 and s[0]=='0':
            try:
                v=float(s[2])
                if v>=0: vq.append(v); Q.append(float(s[5].replace('D','E'))/mu); U.append(float(s[6].replace('D','E'))/mu)
            except: continue
    vz,refl=np.array(vz),np.array(refl); o=np.argsort(vz); vz,refl=vz[o],refl[o]
    vz,iu=np.unique(vz,return_index=True); refl=refl[iu]
    vq=np.array(vq); o2=np.argsort(vq); vq,Q,U=vq[o2],np.array(Q)[o2],np.array(U)[o2]
    vq,iu2=np.unique(vq,return_index=True); Q,U=Q[iu2],U[iu2]
    fI,fQ,fU=PCHIP(vz,refl),PCHIP(vq,Q),PCHIP(vq,U)
    return {vt:(float(fI(vt)),float(fQ(vt)),float(fU(vt))) for vt in VZA_T}

# OCRT aerosol grid (전 raa, 축퇴 fallback)
ocrt={}
for r in csv.DictReader(open("OCRT_item1_grid.csv")):
    if r['atm']=='ray_aer':
        ocrt[(int(r['band']),float(r['aot865']),int(r['sza']),int(r['vza']),int(r['raa']))]=(float(r['rho_I']),float(r['rho_Q']),float(r['rho_U']))
def olk(band,a,sza,vza):
    for raa in (90,0,45,135,180):
        if (band,a,sza,vza,raa) in ocrt: return ocrt[(band,a,sza,vza,raa)]
    return (None,None,None)

aot_filter=float(sys.argv[1]) if len(sys.argv)>1 else None  # chunk by aot865
rows=[]
for band in BANDS:
    for sza in SZA:
        for a in AOT865:
            if aot_filter is not None and a!=aot_filter: continue
            res=run(band,sza,a)
            if not res: sys.stderr.write(f"FAIL {band} {sza} {a}\n"); continue
            ext=parse(res,sza)
            for vt in VZA_T:
                oi,oq,ou=ext[vt]; ci,cq,cu=olk(band,a,sza,vt)
                if ci is None: continue
                rows.append((band,a,sza,vt,ci,oi,cq,oq,cu,ou))
# append 저장
out="cmp_item1_rayaer.csv"; new=not os.path.exists(out)
with open(out,"a",newline="") as f:
    w=csv.writer(f)
    if new: w.writerow(["band","aot865","sza","vza","OCRT_I","OSOAA_I","OCRT_Q","OSOAA_Q","OCRT_U","OSOAA_U"])
    w.writerows(rows)
import statistics
dI=[100*(r[4]/r[5]-1) for r in rows if r[5]]
if dI: print(f"chunk aot={aot_filter}: N={len(dI)} ρ_I MAPE={statistics.mean(abs(x) for x in dI):.2f}% bias={statistics.mean(dI):+.2f}% max={max(abs(x) for x in dI):.1f}%")
print(f"+{len(rows)} rows -> {out}")
