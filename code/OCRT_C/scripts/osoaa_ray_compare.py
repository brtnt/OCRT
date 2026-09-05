#!/usr/bin/env python3
"""item#1 Rayleigh: OSOAA(black Fresnel ocean TOA) ↔ OCRT 대조. raa=90(invariant) 먼저.
- OSOAA: View.Level 1(TOA), MOT=OCRT bodhaine tau_R(matched), YS.Abs440 1000(black water), coxmunk wind3.
- 추출: ρ_I=REFL(vsVZA), ρ_Q=Q/μsun, ρ_U=U/μsun (Adv_UP level0). PCHIP→vza{0,30,60}.
- TOA air-side라 OCRT vza 직접 비교(Snell 불필요).
"""
import subprocess, os, csv, math, sys
import numpy as np
from scipy.interpolate import PchipInterpolator as PCHIP

ROOT="/home/claude/osoaa_build"; EXE=f"{ROOT}/exe/OSOAA_MAIN.exe"; SURF=f"{ROOT}/DATABASE/SURF_MATR"
TAU={412:0.318540221,443:0.2360545301,490:0.1559743814,555:0.09375162019,660:0.04636249595,865:0.01554085494}
WUM={b:f"{b/1000:.5f}" for b in TAU}
BANDS=[412,443,490,555,660,865]; SZA=[0,40,80]; VZA_T=[0,30,60]; RAA=90; PHI=90  # raa90↔Phi90 invariant

def run(band,sza):
    res=f"/tmp/oray/{band}_{sza}"; os.makedirs(res,exist_ok=True)
    cmd=[EXE,"-OSOAA.ResRoot",res,"-OSOAA.Log","Main.Log","-OSOAA.Wa",WUM[band],
         "-ANG.Thetas",f"{sza}.","-ANG.Rad.NbGauss","48","-ANG.Mie.NbGauss","100",
         "-AP.MOT",f"{TAU[band]}","-AP.HR","8.0","-AER.AOTref","0.0",
         "-PHYTO.Chl","0.0","-SED.Csed","0.0","-YS.Abs440","1000.0","-DET.Abs440","0.0",
         "-SEA.Depth","1000.0","-SEA.Ind","1.34","-SEA.Wind","3","-SEA.Dir",SURF,
         "-SEA.SurfAlb","0.0","-SEA.BotType","1","-SEA.BotAlb","0.0",
         "-OSOAA.View.Phi",f"{PHI}.","-OSOAA.View.Level","1",
         "-OSOAA.ResFile.vsVZA","RESLUM_vsVZA.txt","-OSOAA.ResFile.Adv.Up","RESLUM_Adv_UP.txt"]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=250,
                     env={**os.environ,"OSOAA_NO_DIRECT_GLINT":"1"})  # OCRT는 glint-free → OSOAA도 억제
    return res if os.path.exists(f"{res}/Standard_outputs/RESLUM_vsVZA.txt") else None

def parse(res,sza):
    mu=math.cos(math.radians(sza))
    # vsVZA: VZA SCA I REFL POL LPOL REFL_POL  (VZA>0 = azimuth 90 = raa90)
    vz,refl=[],[]
    for ln in open(f"{res}/Standard_outputs/RESLUM_vsVZA.txt"):
        s=ln.split()
        if len(s)==7:
            try: v=float(s[0]); R=float(s[3].replace('D','E'))
            except: continue
            if v>=0: vz.append(v); refl.append(R)
    # Adv_UP level0: LEVEL Z VZA SCA I Q U ...  (Q,U → /mu)
    vq,Q,U=[],[],[]
    for ln in open(f"{res}/Advanced_outputs/RESLUM_Adv_UP.txt"):
        s=ln.split()
        if len(s)>=7 and s[0]=='0':
            try:
                v=float(s[2])
                if v>=0:
                    vq.append(v); Q.append(float(s[5].replace('D','E'))/mu); U.append(float(s[6].replace('D','E'))/mu)
            except: continue
    vz,refl=np.array(vz),np.array(refl); o=np.argsort(vz); vz,refl=vz[o],refl[o]
    vz,iu=np.unique(vz,return_index=True); refl=refl[iu]
    vq=np.array(vq); o2=np.argsort(vq); vq,Q,U=vq[o2],np.array(Q)[o2],np.array(U)[o2]
    vq,iu2=np.unique(vq,return_index=True); Q,U=Q[iu2],U[iu2]
    fI=PCHIP(vz,refl); fQ=PCHIP(vq,Q); fU=PCHIP(vq,U)
    return {vt:(float(fI(vt)),float(fQ(vt)),float(fU(vt))) for vt in VZA_T}

# OCRT Rayleigh grid 로드 (전 raa; 축퇴 fallback용)
ocrt_all={}
for r in csv.DictReader(open("OCRT_item1_grid.csv")):
    if r['atm']=='rayleigh':
        ocrt_all[(int(r['band']),int(r['sza']),int(r['vza']),int(r['raa']))]=(float(r['rho_I']),float(r['rho_Q']),float(r['rho_U']))
def ocrt_lookup(band,sza,vza):
    # raa=90 우선; sza=0 또는 vza=0(축퇴)이면 저장된 대표(raa=0) 사용
    for raa in (90,0,45,135,180):
        if (band,sza,vza,raa) in ocrt_all: return ocrt_all[(band,sza,vza,raa)]
    return (None,None,None)

rows=[]; print(f"{'band':>5}{'sza':>4}{'vza':>4} | {'OCRT_I':>10}{'OSOAA_I':>10}{'dI%':>7} | {'OCRT_Q':>10}{'OSOAA_Q':>10}")
for band in BANDS:
    for sza in SZA:
        res=run(band,sza)
        if not res: print(f"FAIL {band} {sza}"); continue
        ext=parse(res,sza)
        for vt in VZA_T:
            oi,oq,ou=ext[vt]; ci,cq,cu=ocrt_lookup(band,sza,vt)
            if ci is None: continue
            dI=100*(ci/oi-1) if oi else float('nan')
            rows.append((band,sza,vt,ci,oi,dI,cq,oq,cu,ou))
            if sza==40 or (vt==0):  # 대표만 출력
                print(f"{band:>5}{sza:>4}{vt:>4} | {ci:>10.4e}{oi:>10.4e}{dI:>+7.1f} | {cq:>10.3e}{oq:>10.3e}")
# MAPE
import statistics
dIs=[r[5] for r in rows if not math.isnan(r[5])]
print(f"\n=== Rayleigh raa=90 ρ_I: N={len(dIs)}  MAPE={statistics.mean(abs(x) for x in dIs):.2f}%  bias={statistics.mean(dIs):+.2f}%  max|d|={max(abs(x) for x in dIs):.1f}% ===")
# 저장
with open("cmp_item1_rayleigh.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["band","sza","vza","OCRT_I","OSOAA_I","dI%","OCRT_Q","OSOAA_Q","OCRT_U","OSOAA_U"]); w.writerows(rows)
print("saved cmp_item1_rayleigh.csv")
