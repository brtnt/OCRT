#!/usr/bin/env python3
"""item#4 aDOM(CDOM/yellow substance): OCRT <-> OSOAA, TOA IQU + Rrs(0-).
- CDOM = 흡수만(a_CDOM(lambda)=a440*exp(-S*(lambda-440)), b=bb=0). 산란/위상 없음 -> in-water 흡수 응답만 격리.
- parity: OCRT --water-model ocrt --ocrt-adom440/--ocrt-adom-slope == OSOAA -YS.Abs440/-YS.Swa (둘 다 동일 지수식; slope 0.014 명시).
- DET=SED=PHYTO.Chl=0. 물 IOP는 #2와 동일(OCRT 기본 = OSOAA 일치).
- 규약(#2와 동일): REFL/glint-free(raa90)/U부호(OCRT_U=-OSOAA_U)/PCHIP/nadir Q는 |Q|/Rrs는 0-(level27 water-side).
인자: python3 osoaa_adom_compare.py <sza> <adom440>  (둘 다 옵션; 없으면 전체)
"""
import subprocess, os, csv, math, sys
import numpy as np
from scipy.interpolate import PchipInterpolator as PCHIP
ROOT="/home/user/osoaa_build"; EXE=f"{ROOT}/exe/OSOAA_MAIN.exe"; SURF=f"{ROOT}/DATABASE/SURF_MATR"
OCRT="./build/v2_solver_vk"
BANDS=[412,443,490,555,660,865]; WUM={b:f"{b/1000:.5f}" for b in BANDS}
TAU={412:0.318540221,443:0.236054530,490:0.155974381,555:0.093751620,660:0.046362496,865:0.015540855}
SZA=[0,40,80]; VZA_AIR=[0,30,60]; NW=1.34
ADOM=[0.01,0.1,1.0]; SLOPE=0.014
def water_vza(air): return math.degrees(math.asin(min(1.0,math.sin(math.radians(air))/NW)))
os.makedirs("/tmp/adom",exist_ok=True)

def osoaa(band,sza,adom):
    tag=f"{band}_{sza}_{adom:g}"
    res=f"/tmp/adom/{tag}"
    if os.path.exists(f"{res}/Advanced_outputs/Flux.txt"): return res
    os.makedirs(res,exist_ok=True)
    cmd=[EXE,"-OSOAA.ResRoot",res,"-OSOAA.Log","M.Log","-OSOAA.Wa",WUM[band],
         "-ANG.Thetas",f"{sza}.","-ANG.Rad.NbGauss","48","-ANG.Mie.NbGauss","100",
         "-AP.MOT",f"{TAU[band]}","-AP.HR","8.0","-AER.AOTref","0.0",
         "-PHYTO.Chl","0.0","-SED.Csed","0.0",
         "-YS.Abs440",f"{adom}","-YS.Swa",f"{SLOPE}","-DET.Abs440","0.0",
         "-SEA.Depth","200.0","-SEA.Ind","1.34","-SEA.Wind","3","-SEA.Dir",SURF,
         "-SEA.SurfAlb","0.0","-SEA.BotType","1","-SEA.BotAlb","0.0",
         "-OSOAA.View.Phi","90.","-OSOAA.View.Level","1",
         "-OSOAA.ResFile.vsVZA","RESLUM_vsVZA.txt","-OSOAA.ResFile.Adv.Up","RESLUM_Adv_UP.txt"]
    subprocess.run(cmd,capture_output=True,text=True,timeout=250,env={**os.environ,"OSOAA_NO_DIRECT_GLINT":"1"})
    return res if os.path.exists(f"{res}/Advanced_outputs/Flux.txt") else None

def f2(x):
    try: return float(x.replace('D','E'))
    except: return None

def parse_osoaa(res,sza):
    mu=math.cos(math.radians(sza))
    vz,refl=[],[]
    for ln in open(f"{res}/Standard_outputs/RESLUM_vsVZA.txt"):
        s=ln.split()
        if len(s)==7:
            v=f2(s[0]); R=f2(s[3])
            if v is not None and R is not None and v>=0: vz.append(v); refl.append(R)
    vz,refl=np.array(vz),np.array(refl); o=np.argsort(vz); vz,refl=vz[o],refl[o]
    vz,iu=np.unique(vz,return_index=True); refl=refl[iu]; fREFL=PCHIP(vz,refl)
    def adv(level):
        V,I,Q,U=[],[],[],[]
        for ln in open(f"{res}/Advanced_outputs/RESLUM_Adv_UP.txt"):
            s=ln.split()
            if len(s)>=7 and s[0]==str(level):
                v=f2(s[2]); ii=f2(s[4]); q=f2(s[5]); u=f2(s[6])
                if None in (v,ii,q,u): continue
                if v>=0: V.append(v); I.append(ii); Q.append(q); U.append(u)
        V=np.array(V); o=np.argsort(V); V,I,Q,U=V[o],np.array(I)[o],np.array(Q)[o],np.array(U)[o]
        V,iu=np.unique(V,return_index=True); I,Q,U=I[iu],Q[iu],U[iu]
        return V,PCHIP(V,I),PCHIP(V,Q),PCHIP(V,U)
    V0,fI0,fQ0,fU0=adv(0); V27,fI27,fQ27,fU27=adv(27)
    ed={}
    for ln in open(f"{res}/Advanced_outputs/Flux.txt"):
        s=ln.split()
        if s and s[0].isdigit():
            e=f2(s[4]) if len(s)>4 else None
            if e is not None: ed[int(s[0])]=e
    return dict(fREFL=fREFL, fQ0=fQ0,fU0=fU0,mu=mu, fI27=fI27,fQ27=fQ27,fU27=fU27,
               Ed0=ed.get(0,np.nan), Ed27=ed.get(27,np.nan))

def run_ocrt(band,sza,vza,adom):
    out=subprocess.run([OCRT,"--surface","ocean","--decouple-sunglint","--gas-column-h2o","0","--gas-column-o3","0","--gas-column-no2","0","--gas-column-o2","0","--gas-column-co2","0","--gas-column-ch4","0","--water-model","ocrt","--ocrt-chl","0","--ocrt-tsm","0",
        "--ocrt-adom440",str(adom),"--ocrt-adom-slope",str(SLOPE),"--wind-speed","3",
        "--sza",str(sza),"--vza",str(vza),"--raa","90","--wavelength",str(band),
        "--pressure","1013.25","--n-mu-water","8"],
        capture_output=True,text=True,timeout=120).stdout.split('\n')[-2]
    toks=out.split()
    rI,rQ,rU=float(toks[0]),float(toks[1]),float(toks[2])
    d={}
    for t in toks:
        if '=' in t:
            k,v=t.split('=',1)
            try: d[k]=float(v)
            except: pass
    return rI,rQ,rU,d.get('rrs0minus',np.nan),d.get('a_cdom_used',np.nan),d.get('a_total',np.nan)

sza_f=int(sys.argv[1]) if len(sys.argv)>1 else None
adom_f=float(sys.argv[2]) if len(sys.argv)>2 else None
rows=[]
for adom in ADOM:
    if adom_f is not None and abs(adom-adom_f)>1e-9: continue
    for band in BANDS:
        for sza in SZA:
            if sza_f is not None and sza!=sza_f: continue
            res=osoaa(band,sza,adom)
            if not res: sys.stderr.write(f"OSOAA FAIL {band} {sza} {adom}\n"); continue
            P=parse_osoaa(res,sza); mu=P['mu']
            for va in VZA_AIR:
                cI,cQ,cU,crrs,acd,atot=run_ocrt(band,sza,va,adom)
                oI=float(P['fREFL'](va)); oQ=float(P['fQ0'](va))/mu; oU=float(P['fU0'](va))/mu
                vw=water_vza(va); oRrs=float(P['fI27'](vw))/P['Ed27']
                rows.append((adom,band,sza,va,vw,cI,oI,cQ,oQ,cU,oU,crrs,oRrs,acd,atot))
out="cmp_item4_adom.csv"; new=not os.path.exists(out)
with open(out,"a",newline="") as f:
    w=csv.writer(f)
    if new: w.writerow(["adom440","band","sza","vza_air","vza_water","OCRT_TOA_I","OSOAA_TOA_I","OCRT_TOA_Q","OSOAA_TOA_Q","OCRT_TOA_U","OSOAA_TOA_U","OCRT_rrs0minus","OSOAA_rrs0minus","OCRT_a_cdom","OCRT_a_total"])
    w.writerows(rows)
import statistics
dI=[100*(r[5]/r[6]-1) for r in rows if r[6]]
dR=[100*(r[11]/r[12]-1) for r in rows if r[12] and r[3]==0]
print(f"sza={sza_f} adom={adom_f}: TOA_I MAPE={statistics.mean(abs(x) for x in dI):.2f}%" + (f" | Rrs(0-)@nadir bias={statistics.mean(dR):+.2f}%" if dR else ""))
print(f"+{len(rows)} rows -> {out}")
