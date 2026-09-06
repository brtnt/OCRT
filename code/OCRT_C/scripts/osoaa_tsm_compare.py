#!/usr/bin/env python3
"""item#5 TSM(mineral) rrs(0-): FIXEDBULK_SHAREDPHASE_HARNESS §A-§I 정본 규약.
OCRT: fixed-bulk, 무대기(P=0,aod=0,§E), nmw48(§D). OSOAA: HYD.Model 3, thin MOT=0.01(앵커 0.06% 등가).
§C 일관성(블렌드 P11 후방적분=bb/b) 내장. Brown_earth(§H bb/b=0.0099).
인자: <sza> [tsm]"""
import subprocess, os, csv, math, sys, re
import numpy as np
from scipy.interpolate import PchipInterpolator as PCHIP
ROOT="/home/user/ws/SPEEDUP_2026-07-09/osoaa"; SURF=f"{ROOT}/DATABASE/SURF_MATR"
OCRT="/home/user/ws/SPEEDUP_2026-07-09/ocrt_opt/build/v2o_s23"
CONV="/home/user/ws/SPEEDUP_2026-07-09/OCRT_FULLPKG_v1.09_2026-07-08/ocrt/scripts/mie_to_osoaa_extdata.py"
MIE="/home/user/ws/SPEEDUP_2026-07-09/OCRT_FULLPKG_v1.09_2026-07-08/inputs/tsm_ahn/Brown_earth_AHN.mie"
BANDS=[412,443,490,555,660,865]
A_W={412:4.5506e-3,443:7.069e-3,490:1.500e-2,555:5.960e-2,660:4.100e-1,865:4.605}
B_W={412:6.650e-3,443:4.872e-3,490:3.164e-3,555:1.859e-3,660:8.875e-4,865:2.763e-4}
TSM=[0.1,1.0,10.0,50.0]; SZA=[0,40,80]; WORK="/tmp/tsm"; os.makedirs(WORK,exist_ok=True)
DEL=0.039; DD=(1-DEL)/(1+DEL/2)
def p_ray(th): c=np.cos(np.radians(th)); return DD*0.75*(1+c*c)+(1-DD)
def extdata(band):
    out=f"{WORK}/extdata/{band}.extdata"
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out),exist_ok=True)
        subprocess.run(["python3",CONV,MIE,str(band),out],capture_output=True,text=True)
    ext=sca=None; rows=[]
    for ln in open(out):
        s=ln.split()
        if "EXTINCTION_COEF" in ln: ext=float(s[-1])
        elif "SCATTERING_COEF" in ln: sca=float(s[-1])
        elif len(s)>=2:
            try: rows.append((float(s[0]),float(s[1])))
            except ValueError: pass
    a=np.array(rows); return out,ext,sca,a[:,0],a[:,1]
def back_int(th,p):
    t=np.radians(th); s=np.sin(t)
    norm=0.5*np.trapezoid(p*s,t); p=p/norm
    m=th>=90.0
    return 0.5*np.trapezoid(p[m]*s[m],t[m])
def f2(x):
    try: return float(x.replace('D','E'))
    except: return None
def osoaa(band,sza,cs,ext_path,prof):
    res=f"{WORK}/osoaa5/{band}_{sza}_{cs:g}"
    if os.path.exists(f"{res}/Advanced_outputs/Flux.txt"): return res
    os.makedirs(res,exist_ok=True)
    cmd=[f"{ROOT}/exe/OSOAA_MAIN.exe","-OSOAA.ResRoot",res,"-OSOAA.Log","M.Log",
      "-OSOAA.Wa",f"{band/1000:.5f}","-ANG.Thetas",f"{sza}.","-ANG.Rad.NbGauss","48","-ANG.Mie.NbGauss","100",
      "-AP.MOT","0.01","-AP.HR","8.0","-AER.AOTref","0.0",
      "-HYD.Model","3","-HYD.ExtData",ext_path,"-HYD.UserProfile",prof,
      "-PHYTO.Chl","0.0","-SED.Csed","0.0","-YS.Abs440","0.0","-DET.Abs440","0.0",
      "-SEA.Depth","200.0","-SEA.Ind","1.34","-SEA.Wind","3","-SEA.Dir",SURF,
      "-SEA.SurfAlb","0.0","-SEA.BotType","1","-SEA.BotAlb","0.0",
      "-OSOAA.View.Phi","90.","-OSOAA.View.Level","1",
      "-OSOAA.ResFile.vsVZA","RESLUM_vsVZA.txt","-OSOAA.ResFile.Adv.Up","RESLUM_Adv_UP.txt"]
    subprocess.run(cmd,capture_output=True,text=True,timeout=250,env={**os.environ,"OSOAA_NO_DIRECT_GLINT":"1","OSOAA_ROOT":ROOT,"HOME":"/tmp"})
    return res if os.path.exists(f"{res}/Advanced_outputs/Flux.txt") else None
def osoaa_rrs(res):
    V,I=[],[]
    for ln in open(f"{res}/Advanced_outputs/RESLUM_Adv_UP.txt"):
        s=ln.split()
        if len(s)>=7 and s[0]=="27":
            v,ii=f2(s[2]),f2(s[4])
            if v is not None and ii is not None and v>=0: V.append(v); I.append(ii)
    V,I=np.array(V),np.array(I); o=np.argsort(V); V,I=V[o],I[o]
    Lu=float(np.interp(0.0,V,I))
    Ed=None
    for ln in open(f"{res}/Advanced_outputs/Flux.txt"):
        s=ln.split()
        if s and s[0]=="27" and len(s)>4: Ed=f2(s[4]); break
    return Lu/Ed if Ed else float('nan')
def run_ocrt(band,A,B,BB,lut,sza):
    r=subprocess.run([OCRT,"--surface","ocean","--decouple-sunglint","--gas-column-h2o","0","--gas-column-o3","0","--gas-column-no2","0","--gas-column-o2","0","--gas-column-co2","0","--gas-column-ch4","0","--wind-speed","3","--sza",str(sza),"--vza","0",
       "--raa","90","--wavelength",str(band),"--pressure","0","--aod","0",
       "--water-model","iop","--iop-a",f"{A:.6f}","--iop-b",f"{B:.6f}","--iop-bb",f"{BB:.8f}","--iop-phase-lut",lut,
       "--n-mu-water","48"],capture_output=True,text=True,timeout=250,
       env={**os.environ,"OCRT_ADVANCED":"1","OMP_NUM_THREADS":"1"})
    m=re.search(r"rrs0minus=([0-9.eE+-]+)",r.stdout)
    return float(m.group(1)) if m else float('nan')
sza_f=int(sys.argv[1]) if len(sys.argv)>1 else None
tsm_f=float(sys.argv[2]) if len(sys.argv)>2 else None
rows=[]
for cs in TSM:
    if tsm_f is not None and abs(cs-tsm_f)>1e-9: continue
    for band in BANDS:
        ext_path,ext,sca,th,F11=extdata(band)
        astar=ext-sca; bstar=sca
        a_min=cs*astar; b_min=cs*bstar
        A=a_min+A_W[band]; B=b_min+B_W[band]
        pr=p_ray(th)
        blend=(b_min*F11+B_W[band]*pr)/(b_min+B_W[band])
        bbb_blend=back_int(th,blend)
        bbb_min=back_int(th,F11)
        BB_tot=cs*bbb_min*bstar+0.5*B_W[band]   # §C 자기일관: 위상 실측 후방적분(§H 상수 0.0099는 Ahn 문헌값, mie 실측과 파장별 상이 — 기록)
        cons=bbb_blend*B/BB_tot
        lut=f"{WORK}/p11_cs{cs:g}_{band}.csv"
        with open(lut,"w") as f:
            f.write("theta_deg,P11\n")
            for t,p in zip(th,blend): f.write(f"{t:.4f},{p:.6e}\n")
        BB_pass=min(BB_tot,0.49999*B)
        prof=f"{WORK}/prof_cs{cs:g}_{band}.txt"
        with open(prof,"w") as f:
            f.write(f"# mineral Csed={cs} {band}nm\n#\n#\n# depth a b\n#\n0.0 {a_min:.6f} {b_min:.6f}\n200.0 {a_min:.6f} {b_min:.6f}\n")
        for sza in SZA:
            if sza_f is not None and sza!=sza_f: continue
            res=osoaa(band,sza,cs,ext_path,prof)
            if not res: print(f"OSOAA FAIL {band} {sza} {cs}"); continue
            oR=osoaa_rrs(res); cR=run_ocrt(band,A,B,BB_pass,lut,sza)
            rows.append(dict(tsm=cs,band=band,sza=sza,OCRT_rrs0minus=cR,OSOAA_rrs0minus=oR,consistency=round(cons,4)))
import statistics
if rows:
    new=not os.path.exists("cmp_item5_tsm.csv")
    with open("cmp_item5_tsm.csv","a",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        if new: w.writeheader()
        w.writerows(rows)
    dR=[100*(r['OCRT_rrs0minus']/r['OSOAA_rrs0minus']-1) for r in rows if r['OSOAA_rrs0minus']==r['OSOAA_rrs0minus']]
    badc=[r for r in rows if abs(r['consistency']-1)>0.02]
    print(f"sza={sza_f}: rrs MAPE={statistics.mean(abs(x) for x in dR):.2f}% bias={statistics.mean(dR):+.2f}% | §C 위반 {len(badc)}건")
    print(f"+{len(rows)} rows -> cmp_item5_tsm.csv")
