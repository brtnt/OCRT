#!/usr/bin/env python3
"""item#3 chl(phyto) rrs(0-): #5 기계 + blend_chl.mie(value LUT) 재료.
NOTE_item3_chl_setup 규약. OSOAA HYD.Model 3(외부 phase 주입, PM 경로 회피 = value LUT).
DT-HIOM table 커널이 double-cap 우회하는지 검증. 인자: <sza> [chl]"""
import subprocess, os, csv, math, sys, re
import numpy as np
ROOT="/home/claude/ws/SPEEDUP_2026-07-09/osoaa"; SURF=f"{ROOT}/DATABASE/SURF_MATR"
OCRT="/home/claude/ws/SPEEDUP_2026-07-09/ocrt_opt/build/v2o_hiom"
CONV="/home/claude/ws/SPEEDUP_2026-07-09/OCRT_FULLPKG_v1.09_2026-07-08/ocrt/scripts/mie_to_osoaa_extdata.py"
MIEDIR="/home/claude/ws/SPEEDUP_2026-07-09/OCRT_FULLPKG_v1.09_2026-07-08/aux/phase_mie"
BANDS=[412,443,490,555,660,865]
A_W={412:4.5506e-3,443:7.069e-3,490:1.500e-2,555:5.960e-2,660:4.100e-1,865:4.605}
B_W={412:6.650e-3,443:4.872e-3,490:3.164e-3,555:1.859e-3,660:8.875e-4,865:2.763e-4}
CHL=[0.03,0.3,3.0,30.0]; SZA=[0,40,80]; WORK="/tmp/chl"; os.makedirs(WORK,exist_ok=True)
DEL=0.039; DD=(1-DEL)/(1+DEL/2)
BANDCOL={412:1,443:2,490:3,555:4,660:5,865:6}
def read_blend(chl,band):
    mie=f"{MIEDIR}/blend_chl{chl:g}.mie"
    lines=open(mie).read().splitlines()
    ext=sca=g=None
    for ln in lines[2:8]:
        p=ln.split()
        if p and abs(float(p[0])-band/1000)<1e-4:
            ext=float(p[5]); sca=float(p[6]); g=float(p[4]); break
    phs=[i for i,l in enumerate(lines) if 'TETA' in l]
    ph=phs[0]; end=phs[1] if len(phs)>1 else len(lines)
    th=[]; P=[]
    for ln in lines[ph+1:end]:
        p=ln.split()
        if len(p)>=7:
            try: th.append(float(p[0])); P.append(float(p[BANDCOL[band]]))
            except ValueError: pass
    th=np.array(th); P=np.array(P); o=np.argsort(th); th,P=th[o],P[o]
    return ext,sca,g,th,P
def back_int(th,p):
    t=np.radians(th); s=np.sin(t)
    norm=0.5*np.trapezoid(p*s,t); p=p/norm
    m=th>=90.0; return 0.5*np.trapezoid(p[m]*s[m],t[m])
def make_extdata(chl,band):
    # OSOAA용 extdata: blend_chl.mie를 그대로 변환
    out=f"{WORK}/extdata/chl{chl:g}_{band}.extdata"
    if not os.path.exists(out):
        os.makedirs(os.path.dirname(out),exist_ok=True)
        subprocess.run(["python3","/tmp/mkext.py",f"{chl:g}",str(band),out],capture_output=True,text=True)
    return out if os.path.exists(out) else None
def f2(x):
    try: return float(x.replace('D','E'))
    except: return None
def osoaa(chl,band,sza,ext_path,prof):
    res=f"{WORK}/osoaa/chl{chl:g}_{band}_{sza}"
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
    Lu=float(np.interp(0.0,V,I)); Ed=None
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
chl_f=float(sys.argv[2]) if len(sys.argv)>2 else None
rows=[]
for chl in CHL:
    if chl_f is not None and abs(chl-chl_f)>1e-9: continue
    for band in BANDS:
        ext,sca,g,th,P=read_blend(chl,band)
        a_part=ext-sca; b_part=sca
        A=a_part+A_W[band]; B=b_part+B_W[band]
        t=np.radians(th)
        pr=DD*0.75*(1+np.cos(t)**2)+(1-DD)
        norm=0.5*np.trapezoid(P*np.sin(t),t); Pn=P/norm
        blend=(b_part*Pn+B_W[band]*pr)/(b_part+B_W[band])
        bbb=back_int(th,blend)
        BB_tot=bbb*B
        BB_pass=min(BB_tot,0.49999*B)
        lut=f"{WORK}/p11_chl{chl:g}_{band}.csv"
        with open(lut,"w") as f:
            f.write("theta_deg,P11\n")
            for tt,pp in zip(th,blend): f.write(f"{tt:.4f},{pp:.6e}\n")
        ext_path=make_extdata(chl,band)
        prof=f"{WORK}/prof_chl{chl:g}_{band}.txt"
        with open(prof,"w") as f:
            f.write(f"# phyto chl={chl} {band}nm\n#\n#\n# depth a b\n#\n0.0 {a_part:.6f} {b_part:.6f}\n200.0 {a_part:.6f} {b_part:.6f}\n")
        for sza in SZA:
            if sza_f is not None and sza!=sza_f: continue
            res=osoaa(chl,band,sza,ext_path,prof)
            if not res: print(f"OSOAA FAIL chl{chl} {band} {sza}"); continue
            oR=osoaa_rrs(res); cR=run_ocrt(band,A,B,BB_pass,lut,sza)
            rows.append(dict(chl=chl,band=band,sza=sza,OCRT_rrs0minus=cR,OSOAA_rrs0minus=oR))
import statistics
if rows:
    new=not os.path.exists("cmp_item3_chl.csv")
    with open("cmp_item3_chl.csv","a",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        if new: w.writeheader()
        w.writerows(rows)
    dR=[100*(r['OCRT_rrs0minus']/r['OSOAA_rrs0minus']-1) for r in rows if r['OSOAA_rrs0minus']==r['OSOAA_rrs0minus'] and r['OSOAA_rrs0minus']>0]
    print(f"sza={sza_f}: rrs MAPE={statistics.mean(abs(x) for x in dR):.2f}% bias={statistics.mean(dR):+.2f}%")
    print(f"+{len(rows)} rows")
