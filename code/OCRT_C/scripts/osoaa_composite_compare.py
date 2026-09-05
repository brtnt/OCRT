#!/usr/bin/env python3
"""item#6 복합(TSM+aDOM+chl, Rayleigh+aerosol): TOA IQU + rrs(0-).
3성분 IOP 합산 + 산란가중 블렌드 위상(§C). OCRT fixed-bulk + Rayleigh atm + M80C aerosol,
OSOAA HYD.Model 3(ExtData=합성위상, UserProfile=합성 a,b) + AER 동일. 인자: <sza> [case]"""
import subprocess, os, csv, math, sys, re
import numpy as np
from scipy.interpolate import PchipInterpolator as PCHIP
ROOT="/home/claude/ws/SPEEDUP_2026-07-09/osoaa"; SURF=f"{ROOT}/DATABASE/SURF_MATR"
OCRT="/home/claude/ws/SPEEDUP_2026-07-09/ocrt_opt/build/v2o_fbaer"
MIEDIR="/home/claude/ws/SPEEDUP_2026-07-09/OCRT_FULLPKG_v1.09_2026-07-08/aux/phase_mie"
BANDS=[412,443,490,555,660,865]
TAU={412:0.318540221,443:0.236054530,490:0.155974381,555:0.093751620,660:0.046362496,865:0.015540855}
A_W={412:4.5506e-3,443:7.069e-3,490:1.500e-2,555:5.960e-2,660:4.100e-1,865:4.605}
B_W={412:6.650e-3,443:4.872e-3,490:3.164e-3,555:1.859e-3,660:8.875e-4,865:2.763e-4}
# 복합 케이스: (chl mg/m3, TSM g/m3, aCDOM440 m-1)  — 연안/외양/혼탁
CASES={"oceanic":(0.3,0.1,0.01), "coastal":(3.0,1.0,0.1), "turbid":(30.0,10.0,1.0)}
SZA=[0,40,80]; VZA=[0,30,60]; AOT865=0.1; WORK="/tmp/it6"; os.makedirs(WORK,exist_ok=True)
DEL=0.039; DD=(1-DEL)/(1+DEL/2); BANDCOL={412:1,443:2,490:3,555:4,660:5,865:6}
NW=1.34
def _ext(b):
    for ln in open(f"/tmp/M80C_{b}.extdata"):
        if "EXTINCTION_COEF" in ln: return float(ln.split()[-1])
RATIO={b:_ext(b)/_ext(865) for b in BANDS}
def read_mie(name, band, chl_style):
    lines=open(f"{MIEDIR}/{name}.mie").read().splitlines()
    ext=sca=None
    for ln in lines[2:9]:
        p=ln.split()
        if p and len(p)>=7:
            try:
                if abs(float(p[0])-band/1000)<1e-4: ext=float(p[5]); sca=float(p[6]); break
            except ValueError: pass
    hdrs=[i for i,l in enumerate(lines) if 'TETA' in l]
    th=[]; P=[]
    end=hdrs[1] if len(hdrs)>1 else len(lines)
    col=BANDCOL[band] if chl_style else 1
    for ln in lines[hdrs[0]+1:end]:
        p=ln.split()
        if len(p)>=7:
            try: th.append(float(p[0])); P.append(float(p[col]))
            except ValueError: pass
    th=np.array(th); P=np.array(P); o=np.argsort(th); th,P=th[o],P[o]
    t=np.radians(th); nm=0.5*np.trapezoid(P*np.sin(t),t)
    return ext,sca,th,P/nm
def f2(x):
    try: return float(x.replace('D','E'))
    except: return None
def osoaa(case,band,sza,ext_path,prof,aer_ext):
    res=f"{WORK}/osoaa/{case}_{band}_{sza}"
    if os.path.exists(f"{res}/Advanced_outputs/Flux.txt"): return res
    os.makedirs(res,exist_ok=True)
    cmd=[f"{ROOT}/exe/OSOAA_MAIN.exe","-OSOAA.ResRoot",res,"-OSOAA.Log","M.Log",
      "-OSOAA.Wa",f"{band/1000:.5f}","-ANG.Thetas",f"{sza}.","-ANG.Rad.NbGauss","48","-ANG.Mie.NbGauss","100",
      "-AP.MOT",f"{TAU[band]}","-AP.HR","8.0",
      "-AER.Model","4","-AER.ExtData",aer_ext,"-AER.DirMie","/tmp/mie_dir","-AER.Tronca","1",
      "-AER.Waref",f"{band/1000:.5f}","-AER.AOTref",f"{AOT865*RATIO[band]:.5f}","-AP.HA","2.0",
      "-HYD.Model","3","-HYD.ExtData",ext_path,"-HYD.UserProfile",prof,
      "-PHYTO.Chl","0.0","-SED.Csed","0.0","-YS.Abs440","0.0","-DET.Abs440","0.0",
      "-SEA.Depth","200.0","-SEA.Ind","1.34","-SEA.Wind","3","-SEA.Dir",SURF,
      "-SEA.SurfAlb","0.0","-SEA.BotType","1","-SEA.BotAlb","0.0",
      "-OSOAA.View.Phi","90.","-OSOAA.View.Level","1",
      "-OSOAA.ResFile.vsVZA","RESLUM_vsVZA.txt","-OSOAA.ResFile.Adv.Up","RESLUM_Adv_UP.txt"]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=250,env={**os.environ,"OSOAA_NO_DIRECT_GLINT":"1","OSOAA_ROOT":ROOT,"HOME":"/tmp"})
    if not os.path.exists(f"{res}/Advanced_outputs/Flux.txt"):
        print("  OSOAA err:", (r.stdout or "")[-160:].replace("\n"," "))
        return None
    return res
def parse_osoaa(res,sza):
    mu=math.cos(math.radians(sza))
    vz,refl=[],[]
    for ln in open(f"{res}/Standard_outputs/RESLUM_vsVZA.txt"):
        s=ln.split()
        if len(s)==7:
            v=f2(s[0]); R=f2(s[3])
            if v is not None and R is not None and v>=0: vz.append(v); refl.append(R)
    vz,refl=np.array(vz),np.array(refl); o=np.argsort(vz); vz,refl=vz[o],refl[o]
    vz,iu=np.unique(vz,return_index=True); fR=PCHIP(vz,refl[iu])
    def adv(level):
        V,I,Q,U=[],[],[],[]
        for ln in open(f"{res}/Advanced_outputs/RESLUM_Adv_UP.txt"):
            s=ln.split()
            if len(s)>=7 and s[0]==str(level):
                v,ii,q,u=f2(s[2]),f2(s[4]),f2(s[5]),f2(s[6])
                if None in (v,ii,q,u) or v<0: continue
                V.append(v); I.append(ii); Q.append(q); U.append(u)
        V=np.array(V); o=np.argsort(V); V=V[o]
        I=np.array(I)[o]; Q=np.array(Q)[o]; U=np.array(U)[o]
        V,iu=np.unique(V,return_index=True)
        return PCHIP(V,I[iu]),PCHIP(V,Q[iu]),PCHIP(V,U[iu])
    fI0,fQ0,fU0=adv(0); fI27,_,_=adv(27)
    Ed=None
    for ln in open(f"{res}/Advanced_outputs/Flux.txt"):
        s=ln.split()
        if s and s[0]=="27" and len(s)>4: Ed=f2(s[4]); break
    return fR,fQ0,fU0,mu,fI27,Ed
sza_f=int(sys.argv[1]) if len(sys.argv)>1 else None
case_f=sys.argv[2] if len(sys.argv)>2 else None
DONE=set()
if os.path.exists("cmp_item6_composite.csv"):
    for r in csv.DictReader(open("cmp_item6_composite.csv")):
        DONE.add((r['case'],int(r['band']),int(r['sza']),int(float(r['vza']))))
rows=[]
for case,(chl,tsm,acdom) in CASES.items():
    if case_f and case!=case_f: continue
    for band in BANDS:
        e_c,s_c,th,P_chl=read_mie(f"blend_chl{chl:g}",band,True)
        # 광물: #5에서 검증된 extdata(변환기 산출) 재사용 — .mie 원본은 1nm 전파장표라 헤더 구조가 다름
        e_m=s_m=None; thm=[]; Pm=[]
        for ln in open(f"/tmp/tsm/extdata/{band}.extdata"):
            sp=ln.split()
            if "EXTINCTION_COEF" in ln: e_m=float(sp[-1])
            elif "SCATTERING_COEF" in ln: s_m=float(sp[-1])
            elif len(sp)>=2:
                try: thm.append(float(sp[0])); Pm.append(float(sp[1]))
                except ValueError: pass
        thm=np.array(thm); Pm=np.array(Pm)
        tm=np.radians(thm); nmm=0.5*np.trapezoid(Pm*np.sin(tm),tm)
        P_min=np.interp(th, thm, Pm/nmm)
        a_chl=(e_c-s_c); b_chl=s_c
        a_min=tsm*(e_m-s_m); b_min=tsm*s_m
        a_cdom=acdom*math.exp(-0.014*(band-440.0))
        A=a_chl+a_min+a_cdom+A_W[band]
        B=b_chl+b_min+B_W[band]
        t=np.radians(th); pr=DD*0.75*(1+np.cos(t)**2)+(1-DD)
        blend=(b_chl*P_chl+b_min*P_min+B_W[band]*pr)/B
        nb=0.5*np.trapezoid(blend*np.sin(t),t); bl=blend/nb
        m=th>=90; bbb=0.5*np.trapezoid(bl[m]*np.sin(t)[m],t[m])
        BB=min(bbb*B,0.49999*B)
        lut=f"{WORK}/p11_{case}_{band}.csv"
        with open(lut,"w") as f:
            f.write("theta_deg,P11\n")
            for tt,pp in zip(th,blend): f.write(f"{tt:.4f},{pp:.6e}\n")
        # OSOAA hydrosol ExtData: 입자 전용 위상이어야 한다.
        # OSOAA는 순수해수 분자산란을 자체적으로(레일리 위상으로) 더한다. 여기에 물 레일리를
        # 섞어 넣으면 물 후방산란이 이중 계상되어 OSOAA rrs가 과대해진다(실측 -19%@412 ~ -3.5%@865,
        # b_w/b_tot에 비례하는 파장 구조로 확인). OCRT fixed-bulk는 단일매질이라 반대로 물 레일리를
        # 포함한 블렌드가 맞다(§B/§C) — 두 코드의 규약이 다르다.
        b_part_=b_chl+b_min
        blend_p=(b_chl*P_chl+b_min*P_min)/b_part_
        nbp=0.5*np.trapezoid(blend_p*np.sin(t),t); blp=blend_p/nbp
        ext_path=f"{WORK}/hyd_{case}_{band}.extdata"
        with open(ext_path,"w") as f:
            f.write(f"EXTINCTION_COEF : {(a_chl+a_min+a_cdom)+b_part_:.6E}\nSCATTERING_COEF : {b_part_:.6E}\nNB_LINES : {len(th)}\nANGLE  F11  -F12/F11  F22/F11  F33/F11\n")
            for tt,pp in zip(th,blp):
                f.write(f"   {tt:.3f}   {pp:.6E}   0.000000E+00   1.000000E+00   1.000000E+00\n")
        prof=f"{WORK}/prof_{case}_{band}.txt"
        a_part=a_chl+a_min+a_cdom; b_part=b_chl+b_min
        with open(prof,"w") as f:
            f.write(f"# composite {case} {band}nm\n#\n#\n# depth a b\n#\n0.0 {a_part:.6f} {b_part:.6f}\n200.0 {a_part:.6f} {b_part:.6f}\n")
        aer_ext=f"/tmp/M80C_{band}.extdata"
        for sza in SZA:
            if sza_f is not None and sza!=sza_f: continue
            res=osoaa(case,band,sza,ext_path,prof,aer_ext)
            if not res: print(f"OSOAA FAIL {case} {band} {sza}"); continue
            fR,fQ0,fU0,mu,fI27,Ed=parse_osoaa(res,sza)
            for va in VZA:
                if (case,band,sza,va) in DONE: continue
                vw=math.degrees(math.asin(min(1.0,math.sin(math.radians(va))/NW)))
                oI=float(fR(va)); oQ=float(fQ0(va))/mu; oU=float(fU0(va))/mu
                oRrs=float(fI27(vw))/Ed if Ed else float('nan')
                r=subprocess.run([OCRT,"--surface","ocean","--decouple-sunglint","--gas-column-h2o","0","--gas-column-o3","0","--gas-column-no2","0","--gas-column-o2","0","--gas-column-co2","0","--gas-column-ch4","0","--wind-speed","3","--sza",str(sza),"--vza",str(va),
                    "--raa","90","--wavelength",str(band),"--pressure","1013.25",
                    "--mie","inputs/M80C.mie","--aod",f"{AOT865*RATIO[band]:.4f}","--n-layers","400",
                    "--aer-h-km","2.0","--aer-l-max","80",
                    "--water-model","iop","--iop-a",f"{A:.6f}","--iop-b",f"{B:.6f}","--iop-bb",f"{BB:.8f}","--iop-phase-lut",lut,
                    "--n-mu-water","48"],capture_output=True,text=True,timeout=250,
                    env={**os.environ,"OCRT_ADVANCED":"1","OMP_NUM_THREADS":"1"})
                out=r.stdout.strip().split("\n")[-1]; tk=out.split()
                try: cI,cQ,cU=float(tk[0]),float(tk[1]),float(tk[2])
                except (ValueError,IndexError): print(f"OCRT FAIL {case} {band} {sza} {va}"); continue
                mm=re.search(r"rrs0minus=([0-9.eE+-]+)",out)
                cRrs=float(mm.group(1)) if mm else float('nan')
                row=dict(case=case,chl=chl,tsm=tsm,acdom=acdom,band=band,sza=sza,vza=va,
                    OCRT_TOA_I=cI,OSOAA_TOA_I=oI,OCRT_TOA_Q=cQ,OSOAA_TOA_Q=oQ,
                    OCRT_TOA_U=cU,OSOAA_TOA_U=oU,OCRT_rrs0minus=cRrs,OSOAA_rrs0minus=oRrs)
                rows.append(row)
                newf=not os.path.exists("cmp_item6_composite.csv")
                with open("cmp_item6_composite.csv","a",newline="") as fo:   # 셀 단위 즉시 기록(도구 시간 한계 대비)
                    w=csv.DictWriter(fo,fieldnames=list(row.keys()))
                    if newf: w.writeheader()
                    w.writerow(row)
import statistics
if rows:
    dI=[100*(r['OCRT_TOA_I']/r['OSOAA_TOA_I']-1) for r in rows]
    dR=[100*(r['OCRT_rrs0minus']/r['OSOAA_rrs0minus']-1) for r in rows if r['OSOAA_rrs0minus']==r['OSOAA_rrs0minus'] and r['OSOAA_rrs0minus']>0]
    print(f"sza={sza_f} case={case_f}: TOA_I MAPE={statistics.mean(abs(x) for x in dI):.2f}% | Rrs bias={statistics.mean(dR):+.2f}%")
    print(f"+{len(rows)} rows")
