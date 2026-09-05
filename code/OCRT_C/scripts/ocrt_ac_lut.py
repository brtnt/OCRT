#!/usr/bin/env python3
"""
OCRT atmospheric-correction LUT driver (rigorous, coupled, SZA-swept).

Surface (argv[3]): black_fresnel | black   (ocean = separate heavier path, later)
Atmospheres: rayleigh (AOD=0) | aerosol | rayleigh_aerosol(both combined + Rayleigh baseline)
SZA sweep: 0,20,40,60,80.  Grid: VZA 0..80 step10 x RAA 0..330 step30.

Per (atm x model x aod865 x wavelength x sza x grid-point):
  rho_I/Q/U          TOA reflectance (atm+surface), rigorous SOS
  rho_rayleigh_I/Q/U Rayleigh-only TOA baseline (combined rows; reused from rayleigh case, ~0 cost)
  Rrs_I/Q/U          surface remote-sensing reflectance (=0 for black surfaces)
  Tdn_total          COUPLED downward total transmittance = T(sza)  [surface solve, rigorous]
  Ed_down            downward irradiance at surface = Tdn * mu_sun   (per unit F_sun)
  Tup_total          NaN for black/black_fresnel runs: no upward BOA source exists
  Tup_rt_valid        0 for these runs

Tdn is obtained from an explicit m=0 RT surface solve.  The former assignment
Tup(vza)=Tdn(theta=vza) by reciprocity has been removed: directional upward transmission
is source/BRDF dependent and must be taken from an ocean-coupled OCRT solve that reports
T_total_up_view with T_up_rt_valid=1.
"""
import csv, math, os, subprocess, sys, tempfile
TMPDIR = tempfile.gettempdir()

RUNDIR   = os.path.abspath(os.environ.get("OCRT_RUNDIR", "."))
MIE_DIR  = os.environ.get("OCRT_MIE_DIR", "inputs")
WIND     = float(os.environ.get("OCRT_WIND", "3.0"))

def _resolve_bin(b):
    # binary is invoked with cwd=RUNDIR; resolve to an absolute path, trying .exe on Windows.
    for c in ([b] if b.endswith(".exe") else [b, b+".exe"]):
        p = c if os.path.isabs(c) else os.path.join(RUNDIR, c)
        if os.path.exists(p): return os.path.abspath(p)
    return b if os.path.isabs(b) else os.path.abspath(os.path.join(RUNDIR, b))
OCRT_BIN = _resolve_bin(os.environ.get("OCRT_BIN", "./build/v2_solver"))

SURFACE   = sys.argv[3] if len(sys.argv) > 3 else "black_fresnel"
SURF_FLAG = {"black_fresnel": "coxmunk", "black": "black"}[SURFACE]
MODELS    = (sys.argv[2].split(",") if len(sys.argv) > 2 and sys.argv[2] else ["T50","C50","M80C"])

SZA_LIST     = [float(x) for x in os.environ.get("OCRT_SZA_LIST","0,20,40,60,80").split(",")]
WAVELENGTHS  = [float(x) for x in os.environ.get("OCRT_WLS","412,555,865").split(",")]
AOD865_LIST  = [float(x) for x in os.environ.get("OCRT_AODS","0.05,0.1,0.5,1.0").split(",")]
AER_L_MAX, M_MAX, AOD_REF_NM = 80, 16, 865.0

LUT_VZA_STEP, LUT_VZA_MAX, LUT_RAA_STEP = 10.0, 80.0, 30.0
VZA_NODES = [i*LUT_VZA_STEP for i in range(int(LUT_VZA_MAX/LUT_VZA_STEP)+1)]
C_VZA, C_RAA, C_RHO_I, C_RHO_Q, C_RHO_U, C_TIRRAD = 2, 3, 4, 5, 6, 8


def read_mie(model):
    path=os.path.join(RUNDIR,MIE_DIR,model+".mie"); rows=[]; seen=False
    for line in open(path):
        if "Wlgth" in line: seen=True; continue
        if not seen: continue
        p=line.split()
        if len(p)<7: break
        try: wl,ne=float(p[0]),float(p[1])
        except ValueError: break
        if not(0.2<=wl<=5.0): break
        rows.append((wl,ne))
    rows.sort(); return rows

def linterp(xs,ys,x):
    if x<=xs[0]: return ys[0]
    if x>=xs[-1]: return ys[-1]
    for i in range(1,len(xs)):
        if x<=xs[i]:
            t=(x-xs[i-1])/(xs[i]-xs[i-1]); return ys[i-1]+t*(ys[i]-ys[i-1])
    return ys[-1]

def aod_at(tbl,aod865,wl):
    xs=[r[0] for r in tbl]; ys=[r[1] for r in tbl]
    return aod865*linterp(xs,ys,wl/1000.0)/linterp(xs,ys,AOD_REF_NM/1000.0)

def _run(sza,wl,extra,tag,vstep,vmax,rstep):
    out=os.path.join(TMPDIR, f"ocrt_{tag}.csv")
    cmd=[OCRT_BIN,"--surface",SURF_FLAG,"--decouple-sunglint","--gas-column-h2o","0","--gas-column-o3","0","--gas-column-no2","0","--gas-column-o2","0","--gas-column-co2","0","--gas-column-ch4","0","--sza",str(sza),
         "--wavelength",str(wl),"--lut-vza-step",str(vstep),"--lut-vza-max",str(vmax),"--lut-raa-step",str(rstep),
         "--output-full-grid",out]
    if SURF_FLAG=="coxmunk": cmd+=["--wind-speed",str(WIND)]
    else: cmd+=["--wind-speed","0"]  # 2026-07-15 G4: wind 필수화 -- 종전 ocean 경로는 미지정(=0) 물리였으므로 0 명시
    cmd+=extra
    run_env=os.environ.copy()
    # The internal m=0 irradiance solve is an intentional convergence override.
    run_env.setdefault("OCRT_ADVANCED", "1")
    r=subprocess.run(cmd,cwd=RUNDIR,capture_output=True,text=True,env=run_env)
    if r.returncode!=0:
        sys.stderr.write(f"[FAIL] {tag} sza={sza}: {r.stderr.strip()[-160:]}\n"); return None
    path=out if os.path.isabs(out) else os.path.join(RUNDIR,out)
    return list(csv.reader(open(path)))[1:]

def _force_m0(extra):
    """Table needs only m=0 (irradiance_transmittance is the hemispheric=m0 mode). ~15x faster."""
    out=[]; i=0
    while i<len(extra):
        if extra[i]=="--m-max": i+=2; continue
        out.append(extra[i]); i+=1
    return out+["--m-max","0"]

def trans_table(extra,wl,tag):
    """Exact downward total transmittance at every requested SZA.

    The hemispheric irradiance uses only m=0, but each SZA is solved explicitly;
    no zenith-angle interpolation is used for the published Tdn_total values.
    """
    e0=_force_m0(extra)
    T={}
    for th in sorted(set(SZA_LIST)):
        rows=_run(th,wl,e0,f"{tag}_T{int(th)}",80.0,80.0,180.0)
        if rows is None: return None
        T[th]=float(rows[0][C_TIRRAD])
    return T

def refl_grid(extra,wl,sza,tag):
    rows=_run(sza,wl,extra,f"{tag}_R{int(sza)}",LUT_VZA_STEP,LUT_VZA_MAX,LUT_RAA_STEP)
    if rows is None: return None
    return {(round(float(v[C_VZA]),4),round(float(v[C_RAA]),4)):
            (float(v[C_RHO_I]),float(v[C_RHO_Q]),float(v[C_RHO_U])) for v in rows}

FIELDS=["surface","atm","model","aod865","wavelength","sza","vza","raa",
        "rho_I","rho_Q","rho_U","rho_rayleigh_I","rho_rayleigh_Q","rho_rayleigh_U",
        "Rrs_I","Rrs_Q","Rrs_U","Tdn_total","Ed_down",
        "Tup_total","Tup_rt_valid","Tup_method"]

def config_rows(atm,model,aod865,wl,extra,tag,ray_by_sza):
    T=trans_table(extra,wl,tag)
    if T is None: return [],None
    out=[]; my_ray={}
    for sza in SZA_LIST:
        grid=refl_grid(extra,wl,sza,tag)
        if grid is None: continue
        mu=math.cos(math.radians(sza))
        Tdn=T[sza]
        Ed=Tdn*mu
        if atm=="rayleigh": my_ray[sza]={}
        for (vza,raa),(rI,rQ,rU) in sorted(grid.items()):
            # No upward BOA source is prescribed in black/black_fresnel LUT solves.
            # Reporting T(vza) from the downward table was a reciprocity approximation
            # that ignored source angular distribution and BRDF.
            Tup=float("nan")
            if atm=="rayleigh_aerosol" and ray_by_sza:
                rr=ray_by_sza.get(sza,{}).get((vza,raa),("","",""))
            elif atm=="rayleigh":
                rr=(rI,rQ,rU); my_ray[sza][(vza,raa)]=(rI,rQ,rU)
            else: rr=("","","")
            out.append(dict(surface=SURFACE,atm=atm,model=model,aod865=aod865,wavelength=wl,sza=sza,
                vza=vza,raa=raa,rho_I=rI,rho_Q=rQ,rho_U=rU,
                rho_rayleigh_I=rr[0],rho_rayleigh_Q=rr[1],rho_rayleigh_U=rr[2],
                Rrs_I=0.0,Rrs_Q=0.0,Rrs_U=0.0,Tdn_total=Tdn,Ed_down=Ed,
                Tup_total=Tup,Tup_rt_valid=0,Tup_method="undefined_no_upward_BOA_source"))
    return out,(my_ray if atm=="rayleigh" else None)

def main():
    if not os.path.exists(OCRT_BIN):
        sys.stderr.write(f"ERROR: OCRT binary not found at {OCRT_BIN}\n"
                         f"  Build it first (run 'make' in this folder). Expected build/v2_solver(.exe).\n"
                         f"  Or set OCRT_BIN to the binary path.\n")
        sys.exit(1)
    sys.stderr.write(f"OCRT_BIN = {OCRT_BIN}\nRUNDIR   = {RUNDIR}\nsurface  = {SURFACE}  models={MODELS}\n")
    out_csv=os.path.abspath(sys.argv[1] if len(sys.argv)>1 else os.path.join(TMPDIR,"ac_lut.csv"))
    tables={m:read_mie(m) for m in MODELS}
    rows=[]; nsolve=0; NPER=2*len(SZA_LIST)
    ray_by_wl={}
    for wl in WAVELENGTHS:
        r,ray=config_rows("rayleigh","none",0.0,wl,["--pressure","1013.25","--aod","0"],f"ray{int(wl)}",None)
        rows+=r; nsolve+=NPER; ray_by_wl[wl]=ray
        sys.stderr.write(f"  rayleigh wl={wl}: {len(r)} rows\n")
    for atm,pres in [("aerosol","0"),("rayleigh_aerosol","1013.25")]:
        for model in MODELS:
            mie=os.path.join(MIE_DIR,model+".mie")
            for aod865 in AOD865_LIST:
                for wl in WAVELENGTHS:
                    al=aod_at(tables[model],aod865,wl)
                    # coupled 에어로졸 표준설정 3종 필수(누락 시 OSOAA 대비 검증 불성립):
                    #   --aer-h-km 2.0        에어로졸 scale height (OSOAA AP.HA 2.0 정합)
                    #   --n-layers 400        수직 세분(에어로졸 coupled 수렴)
                    #      전방peak 실각 log10-선형 truncation
                    extra=["--pressure",pres,"--mie",mie,"--aod",f"{al:.6f}",
                           "--aer-l-max",str(AER_L_MAX),"--m-max",str(M_MAX),
                           "--aer-h-km","2.0","--n-layers","400"]
                    rb=ray_by_wl.get(wl) if atm=="rayleigh_aerosol" else None
                    r,_=config_rows(atm,model,aod865,wl,extra,f"{atm[:3]}_{model}_{aod865}_{int(wl)}",rb)
                    rows+=r; nsolve+=NPER
            sys.stderr.write(f"  {atm} {model}: done\n")
    with open(out_csv,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    sys.stderr.write(f"TOTAL: {len(rows)} rows, ~{nsolve} solves, surface={SURFACE}, SZA={SZA_LIST} -> {out_csv}\n")

if __name__=="__main__":
    main()
