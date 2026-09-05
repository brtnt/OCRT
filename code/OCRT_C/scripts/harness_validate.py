#!/usr/bin/env python3
"""OCRT↔OSOAA 검증 harness (matrix v1.08).
item #1: Black Fresnel ocean TOA — Rayleigh / aerosol / Ray+aerosol IQU.
OCRT 측 grid 생성 + OSOAA 명령 emit. 확장: items 2-6 (IOP 매핑 주석).
"""
import subprocess, itertools, csv, sys, os

OCRT = "./build/v2_solver_vk"
MIE  = "inputs/M80C.mie"          # #1 aerosol model (maritime M80C)
WIND = 3                           # OSOAA Cox-Munk 사전계산 wind {3,5,7,10}

BANDS = [412, 443, 490, 555, 660, 865]
# 밴드별 AOD = aot865 × ratio (M80C Nor_Ext_Co PCHIP, 2026-06-27 계산)
AOT_RATIO = {412:1.1309, 443:1.1152, 490:1.0937, 555:1.0691, 660:1.0409, 865:1.0000}
AOT865 = [0.05, 0.3, 1.0]

# 기하 grid (25 distinct; raa 축퇴는 LUT가 자동 처리, 후처리서 {0,45,90,135,180}만)
SZA = [0, 40, 80]
VZA_KEEP = [0, 30, 60]
RAA_KEEP = [0, 45, 90, 135, 180]

AER_FLAGS = ["--aer-l-max","80","--aer-h-km","2.0"]
COMMON    = ["--decouple-sunglint","--gas-column-h2o","0","--gas-column-o3","0","--gas-column-no2","0","--gas-column-o2","0","--gas-column-co2","0","--gas-column-ch4","0","--surface","coxmunk","--wind-speed",str(WIND),"--m-max","16","--lut-vza-step","30","--lut-vza-max","60","--lut-raa-step","45"]

def ocrt_cmd(atm, wl, aot, sza, out):
    c = [OCRT,"--sza",str(sza),"--vza","0","--raa","90","--wavelength",str(wl)] + COMMON
    if atm == "rayleigh":
        c += ["--pressure","1013.25","--n-layers","100"]
    elif atm == "aerosol":   # aerosol-only: n-layers 100 수렴 확정(<0.003% vs 400)
        c += ["--pressure","0","--mie",MIE,"--aod",f"{aot:.4f}","--n-layers","100"] + AER_FLAGS
    elif atm == "ray_aer":   # coupled: 400 필요(100은 0.48% off)
        c += ["--pressure","1013.25","--mie",MIE,"--aod",f"{aot:.4f}","--n-layers","400"] + AER_FLAGS
    c += ["--output-full-grid",out]
    return c

def run_grid(atm, wl, aot865):
    aot = aot865 * AOT_RATIO[wl] if atm != "rayleigh" else 0.0
    rows = []
    for sza in SZA:
        out = f"/tmp/h1_{atm}_{wl}_{aot865}_{sza}.csv"
        r = subprocess.run(ocrt_cmd(atm, wl, aot, sza, out), capture_output=True, text=True)
        if not os.path.exists(out):
            sys.stderr.write(f"FAIL {atm} wl{wl} aot{aot865} sza{sza}: {r.stderr[-200:]}\n"); continue
        seen=set()
        for d in csv.DictReader(open(out)):
            vza=round(float(d['vza_deg'])); raa=round(float(d['raa_deg']))
            if vza in VZA_KEEP and raa in RAA_KEEP:
                # raa 축퇴(sza=0 or vza=0): 대표 1개만
                key=(sza,vza, raa if (sza!=0 and vza!=0) else 0)
                if key in seen: continue
                seen.add(key)
                rows.append(dict(item=1, atm=atm, band=wl, aot865=aot865, aot_band=round(aot,4),
                                 sza=sza, vza=vza, raa=raa,
                                 rho_I=float(d['rho_I']), rho_Q=float(d['rho_Q']), rho_U=float(d['rho_U'])))
    return rows

if __name__ == "__main__":
    # 사용: harness_validate.py <atm> <band_group>
    #   atm: rayleigh | aerosol | ray_aer | pilot
    #   band_group: lo(412,443,490) | hi(555,660,865) | all
    # grid 파일에 append (timeout 회피용 chunk 실행). 각 chunk는 단일 command 내 완료.
    out="/home/claude/work/OCRT_item1_grid.csv"
    fields=["item","atm","band","aot865","aot_band","sza","vza","raa","rho_I","rho_Q","rho_U"]
    atm = sys.argv[1] if len(sys.argv)>1 else "pilot"
    bg  = sys.argv[2] if len(sys.argv)>2 else "all"
    band_set = {"lo":[412,443,490], "hi":[555,660,865], "all":BANDS}[bg]

    if atm == "pilot":
        combos=[("rayleigh",443,0.0),("aerosol",443,0.3),("ray_aer",443,0.3)]
    elif atm == "rayleigh":
        combos=[("rayleigh",wl,0.0) for wl in band_set]
    else:  # aerosol or ray_aer
        combos=[(atm,wl,a) for wl in band_set for a in AOT865]

    new = not os.path.exists(out)
    # resume: 이미 grid에 있는 (atm,band,aot865) combo skip
    done=set()
    if not new:
        for d in csv.DictReader(open(out)):
            done.add((d['atm'], int(d['band']), float(d['aot865'])))
    f=open(out,"a",newline="")
    w=csv.DictWriter(f, fieldnames=fields)
    if new: w.writeheader(); f.flush()
    total=0
    for at,wl,a in combos:
        if (at, wl, float(a)) in done:
            sys.stderr.write(f"skip {at} wl{wl} aot{a} (done)\n"); continue
        rs=run_grid(at,wl,a)
        w.writerows(rs); f.flush(); total+=len(rs)
        sys.stderr.write(f"done {at} wl{wl} aot{a}: +{len(rs)} (chunk total {total})\n"); sys.stderr.flush()
    f.close()
    print(f"chunk [{atm}/{bg}] -> +{total} rows in {out}")
