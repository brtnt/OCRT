#!/usr/bin/env python3
"""OCRT-OSOAA 자동 교차 판정기 — 저장된 OSOAA 참조(runs/)와 신규 OCRT fullgrid를 대조.
사용: compare_osoaa.py <pkg_root> <ocrt_pure.csv> <ocrt_tsm.csv> <ocrt_bfo.csv>
임계: 수중 rrs |중앙|<=0.15%, Rrs |중앙|<=0.60%, 대기 rho_I 평균|.|<=0.50%.
주: 풍속 3 m/s 기준 참조(풍<3 비교 제외 규칙은 수중 전용)."""
import sys, os, csv, math
import numpy as np
PKG=os.path.abspath(sys.argv[1]); n=1.34
sys.path.insert(0, os.path.join(PKG,'02_OSOAA','OSOAA_OCRT_REFERENCE'))
import wl_read
f2=lambda x: float(x.replace('D','E'))
def eds(R):
    d={}
    for ln in open(f'{R}/Advanced_outputs/Flux.txt'):
        s=ln.split()
        if s and s[0] in ('26','27') and len(s)>4: d[s[0]]=f2(s[4])
    return d['26'], d['27']
def lvl27(R):
    V=[];I=[]
    for ln in open(f'{R}/Advanced_outputs/RESLUM_Adv_UP.txt'):
        s=ln.split()
        if len(s)>=7 and s[0]=='27':
            try: V.append(f2(s[2])); I.append(f2(s[4]))
            except Exception: pass
    V=np.array(V);I=np.array(I);o=np.argsort(V);V,I=V[o],I[o]
    k=np.concatenate([[True],np.diff(V)>1e-9]); return V[k],I[k]
def load(F):
    O={}
    for r in csv.DictReader(open(F)):
        if float(r['raa_deg'])<=180: O[(float(r['vza_deg']),float(r['raa_deg']))]=r
    return O
def water_case(tag, Ros, F, thr_rrs=0.15, thr_R=0.60):
    e26,e27=eds(Ros)
    wlp=f'{Ros}/wl.txt' if os.path.exists(f'{Ros}/wl.txt') else f'{Ros}/wl_dump.txt'
    wl=wl_read.WaterLeaving(wlp,Ros,e26); V,I27=lvl27(Ros); O=load(F)
    keys=[k for k in sorted(O) if k[1]==90.0 and k[0]>=5.0]
    thw=np.array([math.degrees(math.asin(math.sin(math.radians(k[0]))/n)) for k in keys])
    osr=np.interp(thw,V,I27)/e27
    ocr=np.array([float(O[k]['rrs_I']) for k in keys])
    ocR=np.array([float(O[k]['Rrs_I']) for k in keys])
    osR=np.array([wl.rrs_above(k[0],90.0)[0] for k in keys])
    m_r=float(np.median(100*(ocr/osr-1))); m_R=float(np.median(100*(ocR/osR-1)))
    ok = abs(m_r)<=thr_rrs and abs(m_R)<=thr_R
    print(f"[{tag}] rrs 중앙 {m_r:+.3f}% (임계 {thr_rrs}) | Rrs 중앙 {m_R:+.3f}% (임계 {thr_R}) -> {'PASS' if ok else 'FAIL'}")
    return ok
def atm_case(tag, base, F, thr=0.50):
    O={}
    for r in csv.DictReader(open(F)):
        O[(float(r['vza_deg']),float(r['raa_deg']))]=float(r['rho_I'])
    mu=math.cos(math.radians(40.0)); devs=[]
    for phi in (0,60,90,120,180):
        p=f'{base}_p{phi}/Standard_outputs/RESLUM_vsVZA.txt'
        if not os.path.exists(p): continue
        V=[];RF=[]
        for ln in open(p):
            s=ln.split()
            if len(s)==7:
                try: V.append(float(s[0])); RF.append(float(s[3]))
                except Exception: pass
        V=np.array(V);RF=np.array(RF);o=np.argsort(V);V,RF=V[o],RF[o]
        for side,pe in ((+1,phi),(-1,(phi+180)%360)):
            raa=(180-pe)%360
            if raa>180: continue
            for vza in np.arange(5,61,5.0):
                if (vza,float(raa)) in O:
                    devs.append(100*(O[(vza,float(raa))]/np.interp(side*vza,V,RF)-1))
    m=float(np.mean(np.abs(devs))); ok=m<=thr
    print(f"[{tag}] 대기 rho_I 평균|dev| {m:.3f}% (임계 {thr}) n={len(devs)} -> {'PASS' if ok else 'FAIL'}")
    return ok
ok = True
ok &= water_case('pure_d200', os.path.join(PKG,'05_VALIDATION','runs','nt640_osoaa_pw_d200'), sys.argv[2])
ok &= water_case('red_clay_tsm555', os.path.join(PKG,'05_VALIDATION','runs','nt640_osoaa_tsm555_fr631'), sys.argv[3])
ok &= atm_case('rayleigh_bfo', os.path.join(PKG,'05_VALIDATION','runs','atm_osoaa_bfo'), sys.argv[4])
print('CROSS_CHECK:', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
