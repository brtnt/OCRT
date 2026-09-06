import sys, numpy as np
chl, band, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
MIEDIR="/home/user/ws/SPEEDUP_2026-07-09/OCRT_FULLPKG_v1.09_2026-07-08/aux/phase_mie"
BANDCOL={412:1,443:2,490:3,555:4,660:5,865:6}
lines=open(f"{MIEDIR}/blend_chl{chl}.mie").read().splitlines()
for ln in lines[2:8]:
    p=ln.split()
    if p and abs(float(p[0])-band/1000)<1e-4:
        ext=float(p[5]); sca=float(p[6]); break
# P11(band열), P12/P33 블록
blocks=[i for i,l in enumerate(lines) if 'Phase Function' in l or 'P11' in l or 'P12' in l or 'P33' in l]
# blend_chl.mie: P11 블록만 확인됨. P12/P33도 동일 구조 가정
hdrs=[i for i,l in enumerate(lines) if 'TETA' in l]
def parse(hi):
    th=[]; v=[]
    for ln in lines[hi+1:]:
        p=ln.split()
        if len(p)>=7:
            try: th.append(float(p[0])); v.append(float(p[BANDCOL[band]]))
            except ValueError: break
    return np.array(th), np.array(v)
th,P11=parse(hdrs[0])
P12=np.zeros_like(P11); P33=P11.copy()
if len(hdrs)>=2: _,P12=parse(hdrs[1])
if len(hdrs)>=3: _,P33=parse(hdrs[2])
o=np.argsort(th); th,P11,P12,P33=th[o],P11[o],P12[o],P33[o]
with open(out,'w') as f:
    f.write(f"EXTINCTION_COEF : {ext:.6E}\nSCATTERING_COEF : {sca:.6E}\nNB_LINES : {len(th)}\nANGLE  F11  -F12/F11  F22/F11  F33/F11\n")
    for i in range(len(th)):
        p12r=-P12[i]/P11[i] if P11[i]!=0 else 0.0
        p33r=P33[i]/P11[i] if P11[i]!=0 else 1.0
        f.write(f"   {th[i]:.3f}   {P11[i]:.6E}   {p12r:.6E}   1.000000E+00   {p33r:.6E}\n")
print(f"ext={ext} sca={sca} nang={len(th)} hdrs={len(hdrs)}")
