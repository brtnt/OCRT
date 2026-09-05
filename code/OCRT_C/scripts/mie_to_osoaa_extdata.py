#!/usr/bin/env python3
"""mie_to_osoaa_extdata.py — .mie(OPAC) → OSOAA AER.ExtData(FICEXTDATA, IMOD=4).
포맷: EXTINCTION_COEF / SCATTERING_COEF / NB_LINES / 주석 / ANGLE F11 -F12/F11 F22/F11 F33/F11 (오름차순).
구면 aerosol: F22/F11=1, F33/F11=P33/P11. 부호: col3 = -P12/P11 (OS 규약, F12<0 for Rayleigh).
사용: mie_to_osoaa_extdata.py <mie> <band_nm> <out>
"""
import sys, numpy as np
from scipy.interpolate import PchipInterpolator as PCHIP

mie, band, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
wl_um = band/1000.0
lines = open(mie, encoding='latin1').read().splitlines()

# --- 1) 추출/산란 계수 테이블: 파장행에서 Extinct_Co(col6), Scatter_Co(col7) ---
# 테이블 행 예: "0.4120  1.0560  1.0556  0.9963  0.7708  5.0667E-02  5.0480E-02"
tbl_wl=[]; tbl_ext=[]; tbl_sca=[]
for ln in lines:
    s=ln.split()
    if len(s)==7:
        try:
            w=float(s[0]); e=float(s[5]); sc=float(s[6])
            if 0.2<w<5.0: tbl_wl.append(w); tbl_ext.append(e); tbl_sca.append(sc)
        except: pass
tbl_wl=np.array(tbl_wl); 
ext=float(PCHIP(tbl_wl,np.array(tbl_ext))(wl_um))
sca=float(PCHIP(tbl_wl,np.array(tbl_sca))(wl_um))

# --- 2) phase 블록 3개 (P11/P12/P33) 파싱 ---
def find_blocks():
    hdr=[i for i,l in enumerate(lines) if 'Phase Function' in l]
    return hdr  # [P11_hdr, P12_hdr, P33_hdr]
hdr=find_blocks()
assert len(hdr)>=3, f"phase 블록 3개 필요, 찾음 {len(hdr)}"

def parse_block(hi):
    # hi: 'Phase Function' 헤더 줄 index. 다음 줄 = TETA 헤더(파장 리스트), 그 다음부터 데이터.
    teta=lines[hi+1].split()           # ['TETA','0.3500',...]
    wls=[float(x) for x in teta[1:]]
    ang=[]; rows=[]
    j=hi+2
    while j<len(lines):
        s=lines[j].split()
        if not s or 'Phase Function' in lines[j]: break
        try:
            a=float(s[0]); vals=[float(x) for x in s[1:1+len(wls)]]
            if len(vals)==len(wls): ang.append(a); rows.append(vals)
        except: break
        j+=1
    ang=np.array(ang); rows=np.array(rows)   # [n_ang, n_wl]
    # band 파장으로 보간 (각도별 wl-PCHIP). 정확매치면 그 컬럼.
    wls=np.array(wls)
    if np.any(np.abs(wls-wl_um)<1e-6):
        col=int(np.argmin(np.abs(wls-wl_um))); v=rows[:,col]
    else:
        v=np.array([float(PCHIP(wls,rows[k,:])(wl_um)) for k in range(len(ang))])
    return ang, v

aP,P11=parse_block(hdr[0])
_,P12=parse_block(hdr[1])
_,P33=parse_block(hdr[2])

# --- 3) 오름차순 정렬 + ExtData 작성 ---
o=np.argsort(aP)
ang=aP[o]; P11=P11[o]; P12=P12[o]; P33=P33[o]
with open(out,'w') as f:
    f.write(f"EXTINCTION_COEF : {ext:.6E}\n")
    f.write(f"SCATTERING_COEF : {sca:.6E}\n")
    f.write(f"NB_LINES : {len(ang)}\n")
    f.write("ANGLE  F11  -F12/F11  F22/F11  F33/F11\n")
    for a,p11,p12,p33 in zip(ang,P11,P12,P33):
        f.write(f"{a:8.3f} {p11:14.6E} {(-p12/p11):14.6E} {1.0:14.6E} {(p33/p11):14.6E}\n")
print(f"{out}: band={band}nm ext={ext:.4e} sca={sca:.4e} ssa={sca/ext:.5f} nang={len(ang)} ang[{ang[0]:.1f}→{ang[-1]:.1f}]")
