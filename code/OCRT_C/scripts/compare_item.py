#!/usr/bin/env python3
"""OCRT grid ↔ OSOAA grid 비교 (turnkey). 하네스 §1 규약 적용.
사용: python compare_item.py OCRT_item1_grid.csv OSOAA_item1_grid.csv [item#]
- 매칭 키: (atm,band,aot865,sza,vza,raa)
- RAA 규약(§1.0): RAA=90 및 sza=0/vza=0(축퇴)만 직접 비교. 그 외 raa는 Θ검증 필요 → 'flagged'.
- rho_I primary(MAPE), rho_Q/U는 RMS/I.
- scatter(linear) PNG 생성.
"""
import csv, sys, math, statistics
import matplotlib; matplotlib.use("Agg")
from matplotlib import font_manager as fm
_fp="/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
try: fm.fontManager.addfont(_fp); matplotlib.rcParams["font.family"]=fm.FontProperties(fname=_fp).get_name()
except Exception: pass
matplotlib.rcParams["axes.unicode_minus"]=False
import matplotlib.pyplot as plt

ocrt_f, osoaa_f = sys.argv[1], sys.argv[2]
item = sys.argv[3] if len(sys.argv)>3 else "1"

def load(f):
    d={}
    for r in csv.DictReader(open(f)):
        k=(r['atm'],r['band'],r['aot865'],r['sza'],r['vza'],r['raa'])
        d[k]=(float(r['rho_I']),float(r['rho_Q']),float(r['rho_U']))
    return d
O=load(ocrt_f); S=load(osoaa_f)
keys=sorted(set(O)&set(S))
print(f"matched {len(keys)} cells (OCRT {len(O)}, OSOAA {len(S)})")

def direct(k):  # 직접 비교 가능 기하? (RAA=90 또는 sza0/vza0 축퇴)
    atm,band,a,sza,vza,raa=k
    return raa=='90' or sza=='0' or vza=='0'

# per (atm,band) MAPE on rho_I, direct 기하만
from collections import defaultdict
agg=defaultdict(list); flagged=0; allpts=[]
for k in keys:
    oi,oq,ou=O[k]; si,sq,su=S[k]
    if si==0: continue
    d=100*(oi/si-1)
    if direct(k):
        agg[(k[0],k[1])].append(d); allpts.append((si,oi,k[0]))
    else:
        flagged+=1
print(f"직접비교 {sum(len(v) for v in agg.values())} cells, Θ검증필요(non-90 raa) {flagged} cells\n")
print(f"{'atm':>9} {'band':>5} {'N':>3} {'MAPE_I%':>8} {'maxd%':>7}")
overall=[]
for (atm,band),ds in sorted(agg.items()):
    m=statistics.mean(abs(x) for x in ds); mx=max(abs(x) for x in ds)
    overall+=ds
    print(f"{atm:>9} {band:>5} {len(ds):>3} {m:>8.2f} {mx:>7.1f}")
if overall:
    print(f"\n전체 직접비교 MAPE_I = {statistics.mean(abs(x) for x in overall):.2f}%  "
          f"bias = {statistics.mean(overall):+.2f}%  (목표 ≤1%)")

# scatter (linear), atm 색
fig,ax=plt.subplots(1,2,figsize=(13,5.5))
fig.suptitle(f"OCRT ↔ OSOAA item #{item} (Black Fresnel ocean TOA, 직접비교 기하)",fontweight="bold")
fig.text(0.5,0.01,"약어: rho_I=reflectance I(=π·L/Ed) | 직접비교=RAA90+sza0/vza0 축퇴 | 점선=1:1 | atm: Rayleigh/aerosol/Ray+aerosol",
         ha="center",fontsize=8,bbox=dict(boxstyle="round",fc="#f4f4f4",ec="#999"))
col={'rayleigh':'#1f77b4','aerosol':'#d62728','ray_aer':'#2ca02c'}
for a,c in col.items():
    pts=[(s,o) for s,o,at in allpts if at==a]
    if pts: ax[0].scatter([p[0] for p in pts],[p[1] for p in pts],s=18,c=c,alpha=0.6,label=a,edgecolor='none')
if allpts:
    lo=0; hi=max(max(p[0],p[1]) for p in allpts)*1.05
    ax[0].plot([lo,hi],[lo,hi],'k--',lw=1)
ax[0].set_xlabel("OSOAA rho_I"); ax[0].set_ylabel("OCRT rho_I"); ax[0].set_title("(A) 1:1 scatter")
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)
# d% 분포 (atm별)
for a,c in col.items():
    ds=agg_all=[100*(O[k][0]/S[k][0]-1) for k in keys if direct(k) and k[0]==a and S[k][0]!=0]
    if ds: ax[1].scatter([float(k_) for k_ in range(len(ds))],ds,s=14,c=c,alpha=0.5,label=a,edgecolor='none')
ax[1].axhline(0,color='k',lw=0.8); ax[1].axhline(1,color='gray',ls=':'); ax[1].axhline(-1,color='gray',ls=':')
ax[1].set_xlabel("cell index"); ax[1].set_ylabel("d% = 100·(OCRT/OSOAA−1)"); ax[1].set_title("(B) d% 분포 (점선 ±1%)")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
plt.tight_layout(rect=[0,0.04,1,0.96])
out=f"compare_item{item}.png"; plt.savefig(out,dpi=140,bbox_inches="tight")
print(f"\nscatter -> {out}")
