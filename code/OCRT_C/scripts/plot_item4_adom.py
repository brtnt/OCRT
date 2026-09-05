#!/usr/bin/env python3
"""item#4 aDOM 검증 그림."""
import csv, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
_fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fm.fontManager.addfont(_fp); plt.rcParams['font.family']=fm.FontProperties(fname=_fp).get_name()
plt.rcParams['axes.unicode_minus']=False
import numpy as np

rows=[r for r in csv.DictReader(open('cmp_item4_adom.csv'))]
def fl(r,k): return float(r[k])
ad   =np.array([fl(r,'adom440') for r in rows])
band =np.array([int(r['band']) for r in rows])
va   =np.array([int(r['vza_air']) for r in rows])
ocI  =np.array([fl(r,'OCRT_TOA_I') for r in rows]); osI=np.array([fl(r,'OSOAA_TOA_I') for r in rows])
ocR  =np.array([fl(r,'OCRT_rrs0minus') for r in rows]); osR=np.array([fl(r,'OSOAA_rrs0minus') for r in rows])
atot =np.array([fl(r,'OCRT_a_total') for r in rows])

fig,axes=plt.subplots(2,2,figsize=(13,10))
admk={0.01:'o',0.1:'s',1.0:'^'}; adcol={0.01:'tab:green',0.1:'tab:orange',1.0:'tab:red'}
bandcol={412:'tab:purple',443:'tab:blue',490:'tab:cyan',555:'tab:green',660:'tab:orange',865:'tab:red'}

# (A) Rrs scatter (log-log, 넓은 range 정당화)
ax=axes[0,0]
for a in [0.01,0.1,1.0]:
    m=(ad==a)&(osR>0)
    ax.scatter(osR[m],ocR[m],marker=admk[a],s=26,c=adcol[a],alpha=0.75,label=f'aDOM(440)={a:g}',edgecolors='none')
lim=[min(osR[osR>0].min(),ocR[ocR>0].min())*0.6,max(osR.max(),ocR.max())*1.5]
ax.plot(lim,lim,'k--',lw=1,alpha=0.6,label='1:1'); ax.set_xscale('log');ax.set_yscale('log')
ax.set_xlim(lim);ax.set_ylim(lim);ax.set_xlabel('OSOAA rrs(0$-$)');ax.set_ylabel('OCRT rrs(0$-$)')
ax.set_title('(A) rrs(0$-$) OCRT vs OSOAA (sza=40, 전 band/geometry)\n(log: 강흡수서 rrs 미소)')
ax.legend(fontsize=8,loc='upper left');ax.grid(True,alpha=0.3,which='both')

# (B) Rrs bias vs band, aDOM별 (nadir)
ax=axes[0,1]
bl=[412,443,490,555,660,865]
for a in [0.01,0.1,1.0]:
    bias=[]; 
    for b in bl:
        m=(ad==a)&(band==b)&(va==0)
        bias.append(100*(ocR[m][0]/osR[m][0]-1) if m.any() and osR[m][0] else np.nan)
    ax.plot(bl,bias,marker=admk[a],color=adcol[a],lw=1.8,ms=7,label=f'aDOM(440)={a:g}')
ax.axhline(0,color='k',lw=0.8);ax.axhspan(-1,1,color='gray',alpha=0.15)
ax.set_xlabel('band (nm)');ax.set_ylabel('rrs(0$-$) bias vs OSOAA [%]')
ax.set_title('(B) rrs(0$-$) bias by band & aDOM (nadir)\n흡수 클수록 음의 bias 증가')
ax.legend(fontsize=8);ax.grid(True,alpha=0.3)

# (C) 핵심: Rrs bias vs 총흡수 a_total (흡수 구동 잔차 검증)
ax=axes[1,0]
m=(va==0)&(osR>0)
biasC=100*(ocR/osR-1)
for b in bl:
    mm=m&(band==b)
    ax.scatter(atot[mm],biasC[mm],c=bandcol[b],s=40,label=f'{b} nm',edgecolors='k',linewidths=0.3,zorder=3)
ax.axhline(0,color='k',lw=0.8);ax.axhline(-1,color='gray',ls=':');
ax.set_xscale('log');ax.set_xlabel('OCRT 총흡수 $a_{total}$ [m$^{-1}$] (log)')
ax.set_ylabel('rrs(0$-$) bias vs OSOAA [%]')
ax.set_title('(C) rrs(0$-$) bias vs 총흡수 (nadir, 전 aDOM/band)\n한 곡선으로 모임 = 흡수 구동 잔차(constituent 무관)')
ax.legend(fontsize=7,ncol=2);ax.grid(True,alpha=0.3,which='both')

# (D) TOA scatter
ax=axes[1,1]
for a in [0.01,0.1,1.0]:
    mm=(ad==a)
    ax.scatter(osI[mm],ocI[mm],marker=admk[a],s=26,c=adcol[a],alpha=0.75,label=f'aDOM={a:g}',edgecolors='none')
lim2=[0,max(osI.max(),ocI.max())*1.1];ax.plot(lim2,lim2,'k--',lw=1,alpha=0.6)
ax.set_xlim(lim2);ax.set_ylim(lim2);ax.set_xlabel('OSOAA TOA $\\rho_I$');ax.set_ylabel('OCRT TOA $\\rho_I$')
mape_toa=np.mean(np.abs(ocI/osI-1))*100
ax.set_title(f'(D) TOA $\\rho_I$ OCRT vs OSOAA (sza=40)\nMAPE={mape_toa:.2f}% (전 aDOM)')
ax.legend(fontsize=8,loc='upper left');ax.grid(True,alpha=0.3)

plt.suptitle('item#4 aDOM(CDOM/yellow substance) 검증: OCRT v1.08 vs OSOAA (sza=40, slope=0.014)',fontsize=13,y=1.0)
plt.tight_layout();plt.savefig('item4_adom_validation.png',dpi=130,bbox_inches='tight')
print('saved item4_adom_validation.png')
print(f'sza40 전체: TOA MAPE={mape_toa:.2f}%, Rrs MAPE={np.mean(np.abs(ocR/osR-1))*100:.2f}%')
