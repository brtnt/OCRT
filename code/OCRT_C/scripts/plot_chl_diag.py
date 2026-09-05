import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.font_manager as fm
_fp='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fm.fontManager.addfont(_fp); plt.rcParams['font.family']=fm.FontProperties(fname=_fp).get_name()
plt.rcParams['axes.unicode_minus']=False
import numpy as np
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.2))
# (A) backscatter fraction Bp
labels=['모멘트 진짜값\n(closed-form)','OSOAA\n(rrs 역산)','OCRT\n(rrs 역산)']
Bp=[0.02521,0.0246,0.0122]; col=['tab:gray','tab:blue','tab:red']
b=ax1.bar(labels,Bp,color=col,alpha=0.8,edgecolor='k')
for r,v in zip(b,Bp): ax1.text(r.get_x()+r.get_width()/2,v+0.0006,f'{v:.4f}',ha='center',fontsize=11)
ax1.axhline(0.02521,color='tab:gray',ls='--',lw=1)
ax1.set_ylabel('phyto backscatter fraction $B_p$ (절단 phase)')
ax1.set_title('(A) $B_p$: OSOAA는 모멘트 진짜값과 일치, OCRT는 절반\n(chl=1, 555nm)')
ax1.set_ylim(0,0.030); ax1.grid(True,axis='y',alpha=0.3)
# (B) rrs(0-)
labels2=['OSOAA\n(SOS)','모멘트 예측\n(Gordon)','OCRT\n(moment kernel)']
rrs=[4.012e-3,4.07e-3,2.616e-3]; col2=['tab:blue','tab:gray','tab:red']
b2=ax2.bar(labels2,rrs,color=col2,alpha=0.8,edgecolor='k')
for r,v in zip(b2,rrs): ax2.text(r.get_x()+r.get_width()/2,v+8e-5,f'{v:.3e}',ha='center',fontsize=10)
ax2.axhline(4.012e-3,color='tab:blue',ls='--',lw=1)
ax2.set_ylabel('rrs(0$-$) nadir [sr$^{-1}$]')
ax2.set_title('(B) rrs(0$-$): OCRT −34.8% (backscatter 과소 탓)\nIOP·ω*·truncation은 일치')
ax2.set_ylim(0,4.8e-3); ax2.grid(True,axis='y',alpha=0.3)
plt.suptitle('item#3 chl 진단: OCRT moment-kernel backscatter 2배 과소 버그 (OSOAA-PM phase 경로)',fontsize=13,y=1.0)
plt.tight_layout(); plt.savefig('chl_backscatter_diagnosis.png',dpi=130,bbox_inches='tight')
print('saved chl_backscatter_diagnosis.png')
