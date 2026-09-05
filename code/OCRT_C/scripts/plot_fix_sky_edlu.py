#!/usr/bin/env python3
"""FIX-SKY-EDLU 검증 그림: skylight Ed/Lu sub-cone leak 수정 전후."""
import csv, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
_fp = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fm.fontManager.addfont(_fp)
plt.rcParams['font.family'] = fm.FontProperties(fname=_fp).get_name()
plt.rcParams['axes.unicode_minus'] = False
import numpy as np

# --- sza40 #2 결과 (B 후) ---
rows = [r for r in csv.DictReader(open('cmp_item2_pureocean.csv'))]
oc_rrs = np.array([float(r['OCRT_rrs0minus']) for r in rows])
os_rrs = np.array([float(r['OSOAA_rrs0minus']) for r in rows])
oc_toa = np.array([float(r['OCRT_TOA_I']) for r in rows])
os_toa = np.array([float(r['OSOAA_TOA_I']) for r in rows])
bands  = [r['band'] for r in rows]

# --- 측정값 (하드코딩: 본 세션 런) ---
winds = [0, 3, 12]
drop_before_40 = [-0.23, -0.93, -1.52]   # sza40 원본(버그)
drop_after_40  = [-0.23, -0.14, -0.05]   # sza40 B후
drop_after_80  = [+3.71, +3.59, +3.71]   # sza80 B후 (wind-무관=물리)

band_list = ['412','443','490','555','660','865']
bias_before = {'412':+3.01,'443':-0.02,'490':-1.85,'555':-1.98,'660':-2.82,'865':-2.19}
bias_after  = {'412':+4.27,'443':+1.78,'490':+0.07,'555':-0.30,'660':-1.55,'865':-1.52}

# sub-cone 비중 vs sza (diffuse field에서)
sza_sub = [40, 80]
subfrac = [None, 25.1]  # sza40는 아래서 계산해 채움(없으면 NA)

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# ===== Panel A: OCRT-B vs OSOAA rrs(0-) scatter (log-log, 4桁 range 정당화) =====
ax = axes[0,0]
cmap = {'412':'tab:purple','443':'tab:blue','490':'tab:cyan','555':'tab:green','660':'tab:orange','865':'tab:red'}
for b in band_list:
    m = [i for i in range(len(rows)) if bands[i]==b]
    ax.scatter(os_rrs[m], oc_rrs[m], s=28, c=cmap[b], label=f'{b} nm', alpha=0.8, edgecolors='none')
lim = [min(os_rrs.min(),oc_rrs.min())*0.7, max(os_rrs.max(),oc_rrs.max())*1.4]
ax.plot(lim, lim, 'k--', lw=1, alpha=0.6, label='1:1')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('OSOAA rrs(0$-$)'); ax.set_ylabel('OCRT-B rrs(0$-$)')
ax.set_title('(A) rrs(0$-$) OCRT-B vs OSOAA, sza=40 전 geometry\n(log: rrs 4桁 dynamic range 정당화)')
ax.legend(fontsize=7, loc='upper left'); ax.grid(True, alpha=0.3, which='both')

# ===== Panel B: bias(%) by band 전/후 (선형 — A의 동반) =====
ax = axes[0,1]
x = np.arange(len(band_list)); w=0.38
ax.bar(x-w/2, [bias_before[b] for b in band_list], w, label='수정 전(버그)', color='tab:red', alpha=0.7)
ax.bar(x+w/2, [bias_after[b]  for b in band_list], w, label='수정 후(B)', color='tab:green', alpha=0.7)
ax.axhline(0, color='k', lw=0.8); ax.axhspan(-1,1, color='gray', alpha=0.15, label='±1%')
ax.set_xticks(x); ax.set_xticklabels([f'{b}' for b in band_list])
ax.set_xlabel('band (nm)'); ax.set_ylabel('rrs(0$-$) bias vs OSOAA [%]')
ax.set_title('(B) rrs(0$-$) bias by band (sza=40, nadir)\n490/555 버그 해소; 412/660/865은 노출된 별개 잔차')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis='y')

# ===== Panel C: rrs atm 민감도 wind 의존성 전/후 (선형 — 버그 signature) =====
ax = axes[1,0]
ax.plot(winds, drop_before_40, 'o-', color='tab:red', label='sza40 전(버그)', lw=2, ms=7)
ax.plot(winds, drop_after_40,  's-', color='tab:green', label='sza40 후(B)', lw=2, ms=7)
ax.plot(winds, drop_after_80,  '^--', color='tab:blue', label='sza80 후(B, wind-무관=물리)', lw=2, ms=7)
ax.axhline(0, color='k', lw=0.8)
ax.set_xlabel('wind speed [m/s]'); ax.set_ylabel('rrs(0$-$) atm 민감도 [%]\n(Rayleigh vs no-atm)')
ax.set_title('(C) rrs(0$-$) wind 의존성: 버그 signature\n전(wind 비례 하락) → 후(wind 무관)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# ===== Panel D: TOA OCRT-B vs OSOAA scatter (선형) =====
ax = axes[1,1]
for b in band_list:
    m = [i for i in range(len(rows)) if bands[i]==b]
    ax.scatter(os_toa[m], oc_toa[m], s=28, c=cmap[b], label=f'{b} nm', alpha=0.8, edgecolors='none')
lim2 = [0, max(os_toa.max(),oc_toa.max())*1.1]
ax.plot(lim2, lim2, 'k--', lw=1, alpha=0.6)
ax.set_xlim(lim2); ax.set_ylim(lim2)
ax.set_xlabel('OSOAA TOA $\\rho_I$'); ax.set_ylabel('OCRT-B TOA $\\rho_I$')
ax.set_title('(D) TOA $\\rho_I$ OCRT-B vs OSOAA, sza=40 (선형)\nMAPE 0.84%→0.63% (water-leaving Lu 개선)')
ax.legend(fontsize=7, loc='upper left'); ax.grid(True, alpha=0.3)

plt.suptitle('FIX-SKY-EDLU: skylight Ed/Lu sub-cone leak 수정 검증 (OCRT v1.08 vs OSOAA)', fontsize=13, y=1.00)
plt.tight_layout()
plt.savefig('fix_sky_edlu_validation.png', dpi=130, bbox_inches='tight')
print('saved fix_sky_edlu_validation.png')

# RMS bias 요약
mape_rrs = np.mean(np.abs(oc_rrs/os_rrs - 1))*100
mape_toa = np.mean(np.abs(oc_toa/os_toa - 1))*100
print(f'sza40 전체: rrs MAPE={mape_rrs:.2f}%, TOA MAPE={mape_toa:.2f}%')
