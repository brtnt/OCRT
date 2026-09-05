#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, math
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KEY = ["campaign","chl_mg_m3","tsm_g_m3","adom440_m_inv","band_nm",
       "sza_deg","wind_ms","n_mu_water"]

def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def recompute(df: pd.DataFrame) -> pd.DataFrame:
    d=df[df.ocrt_build != 'v3'].copy()
    rows=[]
    for k,g in d.groupby(KEY, dropna=False, sort=True):
        r=dict(zip(KEY,k)); r['n_directions']=len(g)
        for level in ('Rrs0plus','rrs0minus'):
            for c in 'IQU':
                s=g[f'OSOAA_{level}_{c}'].to_numpy(float)
                o=g[f'OCRT_{level}_{c}'].to_numpy(float)
                n=np.max(np.abs(s))
                e=np.abs(o-s)/n*100 if n>0 else np.zeros_like(s)
                r[f'{level}_{c}_mean_pct']=round(float(e.mean()),4)
                r[f'{level}_{c}_max_pct']=round(float(e.max()),4)
            do=np.hypot(g[f'OCRT_{level}_Q'],g[f'OCRT_{level}_U'])/g[f'OCRT_{level}_I']
            ds=np.hypot(g[f'OSOAA_{level}_Q'],g[f'OSOAA_{level}_U'])/g[f'OSOAA_{level}_I']
            e=np.abs(do-ds).to_numpy(float)*100
            r[f'{level}_DoLP_mean_pp']=round(float(e.mean()),4)
            r[f'{level}_DoLP_max_pp']=round(float(e.max()),4)
        bf=g['black_frac_I'].dropna().to_numpy(float)
        r['black_frac_I_median']=round(float(np.median(bf)),4) if len(bf) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)

def combo_error_mean(g: pd.DataFrame, level: str, c: str) -> float:
    values=[]
    for _,q in g.groupby(KEY, dropna=False):
        s=q[f'OSOAA_{level}_{c}'].to_numpy(float)
        o=q[f'OCRT_{level}_{c}'].to_numpy(float)
        n=np.max(np.abs(s))
        if n>0: values.extend((np.abs(o-s)/n*100).tolist())
    return float(np.mean(values)) if values else math.nan

def scatter(df: pd.DataFrame, out: Path) -> None:
    out.mkdir(parents=True,exist_ok=True)
    chosen={
      'cdom': df[(df.campaign=='cdom') & (df.ocrt_build=='v5') & (df.n_mu_water==64)],
      'tsm': df[(df.campaign=='tsm') & (df.ocrt_build=='v5')],
      'chl': df[(df.campaign=='chl') & (df.ocrt_build=='v11')],
    }
    for camp,g in chosen.items():
        for level,label in [('Rrs0plus','Rrs(0+)'),('rrs0minus','rrs(0−)')]:
            for c in 'IQU':
                x=g[f'OSOAA_{level}_{c}'].to_numpy(float)
                y=g[f'OCRT_{level}_{c}'].to_numpy(float)
                fig,ax=plt.subplots(figsize=(7.8,7.1))
                ax.scatter(x,y,s=9,alpha=0.38)
                vals=np.concatenate([x[np.isfinite(x)],y[np.isfinite(y)]])
                lo,hi=float(vals.min()),float(vals.max())
                pad=max((hi-lo)*0.05,1e-15); lo-=pad; hi+=pad
                if c in 'QU':
                    lim=max(abs(lo),abs(hi)); lo,hi=-lim,lim
                    ax.axhline(0,linewidth=.7); ax.axvline(0,linewidth=.7)
                ax.plot([lo,hi],[lo,hi],'--',linewidth=1.0,label='1:1')
                ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_aspect('equal',adjustable='box')
                ax.set_xlabel(f'OSOAA {label} {c}')
                ax.set_ylabel(f'OCRT {label} {c}')
                ax.set_title(f'{camp.upper()} full-angle validation — {label} {c}')
                ax.grid(True,alpha=.25); ax.legend(loc='best')
                m=combo_error_mean(g,level,c)
                ax.text(.03,.97,f'Combination-normalized mean |Δ{c}| = {m:.3f}%',
                        transform=ax.transAxes,va='top',ha='left',fontsize=9,
                        bbox=dict(boxstyle='round',alpha=.75))
                fig.tight_layout()
                fig.savefig(out/f'{camp}_{level}_{c}_scatter.png',dpi=180)
                plt.close(fig)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('dir',type=Path); a=ap.parse_args()
    valp=a.dir/'OCRT_OSOAA_underwater_validation_2026-07-26.csv'
    sump=a.dir/'OCRT_OSOAA_underwater_summary_2026-07-26.csv'
    df=pd.read_csv(valp); supplied=pd.read_csv(sump)
    rec=recompute(df)
    rec.to_csv(a.dir/'OCRT_OSOAA_underwater_summary_RECOMPUTED_2026-07-26.csv',index=False)
    left=supplied.sort_values(KEY).reset_index(drop=True)
    right=rec.sort_values(KEY).reset_index(drop=True)
    if left.shape != right.shape or list(left.columns)!=list(right.columns):
        raise SystemExit(f'summary schema mismatch {left.shape} {right.shape}')
    maxdiff=0.0
    mismatches=[]
    for col in left.columns:
        if col in KEY or col=='n_directions':
            if not left[col].fillna('').astype(str).equals(right[col].fillna('').astype(str)):
                mismatches.append(col)
        else:
            x=pd.to_numeric(left[col],errors='coerce').to_numpy(float)
            y=pd.to_numeric(right[col],errors='coerce').to_numpy(float)
            d=np.nanmax(np.abs(x-y)) if np.any(np.isfinite(x-y)) else 0.0
            maxdiff=max(maxdiff,float(d))
            if d>5e-5: mismatches.append(col)
    scatter(df,a.dir/'figures_recomputed')
    base=df[(df.ocrt_build!='v3')]
    claims=[]
    for camp in ('cdom','tsm','chl'):
        q=supplied[supplied.campaign==camp]
        claims.append({'campaign':camp,'n_combinations':len(q),
                       'n_direction_rows':int(len(base[base.campaign==camp])),
                       'rrs_I_mean_min_pct':float(q.rrs0minus_I_mean_pct.min()),
                       'rrs_I_mean_max_pct':float(q.rrs0minus_I_mean_pct.max()),
                       'rrs_Q_mean_min_pct':float(q.rrs0minus_Q_mean_pct.min()),
                       'rrs_Q_mean_max_pct':float(q.rrs0minus_Q_mean_pct.max()),
                       'rrs_U_mean_min_pct':float(q.rrs0minus_U_mean_pct.min()),
                       'rrs_U_mean_max_pct':float(q.rrs0minus_U_mean_pct.max())})
    pd.DataFrame(claims).to_csv(a.dir/'UNDERWATER_CAMPAIGN_CLAIM_RANGES.csv',index=False)
    audit={
      'validation_shape':list(df.shape),'summary_shape':list(supplied.shape),
      'recomputed_summary_shape':list(rec.shape),'summary_max_abs_diff':maxdiff,
      'summary_mismatch_columns':mismatches,'validation_sha256':sha256(valp),
      'summary_sha256':sha256(sump),'campaign_counts':df.campaign.value_counts().to_dict(),
      'build_counts':df.ocrt_build.value_counts().to_dict(),
      'pass': not mismatches and maxdiff<=5e-5 and tuple(df.shape)==(11154,37) and tuple(supplied.shape)==(50,26)
    }
    import json
    (a.dir/'UNDERWATER_DATA_INTEGRITY_AUDIT.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(audit,indent=2,ensure_ascii=False))
    if not audit['pass']: raise SystemExit(1)
if __name__=='__main__': main()
