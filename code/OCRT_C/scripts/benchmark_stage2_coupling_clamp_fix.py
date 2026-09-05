#!/usr/bin/env python3
from __future__ import annotations
import csv, json, os, statistics, subprocess, tempfile, time
from pathlib import Path

root=Path(__file__).resolve().parents[1]
val=root/'validation/stage2_coupling_clamp_fix_2026-07-25/runtime'
val.mkdir(parents=True,exist_ok=True)
bins={'before':root/'validation/stage2_coupling_clamp_fix_2026-07-25/runtime/bin/ocrt_before_clamp','after':root/'build/ocrt'}
env=dict(os.environ,OCRT_ADVANCED='1',OCRT_DEBUG='1',OMP_NUM_THREADS='1')
base=[
 '--sza','40','--wavelength','443','--surface','ocean','--wind-speed','3','--pressure','1013.25',
 '--n-mu','24','--n-layers','40','--m-max','2','--sos-max-orders','100','--debug-water-max-orders','100',
 '--n-water','1.34','--water-temperature','20','--water-salinity','38.4','--water-model','ocrt',
 '--ocrt-chl','0','--ocrt-tsm','0','--ocrt-adom440','0','--ocrt-adom-slope','0.014','--water-m-max','4',
 '--gas-column-h2o','0','--gas-column-o3','0','--gas-column-no2','0','--gas-column-o2','0',
 '--gas-column-co2','0','--gas-column-ch4','0','--decouple-sunglint',
 '--lut-vza-step','30','--lut-vza-max','60','--lut-raa-step','45'
]
raw=[]
with tempfile.TemporaryDirectory(prefix='ocrt-clamp-bench-') as td:
 td=Path(td)
 def one(n,variant,tag):
  out=td/f'n{n}_{variant}_{tag}.csv'
  cmd=[str(bins[variant]),*base,'--n-mu-water',str(n),'--output-full-grid',str(out)]
  t0=time.perf_counter()
  p=subprocess.run(cmd,cwd=root,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=240)
  elapsed=time.perf_counter()-t0
  if p.returncode or not out.exists() or out.stat().st_size==0: raise RuntimeError((n,variant,p.returncode))
  return elapsed
 for n in (48,64,96):
  # one warm-up per binary
  for v in ('before','after'):
   t=one(n,v,'warmup'); raw.append({'n_mu_water':n,'variant':v,'phase':'warmup','rep':0,'elapsed_s':t})
   print(f'n={n} {v} warmup {t:.6f}s',flush=True)
  for rep,v in enumerate(('before','after','after','before','before','after'),1):
   t=one(n,v,f'rep{rep}'); raw.append({'n_mu_water':n,'variant':v,'phase':'measure','rep':rep,'elapsed_s':t})
   print(f'n={n} {v} rep={rep} {t:.6f}s',flush=True)

raw_path=val/'COUPLING_CLAMP_RUNTIME_RAW_20260725.csv'
with raw_path.open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(raw[0])); w.writeheader(); w.writerows(raw)
summary=[]
for n in (48,64,96):
 vals={v:[r['elapsed_s'] for r in raw if r['n_mu_water']==n and r['variant']==v and r['phase']=='measure'] for v in ('before','after')}
 bm,am=statistics.mean(vals['before']),statistics.mean(vals['after'])
 bmed,amed=statistics.median(vals['before']),statistics.median(vals['after'])
 summary.append({
  'n_mu_water':n,'repetitions_per_variant':3,
  'before_mean_s':bm,'after_mean_s':am,'mean_change_pct':100*(am/bm-1),
  'before_median_s':bmed,'after_median_s':amed,'median_change_pct':100*(amed/bmed-1),
  'before_min_s':min(vals['before']),'after_min_s':min(vals['after']),
 })
sum_path=val/'COUPLING_CLAMP_RUNTIME_SUMMARY_20260725.csv'
with sum_path.open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0])); w.writeheader(); w.writerows(summary)
decision={
 'profile_required':any(max(r['mean_change_pct'],r['median_change_pct'])>=5.0 for r in summary),
 'threshold_pct':5.0,'summary':summary,
}
(val/'COUPLING_CLAMP_RUNTIME_DECISION_20260725.json').write_text(json.dumps(decision,indent=2)+'\n',encoding='utf-8')
print(json.dumps(decision,indent=2))
