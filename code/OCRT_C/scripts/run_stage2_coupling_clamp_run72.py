from __future__ import annotations
import csv, hashlib, json, os, shlex, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

root=Path('/mnt/data/ocrt_stage2_clamp_fix_work/OCRT_v1.2_STAGE2_COUPLING_CLAMP_FIX_20260725')
binary=root/'build/ocrt'
man=root/'validation/stage2_water_raa_output_fix_2026-07-23/RUN_MANIFEST_72CASES.csv'
base=root/'validation/stage2_surface_pole_limit_2026-07-24/run72/outputs'
outroot=root/'validation/stage2_coupling_clamp_fix_2026-07-25/run72'
outdir=outroot/'outputs'; logdir=outroot/'logs'
outdir.mkdir(parents=True,exist_ok=True); logdir.mkdir(parents=True,exist_ok=True)
with man.open(newline='',encoding='utf-8') as f: rows=list(csv.DictReader(f))
env=dict(os.environ, OCRT_ADVANCED='1', OCRT_DEBUG='1', OMP_NUM_THREADS='1')

def sha(p:Path)->str:
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def run_one(row):
 rid=row['run_id']; out=outdir/f'{rid}.csv'; stdout=logdir/f'{rid}.stdout'; stderr=logdir/f'{rid}.stderr'
 ref=base/f'{rid}.csv'
 if out.exists() and out.stat().st_size and ref.exists() and sha(out)==sha(ref):
  with out.open(newline='', encoding='utf-8') as f: data=list(csv.DictReader(f))
  return {
   'run_id':rid,'case_id':row['case_id'],'band_nm':row['band_nm'],'sza_deg':row['sza_deg'],
   'returncode':0,'elapsed_s':0.0,'rows':len(data),'expected_rows':int(row['expected_grid_rows']),
   'water_converged_all':int(all(int(float(r.get('water_converged','0')))==1 for r in data)),
   'max_water_orders':max(int(float(r.get('water_orders','0'))) for r in data),
   'byte_identical_to_pole_baseline':1,'output_sha256':sha(out),'reference_sha256':sha(ref),
  }
 tmpl=row['command_template'].strip()
 if tmpl.startswith('$ '): tmpl=tmpl[2:]
 tmpl=tmpl.replace('<OCRT_BINARY>',shlex.quote(str(binary))).replace('<OUTPUT_CSV>',shlex.quote(str(out)))
 cmd=shlex.split(tmpl)
 t0=time.perf_counter()
 p=subprocess.run(cmd,cwd=root,env=env,text=True,capture_output=True,timeout=600)
 elapsed=time.perf_counter()-t0
 stdout.write_text(p.stdout,encoding='utf-8'); stderr.write_text(p.stderr,encoding='utf-8')
 nrows=-1; conv=False; maxord=None
 if out.exists() and out.stat().st_size:
  with out.open(newline='',encoding='utf-8') as f: data=list(csv.DictReader(f))
  nrows=len(data)
  try:
   conv=all(int(float(r.get('water_converged','0'))) == 1 for r in data)
   maxord=max(int(float(r.get('water_orders','0'))) for r in data)
  except Exception: pass
 identical=(p.returncode==0 and ref.exists() and out.exists() and sha(ref)==sha(out))
 return {
  'run_id':rid,'case_id':row['case_id'],'band_nm':row['band_nm'],'sza_deg':row['sza_deg'],
  'returncode':p.returncode,'elapsed_s':elapsed,'rows':nrows,'expected_rows':int(row['expected_grid_rows']),
  'water_converged_all':int(conv),'max_water_orders':maxord if maxord is not None else '',
  'byte_identical_to_pole_baseline':int(identical),
  'output_sha256':sha(out) if out.exists() else '',
  'reference_sha256':sha(ref) if ref.exists() else '',
 }

results=[]; start=time.perf_counter()
with ThreadPoolExecutor(max_workers=8) as ex:
 futs={ex.submit(run_one,r):r['run_id'] for r in rows}
 for i,fut in enumerate(as_completed(futs),1):
  res=fut.result(); results.append(res)
  print(f"[{i:02d}/72] {res['run_id']} rc={res['returncode']} rows={res['rows']} identical={res['byte_identical_to_pole_baseline']} elapsed={res['elapsed_s']:.2f}s",flush=True)
results.sort(key=lambda x:x['run_id'])
status=outroot/'RUN72_STATUS_COUPLING_CLAMP_20260725.csv'
with status.open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0])); w.writeheader(); w.writerows(results)
summary={
 'runs':len(results),'success':sum(r['returncode']==0 for r in results),
 'expected_row_count_pass':sum(r['rows']==r['expected_rows'] for r in results),
 'water_converged_pass':sum(r['water_converged_all']==1 for r in results),
 'byte_identical_to_pole_baseline':sum(r['byte_identical_to_pole_baseline']==1 for r in results),
 'max_water_orders':max(int(r['max_water_orders']) for r in results if r['max_water_orders']!=''),
 'wall_s':time.perf_counter()-start,
 'binary_sha256':sha(binary),
}
(outroot/'RUN72_SUMMARY_COUPLING_CLAMP_20260725.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))
if not (summary['success']==72 and summary['expected_row_count_pass']==72 and summary['water_converged_pass']==72 and summary['byte_identical_to_pole_baseline']==72):
 raise SystemExit(1)
