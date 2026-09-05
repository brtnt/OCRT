#!/usr/bin/env python3
from pathlib import Path
import argparse, hashlib, json, os, subprocess, time

def h(b): return hashlib.sha256(b).hexdigest()

def run(binpath, cwd, args):
    env=os.environ.copy(); env['TERM']='xterm'; env['OMP_NUM_THREADS']='1'
    t=time.perf_counter()
    p=subprocess.run([str(binpath),*args],cwd=cwd,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    return p, time.perf_counter()-t

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--old',required=True,type=Path); ap.add_argument('--new',required=True,type=Path)
    ap.add_argument('--root',required=True,type=Path); ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--start',type=int,default=0); ap.add_argument('--end',type=int,default=18)
    a=ap.parse_args(); a.out.parent.mkdir(parents=True,exist_ok=True)
    base=['--surface','ocean','--water-model','ocrt','--wind-speed','3','--sza','40','--vza','30','--raa','90',
          '--pressure','0','--gas-column-h2o','0','--gas-column-o3','0','--gas-column-no2','0',
          '--gas-column-o2','0','--gas-column-co2','0','--gas-column-ch4','0','--decouple-sunglint']
    cases=[]
    for adom in ('0.01','0.1','1.0'):
        for wl in ('412','443','555','660'):
            cases.append((f'cdom_{adom}_{wl}',base+['--ocrt-chl','0','--ocrt-tsm','0','--ocrt-adom440',adom,'--wavelength',wl]))
    for tsm in ('0.5','2.0'):
        for wl in ('443','555','660'):
            cases.append((f'tsm_{tsm}_{wl}',base+['--ocrt-chl','0','--ocrt-tsm',tsm,'--ocrt-adom440','0','--wavelength',wl]))
    cases=cases[a.start:a.end]
    rows=[]; ok=True
    for name,args in cases:
        po,to=run(a.old,a.root,args); pn,tn=run(a.new,a.root,args)
        same=po.returncode==pn.returncode==0 and po.stdout==pn.stdout and po.stderr==pn.stderr
        rows.append({'case':name,'same':same,'old_rc':po.returncode,'new_rc':pn.returncode,
                     'old_stdout_sha256':h(po.stdout),'new_stdout_sha256':h(pn.stdout),
                     'old_stderr_sha256':h(po.stderr),'new_stderr_sha256':h(pn.stderr),
                     'old_seconds':to,'new_seconds':tn})
        print(name,'PASS' if same else 'FAIL',f'{to:.3f}s {tn:.3f}s',flush=True)
        ok &= same
    import csv
    with a.out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    if not ok: raise SystemExit(1)
    print(f'PASS: selected {len(cases)} non-Chl cases byte-identical')
if __name__=='__main__': main()
