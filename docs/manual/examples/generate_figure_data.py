#!/usr/bin/env python3
"""Regenerate the manual's measured Rayleigh/glint illustration (Python standard library)."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ocrt-root',type=Path,required=True)
    parser.add_argument('--exe',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    root=args.ocrt_root.resolve(); exe=args.exe.resolve(); out=args.outdir.resolve()
    if not exe.is_file() or not (root/'inputs').is_dir():
        parser.error('Specify a built C executable and its code/OCRT_C working directory.')
    cases=[('glint_w1',1,False),('glint_w5',5,False),('decoupled_w5',5,True)]
    if any((out/(name+'.csv')).exists() for name,_,_ in cases):
        parser.error('Output already exists. Use a new output directory.')
    out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    # Reproduce default physics and quadrature; do not inherit experimental overrides.
    for key in list(env):
        if key.startswith('OCRT_'): del env[key]
    env['OMP_NUM_THREADS']='1'
    records=[]; rows_by_case={}
    for name,wind,decoupled in cases:
        argv=[str(exe),'--surface','black_fresnel_ocean','--wind-speed',str(wind),
              '--sza','40','--wavelength','555','--pressure','1013.25',
              '--lut-vza-step','40','--lut-vza-max','40','--lut-raa-step','5',
              '--output-full-grid',str(out/(name+'.csv'))]
        for gas in ('h2o','o3','no2','o2','co2','ch4'):
            argv += ['--gas-column-'+gas,'0']
        if decoupled: argv += ['--decouple-sunglint']
        with (out/(name+'.stdout.txt')).open('wb') as stdout, (out/(name+'.stderr.txt')).open('wb') as stderr:
            result=subprocess.run(argv,cwd=root,env=env,stdout=stdout,stderr=stderr,timeout=180)
        if result.returncode: raise SystemExit(f'{name} failed: inspect its stderr log.')
        with (out/(name+'.csv')).open(newline='',encoding='utf-8-sig') as stream:
            rows=list(csv.DictReader(stream))
        if len(rows)!=144 or any(not math.isfinite(float(r[k])) for r in rows for k in ('rho_I','rho_Q','rho_U')):
            raise SystemExit(f'{name}: unexpected rows or nonfinite Stokes reflectances.')
        rows_by_case[name]={int(float(r['raa_deg'])):r for r in rows if float(r['vza_deg'])==40}
        documented_argv=[item.replace(str(exe),'OCRT_EXE').replace(str(out),'FIGURE_OUTDIR') for item in argv]
        records.append({'name':name,'argv_template':documented_argv,'returncode':result.returncode,'grid_rows':len(rows)})
    fields=['raa_deg','rho_I_w1','rho_I_w5','rho_I_decoupled','rho_U_decoupled']
    plot_rows=[]
    for raa in range(0,361,5):
        key=raa%360
        plot_rows.append(dict(zip(fields,[raa,
            rows_by_case['glint_w1'][key]['rho_I'],rows_by_case['glint_w5'][key]['rho_I'],
            rows_by_case['decoupled_w5'][key]['rho_I'],rows_by_case['decoupled_w5'][key]['rho_U']])))
    with (out/'rayleigh_azimuth.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields); writer.writeheader(); writer.writerows(plot_rows)
    dec=rows_by_case['decoupled_w5']
    symmetry={key:max(abs(float(dec[a][key])-sign*float(dec[(-a)%360][key])) for a in dec)
              for key,sign in [('rho_I',1),('rho_Q',1),('rho_U',-1)]}
    peaks={name:max(rows,key=lambda a:float(rows[a]['rho_I'])) for name,rows in rows_by_case.items()}
    version=subprocess.run([str(exe),'--version'],cwd=root,env=env,capture_output=True,text=True,check=True).stdout.strip()
    commit=subprocess.run(['git','rev-parse','HEAD'],cwd=root,capture_output=True,text=True,check=True).stdout.strip()
    meta={'created_utc':datetime.now(timezone.utc).isoformat(),'source_commit':commit,
          'exe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'version_string':version,
          'environment':{'OMP_NUM_THREADS':'1','OCRT_overrides':'removed'},'cases':records,
          'path_placeholders':{'OCRT_EXE':'the --exe argument','FIGURE_OUTDIR':'the --outdir argument'},
          'plot_selection':'VZA=40 only; derived 360-degree closure duplicates RAA=0; native grid excludes 360',
          'symmetry_max_absolute_residual':symmetry,'peak_rho_I_raa_deg':peaks,
          'scope':'Illustration of this configuration; not a production-LUT convergence certificate.'}
    (out/'rayleigh_azimuth_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps({'outdir':str(out),'symmetry':symmetry,'peak_angles':peaks},indent=2))


if __name__=='__main__': main()
