#!/usr/bin/env python3
"""Plan and optionally run small OCRT C LUT campaigns (standard library only)."""
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess


def numbers(text):
    result = [float(x) for x in text.split(',')]
    if not result or not all(math.isfinite(x) for x in result):
        raise argparse.ArgumentTypeError('Use a comma-separated list of finite numbers')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['rayleigh', 'simulation'])
    p.add_argument('--ocrt-root', type=Path, required=True, help='Directory containing src/ and inputs/')
    p.add_argument('--exe', type=Path, required=True, help='Executable; relative paths resolve from caller cwd')
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--wavelengths', type=numbers, default=[443, 555])
    p.add_argument('--szas', type=numbers, default=[30])
    p.add_argument('--pressures', type=numbers, default=[1013.25])
    p.add_argument('--winds', type=numbers, default=[5])
    p.add_argument('--aods', type=numbers, default=[0], help='Simulation AOD at 555 nm')
    p.add_argument('--adom440', type=numbers, default=[0.1], help='Simulation aDOM(440), m^-1; Chl and TSM fixed to zero')
    p.add_argument('--mie', help='Aerosol Mie path, relative to ocrt-root unless absolute')
    p.add_argument('--vza-step', type=float, default=30)
    p.add_argument('--vza-max', type=float, default=60)
    p.add_argument('--raa-step', type=float, default=90)
    p.add_argument('--pssa', action='store_true', help='IPSS for Rayleigh mode only')
    p.add_argument('--include-direct-glint', action='store_true', help='Rayleigh only; default removes direct glint')
    p.add_argument('--gas-off', action='store_true', help='Simulation only; Rayleigh always sets all six gas columns to zero')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--timeout', type=float, default=900, help='Per case seconds')
    p.add_argument('--run', action='store_true', help='Without this switch, only write the plan')
    a = p.parse_args()
    a.ocrt_root = a.ocrt_root.resolve(); a.exe = a.exe.resolve(); a.outdir = a.outdir.resolve()
    if not a.exe.is_file() or not (a.ocrt_root/'src').is_dir() or not (a.ocrt_root/'inputs').is_dir():
        p.error('Check --exe and --ocrt-root; root must be code/OCRT_C')
    if a.threads < 1 or a.timeout <= 0:
        p.error('threads and timeout must be positive')
    if not (0 < a.vza_step <= a.vza_max < 90 and 0 < a.raa_step <= 180):
        p.error('Use 0 < vza-step <= vza-max < 90 and 0 < raa-step <= 180')
    if not math.isclose(360/a.raa_step, round(360/a.raa_step), abs_tol=1e-9):
        p.error('Choose raa-step that divides 360 exactly')
    if not math.isclose(a.vza_max/a.vza_step, round(a.vza_max/a.vza_step), abs_tol=1e-9):
        p.error('Choose vza-step that divides vza-max exactly')
    if any(not 330 <= x <= 1100 for x in a.wavelengths) or any(not 0 <= x < 90 for x in a.szas):
        p.error('Wavelength range is 330..1100 nm; SZA range for this example is 0..<90 degrees')
    if any(x < 0 for seq in [a.pressures,a.winds,a.aods,a.adom440] for x in seq):
        p.error('Pressure, wind, AOD and aDOM must be nonnegative')
    if a.mode == 'rayleigh' and any(x != 0 for x in a.aods):
        p.error('Rayleigh mode requires AOD=0')
    if a.mode == 'simulation' and (a.pssa or a.include_direct_glint):
        p.error('These switches are Rayleigh-only; ocean grid TOA separates direct glint and does not support IPSS')
    if a.mode == 'simulation' and round(360/a.raa_step)>512:
        p.error('Current ocean output arrays support at most 512 RAA samples; use a larger raa-step')
    if any(x > 0 for x in a.aods):
        if not a.mie:
            p.error('Positive AOD requires --mie')
        if not (a.ocrt_root/Path(a.mie)).is_file():
            p.error('Mie file does not exist')
    axes = [a.wavelengths,a.szas,a.pressures,a.winds]
    axes += [a.aods,a.adom440] if a.mode == 'simulation' else [[0],[0]]
    jobs=[]
    for i,(wl,sza,pressure,wind,aod,adom) in enumerate(itertools.product(*axes),1):
        case_id=f'case_{i:06d}'
        result=a.outdir/f'{case_id}.csv'
        args=['--sza',str(sza),'--wavelength',str(wl),'--pressure',str(pressure),
              '--wind-speed',str(wind),'--lut-vza-step',str(a.vza_step),
              '--lut-vza-max',str(a.vza_max),'--lut-raa-step',str(a.raa_step),
              '--output-full-grid',str(result)]
        if a.mode == 'rayleigh':
            args += ['--surface','black_fresnel_ocean']
            if not a.include_direct_glint: args += ['--decouple-sunglint']
            if a.pssa: args += ['--pssa']
        else:
            args += ['--surface','ocean','--water-model','ocrt','--ocrt-chl','0',
                     '--ocrt-tsm','0','--ocrt-adom440',str(adom)]
            if aod > 0: args += ['--aod-555',str(aod),'--mie',a.mie]
        if a.mode == 'rayleigh' or a.gas_off:
            for gas in ['h2o','o3','no2','o2','co2','ch4']:
                args += [f'--gas-column-{gas}','0']
        jobs.append(dict(case_id=case_id,wavelength_nm=wl,sza_deg=sza,pressure_hpa=pressure,
                         wind_m_s=wind,aod_555=aod,adom440_m_1=adom,output=str(result),argv=[str(a.exe),*args]))
    if a.run and any(Path(j['output']).exists() for j in jobs):
        p.error('A result CSV already exists; use a new outdir to prevent accidental overwrite')
    a.outdir.mkdir(parents=True,exist_ok=True)
    env={**os.environ,'OMP_NUM_THREADS':str(a.threads)}
    # Explicitly record inherited OCRT controls rather than silently hiding them.
    version=subprocess.run([str(a.exe),'--version'],cwd=a.ocrt_root,env=env,
                           capture_output=True,text=True,timeout=30,check=True).stdout.strip()
    meta=dict(mode=a.mode,version=version,exe_sha256=hashlib.sha256(a.exe.read_bytes()).hexdigest(),
              cwd=str(a.ocrt_root),OMP_NUM_THREADS=a.threads,
              OCRT_environment={k:v for k,v in env.items() if k.startswith('OCRT_')},
              vza_step=a.vza_step,vza_max=a.vza_max,raa_step=a.raa_step,
              gas_mode='off' if a.mode=='rayleigh' or a.gas_off else 'default_usstd76',
              direct_glint=('separated_from_ocean_TOA_not_exported_in_grid' if a.mode=='simulation'
                            else 'included' if a.include_direct_glint else 'decoupled'),
              ipss=a.pssa,jobs=jobs)
    (a.outdir/'manifest.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
    with (a.outdir/'jobs.csv').open('w',newline='',encoding='utf-8') as fp:
        names=[k for k in jobs[0] if k!='argv']
        w=csv.DictWriter(fp,fieldnames=names,extrasaction='ignore'); w.writeheader(); w.writerows(jobs)
    print(f'{len(jobs)} cases; plan: {a.outdir / "manifest.json"}',flush=True)
    if not a.run: return
    for j in jobs:
        with (a.outdir/f'{j["case_id"]}.stdout.txt').open('w',encoding='utf-8') as so, \
             (a.outdir/f'{j["case_id"]}.stderr.log').open('w',encoding='utf-8') as se:
            r=subprocess.run(j['argv'],cwd=a.ocrt_root,env=env,stdout=so,stderr=se,timeout=a.timeout)
        if r.returncode:
            raise SystemExit(f'{j["case_id"]} failed ({r.returncode}); inspect its stderr log')
        with Path(j['output']).open(newline='',encoding='utf-8') as fp:
            rows=list(csv.DictReader(fp))
        # Match the current C mode-specific rounding, including fractional steps.
        n_vza=math.floor(a.vza_max/a.vza_step+(1e-9 if a.mode=='simulation' else 0))+1
        n_raa=math.floor(360/a.raa_step+(0.5 if a.mode=='simulation' else 0))
        expected=n_vza*n_raa
        required=['sza_deg','wavelength_nm','vza_deg','raa_deg']
        required += ['rho_I','rho_Q','rho_U'] if a.mode=='rayleigh' else ['TOA_rho_I','TOA_rho_Q','TOA_rho_U','Rrs_I','Rrs_Q','Rrs_U','rrs_I','rrs_Q','rrs_U']
        if len(rows)!=expected or any(not math.isfinite(float(row[k])) for row in rows for k in required):
            raise SystemExit(f'{j["case_id"]}: check row count and nonfinite output values')
        if a.mode=='simulation' and any(int(row['water_converged'])!=1 for row in rows):
            raise SystemExit(f'{j["case_id"]}: water solver did not converge')
        print(f'{j["case_id"]}: {len(rows)} rows written',flush=True)


if __name__ == '__main__':
    main()
