#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, re, shlex
from pathlib import Path

ap=argparse.ArgumentParser(description='Convert exact archived OCRT commands into the auditable campaign matrix format.')
ap.add_argument('--commands',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--default-expected-rows',type=int,default=1296)
args=ap.parse_args()
lines=[]
for raw in args.commands.read_text(encoding='utf-8-sig',errors='replace').splitlines():
    s=raw.strip()
    if s and not s.startswith('#'): lines.append(s)
if not lines: raise SystemExit('command file contains no commands')
fields=['run_id','stage','uses_water','expected_rows','output_relpath','command_args_json','required_columns_json','expected_vza_json','expected_raa_json','notes']
rows=[]; audit=[]
for idx,line in enumerate(lines,1):
    # posix=False is friendlier to archived Windows commands; surrounding quotes are stripped below.
    tokens=shlex.split(line,posix=False)
    tokens=[t.strip('"') for t in tokens]
    if not tokens: continue
    executable=tokens[0]; cmd=tokens[1:]
    output=None
    for option in ['--output-full-grid']:
        if option in cmd:
            i=cmd.index(option)
            if i+1>=len(cmd): raise SystemExit(f'line {idx}: missing value after {option}')
            output=cmd[i+1]
            del cmd[i:i+2]
    if '--pssa' not in cmd:
        raise SystemExit(f'line {idx}: command lacks required --pssa; stopped without modification')
    for forbidden in ['--ocrt-mie-truncation','--iop-mie-truncation']:
        if forbidden in cmd:
            raise SystemExit(f'line {idx}: canonical matrix contains forbidden {forbidden}')
    uses_water='--water-model' in cmd
    if uses_water and '--n-mu-water' not in cmd:
        raise SystemExit(f'line {idx}: water command lacks explicit --n-mu-water')
    if output:
        base=Path(output).name
        rel=f'fullgrid/{Path(base).stem}.csv.gz'
        run_id=re.sub(r'[^A-Za-z0-9_.-]+','_',Path(base).stem)
    else:
        digest=hashlib.sha256(line.encode()).hexdigest()[:16]
        run_id=f'run_{idx:05d}_{digest}'; rel=f'fullgrid/{run_id}.csv.gz'
    rows.append({
        'run_id':run_id,'stage':'historical_import','uses_water':str(uses_water).lower(),
        'expected_rows':str(args.default_expected_rows),'output_relpath':rel,
        'command_args_json':json.dumps(cmd,separators=(',',':')),
        'required_columns_json':json.dumps([],separators=(',',':')),
        'expected_vza_json':json.dumps([],separators=(',',':')),
        'expected_raa_json':json.dumps([],separators=(',',':')),
        'notes':f'Imported verbatim from {args.commands.name}; original executable={executable}; original output={output or "none"}'
    })
    audit.append({'line':idx,'run_id':run_id,'uses_water':uses_water,'has_pssa':True,'has_n_mu_water':('--n-mu-water' in cmd),'original':line})
args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
with args.output.with_suffix('.audit.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(audit[0])); w.writeheader(); w.writerows(audit)
print(json.dumps({'commands':len(lines),'matrix_rows':len(rows),'output':str(args.output)},indent=2))
