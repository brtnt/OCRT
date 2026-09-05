#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, hashlib, json, tarfile
from pathlib import Path

TARGET_BASENAMES={
    'runs_summary.csv','aerlib.csv','cases_thin.csv.gz','atm_reference.csv.gz',
    'run_commands.txt','run_manifest.csv','water_states.csv','water_state_table.csv'
}

def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

ap=argparse.ArgumentParser(description='Extract provenance tables from the historical OCRT polarization package without inventing missing state values.')
ap.add_argument('--archive',type=Path,required=True)
ap.add_argument('--output-dir',type=Path,default=Path('historical_import'))
args=ap.parse_args()
archive=args.archive.resolve(); out=args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
found=[]
with tarfile.open(archive,'r:*') as tf:
    for member in tf.getmembers():
        if not member.isfile(): continue
        base=Path(member.name).name
        if base in TARGET_BASENAMES or base.startswith('run_commands') or ('water' in base.lower() and base.lower().endswith('.csv')):
            target=out/member.name
            target.parent.mkdir(parents=True,exist_ok=True)
            src=tf.extractfile(member)
            if src is None: continue
            target.write_bytes(src.read())
            found.append({'member':member.name,'output':str(target),'size_bytes':target.stat().st_size,'sha256':sha256(target)})
summary={'archive':str(archive),'archive_sha256':sha256(archive),'found':found}
# Report actual columns/counts when runs_summary exists; do not infer commands.
for item in found:
    if Path(item['member']).name=='runs_summary.csv':
        path=Path(item['output'])
        with path.open(newline='',encoding='utf-8-sig') as f:
            reader=csv.reader(f); header=next(reader); count=sum(1 for _ in reader)
        summary['runs_summary']={'path':str(path),'columns':header,'rows':count}
(out/'HISTORICAL_PACKAGE_INVENTORY.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
if not found:
    raise SystemExit('no known provenance tables found; exact historical matrix remains unavailable')
