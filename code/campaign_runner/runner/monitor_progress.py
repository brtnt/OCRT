#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json
from collections import Counter
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--output-dir',type=Path,required=True)
a=p.parse_args()
manifest=a.output_dir/'run_manifest.csv'
if not manifest.is_file():
    raise SystemExit(f'no manifest yet: {manifest}')
with manifest.open(newline='',encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
counts=Counter(r['status'] for r in rows)
bytes_total=sum(int(r.get('output_bytes') or 0) for r in rows)
elapsed=sum(float(r.get('elapsed_s') or 0) for r in rows)
print(json.dumps({
    'completed_manifest_rows':len(rows),
    'status_counts':dict(counts),
    'compressed_output_bytes':bytes_total,
    'sum_child_wall_s':elapsed,
    'latest_runs':rows[-10:],
},ensure_ascii=False,indent=2))
