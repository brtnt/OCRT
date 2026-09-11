#!/usr/bin/env python3
"""Plan/run the 10 core C examples in manual Appendix C; Python standard library only.

No solver runs without --run. Paths supplied on the CLI are relative to the
invocation directory; paths in recipe arguments are relative to --ocrt-root.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from datetime import datetime, timezone


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_manifest(path, record):
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


def inspect_grid(path, spec):
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        rows = list(reader)
    if len(rows) != spec['rows']:
        raise ValueError(f'{path.name}: expected {spec["rows"]} rows, got {len(rows)}')
    ocean = spec['kind'] == 'ocean_grid'
    expected_columns = 60 if ocean else 13
    if len(fields) != expected_columns:
        raise ValueError(f'{path.name}: expected {expected_columns} columns, got {len(fields)}')
    keys = ['TOA_rho_I', 'TOA_rho_Q', 'TOA_rho_U', 'Rrs_I', 'rrs_I', 'a_total', 'b_total'] if ocean else ['rho_I', 'rho_Q', 'rho_U']
    for i, row in enumerate(rows, 2):
        for key in keys:
            if not math.isfinite(float(row[key])):
                raise ValueError(f'{path.name}:{i}: nonfinite {key}')
        if ocean and int(row['water_converged']) != 1:
            raise ValueError(f'{path.name}:{i}: water_converged != 1')
    angles = [(float(r['vza_deg']), float(r['raa_deg'])) for r in rows]
    expected = [(v, r) for v in (0, 30, 60) for r in (0, 90, 180, 270)]
    if angles != expected:
        raise ValueError(f'{path.name}: unexpected angle coordinates/order')
    result = {'path': str(path), 'rows': len(rows), 'columns': len(fields), 'sha256': digest(path)}
    if ocean:
        result['invalid_upward_transmittance_rows'] = sum(int(r['T_up_rt_valid']) != 1 for r in rows)
    return result


def inspect_single(path, args):
    values = dict(re.findall(r'([A-Za-z][A-Za-z0-9_]*)=([^\s]+)', path.read_text(encoding='utf-8', errors='replace')))
    keys = ['TOA_rho_I', 'TOA_rho_Q', 'TOA_rho_U']
    if '--water-model' in args:
        keys += ['Rrs0plus_I', 'rrs0minus_I', 'a_total', 'b_total']
    for key in keys:
        if key not in values or not math.isfinite(float(values[key])):
            raise ValueError(f'{path.name}: missing/nonfinite {key}')
    if int(values.get('conv', '0')) != 1:
        raise ValueError(f'{path.name}: conv != 1')
    return {'path': str(path), 'sha256': digest(path), 'observables': {k: float(values[k]) for k in keys},
            'conv': int(values['conv']), 'T_up_rt_valid': int(values.get('T_up_rt_valid', '0')),
            'note': ('Ocean single-text conv reports the returned water solve convergence; it is not an independent atmospheric/full-coupling certificate.'
                     if '--water-model' in args else 'Atmospheric single-text conv aggregates Fourier SOS convergence; it is not an accuracy bound.')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--list', action='store_true', help='list recipe IDs/groups; no paths required')
    p.add_argument('--select', default='all', help='comma-separated recipe IDs, or all (default)')
    p.add_argument('--group', help='comma-separated groups; filters --select')
    p.add_argument('--ocrt-root', type=Path, help='C working directory, normally code/OCRT_C')
    p.add_argument('--exe', type=Path, help='built C executable, relative to invocation directory')
    p.add_argument('--outdir', type=Path, help='new output directory, or an existing unexecuted plan directory')
    p.add_argument('--timeout', type=float, default=900, help='seconds per C invocation (default 900)')
    p.add_argument('--threads', type=int, default=1, help='OMP threads per invocation (default 1)')
    p.add_argument('--run', action='store_true', help='execute serially after preflight')
    opts = p.parse_args()
    catalog_path = Path(__file__).with_name('recipe_catalog.json')
    catalog = json.loads(catalog_path.read_text(encoding='utf-8'))
    recipes = catalog['recipes']
    if opts.list:
        for r in recipes:
            print(f'{r["id"]}  {r["group"]:10}  {r["title_en"]}')
        return 0
    if not all((opts.ocrt_root, opts.exe, opts.outdir)):
        p.error('--ocrt-root, --exe and --outdir are required unless --list is used')
    if opts.threads < 1 or not math.isfinite(opts.timeout) or opts.timeout <= 0:
        p.error('--threads must be positive and --timeout must be finite and positive')
    ids = {r['id'] for r in recipes}
    selected = ids if opts.select == 'all' else {s.strip().upper() for s in opts.select.split(',')}
    if selected - ids:
        p.error('unknown recipe IDs: ' + ', '.join(sorted(selected - ids)))
    groups = {r['group'] for r in recipes}
    wanted_groups = groups if not opts.group else {s.strip() for s in opts.group.split(',')}
    if wanted_groups - groups:
        p.error('unknown groups: ' + ', '.join(sorted(wanted_groups - groups)))
    recipes = [r for r in recipes if r['id'] in selected and r['group'] in wanted_groups]
    if not recipes:
        p.error('selection is empty')
    root, exe, out = opts.ocrt_root.resolve(), opts.exe.resolve(), opts.outdir.resolve()
    if not root.is_dir() or not exe.is_file():
        p.error('C working directory or executable does not exist')
    if any(c in str(out) for c in ',\r\n"') and any(r.get('fixtures') for r in recipes):
        p.error('native batch recipes require an output path without commas, quotes or newlines')
    missing = sorted({str(root / name) for r in recipes for name in r['requires'] if not (root / name).is_file()})
    if missing:
        print('Missing required input data:\n' + '\n'.join(missing), file=sys.stderr)
        return 2
    replace_out = lambda value: value.replace('{out}', out.as_posix())
    commands = []
    for r in recipes:
        for c in r['commands']:
            command = dict(c)
            command['recipe_id'] = r['id']
            if r.get('scientific_caution'):
                command['scientific_caution'] = r['scientific_caution']
            command['argv'] = [str(exe)] + [replace_out(a) for a in c['args']]
            command['environment'] = {'OMP_NUM_THREADS': str(opts.threads), **c['env']}
            command['stdout'] = str(out / c['stdout'])
            command['stderr'] = str(out / c['stderr'])
            command['status'] = 'planned'
            commands.append(command)
    planned_results = [Path(c[k]) for c in commands for k in ('stdout', 'stderr')]
    planned_results += [out / s['path'] for c in commands for s in c['outputs']]
    existing = [str(path) for path in planned_results if path.exists()]
    if existing:
        print('Existing results will not be overwritten; use a new --outdir:\n' + '\n'.join(existing), file=sys.stderr)
        return 2
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / 'manifest.json'
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding='utf-8'))
        if previous.get('status') != 'planned' or [c['argv'] for c in previous.get('commands', [])] != [c['argv'] for c in commands]:
            p.error('existing manifest is not the same unexecuted plan; use a new --outdir')
    fixtures = []
    for r in recipes:
        for item in r.get('fixtures', []):
            path, text = out / item['path'], replace_out(item['text'])
            if path.exists() and path.read_text(encoding='utf-8') != text:
                p.error(f'existing input fixture differs: {path}')
            path.write_text(text, encoding='utf-8', newline='\n')
            fixtures.append({'path': str(path), 'sha256': digest(path)})
    inputs = sorted({name for r in recipes for name in r['requires']})
    record = {'status': 'planned', 'created_utc': datetime.now(timezone.utc).isoformat(),
              'source_commit_documented': catalog['source_commit'], 'catalog_sha256': digest(catalog_path),
              'working_directory': str(root), 'executable': str(exe), 'executable_sha256': digest(exe),
              'timeout_seconds': opts.timeout, 'recipe_ids': [r['id'] for r in recipes],
              'input_files': [{'path': name, 'sha256': digest(root / name)} for name in inputs],
              'fixtures': fixtures, 'commands': commands,
              'environment_policy': 'Inherited OCRT_* variables are removed; only per-case settings below are enabled.',
              'removed_inherited_ocrt_variable_names': sorted(k for k in os.environ if k.startswith('OCRT_'))}
    write_manifest(manifest_path, record)
    print(f'{len(recipes)} recipes / {len(commands)} C invocations; manifest: {manifest_path}', flush=True)
    if not opts.run:
        print('Plan only. Reuse this command with --run to execute. No solver has run.')
        return 0
    record['status'] = 'running'
    clean_env = {k: v for k, v in os.environ.items() if not k.startswith('OCRT_')}
    failed = 0
    for c in commands:
        print(f'Running {c["id"]} ...', flush=True)
        started = time.monotonic()
        try:
            with Path(c['stdout']).open('xb') as stdout, Path(c['stderr']).open('xb') as stderr:
                result = subprocess.run(c['argv'], cwd=root, env={**clean_env, **c['environment']},
                                        stdout=stdout, stderr=stderr, timeout=opts.timeout, check=False)
            c['returncode'] = result.returncode
            if result.returncode:
                raise ValueError(f'C executable returned {result.returncode}')
            if c['outputs']:
                c['checks'] = [inspect_grid(out / spec['path'], spec) for spec in c['outputs']]
            else:
                c['checks'] = [inspect_single(Path(c['stdout']), c['args'])]
            c['status'] = 'passed_basic_checks'
        except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as exc:
            c['status'] = 'failed'
            c['error'] = str(exc)
            failed += 1
            print(f'  FAILED: {exc}', file=sys.stderr, flush=True)
        c['elapsed_seconds'] = round(time.monotonic() - started, 3)
        write_manifest(manifest_path, record)
    record['status'] = 'completed' if not failed else 'completed_with_failures'
    record['failed_commands'] = failed
    record['finished_utc'] = datetime.now(timezone.utc).isoformat()
    write_manifest(manifest_path, record)
    print(f'Finished: {len(commands) - failed}/{len(commands)} commands passed basic checks. See manifest/logs.')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
