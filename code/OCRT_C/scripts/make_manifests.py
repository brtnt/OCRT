#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate the two file manifests of code/OCRT_C.

    MANIFEST_SHA256.csv  relative_path,size_bytes,sha256 — the source package (no Mie tables
                         except the test fixture under tests/fixtures/).
    SHA256SUMS.txt       `sha256sum` format — the source package, MANIFEST_SHA256.csv and the
                         181 Mie tables installed under inputs/ (GitHub Release data-v1).

Conventions (2026-09-13):
  * files are listed relative to code/OCRT_C and sorted in byte order of the path;
  * excluded: build/, _RELEASE/, validation/ipss_2026-09-05/ (cross-check inputs of third-party
    codes), __pycache__/, any path component starting with '.', *.o / *.pyc / *.exe, and the two
    manifests themselves;
  * a Mie table that is not installed locally is listed with the hash recorded in
    ../../scripts/MIE_SHA256SUMS_data-v1.txt, so the output does not depend on whether the
    tables are installed; an installed table whose hash differs from that record is reported
    and the run stops (the table is corrupt or not the data-v1 file).

Usage:  cd code/OCRT_C && python3 scripts/make_manifests.py
Check:  sh verify_package.sh        (sha256sum -c --ignore-missing SHA256SUMS.txt)

한국어: code/OCRT_C 의 매니페스트 두 종을 재생성한다. MANIFEST_SHA256.csv 는 소스 패키지(Mie 표 제외,
tests/fixtures/ 의 시험용 표만 포함), SHA256SUMS.txt 는 소스 패키지 + MANIFEST_SHA256.csv + inputs/ 의
Mie 표 181개다. 제외 규약과 정렬은 위와 같고, 로컬에 없는 Mie 표는 scripts/MIE_SHA256SUMS_data-v1.txt 의
해시로 기재한다. 설치된 표의 해시가 그 기록과 다르면 중단한다.
"""
import csv
import hashlib
import io
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, '..'))
DATA_MANIFEST = os.path.normpath(os.path.join(PKG, '..', '..', 'scripts', 'MIE_SHA256SUMS_data-v1.txt'))
MANIFEST_CSV = 'MANIFEST_SHA256.csv'
SUMS_TXT = 'SHA256SUMS.txt'
EXCLUDED_DIRS = ('build', '_RELEASE', os.path.join('validation', 'ipss_2026-09-05'))
EXCLUDED_SUFFIXES = ('.o', '.pyc', '.exe')


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def is_excluded(rel):
    parts = rel.split('/')
    if any(p.startswith('.') or p == '__pycache__' for p in parts):
        return True
    for d in EXCLUDED_DIRS:
        d = d.replace(os.sep, '/')
        if rel == d or rel.startswith(d + '/'):
            return True
    if rel in (MANIFEST_CSV, SUMS_TXT):
        return True
    return rel.endswith(EXCLUDED_SUFFIXES)


def walk_files():
    out = []
    for root, dirs, files in os.walk(PKG):
        rel_root = os.path.relpath(root, PKG).replace(os.sep, '/')
        rel_root = '' if rel_root == '.' else rel_root + '/'
        # prune excluded directories early
        dirs[:] = sorted(d for d in dirs if not is_excluded(rel_root + d))
        for name in files:
            rel = rel_root + name
            if not is_excluded(rel):
                out.append(rel)
    return sorted(out, key=lambda s: s.encode('utf-8'))


def read_data_manifest():
    ref = {}
    if not os.path.exists(DATA_MANIFEST):
        sys.stderr.write('warning: %s not found; Mie tables not installed locally will be missing from %s\n'
                         % (DATA_MANIFEST, SUMS_TXT))
        return ref
    with io.open(DATA_MANIFEST, encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            digest, _, rel = line.partition('  ')
            rel = rel.strip()
            if rel.startswith('./'):
                rel = rel[2:]
            ref['inputs/' + rel] = digest.strip()
    return ref


def main():
    os.chdir(PKG)
    files = walk_files()
    ref_mie = read_data_manifest()

    source_rows = []   # (rel, size, sha) — everything except Mie tables under inputs/
    local_mie = {}
    for rel in files:
        if rel.endswith('.mie') and rel.startswith('inputs/'):
            local_mie[rel] = sha256_of(rel)
            continue
        source_rows.append((rel, os.path.getsize(rel), sha256_of(rel)))

    # installed tables must match the data-v1 record
    bad = [(rel, h, ref_mie[rel]) for rel, h in local_mie.items() if rel in ref_mie and ref_mie[rel] != h]
    if bad:
        for rel, h, r in bad:
            sys.stderr.write('MISMATCH %s\n  local   %s\n  data-v1 %s\n' % (rel, h, r))
        sys.stderr.write('installed Mie table(s) differ from scripts/MIE_SHA256SUMS_data-v1.txt; nothing written\n')
        return 1
    extra = sorted(set(local_mie) - set(ref_mie))
    if extra:
        sys.stderr.write('note: %d local .mie file(s) not in the data-v1 record are listed with their local hash: %s\n'
                         % (len(extra), ', '.join(extra[:5])))

    # 1) MANIFEST_SHA256.csv
    with io.open(MANIFEST_CSV, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['relative_path', 'size_bytes', 'sha256'])
        for rel, size, sha in source_rows:
            w.writerow([rel, size, sha])

    # 2) SHA256SUMS.txt = source rows + MANIFEST_SHA256.csv + Mie tables (installed or recorded)
    entries = {rel: sha for rel, _, sha in source_rows}
    entries[MANIFEST_CSV] = sha256_of(MANIFEST_CSV)
    mie_entries = dict(ref_mie)
    mie_entries.update(local_mie)
    entries.update(mie_entries)
    with io.open(SUMS_TXT, 'w', encoding='utf-8', newline='') as f:
        for rel in sorted(entries, key=lambda s: s.encode('utf-8')):
            f.write('%s  %s\n' % (entries[rel], rel))

    print('%s: %d files' % (MANIFEST_CSV, len(source_rows)))
    print('%s: %d entries (%d Mie tables, %d of them installed locally)'
          % (SUMS_TXT, len(entries), len(mie_entries), len(local_mie)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
