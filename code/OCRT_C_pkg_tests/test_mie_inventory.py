"""pytest: canonical mie inventory & audit integrity (release blockers #1,#2)."""
import csv, os
ROOT=os.path.join(os.path.dirname(__file__),'..','01_OCRT_C')
def test_unique_mie_names_198():
    names=[l.strip() for l in open(os.path.join(ROOT,'MIE_CANONICAL_INVENTORY_198.txt')) if l.strip()]
    assert len(names)==198 and len(set(names))==198
def test_full_audit_rows_198():
    rows=list(csv.DictReader(open(os.path.join(ROOT,'MIE_FULL_AUDIT_198.csv'))))
    assert len(rows)==198
    for r in rows:
        assert os.path.exists(os.path.join(ROOT,r['path'])), r['path']
