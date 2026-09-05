import json
import os
import numpy as np


def _root():
    return os.path.dirname(os.path.dirname(__file__))


def test_water_table_matches_full_reference_range():
    p=os.path.join(_root(),'data','water_iop','water_coef_z09_1nm.txt')
    rows=[]; data=False
    for line in open(p):
        if line.strip().startswith('/end_header'):
            data=True; continue
        if data and line.strip():
            parts=line.split()
            if len(parts)>=3: rows.append(tuple(map(float,parts[:3])))
    a=np.asarray(rows)
    assert a.shape == (2250,3)
    assert a[0,0] == 200.0 and a[-1,0] == 2449.0
    assert np.all(np.diff(a[:,0]) == 1.0)


def test_spectral_manifest_records_complete_330_1100_contract():
    p=os.path.join(_root(),'data','water_iop','SPECTRAL_SUPPORT_330_1100.json')
    m=json.load(open(p))
    assert m['official_elastic_range_nm'] == [330.0,1100.0]
    assert m['runtime_policy']['silent_endpoint_clamp'] is False
    assert m['phytoplankton_absorption']['range_nm'] == [330.0,1100.0]
    assert m['phytoplankton_absorption']['method'].endswith('800-1100 exact zero')
    d=m['chl_linked_detritus_scattering']
    assert d['active_phase_range_nm'] == [350.0,850.0]
    assert d['status'].startswith('validated legacy phase retained')
    assert d['candidate_range_nm'] == [330.0,1100.0]
    assert d['candidate_decision'].startswith('rejected')
    assert m['atmospheric_aerosol']['mie_files'] == 176
    assert m['absorbing_gases']['range_nm'] == [330.0,1100.0]
