from pathlib import Path
import numpy as np
import pytest

from ocrt_py.aerosol import read_mie
from ocrt_py.spectral_contract import (
    SPECTRAL_MIN_NM, SPECTRAL_MAX_NM, require_wavelength,
    require_table_coverage, require_query_in_table)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'


def _active_mie_paths():
    yield from sorted(DATA.glob('*.mie'))
    yield from sorted((DATA/'aerosol_ahmad2010_paper_mie').glob('*.mie'))
    yield from sorted((DATA/'aerosol_ahmad2010_accurt_mie').glob('*.mie'))
    yield from sorted((DATA/'tsm_ahn').glob('*.mie'))


def test_official_runtime_wavelength_contract_and_out_of_range_failure():
    assert require_wavelength(330.0) == 330.0
    assert require_wavelength(1100.0) == 1100.0
    with pytest.raises(ValueError):
        require_wavelength(329.999)
    with pytest.raises(ValueError):
        require_wavelength(1100.001)


def test_all_official_full_range_mie_files_explicitly_cover_330_1100():
    paths = list(_active_mie_paths())
    assert len(paths) == 180
    for path in paths:
        mie = read_mie(path)
        require_table_coverage(mie.wavelengths * 1000.0,
                               context=f'{path} bulk')
        require_table_coverage(mie.phase_wavelengths * 1000.0,
                               context=f'{path} phase')
        for q in (330.0, 1100.0):
            require_query_in_table(q, mie.wavelengths * 1000.0,
                                   context=f'{path} bulk')
            require_query_in_table(q, mie.phase_wavelengths * 1000.0,
                                   context=f'{path} phase')


def test_six_gas_tables_have_40_levels_and_771_wavelengths():
    for gas in ('h2o','o2','co2','ch4','o3','no2'):
        lines = [x for x in (DATA/'xsec'/f'xsec_{gas}.dat').read_text().splitlines()
                 if x.strip() and not x.lstrip().startswith('#')]
        a = np.asarray([[float(v) for v in x.split()] for x in lines])
        assert a.shape == (41, 771)
        assert a[0,0] == SPECTRAL_MIN_NM
        assert a[0,-1] == SPECTRAL_MAX_NM
        assert np.all(np.isfinite(a))
        assert np.all(a[1:] >= 0.0)


def test_eap_catalog_is_explicitly_excluded_from_official_range():
    eap = sorted((DATA/'water_iop'/'eap').glob('*.mie'))
    assert len(eap) == 17
    for path in eap:
        mie = read_mie(path)
        assert mie.wavelengths[0] * 1000.0 == 350.0
        assert mie.wavelengths[-1] * 1000.0 == 850.0
        with pytest.raises(ValueError):
            require_query_in_table(330.0, mie.phase_wavelengths * 1000.0,
                                   context=str(path))


def test_active_detritus_keeps_validated_native_range_and_fails_outside():
    path = DATA/'water_iop'/'Detritus_Stramski2001.mie'
    mie = read_mie(path)
    phase_nm = mie.phase_wavelengths * 1000.0
    assert phase_nm[0] == 350.0
    assert phase_nm[-1] == 850.0
    for q in (350.0, 443.0, 850.0):
        require_query_in_table(q, phase_nm, context=str(path))
    for q in (330.0, 1100.0):
        with pytest.raises(ValueError):
            require_query_in_table(q, phase_nm, context=str(path))


def test_rejected_detritus_candidate_is_archived_not_active():
    candidate = (DATA/'water_iop'/'candidates'/
                 'Detritus_Stramski2001_330_1100_CANDIDATE_REJECTED_L200.mie')
    assert candidate.is_file()
    mie = read_mie(candidate)
    require_table_coverage(mie.wavelengths * 1000.0, context='candidate bulk')
    require_table_coverage(mie.phase_wavelengths * 1000.0, context='candidate phase')
