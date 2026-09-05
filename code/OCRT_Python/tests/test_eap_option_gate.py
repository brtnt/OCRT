from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pytest
from ocrt_py.constituent import (
    phyto_group_is_default, validate_phyto_group_option)


def test_default_alias_and_canonical():
    assert phyto_group_is_default('micro')
    assert phyto_group_is_default('eap_diatoms_centric')
    assert not phyto_group_is_default('pico')


def test_catalog_gate_when_chl_zero(monkeypatch):
    monkeypatch.delenv('OCRT_ADVANCED', raising=False)
    assert validate_phyto_group_option('micro', chl=0, explicit=True) == 'micro'
    with pytest.raises(ValueError):
        validate_phyto_group_option('eap_synechococcus', chl=0, explicit=True)
    monkeypatch.setenv('OCRT_ADVANCED', '1')
    assert validate_phyto_group_option('eap_synechococcus', chl=0, explicit=True) == 'eap_synechococcus'
