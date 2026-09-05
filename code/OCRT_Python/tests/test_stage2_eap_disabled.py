from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pathlib import Path
import numpy as np
import pytest

from ocrt_py.constituent import (
    OCRTConstituentModel, validate_phyto_group_option,
    DEFAULT_PHYTO_GROUP,
)

ROOT = Path(__file__).resolve().parents[1]


def test_default_chl_is_absorption_only_and_detritus_scattering_remains():
    model = OCRTConstituentModel(str(ROOT / 'data'), phyto_group=DEFAULT_PHYTO_GROUP)
    r = model.evaluate(443.0, chl=0.3, tsm=0.0, acdom440=0.0)
    assert r['a_phyto'] > 0.0
    assert r['b_phyto'] == 0.0
    assert r['bb_phyto'] == 0.0
    assert r['b_det'] > 0.0
    assert r['bb_det'] > 0.0
    assert np.isfinite([r['a'], r['b'], r['bb']]).all()
    # PLOPS-derived class-balanced absorption: a*(443)=0.01838492806522 m2/mg.
    assert r['a_phyto'] == pytest.approx(0.005515478419566, abs=2e-15)


def test_batch_chl_is_absorption_only():
    model = OCRTConstituentModel(str(ROOT / 'data'), phyto_group=DEFAULT_PHYTO_GROUP)
    r = model.evaluate_batch(443.0, np.array([0.0, 0.3, 1.0]),
                             np.zeros(3), np.zeros(3))
    assert np.all(np.asarray(r['b_phyto']) == 0.0)
    assert np.all(np.asarray(r['bb_phyto']) == 0.0)
    assert float(np.asarray(r['a_phyto'])[1]) == pytest.approx(0.005515478419566, abs=2e-15)
    assert float(np.asarray(r['b_det'])[1]) > 0.0


def test_explicit_species_rejected_for_positive_chl(monkeypatch):
    monkeypatch.setenv('OCRT_ADVANCED', '1')
    for group in ('micro', 'eap_synechococcus', 'eap_diatoms_pennate'):
        with pytest.raises(ValueError, match='not supported yet'):
            validate_phyto_group_option(group, chl=0.3, explicit=True)


def test_implicit_default_resolves():
    assert validate_phyto_group_option(None, chl=1.0, explicit=False) == 'micro'
