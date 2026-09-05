import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ocrt_py.coupling import apply_diffuse_top_basis_sign_inplace


def test_dtpsign_flips_only_odd_modes_and_is_involution():
    base = np.arange(5 * 3, dtype=float).reshape(5, 3) + 0.25
    a = base.copy()
    b = (10.0 + base).copy()
    c = (-3.0 - base).copy()
    apply_diffuse_top_basis_sign_inplace(a, b, c, mode_axis=0)
    np.testing.assert_array_equal(a[0::2], base[0::2])
    np.testing.assert_array_equal(a[1::2], -base[1::2])
    np.testing.assert_array_equal(b[0::2], (10.0 + base)[0::2])
    np.testing.assert_array_equal(b[1::2], -(10.0 + base)[1::2])
    np.testing.assert_array_equal(c[0::2], (-3.0 - base)[0::2])
    np.testing.assert_array_equal(c[1::2], -(-3.0 - base)[1::2])
    apply_diffuse_top_basis_sign_inplace(a, b, c, mode_axis=0)
    np.testing.assert_array_equal(a, base)
    np.testing.assert_array_equal(b, 10.0 + base)
    np.testing.assert_array_equal(c, -3.0 - base)


def test_dtpsign_batch_mode_axis():
    base = np.arange(2 * 4 * 3, dtype=float).reshape(2, 4, 3)
    a = base.copy()
    apply_diffuse_top_basis_sign_inplace(a, mode_axis=1)
    np.testing.assert_array_equal(a[:, 0::2, :], base[:, 0::2, :])
    np.testing.assert_array_equal(a[:, 1::2, :], -base[:, 1::2, :])


def test_single_tsm_component_phase_diagnostic(monkeypatch, capsys):
    import ocrt_solve

    monkeypatch.setenv('OCRT_DEBUG', '1')
    monkeypatch.setenv('OCRT_DUMP_IOP', '1')
    cm = ocrt_solve._constituent_model('micro', 'red_clay')
    result, constituent = ocrt_solve._constituent_greek(
        cm, 555.0, chl=0.0, tsm=5.0, acdom440=0.0, L=32, nmg=64)
    err = capsys.readouterr().err
    assert result['b_min'] > 0.0
    assert constituent['b_particle'] > 0.0
    assert 'OCRT_PHASE_COMPONENT name=tsm ' in err
    assert 'weight_b=1' in err
    assert 'OCRT_PHASE_MIX components=1 ' in err
