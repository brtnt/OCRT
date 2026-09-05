import warnings
from ocrt_py import water


def test_exact_nadir_warns_without_substitution():
    water._VALUE_POL_NADIR_WARNED = False
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        water._warn_value_kernel_exact_nadir([0.0, 30.0])
    assert len(rec) == 1
    msg = str(rec[0].message)
    assert 'VZA=0 deg' in msg
    assert 'VZA=0.001 deg' in msg
    assert 'not modified automatically' in msg


def test_validated_near_nadir_workaround_does_not_warn():
    water._VALUE_POL_NADIR_WARNED = False
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        water._warn_value_kernel_exact_nadir([0.001])
    assert not rec
