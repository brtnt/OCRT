import os
from ocrt_py.constituent import PhytoAbsorptionLUT


def _lut():
    root = os.path.dirname(os.path.dirname(__file__))
    return PhytoAbsorptionLUT(
        os.path.join(root, 'data', 'water_iop', 'phyto_absorption_default.csv'))


def test_phyto_absorption_source_taper_and_explicit_zero_from_800_to_1100():
    lut = _lut()
    assert lut.wl[0] <= 330.0
    assert lut.wl[-1] >= 1100.0
    assert lut.eval(750.0) > 0.0
    assert lut.eval(775.0) > 0.0
    assert lut.eval(799.0) > 0.0
    assert lut.eval(800.0) == 0.0
    assert lut.eval(850.0) == 0.0
    assert lut.eval(865.0) == 0.0
    assert lut.eval(1100.0) == 0.0
