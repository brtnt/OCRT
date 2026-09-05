from pathlib import Path
import numpy as np

from ocrt_py.constituent import WaterCoefLUT, pure_water


def _lut():
    root = Path(__file__).resolve().parents[1]
    return WaterCoefLUT(str(root / 'data' / 'water_iop' / 'water_coef_z09_1nm.txt'))


def test_table_shape_range_and_monotonicity():
    lut = _lut()
    assert len(lut.wl) == 2250
    assert lut.wl[0] == 200.0
    assert lut.wl[-1] == 2449.0
    assert np.all(np.diff(lut.wl) == 1.0)
    assert np.all(np.isfinite(lut.aw)) and np.all(lut.aw > 0.0)
    assert np.all(np.isfinite(lut.bw)) and np.all(lut.bw > 0.0)


def test_key_band_values_and_backscatter_rule():
    lut = _lut()
    refs = {
        380.0: (1.155917e-02, 8.383890e-03),
        412.0: (4.553622e-03, 5.909969e-03),
        443.0: (7.067186e-03, 4.331984e-03),
        490.0: (1.502097e-02, 2.824033e-03),
        510.0: (3.250344e-02, 2.385585e-03),
        555.0: (5.963466e-02, 1.672631e-03),
        620.0: (2.755385e-01, 1.053406e-03),
        680.0: (4.653764e-01, 7.176552e-04),
        745.0: (2.835900e+00, 4.917794e-04),
        865.0: (4.605200e+00, 2.656188e-04),
        970.0: (4.799918e+01, 1.658642e-04),
        1100.0: (1.908413e+01, 9.904724e-05),
    }
    for wl, (aw_ref, bw_ref) in refs.items():
        aw, bw, bbw = pure_water(lut, None, wl, 20.0)
        assert aw == aw_ref
        assert bw == bw_ref
        assert bbw == 0.5 * bw_ref
