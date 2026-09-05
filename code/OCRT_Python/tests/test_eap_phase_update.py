from pathlib import Path
import hashlib
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ocrt_py.aerosol import read_mie, water_component_phase_moments_gauss
from ocrt_py.constituent import OCRTConstituentModel, PHYTO_FILE, PHYTO_GROUP_CHOICES


def test_eap_catalog_files_and_legacy_aliases():
    root = Path(__file__).resolve().parents[1]
    data = root / "data" / "water_iop"
    files = sorted((data / "eap").glob("*.mie"))
    assert len(files) == 17
    assert not (data / "pico_Synechococcus_EAP.mie").exists()
    assert not (data / "nano_Haptophytes_EAP.mie").exists()
    assert not (data / "Diatoms_centric_EAP.mie").exists()
    assert PHYTO_FILE["pico"] == PHYTO_FILE["eap_synechococcus"]
    assert PHYTO_FILE["nano"] == PHYTO_FILE["eap_hapto_prymnesiaceae"]
    assert PHYTO_FILE["micro"] == PHYTO_FILE["eap_diatoms_centric"]
    assert len(PHYTO_GROUP_CHOICES) == 20


def test_eap_file_shape_normalization_and_physicality():
    root = Path(__file__).resolve().parents[1]
    for path in sorted((root / "data" / "water_iop" / "eap").glob("*.mie")):
        mie = read_mie(path)
        assert (mie.n_wl, mie.n_phase_wl, mie.n_ang) == (101, 12, 361)
        assert np.isfinite(mie.spectral).all()
        assert np.isfinite(mie.P11).all()
        assert np.isfinite(mie.P12).all()
        assert np.isfinite(mie.P33).all()
        theta = np.deg2rad(mie.angles)
        norm = np.abs(0.5 * np.trapezoid(mie.P11 * np.sin(theta)[:, None], theta, axis=0))
        assert np.max(np.abs(norm - 1.0)) < 2.0e-5
        assert np.max(np.abs(mie.P12) - mie.P11) <= 1.0e-12
        assert np.max(np.abs(mie.P33) - mie.P11) <= 1.0e-12


def test_eap_alias_phase_moments_are_identical():
    root = Path(__file__).resolve().parents[1]
    data = root / "data" / "water_iop"
    for alias, canonical in (("pico", "eap_synechococcus"),
                             ("nano", "eap_hapto_prymnesiaceae"),
                             ("micro", "eap_diatoms_centric")):
        a = read_mie(data / PHYTO_FILE[alias])
        b = read_mie(data / PHYTO_FILE[canonical])
        assert hashlib.sha256(Path(a.source_path).read_bytes()).digest() == \
               hashlib.sha256(Path(b.source_path).read_bytes()).digest()
        am = water_component_phase_moments_gauss(a, 555.0, L_max=24, mie_n_mu=40)
        bm = water_component_phase_moments_gauss(b, 555.0, L_max=24, mie_n_mu=40)
        for x, y in zip(am[:4], bm[:4]):
            assert np.array_equal(x, y)


def test_all_eap_groups_produce_finite_iops():
    for group in PHYTO_GROUP_CHOICES:
        model = OCRTConstituentModel(str(ROOT / "data"), phyto_group=group)
        iop = model.evaluate(555.0, chl=1.0, tsm=0.0, acdom440=0.0)
        values = np.asarray([iop["a"], iop["b"], iop["bb"]], float)
        assert np.isfinite(values).all()
        assert np.all(values > 0.0)
        assert iop["a_phyto"] > 0.0
        assert iop["b_phyto"] == 0.0
        assert iop["bb_phyto"] == 0.0
        assert iop["b_det"] > 0.0 and iop["bb_det"] > 0.0
