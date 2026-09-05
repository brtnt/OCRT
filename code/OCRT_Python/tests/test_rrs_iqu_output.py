import os
import numpy as np


def test_batch_rrs_iqu_reuses_existing_vector_solution(monkeypatch):
    from ocrt_py import aerosol as AER
    from ocrt_py import batch_driver as BD
    from ocrt_py import atmos_batch as AB
    from ocrt_py.constituent import OCRTConstituentModel

    root = os.path.dirname(os.path.dirname(__file__))
    data = os.path.join(root, 'data')
    mie_cache = {'C50': AER.read_mie(os.path.join(data, 'C50.mie'))}
    cm = OCRTConstituentModel(data, 'micro', 'red_clay')
    cases = [dict(sza=30.0, vza=30.0, raa=90.0, wind=3.0,
                  aer='C50', wl=490.0, aod=0.2,
                  chl=1.2, tsm=3.5, ad=0.04)]

    calls = {'atm': 0, 'water': 0}
    captured = {}
    atm0 = AB.solve_atm_batch
    wat0 = AB.solve_water_batch_r1

    def atm_count(*a, **k):
        calls['atm'] += 1
        return atm0(*a, **k)

    def water_count(*a, **k):
        calls['water'] += 1
        out = wat0(*a, **k)
        captured.update(out)
        return out

    monkeypatch.setattr(AB, 'solve_atm_batch', atm_count)
    monkeypatch.setattr(AB, 'solve_water_batch_r1', water_count)
    _, legacy_Rrs_I, tr, _ = BD.solve_r1_grid(
        cases, mie_cache, cm, L_max=16, Lmix=24,
        pressure_hpa=1013.25, n_mu_water=4, nt_atm=12,
        fourier_m_max=3, max_it_atm=6, tol_atm=1e-4,
        n_layers_water_min=12, max_it_water=8, tol_water=1e-4,
        n_phi_quad=64, n_phi_taw=32, absorption=None)

    # No extra RT pass was added: pass-1 atmosphere, one water solve,
    # and pass-2 atmosphere are exactly the pre-existing call graph.
    assert calls == {'atm': 2, 'water': 1}
    assert tr['Rrs_I'][0] == legacy_Rrs_I[0]
    assert np.isfinite([tr[k][0] for k in
                        ('Rrs_I', 'Rrs_Q', 'Rrs_U',
                         'rrs_I', 'rrs_Q', 'rrs_U')]).all()

    # A common irradiance denominator is applied at each interface, so the
    # Stokes reflectance ratios must reproduce the already-solved radiance
    # ratios without a second solver invocation.
    assert np.isclose(tr['Rrs_Q'][0] / tr['Rrs_I'][0],
                      float(captured['Qu_0plus'][0] /
                            captured['Lu_0plus'][0]),
                      rtol=2e-13, atol=2e-15)
    assert np.isclose(tr['Rrs_U'][0] / tr['Rrs_I'][0],
                      float(captured['Uu_0plus'][0] /
                            captured['Lu_0plus'][0]),
                      rtol=2e-13, atol=2e-15)
    assert np.isclose(tr['rrs_Q'][0] / tr['rrs_I'][0],
                      float(captured['Qu_0minus'][0] /
                            captured['Lu_0minus'][0]),
                      rtol=2e-13, atol=2e-15)
    assert np.isclose(tr['rrs_U'][0] / tr['rrs_I'][0],
                      float(captured['Uu_0minus'][0] /
                            captured['Lu_0minus'][0]),
                      rtol=2e-13, atol=2e-15)
