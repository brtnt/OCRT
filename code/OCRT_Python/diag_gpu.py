"""GPU 진단: 백엔드가 실제 CuPy인지, R1 1케이스 정밀조건이 몇 초인지 확인.

실행:  set OCRT_PY_GPU=1  &&  python diag_gpu.py
실행 중 다른 창에서 nvidia-smi 를 한 번 더 찍어 GPU-Util 을 관찰한다.
"""
import os
os.environ.setdefault('OCRT_PY_GPU', '1')
import sys
sys.path.insert(0, '.')
import time

from ocrt_py import backend

print('=' * 55)
print('backend.GPU =', backend.GPU, '  (True 여야 GPU 사용)')
print('backend.xp  =', backend.xp.__name__, '  (cupy 여야 GPU 사용)')
xp = backend.xp
# CuPy 가 실제 GPU 연산을 하는지 (이 줄에서 GPU 컨텍스트가 생김)
t = time.perf_counter()
a = xp.arange(5_000_000.0)
s = float((a * a).sum())
print(f'GPU 워밍업 연산: sum={s:.3e}  ({(time.perf_counter()-t)*1000:.0f}ms)')
print('=' * 55)

from ocrt_py import batch_driver as BD, aerosol as AER
from ocrt_py.constituent import OCRTConstituentModel

cm = OCRTConstituentModel('data', 'micro', 'red_clay')
mie = {'C50': AER.read_mie('data/C50.mie')}
c = dict(sza=32, vza=28, raa=95, wl=555, wind=3, aer='C50', aod=0.15,
         chl=0.3, tsm=0.5, ad=0.03)

# ---- CPU 부분(위상 계수)만 따로 계측 -------------------------------------
from ocrt_py.aerosol import water_component_phase_moments_gauss
r = cm.evaluate(555, 0.3, 0.5, 0.03)
t = time.perf_counter()
for mc in [cm.det_mie, cm.min_mie]:
    water_component_phase_moments_gauss(mc, 555, 200, 400)
print(f'[CPU] 수중 위상 계수 2 산란성분(nmg=400): {(time.perf_counter()-t)*1000:.0f}ms')

# ---- R1 1케이스 정밀조건 전체 --------------------------------------------
print('R1 1케이스 정밀조건 실행 중... (지금 nvidia-smi 관찰)')
t = time.perf_counter()
rho, rrs = BD.solve_r1_grid(
    [dict(c)], mie, cm, L_max=80, Lmix=200,
    n_mu_water=24, nt_atm=400, fourier_m_max=16,
    max_it_atm=100, tol_atm=1e-7, max_it_water=500, tol_water=1e-7,
    n_water=1.34, pressure_hpa=1013.25, F_sun_TOA=1.0)
dt = time.perf_counter() - t
print(f'[전체] R1 1케이스: {dt:.1f}s   rho={rho[0]:.6e}')
print('=' * 55)
print(f'추정: 12000행(1000x12밴드) x 3회(R1/R2/R3) ~ {dt*12000*3/3600:.1f}시간 (R1 기준 상한)')
