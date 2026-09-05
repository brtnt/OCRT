#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""설치 검증: 완전 ρ_TOA(축소조건)를 산출해 C 기준값과 대조한다."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ocrt_solve import solve_reflectances

C_REF = 2.2748060400e-01   # C --water-model ocrt, 축소조건 n_mu_water=16 nt=100 m_max=3
print("설치 검증 실행(축소조건, ~1분)...")
t = time.time()
r = solve_reflectances(555.0, 30.0, 30.0, 90.0, 3.0, 1.2, 3.5, 0.04, 0.18,
                       'C50', which=('R1',), n_mu_water=16, fourier_m_max=3,
                       nt_atm=100, max_it_water=500)
dt = time.time() - t
rel = abs(r['rho_TOA_I'] - C_REF) / C_REF
print(f"rho_TOA_I = {r['rho_TOA_I']:.10e}  (C {C_REF:.10e}, rel {rel:.2e})")
print(f"orders={r['orders']} conv={r['conv']}  [{dt:.0f}s]")
if rel < 1e-8 and r['conv'] == 1:
    print("PASS: C와 일치, 정상 설치.")
else:
    print("FAIL: 기준값 불일치 또는 미수렴. 환경/데이터 확인 필요.")
    sys.exit(1)
