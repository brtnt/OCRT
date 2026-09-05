#!/usr/bin/env python3
"""RTSOS 교차검증 케이스 생성.
대기: 보존 Rayleigh, 지수 8 km, tau_total 0.0935, TOA 100 km, 등-Δτ 40층 (OCRT 방식),
      흑면 (Lambertian albedo 0).  위상행렬은 동봉 ray.pmtx (κ 는 위상함수 무관).
각도: SZA {0,40,60,70,80,85}, VZA {0,20,40,55,70}, PHI {0,90,180}.
같은 층 프로파일을 profile.csv 로 내보내 OCRT rt_ipss 측 κ 계산에 쓴다."""
import numpy as np, os, shutil

TAU, HS, HTOA, N = 0.0935, 8.0, 100.0, 40
SZA = [0.0, 40.0, 60.0, 70.0, 80.0, 85.0]
VZA = [0.0, 20.0, 40.0, 55.0, 70.0]
PHI = [0.0, 90.0, 180.0]
D = "/root/rtsos/xcheck/case"
os.makedirs(D, exist_ok=True)

# 등-Δτ 층 경계 (지수 프로파일의 누적 질량 등분)
F = 1.0 - np.exp(-HTOA / HS)
edges = [-HS * np.log(1.0 - (i / N) * F) for i in range(N + 1)]   # 0 .. 100 (상향)
edges[-1] = HTOA
tops = edges[1:][::-1]                                            # 층 상단, 위에서 아래로
bots = edges[:-1][::-1]

with open(f"{D}/sosi.dat", "w") as f:
    f.write("-1000  0.0  0.0  0.0  0.0  0.0  0.0\n")               # TOA 검출기
    for t, b in zip(tops, bots):
        f.write(f"1  {TAU/N:.10e}  1.0  0.0  0.0  0.0  {t:.6f}\n")   # IPT>0: TAU LBDOM wl T EFLRSC ALH
    f.write("-200  0.0  1.0  0.0  0.0  0.0  0.0\n")               # 흑면 Lambertian
    f.write("40  30  62\n1  3  3\nray.pmtx\n")
with open(f"{D}/sosi.amu", "w") as f:
    f.write(f"{len(SZA)}\n"); [f.write(f"{s:.4f}\n") for s in SZA]
    f.write(f"{len(VZA)}  {len(PHI)}\n"); [f.write(f"{v:.4f}\n") for v in VZA]
    [f.write(f"{p:.1f}\n") for p in PHI]
shutil.copy("/root/rtsos/RTSOS/validation/benchmark/Coulson_thick/ray.pmtx", D)
shutil.copy("/root/rtsos/RTSOS/validation/benchmark/Coulson_thick/auxiliary_directory", D)
with open(f"{D}/profile.csv", "w") as f:
    f.write("h_bot_km,h_top_km,tau,omega\n")
    for t, b in zip(tops, bots):
        f.write(f"{b:.6f},{t:.6f},{TAU/N:.10e},1.0\n")
print("layers", N, "bottom-layer top", tops[-1], "km;  top layer", bots[0], "-", tops[0])
