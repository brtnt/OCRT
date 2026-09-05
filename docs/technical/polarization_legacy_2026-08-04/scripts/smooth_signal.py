#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
매끄러운 최종 지표.

앞선 계산이 방위각을 따라 튀었던 이유는 광학두께를 격자로 훑고 그 위에서
최솟값을 고른 탓이다. 격자가 조금만 어긋나도 최솟값이 다른 점으로 건너뛴다.

여기서는 격자에서 고르지 않고 풀어서 구한다.

  두 모델의 748 nm 레일리 보정 반사도가 같아지도록 각자의 광학두께를
  맞춘다. 이것은 한 변수를 다른 변수로 나타내는 관계이므로 곡선 하나가
  된다. 그 곡선을 따라가며 865 nm 세기 차이와 선형편광도 차이를 잰다.

  e_I : 곡선 위에서 남는 865 nm 세기 차이의 최솟값
  e_P : 그 지점에서의 선형편광도 차이

두 값 모두 연속량이라 기하가 조금 변하면 조금만 변한다.
잡음은 가정하지 않는다. 필요한 정확도는 e/3 으로 읽으면 된다.
"""
import glob
import os
import numpy as np
import pandas as pd

NIR = [748, 865]
WIND = 5
FFS = [20, 30, 50, 80]
RHS = [30, 50, 80, 95]
NFINE = 1200
REP = ("r30f50v01", "r95f50v01")

lib = pd.read_csv("/tmp/aerlib.csv").set_index("model")


def load_ray():
    out = []
    for f in sorted(glob.glob("/tmp/ray/R_*.csv")):
        b, s = os.path.basename(f)[2:-4].split("_")
        d = pd.read_csv(f)
        d["band_nm"] = int(b); d["sza_deg"] = int(s)
        out.append(d[["band_nm", "sza_deg", "vza_deg", "raa_deg",
                      "rho_I", "rho_Q", "rho_U"]])
    r = pd.concat(out, ignore_index=True)
    r.columns = ["band_nm", "sza_deg", "vza_deg", "raa_deg", "rI", "rQ", "rU"]
    return r


R = load_ray()
A = pd.read_csv("atm_reference.csv.gz")
A = A[A.band_nm.isin(NIR) & (A.wind_ms == WIND)]
A["ff"] = A.aer_model.map(lib.ff); A["rh"] = A.aer_model.map(lib.rh)
D = A[A.ff.isin(FFS) & A.rh.isin(RHS)].copy()

MODELS = sorted(D.aer_model.unique())
NM = len(MODELS)
MI = {m: i for i, m in enumerate(MODELS)}
AODS = np.array(sorted(D.aod865.unique()))
LA = np.log(AODS)
LAD = np.linspace(LA.min(), LA.max(), NFINE)
FF = np.array([lib.ff[m] for m in MODELS])
SSA = np.array([lib.ssa865[m] for m in MODELS])
PAIRS = [(i, j) for i in range(NM) for j in range(i + 1, NM) if FF[i] == FF[j]]
DSSA = np.array([abs(SSA[i] - SSA[j]) for (i, j) in PAIRS])
REP_K = [k for k, (i, j) in enumerate(PAIRS)
         if {MODELS[i], MODELS[j]} == set(REP)][0]
print("모델 %d종, 흡수만 다른 짝 %d개, 대표 짝 번호 %d" % (NM, len(PAIRS), REP_K))

Rk = R.set_index(["band_nm", "sza_deg", "vza_deg", "raa_deg"]).rI
rows = []
for gi, ((sza, vza, raa), g) in enumerate(D.groupby(["sza_deg", "vza_deg",
                                                     "raa_deg"])):
    try:
        rI = np.array([Rk[(b, sza, vza, raa)] for b in NIR])
    except KeyError:
        continue
    piv = g.pivot_table(index="aer_model", columns=["aod865", "band_nm"],
                        values=["atm_I", "atm_Q", "atm_U"])
    if len(piv) != NM:
        continue
    Y = np.empty((NM, len(AODS), 4))
    ok = True
    for ai, a in enumerate(AODS):
        try:
            I = np.column_stack([piv[("atm_I", a, b)] for b in NIR])
            Qv = np.column_stack([piv[("atm_Q", a, b)] for b in NIR])
            Uv = np.column_stack([piv[("atm_U", a, b)] for b in NIR])
        except KeyError:
            ok = False
            break
        Y[:, ai, 0:2] = I - rI
        Y[:, ai, 2:4] = np.hypot(Qv, Uv) / I
    if not ok:
        continue
    Yd = np.empty((NM, NFINE, 4))
    for m in range(NM):
        for k in range(4):
            Yd[m, :, k] = np.interp(LAD, LA, Y[m, :, k])

    eI, eP = np.full(len(PAIRS), np.nan), np.full(len(PAIRS), np.nan)
    for k, (i, j) in enumerate(PAIRS):
        xj = Yd[j][:, 0]
        if xj[-1] <= xj[0]:
            continue
        # 748 nm 를 같게 만드는 j 의 광학두께 (한 값으로 정해진다)
        lj = np.interp(Yd[i][:, 0], xj, LAD, left=np.nan, right=np.nan)
        good = np.isfinite(lj)
        if good.sum() < 5:
            continue
        r865 = np.abs(Yd[i][:, 1] -
                      np.interp(lj, LAD, Yd[j][:, 1], left=np.nan, right=np.nan))
        r865 = np.where(good, r865, np.nan)
        t = int(np.nanargmin(r865))
        eI[k] = r865[t]
        pj = np.array([np.interp(lj[t], LAD, Yd[j][:, c]) for c in (2, 3)])
        eP[k] = float(np.sqrt(((Yd[i][t, 2:4] - pj) ** 2).sum()))
    if not np.isfinite(eP).any():
        continue
    rows.append(dict(sza_deg=sza, vza_deg=vza, raa_deg=raa,
                     rep_eI=eI[REP_K], rep_eP=eP[REP_K],
                     med_eI=float(np.nanmedian(eI)),
                     med_eP=float(np.nanmedian(eP)),
                     q25_eP=float(np.nanpercentile(eP, 25)),
                     slope=float(np.nanmedian(eP / DSSA))))
    if (gi + 1) % 1000 == 0:
        print("  %d 기하" % (gi + 1), flush=True)

Z = pd.DataFrame(rows)
Z.to_csv("/tmp/smooth.csv", index=False)
print("기하 %d개 완료" % len(Z))
print()
print("=== 방위각 변화의 매끄러움 (태양천정각 50도, 관측천정각 40도) ===")
d = Z[(Z.sza_deg == 50) & (Z.vza_deg == 40)].sort_values("raa_deg")
for c, nm in [("rep_eP", "대표 짝 편광 신호"), ("rep_eI", "대표 짝 세기 잔여"),
              ("med_eP", "짝 전체 중앙값")]:
    v = d[c].values
    F = np.abs(np.fft.rfft(v - v.mean()))
    print("  %-16s 9차 이상 %.1f%%, 이웃간 변화 중앙 %.1f%%"
          % (nm, 100 * F[9:].sum() / F[1:].sum(),
             100 * np.median(np.abs(np.diff(v)) / np.abs(v[:-1]))))
print()
print("=== 세기를 맞췄을 때 남는 신호 (전 기하 중앙값) ===")
print("  대표 짝: 세기 잔여 %.2e, 편광 차이 %.2e"
      % (Z.rep_eI.median(), Z.rep_eP.median()))
print("  짝 전체: 세기 잔여 %.2e, 편광 차이 %.2e"
      % (Z.med_eI.median(), Z.med_eP.median()))
print()
print("=== 3시그마로 가르는 데 필요한 정확도 (전 기하 중앙값) ===")
print("%-12s %14s %14s %14s" % ("태양천정각", "세기 필요", "편광 필요", "비"))
for s in [25, 50, 75]:
    d2 = Z[Z.sza_deg == s]
    a = d2.med_eI.median() / 3
    b = d2.med_eP.median() / 3
    print("%-12d %14.2e %14.2e %13.0f배" % (s, a, b, b / a))
print()
print("=== 방위각별 필요 편광 정확도 (관측천정각 30~60도, 짝 전체 중앙값) ===")
sub = Z[(Z.vza_deg >= 30) & (Z.vza_deg <= 60)].copy()
sub["raa180"] = np.minimum(sub.raa_deg, 360 - sub.raa_deg)
sub["bin"] = (sub.raa180 // 30 * 30).astype(int)
p = sub.pivot_table(index="bin", columns="sza_deg", values="med_eP") / 3
print(p.map(lambda x: "%.2e" % x).to_string())
