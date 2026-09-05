#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""최종 그림 — 잡음 가정 없이, 선형화 없이."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

rcParams["font.family"] = "Noto Sans CJK JP"
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 10.5
rcParams["axes.grid"] = True
rcParams["grid.alpha"] = 0.25
rcParams["figure.facecolor"] = "white"
rcParams["axes.axisbelow"] = True

OUT = "/mnt/user-data/outputs/"
SZAS = [25, 50, 75]
COL = {25: "#3B4CC0", 50: "#2E8B57", 75: "#D9534F"}

lib = pd.read_csv("/tmp/aerlib.csv")
Z = pd.read_csv("/tmp/smooth.csv")
A = pd.read_csv("atm_reference.csv.gz")
Z["raa180"] = np.minimum(Z.raa_deg, 360 - Z.raa_deg)


# =============================================================== 그림 1
def figA():
    """축퇴 구조와 대표 짝."""
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0))
    G = lib[lib.ff.isin([20, 30, 50, 80]) & lib.rh.isin([30, 50, 80, 95])]
    ax = axes[0]
    ax.scatter(lib.ang748_865, lib.ssa865, s=20, c="#D5D5D5",
               edgecolors="none", zorder=1, label="라이브러리 80종")
    cmap = plt.get_cmap("plasma")
    for k, ff in enumerate([20, 30, 50, 80]):
        d = G[G.ff == ff].sort_values("rh")
        ax.plot(d.ang748_865, d.ssa865, "o-", color=cmap(k / 3.4), ms=9,
                lw=2.2, zorder=3, label="미세입자 %d%%" % ff)
        for _, r in d.iterrows():
            ax.annotate("%d" % r.rh, (r.ang748_865, r.ssa865),
                        textcoords="offset points", xytext=(7, -3),
                        fontsize=7.2, color=cmap(k / 3.4))
    rep = G[G.model.isin(["r30f50v01", "r95f50v01"])]
    ax.plot(rep.ang748_865, rep.ssa865, "o", ms=17, mfc="none",
            mec="k", mew=2.0, zorder=5)
    ax.annotate("대표 짝", (rep.ang748_865.mean(), rep.ssa865.mean()),
                textcoords="offset points", xytext=(-58, 0), fontsize=9.5,
                weight="bold")
    ax.set_xlabel("옹스트롬 지수 (748–865 nm)  ← 입자 크기")
    ax.set_ylabel("단일산란알베도 (865 nm)  ← 흡수")
    ax.set_title("축퇴 격자 16종.  점 옆 숫자는 상대습도", fontsize=11)
    ax.legend(fontsize=8, loc="lower left", ncol=2)

    # 대표 짝의 관측량을 방위각을 따라
    ax = axes[1]
    for m, lab, c in [("r30f50v01", "흡수 강함 (알베도 0.947)", "#C0392B"),
                      ("r95f50v01", "흡수 약함 (알베도 0.984)", "#1A5FB4")]:
        d = A[(A.band_nm == 865) & (A.sza_deg == 50) & (A.vza_deg == 40) &
              (A.wind_ms == 5) & (A.aer_model == m) &
              (A.aod865 == 0.1)].sort_values("raa_deg")
        ax.plot(d.raa_deg, np.hypot(d.atm_Q, d.atm_U) / d.atm_I, "-",
                color=c, lw=2.4, label=lab)
    ax.set_xlim(0, 355)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel("상대방위각 (도)")
    ax.set_ylabel("선형편광도 865 nm")
    ax.set_title("두 모델은 세기 기울기가 같지만\n편광은 방위각 전체에서 갈린다",
                 fontsize=11)
    ax.legend(fontsize=8.8, loc="upper left")
    fig.suptitle("입자 크기가 같고 흡수만 다른 축퇴 짝  "
                 "(태양천정각 50°, 관측천정각 40°, 광학두께 0.1)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT + "F1_축퇴격자와_대표짝.png", dpi=145, bbox_inches="tight")
    plt.close(fig)


# =============================================================== 그림 2
def figB():
    """필요 정확도 — 세기와 편광."""
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
    ax = axes[0]
    x = np.arange(len(SZAS)); w = 0.36
    a = [Z[Z.sza_deg == s].med_eI.median() / 3 for s in SZAS]
    b = [Z[Z.sza_deg == s].med_eP.median() / 3 for s in SZAS]
    ax.bar(x - w / 2, a, w, color="#9E9E9E", label="세기 반사도에 필요")
    ax.bar(x + w / 2, b, w, color="#1A5FB4", label="선형편광도에 필요")
    ax.axhline(5e-4, color="#C62828", lw=1.6, ls="--")
    ax.text(2.44, 5.6e-4, "세기 반사도의 현실적 정확도\n(보정 2 %, 약 5e-4)",
            fontsize=8, color="#C62828", ha="right", va="bottom")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(["%d°" % s for s in SZAS])
    ax.set_xlabel("태양천정각")
    ax.set_ylabel("3시그마로 가르는 데 필요한 정확도")
    ax.set_title("흡수만 다른 짝을 가르려면\n얼마나 정확해야 하는가", fontsize=11)
    ax.legend(fontsize=8.8, loc="center left")

    ax = axes[1]
    r = [Z[Z.sza_deg == s].med_eP.median() / Z[Z.sza_deg == s].med_eI.median()
         for s in SZAS]
    ax.bar(x, r, 0.5, color="#2E8B57")
    for xi, v in zip(x, r):
        ax.text(xi, v * 1.04, "%.0f배" % v, ha="center", fontsize=10)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(["%d°" % s for s in SZAS])
    ax.set_xlabel("태양천정각")
    ax.set_ylabel("요구 정확도가 느슨해지는 배수")
    ax.set_title("편광을 쓰면 요구 사양이\n몇 배 느슨해지는가", fontsize=11)
    fig.suptitle("잡음을 가정하지 않고, 필요한 정확도를 거꾸로 구한 결과  "
                 "(3,888개 관측 방향의 중앙값)", fontsize=12, y=1.03)
    fig.tight_layout()
    fig.savefig(OUT + "F2_필요정확도.png", dpi=145, bbox_inches="tight")
    plt.close(fig)


# =============================================================== 그림 3
def figC():
    """방위각 의존."""
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    sub = Z[(Z.vza_deg >= 30) & (Z.vza_deg <= 60)]
    for ax, (c, lab) in zip(axes[:2],
                            [("med_eP", "짝 전체의 중앙값"),
                             ("q25_eP", "어려운 쪽 짝 (하위 사분위)")]):
        for s in SZAS:
            d = sub[sub.sza_deg == s].groupby("raa180")[c].median()
            ax.plot(d.index, d.values / 3, "-", color=COL[s], lw=2.4,
                    label="태양천정각 %d°" % s)
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 45, 90, 135, 180])
        ax.set_xlabel("상대방위각 (도)")
        ax.set_ylabel("필요한 선형편광도 정확도")
        ax.set_title("%s\n(위로 갈수록 잡음에 관대하다)" % lab, fontsize=10.8)
    axes[0].legend(fontsize=8.8, loc="upper left")
    ax = axes[2]
    for s in SZAS:
        d = sub[sub.sza_deg == s].groupby("raa180").med_eI.median()
        ax.plot(d.index, d.values / 3, "-", color=COL[s], lw=2.4)
    ax.axhline(5e-4, color="#C62828", lw=1.6, ls="--")
    ax.text(178, 5.6e-4, "현실적 정확도", fontsize=8, color="#C62828",
            ha="right", va="bottom")
    ax.set_yscale("log")
    ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180])
    ax.set_xlabel("상대방위각 (도)")
    ax.set_ylabel("필요한 세기 반사도 정확도")
    ax.set_title("세기 쪽은 어느 방위각에서도\n현실적 정확도에 못 미친다",
                 fontsize=10.8)
    fig.suptitle("어느 방위각이 유리한가  (관측천정각 30~60°)   "
                 "0°는 태양과 같은 쪽, 180°는 태양 반대쪽",
                 fontsize=12, y=1.03)
    fig.tight_layout()
    fig.savefig(OUT + "F3_방위각의존.png", dpi=145, bbox_inches="tight")
    plt.close(fig)


# F4 와 F5(관측 방향 지도)는 make_map.py 가 면 형식으로 만든다.
# 여기서 점 형식으로 다시 그리면 그 결과를 덮어쓰므로 넣지 않는다.
for f in (figA, figB, figC):
    f()
    print("완료", f.__name__)
