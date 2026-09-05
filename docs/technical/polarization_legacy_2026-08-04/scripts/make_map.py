#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""방향지도를 점이 아닌 연속된 면으로 다시 그린다.

자료는 관측천정각 18개 x 방위각 72개의 완전한 격자다. 빠진 점이 없으므로
격자점 사이를 선형으로 이어 면을 채우면 된다(가우로 음영). 없는 값을
만들어 내는 것이 아니라, 이미 있는 값 사이를 잇는 것이다.
"""
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
rcParams["figure.facecolor"] = "white"

OUT = "/mnt/user-data/outputs/"
SZAS = [25, 50, 75]
Z = pd.read_csv("/tmp/smooth.csv")


def surface(ax, d, col, norm, cmap="RdYlGn", levels=None):
    """격자 자료를 극좌표 면으로 그린다. 방위각은 한 바퀴 닫는다."""
    p = d.pivot_table(index="vza_deg", columns="raa_deg", values=col)
    r = p.index.to_numpy(float)
    th = p.columns.to_numpy(float)
    C = p.to_numpy(float)
    # 355도와 0도를 이어 붙여 빈틈을 없앤다
    th = np.append(th, 360.0)
    C = np.hstack([C, C[:, :1]])
    TH, RR = np.meshgrid(np.radians(th), r)
    m = ax.pcolormesh(TH, RR, C, shading="gouraud", cmap=cmap, norm=norm)
    if levels is not None:
        ax.contour(TH, RR, C, levels=levels, colors="k", linewidths=0.55,
                   alpha=0.45)
    ax.set_theta_zero_location("N")
    ax.set_rmax(85)
    ax.set_rticks([30, 60])
    ax.set_rlabel_position(112)
    ax.tick_params(labelsize=7.5)
    ax.grid(color="k", alpha=0.18, lw=0.6)
    return m


def fig_map():
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 5.1),
                             subplot_kw=dict(projection="polar"))
    norm = mcolors.LogNorm(vmin=1e-3, vmax=2e-2)
    lv = [2e-3, 3e-3, 5e-3, 8e-3, 1.2e-2]
    for ax, sza in zip(axes, SZAS):
        d = Z[Z.sza_deg == sza].copy()
        d["val"] = d.med_eP / 3.0
        m = surface(ax, d, "val", norm, levels=lv)
        ax.set_title("태양천정각 %d°" % sza, fontsize=11.5, pad=14)
    cb = fig.colorbar(m, ax=axes, fraction=0.022, pad=0.04)
    cb.set_label("필요한 선형편광도 정확도")
    cb.set_ticks([1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
    cb.ax.set_yticklabels(["0.001", "0.002", "0.005", "0.010", "0.020"])
    fig.suptitle("관측 방향별로 본 요구 정확도.  초록일수록 잡음에 관대해 유리하다\n"
                 "반경 = 관측천정각 (0~85°), 각도 = 상대방위각.  "
                 "관측천정각 18개 x 방위각 72개 실측값을 면으로 이었다",
                 fontsize=12, y=1.06)
    fig.savefig(OUT + "F4_방향지도.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_map_hard():
    """어려운 쪽 짝 기준 지도도 함께."""
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 9.4),
                             subplot_kw=dict(projection="polar"))
    specs = [("med_eP", "짝 절반 기준", mcolors.LogNorm(1e-3, 2e-2),
              [2e-3, 3e-3, 5e-3, 8e-3, 1.2e-2]),
             ("q25_eP", "어려운 쪽 짝 기준", mcolors.LogNorm(3e-4, 6e-3),
              [5e-4, 1e-3, 2e-3, 4e-3])]
    for i, (col, lab, norm, lv) in enumerate(specs):
        for j, sza in enumerate(SZAS):
            ax = axes[i, j]
            d = Z[Z.sza_deg == sza].copy()
            d["val"] = d[col] / 3.0
            m = surface(ax, d, "val", norm, levels=lv)
            if i == 0:
                ax.set_title("태양천정각 %d°" % sza, fontsize=11.5, pad=14)
            if j == 0:
                ax.text(-0.30, 0.5, lab, transform=ax.transAxes, rotation=90,
                        va="center", fontsize=12)
        cb = fig.colorbar(m, ax=axes[i, :], fraction=0.022, pad=0.04)
        cb.set_label("필요한 선형편광도 정확도")
    fig.suptitle("관측 방향별 요구 정확도 — 위는 짝의 절반을 가르는 기준, "
                 "아래는 가장 어려운 짝까지 가르는 기준\n"
                 "초록일수록 유리하다.  반경 = 관측천정각, 각도 = 상대방위각",
                 fontsize=12, y=0.98)
    fig.savefig(OUT + "F5_방향지도_두기준.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


fig_map()
print("F4 완료")
fig_map_hard()
print("F5 완료")
