#!/usr/bin/env python3
"""IPSS 교체 전후 비교 — 사용자 지정 구성 (v1.11.1: 지상(화소) 앵커 + 위상함수 가중 κ).

고정: VZA 55, RAA 90, M80C AOD(555)=0.1, black_fresnel_ocean, 풍속 5 m/s,
      US62, 400층, 기체흡수 ON.  스윕: SZA 0..85.
기준선(0 %) = IPSS.  비교 대상 = 무보정 평면평행, 그리고 v1.11 에서 삭제된
구 평균할선 Chapman PSSA(교체 전 바이너리로 산출), 그리고 v1.11 의 위상함수 미가중 κ
(2026-09-05 재조사에서 발견·수정; sza_sweep_v1p11_phasefree/).
상단 = 전 구간, 하단 = ±1 % 확대(저·중 SZA 에서 보정량이 0 에 수렴함을 확인)."""
import math
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GRAY, SURF = "#2a78d6", "#eb6834", "#b9b7ae", "#fcfcfb"
D = "/root/ocrt/ipss_bench"
NEW, OLD, PF = f"{D}/sza_sweep", f"{D}/sza_sweep_toaanchor", f"{D}/sza_sweep_v1p11_phasefree"
SZA = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85]
WLS = [(412, 0.3186, 0.1058), (555, 0.0935, 0.1000), (865, 0.0155, 0.0935)]


def rho(d, wl, sza, mode):
    f = f"{d}/wl{wl}_sza{sza}_{mode}.txt"
    m = re.search(r"TOA_rho_I=(-?[0-9.eE+-]+)", open(f).read())
    return float(m.group(1))


fig = plt.figure(figsize=(13.0, 7.0), facecolor=SURF)
gs = fig.add_gridspec(2, 3, height_ratios=[1.3, 1.0], hspace=0.33, wspace=0.24,
                      left=0.065, right=0.985, top=0.845, bottom=0.085)
x = np.array(SZA, float)

for j, (wl, tR, aod) in enumerate(WLS):
    ip = np.array([rho(NEW, wl, s, "ipss") for s in SZA])
    pp = np.array([rho(NEW, wl, s, "pp") for s in SZA])
    lg = np.array([rho(OLD, wl, s, "legacy") for s in SZA])
    pf = np.array([rho(PF, wl, s, "ipss") for s in SZA])
    dp, dl, df = 100 * (pp / ip - 1), 100 * (lg / ip - 1), 100 * (pf / ip - 1)

    for row in (0, 1):
        ax = fig.add_subplot(gs[row, j])
        ax.set_facecolor(SURF)
        ax.axhline(0, color=BLUE, lw=2.0, zorder=2)
        ax.plot(x, dp, "--^", color=GRAY, mfc="none", ms=5, lw=1.5,
                label="no correction (plane-parallel)")
        ax.plot(x, dl, "-o", color=ORANGE, mfc="none", ms=5, lw=1.8,
                label="legacy Chapman PSSA (deleted in v1.11)")
        ax.plot(x, df, ":s", color="#6a5acd", mfc="none", ms=4, lw=1.4,
                label=r"IPSS v1.11, phase-free $\kappa$ (fixed in v1.11.1)")
        ax.set_xlim(-3, 88)
        ax.set_xticks(range(0, 90, 15))
        ax.tick_params(labelsize=8)
        ax.set_xlabel("SZA (deg)", fontsize=9)
        if row == 0:
            lo = min(dp.min(), dl.min(), df.min())
            ax.set_ylim(lo - 0.12 * abs(lo), max(0.9, 0.22 * abs(lo)))
            ax.set_title(f"{wl} nm   " r"$\tau_R$=" f"{tR:.4f}, "
                         r"$\tau_a$=" f"{aod:.4f}", fontsize=10)
            ax.text(0.965, 0.90, "IPSS = 0 % (reference)", transform=ax.transAxes,
                    color=BLUE, fontsize=8.0, fontweight="bold", ha="right", va="top")
            if j == 0:
                ax.set_ylabel(r"$\rho_I$  deviation from IPSS (%)", fontsize=9.5)
                ax.legend(fontsize=7.4, frameon=False, loc="lower left",
                          bbox_to_anchor=(-0.012, 0.02))
            ax.annotate(f"{dp[-1]:+.1f} %", xy=(85, dp[-1]),
                        xytext=(0.52, 0.30), textcoords="axes fraction",
                        fontsize=7.6, color="0.35", ha="left", va="center",
                        arrowprops=dict(arrowstyle="->", color="0.45", lw=0.9))
        else:
            ax.set_ylim(-1.0, 1.0)
            ax.axhspan(-0.2, 0.2, color=BLUE, alpha=0.06, zorder=0)
            ax.set_title("same data, $\\pm$1 % zoom", fontsize=8.5, color="0.35")
            if j == 0:
                ax.set_ylabel("deviation (%)  [zoom]", fontsize=9.5)
                ax.text(0.035, 0.075,
                        "shaded band = $\\pm$0.2 %\nno correction needed below SZA $\\approx$ 70",
                        transform=ax.transAxes, fontsize=7.0, color="0.35")

fig.suptitle(
    "Spherical correction, before / after — VZA 55$^\\circ$ (pixel-referenced), RAA 90$^\\circ$, "
    "M80C AOD(555)=0.1, black Fresnel ocean, wind 5 m s$^{-1}$\n"
    "US62 profile, 400 layers, gas absorption on;  IPSS (Zhai & Hu 2022) v1.11.1 — pixel-anchored VZA, "
    r"phase-weighted $\kappa$ — is the 0 % reference", fontsize=10.5, y=0.985)
fig.savefig(f"{D}/figI3_sza_diff.png", dpi=160, facecolor=SURF, bbox_inches="tight")
print("saved figI3_sza_diff.png")

print(f"\n{'wl':>5} {'pp@65':>8} {'pp@70':>8} {'pp@75':>8} {'pp@80':>8} {'pp@85':>8}"
      f" | {'lg@65':>8} {'lg@85':>8}")
for wl, _, _ in WLS:
    g = lambda m, s, d: 100 * (rho(d, wl, s, m) / rho(NEW, wl, s, "ipss") - 1)
    print(f"{wl:5d} " + " ".join(f"{g('pp',s,NEW):+8.3f}" for s in (65, 70, 75, 80, 85))
          + f" | {g('legacy',65,OLD):+8.3f} {g('legacy',85,OLD):+8.3f}")
