#!/usr/bin/env python3
"""IPSS 개선 전후 비교 그림.
figI1: 1:1 산포도 — x=IPSS(개선 후, 기준), y=legacy PSSA / 무보정 pp.
figI2: 오차 구조 — VZA·방위별 legacy 편차.
프로젝트 표준 팔레트/레이아웃(figS1~S6과 동일)."""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

D = "/root/ocrt/ipss_bench"
BLUE, ORANGE, GRAY, SURF = "#2a78d6", "#eb6834", "#b9b7ae", "#fcfcfb"


def load(tag):
    return {(round(float(r["vza_deg"]), 3), round(float(r["raa_deg"]), 3)):
            float(r["rho_I"]) for r in csv.DictReader(open(f"{D}/case_{tag}.csv"))}


CASES = [("0.25", "025", "84.26", "84p26"), ("1.00", "100", "84.26", "84p26"),
         ("0.25", "025", "70.47", "70p47")]

# ---------------------------------------------------------------- figI1 (1:1)
fig = plt.figure(figsize=(13.2, 5.6), facecolor=SURF)
gs = fig.add_gridspec(2, 3, height_ratios=[3.2, 1.0], hspace=0.36, wspace=0.30,
                      left=0.06, right=0.985, top=0.85, bottom=0.10)
for j, (tau, tg, sza, sg) in enumerate(CASES):
    pp, lg, ip = load(f"t{tg}_s{sg}_pp"), load(f"t{tg}_s{sg}_legacy"), load(f"t{tg}_s{sg}_ipss")
    ks = sorted(k for k in ip if k[0] <= 70.001 and ip[k] > 0)
    x = np.array([ip[k] for k in ks])
    yl = np.array([lg[k] for k in ks])
    yp = np.array([pp[k] for k in ks])
    ax = fig.add_subplot(gs[0, j]); ax.set_facecolor(SURF)
    rx = fig.add_subplot(gs[1, j]); rx.set_facecolor(SURF)
    ax.scatter(x, yp, s=16, marker="^", facecolors="none", edgecolors=GRAY,
               linewidths=0.8, label="no correction (plane-parallel)")
    ax.scatter(x, yl, s=22, marker="o", facecolors="none", edgecolors=ORANGE,
               linewidths=1.0, label="legacy PSSA (before)")
    lo, hi = min(x.min(), yl.min(), yp.min()), max(x.max(), yl.max(), yp.max())
    pad = 0.06 * (hi - lo); lo -= pad; hi += pad
    ax.plot([lo, hi], [lo, hi], color="0.25", lw=0.9)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal", "box")
    dl = 100 * (yl / x - 1)
    ax.set_title(f"tau={tau}, SZA={sza}$^\\circ$\nlegacy dev  median {np.median(dl):+.2f}%"
                 f"   range {dl.min():+.1f} .. {dl.max():+.1f}%", fontsize=9.5)
    ax.set_xlabel("IPSS  rho_I  (after)", fontsize=8.5)
    if j == 0:
        ax.set_ylabel("legacy / uncorrected  rho_I", fontsize=8.5)
        ax.legend(fontsize=6.8, frameon=False, handletextpad=0.2, loc="upper left")
    ax.tick_params(labelsize=7.5)
    ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(5))
    dp = 100 * (yp / x - 1)
    rx.axhline(0, color="0.25", lw=0.8)
    rx.scatter(x, dp, s=10, marker="^", facecolors="none", edgecolors=GRAY, linewidths=0.7)
    rx.scatter(x, dl, s=12, marker="o", facecolors="none", edgecolors=ORANGE, linewidths=0.8)
    rx.set_xlim(lo, hi)
    rx.set_xlabel("IPSS rho_I", fontsize=8); rx.set_ylabel("dev vs IPSS (%)", fontsize=7.5)
    rx.tick_params(labelsize=7); rx.xaxis.set_major_locator(MaxNLocator(5))
fig.suptitle("PSSA replacement: legacy average-secant Chapman vs IPSS (Zhai & Hu 2022)\n"
             "conservative Rayleigh, black surface, US62 profile, 400 layers, pixel-referenced VZA<=70 x RAA 0-355",
             fontsize=10.5, y=0.995)
fig.savefig(f"{D}/figI1_ipss_1to1.png", dpi=150, facecolor=SURF, bbox_inches="tight")
print("saved figI1")

# ------------------------------------------------------------- figI2 (구조)
fig2 = plt.figure(figsize=(12.4, 4.0), facecolor=SURF)
gs2 = fig2.add_gridspec(1, 3, wspace=0.26, left=0.06, right=0.985, top=0.80, bottom=0.15)
STY = {0.0: (BLUE, "-", "RAA 0 (anti-glint half-plane)"),
       90.0: (GRAY, "--", "RAA 90"),
       180.0: (ORANGE, "-.", "RAA 180 (glint half-plane)")}
for j, (tau, tg, sza, sg) in enumerate(CASES):
    lg, ip = load(f"t{tg}_s{sg}_legacy"), load(f"t{tg}_s{sg}_ipss")
    ax = fig2.add_subplot(gs2[0, j]); ax.set_facecolor(SURF)
    for raa, (c, ls, lab) in STY.items():
        v = sorted({k[0] for k in ip if k[0] <= 70.001})
        y = [100 * (lg[(vz, raa)] / ip[(vz, raa)] - 1) for vz in v]
        ax.plot(v, y, ls, color=c, lw=1.6, label=lab)
    ax.axhline(0, color="0.25", lw=0.8)
    ax.set_title(f"tau={tau}, SZA={sza}$^\\circ$", fontsize=10)
    ax.set_xlabel("VZA at the pixel (deg)", fontsize=8.5)
    if j == 0:
        ax.set_ylabel("legacy PSSA deviation from IPSS (%)", fontsize=8.5)
        ax.legend(fontsize=7.2, frameon=False)
    ax.tick_params(labelsize=8)
fig2.suptitle("Legacy PSSA residual removed by IPSS — small in size but structured: grows with VZA\n"
              "and changes sign across the principal plane, a degree of freedom legacy cannot have",
              fontsize=10.0, y=1.00)
fig2.savefig(f"{D}/figI2_ipss_structure.png", dpi=150, facecolor=SURF, bbox_inches="tight")
print("saved figI2")
