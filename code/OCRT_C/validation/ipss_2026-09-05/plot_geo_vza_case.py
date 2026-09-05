#!/usr/bin/env python3
"""실제 기하 계산 그림 (2026-09-05).

표적 36.0N 128.25E, 센서 0.0N 128.25E (같은 경도, 적도).
(a) Delta theta 는 theta_surface 만의 함수 — 700 km / 2000 / 5000 / GEO 가 한 점에 겹친다.
(b) gamma = 36 deg 고정 시 센서 고도를 올리면 theta_surface 가 gamma 로 수렴하고
    Delta theta 도 9.18 -> 0.64 deg 로 줄어든다 (0 으로 가지는 않는다).
프로젝트 표준 팔레트 (dataviz 검증: #2a78d6/#eb6834 CVD dE 24.7, normal 33.6 — PASS).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

_KO = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
fm.fontManager.addfont(_KO)
matplotlib.rcParams["font.family"] = [fm.FontProperties(fname=_KO).get_name(), "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

BLUE, ORANGE, GRAY, SURF = "#2a78d6", "#eb6834", "#b9b7ae", "#fcfcfb"
INK, INK2, MUTED = "#1c1c1a", "#4a4a46", "#8a8a84"
RE, H, GAMMA = 6371.0, 100.0, 36.0
D, R = np.degrees, np.radians


def th_surf(gamma, h_sat):
    g = R(gamma); Rs = RE + h_sat
    d = np.sqrt(RE ** 2 + Rs ** 2 - 2 * RE * Rs * np.cos(g))
    return D(np.arccos(np.clip((Rs * np.cos(g) - RE) / d, -1.0, 1.0)))


def th_toa(t):
    return D(np.arcsin(np.clip(RE / (RE + H) * np.sin(R(t)), -1.0, 1.0)))


def gamma_for(t_target, h_sat):
    lo, hi = 1e-9, D(np.arccos(RE / (RE + h_sat))) - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if th_surf(mid, h_sat) < t_target: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


T_GEO = th_surf(GAMMA, 35786.0)

fig = plt.figure(figsize=(12.8, 5.0), facecolor=SURF)
gs = fig.add_gridspec(1, 2, wspace=0.24, left=0.055, right=0.985, top=0.80, bottom=0.16)

# ---------------------------------------------------------------- (a)
ax = fig.add_subplot(gs[0, 0]); ax.set_facecolor(SURF)
t = np.linspace(0, 85, 500)
ax.plot(t, t - th_toa(t), "-", color=BLUE, lw=2.4)
for h_sat, lab in ((700.0, "700 km"), (2000.0, "2 000 km"),
                   (5000.0, "5 000 km"), (35786.0, "GEO 35 786 km")):
    g = gamma_for(T_GEO, h_sat)
    ax.plot(T_GEO, T_GEO - th_toa(T_GEO), "o", color=ORANGE, ms=11,
            mfc="none", mew=2.0, zorder=5)
ax.annotate("센서 고도 700 km / 2 000 / 5 000 / GEO —\n"
            "네 경우가 정확히 같은 한 점에 겹친다\n"
            r"($\theta_{\rm surface}$=41.779$^\circ$ $\Rightarrow$ "
            r"$\Delta\theta$=0.786$^\circ$)",
            xy=(T_GEO, T_GEO - th_toa(T_GEO)), xytext=(0.10, 0.62),
            textcoords="axes fraction", fontsize=8.4, color=ORANGE,
            fontweight="bold", ha="left",
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2))
ax.set_xlim(0, 86); ax.set_ylim(0, 6.0)
ax.set_xlabel(r"$\theta_{\rm surface}$  =  L1B senz  (deg)", fontsize=9.5)
ax.set_ylabel(r"$\Delta\theta=\theta_{\rm surface}-\theta_{\rm TOA}$  (deg)", fontsize=9.5)
ax.tick_params(labelsize=8.4)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.set_title(r"a) $\Delta\theta$ 는 $\theta_{\rm surface}$ 만의 함수다 (H=100 km)",
             fontsize=10.0, color=INK, fontweight="bold", loc="left", pad=10)
ax.text(0.0, -0.175, "센서 고도는 어디에도 들어가지 않는다. 같은 senz 면 700 km 든 GEO 든\n"
        r"$\Delta\theta$ 가 완전히 동일하다 — 이것이 앵커 규약이 고도와 무관한 이유다.",
        transform=ax.transAxes, fontsize=7.9, color=INK2, va="top", linespacing=1.5)

# ---------------------------------------------------------------- (b)
bx = fig.add_subplot(gs[0, 1]); bx.set_facecolor(SURF)
hmin = RE / np.cos(R(GAMMA)) - RE            # 1504 km
hs = np.geomspace(hmin * 1.02, 3.0e5, 600)
ts = np.array([th_surf(GAMMA, h) for h in hs])
bx.plot(hs, ts, "-", color=BLUE, lw=2.4, label=r"$\theta_{\rm surface}$ (senz)")
bx.plot(hs, ts - th_toa(ts), "-", color=ORANGE, lw=2.4, label=r"$\Delta\theta$")
bx.axhline(GAMMA, color=BLUE, lw=1.1, dashes=(5, 3))
bx.annotate(r"$h\to\infty$ 극한: $\theta_{\rm surface}\to\gamma$=36.0$^\circ$",
            xy=(1.6e5, GAMMA), xytext=(0.44, 0.56), textcoords="axes fraction",
            fontsize=8.0, color=BLUE, ha="left",
            arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))
bx.axhline(0.6412, color=ORANGE, lw=1.1, dashes=(5, 3))
bx.annotate(r"$\Delta\theta$ 는 0.641$^\circ$ 로 수렴 — 0 이 되지 않는다",
            xy=(1.2e5, 0.6412), xytext=(0.16, 0.135), textcoords="axes fraction",
            fontsize=8.0, color=ORANGE, ha="left",
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9))
bx.axvspan(hs[0] / 1.6, hmin, color=GRAY, alpha=0.32, lw=0)
bx.annotate("700 km 는 이 영역 —\n36.0N 이 지평선 아래\n(최소 고도 1 504 km)",
            xy=(hmin * 0.72, 55), xytext=(0.035, 0.52), textcoords="axes fraction",
            fontsize=7.8, color=INK2, ha="left")
for h_sat, lab in ((35786.0, "GEO"), (5000.0, "5 000 km"), (2000.0, "2 000 km")):
    tt = th_surf(GAMMA, h_sat)
    bx.plot([h_sat, h_sat], [0, tt], color=MUTED, lw=0.8, dashes=(2, 3), zorder=0)
    bx.plot(h_sat, tt, "o", color=BLUE, ms=6, zorder=5)
    bx.plot(h_sat, tt - th_toa(tt), "o", color=ORANGE, ms=6, zorder=5)
    bx.annotate(lab, xy=(h_sat, tt), xytext=(h_sat * 1.06, tt + 3.0), fontsize=7.4,
                color=INK2, ha="left")
bx.set_xscale("log"); bx.set_xlim(hs[0] / 1.6, 3.0e5); bx.set_ylim(0, 95)
bx.set_xlabel("센서 고도 $h_{\\rm sat}$ (km, log)", fontsize=9.5)
bx.set_ylabel("angle (deg)", fontsize=9.5)
bx.tick_params(labelsize=8.4)
for sp in ("top", "right"): bx.spines[sp].set_visible(False)
bx.legend(fontsize=8.2, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.98))
bx.set_title(r"b) 표적·부성점 고정($\gamma$=36$^\circ$) 시에는 고도가 senz 를 바꾼다",
             fontsize=10.0, color=INK, fontweight="bold", loc="left", pad=10)
bx.text(0.0, -0.175, "센서를 올리면 senz 가 작아지고, 그래서 $\\Delta\\theta$ 도 9.18$^\\circ$(1 600 km) →\n"
        "0.79$^\\circ$(GEO) → 0.64$^\\circ$(무한)로 줄어든다. 직관이 맞는 부분은 이 경로다.",
        transform=bx.transAxes, fontsize=7.9, color=INK2, va="top", linespacing=1.5)

fig.suptitle("표적 36.0$^\\circ$N 128.25$^\\circ$E, 센서 0.0$^\\circ$N 128.25$^\\circ$E — "
             "표층 VZA 와 TOA VZA (구면 $R_e$=6371 km, 모델 TOA H=100 km)",
             fontsize=11.0, y=0.955, color=INK)
fig.savefig("/root/ocrt/ipss_bench/figI5_geo_vza_case.png", dpi=160,
            facecolor=SURF, bbox_inches="tight")
print("saved figI5_geo_vza_case.png")
