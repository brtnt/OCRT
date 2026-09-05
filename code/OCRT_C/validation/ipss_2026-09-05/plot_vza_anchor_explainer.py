#!/usr/bin/env python3
"""VZA 앵커 설명 그림 (2026-09-05).

(a) 직선 광선의 국소 천정각이 고도에 따라 변하는 이유 — 모식도(곡률 과장).
(b) 센서 거리는 무관하다 — 광선은 하나의 같은 직선, 변하는 건 국소 연직뿐.
(c) 평면평행 코드에는 이 문제가 없다 (OSOAA 등: mu = 상수).
(d) 그래서 앵커는 화소여야 한다 — 소멸계수가 몰린 고도에서의 각 불일치.

프로젝트 표준 팔레트(figS1~S6, figI1~I3 동일).  dataviz 검증:
  #2a78d6 / #eb6834 — CVD ΔE 24.7 (protan), normal ΔE 33.6 — 전 항목 PASS.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import Arc

_KO = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
fm.fontManager.addfont(_KO)
matplotlib.rcParams["font.family"] = [fm.FontProperties(fname=_KO).get_name(), "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

BLUE, ORANGE, GRAY, SURF = "#2a78d6", "#eb6834", "#b9b7ae", "#fcfcfb"
INK, INK2, MUTED = "#1c1c1a", "#4a4a46", "#8a8a84"
D = "/root/ocrt/ipss_bench"
RE, HATM, TH_G = 6371.0, 100.0, 55.0


def theta_at(h, p):
    return np.degrees(np.arcsin(np.clip(p / (RE + h), -1.0, 1.0)))


P_SURF = RE * np.sin(np.radians(TH_G))
P_TOA = (RE + HATM) * np.sin(np.radians(TH_G))

fig = plt.figure(figsize=(13.4, 9.6), facecolor=SURF)
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.16,
                      left=0.045, right=0.985, top=0.855, bottom=0.085)

RES, HS = 260.0, 100.0
tv = np.radians(TH_G)
u = np.array([np.sin(tv), np.cos(tv)])
Pp = np.array([0.0, RES])
s_toa = -RES * np.cos(tv) + np.sqrt((RES * np.cos(tv)) ** 2 + (RES + HS) ** 2 - RES ** 2)
Q = Pp + s_toa * u


def head(ax, txt):
    ax.set_title(txt, fontsize=9.6, color=INK, fontweight="bold",
                 loc="left", pad=10, linespacing=1.45)


def foot(ax, txt, y=-0.035):
    ax.text(0.0, y, txt, transform=ax.transAxes, fontsize=7.9,
            color=INK2, va="top", ha="left", linespacing=1.5)


def shells(ax, span=(58, 122)):
    a = np.linspace(np.radians(span[0]), np.radians(span[1]), 400)
    ax.fill_between(np.concatenate([RES * np.cos(a), (RES + HS) * np.cos(a[::-1])]),
                    np.concatenate([RES * np.sin(a), (RES + HS) * np.sin(a[::-1])]),
                    color=GRAY, alpha=0.18, lw=0, zorder=1)
    ax.plot(RES * np.cos(a), RES * np.sin(a), "-", color=INK, lw=2.0, zorder=3)
    ax.plot((RES + HS) * np.cos(a), (RES + HS) * np.sin(a), "--", color=MUTED,
            lw=1.4, zorder=2)


def vertical(ax, pt, length, dashes=(4, 3)):
    n = pt / np.linalg.norm(pt)
    ax.plot([pt[0], pt[0] + length * n[0]], [pt[1], pt[1] + length * n[1]],
            color=ORANGE, lw=1.3, dashes=dashes, zorder=4)
    return n


def arc(ax, pt, n, rad, text="", fs=9.5):
    a0 = np.degrees(np.arctan2(n[1], n[0]))
    a1 = np.degrees(np.arctan2(u[1], u[0]))
    ax.add_patch(Arc(pt, 2 * rad, 2 * rad, theta1=min(a0, a1), theta2=max(a0, a1),
                     color=INK, lw=1.4, zorder=5))
    if text:
        m = np.radians(0.5 * (a0 + a1))
        ax.annotate(text, xy=(pt[0] + 1.42 * rad * np.cos(m),
                              pt[1] + 1.42 * rad * np.sin(m)),
                    fontsize=fs, color=INK, ha="center", va="center",
                    fontweight="bold")


# ------------------------------------------------------------------ (a)
axA = fig.add_subplot(gs[0, 0]); axA.set_facecolor(SURF)
shells(axA)
end = Pp + 1.62 * s_toa * u
axA.annotate("sensor: anywhere above H\n(does not enter the model)", xy=end,
             xytext=(end[0] - 118, end[1] + 44), fontsize=7.4, color=MUTED,
             ha="left", arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
axA.annotate("", xy=end, xytext=Pp,
             arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.6,
                             shrinkA=0, shrinkB=0), zorder=6)
nP = vertical(axA, Pp, 108)
nQ = vertical(axA, Q, 108)
arc(axA, Pp, nP, 50, r"$\theta_{\rm pixel}$")
arc(axA, Q, nQ, 50, r"$\theta_{\rm TOA}$")
axA.plot(*Pp, "o", color=BLUE, ms=8, zorder=7)
axA.plot(*Q, "o", color=MUTED, ms=6, zorder=7)
axA.annotate("target pixel P", xy=Pp, xytext=(30, RES - 62), fontsize=8.4,
             color=BLUE, fontweight="bold",
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=0.9))
axA.annotate("local vertical\nat the pixel", xy=(0, RES + 100),
             xytext=(-124, RES + 146), fontsize=7.6, color=ORANGE,
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.8))
axA.annotate("local vertical\nat 100 km", xy=Q + 108 * nQ,
             xytext=(Q[0] + 118, Q[1] + 60), fontsize=7.6, color=ORANGE)
axA.annotate("TOA, 100 km", xy=(-100, RES + HS - 20), xytext=(-124, RES + 60),
             fontsize=7.6, color=MUTED)
axA.annotate("ground", xy=(-95, RES - 18), xytext=(-124, RES - 62),
             fontsize=7.6, color=INK2)
axA.set_xlim(-125, 360); axA.set_ylim(180, 480)
axA.set_aspect("equal", adjustable="datalim"); axA.axis("off")
head(axA, "a) 직선 광선인데도 국소 천정각은 고도에 따라 변한다\n"
          r"     직선의 불변량 $r\,\sin\theta(r)=$const  $\Rightarrow$  위로 갈수록 $\theta$ 가 작아진다")
foot(axA, "모식도 — 곡률 과장(그림은 H/R$_e$=1/2.6, 실제는 100/6371=1/63.7).\n"
          "실제 값:  $\\theta_{\\rm pixel}$ = 55.000$^\\circ$  →  $\\theta_{\\rm TOA}$ = 53.755$^\\circ$   (차 1.245$^\\circ$).\n"
          "모델에 등장하는 고도는 이 둘뿐이다. 센서는 h > H 어디에 있어도 무관하다 —\n"
          "H 위에는 소광이 없어 복사휘도가 변하지 않고, 모델은 h = H 에서 값을 낸다.")

# ------------------------------------------------------------------ (b)
axB = fig.add_subplot(gs[0, 1]); axB.set_facecolor(SURF)
hs = np.genfromtxt(f"{D}/hscan_sza0.csv", delimiter=",", names=True)
axB.axhline(0, color=INK2, lw=1.2, zorder=1)
axB.plot(hs["H_km"], hs["dev_surf_pct"], "-", color=BLUE, lw=2.4,
         label="surface (pixel) anchor — v1.11")
axB.plot(hs["H_km"], hs["dev_toa_pct"], "-", color=ORANGE, lw=2.4,
         label="TOA anchor — 초판 결함")
axB.axvspan(20, 30, color=GRAY, alpha=0.22, lw=0)
axB.annotate("소광이 실제로\n있는 구간", xy=(28, 1.6), xytext=(38, 2.6),
             fontsize=7.4, color=INK2, ha="left",
             arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
for H, dth in ((100, -1.25), (180, -2.20)):
    axB.annotate(f"H={H}: $\\Delta\\theta$={dth:.2f}$^\\circ$",
                 xy=(H, float(np.interp(H, hs["H_km"], hs["dev_toa_pct"]))),
                 xytext=(H + 12, float(np.interp(H, hs["H_km"], hs["dev_toa_pct"])) - 0.70),
                 fontsize=7.2, color=ORANGE, ha="left",
                 arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.7))
axB.annotate("H 를 바꿔도 답이 거의 안 변한다\n(30 km 위는 비어 있으므로 — 물리적으로 옳다)",
             xy=(150, float(np.interp(150, hs["H_km"], hs["dev_surf_pct"]))),
             xytext=(0.30, 0.16), textcoords="axes fraction", fontsize=7.8,
             color=BLUE, fontweight="bold", ha="left",
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.0))
axB.set_xlim(18, 208); axB.set_ylim(-0.9, 7.2)
axB.set_xlabel("모델 TOA 높이 H (km)", fontsize=9.2)
axB.set_ylabel(r"(plane-parallel / spherical $-$ 1)  (%)", fontsize=9.2)
axB.tick_params(labelsize=8.2)
for sp in ("top", "right"):
    axB.spines[sp].set_visible(False)
axB.legend(fontsize=7.8, frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.99))
head(axB, "b) 차이를 만드는 것은 모델 TOA 높이 H 이고, 센서 높이가 아니다\n"
          "     VZA 55$^\\circ$(화소 기준), SZA 0, 지수 8 km 프로파일, "
          "$\\tau$=0.0935")
foot(axB, "$\\Delta\\theta=\\theta_{\\rm TOA}-\\theta_{\\rm pixel}$ 은 H 와 함께 커진다 — 여기까지는 맞다. 그런데 지상 앵커는\n"
          "H=60 km 부터 값이 고정되고(−0.24 %), TOA 앵커는 H 에 비례해 계속 커진다(+0.7→+2.8→+6.3 %).\n"
          "모델을 어디서 잘랐는지가 답을 바꾸면 그 앵커는 틀렸다. 센서 고도는 어느 쪽에도 들어가지 않는다.",
          y=-0.155)

# ------------------------------------------------------------------ (c)
axC = fig.add_subplot(gs[1, 0]); axC.set_facecolor(SURF)
axC.fill_between([-25, 250], 0, 100, color=GRAY, alpha=0.18, lw=0)
axC.plot([-25, 250], [0, 0], color=INK, lw=2.0)
axC.plot([-25, 250], [100, 100], "--", color=MUTED, lw=1.4)
xe = 100.0 * np.tan(np.radians(TH_G))
axC.annotate("", xy=(1.30 * xe, 130), xytext=(0, 0),
             arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.6,
                             shrinkA=0, shrinkB=0), zorder=5)
for hh in (0, 25, 50, 75, 100):
    xh = hh * np.tan(np.radians(TH_G))
    axC.plot([xh, xh], [hh, hh + 30], color=ORANGE, lw=1.3, dashes=(4, 3), zorder=4)
    axC.add_patch(Arc((xh, hh), 30, 30, theta1=55.0, theta2=90.0,
                      color=INK, lw=1.2, zorder=5))
    axC.annotate("55$^\\circ$", xy=(xh + 19, hh + 12), fontsize=7.8, color=INK)
axC.plot(0, 0, "o", color=BLUE, ms=8, zorder=6)
axC.annotate("TOA", xy=(-20, 104), fontsize=7.6, color=MUTED)
axC.annotate("ground", xy=(-20, -12), fontsize=7.6, color=INK2)
axC.set_xlim(-25, 250); axC.set_ylim(-24, 150)
axC.set_aspect("equal", adjustable="datalim"); axC.axis("off")
head(axC, "c) 평면평행 코드에는 앵커 문제가 아예 없다\n"
          r"     국소 연직이 모든 고도에서 평행  $\Rightarrow$  $\theta$ 는 고도 무관 상수")
foot(axC, "확인: OSOAA 소스 전체에 spherical / Chapman / Earth radius 가 없다.\n"
          "투과항은 그냥  DEXP(-DTAU/RMU(K))  —  RMU 는 상수 방향코사인.\n"
          "→ \"RT 코드에서 VZA·SZA 는 TOA 와 표면이 같다\" 는 이해는 평면평행 코드에서 맞다.")

# ------------------------------------------------------------------ (d)
gsD = gs[1, 1].subgridspec(1, 2, width_ratios=[0.30, 0.70], wspace=0.06)
h = np.linspace(0, 40, 500)
beta = np.exp(-h / 8.0)

axD0 = fig.add_subplot(gsD[0, 0]); axD0.set_facecolor(SURF)
axD0.fill_betweenx(h, 0, beta, color=GRAY, alpha=0.45, lw=0)
axD0.plot(beta, h, "-", color=INK2, lw=1.2)
for hh, lab in ((10.0, "71 %"), (20.0, "92 %")):
    axD0.axhline(hh, color=MUTED, lw=0.9, dashes=(3, 3))
    axD0.annotate(f"{lab} of mass\nbelow {hh:.0f} km", xy=(0.96, hh + 1.0),
                  fontsize=7.2, color=INK2, ha="left", va="bottom")
axD0.set_xlim(0, 1.06); axD0.set_ylim(0, 40)
axD0.invert_xaxis()
axD0.set_ylabel("altitude (km)", fontsize=9.2)
axD0.set_xlabel(r"$\beta(h)$  (norm.)", fontsize=8.6)
axD0.tick_params(labelsize=8.2); axD0.set_xticks([0, 0.5, 1.0])
for sp in ("top", "left"):
    axD0.spines[sp].set_visible(False)

axD = fig.add_subplot(gsD[0, 1], sharey=axD0); axD.set_facecolor(SURF)
axD.plot(theta_at(h, P_SURF) - TH_G, h, "-", color=BLUE, lw=2.4,
         label="surface (pixel) anchor — v1.11")
axD.plot(theta_at(h, P_TOA) - TH_G, h, "-", color=ORANGE, lw=2.4,
         label="TOA anchor — 초판 결함")
axD.axvline(0, color=INK2, lw=1.5, zorder=1)
axD.axhline(10.0, color=MUTED, lw=0.9, dashes=(3, 3), zorder=0)
axD.axhline(20.0, color=MUTED, lw=0.9, dashes=(3, 3), zorder=0)
axD.annotate("plane-parallel model\nuses 55$^\\circ$ at every h",
             xy=(0, 31.0), xytext=(0.04, 0.70), textcoords="axes fraction",
             fontsize=7.6, color=INK2, ha="left",
             arrowprops=dict(arrowstyle="->", color=INK2, lw=0.9))
axD.annotate("질량이 있는 곳에서 일치\n$\\Rightarrow$ 보정량 $\\to$ 0",
             xy=(-0.055, 2.5), xytext=(0.03, 0.115), textcoords="axes fraction",
             fontsize=8.0, color=BLUE, fontweight="bold", ha="left",
             arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.1,
                             connectionstyle="arc3,rad=-0.25"))
axD.annotate("질량이 있는 곳에서 1.31$^\\circ$ 어긋남\n"
             "$\\Rightarrow$ 경로 +3.4 % $\\Rightarrow$ $\\kappa$ 가 편향된다",
             xy=(1.29, 4.0), xytext=(0.24, 0.44), textcoords="axes fraction",
             fontsize=8.0, color=ORANGE, fontweight="bold", ha="left",
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.1))
axD.set_xlim(-0.72, 1.52)
axD.set_xlabel(r"local zenith angle $-$ 55$^\circ$   (deg)", fontsize=9.2)
axD.tick_params(labelsize=8.2, labelleft=False)
for sp in ("top", "right"):
    axD.spines[sp].set_visible(False)
axD.legend(fontsize=7.8, frameon=False, loc="upper left",
           bbox_to_anchor=(0.02, 1.0))
head(axD0, "d) 소멸계수는 하층에 몰려 있다 — 그래서 앵커가 화소여야 한다")
axD0.text(0.0, -0.135,
          "척도높이를 0 으로 줄이면 지상앵커 편차 $-$0.24 → $-$0.004 %,\n"
          "TOA 앵커는 +2.85 → +3.11 % 로 남는다 (게이트 G-I3).",
          transform=axD0.transAxes, fontsize=7.9, color=INK2, va="top",
          ha="left", linespacing=1.5)

fig.suptitle("VZA 앵커 — 평면평행 해에 구면 보정을 곱할 때 각을 어느 고도에 고정할 것인가"
             "   |   OCRT v1.11: 화소(지상) 기준 확정   |   예시 VZA 55$^\\circ$, RAA 90$^\\circ$",
             fontsize=11.0, y=0.965, color=INK)
fig.savefig(f"{D}/figI4_vza_anchor.png", dpi=160, facecolor=SURF, bbox_inches="tight")
print("saved figI4_vza_anchor.png")

for lab, p in (("surface anchor", P_SURF), ("TOA anchor", P_TOA)):
    print(f"{lab:16s} theta(0)={theta_at(0.0,p):8.4f}  theta(100)={theta_at(HATM,p):8.4f}")
for r, lab in ((RE + 100, "100 km"), (RE + 800, "800 km"), (42164.0, "GEO")):
    print(f"  sensor {lab:8s} local zenith = {np.degrees(np.arcsin(P_SURF/r)):7.3f} deg")
print(f"  mass below 10/20 km = {(1-np.exp(-10/8))*100:.1f} / {(1-np.exp(-20/8))*100:.1f} %")
print(f"  sec(56.306)/sec(55)-1 = {np.cos(np.radians(55))/np.cos(np.radians(56.3057))-1:+.4f}")
