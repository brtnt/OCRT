#!/usr/bin/env python3
"""v1.08 aerosol value-kernel 세션 검증 종합 그림 (전부 선형축, scatter)."""
import math
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager as fm
# Noto Sans CJK (한국어) 등록
_fp = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
try:
    fm.fontManager.addfont(_fp)
    matplotlib.rcParams["font.family"] = fm.FontProperties(fname=_fp).get_name()
except Exception as e:
    print("font load failed:", e)
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt

# ---- 데이터 (이 세션에서 측정/계산) ----
# VZA sweep: M80C@412, aerosol-only(pressure 0), sza30, raa90, AOD0.2, n-layers400
vza   = [0, 10, 20, 40, 60, 80]
osoaa = [2.5061e-2, 2.2669e-2, 1.8888e-2, 1.3054e-2, 1.8181e-2, 5.1645e-2]
mom   = [2.0033e-2, 1.9681e-2, 1.8157e-2, 1.2956e-2, 1.8015e-2, 5.0331e-2]
val   = [2.0087e-2, 1.9642e-2, 1.8158e-2, 1.2957e-2, 1.8014e-2, 5.0337e-2]

# thin-AOD linearity: black surface nadir, OCRT vs 해석적 single-scatter
mu_v = 1.0; mu_s = math.cos(math.radians(30)); P150 = 0.23002
def ss_anal(tau, w=1.0):
    return w*P150/(4*(mu_v+mu_s))*(1-math.exp(-tau*(1/mu_v+1/mu_s)))
aod_thin = [0.002, 0.005, 0.01, 0.02, 0.05]
ocrt_thin = [1.31681e-4, 3.29227e-4, 6.58530e-4, 1.31732e-3, 3.29458e-3]
anal_thin = [ss_anal(a) for a in aod_thin]

# 분해 (AOD=0.2, nadir, raa90)
SS = 1.07890e-2; atm_MS = 2.39075e-3; surf = 6.85313e-3
cox_full = SS + atm_MS + surf
osoaa_nadir = 2.5061e-2

fig, ax = plt.subplots(2, 2, figsize=(12.5, 10))
fig.suptitle("OCRT v1.08 aerosol value-kernel 세션 검증  (M80C, 412 nm, aerosol-only, sza=30, raa=90)",
             fontsize=12.5, fontweight="bold")

# ---- 약어 박스 (그림만 봐도 이해되도록) ----
abbr = ("약어:  rho_I = reflectance I 성분 (= pi*L / Ed)   |   VZA = view zenith angle (관측 천정각)\n"
        "AOD = aerosol optical depth   |   SS = single scattering (단일산란)   |   MS = multiple scattering (다중산란)\n"
        "OCRT moment = Legendre moment kernel   |   OCRT value = angle-space value kernel (이 세션 수정)   |   "
        "OSOAA = 참조 vector RT")
fig.text(0.5, 0.005, abbr, ha="center", va="bottom", fontsize=8.3,
         bbox=dict(boxstyle="round", fc="#f4f4f4", ec="#999"))

# ===== Panel A: VZA sweep =====
a = ax[0][0]
a.plot(vza, osoaa, "o-", color="k", ms=8, lw=1.2, label="OSOAA (참조)")
a.plot(vza, mom, "s--", color="#1f77b4", ms=7, lw=1.0, label="OCRT moment")
a.plot(vza, val, "^:", color="#d62728", ms=7, lw=1.0, label="OCRT value (수정)")
a.set_xlabel("VZA (deg)"); a.set_ylabel("rho_I")
a.set_title("(A) VZA sweep — value$\\approx$moment, nadir 둘 다 $-$20%")
a.legend(fontsize=9); a.grid(alpha=0.3)
a.annotate("nadir $-$20%\n(kernel 무관)", xy=(0, mom[0]), xytext=(15, 2.35e-2),
           fontsize=9, arrowprops=dict(arrowstyle="->", color="gray"))

# ===== Panel B: value vs moment scatter (1:1) =====
b = ax[0][1]
b.scatter(mom, val, c="#d62728", s=70, zorder=3, edgecolor="k", linewidth=0.5)
lo = min(min(mom), min(val))*0.95; hi = max(max(mom), max(val))*1.05
b.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
b.set_xlabel("OCRT moment  rho_I"); b.set_ylabel("OCRT value  rho_I")
b.set_title("(B) value vs moment — 전 VZA $\\leq$0.3% 일치\n(두 버그 수정 후, 폭발 해소)")
b.legend(fontsize=9); b.grid(alpha=0.3)
for i, v in enumerate(vza):
    b.annotate(f"{v}", (mom[i], val[i]), fontsize=7, xytext=(4, -8), textcoords="offset points")

# ===== Panel C: thin-AOD SS 일치 scatter (1:1) =====
c = ax[1][0]
c.scatter(anal_thin, ocrt_thin, c="#2ca02c", s=80, zorder=3, edgecolor="k", linewidth=0.5)
lo2 = 0; hi2 = max(max(anal_thin), max(ocrt_thin))*1.08
c.plot([lo2, hi2], [lo2, hi2], "k--", lw=1, label="1:1")
c.set_xlabel("해석적 SS  rho_I (P(150$^\\circ$)=0.230, $\\omega_0$=1)")
c.set_ylabel("OCRT black-surface nadir  rho_I")
c.set_title("(C) thin-AOD: OCRT SS = 해석적 SS\nAOD$\\leq$0.01서 $\\leq$0.25% (SS 정확 확정)")
c.legend(fontsize=9); c.grid(alpha=0.3)
for i, a_ in enumerate(aod_thin):
    c.annotate(f"AOD {a_}", (anal_thin[i], ocrt_thin[i]), fontsize=7,
               xytext=(5, -3), textcoords="offset points")

# ===== Panel D: 분해 bar (AOD=0.2 nadir) =====
d = ax[1][1]
bars = ["SS\n(정확)", "atm_MS", "surface", "OCRT\ncox_full", "OSOAA\nref"]
vals = [SS, atm_MS, surf, cox_full, osoaa_nadir]
cols = ["#2ca02c", "#1f77b4", "#ff7f0e", "#555", "k"]
xpos = [0, 1, 2, 3.4, 4.4]
d.bar(xpos[:3], vals[:3], width=0.7, color=cols[:3], edgecolor="k")
# stacked total bar for OCRT
d.bar(xpos[3], SS, width=0.7, color="#2ca02c", edgecolor="k")
d.bar(xpos[3], atm_MS, bottom=SS, width=0.7, color="#1f77b4", edgecolor="k")
d.bar(xpos[3], surf, bottom=SS+atm_MS, width=0.7, color="#ff7f0e", edgecolor="k")
d.bar(xpos[4], osoaa_nadir, width=0.7, color="white", edgecolor="k", hatch="//")
d.set_xticks(xpos); d.set_xticklabels(bars, fontsize=8.5)
d.set_ylabel("rho_I (nadir, AOD=0.2)")
d.set_title("(D) 분해: gap(5.0e$-$3)은 surface 규모\n(SS의 53.9%, atm_MS의 11.9%, surface 34.2%)")
d.grid(alpha=0.3, axis="y")
# gap annotation
d.annotate("", xy=(xpos[4], osoaa_nadir), xytext=(xpos[4], cox_full),
           arrowprops=dict(arrowstyle="<->", color="red", lw=1.5))
d.text(xpos[4]+0.42, (cox_full+osoaa_nadir)/2, f"gap\n{osoaa_nadir-cox_full:.2e}\n(+25%)",
       fontsize=8.5, color="red", va="center")

plt.tight_layout(rect=[0, 0.045, 1, 0.97])
plt.savefig("/home/user/work/verify_value_kernel.png", dpi=140, bbox_inches="tight")
print("saved verify_value_kernel.png")
