#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""7월 tsm0.5/555 구성 회귀 파일럿: 현재 코드(2모드) 대 보관 OSOAA."""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/user/work/diag_out/OCRT_JULY_CONFIG_REGRESSION_PILOT_TSM555_2026-08-18"
for d in ("tables", "figures", "raw"):
    os.makedirs(f"{OUT}/{d}", exist_ok=True)

ARCH = ("/home/user/work/OCRT_MIE_COMPLETE_HANDOFF_COMPACT_2026-08-17/06_VALIDATION_REFERENCE/"
        "OCRT_OSOAA_STEP_A_RAW_MATCHUP_LOCATIONS_2026-08-16/OCRT_OSOAA_underwater_validation_2026-07-26.csv")

july = {}
for r in csv.DictReader(open(ARCH)):
    if r["campaign"] == "tsm" and r["note"] == "총부유물 0.5 g/m3" and r["band_nm"] == "555":
        key = (round(float(r["vza_deg"]), 3), round(float(r["raa_deg"]), 3))
        july[key] = {f"{s}_{lv}_{c}": float(r[f"{s}_{lv}_{c}"])
                     for s in ("OCRT", "OSOAA") for lv in ("rrs0minus", "Rrs0plus")
                     for c in ("I", "Q", "U")}
print(f"7월 보관 555 셀: {len(july)}")

def load_grid(path):
    out = {}
    for r in csv.DictReader(open(path)):
        raa = float(r["raa_deg"])
        if raa > 180.0 + 1e-9:
            continue
        key = (round(float(r["vza_deg"]), 3), round(raa, 3))
        out[key] = {"rrs0minus_I": float(r["rrs_I"]), "rrs0minus_Q": float(r["rrs_Q"]),
                    "rrs0minus_U": float(r["rrs_U"]), "Rrs0plus_I": float(r["Rrs_I"]),
                    "Rrs0plus_Q": float(r["Rrs_Q"]), "Rrs0plus_U": float(r["Rrs_U"])}
    return out

cur = {m: load_grid(f"/home/user/work/runs/regpilot_{m}/fullgrid.csv")
       for m in ("moment", "direct")}
keys = sorted(july.keys())
assert all(k in cur["moment"] and k in cur["direct"] for k in keys), "기하 매칭 실패"
print(f"기하 매칭: {len(keys)}/169")

def dolp(I, Q, U):
    return math.sqrt(Q * Q + U * U) / I if abs(I) > 0 else float("nan")

def metrics(name, get_test, get_ref, level):
    """관례: I 상대%, Q/U 차이/max|I_ref|%, DoLP pp."""
    maxI = max(abs(get_ref(k, "I")) for k in keys)
    out = []
    dI = [100 * (get_test(k, "I") - get_ref(k, "I")) / get_ref(k, "I") for k in keys]
    out.append([name, level, "I_relative_pct", len(keys),
                f"{np.mean(np.abs(dI)):.4f}", f"{np.max(np.abs(dI)):.4f}",
                f"{np.median([get_test(k,'I')/get_ref(k,'I') for k in keys]):.6f}"])
    for c in ("Q", "U"):
        d = [100 * (get_test(k, c) - get_ref(k, c)) / maxI for k in keys]
        out.append([name, level, f"{c}_diff_over_maxI_pct", len(keys),
                    f"{np.mean(np.abs(d)):.4f}", f"{np.max(np.abs(d)):.4f}", ""])
    dd = [100 * (dolp(*[get_test(k, c) for c in "IQU"]) -
                 dolp(*[get_ref(k, c) for c in "IQU"])) for k in keys]
    out.append([name, level, "DoLP_pp", len(keys),
                f"{np.mean(np.abs(dd)):.4f}", f"{np.max(np.abs(dd)):.4f}", ""])
    return out

rows = []
for level in ("rrs0minus", "Rrs0plus"):
    J_O = lambda k, c, lv=level: july[k][f"OCRT_{lv}_{c}"]
    J_S = lambda k, c, lv=level: july[k][f"OSOAA_{lv}_{c}"]
    M   = lambda k, c, lv=level: cur["moment"][k][f"{lv}_{c}"]
    D   = lambda k, c, lv=level: cur["direct"][k][f"{lv}_{c}"]
    rows += metrics("기준선: 7월OCRT vs 7월OSOAA", J_O, J_S, level)
    rows += metrics("현재(moment) vs 7월OSOAA",     M,   J_S, level)
    rows += metrics("현재(direct) vs 7월OSOAA",     D,   J_S, level)
    rows += metrics("현재(moment) vs 7월OCRT",      M,   J_O, level)
    rows += metrics("현재(direct) vs 7월OCRT",      D,   J_O, level)

with open(f"{OUT}/tables/regression_pilot_metrics.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["comparison", "level", "metric", "n", "mean_abs", "max_abs", "median_ratio"])
    w.writerows(rows)
for r in rows:
    if r[1] == "rrs0minus":
        print(f"{r[0]:<28} {r[2]:<22} mean|{r[4]}| max|{r[5]}| {('ratio '+r[6]) if r[6] else ''}")

# 원자료 보존
with open(f"{OUT}/raw/per_cell_values.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["vza", "raa"] +
               [f"{s}_{lv}_{c}" for s in ("julyOCRT", "julyOSOAA", "curMoment", "curDirect")
                for lv in ("rrs0minus", "Rrs0plus") for c in "IQU"])
    for k in keys:
        row = [k[0], k[1]]
        for src in (lambda lv, c: july[k][f"OCRT_{lv}_{c}"],
                    lambda lv, c: july[k][f"OSOAA_{lv}_{c}"],
                    lambda lv, c: cur["moment"][k][f"{lv}_{c}"],
                    lambda lv, c: cur["direct"][k][f"{lv}_{c}"]):
            row += [f"{src(lv, c):.10e}" for lv in ("rrs0minus", "Rrs0plus") for c in "IQU"]
        w.writerow(row)

# 그림: rrs 1:1 + 잔차 (I/Q/U). 파랑=현재(moment), 회색=7월OCRT 기준선.
fig, axes = plt.subplots(2, 3, figsize=(12.8, 7.0),
                         gridspec_kw={"height_ratios": [1.6, 1]})
maxI = max(abs(july[k]["OSOAA_rrs0minus_I"]) for k in keys)
for ci, c in enumerate("IQU"):
    ref = np.array([july[k][f"OSOAA_rrs0minus_{c}"] for k in keys])
    jul = np.array([july[k][f"OCRT_rrs0minus_{c}"] for k in keys])
    now = np.array([cur["moment"][k][f"rrs0minus_{c}"] for k in keys])
    ax = axes[0, ci]
    lo = min(ref.min(), now.min()); hi = max(ref.max(), now.max())
    pad = 0.06 * (hi - lo + 1e-12); lo -= pad; hi += pad
    ax.plot([lo, hi], [lo, hi], "k-", lw=0.8, label="1:1")
    for b in (0.2, 0.5):
        d = b / 100 * maxI
        ax.plot([lo, hi], [lo + d, hi + d], "--", color="0.6", lw=0.6)
        ax.plot([lo, hi], [lo - d, hi - d], "--", color="0.6", lw=0.6)
    ax.scatter(ref, jul, s=13, c="0.6", label="July OCRT (baseline)")
    ax.scatter(ref, now, s=10, c="#1f77b4", label="current (moment)")
    ax.set_title(f"rrs {c}", fontsize=10); ax.grid(alpha=0.25, lw=0.4)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    if ci == 0:
        ax.set_ylabel("OCRT rrs(0-)"); ax.legend(fontsize=7)
    axr = axes[1, ci]
    axr.scatter(ref, 100 * (jul - ref) / maxI, s=11, c="0.6")
    axr.scatter(ref, 100 * (now - ref) / maxI, s=9, c="#1f77b4")
    axr.axhline(0, color="k", lw=0.8)
    for b in (0.2, 0.5):
        axr.axhline(+b, ls="--", color="0.6", lw=0.6)
        axr.axhline(-b, ls="--", color="0.6", lw=0.6)
    axr.set_xlabel(f"July OSOAA rrs {c}")
    axr.grid(alpha=0.25, lw=0.4)
    if ci == 0:
        axr.set_ylabel("residual (% of max|I_OSOAA|)")
fig.suptitle("July-config regression pilot: Red_clay TSM 0.5, 555 nm, 169 geometries "
             "(bands ±0.2/±0.5% of max|I|)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{OUT}/figures/R1_rrs_current_vs_julyOSOAA.png", dpi=160)
plt.close(fig)
print("분석 완료 →", OUT)
