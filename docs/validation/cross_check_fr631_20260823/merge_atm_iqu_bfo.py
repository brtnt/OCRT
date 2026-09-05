#!/usr/bin/env python3
"""BFO Rayleigh IQU 비교 CSV 병합 (저장본만 사용, 신규 RT 실행 없음). 2026-08-25.

소스:
  OCRT : 05_VALIDATION/runs/atm_ocrt_bfo/fullgrid.csv  (rho_I/Q/U; VZA 0-60x5, RAA 0-345x15)
  OSOAA: 05_VALIDATION/runs/atm_osoaa_bfo_p{0,60,90,120,180}/Advanced_outputs/RESLUM_Adv_UP.txt
         level 0(TOA) I,Q,U (pi*L/Esun) -> /cos(SZA) 로 rho 환산, 부호 VZA축 선형보간
검증 게이트:
  G1: osoaa_rho_I 가 기존 atm_rhoI_bfo.csv 와 일치 (동일 보간 프로토콜 재현; rel<=2e-4
      허용 - 기존 CSV는 Standard REFL 열, 본 병합은 Adv IQU/cos(SZA): 소스 파일 유효자리 차)
  G2: ocrt_rho_I 가 기존 CSV와 완전 일치 (거울 raa 노드 대칭 포함; rel<=1e-9)
  G3: U 분기 판정 - 거울쌍 {180-phi, 180+phi} 중 U 잔차가 작은 쪽 채택, 반대쪽 기각 통계 출력
  G4: 주평면(플레인 0/180) U ~ 0 (양 코드)
"""
import csv, math, os
import numpy as np

V   = "/root/ocrt/runtime/MIGRATION_PKG_2026-08-19/05_VALIDATION"
OUT = "/root/ocrt/xcheck_fig/atm_rhoIQU_bfo.csv"
SZA = 40.0
MU0 = math.cos(math.radians(SZA))
PLANES = (0, 60, 90, 120, 180)
VZAS = np.arange(5.0, 61.0, 5.0)

# ---------- OCRT fullgrid ----------
O = {}
for r in csv.DictReader(open(f"{V}/runs/atm_ocrt_bfo/fullgrid.csv")):
    O[(float(r["vza_deg"]), float(r["raa_deg"]))] = (
        float(r["rho_I"]), float(r["rho_Q"]), float(r["rho_U"]))

# ---------- OSOAA Adv_UP level-0 파서 ----------
def read_adv_toa(plane):
    p = f"{V}/runs/atm_osoaa_bfo_p{plane}/Advanced_outputs/RESLUM_Adv_UP.txt"
    vza, I, Q, U = [], [], [], []
    for ln in open(p):
        s = ln.split()
        if len(s) == 10:
            try:
                lvl = int(s[0])
            except ValueError:
                continue
            if lvl != 0:
                continue
            vza.append(float(s[2])); I.append(float(s[4]))
            Q.append(float(s[5]));  U.append(float(s[6]))
    a = np.argsort(np.array(vza))
    return (np.array(vza)[a], np.array(I)[a], np.array(Q)[a], np.array(U)[a])

ADV = {p: read_adv_toa(p) for p in PLANES}

def osoaa_rho(plane, vza_signed):
    v, I, Q, U = ADV[plane]
    return tuple(np.interp(vza_signed, v, x) / MU0 for x in (I, Q, U))

# ---------- 기존 I-only CSV (게이트 기준) ----------
LEG = {}
for r in csv.DictReader(open("/root/ocrt/xcheck_fig/atm_rhoI_bfo.csv")):
    LEG[(int(r["plane_phi"]), float(r["vza_signed"]))] = (
        float(r["raa_deg"]), float(r["ocrt_rho_I"]), float(r["osoaa_rho_I"]))

# ---------- 행 구조 (기존 84행과 동일한 plane/side 집합) ----------
rows_spec = []           # (plane, side) ; side=+1 -> phi_view=plane, -1 -> plane+180
for p in PLANES:
    rows_spec.append((p, +1))
for p in (0, 180):       # 주평면은 양쪽 분기 모두 (기존 84행 구조)
    rows_spec.append((p, -1))

# ---------- U 분기 판정: raa 후보 A=(180-phi)%360, B=(phi+180)%360 ----------
def build(rows_spec, pick):   # pick: 'A' or 'B'
    out = []
    for plane, side in rows_spec:
        phi_view = plane if side > 0 else (plane + 180) % 360
        raaA = (180.0 - phi_view) % 360.0
        raaB = (phi_view + 180.0) % 360.0
        raa = raaA if pick == "A" else raaB
        for vz in VZAS:
            oc = O[(vz, raa)]
            os_ = osoaa_rho(plane, side * vz)
            out.append(dict(plane=plane, side=side, phi_view=phi_view,
                            vza_signed=side * vz, raa=raa, oc=oc, os=os_))
    return out

resA, resB = build(rows_spec, "A"), build(rows_spec, "B")
maxI = max(abs(r["os"][0]) for r in resA)

def ustat(res):
    d = [100 * (r["oc"][2] - r["os"][2]) / maxI for r in res
         if r["plane"] in (60, 90, 120)]
    d = np.array(d)
    return float(np.abs(d).max()), float(np.median(d))

uA, uB = ustat(resA), ustat(resB)
pick, res, ustats, urej = ("A", resA, uA, uB) if uA[0] < uB[0] else ("B", resB, uB, uA)
print(f"[G3] U-branch: A(raa=180-phi) max|dev| {uA[0]:.3f} | B(raa=phi+180) max|dev| {uB[0]:.3f}"
      f"  -> 채택 {pick} (기각측 {urej[0]:.3f})")

# ---------- 게이트 G1/G2 ----------
g1 = g2 = 0.0
for r in res:
    key = (r["plane"], r["vza_signed"])
    if key in LEG:
        _, leg_oc, leg_os = LEG[key]
        g2 = max(g2, abs(r["oc"][0] / leg_oc - 1))
        g1 = max(g1, abs(r["os"][0] / leg_os - 1))
n_leg = sum(1 for r in res if (r["plane"], r["vza_signed"]) in LEG)
print(f"[G1] osoaa_rho_I vs 기존CSV: max rel {g1:.2e} (n={n_leg}/84) -> {'PASS' if g1 <= 2e-4 else 'FAIL'}")
print(f"[G2] ocrt_rho_I  vs 기존CSV: max rel {g2:.2e} -> {'PASS' if g2 <= 5e-9 else 'FAIL'}"
      "  (기존 CSV 9유효자리 인쇄 반올림 한계 5e-9; 거울 raa 노드 대칭 동시 검증)")

# ---------- 게이트 G4: 주평면 U ----------
u_pp = max(max(abs(r["oc"][2]), abs(r["os"][2])) for r in res if r["plane"] in (0, 180))
print(f"[G4] 주평면 |U| max {u_pp:.2e} (rho 단위) -> {'PASS' if u_pp < 1e-6 else 'FAIL'}")

# ---------- 통계 ----------
for comp, idx in (("I", 0), ("Q", 1), ("U", 2)):
    if comp == "I":
        d = np.array([100 * (r["oc"][0] / r["os"][0] - 1) for r in res])
        print(f"rho_I : median {np.median(d):+.3f}%  max|.| {np.abs(d).max():.3f}%")
    else:
        d = np.array([100 * (r["oc"][idx] - r["os"][idx]) / maxI for r in res])
        print(f"rho_{comp} : median {np.median(d):+.3f}  max|.| {np.abs(d).max():.3f}  (%/max|I|, maxI={maxI:.6f})")

# ---------- CSV ----------
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["plane_phi", "osoaa_phi_view", "vza_signed", "ocrt_raa_deg",
                "ocrt_rho_I", "ocrt_rho_Q", "ocrt_rho_U",
                "osoaa_rho_I", "osoaa_rho_Q", "osoaa_rho_U",
                "dev_I_pct", "dev_Q_pctmaxI", "dev_U_pctmaxI"])
    for r in sorted(res, key=lambda x: (x["plane"], -x["side"], abs(x["vza_signed"]))):
        oc, os_ = r["oc"], r["os"]
        w.writerow([r["plane"], r["phi_view"], f'{r["vza_signed"]:.1f}', f'{r["raa"]:.1f}',
                    *[f"{x:.8e}" for x in oc], *[f"{x:.8e}" for x in os_],
                    f"{100*(oc[0]/os_[0]-1):.4f}",
                    f"{100*(oc[1]-os_[1])/maxI:.4f}",
                    f"{100*(oc[2]-os_[2])/maxI:.4f}"])
print(f"WROTE {OUT} rows={len(res)}")
