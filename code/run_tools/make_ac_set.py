#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
대기보정 검증용 무작위 표본 생성기 (GOCI-II 12밴드).

해수 성분, 관측 기하, 에어로졸 종류와 양을 무작위로 뽑아 대기보정 검증에
쓰는 항을 모두 만든다.

한 사례·한 밴드마다 다섯 번 실행한다. 코드가 이 항들을 한 번에 내주지
않기 때문이다.

  ① 해수 결합 + 에어로졸, 관측 방향   → ρ_TOA, Td_s, Td_v, Rrs(관측 방향)
  ② 흑색 프레넬 + 에어로졸            → ρ_(R+A)
  ③ 흑색 프레넬, 에어로졸 없음        → ρ_R
  ④ 해수 결합, 에어로졸 없음          → Td_R_s, Td_R_v
  ⑤ 해수 결합 + 에어로졸, 천저 방향   → Rrs(천저)

유도

  ρ_C   = ① − ③                       레일리 보정 반사도
  ρ_Am  = ② − ③                       에어로졸 + 레일리·에어로졸 상호작용
  Td_s  = ①의 하향 전투과율            해수면 ρ_w 를 TOA 로 옮기는 하향 몫
  Td_v  = ①의 상향 전투과율            같은 것의 상향 몫
  Td_R_* = ④에서 얻은 레일리만의 값
  Td_Am_* = Td_* ÷ Td_R_*              곱셈 분해

Td 의 정의가 코드 출력과 맞는지는 검산으로 확인했다.
수면 위 ρ_w = π·Rrs 이고, TOA 의 ρ_w 를 그것으로 나눈 값이 Td_s×Td_v 와
상대차 1e-7 로 같다. 결과표에 이 검산 값을 열로 함께 남긴다.

쓰는 법
  python make_ac_set.py --dry-run          계획과 예상 시간만 본다
  python make_ac_set.py                    실행 (사례 1000개, 24코어)
  python make_ac_set.py --n 50             적은 수로 먼저 시험
"""

import argparse
import csv
import glob
import math
import os
import random
import subprocess
import sys
import time
from datetime import datetime, timedelta

if os.name == "nt":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
        ctypes.windll.kernel32.SetConsoleCP(65001)
    except Exception:
        pass
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


# ===========================================================================
# 1. 설계 상수 — 여기만 고치면 표본의 성격이 바뀐다
# ===========================================================================

# GOCI-II 12밴드
BANDS_NM = [380, 412, 443, 490, 510, 555, 620, 660, 680, 709, 745, 865]

SZA_RANGE = (0.0, 60.0)
VZA_RANGE = (0.0, 50.0)
RAA_RANGE = (0.0, 180.0)

CHL_RANGE = (0.1, 5.0)                     # mg/m^3
TSM_RANGE = (0.1, 20.0)                    # g/m^3
CDOM_RANGE = (0.01, 0.1)                   # 1/m at 440 nm

AOD_RANGE = (0.05, 0.30)                   # 865 nm
WIND_RANGE = (2.0, 10.0)                   # m/s

AER_DIR = os.path.join("inputs", "aerosol_ahmad2010_accurt_mie")

LOG_SAMPLE = True          # 해수 세 값은 로그 공간에서 균등하게 뽑는다
COVARY_SCATTER = 0.35      # --covary 를 켤 때의 흩뿌림 폭

# 한 사례·한 밴드에 드는 시간 (1코어 실측 기준)
NOMINAL_SEC = 21.0 + 3.6 + 0.5 + 11.0 + 14.0


# ===========================================================================
# 2. 설계표
# ===========================================================================

def loguni(rng, lo, hi):
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def glint_tilt_deg(sza, vza, raa):
    """직달 선글린트가 생기려면 수면이 기울어야 하는 각도.

    0 이면 잔잔한 수면에서 곧바로 선글린트가 생기는 방향이다.
    """
    s, v, r = math.radians(sza), math.radians(vza), math.radians(raa)
    sv = math.cos(s) * math.cos(v) + math.sin(s) * math.sin(v) * math.cos(r)
    den = math.sqrt(max(2.0 + 2.0 * sv, 1e-300))
    return math.degrees(math.acos(min(max((math.cos(s) + math.cos(v)) / den,
                                          -1.0), 1.0)))


def scattering_angle_deg(sza, vza, raa):
    s, v, r = math.radians(sza), math.radians(vza), math.radians(raa)
    c = -math.cos(s) * math.cos(v) - math.sin(s) * math.sin(v) * math.cos(r)
    return math.degrees(math.acos(min(max(c, -1.0), 1.0)))


def make_design(n, seed, aer_models, covary):
    rng = random.Random(seed)
    rows = []
    for i in range(1, n + 1):
        if covary:
            t = rng.random()
            def draw(lo, hi):
                base = math.log(lo) + t * (math.log(hi) - math.log(lo))
                return min(max(math.exp(base + rng.gauss(0.0, COVARY_SCATTER)),
                               lo), hi)
            chl, tsm, cdom = draw(*CHL_RANGE), draw(*TSM_RANGE), draw(*CDOM_RANGE)
        elif LOG_SAMPLE:
            chl = loguni(rng, *CHL_RANGE)
            tsm = loguni(rng, *TSM_RANGE)
            cdom = loguni(rng, *CDOM_RANGE)
        else:
            chl = rng.uniform(*CHL_RANGE)
            tsm = rng.uniform(*TSM_RANGE)
            cdom = rng.uniform(*CDOM_RANGE)
        sza = round(rng.uniform(*SZA_RANGE), 4)
        vza = round(rng.uniform(*VZA_RANGE), 4)
        raa = round(rng.uniform(*RAA_RANGE), 4)
        rows.append(dict(
            case_id=i, sza_deg=sza, vza_deg=vza, raa_deg=raa,
            scat_angle_deg=round(scattering_angle_deg(sza, vza, raa), 4),
            glint_tilt_deg=round(glint_tilt_deg(sza, vza, raa), 4),
            wind_ms=round(rng.uniform(*WIND_RANGE), 4),
            chl_mg_m3=round(chl, 6),
            tsm_g_m3=round(tsm, 6),
            acdom440_m_inv=round(cdom, 6),
            aer_model=rng.choice(aer_models),
            aod865=round(rng.uniform(*AOD_RANGE), 5)))
    return rows


# ===========================================================================
# 3. 실행 명령 다섯 가지
# ===========================================================================

def find_root(explicit):
    if explicit:
        return os.path.abspath(explicit)
    here = os.path.dirname(os.path.abspath(__file__))
    for c in [here, os.path.dirname(here),
              os.path.dirname(os.path.dirname(here)),
              os.getcwd(), os.path.dirname(os.getcwd())]:
        if c and os.path.isdir(os.path.join(c, "inputs")) \
              and os.path.isdir(os.path.join(c, "src")):
            return c
    return None


def find_exe(explicit, root):
    if explicit:
        return os.path.abspath(explicit)
    names = (["ocrt_v1.2.exe", "ocrt.exe", "ocrt_v1.2", "ocrt"] if os.name == "nt"
             else ["ocrt_v1.2", "ocrt", "ocrt_v1.2.exe", "ocrt.exe"])
    here = os.path.dirname(os.path.abspath(__file__))
    for d in [os.path.join(root, "build"), here, root]:
        for nm in names:
            p = os.path.join(d, nm)
            if os.path.isfile(p):
                return p
    return None


def base_geom(r, band, vza, raa):
    return ["--sza", "%g" % r["sza_deg"],
            "--vza", "%g" % vza,
            "--raa", "%g" % raa,
            "--wavelength", "%g" % band,
            "--wind-speed", "%g" % r["wind_ms"],
            "--decouple-sunglint",
            "--pssa"]


def water_opts(r):
    return ["--water-model", "ocrt",
            "--ocrt-chl", "%g" % r["chl_mg_m3"],
            "--ocrt-tsm", "%g" % r["tsm_g_m3"],
            "--ocrt-adom440", "%g" % r["acdom440_m_inv"]]


def aer_opts(r, root):
    return ["--aod-865", "%g" % r["aod865"],
            "--mie", os.path.join(AER_DIR, r["aer_model"] + ".mie")]


def build_cmds(exe, root, r, band):
    g = base_geom(r, band, r["vza_deg"], r["raa_deg"])
    gn = base_geom(r, band, 0.0, 0.0)          # 천저 방향
    return {
        "ocean_aer":   [exe] + g + ["--surface", "ocean"] + aer_opts(r, root)
                       + water_opts(r),
        "black_aer":   [exe] + g + ["--surface", "black_fresnel_ocean"]
                       + aer_opts(r, root),
        "black_ray":   [exe] + g + ["--surface", "black_fresnel_ocean"],
        "ocean_ray":   [exe] + g + ["--surface", "ocean"] + water_opts(r),
        "ocean_nadir": [exe] + gn + ["--surface", "ocean"] + aer_opts(r, root)
                       + water_opts(r),
    }


def parse_tail(text):
    out = {}
    lines = [l for l in text.splitlines() if l.strip()]
    if not lines:
        return out
    for tok in lines[-1].split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def run_one(args):
    exe, root, r, band, timeout = args
    rec = dict(case_id=r["case_id"], band_nm=band)
    t0 = time.time()
    cmds = build_cmds(exe, root, r, band)
    got = {}
    for key, cmd in cmds.items():
        try:
            p = subprocess.run(cmd, cwd=root, capture_output=True, text=True,
                               encoding="utf-8", errors="replace",
                               timeout=timeout)
        except Exception as e:
            rec["status"] = "failed"
            rec["note"] = "%s: %s" % (key, str(e)[:120])
            return rec
        if p.returncode != 0:
            rec["status"] = "failed"
            rec["note"] = "%s: %s" % (key, " ".join((p.stderr or "").split())[:180])
            return rec
        d = parse_tail(p.stdout)
        if "TOA_rho_I" not in d:
            rec["status"] = "failed"
            rec["note"] = "%s: 출력에서 값을 읽지 못함" % key
            return rec
        got[key] = d

    oa, ba, br, orr, on = (got["ocean_aer"], got["black_aer"], got["black_ray"],
                           got["ocean_ray"], got["ocean_nadir"])

    rho_TOA = oa["TOA_rho_I"]
    rho_R = br["TOA_rho_I"]
    rho_RA = ba["TOA_rho_I"]
    Td_s = oa["T_total_dn_hemi"]
    Td_v = oa["T_total_up_view"]
    Td_R_s = orr["T_total_dn_hemi"]
    Td_R_v = orr["T_total_up_view"]

    rec.update(
        Rho_TOA=rho_TOA,
        Rho_R=rho_R,
        Rho_C=rho_TOA - rho_R,
        Rho_Am=rho_RA - rho_R,
        Td_R_s=Td_R_s,
        Td_R_v=Td_R_v,
        Td_Am_s=Td_s / Td_R_s if Td_R_s else float("nan"),
        Td_Am_v=Td_v / Td_R_v if Td_R_v else float("nan"),
        Td_s=Td_s,
        Td_v=Td_v,
        Rrs=on["Rrs0plus_I"],                      # 천저 방향
        Rrs_view=oa["Rrs0plus_I"],                 # 관측 방향
        rrs0minus=on.get("rrs0minus_I"),
        t_rho_w=(rho_TOA - rho_RA),                # 차감으로 얻은 수출광 항
        AOD_band=oa.get("AOD_band"),
        orders=oa.get("orders"),
    )
    # 검산: 수면 위 rho_w 를 Td 로 옮긴 값이 TOA 수출광 항과 맞는가
    rho_w_surf = math.pi * oa["Rrs0plus_I"]
    if rho_w_surf > 0:
        rec["Td_check"] = (oa.get("TOA_rho_water_total_I", float("nan"))
                           / rho_w_surf) / (Td_s * Td_v)
    conv = int(oa.get("conv", 1)) * int(orr.get("conv", 1)) * \
        int(on.get("conv", 1))
    rec["status"] = "ok" if conv == 1 else "unconverged"
    rec["wall_s"] = round(time.time() - t0, 2)
    return rec


def hms(sec):
    if sec is None or sec != sec or sec < 0:
        return "--:--:--"
    sec = int(sec)
    return "%02d:%02d:%02d" % (sec // 3600, (sec % 3600) // 60, sec % 60)


OUT_FIELDS = ["case_id", "band_nm",
              "sza_deg", "vza_deg", "raa_deg", "scat_angle_deg",
              "glint_tilt_deg", "wind_ms",
              "chl_mg_m3", "tsm_g_m3", "acdom440_m_inv",
              "aer_model", "aod865", "AOD_band",
              "Rho_TOA", "Rho_R", "Rho_C", "Rho_Am",
              "Td_R_s", "Td_R_v", "Td_Am_s", "Td_Am_v",
              "Td_s", "Td_v", "Rrs", "Rrs_view", "rrs0minus",
              "t_rho_w", "Td_check", "orders", "status", "wall_s", "note"]


def main():
    ap = argparse.ArgumentParser(description="대기보정 검증용 표본 생성기")
    ap.add_argument("--root", default=None)
    ap.add_argument("--exe", default=None)
    ap.add_argument("--out", default=None, help="산출물 폴더 (기본 <root>/ac_set)")
    ap.add_argument("--n", type=int, default=1000, help="사례 수")
    ap.add_argument("--seed", type=int, default=20260804)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--covary", action="store_true",
                    help="해수 세 값을 실제 바다처럼 함께 움직이게 뽑는다")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = find_root(args.root)
    if not root:
        sys.exit("오류: OCRT 패키지 최상위를 찾지 못했다. --root 로 지정한다.")
    exe = find_exe(args.exe, root)
    if not exe:
        sys.exit("오류: 실행 파일을 찾지 못했다. --exe 로 지정한다.")
    out_dir = os.path.abspath(args.out or os.path.join(root, "ac_set"))
    os.makedirs(out_dir, exist_ok=True)

    aer_models = sorted(os.path.splitext(os.path.basename(p))[0]
                        for p in glob.glob(os.path.join(root, AER_DIR, "*.mie")))
    if not aer_models:
        sys.exit("오류: 에어로졸 파일이 없다: %s" % AER_DIR)

    design = make_design(args.n, args.seed, aer_models, args.covary)
    dpath = os.path.join(out_dir, "design.csv")
    with open(dpath, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(design[0].keys()))
        w.writeheader(); w.writerows(design)

    total = len(design) * len(BANDS_NM)
    print("대기보정 검증용 표본 생성 계획")
    print("  패키지 폴더  : %s" % root)
    print("  산출물 폴더  : %s" % out_dir)
    print("  사례 %d개 x 밴드 %d개 = %d 조합, 조합마다 실행 5회"
          % (len(design), len(BANDS_NM), total))
    print("  밴드         : %s nm" % ", ".join(str(b) for b in BANDS_NM))
    print("  에어로졸     : %d종 무작위, 광학두께 %g~%g" % (len(aer_models), *AOD_RANGE))
    print("  기하         : 태양천정각 %g~%g, 관측천정각 %g~%g, 방위각 %g~%g도"
          % (*SZA_RANGE, *VZA_RANGE, *RAA_RANGE))
    print("  해수         : 엽록소 %g~%g, 총부유물 %g~%g, CDOM %g~%g%s"
          % (*CHL_RANGE, *TSM_RANGE, *CDOM_RANGE,
             " (함께 움직임)" if args.covary else " (각각 독립, 로그 균등)"))
    print("  풍속         : %g~%g m/s" % WIND_RANGE)
    print("  Rrs          : 천저 방향 (관측 방향 값도 함께 기록)")
    print("  설계표       : %s" % dpath)
    est = total * NOMINAL_SEC / max(1, args.workers)
    print("  예상 소요    : 약 %s (조합당 %.0f초)" % (hms(est), NOMINAL_SEC))

    if args.dry_run:
        print("\n--dry-run 이므로 실행하지 않는다.")
        return

    res_path = os.path.join(out_dir, "results.csv")
    done = set()
    if os.path.isfile(res_path):
        with open(res_path, "r", encoding="utf-8", errors="replace",
                  newline="") as f:
            for row in csv.DictReader(f):
                if row.get("status") in ("ok", "unconverged"):
                    done.add((int(float(row["case_id"])),
                              int(float(row["band_nm"]))))
        print("  이미 끝난 것 : %d 조합" % len(done))

    jobs = [(exe, root, r, b, args.timeout)
            for r in design for b in BANDS_NM
            if (r["case_id"], b) not in done]
    print("  이번에 돌릴 것: %d 조합\n" % len(jobs))
    if not jobs:
        print("돌릴 것이 없다.")
        return

    new_file = not os.path.isfile(res_path)
    fh = open(res_path, "a", encoding="utf-8", newline="")
    wr = csv.DictWriter(fh, fieldnames=OUT_FIELDS, extrasaction="ignore")
    if new_file:
        wr.writeheader(); fh.flush()
    dmap = {r["case_id"]: r for r in design}

    t0 = time.time(); durs = []; nfail = nunconv = 0; width = 0
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for k, res in enumerate(ex.map(run_one, jobs, chunksize=1), 1):
            row = dict(dmap[res["case_id"]]); row.update(res)
            wr.writerow(row)
            fh.flush()          # 한 건마다 기록한다. 중간에 멈춰도 남는다.
            if res.get("status") == "failed":
                nfail += 1
            elif res.get("status") == "unconverged":
                nunconv += 1
            if res.get("wall_s"):
                durs.append(res["wall_s"])
            avg = (sum(durs) / len(durs)) if durs else NOMINAL_SEC
            eta = (len(jobs) - k) * avg / max(1, args.workers)
            pct = 100.0 * k / len(jobs)
            nb = int(round(pct / 10))
            line = ("[%s] %d/%d %5.1f%% | 남은 %s | 완료 %s | 평균 %.0f초"
                    % ("#" * nb + "." * (10 - nb), k, len(jobs), pct, hms(eta),
                       (datetime.now() + timedelta(seconds=eta)).strftime("%m-%d %H:%M"),
                       avg))
            if nfail:
                line += " | 실패 %d" % nfail
            if nunconv:
                line += " | 미수렴 %d" % nunconv
            sys.stdout.write("\r" + line + " " * max(0, width - len(line)))
            sys.stdout.flush(); width = len(line)
    fh.close()
    print("\n")
    print("끝났다.")
    print("  결과 파일 : %s" % res_path)
    print("  완료 %d조합, 실패 %d, 미수렴 %d" % (len(jobs), nfail, nunconv))
    print("  총 경과 %s" % hms(time.time() - t0))
    if nfail:
        print("  같은 명령을 다시 돌리면 실패한 것만 재실행된다.")


if __name__ == "__main__":
    main()
