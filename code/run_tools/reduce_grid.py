#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OCRT 편광 민감도 결과 축약기

수천 개의 격자 결과 파일(전체 10 GB 안팎)에서 분석에 필요한 것만 뽑아
수십 MB 로 줄인다. 결과 파일 자체는 옮기지 않아도 된다.

만드는 것은 세 가지다.

  metrics.csv.gz     주 산출물. 두 측도 N 과 S 를 각도마다 계산해 둔 표.
                     한 행 = (밴드, 태양천정각, 풍속, N가지 광학두께,
                              관측천정각, 방위각)
  cases_thin.csv.gz  성긴 각도에서의 원값. 산포도를 그리거나 계산을
                     다시 확인할 때 쓴다.
  runs_summary.csv   실행 한 건당 요약 한 줄. 품질 점검용.

쓰는 법
  python reduce_grid.py --run-dir D:\\temp\\star_run --out D:\\temp\\reduced
  python reduce_grid.py --run-dir D:\\temp\\star_run --out D:\\temp\\reduced --workers 12
"""

import argparse
import csv
import glob
import gzip
import os
import sys
import time
from datetime import datetime

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
if os.name == "nt":
    try:
        import ctypes
        ctypes.windll.kernel32.SetConsoleOutputCP(65001)
    except Exception:
        pass

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("오류: numpy 와 pandas 가 필요하다. 아나콘다 프롬프트에서 실행한다.")


# 격자 파일에서 실제로 읽는 열. 60개 중 이것만 읽으면 된다.
USECOLS = ["vza_deg", "raa_deg",
           "TOA_rho_I", "TOA_rho_Q", "TOA_rho_U",
           "rho_I", "rho_Q", "rho_U",           # 해수 없는 기준 대기의 열 이름
           "TOA_water_signal_I", "Rrs_I", "Ed0plus_air"]

# 설계 축. 기록 파일에서 가져온다.
AXES = ["band_nm", "sza_deg", "wind_ms", "aod865", "aer_model", "ssa443",
        "chl", "tsm", "cdom", "branch", "surface"]


# ===========================================================================
# 기하 보조
# ===========================================================================

def scattering_angle_deg(sza, vza, raa):
    """OCRT 방위각 규약에서의 산란각.

    cos(theta) = -cos(sza)cos(vza) - sin(sza)sin(vza)cos(raa)
    레일리 대기에서 선형편광도 최댓값이 정확히 90도에 오는 것으로 확인했다.
    """
    s, v, r = np.radians(sza), np.radians(vza), np.radians(raa)
    c = -np.cos(s) * np.cos(v) - np.sin(s) * np.sin(v) * np.cos(r)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def glint_tilt_deg(sza, vza, raa):
    """직달 선글린트를 만들려면 해수면이 얼마나 기울어야 하는가(도).

    0 이면 잔잔한 수면에서 곧바로 선글린트가 생기는 방향이다.
    태양천정각 40도에서 이 값이 0 인 지점이 관측천정각 40도·방위각 180도로
    나오고, 선글린트를 켠 실행의 세기 최대 지점과 일치함을 확인했다.
    """
    s, v, r = np.radians(sza), np.radians(vza), np.radians(raa)
    sv = np.cos(s) * np.cos(v) + np.sin(s) * np.sin(v) * np.cos(r)
    num = np.cos(s) + np.cos(v)
    den = np.sqrt(np.maximum(2.0 + 2.0 * sv, 1e-300))
    return np.degrees(np.arccos(np.clip(num / den, -1.0, 1.0)))


# ===========================================================================
# 기록 파일 읽기
# ===========================================================================

def hhmm(sec):
    sec = int(max(sec, 0))
    return "%d분 %02d초" % (sec // 60, sec % 60)


def load_index(run_dir):
    """모든 기록 파일을 읽어 하나로 만든다."""
    pats = ["index.csv", "index_*.csv"]
    files = []
    for p in pats:
        files += glob.glob(os.path.join(run_dir, p))
    files = sorted(set(files))
    if not files:
        sys.exit("오류: 기록 파일(index*.csv)이 없다: %s" % run_dir)

    rows = []
    for f in files:
        with open(f, "r", encoding="utf-8", errors="replace", newline="") as fh:
            sub = list(csv.DictReader(fh))
        rows.extend(sub)
        print("  %s (%d 행)" % (os.path.basename(f), len(sub)))

    df = pd.DataFrame(rows)
    # 결과 파일 이름만 뽑는다. 윈도우 경로와 리눅스 경로를 모두 받는다.
    df["fname"] = (df["out_path"].astype(str)
                   .str.replace("\\", "/", regex=False)
                   .str.rsplit("/", n=1).str[-1])

    # 같은 실행이 여러 번 기록돼 있으면 성공한 것을 남긴다.
    rank = {"ok": 0, "unconverged": 1, "failed": 2}
    df["_rank"] = df["status"].map(lambda x: rank.get(x, 3))
    df = df.sort_values("_rank").drop_duplicates("fname", keep="first")

    if "surface" not in df.columns:
        df["surface"] = "ocean"
    df["surface"] = df["surface"].fillna("ocean").replace("", "ocean")
    for c in ["band_nm", "sza_deg", "wind_ms", "aod865", "ssa443",
              "chl", "tsm", "cdom"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.reset_index(drop=True)


def resolve_paths(df, run_dir):
    """결과 파일의 실제 위치를 찾는다."""
    base = os.path.join(run_dir, "out", "ocn")
    df["path"] = df["fname"].map(lambda n: os.path.join(base, n))
    have = df["path"].map(os.path.isfile)
    if not have.all():
        missing = int((~have).sum())
        print("")
        print("경고: 결과 파일 %d 개를 찾지 못했다. 그 실행은 건너뛴다." % missing)
        for n in df.loc[~have, "fname"].head(5):
            print("       %s" % n)
    return df[have].reset_index(drop=True)


# ===========================================================================
# 격자 파일 읽기
# ===========================================================================

def read_grid(path):
    """격자 파일에서 필요한 열만 읽는다."""
    d = pd.read_csv(path, usecols=lambda c: c.strip() in USECOLS)
    d.columns = [c.strip() for c in d.columns]
    return d


def _worker(args):
    """한 파일을 읽어 필요한 배열만 돌려준다."""
    i, path = args
    try:
        d = read_grid(path)
    except Exception as e:
        return i, None, str(e)
    pre = "TOA_rho_I" if "TOA_rho_I" in d.columns else "rho_I"
    tag = "TOA_rho_" if pre == "TOA_rho_I" else "rho_"
    out = dict(
        vza=d["vza_deg"].to_numpy(np.float64),
        raa=d["raa_deg"].to_numpy(np.float64),
        I=d[tag + "I"].to_numpy(np.float64),
        Q=d[tag + "Q"].to_numpy(np.float64),
        U=d[tag + "U"].to_numpy(np.float64),
    )
    for extra, key in [("TOA_water_signal_I", "wI"), ("Rrs_I", "rrs"),
                       ("Ed0plus_air", "ed")]:
        if extra in d.columns:
            out[key] = d[extra].to_numpy(np.float64)
    return i, out, None


# ===========================================================================
# 두 측도 계산
# ===========================================================================

def compute_N(I, P):
    """해수를 훑을 때 세기로 설명되지 않는 편광 변화.

    I : (해수종수, 각도수)  세기
    P : (해수종수, 각도수)  편광 반사도 sqrt(Q^2+U^2)

    각 각도에서 P 를 I 의 일차식에 맞추고, 남는 흩어짐을 표준편차로 잰다.
    함께 결정계수와 기울기도 낸다.
    """
    n, m = I.shape
    resid = np.full(m, np.nan)
    r2 = np.full(m, np.nan)
    slope = np.full(m, np.nan)
    if n < 3:
        return resid, r2, slope
    Ic = I - I.mean(0)
    Pc = P - P.mean(0)
    sxx = (Ic * Ic).sum(0)
    sxy = (Ic * Pc).sum(0)
    syy = (Pc * Pc).sum(0)
    good = sxx > 0
    b = np.zeros(m)
    b[good] = sxy[good] / sxx[good]
    ss = syy - b * sxy                      # 잔차 제곱합
    ss = np.maximum(ss, 0.0)
    resid = np.sqrt(ss / max(n - 2, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where(syy > 0, 1.0 - ss / syy, np.nan)
    slope = b
    return resid, r2, slope


def compute_S(I, P, model_ids, aods):
    """세기가 같은 지점에서 에어로졸 모델 사이에 편광이 얼마나 벌어지는가.

    모델마다 광학두께를 훑어 (세기, 편광) 곡선을 만들고, 모든 모델이
    함께 덮는 세기 구간의 한가운데에서 편광을 읽어 최대와 최소의 차를 잰다.
    겹치는 구간이 없으면 그 각도는 결측으로 둔다.
    """
    models = sorted(set(model_ids))
    m = I.shape[1]
    spread = np.full(m, np.nan)
    overlap = np.full(m, np.nan)
    if len(models) < 2:
        return spread, overlap

    # 모델별로 광학두께 순서대로 정렬해 둔다.
    curves = []
    for mo in models:
        sel = [k for k, x in enumerate(model_ids) if x == mo]
        order = np.argsort([aods[k] for k in sel])
        idx = [sel[k] for k in order]
        curves.append((I[idx, :], P[idx, :]))

    lo = np.maximum.reduce([c[0].min(0) for c in curves])
    hi = np.minimum.reduce([c[0].max(0) for c in curves])
    span_all = np.maximum.reduce([c[0].max(0) for c in curves]) - \
        np.minimum.reduce([c[0].min(0) for c in curves])
    ok = hi > lo
    mid = 0.5 * (lo + hi)

    vals = np.full((len(curves), m), np.nan)
    for ci, (Ic, Pc) in enumerate(curves):
        for k in range(m):
            if not ok[k]:
                continue
            x, y = Ic[:, k], Pc[:, k]
            o = np.argsort(x)
            vals[ci, k] = np.interp(mid[k], x[o], y[o])
    with np.errstate(invalid="ignore"):
        spread = np.nanmax(vals, 0) - np.nanmin(vals, 0)
        overlap = np.where(span_all > 0, (hi - lo) / span_all, np.nan)
    spread[~ok] = np.nan
    return spread, overlap


def compute_S_pairwise(I, P, model_ids, aods, ssas):
    """모델을 둘씩 짝지어 세기가 겹치는 구간에서 편광 차를 잰다.

    모든 모델이 함께 덮는 구간을 요구하지 않으므로 모델 수가 늘어도
    측정이 가능하다. 짝마다의 최댓값과, 단일산란알베도 차가 가장 큰
    짝의 값을 함께 낸다.
    """
    models = sorted(set(model_ids))
    m = I.shape[1]
    if len(models) < 2:
        return np.full(m, np.nan), np.full(m, np.nan), np.full(m, np.nan)

    curves = {}
    ssa_of = {}
    for mo in models:
        sel = [k for k, x in enumerate(model_ids) if x == mo]
        order = np.argsort([aods[k] for k in sel])
        idx = [sel[k] for k in order]
        curves[mo] = (I[idx, :], P[idx, :])
        ssa_of[mo] = ssas[sel[0]]

    best = np.zeros(m)
    npair = np.zeros(m)
    extreme = np.full(m, np.nan)
    lo_m = min(models, key=lambda x: ssa_of[x])
    hi_m = max(models, key=lambda x: ssa_of[x])

    for a in range(len(models)):
        for b in range(a + 1, len(models)):
            ma, mb = models[a], models[b]
            Ia, Pa = curves[ma]
            Ib, Pb = curves[mb]
            lo = np.maximum(Ia.min(0), Ib.min(0))
            hi = np.minimum(Ia.max(0), Ib.max(0))
            ok = hi > lo
            if not ok.any():
                continue
            mid = 0.5 * (lo + hi)
            d = np.full(m, np.nan)
            for k in np.where(ok)[0]:
                oa = np.argsort(Ia[:, k])
                ob = np.argsort(Ib[:, k])
                d[k] = abs(np.interp(mid[k], Ia[oa, k], Pa[oa, k]) -
                           np.interp(mid[k], Ib[ob, k], Pb[ob, k]))
            npair += ok
            best = np.where(ok & (d > best), d, best)
            if {ma, mb} == {lo_m, hi_m}:
                extreme = d
    best[npair == 0] = np.nan
    return best, npair, extreme


def compute_S_dolp(I, P, model_ids, aods, ref_aod):
    """광학두께를 한 값으로 고정했을 때 모델 사이 선형편광도의 벌어짐.

    세기를 맞출 필요가 없어 언제나 정의된다. 위의 S 를 보완한다.
    """
    m = I.shape[1]
    picks = [k for k, a in enumerate(aods) if abs(a - ref_aod) < 1e-9]
    if len(picks) < 2:
        return np.full(m, np.nan)
    d = P[picks, :] / np.maximum(I[picks, :], 1e-300)
    return d.max(0) - d.min(0)


# ===========================================================================
# 본체
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="OCRT 편광 민감도 결과 축약기")
    ap.add_argument("--run-dir", required=True,
                    help="실행 산출물 폴더 (index*.csv 와 out/ocn 이 있는 곳)")
    ap.add_argument("--out", required=True, help="축약 결과를 둘 폴더")
    ap.add_argument("--workers", type=int, default=0,
                    help="파일 읽기 병렬 수 (기본: 코어 수의 절반)")
    ap.add_argument("--thin-vza", type=float, default=15.0,
                    help="원값 파일에 남길 관측천정각 간격(도). 기본 15")
    ap.add_argument("--thin-raa", type=float, default=45.0,
                    help="원값 파일에 남길 방위각 간격(도). 기본 45")
    ap.add_argument("--no-thin", action="store_true",
                    help="원값 파일을 만들지 않는다")
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    print("기록 파일 읽기")
    idx = load_index(run_dir)
    idx = resolve_paths(idx, run_dir)
    print("")
    print("  실행 %d 건" % len(idx))
    ok = (idx["status"] == "ok").sum()
    if ok != len(idx):
        print("  그중 성공 %d 건, 나머지는 그대로 포함한다" % ok)

    nw = args.workers or max(1, (os.cpu_count() or 4) // 2)

    # --- 격자 파일을 모두 읽는다 ------------------------------------------
    print("")
    print("격자 파일 읽기 (동시 %d개)" % nw)
    data = [None] * len(idx)
    errs = []
    jobs = list(enumerate(idx["path"].tolist()))

    def report(done, total):
        el = time.time() - t0
        eta = el / max(done, 1) * (total - done)
        sys.stdout.write("\r  %d/%d  경과 %s  남은 %s      "
                         % (done, total, hhmm(el), hhmm(eta)))
        sys.stdout.flush()

    if nw > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=nw) as ex:
            for n, (i, d, e) in enumerate(ex.map(_worker, jobs, chunksize=8), 1):
                data[i] = d
                if e:
                    errs.append((idx["fname"].iloc[i], e))
                if n % 25 == 0 or n == len(jobs):
                    report(n, len(jobs))
    else:
        for n, j in enumerate(jobs, 1):
            i, d, e = _worker(j)
            data[i] = d
            if e:
                errs.append((idx["fname"].iloc[i], e))
            if n % 25 == 0 or n == len(jobs):
                report(n, len(jobs))
    print("")
    if errs:
        print("  읽기 실패 %d 건" % len(errs))
        for n, e in errs[:5]:
            print("    %s : %s" % (n, e[:120]))

    keep = [i for i in range(len(idx)) if data[i] is not None]
    idx = idx.iloc[keep].reset_index(drop=True)
    data = [data[i] for i in keep]

    # 각도 격자는 모든 실행에서 같아야 한다.
    vza = data[0]["vza"]
    raa = data[0]["raa"]
    same = all(len(d["vza"]) == len(vza) for d in data)
    if not same:
        sys.exit("오류: 실행마다 각도 격자가 다르다. 같은 격자로 다시 돌려야 한다.")

    # --- 실행 한 건당 요약 -------------------------------------------------
    print("")
    print("실행별 요약 만들기")
    P_all = [np.hypot(d["Q"], d["U"]) for d in data]
    summ = idx[[c for c in AXES if c in idx.columns] +
               ["fname", "status", "wall_s", "n_rows", "n_unconverged",
                "max_water_orders", "code_version", "binary_sha256"]].copy()
    summ["I_mean"] = [d["I"].mean() for d in data]
    summ["I_max"] = [d["I"].max() for d in data]
    summ["P_mean"] = [p.mean() for p in P_all]
    summ["P_max"] = [p.max() for p in P_all]
    summ["dolp_mean"] = [(p / np.maximum(d["I"], 1e-300)).mean()
                         for p, d in zip(P_all, data)]
    summ["dolp_max"] = [(p / np.maximum(d["I"], 1e-300)).max()
                        for p, d in zip(P_all, data)]
    if any("wI" in d for d in data):
        summ["water_frac_I_mean"] = [
            (d["wI"] / np.maximum(d["I"], 1e-300)).mean() if "wI" in d else np.nan
            for d in data]
    summ.to_csv(os.path.join(out_dir, "runs_summary.csv"), index=False)
    print("  runs_summary.csv  (%d 행)" % len(summ))

    # --- 두 측도 계산 ------------------------------------------------------
    print("")
    print("측도 N 과 S 계산")
    idx["_row"] = np.arange(len(idx))
    ocean = idx["surface"].astype(str).eq("ocean")
    is_N = ocean & idx["branch"].astype(str).str.contains("N")
    is_S = ocean & idx["branch"].astype(str).str.contains("S")
    is_A = idx["surface"].astype(str).eq("black")

    # S 는 (밴드, 태양천정각, 풍속) 마다 하나씩 나온다.
    S_map = {}
    for key, g in idx[is_S].groupby(["band_nm", "sza_deg", "wind_ms"]):
        rows = g["_row"].to_numpy()
        I = np.vstack([data[r]["I"] for r in rows])
        P = np.vstack([P_all[r] for r in rows])
        mids = g["aer_model"].tolist()
        aods = g["aod865"].tolist()
        sp, ov = compute_S(I, P, mids, aods)
        ref = 0.10 if 0.10 in set(aods) else sorted(aods)[len(aods) // 2]
        sd = compute_S_dolp(I, P, mids, aods, ref)
        ssas = g["ssa443"].tolist()
        spair, npair, sext = compute_S_pairwise(I, P, mids, aods, ssas)
        S_map[key] = dict(S_spread=sp, S_overlap=ov, S_dolp=sd,
                          S_pair_max=spair, S_pair_n=npair, S_pair_extreme=sext,
                          n_aer=len(set(mids)), n_aer_runs=len(rows),
                          S_ref_aod=ref)

    # N 은 (밴드, 태양천정각, 풍속, 광학두께) 마다 하나씩 나온다.
    out_rows = []
    for key, g in idx[is_N].groupby(["band_nm", "sza_deg", "wind_ms", "aod865"]):
        band, sza, wind, aod = key
        rows = g["_row"].to_numpy()
        I = np.vstack([data[r]["I"] for r in rows])
        P = np.vstack([P_all[r] for r in rows])
        resid, r2, slope = compute_N(I, P)

        s = S_map.get((band, sza, wind))
        m = len(vza)
        blk = pd.DataFrame(dict(
            band_nm=band, sza_deg=sza, wind_ms=wind, n_aod865=aod,
            vza_deg=vza, raa_deg=raa,
            theta_deg=scattering_angle_deg(sza, vza, raa),
            glint_tilt_deg=glint_tilt_deg(sza, vza, raa),
            N_resid=resid, N_r2=r2, N_slope=slope,
            I_mean=I.mean(0), I_span=I.max(0) - I.min(0),
            P_mean=P.mean(0), P_span=P.max(0) - P.min(0),
            dolp_mean=(P / np.maximum(I, 1e-300)).mean(0),
            n_water=len(rows),
        ))
        if s is not None:
            for c in ["S_spread", "S_overlap", "S_dolp",
                      "S_pair_max", "S_pair_n", "S_pair_extreme"]:
                blk[c] = s[c]
            blk["n_aer"] = s["n_aer"]
            blk["S_ref_aod"] = s["S_ref_aod"]
        else:
            for c in ["S_spread", "S_overlap", "S_dolp", "S_pair_max",
                      "S_pair_n", "S_pair_extreme", "n_aer", "S_ref_aod"]:
                blk[c] = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            blk["R"] = blk["S_spread"] / blk["N_resid"]
            blk["R_pair"] = blk["S_pair_max"] / blk["N_resid"]
        out_rows.append(blk)

    met = pd.concat(out_rows, ignore_index=True)
    thin_fam = met[met.n_water < 8][["band_nm", "sza_deg", "wind_ms",
                                     "n_aod865", "n_water"]].drop_duplicates()
    if len(thin_fam):
        print("")
        print("  주의: 해수 표본이 8종 미만인 무리가 %d개 있다." % len(thin_fam))
        print("        중간에 멈춘 실행이 섞이면 생긴다. N 측도를 믿기 어려우므로")
        print("        분석에서 n_water 열로 걸러 쓴다.")
        for _, r in thin_fam.head(6).iterrows():
            print("        밴드 %d, 태양천정각 %g도, 풍속 %g m/s, 광학두께 %g → 해수 %d종"
                  % (r.band_nm, r.sza_deg, r.wind_ms, r.n_aod865, r.n_water))
    mp = os.path.join(out_dir, "metrics.csv.gz")
    met.to_csv(mp, index=False, compression="gzip", float_format="%.8g")
    print("  metrics.csv.gz    (%d 행, %.1f MB)"
          % (len(met), os.path.getsize(mp) / 1e6))

    # --- 해수 없는 기준 대기 (있을 때만) ------------------------------------
    if is_A.any():
        print("")
        print("해수 없는 기준 대기 정리")
        parts = []
        for r in idx[is_A]["_row"].to_numpy():
            d = data[r]
            n = len(d["vza"])
            rec = {c: np.repeat(idx[c].iloc[r], n)
                   for c in ["band_nm", "sza_deg", "wind_ms", "aod865",
                             "aer_model", "ssa443"]}
            rec["vza_deg"] = d["vza"]
            rec["raa_deg"] = d["raa"]
            rec["atm_I"] = d["I"]
            rec["atm_Q"] = d["Q"]
            rec["atm_U"] = d["U"]
            parts.append(pd.DataFrame(rec))
        atm = pd.concat(parts, ignore_index=True)
        ap = os.path.join(out_dir, "atm_reference.csv.gz")
        atm.to_csv(ap, index=False, compression="gzip", float_format="%.8g")
        print("  atm_reference.csv.gz (%d 행, %.1f MB)"
              % (len(atm), os.path.getsize(ap) / 1e6))

    # --- 성긴 각도의 원값 ---------------------------------------------------
    if not args.no_thin:
        print("")
        print("성긴 각도 원값 만들기")
        pick = (np.isclose(vza % args.thin_vza, 0) |
                np.isclose(vza % args.thin_vza, args.thin_vza)) & \
               (np.isclose(raa % args.thin_raa, 0) |
                np.isclose(raa % args.thin_raa, args.thin_raa))
        pick = np.where(pick)[0]
        print("  각도 %d개 중 %d개를 남긴다" % (len(vza), len(pick)))

        cols = [c for c in AXES if c in idx.columns]
        parts = []
        thin_rows = np.where(ocean.to_numpy())[0]
        for i in thin_rows:
            d = data[i]
            n = len(pick)
            rec = {c: np.repeat(idx[c].iloc[i], n) for c in cols}
            rec["fname"] = np.repeat(idx["fname"].iloc[i], n)
            rec["vza_deg"] = d["vza"][pick]
            rec["raa_deg"] = d["raa"][pick]
            rec["TOA_rho_I"] = d["I"][pick]
            rec["TOA_rho_Q"] = d["Q"][pick]
            rec["TOA_rho_U"] = d["U"][pick]
            for k, name in [("wI", "TOA_water_signal_I"), ("rrs", "Rrs_I")]:
                if k in d:
                    rec[name] = d[k][pick]
            parts.append(pd.DataFrame(rec))
        thin = pd.concat(parts, ignore_index=True)
        thin["theta_deg"] = scattering_angle_deg(
            thin["sza_deg"].to_numpy(np.float64),
            thin["vza_deg"].to_numpy(np.float64),
            thin["raa_deg"].to_numpy(np.float64))
        thin["glint_tilt_deg"] = glint_tilt_deg(
            thin["sza_deg"].to_numpy(np.float64),
            thin["vza_deg"].to_numpy(np.float64),
            thin["raa_deg"].to_numpy(np.float64))
        tp = os.path.join(out_dir, "cases_thin.csv.gz")
        thin.to_csv(tp, index=False, compression="gzip", float_format="%.8g")
        print("  cases_thin.csv.gz (%d 행, %.1f MB)"
              % (len(thin), os.path.getsize(tp) / 1e6))

    # --- 안내문 -------------------------------------------------------------
    note = os.path.join(out_dir, "README_reduced.md")
    total = sum(os.path.getsize(os.path.join(out_dir, f))
                for f in os.listdir(out_dir))
    with open(note, "w", encoding="utf-8") as f:
        f.write(NOTE_TEMPLATE % dict(
            when=datetime.now().strftime("%Y-%m-%d %H:%M"),
            run_dir=run_dir, n_runs=len(idx),
            n_angles=len(vza), thin_v=args.thin_vza, thin_r=args.thin_raa,
            n_met=len(met), total_mb=total / 1e6))

    print("")
    print("끝났다. 총 %.1f MB (원자료 대비 크게 줄었다)" % (total / 1e6))
    print("  폴더: %s" % out_dir)
    print("  이 폴더만 올리면 분석을 이어갈 수 있다.")


NOTE_TEMPLATE = """# 축약 결과 안내

- 만든 때: %(when)s
- 원본 폴더: %(run_dir)s
- 실행 %(n_runs)d 건, 각 실행의 각도 격자 %(n_angles)d 방향

## metrics.csv.gz

주 산출물이다. 한 행이 (밴드, 태양천정각, 풍속, N가지 광학두께, 관측천정각,
방위각) 하나에 해당한다. 모두 %(n_met)d 행이다.

| 열 | 뜻 |
|---|---|
| `theta_deg` | 산란각. cos(theta) = -cos(sza)cos(vza) - sin(sza)sin(vza)cos(raa) |
| `glint_tilt_deg` | 직달 선글린트가 생기려면 수면이 기울어야 하는 각도. 0 이면 잔잔한 수면에서 곧바로 선글린트가 생기는 방향이다 |
| `N_resid` | 해수를 훑을 때 세기로 설명되지 않는 편광 변화 (표준편차) |
| `N_r2` | 세기 하나로 편광을 설명하는 정도 |
| `N_slope` | 세기에 대한 편광의 기울기 |
| `S_spread` | 세기가 같은 지점에서 에어로졸 모델 사이 편광 벌어짐 |
| `S_overlap` | 모델들이 함께 덮는 세기 구간의 비율. 낮으면 `S_spread` 를 믿기 어렵다 |
| `S_pair_max` | 모델을 둘씩 짝지어 세기가 겹치는 구간에서 잰 편광 차의 최댓값. 모든 모델의 공통 구간을 요구하지 않아 언제나 잴 수 있다 |
| `S_pair_n` | 겹침이 있어 계산에 쓰인 짝의 수 |
| `S_pair_extreme` | 단일산란알베도가 가장 낮은 모델과 가장 높은 모델의 짝에서 잰 값 |
| `S_dolp` | 광학두께를 고정했을 때 모델 사이 선형편광도 벌어짐 |
| `R` | 유효도 = `S_spread` / `N_resid` |
| `R_pair` | 유효도 = `S_pair_max` / `N_resid`. 결측이 적어 지도 그리기에 알맞다 |
| `I_mean` `P_mean` `dolp_mean` | 해수 훑기의 평균 세기, 편광 반사도, 선형편광도 |
| `I_span` `P_span` | 해수 훑기에서의 변동폭 |
| `n_water` `n_aer` | 그 계산에 쓰인 해수 종수와 에어로졸 모델 수 |

## atm_reference.csv.gz

해수가 없는 기준 대기(흑색 프레넬 해면)의 전 각도 반사도다. 같은 조건의
결합 실행에서 이 값을 빼면 해수가 더한 몫이 나온다. 편광에서도 세기와
같은 방식으로 뺄 수 있다. 열은 `atm_I`, `atm_Q`, `atm_U` 이며 설계 축과
각도로 맞추면 된다.

편광 반사도는 sqrt(Q^2 + U^2), 선형편광도는 그 값을 I 로 나눈 것이다.

## cases_thin.csv.gz

성긴 각도(관측천정각 %(thin_v)g도, 방위각 %(thin_r)g도 간격)에서의 원값이다.
설계 축과 함께 TOA 반사도 I, Q, U 가 들어 있다. 산포도를 그리거나 위의
계산을 다시 확인할 때 쓴다.

## runs_summary.csv

실행 한 건에 한 줄이다. 상태, 소요, 미수렴 행 수, 산란차수, 코드 판번호와
실행 파일 해시가 들어 있다. 품질 점검용이다.

## 총 용량

%(total_mb).1f MB
"""


if __name__ == "__main__":
    main()
