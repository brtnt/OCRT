#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
조각 나눠 돌린 결과를 합쳐 하나로 만든다.

세 조각(또는 몇 조각이든)을 각기 다른 폴더나 기계에서 돌린 뒤, 그 폴더들을
인자로 주면 다음을 한다.

  1. 조각별 기록 파일(index_part*.csv)을 하나로 합친다.
  2. 실패, 미수렴, 빠진 실행, 중복 실행을 찾아 알린다.
  3. 결과 파일을 한 폴더로 모은다(원하는 경우).
  4. 분석에 바로 쓸 수 있는 하나의 index.csv 를 만든다.

쓰는 법
  python merge_parts.py D:\\run_a D:\\run_b D:\\run_c --out D:\\merged
  python merge_parts.py D:\\run_a D:\\run_b D:\\run_c --out D:\\merged --copy
"""

import argparse
import csv
import glob
import os
import shutil
import sys

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


def read_index_files(folders):
    """조각 폴더들에서 기록 파일을 모두 읽는다."""
    rows, files = [], []
    for d in folders:
        d = os.path.abspath(d)
        if not os.path.isdir(d):
            print("경고: 폴더가 없다: %s" % d)
            continue
        found = sorted(glob.glob(os.path.join(d, "index_part*.csv"))) \
            + sorted(glob.glob(os.path.join(d, "index.csv")))
        if not found:
            print("경고: 기록 파일이 없다: %s" % d)
            continue
        for f in found:
            with open(f, "r", encoding="utf-8", errors="replace", newline="") as fh:
                sub = list(csv.DictReader(fh))
            for r in sub:
                r["_src_folder"] = d
                r["_src_index"] = os.path.basename(f)
            rows.extend(sub)
            files.append((f, len(sub)))
    return rows, files


def main():
    ap = argparse.ArgumentParser(description="조각 결과 합치기")
    ap.add_argument("folders", nargs="+", help="조각별 산출물 폴더")
    ap.add_argument("--out", required=True, help="합친 결과를 둘 폴더")
    ap.add_argument("--copy", action="store_true",
                    help="결과 CSV 파일도 한 폴더로 복사한다")
    ap.add_argument("--expect", type=int, default=615,
                    help="전체 실행 수 (기본 615)")
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)

    rows, files = read_index_files(args.folders)
    if not rows:
        sys.exit("합칠 기록이 하나도 없다.")

    print("읽은 기록 파일")
    for f, n in files:
        print("  %s  (%d 행)" % (f, n))
    print("")

    # --- 중복 정리: 같은 실행이 여러 번 있으면 성공한 것을 남긴다 -----------
    best = {}
    dup = 0
    for r in rows:
        key = os.path.basename(r.get("out_path", ""))
        if not key:
            continue
        if key in best:
            dup += 1
            old = best[key]
            rank = {"ok": 0, "unconverged": 1, "failed": 2}
            if rank.get(r.get("status"), 3) < rank.get(old.get("status"), 3):
                best[key] = r
        else:
            best[key] = r

    merged = list(best.values())
    merged.sort(key=lambda r: (int(float(r.get("band_nm", 0))),
                               float(r.get("sza_deg", 0)),
                               r.get("out_path", "")))

    ok = [r for r in merged if r.get("status") == "ok"]
    unconv = [r for r in merged if r.get("status") == "unconverged"]
    failed = [r for r in merged if r.get("status") == "failed"]

    # --- 결과 파일이 실제로 있는지 확인 -------------------------------------
    missing_file = []
    for r in merged:
        p = r.get("out_path", "")
        if not p:
            continue
        if not os.path.isfile(p):
            # 폴더가 옮겨졌을 수 있으니 조각 폴더 기준으로도 찾아 본다
            alt = os.path.join(r["_src_folder"], "out", "ocn",
                               os.path.basename(p))
            if os.path.isfile(alt):
                r["out_path"] = alt
            else:
                missing_file.append(os.path.basename(p))

    # --- 요약 ---------------------------------------------------------------
    print("합친 결과")
    print("  전체 기록      : %d 건" % len(merged))
    print("  성공           : %d 건" % len(ok))
    if unconv:
        print("  미수렴         : %d 건" % len(unconv))
    if failed:
        print("  실패           : %d 건" % len(failed))
    if dup:
        print("  중복 기록      : %d 건 (성공한 쪽을 남겼다)" % dup)
    if missing_file:
        print("  결과 파일 없음 : %d 건" % len(missing_file))

    short = args.expect - len(merged)
    if short > 0:
        print("  빠진 실행      : %d 건 (기대 %d 건)" % (short, args.expect))
    elif short < 0:
        print("  기대보다 많다  : %d 건 (기대 %d 건)" % (-short, args.expect))

    # --- 판본 일치 확인 -----------------------------------------------------
    versions = sorted({r.get("code_version", "") for r in merged})
    hashes = sorted({r.get("binary_sha256", "") for r in merged})
    if len(versions) > 1 or len(hashes) > 1:
        print("")
        print("경고: 조각마다 실행 파일이 다르다. 결과를 섞으면 안 된다.")
        for v in versions:
            print("  판번호: %s" % v[:70])
        for h in hashes:
            print("  해시  : %s" % h[:16])
    else:
        print("  실행 파일      : 모든 조각이 같다 (%s)" % hashes[0][:16])

    grids = sorted({(r.get("vza_step_deg", ""), r.get("raa_step_deg", ""))
                    for r in merged})
    if len(grids) > 1:
        print("")
        print("경고: 조각마다 각도 격자가 다르다.")
        for g in grids:
            print("  관측천정각 간격 %s도 / 방위각 간격 %s도" % g)
    else:
        print("  각도 격자      : 모든 조각이 같다 (관측천정각 %s도, 방위각 %s도)"
              % grids[0])

    # --- 합친 기록 쓰기 -----------------------------------------------------
    fields = [k for k in merged[0].keys() if not k.startswith("_")]
    idx_out = os.path.join(out_dir, "index.csv")
    with open(idx_out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in merged:
            w.writerow(r)
    print("")
    print("합친 기록 : %s" % idx_out)

    # --- 문제 목록 ----------------------------------------------------------
    problems = failed + unconv
    if problems or missing_file:
        prob_out = os.path.join(out_dir, "problems.csv")
        with open(prob_out, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in problems:
                w.writerow(r)
        print("문제 목록 : %s" % prob_out)
        print("            실패한 건은 그 조각 폴더에서 같은 명령을 다시 돌리면 된다.")
        print("            미수렴 건은 결과 파일을 지우고 --max-orders 를 올려 다시 돌린다.")

    # --- 결과 파일 복사 -----------------------------------------------------
    if args.copy:
        dst_dir = os.path.join(out_dir, "out", "ocn")
        os.makedirs(dst_dir, exist_ok=True)
        n_cp = 0
        for r in ok:
            p = r.get("out_path", "")
            if p and os.path.isfile(p):
                dst = os.path.join(dst_dir, os.path.basename(p))
                if not os.path.isfile(dst):
                    shutil.copy2(p, dst)
                r["out_path"] = dst
                n_cp += 1
        # 경로가 바뀌었으므로 기록을 다시 쓴다
        with open(idx_out, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in merged:
                w.writerow(r)
        print("결과 파일 %d 개를 %s 로 모았다." % (n_cp, dst_dir))

    print("")
    if not problems and not missing_file and short == 0:
        print("빠짐도 문제도 없다. 분석으로 넘어가면 된다.")
    else:
        print("위의 문제를 정리한 뒤 분석으로 넘어가는 것이 안전하다.")


if __name__ == "__main__":
    main()
