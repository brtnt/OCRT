#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OCRT 편광 민감도 실행기 (별 모양 설계)

두 가지로 뻗는 실행 목록을 만들어 결합 격자 계산을 돌린다.
  N 가지: 대기를 중심에 고정하고 해수 27종(865 nm는 9종)을 훑는다.
  S 가지: 해수를 중심에 고정하고 대기 15종(3모델 x 광학두께 5단)을 훑는다.
두 가지가 공유하는 중심 실행은 한 번만 돌린다.

특징
  - 한 행이 끝날 때마다 진행 숫자와 남은 시간이 갱신된다.
  - 이미 끝난 행은 건너뛴다(재개 가능).
  - 모든 실행 기록을 index.csv 에 한 행씩 덧붙인다. 나중에 광학두께를
    더해도 기존 파일을 덮지 않으며 행만 늘어난다.
"""

import argparse
import csv
import hashlib
import itertools
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta

# 윈도우 콘솔에서 한글이 깨지지 않도록 한다.
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
# 1. 실험 설계 상수 — 여기만 고치면 설계가 바뀐다
# ===========================================================================

BANDS_NM = [412, 555, 667, 748, 865]        # GOCI-III 편광 밴드
VZA_STEP_DEG = 5                             # 관측천정각 간격 (0~85도)
RAA_STEP_DEG = 5                             # 방위각 간격 (0~355도)

# 엽록소 파장 상한은 2026-07-28 판본에서 1100 nm 로 넓어졌다. 865 nm 특례가 없다.
CHL_MAX_WL_NM = 1100.0

# 중심 해수와 중심 대기. 별 모양 설계의 한가운데다.
CENTER_WATER = dict(chl=1.0, tsm=1.5, cdom=0.035)
CENTER_AER = "r50f05v01"
CENTER_AOD = 0.10
CENTER_WIND = 5.0

AER_DIR = os.path.join("inputs", "aerosol_ahmad2010_accurt_mie")

# 에어로졸 목록. 443 nm 단일산란알베도를 함께 적는다.
AER_3 = [
    ("r30f95v01", 0.9503),
    ("r50f05v01", 0.9841),
    ("r95f00v01", 1.0000),
]
AER_8 = [
    ("r30f95v01", 0.9503),   # 비대칭443=0.649, SSA865=0.9216   1단계부터 사용
    ("r70f80v01", 0.9568),   # 비대칭443=0.669, SSA865=0.9360
    ("r75f50v01", 0.9636),   # 비대칭443=0.692, SSA865=0.9541
    ("r70f20v01", 0.9679),   # 비대칭443=0.698, SSA865=0.9741
    ("r50f10v01", 0.9754),   # 비대칭443=0.702, SSA865=0.9860
    ("r50f05v01", 0.9841),   # 비대칭443=0.719, SSA865=0.9927   1단계부터 사용
    ("r80f02v01", 0.9937),   # 비대칭443=0.790, SSA865=0.9972
    ("r95f00v01", 1.0000),   # 비대칭443=0.836, SSA865=1.0000   1단계부터 사용
]

# 축퇴 격자. 미세입자 비율(크기 축) 4수준 x 상대습도(흡수 축) 4수준.
# 같은 행(미세입자 비율 고정)에서는 세기 분광 기울기가 거의 같은데
# 단일산란알베도만 크게 변한다. 세기만으로는 가려낼 수 없는 상황을
# 일부러 만들어 편광이 이를 푸는지 본다.
AER_DEG = [
    ("r30f80v01", 0.9513),   # 습도30 미세80%, 옹스트롬 2.39, 알베도865 0.9288
    ("r50f80v01", 0.9527),   # 습도50 미세80%, 옹스트롬 2.38, 알베도865 0.9308
    ("r30f50v01", 0.9549),   # 습도30 미세50%, 옹스트롬 1.82, 알베도865 0.9471
    ("r50f50v01", 0.9562),   # 습도50 미세50%, 옹스트롬 1.81, 알베도865 0.9485
    ("r30f30v01", 0.9602),   # 습도30 미세30%, 옹스트롬 1.25, 알베도865 0.9636
    ("r50f30v01", 0.9612),   # 습도50 미세30%, 옹스트롬 1.24, 알베도865 0.9645
    ("r30f20v01", 0.9653),   # 습도30 미세20%, 옹스트롬 0.88, 알베도865 0.9738
    ("r50f20v01", 0.9661),   # 습도50 미세20%, 옹스트롬 0.87, 알베도865 0.9744
    ("r80f80v01", 0.9709),   # 습도80 미세80%, 옹스트롬 2.28, 알베도865 0.9574
    ("r80f50v01", 0.9724),   # 습도80 미세50%, 옹스트롬 1.91, 알베도865 0.9651
    ("r80f30v01", 0.9748),   # 습도80 미세30%, 옹스트롬 1.46, 알베도865 0.9735
    ("r80f20v01", 0.9773),   # 습도80 미세20%, 옹스트롬 1.12, 알베도865 0.9796
    ("r95f80v01", 0.9858),   # 습도95 미세80%, 옹스트롬 2.06, 알베도865 0.9812
    ("r95f50v01", 0.9863),   # 습도95 미세50%, 옹스트롬 1.82, 알베도865 0.9835
    ("r95f30v01", 0.9873),   # 습도95 미세30%, 옹스트롬 1.51, 알베도865 0.9865
    ("r95f20v01", 0.9883),   # 습도95 미세20%, 옹스트롬 1.24, 알베도865 0.9890
]
# 기존 8종을 버리지 않고 합친다. 이미 돌린 실행이 그대로 살아난다.
AER_24 = sorted(set(AER_8) | set(AER_DEG), key=lambda x: x[1])

# 해수 수준. 4수준은 3수준에 한 점씩 더한 것이라 기존 실행을 그대로 재사용한다.
CHL_3 = [0.15, 1.0, 4.0]
CHL_4 = [0.15, 0.4, 1.0, 4.0]
TSM_3 = [0.15, 1.5, 15.0]
TSM_4 = [0.15, 1.5, 5.0, 15.0]
CDOM_3 = [0.012, 0.035, 0.09]
CDOM_4 = [0.012, 0.02, 0.035, 0.09]

S_AOD_5 = [0.05, 0.10, 0.20, 0.30, 0.40]           # S 가지 광학두께 (1차 실행)
S_AOD_7 = [0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.60]
# 0.02 와 0.60 을 더한 이유. S 를 재려면 모델들이 같은 세기를 내는 구간이
# 있어야 하는데, 5단에서는 그 구간이 전체의 8.5% 뿐이었다. 양끝을 넓히면
# 30.9% 로 늘어난다. 555 nm, 태양천정각 50도에서 실측해 확인한 값이다.

# ---------------------------------------------------------------------------
# 실행 단계. 앞 단계에서 만든 파일은 그대로 재사용되므로 새것만 계산한다.
#
#   base : 1차 실행 (2026-08-02 완료분)
#   p4   : 에어로졸 3->8, 해수 27->64, N 가지 광학두께 3단
#   p5   : 풍속 3단, 태양천정각 6단
# ---------------------------------------------------------------------------
STAGES = {
    "base": dict(
        aer=AER_3, chl=CHL_3, tsm=TSM_3, cdom=CDOM_3,
        n_aod=[0.10], s_aod=S_AOD_5, sza=[25, 50, 75], wind=[5.0],
        label="1차 실행"),
    "p4a": dict(
        aer=AER_8, chl=CHL_3, tsm=TSM_3, cdom=CDOM_3,
        n_aod=[0.10], s_aod=S_AOD_5, sza=[25, 50, 75], wind=[5.0],
        label="4-1 에어로졸 3종 -> 8종"),
    "p4b": dict(
        aer=AER_8, chl=CHL_4, tsm=TSM_4, cdom=CDOM_4,
        n_aod=[0.10], s_aod=S_AOD_5, sza=[25, 50, 75], wind=[5.0],
        label="4-2 해수 27종 -> 64종"),
    "p4": dict(
        aer=AER_8, chl=CHL_4, tsm=TSM_4, cdom=CDOM_4,
        n_aod=[0.05, 0.10, 0.30], s_aod=S_AOD_5, sza=[25, 50, 75], wind=[5.0],
        label="4-3 N 가지 광학두께 3단 (4단계 전체)"),
    "p4c": dict(
        aer=AER_8, chl=CHL_4, tsm=TSM_4, cdom=CDOM_4,
        n_aod=[0.05, 0.10, 0.30], s_aod=S_AOD_7, sza=[25, 50, 75], wind=[5.0],
        label="4-4 S 가지 광학두께 0.02 와 0.60 추가"),
    "p6": dict(
        aer=AER_24, chl=CHL_4, tsm=TSM_4, cdom=CDOM_4,
        n_aod=[0.05, 0.10, 0.30], s_aod=S_AOD_7, sza=[25, 50, 75], wind=[5.0],
        label="6-1 축퇴 격자 16종 추가 (에어로졸 24종)"),
    "p6a": dict(
        aer=AER_24, chl=CHL_4, tsm=TSM_4, cdom=CDOM_4,
        n_aod=[0.05, 0.10, 0.30], s_aod=S_AOD_7, sza=[25, 50, 75], wind=[5.0],
        surface="black",
        label="6-2 해수 없는 기준 대기 (흑색 프레넬)"),
    "p5a": dict(
        aer=AER_8, chl=CHL_3, tsm=TSM_3, cdom=CDOM_3,
        n_aod=[0.10], s_aod=S_AOD_7, sza=[25, 50, 75], wind=[3.0, 5.0, 9.0],
        label="5-1 풍속 3, 5, 9 m/s"),
    "p5b": dict(
        aer=AER_8, chl=CHL_3, tsm=TSM_3, cdom=CDOM_3,
        n_aod=[0.10], s_aod=S_AOD_7,
        sza=[15, 25, 35, 50, 75], wind=[3.0, 5.0, 9.0],
        label="5-2 태양천정각 15도와 35도 추가"),
    "p5": dict(
        aer=AER_8, chl=CHL_3, tsm=TSM_3, cdom=CDOM_3,
        n_aod=[0.10], s_aod=S_AOD_7,
        sza=[15, 25, 35, 50, 65, 75], wind=[3.0, 5.0, 9.0],
        label="5-3 태양천정각 65도 추가 (5단계 전체)"),
}

NOMINAL_SEC = 220.0    # 1차 실행 615건 실측 평균 (1코어)


# ===========================================================================
# 2. 실행 목록 만들기
# ===========================================================================

def fmt_num(x):
    """파일 이름에 쓰기 위해 소수점을 p 로 바꾼다."""
    s = ("%g" % x)
    return s.replace(".", "p").replace("-", "m")


def make_case(branch, band, sza, chl, tsm, cdom, model, ssa, aod, wind,
              surface="ocean"):
    """실행 한 건을 사전으로 만든다."""
    if band > CHL_MAX_WL_NM:
        chl = 0.0
    if surface == "black":
        # 해수가 없으므로 이름에 해수 값을 넣지 않는다.
        chl = tsm = cdom = 0.0
        name = ("A_b%d_s%d_m%s_a%s" % (band, sza, model, fmt_num(aod)))
    else:
        name = ("O_b%d_s%d_chl%s_tsm%s_cdm%s_m%s_a%s"
                % (band, sza, fmt_num(chl), fmt_num(tsm), fmt_num(cdom),
                   model, fmt_num(aod)))
    # 풍속은 1차 실행 때 파일 이름에 없었다. 기본값이 아닐 때만 뒤에 붙여
    # 이미 만들어 둔 파일을 그대로 쓸 수 있게 한다.
    if abs(wind - CENTER_WIND) > 1e-9:
        name += "_w" + fmt_num(wind)
    return dict(branch=branch, band_nm=band, sza_deg=sza,
                chl=chl, tsm=tsm, cdom=cdom,
                aer_model=model, ssa443=ssa, aod865=aod,
                wind_ms=wind, surface=surface, name=name)


def build_runlist(stage):
    """단계 정의에서 N 가지와 S 가지를 만들고 중복을 없앤다."""
    st = STAGES[stage]
    center_ssa = dict(st["aer"])[CENTER_AER]
    cases = []

    if st.get("surface") == "black":
        # 해수 없는 기준 대기. 해수를 훑을 필요가 없으므로 대기 조합만 만든다.
        for band, sza, wind in itertools.product(
                BANDS_NM, st["sza"], st["wind"]):
            for (model, ssa), aod in itertools.product(st["aer"], st["s_aod"]):
                cases.append(make_case("A", band, sza, 0, 0, 0,
                                       model, ssa, aod, wind, surface="black"))
        uniq = {}
        for c in cases:
            uniq.setdefault(c["name"], c)
        out = sorted(uniq.values(),
                     key=lambda c: (c["band_nm"], c["sza_deg"], c["name"]))
        for i, c in enumerate(out, 1):
            c["run_id"] = i
        return out

    # --- N 가지: 해수를 훑는다 (대기는 중심 고정) -------------------------
    for band, sza, wind, aod in itertools.product(
            BANDS_NM, st["sza"], st["wind"], st["n_aod"]):
        for chl, tsm, cdom in itertools.product(st["chl"], st["tsm"], st["cdom"]):
            cases.append(make_case("N", band, sza, chl, tsm, cdom,
                                   CENTER_AER, center_ssa, aod, wind))

    # --- S 가지: 대기를 훑는다 (해수는 중심 고정) -------------------------
    for band, sza, wind in itertools.product(BANDS_NM, st["sza"], st["wind"]):
        for (model, ssa), aod in itertools.product(st["aer"], st["s_aod"]):
            cases.append(make_case("S", band, sza,
                                   CENTER_WATER["chl"], CENTER_WATER["tsm"],
                                   CENTER_WATER["cdom"], model, ssa, aod, wind))

    # --- 중복 제거 (중심 실행은 두 가지에 모두 들어 있다) -----------------
    seen, uniq = {}, []
    for c in cases:
        if c["name"] in seen:
            seen[c["name"]]["branch"] = "NS"
            continue
        seen[c["name"]] = c
        uniq.append(c)

    uniq.sort(key=lambda c: (c["band_nm"], c["sza_deg"], c["name"]))
    for i, c in enumerate(uniq, 1):
        c["run_id"] = i
    return uniq


# ===========================================================================
# 3. 명령 만들기와 결과 점검
# ===========================================================================

def build_cmd(exe, case, out_tmp, max_orders, vza_step, raa_step):
    if case.get("surface") == "black":
        mie = os.path.join(AER_DIR, case["aer_model"] + ".mie")
        cmd = [
            exe,
            "--sza", "%g" % case["sza_deg"],
            "--wavelength", "%g" % case["band_nm"],
            "--surface", "black_fresnel_ocean",
            "--wind-speed", "%g" % case["wind_ms"],
            "--decouple-sunglint",
            "--pssa",
            "--aod-865", "%g" % case["aod865"],
            "--mie", mie,
            "--lut-vza-step", "%g" % vza_step,
            "--lut-raa-step", "%g" % raa_step,
            "--output-full-grid", out_tmp,
        ]
        return cmd
    mie = os.path.join(AER_DIR, case["aer_model"] + ".mie")
    cmd = [
        exe,
        "--sza", "%g" % case["sza_deg"],
        "--wavelength", "%g" % case["band_nm"],
        "--surface", "ocean",
        "--wind-speed", "%g" % case["wind_ms"],
        "--decouple-sunglint",
        "--pssa",
        "--aod-865", "%g" % case["aod865"],
        "--mie", mie,
        "--water-model", "ocrt",
        "--ocrt-chl", "%g" % case["chl"],
        "--ocrt-tsm", "%g" % case["tsm"],
        "--ocrt-adom440", "%g" % case["cdom"],
        "--lut-vza-step", "%g" % vza_step,
        "--lut-raa-step", "%g" % raa_step,
        "--output-full-grid", out_tmp,
    ]
    if max_orders:
        cmd += ["--max-orders", str(max_orders)]
    return cmd


def inspect_output(path):
    """출력 파일의 행 수와 수렴 여부를 살핀다."""
    n_rows, n_bad, max_orders = 0, 0, 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
            rd = csv.DictReader(f)
            for row in rd:
                n_rows += 1
                try:
                    if int(float(row.get("water_converged", "1"))) == 0:
                        n_bad += 1
                    max_orders = max(max_orders, int(float(row.get("water_orders", "0"))))
                except (TypeError, ValueError):
                    pass
    except OSError:
        return 0, 0, 0
    return n_rows, n_bad, max_orders


def sha256_of(path):
    h = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return ""
    return h.hexdigest()


# ===========================================================================
# 4. 진행 표시
# ===========================================================================

def hms(sec):
    if sec is None or sec != sec or sec < 0:
        return "--:--:--"
    sec = int(sec)
    return "%02d:%02d:%02d" % (sec // 3600, (sec % 3600) // 60, sec % 60)


class Progress:
    """한 줄짜리 진행 표시. 행이 끝날 때마다 숫자가 갱신된다."""

    def __init__(self, total, skipped, workers):
        self.total = total
        self.skipped = skipped
        self.workers = workers
        self.done = 0
        self.failed = 0
        self.unconv = 0
        self.durations = []
        self.t0 = time.time()
        self.width = 0

    def avg(self):
        if self.durations:
            return sum(self.durations) / len(self.durations)
        return NOMINAL_SEC

    def eta(self, running_elapsed):
        left = self.total - self.done
        if left <= 0:
            return 0.0
        work = left * self.avg() - sum(running_elapsed)
        if work < 0:
            work = 0.0
        return work / max(1, min(self.workers, left))

    def draw(self, running_elapsed):
        """한 줄로 진행 상황을 갱신한다. 콘솔 80칸 안에 들어가도록 압축했다."""
        done_all = self.done + self.skipped
        total_all = self.total + self.skipped
        pct = 100.0 * done_all / total_all if total_all else 100.0
        nfill = int(round(pct / 10.0))
        bar = "#" * nfill + "." * (10 - nfill)
        eta = self.eta(running_elapsed)
        finish = datetime.now() + timedelta(seconds=eta)
        line = ("[%s] %d/%d %5.1f%% | 남은 %s | 완료 %s | 평균 %.1f분 | 동시 %d"
                % (bar, done_all, total_all, pct, hms(eta),
                   finish.strftime("%H:%M"), self.avg() / 60.0,
                   len(running_elapsed)))
        if self.failed:
            line += " | 실패 %d" % self.failed
        if self.unconv:
            line += " | 미수렴 %d" % self.unconv
        pad = " " * max(0, self.width - len(line))
        self.width = len(line)
        sys.stdout.write("\r" + line + pad)
        sys.stdout.flush()

    def finish_line(self):
        sys.stdout.write("\n")
        sys.stdout.flush()


# ===========================================================================
# 5. 본체
# ===========================================================================

# 실행 파일이 제대로 동작하는지 확인하는 기준 사례.
# 판본 OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic 에서 낸 값이다.
# 홀수 방위각 모드 부호 수정으로 결합 사례의 U 성분이 이전 판본과 다르다.
VERIFY_CASES = [
    ("대기만 (흑색 프레넬, 555 nm)",
     ["--sza", "25", "--vza", "30", "--raa", "90", "--wavelength", "555",
      "--surface", "black_fresnel_ocean", "--wind-speed", "5",
      "--decouple-sunglint", "--pssa", "--aod-865", "0.1"],
     (4.4991861596e-02, -9.0983939253e-04, -7.5031808484e-03)),
    ("해수 결합 (중심 조건, 555 nm)",
     ["--sza", "25", "--vza", "30", "--raa", "90", "--wavelength", "555",
      "--surface", "ocean", "--wind-speed", "5",
      "--decouple-sunglint", "--pssa", "--aod-865", "0.1",
      "--water-model", "ocrt", "--ocrt-chl", "1.0",
      "--ocrt-tsm", "1.5", "--ocrt-adom440", "0.035"],
     (1.2305024573e-01, -6.3070476021e-04, -1.1478355407e-02)),
]


# 판본 OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic 의 주요 자료 파일 해시.
# 실행 파일만 갈아 끼우고 자료를 옛것으로 두면 결과가 조용히 틀어지므로 함께 대조한다.
DATA_FILES = [
    ("inputs/water_iop/water_coef_z09_1nm.txt",
     "f173e1a4514b653973def109a5415678afbe4f39564a7b19f94a71ef55ea0e46",
     "순수해수 흡수·산란표"),
]


def check_data_files(root):
    """주요 자료 파일이 이 판본의 것인지 확인한다."""
    ok = True
    print("자료 파일 확인")
    for rel, want, label in DATA_FILES:
        path = os.path.join(root, rel.replace("/", os.sep))
        if not os.path.isfile(path):
            print("  [없음] %s" % rel)
            print("         %s 가 있어야 한다." % label)
            ok = False
            continue
        got = sha256_of(path)
        if got == want:
            print("  [통과] %s" % rel)
        else:
            print("  [불일치] %s" % rel)
            print("           %s 가 이 판본의 것이 아니다." % label)
            print("           나온 해시 : %s" % got[:32])
            print("           기준 해시 : %s" % want[:32])
            ok = False
    if not ok:
        print("")
        print("  자료가 옛 판본이면 실행 파일이 맞아도 결과가 틀어진다.")
        print("  최신 패키지의 자료로 바꾼 뒤 다시 확인한다.")
    print("")
    return ok


def run_verify(exe, root):
    """기준 사례를 돌려 결과를 대조한다."""
    mie = os.path.join(AER_DIR, CENTER_AER + ".mie")
    ok_all = check_data_files(root)
    print("실행 파일 확인")
    print("  대상: %s" % exe)
    print("")
    for name, extra, expect in VERIFY_CASES:
        cmd = [exe] + extra + ["--mie", mie]
        t0 = time.time()
        try:
            r = subprocess.run(cmd, cwd=root, capture_output=True, text=True,
                               encoding="utf-8", errors="replace", timeout=1800)
        except Exception as e:
            print("  [실패] %s : 실행 자체가 되지 않는다 (%s)" % (name, e))
            ok_all = False
            continue
        dur = time.time() - t0
        got = None
        for ln in reversed((r.stdout or "").splitlines()):
            t = ln.split()
            if len(t) >= 3:
                try:
                    got = (float(t[0]), float(t[1]), float(t[2]))
                    break
                except ValueError:
                    continue
        if got is None:
            print("  [실패] %s : 결과를 읽지 못했다. 종료코드 %s" % (name, r.returncode))
            err = (r.stderr or "").strip().splitlines()
            if err:
                print("         %s" % err[-1][:200])
            ok_all = False
            continue
        worst = 0.0
        for g, e in zip(got, expect):
            denom = abs(e) if abs(e) > 0 else 1.0
            worst = max(worst, abs(g - e) / denom)
        mark = "통과" if worst < 1e-6 else "불일치"
        if worst >= 1e-6:
            ok_all = False
        print("  [%s] %s  (%.1f초, 최대 상대차 %.2e)" % (mark, name, dur, worst))
        print("         나온 값 : %.10e %.10e %.10e" % got)
        print("         기준 값 : %.10e %.10e %.10e" % expect)
    print("")
    if ok_all:
        print("두 사례 모두 기준값과 일치한다. 본 실행을 진행해도 된다.")
    else:
        print("일치하지 않는다. 본 실행을 진행하면 안 된다.")
        print("실행 파일이나 inputs 폴더 구성을 먼저 확인해야 한다.")
    return 0 if ok_all else 1


def find_root(explicit):
    """패키지 최상위 폴더(inputs 와 src 가 있는 곳)를 찾는다."""
    if explicit:
        return os.path.abspath(explicit)
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [here,
                  os.path.dirname(here),
                  os.path.dirname(os.path.dirname(here)),
                  os.getcwd(),
                  os.path.dirname(os.getcwd())]
    for c in candidates:
        if c and os.path.isdir(os.path.join(c, "inputs")) \
              and os.path.isdir(os.path.join(c, "src")):
            return c
    return None


def find_exe(explicit, root):
    """OCRT 실행 파일을 찾는다."""
    if explicit:
        return os.path.abspath(explicit)
    here = os.path.dirname(os.path.abspath(__file__))
    if os.name == "nt":
        names = ["ocrt_v1.2.exe", "ocrt.exe", "ocrt_v1.2", "ocrt"]
    else:
        names = ["ocrt_v1.2", "ocrt", "ocrt_v1.2.exe", "ocrt.exe"]
    places = [os.path.join(root, "build"), here, root]
    for d in places:
        for n in names:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return p
    return None


INDEX_FIELDS = ["run_id", "branch", "surface", "band_nm", "sza_deg",
                "chl", "tsm", "cdom",
                "aer_model", "ssa443", "aod865", "wind_ms", "pssa", "glint",
                "vza_step_deg", "raa_step_deg", "out_path", "status", "wall_s", "n_rows",
                "n_unconverged", "max_water_orders", "code_version",
                "binary_sha256", "run_tag", "finished_at"]


def main():
    ap = argparse.ArgumentParser(description="OCRT 편광 민감도 배치 실행기")
    ap.add_argument("--root", default=None,
                    help="OCRT 패키지 최상위 폴더 (생략하면 스스로 찾는다)")
    ap.add_argument("--exe", default=None,
                    help="실행 파일 경로 (기본: <root>/build/ocrt_v1.2.exe)")
    ap.add_argument("--out", default=None,
                    help="산출물 폴더 (생략하면 <root>/star_run)")
    ap.add_argument("--workers", type=int, default=24, help="동시 실행 수 (기본 24)")
    ap.add_argument("--stage", default="base", choices=sorted(STAGES.keys()),
                    help="실행 단계. 앞 단계 결과는 그대로 재사용된다")
    ap.add_argument("--max-orders", type=int, default=0,
                    help="수중 산란차수 상한 (0 이면 코드 기본값)")
    ap.add_argument("--tag", default=None, help="이번 실행 묶음 이름")
    ap.add_argument("--vza-step", type=float, default=VZA_STEP_DEG,
                    help="관측천정각 간격(도). 기본 %g" % VZA_STEP_DEG)
    ap.add_argument("--raa-step", type=float, default=RAA_STEP_DEG,
                    help="방위각 간격(도). 기본 %g" % RAA_STEP_DEG)
    ap.add_argument("--part", default=None, metavar="K/N",
                    help="실행 목록을 N 조각으로 나눠 그중 K 번째만 돌린다 (예: 1/3)")
    ap.add_argument("--limit", type=int, default=0,
                    help="이번에 돌릴 건수를 이 값으로 제한한다 (속도 측정용)")
    ap.add_argument("--verify", action="store_true",
                    help="기준 사례 두 건을 돌려 실행 파일이 제대로 동작하는지 확인한다")
    ap.add_argument("--dry-run", action="store_true",
                    help="실행하지 않고 목록과 예상 시간만 낸다")
    args = ap.parse_args()

    root = find_root(args.root)
    if not root:
        print("오류: OCRT 패키지 최상위 폴더를 찾지 못했다.")
        print("      inputs 폴더와 src 폴더가 함께 있는 곳이 최상위다.")
        print("      이 스크립트를 그 폴더나 그 아래 build 폴더에 두거나,")
        print("      --root 로 경로를 직접 지정한다.")
        print("      지금 스크립트 위치: %s" % os.path.dirname(os.path.abspath(__file__)))
        sys.exit(2)

    exe = find_exe(args.exe, root)
    if not exe:
        print("오류: OCRT 실행 파일을 찾지 못했다.")
        print("      %s 에 ocrt_v1.2.exe 를 두거나 --exe 로 지정한다."
              % os.path.join(root, "build"))
        sys.exit(2)

    out_dir = os.path.abspath(args.out or os.path.join(root, "star_run"))
    csv_dir = os.path.join(out_dir, "out", "ocn")
    tmp_dir = os.path.join(csv_dir, "_tmp")
    st = STAGES[args.stage]
    tag = args.tag or (datetime.now().strftime("%Y-%m-%d") + "_" + args.stage)
    if args.part:
        tag += "_part" + args.part.replace("/", "of")

    if args.verify:
        sys.exit(run_verify(exe, root))

    os.makedirs(csv_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    runlist = build_runlist(args.stage)

    # 실행 목록을 파일로 남긴다.
    rl_path = os.path.join(out_dir, "runlist_%s.csv" % tag)
    with open(rl_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["run_id", "branch", "band_nm", "sza_deg",
                                          "chl", "tsm", "cdom", "aer_model",
                                          "ssa443", "aod865", "wind_ms", "name"])
        w.writeheader()
        for c in runlist:
            w.writerow({k: c[k] for k in w.fieldnames})

    # 이미 끝난 행을 가려낸다.
    todo, skipped = [], 0
    for c in runlist:
        c["out_path"] = os.path.join(csv_dir, c["name"] + ".csv")
        if os.path.isfile(c["out_path"]) and os.path.getsize(c["out_path"]) > 0:
            skipped += 1
        else:
            todo.append(c)

    part_k = part_n = 0
    if args.part:
        try:
            a, b = args.part.split("/")
            part_k, part_n = int(a), int(b)
        except Exception:
            sys.exit("오류: --part 는 1/3 처럼 적는다. 받은 값: %s" % args.part)
        if part_n < 1 or part_k < 1 or part_k > part_n:
            sys.exit("오류: --part 값이 올바르지 않다. 받은 값: %s" % args.part)
        # 조각을 나누는 기준은 전체 실행 목록이다. 이미 끝난 것과 무관하게
        # 같은 조각은 언제 돌려도 같은 실행을 맡는다.
        mine = {c["name"] for c in runlist if (c["run_id"] - 1) % part_n == part_k - 1}
        todo = [c for c in todo if c["name"] in mine]
        skipped = sum(1 for c in runlist
                      if c["name"] in mine and c["name"] not in {t["name"] for t in todo})

    if args.limit and len(todo) > args.limit:
        todo = todo[:args.limit]
        limited = True
    else:
        limited = False

    n_branch = sum(1 for c in runlist if "N" in c["branch"])
    s_branch = sum(1 for c in runlist if "S" in c["branch"])
    print("실행 계획")
    print("  패키지 폴더    : %s" % root)
    print("  산출물 폴더    : %s" % out_dir)
    print("  실행 묶음 이름 : %s" % tag)
    print("  전체 실행 수   : %d 건 (N 가지 %d, S 가지 %d, 중심 공유분 제외)"
          % (len(runlist), n_branch, s_branch))
    if part_n:
        n_mine = sum(1 for c in runlist if (c["run_id"] - 1) % part_n == part_k - 1)
        print("  이 조각         : %d / %d 번째, 맡은 실행 %d 건"
              % (part_k, part_n, n_mine))
    print("  이미 끝난 것   : %d 건" % skipped)
    print("  이번에 돌릴 것 : %d 건%s"
          % (len(todo), " (--limit 로 제한됨)" if limited else ""))
    print("  동시 실행 수   : %d" % args.workers)
    print("  단계           : %s (%s)" % (args.stage, st["label"]))
    print("  에어로졸       : %d종 | 해수 %d종 | N 광학두께 %s"
          % (len(st["aer"]), len(st["chl"]) * len(st["tsm"]) * len(st["cdom"]),
             ", ".join("%g" % a for a in st["n_aod"])))
    print("  태양천정각     : %s도 | 풍속 %s m/s"
          % (", ".join(str(x) for x in st["sza"]),
             ", ".join("%g" % w for w in st["wind"])))
    nv = int(85.0 / args.vza_step) + 1
    nr = int(round(360.0 / args.raa_step))
    print("  각도 격자      : 관측천정각 %g도 간격 %d개 x 방위각 %g도 간격 %d개 = %d점"
          % (args.vza_step, nv, args.raa_step, nr, nv * nr))
    est = len(todo) * NOMINAL_SEC / max(1, args.workers)
    print("  예상 소요      : 약 %s (1건 %g분 가정, 실측되면 갱신된다)"
          % (hms(est), NOMINAL_SEC / 60.0))
    print("  실행 목록 파일 : %s" % rl_path)

    if args.dry_run:
        print("\n--dry-run 이므로 실행하지 않는다.")
        return
    if not todo:
        print("\n돌릴 것이 없다. 모두 끝나 있다.")
        return

    if not check_data_files(root):
        sys.exit("자료 파일이 이 판본의 것이 아니다. 본 실행을 중단한다.")

    binary_hash = sha256_of(exe)
    code_version = ""
    try:
        r = subprocess.run([exe, "--help"], cwd=root, capture_output=True,
                           text=True, encoding="utf-8", errors="replace",
                           timeout=60)
        first = (r.stdout or "").splitlines()
        if first:
            code_version = first[0].strip()
    except Exception:
        pass

    # 실행 파일 사본을 남겨 나중 합치기에서 대조할 수 있게 한다.
    bin_dir = os.path.join(out_dir, "bin")
    os.makedirs(bin_dir, exist_ok=True)
    try:
        dst = os.path.join(bin_dir, os.path.basename(exe))
        if not os.path.isfile(dst):
            shutil.copy2(exe, dst)
        with open(os.path.join(bin_dir, "sha256.txt"), "a", encoding="utf-8") as f:
            f.write("%s  %s  %s\n" % (binary_hash, os.path.basename(exe), tag))
    except OSError:
        pass

    if args.stage != "base":
        base_idx = "index_%s" % args.stage
    else:
        base_idx = "index"
    if part_n:
        index_path = os.path.join(out_dir,
                                  "%s_part%dof%d.csv" % (base_idx, part_k, part_n))
    else:
        index_path = os.path.join(out_dir, base_idx + ".csv")
    new_index = not os.path.isfile(index_path)
    idx_f = open(index_path, "a", encoding="utf-8", newline="")
    idx = csv.DictWriter(idx_f, fieldnames=INDEX_FIELDS)
    if new_index:
        idx.writeheader()
        idx_f.flush()

    log_path = os.path.join(out_dir, "run_%s.log" % tag)
    log_f = open(log_path, "a", encoding="utf-8")

    print("  기록 파일      : %s" % index_path)
    print("  자세한 기록    : %s" % log_path)
    print("")

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"     # 한 건은 1코어, 병렬은 이 스크립트가 맡는다

    prog = Progress(total=len(todo), skipped=skipped, workers=args.workers)
    queue = list(todo)
    running = []          # (Popen, case, t0, tmp_path, devnull)
    interrupted = False

    try:
        while queue or running:
            # 빈 자리를 채운다.
            while queue and len(running) < args.workers and not interrupted:
                c = queue.pop(0)
                tmp = os.path.join(tmp_dir, c["name"] + ".csv")
                cmd = build_cmd(exe, c, tmp, args.max_orders,
                                args.vza_step, args.raa_step)
                dn = open(os.devnull, "w")
                p = subprocess.Popen(cmd, cwd=root, env=env,
                                     stdout=dn, stderr=subprocess.PIPE, text=True,
                                     encoding="utf-8", errors="replace")
                running.append([p, c, time.time(), tmp, dn])

            # 끝난 것을 거둔다.
            still = []
            for item in running:
                p, c, t0, tmp, dn = item
                if p.poll() is None:
                    still.append(item)
                    continue
                dur = time.time() - t0
                err = ""
                try:
                    raw = (p.stderr.read() or "")
                    err = " | ".join(ln.strip() for ln in raw.splitlines()
                                     if ln.strip() and not ln.lstrip().startswith("#"))
                except Exception:
                    pass
                dn.close()

                if p.returncode == 0 and os.path.isfile(tmp):
                    n_rows, n_bad, mo = inspect_output(tmp)
                    os.replace(tmp, c["out_path"])
                    status = "ok" if n_bad == 0 else "unconverged"
                    if n_bad:
                        prog.unconv += 1
                else:
                    n_rows, n_bad, mo = 0, 0, 0
                    status = "failed"
                    prog.failed += 1
                    if os.path.isfile(tmp):
                        try:
                            os.remove(tmp)
                        except OSError:
                            pass

                prog.done += 1
                prog.durations.append(dur)

                idx.writerow(dict(
                    run_id=c["run_id"], branch=c["branch"],
                    surface=c.get("surface", "ocean"), band_nm=c["band_nm"],
                    sza_deg=c["sza_deg"], chl=c["chl"], tsm=c["tsm"], cdom=c["cdom"],
                    aer_model=c["aer_model"], ssa443=c["ssa443"], aod865=c["aod865"],
                    wind_ms=c["wind_ms"], pssa=1, glint="decoupled",
                    vza_step_deg=args.vza_step, raa_step_deg=args.raa_step,
                    out_path=c["out_path"],
                    status=status, wall_s="%.1f" % dur, n_rows=n_rows,
                    n_unconverged=n_bad, max_water_orders=mo,
                    code_version=code_version, binary_sha256=binary_hash,
                    run_tag=tag,
                    finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
                idx_f.flush()

                log_f.write("[%s] %s %s %.1fs rows=%d unconv=%d orders=%d %s\n"
                            % (datetime.now().strftime("%H:%M:%S"), status,
                               c["name"], dur, n_rows, n_bad, mo, err[:600]))
                log_f.flush()

            running = still
            prog.draw([time.time() - it[2] for it in running])
            time.sleep(0.25)

    except KeyboardInterrupt:
        interrupted = True
        prog.finish_line()
        print("중단 신호를 받았다. 실행 중인 %d 건이 끝날 때까지 기다린다." % len(running))
        for p, c, t0, tmp, dn in running:
            try:
                p.wait()
            except Exception:
                pass
            dn.close()
            if os.path.isfile(tmp):
                try:
                    os.remove(tmp)
                except OSError:
                    pass

    prog.finish_line()
    idx_f.close()
    log_f.close()

    print("")
    print("끝났다.")
    print("  완료 %d 건, 실패 %d 건, 미수렴 %d 건" % (prog.done, prog.failed, prog.unconv))
    if prog.durations:
        print("  1건 평균 %.1f분, 최장 %.1f분"
              % (prog.avg() / 60.0, max(prog.durations) / 60.0))
    print("  총 경과 %s" % hms(time.time() - prog.t0))
    if prog.failed:
        print("  실패한 건은 %s 에서 사유를 확인한다. 같은 명령을 다시 돌리면"
              " 실패분만 재실행된다." % log_path)
    if prog.unconv:
        print("  미수렴 건이 있다. --max-orders 를 올려 다시 돌려야 한다."
              " 해당 출력 파일을 지운 뒤 재실행하면 된다.")


if __name__ == "__main__":
    main()
