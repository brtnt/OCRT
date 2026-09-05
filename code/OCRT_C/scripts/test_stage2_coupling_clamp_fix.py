#!/usr/bin/env python3
"""Regression for the water-to-air coupling interpolation clamp correction.

The test is intentionally local to the corrected module boundary.  It does not
re-open the validated atmosphere, FIX1--FIX4, exact-pole, or in-water physics.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
from typing import Dict, Iterable, List

EXPECTED_VERSION = "OCRT-v1.2-2026-08-02-KST-dtpsign-phase-diagnostic"
NODES = (48, 64, 96)


def run(cmd: List[str], cwd: Path, env: Dict[str, str], timeout: int = 360) -> subprocess.CompletedProcess[str]:
    p = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True, timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError(
            f"command failed rc={p.returncode}: {shlex.join(cmd)}\n"
            f"stdout:\n{p.stdout[-2000:]}\nstderr:\n{p.stderr[-4000:]}"
        )
    return p


def common_args() -> List[str]:
    return [
        "--sza", "40", "--wavelength", "443", "--surface", "ocean",
        "--wind-speed", "3", "--pressure", "1013.25",
        "--n-mu", "24", "--n-layers", "40", "--m-max", "2",
        "--sos-max-orders", "100", "--debug-water-max-orders", "100",
        "--n-water", "1.34", "--water-temperature", "20",
        "--water-salinity", "38.4", "--water-model", "ocrt",
        "--ocrt-chl", "0", "--ocrt-tsm", "0", "--ocrt-adom440", "0",
        "--ocrt-adom-slope", "0.014", "--water-m-max", "4",
        "--gas-column-h2o", "0", "--gas-column-o3", "0",
        "--gas-column-no2", "0", "--gas-column-o2", "0",
        "--gas-column-co2", "0", "--gas-column-ch4", "0",
        "--decouple-sunglint",
    ]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(row: Dict[str, str], key: str) -> float:
    return float(row[key])


def find_row(rows: Iterable[Dict[str, str]], vza: float, raa: float) -> Dict[str, str]:
    for row in rows:
        if abs(f(row, "vza_deg") - vza) < 1e-12 and abs(f(row, "raa_deg") - raa) < 1e-12:
            return row
    raise AssertionError(f"missing full-grid row vza={vza}, raa={raa}")


def rrs_component(row: Dict[str, str], component: str) -> float:
    if component == "I":
        return f(row, "rrs_I")
    return f(row, f"Lu0minus_{component}") / f(row, "Ed0minus_water")


def compare_reference(got: Path, reference: Path, atol: float = 5e-13, rtol: float = 5e-13) -> None:
    g = read_csv(got)
    r = read_csv(reference)
    if len(g) != len(r):
        raise AssertionError(f"row-count mismatch: {len(g)} != {len(r)}")
    if not g:
        raise AssertionError("empty regression CSV")
    reference_keys = list(r[0])
    missing = [key for key in reference_keys if key not in g[0]]
    if missing:
        raise AssertionError(f"CSV is missing reference columns: {missing}")
    for i, (ga, ra) in enumerate(zip(g, r), start=2):
        for key in reference_keys:
            try:
                x, y = float(ga[key]), float(ra[key])
            except ValueError:
                if ga[key] != ra[key]:
                    raise AssertionError(f"text mismatch row={i} column={key}: {ga[key]!r} != {ra[key]!r}")
                continue
            if not math.isclose(x, y, rel_tol=rtol, abs_tol=atol):
                raise AssertionError(f"numeric mismatch row={i} column={key}: {x:.17g} != {y:.17g}")


def source_audit(root: Path) -> None:
    text = (root / "src/rt_air_water_coupling.c").read_text(encoding="utf-8")
    forbidden = ["int ord[64]", "double ms[64]", "if (n > 64) n = 64"]
    for needle in forbidden:
        if needle in text:
            raise AssertionError(f"legacy silent clamp remains in source: {needle}")
    required = ["AW_INTERP_STACK_N = 320", "failed to allocate workspace", "abort();"]
    for needle in required:
        if needle not in text:
            raise AssertionError(f"required fail-loud interpolation guard missing: {needle}")


def check_nadir_spin2(rows: List[Dict[str, str]], n_mu_water: int) -> None:
    r0 = find_row(rows, 0.0, 0.0)
    r45 = find_row(rows, 0.0, 45.0)
    r90 = find_row(rows, 0.0, 90.0)
    for label, getter in (
        ("Rrs_Q", lambda r: f(r, "Rrs_Q")),
        ("TOA_rho_Q", lambda r: f(r, "TOA_rho_Q")),
        ("rrs_Q", lambda r: rrs_component(r, "Q")),
    ):
        q0, q90 = getter(r0), getter(r90)
        m0 = 0.5 * (q0 + q90)
        scale = max(abs(q0), abs(q90), 1e-300)
        if abs(m0) > max(1e-14, 1e-12 * scale):
            raise AssertionError(f"n_mu_water={n_mu_water} forbidden nadir m=0 in {label}: {m0:.17g}")
    if abs(f(r45, "Rrs_Q")) > 1e-13:
        raise AssertionError(f"n_mu_water={n_mu_water} nadir Rrs_Q at RAA=45 is not zero")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--binary", type=Path, default=Path("build/ocrt"))
    ap.add_argument("--reference-dir", type=Path, default=Path("validation/stage2_coupling_clamp_fix_2026-08-02_dtpsign/reference"))
    args = ap.parse_args()

    root = args.root.resolve()
    binary = args.binary if args.binary.is_absolute() else root / args.binary
    refdir = args.reference_dir if args.reference_dir.is_absolute() else root / args.reference_dir
    env = dict(os.environ, OCRT_ADVANCED="1", OCRT_DEBUG="1", OMP_NUM_THREADS="1")

    source_audit(root)
    version = run([str(binary), "--version"], root, env, timeout=30).stdout.strip()
    if EXPECTED_VERSION not in version:
        raise AssertionError(f"unexpected version string: {version}")

    with tempfile.TemporaryDirectory(prefix="ocrt-clamp-reg-") as td:
        tmp = Path(td)
        generated: Dict[int, Path] = {}
        for n in NODES:
            print(f"RUN full-grid n_mu_water={n}", flush=True)
            path = tmp / f"fullgrid_n{n}.csv"
            cmd = [
                str(binary), *common_args(), "--n-mu-water", str(n),
                "--lut-vza-step", "30", "--lut-vza-max", "60",
                "--lut-raa-step", "45", "--output-full-grid", str(path),
            ]
            run(cmd, root, env)
            generated[n] = path
            rows = read_csv(path)
            if len(rows) != 24:
                raise AssertionError(f"n_mu_water={n}: expected 24 full-grid rows, got {len(rows)}")
            compare_reference(path, refdir / f"fullgrid_n{n}.csv")
            check_nadir_spin2(rows, n)
            print(f"PASS full-grid n_mu_water={n}", flush=True)


    print("PASS: coupling interpolation retains all water directions")
    print("PASS: n_mu_water=48/64/96 full-grid references")
    print("PASS: forbidden nadir m=0 polarization is zero")
    print("NOTE: single/full parity is checked by test_stage2_coupling_clamp_single_parity.sh")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise
