#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

MODELS = (
    "T50", "T80", "T90", "T95",
    "C50", "C70", "C80", "C90", "C95",
    "M50C", "M70C", "M80C", "M90C", "M95C", "M98C",
    "O99",
)
WAVELENGTHS_NM = (412, 550, 860)
KV = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")
GAS_OFF = [
    "--gas-column-h2o", "0",
    "--gas-column-o3", "0",
    "--gas-column-no2", "0",
    "--gas-column-o2", "0",
    "--gas-column-co2", "0",
    "--gas-column-ch4", "0",
]


def infer_runtime_cwd(exe: Path) -> Path:
    exe = exe.resolve()
    if exe.parent.name == "build":
        return exe.parent.parent
    if exe.parent.name == "bin" and (exe.parent.parent / "ocrt").is_dir():
        return exe.parent.parent / "ocrt"
    return exe.parent


def run_case(job: tuple[Path, Path, int, Path, int]) -> dict[str, object]:
    exe, mie, wavelength_nm, cwd, timeout_s = job
    cmd = [
        str(exe),
        "--surface", "black",
        "--sza", "30",
        "--vza", "20",
        "--raa", "90",
        "--wavelength", str(wavelength_nm),
        "--pressure", "0",
        "--mie", str(mie),
        "--aod-555", "0.15",
        *GAS_OFF,
    ]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    row: dict[str, object] = {
        "model": mie.stem,
        "wavelength_nm": wavelength_nm,
    }
    try:
        p = subprocess.run(
            cmd,
            text=True,
            capture_output=True,
            env=env,
            timeout=timeout_s,
            cwd=cwd,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        row.update(returncode=124, pass_=False, stderr=f"TIMEOUT after {timeout_s}s: {exc}")
        return row

    row["returncode"] = p.returncode
    row["stderr"] = p.stderr.strip().replace("\n", " | ")
    if p.stdout.strip():
        toks = p.stdout.strip().split()
        try:
            row.update(rho_I=float(toks[0]), rho_Q=float(toks[1]), rho_U=float(toks[2]))
        except (ValueError, IndexError):
            pass
        for key, value in KV.findall(p.stdout):
            try:
                row[key] = float(value)
            except ValueError:
                row[key] = value
    expected_policy_rejection = (
        mie.stem == "M50C"
        and p.returncode == 2
        and "diagnostic-only" in p.stderr
        and "M50C_POLICY_v0_9" in p.stderr
    )
    passed = expected_policy_rejection or (
        p.returncode == 0
        and int(float(row.get("conv", 0))) == 1
        and all(math.isfinite(float(row.get(k, float("nan")))) for k in ("rho_I", "rho_Q", "rho_U"))
    )
    row["expected_policy_rejection"] = expected_policy_rejection
    row["pass"] = passed
    return row


def main() -> None:
    ap = argparse.ArgumentParser(
        description=("Parallel OCRT smoke test for the 15 production SnF "
                     "models plus the required M50C diagnostic-only rejection."))
    ap.add_argument("--exe", type=Path, required=True)
    ap.add_argument("--mie-dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--timeout", type=int, default=60, help="per-case timeout in seconds")
    args = ap.parse_args()

    exe = args.exe.resolve()
    mie_dir = args.mie_dir.resolve()
    cwd = infer_runtime_cwd(exe)
    missing = [name for name in MODELS if not (mie_dir / f"{name}.mie").is_file()]
    if missing:
        raise SystemExit(f"missing canonical Mie files: {', '.join(missing)}")

    jobs = [
        (exe, mie_dir / f"{model}.mie", wavelength_nm, cwd, args.timeout)
        for model in MODELS
        for wavelength_nm in WAVELENGTHS_NM
    ]
    rows: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [pool.submit(run_case, job) for job in jobs]
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda r: (str(r["model"]), int(r["wavelength_nm"])))

    preferred = [
        "model", "wavelength_nm", "returncode", "pass",
        "expected_policy_rejection",
        "rho_I", "rho_Q", "rho_U",
        "AOD_ref", "AOD_ref_nm", "AOD_band", "AOD_ext_ratio",
        "tau_R", "orders", "conv", "stderr",
    ]
    fields = [key for key in preferred if any(key in row for row in rows)]
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    passed = sum(bool(row.get("pass")) for row in rows)
    failed = len(rows) - passed
    print(f"SNF_SMOKE pass={passed} fail={failed} cases={len(rows)}")
    for row in rows:
        if not bool(row.get("pass")):
            print(row)
    raise SystemExit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
