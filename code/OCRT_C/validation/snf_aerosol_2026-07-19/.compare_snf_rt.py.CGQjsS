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

KV = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")
GAS_OFF = [
    "--gas-column-h2o", "0", "--gas-column-o3", "0", "--gas-column-no2", "0",
    "--gas-column-o2", "0", "--gas-column-co2", "0", "--gas-column-ch4", "0",
]


def run_one(job: tuple[str, Path, int, float, Path, Path, int]) -> tuple[str, int, float, str, dict[str, float]]:
    tag, mie, wl, pressure, exe, cwd, timeout_s = job
    cmd = [
        str(exe), "--surface", "black", "--sza", "30", "--vza", "20", "--raa", "90",
        "--wavelength", str(wl), "--pressure", str(pressure), "--mie", str(mie),
        "--aod-555", "0.15", "--sos-max-orders", "40", *GAS_OFF,
    ]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OCRT_ADVANCED"] = "1"
    p = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True,
                       timeout=timeout_s, check=False)
    if p.returncode != 0:
        raise RuntimeError(f"{tag} {mie} wl={wl} p={pressure}: rc={p.returncode}: {p.stderr}")
    tokens = p.stdout.split()
    if len(tokens) < 3:
        raise RuntimeError(f"malformed output for {tag} {mie}: {p.stdout!r}")
    values: dict[str, float] = {
        "rho_I": float(tokens[0]), "rho_Q": float(tokens[1]), "rho_U": float(tokens[2])
    }
    for key, value in KV.findall(p.stdout):
        try:
            values[key] = float(value)
        except ValueError:
            continue
    if int(values.get("conv", 0)) != 1:
        raise RuntimeError(f"not converged for {tag} {mie}: {p.stdout}")
    return tag, wl, pressure, mie.stem, values


def pct_delta(new: float, old: float) -> float:
    return math.nan if old == 0.0 else 100.0 * (new - old) / old


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="RT identity and legacy-replacement comparisons for SnF Mie caches.")
    ap.add_argument("--exe", type=Path, required=True)
    ap.add_argument("--canonical-dir", type=Path, required=True)
    ap.add_argument("--generated-dir", type=Path, required=True)
    ap.add_argument("--legacy-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--timeout", type=int, default=60)
    args = ap.parse_args()

    exe = args.exe.resolve()
    cwd = exe.parent.parent if exe.parent.name == "build" else exe.parent
    pressures = (0.0, 1013.25)
    wavelengths = (412, 550, 860)

    jobs: list[tuple[str, Path, int, float, Path, Path, int]] = []
    baseline_models = ("T50", "C50", "M80C", "O99")
    for model in baseline_models:
        for pressure in pressures:
            for wl in wavelengths:
                jobs.append((f"baseline:canonical:{model}", args.canonical_dir / f"{model}.mie", wl, pressure, exe, cwd, args.timeout))
                jobs.append((f"baseline:generated:{model}", args.generated_dir / f"{model}.mie", wl, pressure, exe, cwd, args.timeout))

    legacy_models = ("M50C", "M95C", "M98C")
    for model in legacy_models:
        for pressure in pressures:
            for wl in wavelengths:
                jobs.append((f"legacy:old:{model}", args.legacy_dir / f"{model}_legacy83.mie", wl, pressure, exe, cwd, args.timeout))
                jobs.append((f"legacy:new:{model}", args.canonical_dir / f"{model}.mie", wl, pressure, exe, cwd, args.timeout))

    results: dict[tuple[str, int, float], dict[str, float]] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [pool.submit(run_one, job) for job in jobs]
        for future in as_completed(futures):
            tag, wl, pressure, _stem, values = future.result()
            results[(tag, wl, pressure)] = values

    baseline_rows: list[dict[str, object]] = []
    for model in baseline_models:
        for pressure in pressures:
            for wl in wavelengths:
                a = results[(f"baseline:canonical:{model}", wl, pressure)]
                b = results[(f"baseline:generated:{model}", wl, pressure)]
                row: dict[str, object] = {
                    "model": model, "pressure_hpa": pressure, "wavelength_nm": wl,
                }
                exact = True
                for comp in ("I", "Q", "U"):
                    old = a[f"rho_{comp}"]
                    new = b[f"rho_{comp}"]
                    delta = new - old
                    exact = exact and new == old
                    row[f"canonical_rho_{comp}"] = old
                    row[f"generated_rho_{comp}"] = new
                    row[f"delta_rho_{comp}"] = delta
                row["exact_output"] = exact
                row["canonical_AOD_band"] = a.get("AOD_band")
                row["generated_AOD_band"] = b.get("AOD_band")
                row["delta_AOD_band"] = b.get("AOD_band", math.nan) - a.get("AOD_band", math.nan)
                baseline_rows.append(row)

    legacy_rows: list[dict[str, object]] = []
    for model in legacy_models:
        for pressure in pressures:
            for wl in wavelengths:
                old_values = results[(f"legacy:old:{model}", wl, pressure)]
                new_values = results[(f"legacy:new:{model}", wl, pressure)]
                row = {
                    "model": model, "pressure_hpa": pressure, "wavelength_nm": wl,
                    "legacy_AOD_band": old_values.get("AOD_band"),
                    "canonical_AOD_band": new_values.get("AOD_band"),
                    "AOD_band_delta_pct": pct_delta(new_values.get("AOD_band", math.nan), old_values.get("AOD_band", math.nan)),
                }
                for comp in ("I", "Q", "U"):
                    old = old_values[f"rho_{comp}"]
                    new = new_values[f"rho_{comp}"]
                    row[f"legacy_rho_{comp}"] = old
                    row[f"canonical_rho_{comp}"] = new
                    row[f"delta_rho_{comp}"] = new - old
                    row[f"delta_rho_{comp}_pct"] = pct_delta(new, old)
                legacy_rows.append(row)

    write_rows(args.out_dir / "baseline_rt_identity_24cases.csv", baseline_rows)
    write_rows(args.out_dir / "legacy83_vs_canonical_rt_18cases.csv", legacy_rows)

    exact_count = sum(bool(row["exact_output"]) for row in baseline_rows)
    max_baseline = max(abs(float(row[key])) for row in baseline_rows for key in ("delta_rho_I", "delta_rho_Q", "delta_rho_U"))
    max_legacy_i_pct = max(abs(float(row["delta_rho_I_pct"])) for row in legacy_rows if math.isfinite(float(row["delta_rho_I_pct"])))
    print(f"BASELINE_RT_IDENTITY exact={exact_count}/{len(baseline_rows)} max_abs_delta={max_baseline:.3e}")
    print(f"LEGACY_REPLACEMENT cases={len(legacy_rows)} max_abs_delta_I_pct={max_legacy_i_pct:.3f}%")


if __name__ == "__main__":
    main()
