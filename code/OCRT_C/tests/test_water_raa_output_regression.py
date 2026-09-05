#!/usr/bin/env python3
"""Runtime regression for the stage-2 water RAA output fix.

The test is reference-code independent.  It checks that the reported rrs(0-)
and the TOA water-leaving increment vary with the same public OCRT RAA, that
principal-plane U vanishes, that full-grid and single-geometry outputs agree,
and that the separate wind=0 branch remains operational.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

GAS_OFF = [
    "--gas-column-h2o", "0", "--gas-column-o3", "0",
    "--gas-column-no2", "0", "--gas-column-o2", "0",
    "--gas-column-co2", "0", "--gas-column-ch4", "0",
]
RAAS = (0.0, 45.0, 90.0, 135.0, 180.0)


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env.update(
        OCRT_ADVANCED="1", OCRT_DEBUG="1", OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C",
    )
    cp = subprocess.run(
        cmd, cwd=cwd, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=180,
    )
    if cp.returncode != 0:
        raise RuntimeError(
            f"command failed rc={cp.returncode}\n"
            f"{' '.join(cmd)}\nSTDERR:\n{cp.stderr[-3000:]}"
        )
    return cp


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as fp:
        rows = []
        for row in csv.DictReader(fp):
            numeric: dict[str, float] = {}
            for key, value in row.items():
                if value in (None, ""):
                    continue
                try:
                    numeric[key] = float(value)
                except ValueError:
                    # Full-grid CSVs also carry diagnostic text columns
                    # (for example output_mode="ocean_rrs_grid").
                    continue
            rows.append(numeric)
    return rows


def by_geometry(rows: list[dict[str, float]], vza: float, raa: float) -> dict[str, float]:
    hit = [
        r for r in rows
        if abs(r["vza_deg"] - vza) < 1.0e-12
        and abs(r["raa_deg"] - raa) < 1.0e-12
    ]
    if len(hit) != 1:
        raise RuntimeError(f"geometry vza={vza}, raa={raa}: rows={len(hit)}")
    return hit[0]


def parse_pairs(stdout: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for token in stdout.replace("\n", " ").split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        try:
            out[key] = float(value)
        except ValueError:
            pass
    return out


def ocean_common(binary: Path, wind: float) -> list[str]:
    return [
        str(binary), "--surface", "ocean", "--decouple-sunglint", *GAS_OFF,
        "--wind-speed", f"{wind:g}",
        "--sza", "40", "--wavelength", "443", "--pressure", "1013.25",
        "--aod", "0", "--n-water", "1.34",
        "--n-mu", "48", "--n-layers", "40", "--m-max", "2",
        "--sos-max-orders", "100", "--debug-water-max-orders", "100",
        "--water-temperature", "20", "--water-salinity", "38.4",
        "--water-shared-grid", "--water-model", "ocrt",
        "--ocrt-chl", "0", "--ocrt-tsm", "0",
        "--ocrt-adom440", "0", "--ocrt-adom-slope", "0.014",
    ]


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    binary = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else root / "build/ocrt"
    if not binary.is_file():
        raise RuntimeError(f"binary not found: {binary}")

    with tempfile.TemporaryDirectory(prefix="ocrt_raa_reg_") as td:
        work = Path(td)
        ocean_csv = work / "ocean.csv"
        black_csv = work / "black.csv"
        wind0_csv = work / "wind0.csv"

        grid_args = [
            "--lut-vza-step", "30", "--lut-vza-max", "60",
            "--lut-raa-step", "45", "--output-full-grid",
        ]
        run(ocean_common(binary, 3.0) + grid_args + [str(ocean_csv)], root)

        black_cmd = [
            str(binary), "--surface", "black_fresnel_ocean",
            "--decouple-sunglint", *GAS_OFF,
            "--wind-speed", "3", "--sza", "40", "--wavelength", "443",
            "--pressure", "1013.25", "--aod", "0", "--n-water", "1.34",
            "--n-mu", "48", "--n-layers", "40", "--m-max", "2",
            "--sos-max-orders", "100",
            *grid_args, str(black_csv),
        ]
        run(black_cmd, root)

        ocean = read_rows(ocean_csv)
        black = read_rows(black_csv)
        ratios: list[float] = []
        for raa in RAAS:
            ow = by_geometry(ocean, 60.0, raa)
            bk = by_geometry(black, 60.0, raa)
            toa_water = ow["TOA_rho_I"] - bk["rho_I"]
            if toa_water <= 0.0:
                raise RuntimeError(f"non-positive TOA water increment at RAA={raa}: {toa_water}")
            ratios.append(ow["rrs_I"] / toa_water)

            if raa in (0.0, 180.0):
                rrs_u = ow["Lu0minus_U"] / ow["Ed0minus_water"]
                if abs(rrs_u) / max(abs(ow["rrs_I"]), 1.0e-30) > 1.0e-10:
                    raise RuntimeError(f"principal-plane rrs U nonzero at RAA={raa}: {rrs_u}")
                if abs(ow["Rrs_U"]) / max(abs(ow["Rrs_I"]), 1.0e-30) > 1.0e-10:
                    raise RuntimeError(f"principal-plane Rrs U nonzero at RAA={raa}: {ow['Rrs_U']}")

        spread = max(ratios) / min(ratios) - 1.0
        if spread > 0.08:
            raise RuntimeError(
                f"water/TOA azimuth consistency failed: spread={spread:.6f}; ratios={ratios}"
            )

        # Full-grid versus single-geometry parity at an off-principal azimuth.
        single = run(
            ocean_common(binary, 3.0)
            + ["--vza", "60", "--raa", "45"],
            root,
        )
        s = parse_pairs(single.stdout)
        g = by_geometry(ocean, 60.0, 45.0)
        pairs = {
            "Rrs0plus_I": g["Rrs_I"],
            "Rrs0plus_Q": g["Rrs_Q"],
            "Rrs0plus_U": g["Rrs_U"],
            "rrs0minus_I": g["rrs_I"],
            "rrs0minus_Q": g["Lu0minus_Q"] / g["Ed0minus_water"],
            "rrs0minus_U": g["Lu0minus_U"] / g["Ed0minus_water"],
        }
        for key, expected in pairs.items():
            if key not in s:
                raise RuntimeError(f"single output missing {key}")
            if abs(s[key] - expected) > 1.0e-7 * max(1.0, abs(expected)):
                raise RuntimeError(
                    f"full-grid/single mismatch {key}: single={s[key]:.12e}, "
                    f"grid={expected:.12e}"
                )

        # Separate flat/wind=0 branch: execution, convergence, and principal U.
        run(ocean_common(binary, 0.0) + grid_args + [str(wind0_csv)], root)
        wind0 = read_rows(wind0_csv)
        for raa in (0.0, 180.0):
            w0 = by_geometry(wind0, 60.0, raa)
            if int(round(w0["water_converged"])) != 1:
                raise RuntimeError(f"wind=0 did not converge at RAA={raa}")
            rrs_u = w0["Lu0minus_U"] / w0["Ed0minus_water"]
            if abs(rrs_u) / max(abs(w0["rrs_I"]), 1.0e-30) > 1.0e-10:
                raise RuntimeError(f"wind=0 principal-plane rrs U nonzero: {rrs_u}")

        print(
            "PASS: water RAA output regression "
            f"(water/TOA ratio spread={100.0*spread:.3f}%, "
            "full-grid parity and wind=0 branch passed)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
