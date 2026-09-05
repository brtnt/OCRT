#!/usr/bin/env python3
"""Decode an OSOAA TAW Fourier-matrix file and compare it with OCRT's
surface_T_aw_coxmunk_fourier_kernel on the identical angular grid.

The OSOAA matrix is converted to the kernel normalization used by OCRT from
OSOAA_SOS_CORE.F step 7:

    S_w(K) = 2/mu_w(K) * sum_J w_J TAW(J,K) S_a(J)

and OCRT rt_air_water_couple_atm_to_water:

    S_w(K) = C_m * sum_J mu_a(J) w_J K_m(K,J) S_a(J),
    C_0=2*pi, C_m=pi for m>0.

Thus K_OSOAA = 2*TAW/(C_m*mu_w*mu_a).  The optional (-1)^m factor maps the
OSOAA public azimuth to the OCRT public RAA convention by a 180-degree shift.
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
from typing import Iterable

import numpy as np

COMPONENTS = ("11", "12", "13", "21", "22", "23", "31", "32", "33")
COS_ENTRY = np.array([True, True, False, True, True, False, False, False, True])
MARKERS = ("o", "s", "^", "D", "x", "+", "v", "<", ">")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_angles(path: Path) -> tuple[np.ndarray, np.ndarray]:
    pat = re.compile(
        r"\s*(\d+)\s+([+-]?\d+\.\d+D[+-]\d+)\s+"
        r"([+-]?\d+\.\d+D[+-]\d+)\s+(\d+)"
    )
    mu: list[float] = []
    wt: list[float] = []
    for line in path.read_text(errors="replace").splitlines():
        m = pat.match(line)
        if m:
            mu.append(float(m.group(2).replace("D", "E")))
            wt.append(float(m.group(3).replace("D", "E")))
    if not mu:
        raise ValueError(f"No angle rows parsed from {path}")
    return np.asarray(mu, dtype=np.float64), np.asarray(wt, dtype=np.float64)


def read_osoaa_taw(path: Path, n_mu: int, m_max: int) -> np.ndarray:
    payload_len = 9 * n_mu * n_mu * 4
    records: list[np.ndarray] = []
    with path.open("rb") as f:
        for m in range(m_max + 1):
            head = f.read(4)
            if len(head) != 4:
                raise EOFError(f"Missing record {m} in {path}")
            (nbyte,) = struct.unpack("<i", head)
            if nbyte != payload_len:
                raise ValueError(
                    f"Unexpected record payload at m={m}: {nbyte}; expected {payload_len}"
                )
            payload = f.read(nbyte)
            tail = f.read(4)
            if len(payload) != nbyte or len(tail) != 4:
                raise EOFError(f"Truncated record {m} in {path}")
            (nbyte2,) = struct.unpack("<i", tail)
            if nbyte2 != nbyte:
                raise ValueError(f"Fortran record marker mismatch at m={m}")
            values = np.frombuffer(payload, dtype="<f4")
            mats = np.stack(
                [
                    values[k * n_mu * n_mu : (k + 1) * n_mu * n_mu].reshape(
                        (n_mu, n_mu), order="F"
                    )
                    for k in range(9)
                ]
            )
            records.append(mats)
    return np.stack(records).astype(np.float64)


def compile_surface_library(ocrt_root: Path, out_so: Path) -> None:
    shared = ocrt_root / "src" / "shared"
    surface = shared / "surface.c"
    multibounce = shared / "surface_multibounce.c"
    for path in (surface, multibounce, shared / "surface.h", shared / "mat3.h"):
        if not path.exists():
            raise FileNotFoundError(path)
    cmd = [
        "gcc",
        "-O3",
        "-fPIC",
        "-shared",
        "-DOCRT_FAST_KERNELS",
        f"-I{shared}",
        str(surface),
        str(multibounce),
        "-lm",
        "-o",
        str(out_so),
    ]
    subprocess.run(cmd, check=True)


def load_ocrt_kernel(so_path: Path):
    lib = ctypes.CDLL(str(so_path))
    fn = lib.surface_T_aw_coxmunk_fourier_kernel
    fn.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double),
    ]
    fn.restype = ctypes.c_int
    return fn


def call_ocrt_kernel(
    fn,
    mu: np.ndarray,
    mode: int,
    n_phi: int,
    wind: float,
    sigma_type: int,
    n_water: float,
    q_convention: int,
) -> np.ndarray:
    out = np.zeros((mu.size, mu.size, 9), dtype=np.float64)
    rc = fn(
        mu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        int(mu.size),
        mu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        int(mu.size),
        int(mode),
        int(n_phi),
        float(wind),
        int(sigma_type),
        float(n_water),
        int(q_convention),
        out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )
    if rc != 0:
        raise RuntimeError(f"surface_T_aw_coxmunk_fourier_kernel rc={rc}")
    return out


def osoaa_to_ocrt_kernel(
    osoaa_mode: np.ndarray, mu: np.ndarray, mode: int, public_raa_shift: bool
) -> np.ndarray:
    # OSOAA array: [component, incident-air I, transmitted-water J].
    # OCRT array:  [transmitted-water J, incident-air I, component].
    out = osoaa_mode.transpose(2, 1, 0).copy()
    c_m = 2.0 * math.pi if mode == 0 else math.pi
    out *= (2.0 / (c_m * mu[:, None] * mu[None, :]))[:, :, None]
    if public_raa_shift:
        out *= (-1.0) ** mode
    return out


def metrics(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    max_ref = float(np.max(np.abs(x)))
    if max_ref <= 1e-30:
        return {
            "best_fit_ocrt_over_osoaa": float("nan"),
            "correlation": float("nan"),
            "mean_abs_over_max_ref_pct": float("nan"),
            "rmse_over_max_ref_pct": float("nan"),
            "max_abs_over_max_ref_pct": float("nan"),
        }
    best = float(np.dot(x, y) / np.dot(x, x)) if np.dot(x, x) > 0 else float("nan")
    corr = float(np.corrcoef(x, y)[0, 1]) if x.size > 2 else float("nan")
    d = y - x
    return {
        "best_fit_ocrt_over_osoaa": best,
        "correlation": corr,
        "mean_abs_over_max_ref_pct": 100.0 * float(np.mean(np.abs(d))) / max_ref,
        "rmse_over_max_ref_pct": 100.0 * float(np.sqrt(np.mean(d * d))) / max_ref,
        "max_abs_over_max_ref_pct": 100.0 * float(np.max(np.abs(d))) / max_ref,
    }


def write_csv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(
    outdir: Path,
    long_rows: list[dict[str, object]],
    m_max: int,
) -> list[Path]:
    import matplotlib.pyplot as plt

    files: list[Path] = []
    for component in COMPONENTS:
        rows = [r for r in long_rows if r["component"] == component]
        if not rows:
            continue
        x_all = np.asarray([float(r["osoaa_kernel"]) for r in rows])
        y_all = np.asarray([float(r["ocrt_kernel"]) for r in rows])
        lim = max(float(np.max(np.abs(x_all))), float(np.max(np.abs(y_all))), 1e-15)
        fig, ax = plt.subplots(figsize=(7.6, 7.0))
        for mode in range(m_max + 1):
            rr = [r for r in rows if int(r["mode"]) == mode]
            if not rr:
                continue
            x = [float(r["osoaa_kernel"]) for r in rr]
            y = [float(r["ocrt_kernel"]) for r in rr]
            ax.scatter(x, y, s=12, alpha=0.45, marker=MARKERS[mode % len(MARKERS)], label=f"m={mode}")
        ax.plot([-lim, lim], [-lim, lim], linestyle="--", linewidth=1.1, label="1:1")
        ax.axhline(0.0, linewidth=0.7)
        ax.axvline(0.0, linewidth=0.7)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(f"OSOAA comparable TAW M{component}")
        ax.set_ylabel(f"OCRT TAW M{component}")
        ax.set_title(f"Air-to-water Fourier operator audit — M{component}")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=8)
        fig.tight_layout()
        path = outdir / f"TAW_M{component}_scatter.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        files.append(path)
    return files


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--ocrt-root", type=Path, required=True)
    p.add_argument("--osoaa-taw", type=Path, required=True)
    p.add_argument("--angles", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--m-max", type=int, default=4)
    p.add_argument("--n-phi", type=int, default=1024)
    p.add_argument("--wind", type=float, default=3.0)
    p.add_argument("--n-water", type=float, default=1.34)
    p.add_argument("--sigma-type", type=int, default=1)
    p.add_argument("--q-convention", type=int, default=1)
    p.add_argument("--include-zero-weight", action="store_true")
    p.add_argument("--native-azimuth", action="store_true", help="Do not apply the OCRT/OSOAA 180-degree public RAA shift")
    p.add_argument("--no-plots", action="store_true")
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    mu, weight = parse_angles(args.angles)
    osoaa = read_osoaa_taw(args.osoaa_taw, int(mu.size), args.m_max)
    active = np.ones(mu.size, dtype=bool) if args.include_zero_weight else weight > 0.0

    so_path = args.outdir / "libocrt_surface_audit.so"
    compile_surface_library(args.ocrt_root, so_path)
    fn = load_ocrt_kernel(so_path)

    summary_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []
    for mode in range(args.m_max + 1):
        ocrt = call_ocrt_kernel(
            fn,
            mu,
            mode,
            args.n_phi,
            args.wind,
            args.sigma_type,
            args.n_water,
            args.q_convention,
        )
        ref = osoaa_to_ocrt_kernel(osoaa[mode], mu, mode, not args.native_azimuth)
        for c, component in enumerate(COMPONENTS):
            x = ref[:, :, c][np.ix_(active, active)].ravel()
            y = ocrt[:, :, c][np.ix_(active, active)].ravel()
            mm = metrics(x, y)
            summary_rows.append(
                {
                    "mode": mode,
                    "component": component,
                    "n_pairs": int(x.size),
                    **mm,
                }
            )
        active_idx = np.where(active)[0]
        for jw in active_idx:
            for ia in active_idx:
                for c, component in enumerate(COMPONENTS):
                    long_rows.append(
                        {
                            "mode": mode,
                            "component": component,
                            "water_index_1based": int(jw + 1),
                            "air_index_1based": int(ia + 1),
                            "mu_water": float(mu[jw]),
                            "mu_air": float(mu[ia]),
                            "osoaa_kernel": float(ref[jw, ia, c]),
                            "ocrt_kernel": float(ocrt[jw, ia, c]),
                            "delta": float(ocrt[jw, ia, c] - ref[jw, ia, c]),
                        }
                    )

    summary_path = args.outdir / "TAW_OPERATOR_AUDIT_SUMMARY.csv"
    write_csv(
        summary_path,
        summary_rows,
        [
            "mode",
            "component",
            "n_pairs",
            "best_fit_ocrt_over_osoaa",
            "correlation",
            "mean_abs_over_max_ref_pct",
            "rmse_over_max_ref_pct",
            "max_abs_over_max_ref_pct",
        ],
    )
    long_path = args.outdir / "TAW_OPERATOR_AUDIT_LONG.csv"
    write_csv(
        long_path,
        long_rows,
        [
            "mode",
            "component",
            "water_index_1based",
            "air_index_1based",
            "mu_water",
            "mu_air",
            "osoaa_kernel",
            "ocrt_kernel",
            "delta",
        ],
    )
    plot_files: list[Path] = []
    if not args.no_plots:
        plot_files = make_plots(args.outdir, long_rows, args.m_max)

    manifest = {
        "ocrt_root": str(args.ocrt_root.resolve()),
        "ocrt_surface_c_sha256": sha256(args.ocrt_root / "src" / "shared" / "surface.c"),
        "osoaa_taw": str(args.osoaa_taw.resolve()),
        "osoaa_taw_sha256": sha256(args.osoaa_taw),
        "angles": str(args.angles.resolve()),
        "angles_sha256": sha256(args.angles),
        "n_total_angles": int(mu.size),
        "n_active_gauss_angles": int(np.count_nonzero(active)),
        "m_max": args.m_max,
        "n_phi": args.n_phi,
        "wind": args.wind,
        "n_water": args.n_water,
        "sigma_type": args.sigma_type,
        "q_convention": args.q_convention,
        "public_raa_shift_applied": not args.native_azimuth,
        "normalization": "K_OSOAA=2*TAW/(C_m*mu_water*mu_air), C0=2pi, Cm=pi; optional (-1)^m public-RAA shift",
        "files": [summary_path.name, long_path.name] + [f.name for f in plot_files],
    }
    (args.outdir / "TAW_OPERATOR_AUDIT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    # Compact console result: scalar anchor and polarization diagnostics.
    by_key = {(int(r["mode"]), str(r["component"])): r for r in summary_rows}
    print(f"active Gauss nodes: {np.count_nonzero(active)}/{mu.size}")
    for mode in range(args.m_max + 1):
        r11 = by_key[(mode, "11")]
        print(
            f"m={mode}: M11 best-scale={r11['best_fit_ocrt_over_osoaa']:.8f}, "
            f"max/maxref={r11['max_abs_over_max_ref_pct']:.6f}%"
        )
    print(f"wrote {summary_path}")
    print(f"wrote {long_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
