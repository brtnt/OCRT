#!/usr/bin/env python3
"""Remap validated legacy AHN TSM phase matrices to the OCRT FR631 grid.

This is deliberately a *phase-preserving remap*, not a new Mie regeneration.
The historical AHN generator/.inp contract is not present in the handoff.

Contract:
- scalar bulk block: retain the current 330--1100 nm spectral extension;
- Asymm_Para: replace by interpolation of the validated legacy bulk g table;
- phase wavelengths: 0.330, validated legacy nodes through 0.860, and 1.100 um;
- 0.330 um phase: conservative endpoint hold of the legacy 0.350 um phase;
- 1.100 um phase: linear interpolation between legacy 0.860 and 1.240 um;
- phase-angle grid: fixed FR631, evaluated by theta-linear interpolation;
- P11/P12/P33 use identical wavelength and angular weights.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

TARGET_PHASE_WL_UM = np.array(
    [0.330, 0.350, 0.400, 0.412, 0.443, 0.470, 0.488, 0.515,
     0.550, 0.590, 0.633, 0.670, 0.694, 0.760, 0.860, 1.100],
    dtype=np.float64,
)


def fr631_angles_ascending() -> np.ndarray:
    parts = [
        np.arange(0.0, 0.2 + 0.5 * 0.005, 0.005),
        np.arange(0.2, 1.0 + 0.5 * 0.02, 0.02),
        np.arange(1.0, 5.0 + 0.5 * 0.05, 0.05),
        np.arange(5.0, 20.0 + 0.5 * 0.1, 0.1),
        np.arange(20.0, 180.0 + 0.5 * 0.5, 0.5),
    ]
    grid = np.concatenate([parts[0], *(p[1:] for p in parts[1:])])
    if grid.shape != (631,) or not np.all(np.diff(grid) > 0.0):
        raise RuntimeError(f"invalid FR631 grid: shape={grid.shape}")
    if abs(grid[0]) > 1e-14 or abs(grid[-1] - 180.0) > 1e-12:
        raise RuntimeError("invalid FR631 endpoints")
    return grid


@dataclass(frozen=True)
class MieText:
    spectral_header: str
    spectral: np.ndarray  # (N, 7)
    phase_wl_um: np.ndarray  # (W,)
    angles_deg: np.ndarray  # (A,), file order
    p11: np.ndarray  # (A, W)
    p12: np.ndarray
    p33: np.ndarray


def _parse_phase_block(lines: list[str], marker_index: int) -> tuple[np.ndarray, np.ndarray, int]:
    header = lines[marker_index + 1].split()
    if not header or header[0].upper() != "TETA":
        raise ValueError(f"missing TETA header after line {marker_index + 1}")
    wl = np.asarray([float(x) for x in header[1:]], dtype=np.float64)
    rows: list[list[float]] = []
    idx = marker_index + 2
    while idx < len(lines) and "Phase Function" not in lines[idx]:
        parts = lines[idx].split()
        if len(parts) == 1 + wl.size:
            rows.append([float(x) for x in parts])
        idx += 1
    arr = np.asarray(rows, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 1 + wl.size:
        raise ValueError("malformed phase block")
    return wl, arr, idx


def read_mie_text(path: Path) -> MieText:
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) < 10:
        raise ValueError(f"too short: {path}")
    n_ang = int(lines[0].strip())
    spectral_header = lines[1]
    p11_idx = next(i for i, line in enumerate(lines) if "Phase Function (P11)" in line)
    spectral_rows: list[list[float]] = []
    for line in lines[2:p11_idx]:
        parts = line.split()
        if len(parts) == 7:
            spectral_rows.append([float(x) for x in parts])
    spectral = np.asarray(spectral_rows, dtype=np.float64)
    if spectral.ndim != 2 or spectral.shape[1] != 7:
        raise ValueError(f"malformed spectral block: {path}")

    wl11, block11, next_idx = _parse_phase_block(lines, p11_idx)
    p12_idx = next(i for i in range(next_idx, len(lines)) if "Phase Function (P12" in lines[i])
    wl12, block12, next_idx = _parse_phase_block(lines, p12_idx)
    p33_idx = next(i for i in range(next_idx, len(lines)) if "Phase Function (P33" in lines[i])
    wl33, block33, _ = _parse_phase_block(lines, p33_idx)
    if not (np.array_equal(wl11, wl12) and np.array_equal(wl11, wl33)):
        raise ValueError("P11/P12/P33 wavelength headers differ")
    if not (np.array_equal(block11[:, 0], block12[:, 0]) and
            np.array_equal(block11[:, 0], block33[:, 0])):
        raise ValueError("P11/P12/P33 angle grids differ")
    if block11.shape[0] != n_ang:
        raise ValueError(f"declared n_ang={n_ang}, parsed={block11.shape[0]}")
    return MieText(
        spectral_header=spectral_header,
        spectral=spectral,
        phase_wl_um=wl11,
        angles_deg=block11[:, 0],
        p11=block11[:, 1:],
        p12=block12[:, 1:],
        p33=block33[:, 1:],
    )


def interp_wavelength_common(values: np.ndarray, src_wl: np.ndarray,
                             target_wl: np.ndarray) -> np.ndarray:
    """Common-weight linear interpolation; 0.330 uses 0.350 endpoint hold."""
    if values.shape[1] != src_wl.size:
        raise ValueError("phase wavelength dimension mismatch")
    result = np.empty((values.shape[0], target_wl.size), dtype=np.float64)
    for j, wl in enumerate(target_wl):
        if wl <= src_wl[0]:
            # Explicit science contract: 330--349 nm uses the validated 350-nm shape.
            result[:, j] = values[:, 0]
            continue
        if wl >= src_wl[-1]:
            result[:, j] = values[:, -1]
            continue
        hi = int(np.searchsorted(src_wl, wl, side="right"))
        lo = hi - 1
        if abs(wl - src_wl[lo]) < 1e-14:
            result[:, j] = values[:, lo]
            continue
        if abs(wl - src_wl[hi]) < 1e-14:
            result[:, j] = values[:, hi]
            continue
        weight = (wl - src_wl[lo]) / (src_wl[hi] - src_wl[lo])
        result[:, j] = values[:, lo] + weight * (values[:, hi] - values[:, lo])
    return result


def remap_phase_theta(values: np.ndarray, src_angle: np.ndarray,
                      target_angle_asc: np.ndarray) -> np.ndarray:
    order = np.argsort(src_angle)
    x = src_angle[order]
    y = values[order, :]
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("legacy phase angle grid is not strictly monotone")
    out = np.empty((target_angle_asc.size, y.shape[1]), dtype=np.float64)
    for j in range(y.shape[1]):
        out[:, j] = np.interp(target_angle_asc, x, y[:, j])
    return out


def legacy_g_at_bulk_wavelengths(legacy: MieText, wl_um: np.ndarray) -> np.ndarray:
    src_wl = legacy.spectral[:, 0]
    src_g = legacy.spectral[:, 4]
    if not np.all(np.diff(src_wl) > 0.0):
        raise ValueError("legacy bulk wavelengths are not monotone")
    # Explicit endpoint hold only below the validated 350-nm lower endpoint.
    return np.interp(wl_um, src_wl, src_g, left=src_g[0], right=src_g[-1])


def write_mie(path: Path, spectral_header: str, spectral: np.ndarray,
              phase_wl_um: np.ndarray, angle_desc: np.ndarray,
              p11_desc: np.ndarray, p12_desc: np.ndarray,
              p33_desc: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write(f"{angle_desc.size:4d}\n")
        f.write("   Wlgth  Nor_Ext_Co  Nor_Sca_Co  Sg_Sca_Alb  Asymm_Para  Extinct_Co  Scatter_Co\n")
        for row in spectral:
            f.write(
                f"    {row[0]:0.4f}     {row[1]:.10f}        {row[2]:.10f}"
                f"        {row[3]:.10f}        {row[4]:.10f}"
                f"     {row[5]:.12E}  {row[6]:.12E}\n"
            )
        headers = (
            ("                    Phase Function (P11)", p11_desc),
            ("                    Phase Function (P12 / Q-polarization)", p12_desc),
            ("                    Phase Function (P33 / U-polarization)", p33_desc),
        )
        wl_header = "   TETA" + "".join(f"    {wl:0.4f}" for wl in phase_wl_um)
        for marker, values in headers:
            f.write(marker + "\n")
            f.write(wl_header + "\n")
            for theta, row in zip(angle_desc, values, strict=True):
                f.write(f"{theta:9.6f}" + "".join(f"  {v:+.8E}" for v in row) + "\n")


def validate_output(mie: MieText) -> None:
    expected = fr631_angles_ascending()[::-1]
    if mie.angles_deg.shape != expected.shape or np.max(np.abs(mie.angles_deg - expected)) > 5e-7:
        raise ValueError("output does not use FR631")
    if np.max(np.abs(mie.phase_wl_um - TARGET_PHASE_WL_UM)) > 5e-12:
        raise ValueError("unexpected output phase wavelengths")
    if not np.isfinite(mie.p11).all() or not np.isfinite(mie.p12).all() or not np.isfinite(mie.p33).all():
        raise ValueError("non-finite phase value")
    tol = 5e-10
    if np.min(mie.p11) < -tol:
        raise ValueError(f"negative P11: {np.min(mie.p11)}")
    if np.max(np.abs(mie.p12) - mie.p11) > tol:
        raise ValueError("|P12| > P11")
    if np.max(np.abs(mie.p33) - mie.p11) > tol:
        raise ValueError("|P33| > P11")
    if mie.spectral[0, 0] > 0.3300001 or mie.spectral[-1, 0] < 1.0999999:
        raise ValueError("bulk spectral table does not cover 330--1100 nm")


def process_one(legacy_path: Path, current_path: Path, output_path: Path) -> dict[str, float | int | str]:
    legacy = read_mie_text(legacy_path)
    current = read_mie_text(current_path)
    spectral = current.spectral.copy()
    spectral[:, 4] = legacy_g_at_bulk_wavelengths(legacy, spectral[:, 0])

    p11_w = interp_wavelength_common(legacy.p11, legacy.phase_wl_um, TARGET_PHASE_WL_UM)
    p12_w = interp_wavelength_common(legacy.p12, legacy.phase_wl_um, TARGET_PHASE_WL_UM)
    p33_w = interp_wavelength_common(legacy.p33, legacy.phase_wl_um, TARGET_PHASE_WL_UM)

    target_asc = fr631_angles_ascending()
    p11_asc = remap_phase_theta(p11_w, legacy.angles_deg, target_asc)
    p12_asc = remap_phase_theta(p12_w, legacy.angles_deg, target_asc)
    p33_asc = remap_phase_theta(p33_w, legacy.angles_deg, target_asc)

    write_mie(
        output_path,
        current.spectral_header,
        spectral,
        TARGET_PHASE_WL_UM,
        target_asc[::-1],
        p11_asc[::-1, :],
        p12_asc[::-1, :],
        p33_asc[::-1, :],
    )
    out = read_mie_text(output_path)
    validate_output(out)

    # Exact overlap at historical 0.5-degree nodes and common phase wavelengths.
    common_wl = [0.350, 0.400, 0.412, 0.443, 0.470, 0.488, 0.515,
                 0.550, 0.590, 0.633, 0.670, 0.694, 0.760, 0.860]
    max_overlap = 0.0
    old_angle_index = {round(float(v), 8): i for i, v in enumerate(legacy.angles_deg)}
    new_angle_index = {round(float(v), 8): i for i, v in enumerate(out.angles_deg)}
    old_wl_index = {round(float(v), 8): i for i, v in enumerate(legacy.phase_wl_um)}
    new_wl_index = {round(float(v), 8): i for i, v in enumerate(out.phase_wl_um)}
    for theta_key, oi in old_angle_index.items():
        ni = new_angle_index[theta_key]
        for wl in common_wl:
            oj = old_wl_index[round(wl, 8)]
            nj = new_wl_index[round(wl, 8)]
            for old_mat, new_mat in ((legacy.p11, out.p11), (legacy.p12, out.p12), (legacy.p33, out.p33)):
                max_overlap = max(max_overlap, abs(float(old_mat[oi, oj] - new_mat[ni, nj])))
    return {
        "file": output_path.name,
        "n_angles": int(out.angles_deg.size),
        "n_phase_wavelengths": int(out.phase_wl_um.size),
        "phase_min_um": float(out.phase_wl_um[0]),
        "phase_max_um": float(out.phase_wl_um[-1]),
        "max_common_node_abs_difference": max_overlap,
        "min_p11": float(np.min(out.p11)),
        "max_p12_bound_excess": float(np.max(np.abs(out.p12) - out.p11)),
        "max_p33_bound_excess": float(np.max(np.abs(out.p33) - out.p11)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-dir", type=Path, required=True)
    parser.add_argument("--current-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    names = [
        "Red_clay_AHN.mie",
        "Brown_earth_AHN.mie",
        "Yellow_clay_AHN.mie",
        "Calcareous_sand_AHN.mie",
    ]
    summaries = [
        process_one(args.legacy_dir / name, args.current_dir / name, args.output_dir / name)
        for name in names
    ]
    import json
    text = json.dumps({"contract": "AHN-LEGACY-PHASE-FR631-330HOLD-2026-08-16", "files": summaries}, indent=2)
    if args.summary:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
