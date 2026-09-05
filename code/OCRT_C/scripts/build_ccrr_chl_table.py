#!/usr/bin/env python3
"""Build the OCRT CCRR Chl absorption table from a normalized Morel spectrum.

The source spectrum A_chl(lambda) is dimensionless and normalized to one at
440 nm.  The classic Case-1 closure used here is

    a_p(lambda) = 0.06 * A_chl(lambda) * Chl**0.65  [m^-1]

The OCRT legacy CCRR loader expects five numeric columns:

    wavelength_nm  Ap  Ep  Aphi  Ephi

Only Aphi and Ephi are consumed by the current loader.  Ap and Ep are written
as explicit zero placeholders rather than implying unsupported semantics.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

SCALE = 0.06
EXPONENT = 0.65


def read_normalized_spectrum(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    for raw in path.read_text(encoding="utf-8", errors="strict").splitlines():
        fields = raw.split()
        if len(fields) < 2:
            continue
        try:
            wavelength = float(fields[0])
            value = float(fields[1])
        except ValueError:
            continue
        if wavelength < 0.0 and value < 0.0:
            break
        if wavelength <= 0.0 or value < 0.0 or not (math.isfinite(wavelength) and math.isfinite(value)):
            raise ValueError(f"invalid data row: {raw!r}")
        rows.append((wavelength, value))

    if len(rows) < 2:
        raise ValueError("fewer than two valid spectral rows")
    for (wl0, _), (wl1, _) in zip(rows, rows[1:]):
        if not wl1 > wl0:
            raise ValueError(f"wavelengths are not strictly increasing: {wl0} -> {wl1}")

    value_440 = next((value for wavelength, value in rows if abs(wavelength - 440.0) < 1e-12), None)
    if value_440 is None or abs(value_440 - 1.0) > 1e-12:
        raise ValueError(f"source must contain A_chl(440 nm)=1; got {value_440!r}")
    return rows


def render(rows: list[tuple[float, float]], source_name: str) -> str:
    lines = [
        "! OCRT CCRR chlorophyll absorption coefficient table",
        "! Canonical model: Morel (1988) normalized A_chl(lambda) with the",
        "! Morel-Maritorena / HydroLight classic Case-1 closure:",
        "!   a_p(lambda) = Aphi(lambda) * Chl^Ephi(lambda)",
        "!   Aphi(lambda) = 0.06 * A_chl_normalized(lambda)  [m^-1 at Chl=1]",
        "!   Ephi(lambda) = 0.65",
        "! Input Chl unit: mg m^-3. Output a_p unit: m^-1.",
        f"! Source spectrum: {source_name}",
        "! The source is normalized to 1 at 440 nm and identifies the spectral",
        "! shape as Morel (1988, Fig. 10c; Prieur and Sathyendranath 1981).",
        "! Source warning retained: 300-350 and 700-1000 nm are extrapolated",
        "! and may be unrealistic. OCRT uses linear interpolation and nearest-edge",
        "! extrapolation outside the tabulated range.",
        "!",
        "! Legacy OCRT five-column interface:",
        "! wavelength_nm  Ap_unused  Ep_unused  Aphi  Ephi",
        "! Ap_unused and Ep_unused are explicit zero placeholders; current OCRT",
        "! consumes only Aphi and Ephi in rt_iop_ccrr_pigment_eval().",
    ]
    for wavelength, normalized in rows:
        aphi = SCALE * normalized
        lines.append(f"{wavelength:7.1f}  0.000000000  0.000000000  {aphi:.9f}  {EXPONENT:.9f}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    rows = read_normalized_spectrum(args.source)
    text = render(rows, args.source.name)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
