#!/usr/bin/env python3
"""Compare one single-target solve against the matching full-grid row."""
from __future__ import annotations
import argparse, csv, math
from pathlib import Path


def parse_stdout(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    for token in path.read_text(encoding="utf-8").replace("\n", " ").split():
        if "=" not in token:
            continue
        k, v = token.split("=", 1)
        try:
            out[k] = float(v)
        except ValueError:
            pass
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--single", type=Path, required=True)
    ap.add_argument("--full-grid", type=Path, required=True)
    ap.add_argument("--vza", type=float, default=30.0)
    ap.add_argument("--raa", type=float, default=45.0)
    args = ap.parse_args()

    got = parse_stdout(args.single)
    with args.full_grid.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    ref = next((r for r in rows if abs(float(r["vza_deg"])-args.vza)<1e-12 and abs(float(r["raa_deg"])-args.raa)<1e-12), None)
    if ref is None:
        raise SystemExit("matching full-grid row not found")
    ed = float(ref["Ed0minus_water"])
    expected = {
        "TOA_rho_I": float(ref["TOA_rho_I"]),
        "TOA_rho_Q": float(ref["TOA_rho_Q"]),
        "TOA_rho_U": float(ref["TOA_rho_U"]),
        "Rrs0plus_I": float(ref["Rrs_I"]),
        "Rrs0plus_Q": float(ref["Rrs_Q"]),
        "Rrs0plus_U": float(ref["Rrs_U"]),
        "rrs0minus_I": float(ref["rrs_I"]),
        "rrs0minus_Q": float(ref["Lu0minus_Q"])/ed,
        "rrs0minus_U": float(ref["Lu0minus_U"])/ed,
    }
    for key, exp in expected.items():
        if key not in got:
            raise SystemExit(f"single output missing {key}")
        tol = 2e-6 if key.startswith("TOA_") else 2e-8
        if not math.isclose(got[key], exp, rel_tol=0.0, abs_tol=tol):
            raise SystemExit(f"mismatch {key}: {got[key]:.17g} vs {exp:.17g}")
    print("PASS: single-target/full-grid consistency at vza=30 raa=45 n_mu_water=64")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
