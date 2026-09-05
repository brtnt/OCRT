#!/usr/bin/env python3
"""Platform-neutral prebuild helper for OCRT persistent surface operators.

The helper performs one explicit single-prebuilder run, then writes a SHA-256
manifest for the immutable cache files.  Simulation rows must subsequently use
OCRT_SURFACE_PERSIST_CACHE_MODE=required or readonly; they never build or wait.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys
from typing import Iterable


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(cache_dir: Path) -> Path:
    rows = []
    for path in sorted(cache_dir.glob("ocrt_surface_v*.bin"), key=lambda p: p.name):
        rows.append({
            "file": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    manifest = cache_dir / "SURFACE_CACHE_SHA256.csv"
    with manifest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["file", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Prebuild OCRT persistent surface cache without OS-specific launch logic."
    )
    parser.add_argument("--ocrt", required=True, type=Path, help="OCRT executable")
    parser.add_argument("--cache-dir", required=True, type=Path, help="Cache directory")
    parser.add_argument(
        "ocrt_args",
        nargs=argparse.REMAINDER,
        help="Arguments passed to OCRT; prefix with -- to end helper options",
    )
    ns = parser.parse_args(list(argv) if argv is not None else None)

    exe = ns.ocrt.expanduser().resolve()
    if not exe.is_file():
        parser.error(f"OCRT executable not found: {exe}")

    cache_dir = ns.cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    args = list(ns.ocrt_args)
    if args and args[0] == "--":
        args = args[1:]

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OCRT_SURFACE_PERSIST_CACHE_MODE"] = "build"
    env["OCRT_SURFACE_PERSIST_CACHE_DIR"] = str(cache_dir)

    proc = subprocess.run([str(exe), *args], env=env, check=False)
    if proc.returncode != 0:
        print(
            f"error: OCRT surface-cache prebuild failed with exit code {proc.returncode}",
            file=sys.stderr,
        )
        return proc.returncode or 1

    manifest = write_manifest(cache_dir)
    print(f"Persistent surface-cache prebuild complete: {cache_dir}")
    print(f"Manifest: {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
