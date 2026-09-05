#!/usr/bin/env python3
"""Compare two CSV files without pandas numeric parsing.

Numeric cells are parsed with Python's built-in float() and compared by their
IEEE-754 binary64 bit pattern. Non-numeric cells are compared as exact text.
Use --raw to require byte-for-byte file identity.
"""
from __future__ import annotations

import argparse
import csv
import math
import struct
import sys
from pathlib import Path


def float_bits(text: str) -> bytes | None:
    try:
        value = float(text)
    except ValueError:
        return None
    if math.isnan(value):
        return b"nan"
    return struct.pack(">d", value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--raw", action="store_true", help="require byte identity")
    parser.add_argument("--ignore-column", action="append", default=[])
    args = parser.parse_args()

    if args.raw:
        same = args.left.read_bytes() == args.right.read_bytes()
        print("PASS: raw byte identity" if same else "FAIL: raw bytes differ")
        return 0 if same else 1

    with args.left.open(newline="") as fa, args.right.open(newline="") as fb:
        ra = csv.DictReader(fa)
        rb = csv.DictReader(fb)
        if ra.fieldnames != rb.fieldnames:
            print("FAIL: CSV headers differ", file=sys.stderr)
            return 1
        ignored = set(args.ignore_column)
        for row_index, (a, b) in enumerate(zip(ra, rb), start=2):
            for key in ra.fieldnames or []:
                if key in ignored:
                    continue
                av, bv = a[key], b[key]
                ab, bb = float_bits(av), float_bits(bv)
                equal = (ab == bb) if ab is not None and bb is not None else (av == bv)
                if not equal:
                    print(
                        f"FAIL: row={row_index} column={key} left={av!r} right={bv!r}",
                        file=sys.stderr,
                    )
                    return 1
        try:
            next(ra)
            print("FAIL: left has extra rows", file=sys.stderr)
            return 1
        except StopIteration:
            pass
        try:
            next(rb)
            print("FAIL: right has extra rows", file=sys.stderr)
            return 1
        except StopIteration:
            pass
    print("PASS: CSV values are bit-identical under built-in float parsing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
