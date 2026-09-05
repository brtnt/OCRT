#!/usr/bin/env python3
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PACKAGE_ROOT / "config" / "polarization_campaign_pssa_16core.json"
MANIFEST_FIELDS = [
    "run_id", "stage", "status", "returncode", "elapsed_s", "rows",
    "output_relpath", "output_bytes", "output_sha256", "binary_sha256",
    "command", "stdout_log", "stderr_log", "message"
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def find_binary(runtime_root: Path, explicit: Path | None) -> Path:
    if explicit:
        path = explicit.resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    croot = runtime_root / "01_OCRT_C"
    candidates = [
        croot / "build" / "ocrt_polarization",
        croot / "build" / "ocrt_polarization.exe",
        croot / "build" / "ocrt",
        croot / "build" / "ocrt.exe",
        croot / "build" / "ocrt_v1.2",
        croot / "build" / "ocrt_v1.2.exe",
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"OCRT executable not found; searched {candidates}")


def read_matrix(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as f:
        raw = list(csv.DictReader(f))
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for line_no, row in enumerate(raw, start=2):
        run_id = row.get("run_id", "").strip()
        if not run_id:
            continue
        if run_id in seen:
            raise ValueError(f"duplicate run_id at line {line_no}: {run_id}")
        seen.add(run_id)
        try:
            args = json.loads(row["command_args_json"])
            required_columns = json.loads(row.get("required_columns_json") or "[]")
            expected_vza = json.loads(row.get("expected_vza_json") or "[]")
            expected_raa = json.loads(row.get("expected_raa_json") or "[]")
        except Exception as exc:
            raise ValueError(f"matrix JSON parse error line {line_no}: {exc}") from exc
        if not isinstance(args, list) or not all(isinstance(x, str) for x in args):
            raise ValueError(f"command_args_json must be a JSON string list at line {line_no}")
        rows.append({
            "run_id": run_id,
            "stage": row.get("stage", "").strip(),
            "uses_water": parse_bool(row.get("uses_water", "false")),
            "expected_rows": int(row.get("expected_rows") or 0),
            "output_relpath": row.get("output_relpath", "").strip(),
            "args": args,
            "required_columns": required_columns,
            "expected_vza": expected_vza,
            "expected_raa": expected_raa,
            "notes": row.get("notes", ""),
        })
    if not rows:
        raise RuntimeError(
            f"matrix contains no executable rows: {path}. The historical 8,045-run matrix is intentionally not fabricated."
        )
    return rows


def option_present(args: list[str], option: str) -> bool:
    return option in args or any(x.startswith(option + "=") for x in args)


def value_after(args: list[str], option: str) -> str | None:
    for i, item in enumerate(args):
        if item == option:
            return args[i + 1] if i + 1 < len(args) else None
        if item.startswith(option + "="):
            return item.split("=", 1)[1]
    return None


def validate_row_policy(row: dict[str, Any], config: dict[str, Any]) -> list[str]:
    args = list(row["args"])
    if config.get("enforce_pssa", True) and not option_present(args, "--pssa"):
        args.append("--pssa")
    for forbidden in config.get("forbidden_cli_options", []):
        if option_present(args, forbidden):
            raise RuntimeError(f"{row['run_id']}: forbidden canonical-mode option present: {forbidden}")
    for controlled in ["--output-full-grid", "--batch-full-grid", "--output"]:
        if option_present(args, controlled):
            raise RuntimeError(f"{row['run_id']}: output option is controlled by runner: {controlled}")
    if row["uses_water"] and not option_present(args, "--n-mu-water"):
        raise RuntimeError(
            f"{row['run_id']}: water row lacks explicit advanced option --n-mu-water. "
            "The runner does not invent a numerical default."
        )
    if option_present(args, "--n-mu-water"):
        value = value_after(args, "--n-mu-water")
        if value is None or int(value) <= 0:
            raise RuntimeError(f"{row['run_id']}: invalid --n-mu-water value: {value}")
    row["args"] = args
    return args


def build_child_env(config: dict[str, Any]) -> dict[str, str]:
    env = os.environ.copy()
    for name in config.get("forbidden_environment_variables", []):
        if env.get(name):
            raise RuntimeError(
                f"forbidden environment variable is set in parent shell: {name}={env[name]!r}. "
                "Unset it before running the canonical campaign."
            )
        env.pop(name, None)
    for key, value in config.get("child_environment", {}).items():
        env[str(key)] = str(value)
    if env.get("OMP_NUM_THREADS") != "1":
        raise RuntimeError("canonical 16-worker policy requires OMP_NUM_THREADS=1 per child")
    return env


def open_csv_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", newline="", encoding="utf-8")
    return path.open("r", newline="", encoding="utf-8")


def validate_csv(
    path: Path,
    expected_rows: int,
    required_columns: list[str],
    expected_vza: list[Any],
    expected_raa: list[Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    rows = 0
    vzas: set[float] = set()
    raas: set[float] = set()
    with open_csv_text(path) as f:
        reader = csv.DictReader(f)
        fieldnames = set(reader.fieldnames or [])
        missing = set(required_columns) - fieldnames
        if missing:
            raise RuntimeError(f"{path}: missing columns {sorted(missing)}")
        numeric_candidates = [
            c for c in ["rho_I", "rho_Q", "rho_U", "Rrs_I", "Rrs_Q", "Rrs_U", "rrs_I", "rrs_Q", "rrs_U"]
            if c in fieldnames
        ]
        for record in reader:
            rows += 1
            if "vza_deg" in record and record["vza_deg"] != "":
                vzas.add(round(float(record["vza_deg"]), 10))
            if "raa_deg" in record and record["raa_deg"] != "":
                raas.add(round(float(record["raa_deg"]), 10))
            for col in numeric_candidates:
                value = float(record[col])
                if config["validation"].get("fail_on_nonfinite", True) and not math.isfinite(value):
                    raise RuntimeError(f"{path}: nonfinite {col} at data row {rows}")
            if (
                config["validation"].get("fail_on_nonpositive_I_when_rho_I_present", True)
                and "rho_I" in fieldnames
                and float(record["rho_I"]) <= 0.0
            ):
                raise RuntimeError(f"{path}: nonpositive rho_I at data row {rows}")
    if expected_rows and rows != expected_rows:
        raise RuntimeError(f"{path}: row count {rows} != {expected_rows}")
    if expected_vza:
        expected = {round(float(x), 10) for x in expected_vza}
        if vzas != expected:
            raise RuntimeError(f"{path}: VZA grid mismatch; got={sorted(vzas)} expected={sorted(expected)}")
    if expected_raa:
        expected = {round(float(x), 10) for x in expected_raa}
        if raas != expected:
            raise RuntimeError(f"{path}: RAA grid mismatch; got={sorted(raas)} expected={sorted(expected)}")
    return {"rows": rows, "vza_count": len(vzas), "raa_count": len(raas)}


def gzip_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_name(destination.name + ".tmp")
    with source.open("rb") as src, gzip.GzipFile(
        filename=str(temp), mode="wb", compresslevel=6, mtime=0
    ) as dst:
        shutil_copyfileobj(src, dst)
    os.replace(temp, destination)


def shutil_copyfileobj(src, dst, length: int = 1024 * 1024) -> None:
    while True:
        block = src.read(length)
        if not block:
            return
        dst.write(block)


def write_manifest_atomic(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in MANIFEST_FIELDS})
    os.replace(temp, path)


def run_one(
    row: dict[str, Any],
    binary: Path,
    croot: Path,
    output_dir: Path,
    logs_dir: Path,
    env: dict[str, str],
    config: dict[str, Any],
    keep_raw: bool,
) -> dict[str, Any]:
    run_id = row["run_id"]
    final = output_dir / row["output_relpath"]
    if final.suffix.lower() != ".gz":
        final = final.with_suffix(final.suffix + ".gz")
    final.parent.mkdir(parents=True, exist_ok=True)
    stdout_log = logs_dir / f"{run_id}.stdout.log"
    stderr_log = logs_dir / f"{run_id}.stderr.log"

    if final.is_file():
        try:
            info = validate_csv(
                final, row["expected_rows"], row["required_columns"],
                row["expected_vza"], row["expected_raa"], config,
            )
            return {
                "run_id": run_id, "stage": row["stage"], "status": "skipped_valid",
                "returncode": 0, "elapsed_s": 0.0, "rows": info["rows"],
                "output_relpath": str(final.relative_to(output_dir)),
                "output_bytes": final.stat().st_size, "output_sha256": sha256_file(final),
                "binary_sha256": sha256_file(binary), "command": "",
                "stdout_log": str(stdout_log.relative_to(output_dir)),
                "stderr_log": str(stderr_log.relative_to(output_dir)),
                "message": "existing validated output reused",
            }
        except Exception:
            quarantine = final.with_name(final.name + f".invalid_{int(time.time())}")
            os.replace(final, quarantine)

    work_dir = output_dir / ".work"
    work_dir.mkdir(parents=True, exist_ok=True)
    raw = work_dir / f"{run_id}.csv"
    raw.unlink(missing_ok=True)
    args = list(row["args"])
    cmd = [str(binary), *args, "--output-full-grid", str(raw)]
    command_text = shlex.join(cmd)
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=croot, env=env, text=True, capture_output=True)
    elapsed = time.time() - t0
    stdout_log.write_text(proc.stdout or "", encoding="utf-8")
    stderr_log.write_text(proc.stderr or "", encoding="utf-8")
    if proc.returncode != 0:
        return {
            "run_id": run_id, "stage": row["stage"], "status": "failed",
            "returncode": proc.returncode, "elapsed_s": elapsed, "rows": 0,
            "output_relpath": str(final.relative_to(output_dir)), "output_bytes": 0,
            "output_sha256": "", "binary_sha256": sha256_file(binary),
            "command": command_text,
            "stdout_log": str(stdout_log.relative_to(output_dir)),
            "stderr_log": str(stderr_log.relative_to(output_dir)),
            "message": "OCRT returned nonzero",
        }
    try:
        info = validate_csv(
            raw, row["expected_rows"], row["required_columns"],
            row["expected_vza"], row["expected_raa"], config,
        )
        gzip_atomic(raw, final)
        # Re-open the compressed authoritative artifact and validate once more.
        validate_csv(
            final, row["expected_rows"], row["required_columns"],
            row["expected_vza"], row["expected_raa"], config,
        )
        if not keep_raw:
            raw.unlink(missing_ok=True)
        return {
            "run_id": run_id, "stage": row["stage"], "status": "ok",
            "returncode": 0, "elapsed_s": elapsed, "rows": info["rows"],
            "output_relpath": str(final.relative_to(output_dir)),
            "output_bytes": final.stat().st_size, "output_sha256": sha256_file(final),
            "binary_sha256": sha256_file(binary), "command": command_text,
            "stdout_log": str(stdout_log.relative_to(output_dir)),
            "stderr_log": str(stderr_log.relative_to(output_dir)), "message": "",
        }
    except Exception as exc:
        return {
            "run_id": run_id, "stage": row["stage"], "status": "failed_validation",
            "returncode": 0, "elapsed_s": elapsed, "rows": 0,
            "output_relpath": str(final.relative_to(output_dir)), "output_bytes": 0,
            "output_sha256": "", "binary_sha256": sha256_file(binary),
            "command": command_text,
            "stdout_log": str(stdout_log.relative_to(output_dir)),
            "stderr_log": str(stderr_log.relative_to(output_dir)),
            "message": str(exc),
        }


def make_pilot_estimate(results: list[dict[str, Any]], matrix_count: int) -> dict[str, Any]:
    ok = [r for r in results if r["status"] in {"ok", "skipped_valid"}]
    if not ok:
        return {"status": "NO_VALID_RUNS"}
    active = [r for r in ok if float(r["elapsed_s"]) > 0]
    elapsed_mean = sum(float(r["elapsed_s"]) for r in active) / len(active) if active else 0.0
    bytes_mean = sum(int(r["output_bytes"]) for r in ok) / len(ok)
    return {
        "status": "ESTIMATE",
        "sample_runs": len(ok),
        "matrix_runs": matrix_count,
        "mean_compressed_bytes_per_run": bytes_mean,
        "estimated_compressed_output_bytes": bytes_mean * matrix_count,
        "mean_single_run_wall_s": elapsed_mean,
        "estimated_16_worker_wall_s_ideal": (elapsed_mean * matrix_count / 16.0) if elapsed_mean else None,
        "warning": "Wall-time estimate is idealized; I/O, case complexity, and convergence vary by state.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run OCRT full-grid cases with PSSA ON and 16 independent workers.")
    parser.add_argument("--runtime-root", type=Path, required=True, help="MIGRATION_PKG_2026-08-19 directory")
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--exe", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--pilot", type=int, default=0, help="Run only first N rows and estimate storage/time")
    parser.add_argument("--max-runs", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep-raw", action="store_true")
    args = parser.parse_args()

    config = load_json(args.config.resolve())
    runtime_root = args.runtime_root.resolve()
    croot = runtime_root / "01_OCRT_C"
    binary = find_binary(runtime_root, args.exe)
    matrix = read_matrix(args.matrix.resolve())
    for row in matrix:
        validate_row_policy(row, config)
    if args.pilot:
        matrix = matrix[: args.pilot]
    if args.max_runs:
        matrix = matrix[: args.max_runs]
    workers = args.workers or int(config["workers"])
    if workers <= 0 or workers > 16:
        raise RuntimeError("this package policy permits 1..16 workers; default is 16")
    env = build_child_env(config)
    output_dir = args.output_dir.resolve()
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    run_config = {
        "campaign_config": config,
        "matrix_file": str(args.matrix.resolve()),
        "matrix_sha256": sha256_file(args.matrix.resolve()),
        "selected_rows": len(matrix),
        "workers": workers,
        "binary": str(binary),
        "binary_sha256": sha256_file(binary),
        "c_root": str(croot),
        "actual_child_environment": {k: env[k] for k in config["child_environment"]},
        "forbidden_environment_checked": config.get("forbidden_environment_variables", []),
        "pssa": "ON (explicit --pssa in every actual command)",
        "water_truncation": "OFF (forbidden options audited)",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_config.json").write_text(
        json.dumps(run_config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    commands = []
    for row in matrix:
        final = output_dir / row["output_relpath"]
        raw = output_dir / ".work" / f"{row['run_id']}.csv"
        commands.append(shlex.join([str(binary), *row["args"], "--output-full-grid", str(raw)]))
    (output_dir / "run_commands.txt").write_text("\n".join(commands) + "\n", encoding="utf-8")
    if args.dry_run:
        print(json.dumps(run_config, ensure_ascii=False, indent=2))
        print(f"dry-run rows={len(matrix)}; commands={output_dir / 'run_commands.txt'}")
        return 0

    results: list[dict[str, Any]] = []
    manifest_path = output_dir / "run_manifest.csv"
    lock = threading.Lock()
    print(f"starting {len(matrix)} runs with workers={workers}, OMP_NUM_THREADS=1, PSSA=ON", flush=True)
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="ocrt") as executor:
        futures = {
            executor.submit(
                run_one, row, binary, croot, output_dir, logs_dir, env, config, args.keep_raw
            ): row
            for row in matrix
        }
        for future in as_completed(futures):
            row = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {
                    "run_id": row["run_id"], "stage": row["stage"], "status": "runner_exception",
                    "returncode": -1, "elapsed_s": 0, "rows": 0,
                    "output_relpath": row["output_relpath"], "output_bytes": 0,
                    "output_sha256": "", "binary_sha256": sha256_file(binary),
                    "command": "", "stdout_log": "", "stderr_log": "", "message": repr(exc),
                }
            with lock:
                results.append(result)
                write_manifest_atomic(manifest_path, results)
            print(
                f"[{len(results)}/{len(matrix)}] {result['run_id']} {result['status']} "
                f"{float(result['elapsed_s']):.2f}s {result.get('message','')}",
                flush=True,
            )

    results.sort(key=lambda r: r["run_id"])
    write_manifest_atomic(manifest_path, results)
    estimate = make_pilot_estimate(results, len(read_matrix(args.matrix.resolve())))
    (output_dir / "PILOT_OR_RUN_ESTIMATE.json").write_text(
        json.dumps(estimate, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    failures = [r for r in results if r["status"] not in {"ok", "skipped_valid"}]
    summary = {
        "selected_runs": len(results),
        "ok_or_reused": len(results) - len(failures),
        "failures": len(failures),
        "workers": workers,
        "pssa": True,
        "broad_truncation": False,
        "output_bytes": sum(int(r.get("output_bytes") or 0) for r in results),
        "status": "PASS" if not failures else "FAIL",
    }
    (output_dir / "RUN_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
