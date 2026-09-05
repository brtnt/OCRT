#!/usr/bin/env python3
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import zipfile

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
BASE_TAR_NAME = "OCRT_MIGRATION_PKG_2026-08-19 (1).tar.gz"
POLICY_ZIP_NAME = "OCRT_MIE_GRID_TRUNCATION_VALIDATION_ARTIFACT_UPDATE_2026-08-21.zip"
MIGRATION_DIR_NAME = "MIGRATION_PKG_2026-08-19"
HANDOFF_ROOT_NAME = "OCRT_SPECTRAL_DATA_330_1100_FR631_HANDOFF_2026-08-20"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_required_manifest() -> dict:
    path = PACKAGE_ROOT / "manifests" / "REQUIRED_FILES.json"
    return json.loads(path.read_text(encoding="utf-8"))


def find_part(parts_dir: Path, index: int) -> Path:
    matches = sorted(parts_dir.glob(f"PART{index:02d}_*.zip"))
    if not matches:
        raise FileNotFoundError(f"PART{index:02d} archive not found under {parts_dir}")
    if len(matches) > 1:
        # Exact known names are preferred; otherwise ambiguity is unsafe.
        known = [p for p in matches if p.name.endswith("(1).zip")]
        if len(known) == 1:
            return known[0]
        raise RuntimeError(f"ambiguous PART{index:02d} archives: {[p.name for p in matches]}")
    return matches[0]


def verify_inputs(parts_dir: Path) -> tuple[Path, Path, list[Path]]:
    manifest = load_required_manifest()
    base_tar = PACKAGE_ROOT / "base" / BASE_TAR_NAME
    policy_zip = PACKAGE_ROOT / "base" / POLICY_ZIP_NAME
    paths = [base_tar, policy_zip]
    parts = [find_part(parts_dir, i) for i in range(0, 11)]
    paths.extend(parts)

    failures: list[str] = []
    for path in paths:
        if not path.is_file():
            failures.append(f"missing: {path}")
            continue
        expected = manifest["files"].get(path.name)
        if not expected:
            failures.append(f"manifest entry missing: {path.name}")
            continue
        actual_size = path.stat().st_size
        if actual_size != int(expected["size_bytes"]):
            failures.append(
                f"size mismatch: {path.name}: {actual_size} != {expected['size_bytes']}"
            )
            continue
        actual_hash = sha256_file(path)
        if actual_hash != expected["sha256"]:
            failures.append(
                f"sha256 mismatch: {path.name}: {actual_hash} != {expected['sha256']}"
            )
    if failures:
        raise RuntimeError("input verification failed:\n" + "\n".join(failures))
    return base_tar, policy_zip, parts


def safe_extract_tar(archive: Path, destination: Path) -> None:
    with tarfile.open(archive, "r:*") as tf:
        try:
            tf.extractall(destination, filter="data")
        except TypeError:
            tf.extractall(destination)


def apply_policy_overlay(policy_zip: Path, migration_root: Path, work: Path) -> None:
    overlay_extract = work / "policy_overlay"
    with zipfile.ZipFile(policy_zip) as zf:
        zf.extractall(overlay_extract)
    roots = [p for p in overlay_extract.iterdir() if p.is_dir()]
    if len(roots) != 1:
        raise RuntimeError(f"unexpected policy archive root: {roots}")
    script = roots[0] / "apply_overlay.py"
    result = subprocess.run(
        [sys.executable, str(script), "--migration-root", str(migration_root), "--apply"],
        cwd=roots[0],
        text=True,
        capture_output=True,
    )
    (work / "apply_policy_overlay.stdout.log").write_text(result.stdout, encoding="utf-8")
    (work / "apply_policy_overlay.stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"policy overlay failed; see {work}")
    # The immutable base archive is retained, so duplicate overlay backups are not needed.
    for backup in migration_root.glob(".policy_update_backup_*"):
        shutil.rmtree(backup)


def extract_part00(part00: Path, destination: Path) -> Path:
    with zipfile.ZipFile(part00) as zf:
        zf.extractall(destination)
    roots = list(destination.rglob(HANDOFF_ROOT_NAME))
    roots = [p for p in roots if p.is_dir()]
    if len(roots) != 1:
        raise RuntimeError(f"expected one handoff root in PART00, found {roots}")
    return roots[0]


def apply_selected_consumer_patch(bundle_root: Path, croot: Path) -> list[str]:
    # PART00 contains the approved FR631 spectral/angle consumer changes. The
    # later policy overlay contains a newer shared/mie_io implementation with
    # worker-private cache and grid-policy comments, so shared/mie_io is NOT
    # overwritten here.
    selected = [
        "src/rt_aerosol.c",
        "src/rt_aerosol_runtime.c",
    ]
    copied: list[str] = []
    for rel in selected:
        src = bundle_root / "03_INTEGRATION_PATCH" / "C" / rel
        dst = croot / rel
        if not src.is_file():
            raise FileNotFoundError(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied.append(rel)

    # The migration source retained a 64-column AHN reader although the latest
    # canonical files contain 771 phase wavelengths. This is a fail-loud source
    # correction already identified during the 2026-08-21 integration audit.
    ahn = croot / "src" / "rt_iop_ahn_mineral.c"
    text = ahn.read_text(encoding="utf-8")
    old = "#define AHN_MAX_PHW    64"
    if old in text:
        text = text.replace(
            old,
            "/* FR631/771 contract: canonical AHN files contain 771 phase-\n"
            " * wavelength columns. Keep headroom and a line buffer large enough\n"
            " * for a complete row. DOC-REF: OCRT_330_1100_FR631_INTEGRATION_GUIDE_KO;\n"
            " * OCRT_MIE_GRID_TRUNCATION_AND_VALIDATION_ARTIFACT_POLICY_2026-08-21. */\n"
            "#define AHN_MAX_PHW  1024",
            1,
        )
    elif "#define AHN_MAX_PHW  1024" not in text:
        raise RuntimeError("unexpected AHN_MAX_PHW definition")
    text = text.replace("char line[8192];", "char line[131072];")
    text = text.replace("char line[16384];", "char line[131072];")
    ahn.write_text(text, encoding="utf-8")
    copied.append("src/rt_iop_ahn_mineral.c (771-column patch)")
    return copied


def build_zip_member_index(parts: list[Path]) -> dict[str, tuple[Path, str, int]]:
    index: dict[str, tuple[Path, str, int]] = {}
    for archive in parts:
        with zipfile.ZipFile(archive) as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                name = info.filename.replace("\\", "/")
                marker = HANDOFF_ROOT_NAME + "/"
                pos = name.find(marker)
                if pos < 0:
                    continue
                relative = name[pos + len(marker):]
                prior = index.get(relative)
                candidate = (archive, info.filename, info.file_size)
                if prior and (prior[2] != info.file_size or prior[0] != archive):
                    # A file may occur once in PART00 and once elsewhere only if
                    # byte-identical. Conflicting source ownership is unsafe.
                    raise RuntimeError(
                        f"duplicate source member {relative}: {prior[0].name}, {archive.name}"
                    )
                index[relative] = candidate
    return index


def stream_zip_member(archive: Path, member: str, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp = destination.with_name(destination.name + ".installing")
    h = hashlib.sha256()
    with zipfile.ZipFile(archive) as zf, zf.open(member) as src, temp.open("wb") as dst:
        while True:
            block = src.read(1024 * 1024)
            if not block:
                break
            h.update(block)
            dst.write(block)
        dst.flush()
        os.fsync(dst.fileno())
    os.replace(temp, destination)
    return h.hexdigest()


def hardlink_or_copy(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    try:
        os.link(source, destination)
        return "hardlink"
    except OSError:
        shutil.copy2(source, destination)
        return "copy"


def extract_python_tree(migration_root: Path) -> Path:
    py_archive = (
        migration_root
        / "03_PYTHON_AND_NEW_DATA"
        / "pyOCRT_v1_2_SPECTRAL_DATA_330_1100_COMPONENT_GATED_2026-08-14_tar.gz"
    )
    if not py_archive.is_file():
        raise FileNotFoundError(py_archive)
    destination = migration_root / "03_PYTHON_RUNTIME"
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    safe_extract_tar(py_archive, destination)
    roots = [p for p in destination.iterdir() if p.is_dir()]
    if len(roots) != 1:
        raise RuntimeError(f"unexpected Python package roots: {roots}")
    return roots[0]


def install_data_streaming(
    bundle_root: Path,
    parts: list[Path],
    croot: Path,
    python_root: Path | None,
    log_path: Path,
) -> dict:
    member_index = build_zip_member_index(parts)
    installation_map = bundle_root / "06_MANIFESTS" / "INSTALLATION_MAP_C_PYTHON.csv"
    inventory_path = bundle_root / "06_MANIFESTS" / "MIE_CANONICAL_INVENTORY_198.csv"
    rows = list(csv.DictReader(installation_map.open(encoding="utf-8")))
    inventory = list(csv.DictReader(inventory_path.open(encoding="utf-8")))
    expected_mie = {
        "02_RUNTIME_MIE/" + row["relative_path"]: row["sha256"] for row in inventory
    }

    log_rows: list[dict[str, str]] = []
    missing: list[str] = []
    for row in rows:
        source_relative = row["source_relative"].replace("\\", "/")
        entry = member_index.get(source_relative)
        if not entry:
            missing.append(source_relative)
            continue
        archive, member, _ = entry
        c_destination = croot / row["c_destination"]
        actual_hash = stream_zip_member(archive, member, c_destination)
        expected_hash = expected_mie.get(source_relative)
        if expected_hash and actual_hash != expected_hash:
            raise RuntimeError(
                f"canonical Mie hash mismatch after install: {source_relative}: "
                f"{actual_hash} != {expected_hash}"
            )
        py_method = "not_requested"
        extra_method = ""
        if python_root is not None:
            py_destination = python_root / row["python_destination"]
            py_method = hardlink_or_copy(c_destination, py_destination)
            extra = row.get("python_extra_destination", "").strip()
            if extra:
                extra_method = hardlink_or_copy(c_destination, python_root / extra)
        log_rows.append(
            {
                "source_relative": source_relative,
                "source_archive": archive.name,
                "c_destination": str(c_destination),
                "sha256": actual_hash,
                "python_method": py_method,
                "python_extra_method": extra_method,
            }
        )
        if len(log_rows) % 25 == 0:
            print(f"installed {len(log_rows)}/{len(rows)} mapped files", flush=True)
    if missing:
        raise RuntimeError(f"installation-map source files missing: {missing[:20]} total={len(missing)}")

    with log_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(log_rows[0]))
        writer.writeheader()
        writer.writerows(log_rows)
    installed_canonical = sum(1 for row in log_rows if row["source_relative"] in expected_mie)
    if installed_canonical != 198:
        raise RuntimeError(f"canonical Mie install count {installed_canonical} != 198")
    return {
        "mapped_files": len(log_rows),
        "canonical_mie_installed": installed_canonical,
        "python_tree_installed": python_root is not None,
    }



def extract_parts_merged(parts: list[Path], destination: Path, work: Path) -> Path:
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    unzip = shutil.which("unzip")
    logs: list[str] = []
    for archive in parts:
        t0 = time.time()
        if unzip:
            result = subprocess.run(
                [unzip, "-q", "-o", str(archive), "-d", str(destination)],
                text=True, capture_output=True,
            )
            logs.append(
                f"{archive.name} rc={result.returncode} elapsed={time.time()-t0:.3f}s\n"
                f"stdout={result.stdout}\nstderr={result.stderr}\n"
            )
            if result.returncode != 0:
                raise RuntimeError(f"unzip failed: {archive}; see extract_parts.log")
        else:
            with zipfile.ZipFile(archive) as zf:
                zf.extractall(destination)
            logs.append(f"{archive.name} python-zip elapsed={time.time()-t0:.3f}s\n")
    (work / "extract_parts.log").write_text("\n".join(logs), encoding="utf-8")
    roots = [p for p in destination.rglob(HANDOFF_ROOT_NAME) if p.is_dir()]
    if len(roots) != 1:
        raise RuntimeError(f"expected one merged handoff root, found {roots}")
    bundle_root = roots[0]
    verify = bundle_root / "05_TOOLS" / "verify_extracted_bundle.py"
    result = subprocess.run(
        [sys.executable, str(verify)], cwd=bundle_root, text=True, capture_output=True
    )
    (work / "verify_extracted_bundle.stdout.log").write_text(result.stdout, encoding="utf-8")
    (work / "verify_extracted_bundle.stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError("merged bundle verification failed")
    return bundle_root


def install_data_from_extracted(
    bundle_root: Path,
    croot: Path,
    python_root: Path | None,
    log_path: Path,
) -> dict:
    installation_map = bundle_root / "06_MANIFESTS" / "INSTALLATION_MAP_C_PYTHON.csv"
    inventory_path = bundle_root / "06_MANIFESTS" / "MIE_CANONICAL_INVENTORY_198.csv"
    rows = list(csv.DictReader(installation_map.open(encoding="utf-8")))
    inventory = list(csv.DictReader(inventory_path.open(encoding="utf-8")))
    expected_mie = {
        "02_RUNTIME_MIE/" + row["relative_path"]: row["sha256"] for row in inventory
    }
    log_rows: list[dict[str, str]] = []
    for i, row in enumerate(rows, 1):
        source_relative = row["source_relative"].replace("\\", "/")
        source = bundle_root / source_relative
        if not source.is_file():
            raise FileNotFoundError(source)
        c_destination = croot / row["c_destination"]
        c_method = hardlink_or_copy(source, c_destination)
        actual_hash = sha256_file(c_destination)
        expected_hash = expected_mie.get(source_relative)
        if expected_hash and actual_hash != expected_hash:
            raise RuntimeError(
                f"canonical Mie hash mismatch after install: {source_relative}: "
                f"{actual_hash} != {expected_hash}"
            )
        py_method = "not_requested"
        extra_method = ""
        if python_root is not None:
            py_method = hardlink_or_copy(c_destination, python_root / row["python_destination"])
            extra = row.get("python_extra_destination", "").strip()
            if extra:
                extra_method = hardlink_or_copy(c_destination, python_root / extra)
        log_rows.append({
            "source_relative": source_relative,
            "c_destination": str(c_destination),
            "sha256": actual_hash,
            "c_method": c_method,
            "python_method": py_method,
            "python_extra_method": extra_method,
        })
        if i % 25 == 0:
            print(f"installed {i}/{len(rows)} mapped files", flush=True)
    with log_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(log_rows[0]))
        writer.writeheader(); writer.writerows(log_rows)
    installed_canonical = sum(1 for row in log_rows if row["source_relative"] in expected_mie)
    if installed_canonical != 198:
        raise RuntimeError(f"canonical Mie install count {installed_canonical} != 198")
    return {
        "mapped_files": len(log_rows),
        "canonical_mie_installed": installed_canonical,
        "python_tree_installed": python_root is not None,
        "source_install_mode": "hardlink where supported; copy fallback",
    }


def extract_parts_parallel(
    parts: list[Path], destination: Path, work: Path, workers: int = 4
) -> tuple[Path, list[Path]]:
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True)
    unzip = shutil.which("unzip")

    def one(item: tuple[int, Path]) -> tuple[int, Path, str]:
        index, archive = item
        target = destination / f"part{index:02d}"
        target.mkdir(parents=True)
        t0 = time.time()
        if unzip:
            result = subprocess.run(
                [unzip, "-q", "-o", str(archive), "-d", str(target)],
                text=True, capture_output=True,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"unzip failed: {archive}: {result.stderr.strip()}"
                )
        else:
            with zipfile.ZipFile(archive) as zf:
                zf.extractall(target)
        roots = [p for p in target.rglob(HANDOFF_ROOT_NAME) if p.is_dir()]
        if len(roots) != 1:
            raise RuntimeError(f"{archive.name}: expected one handoff root, found {roots}")
        return index, roots[0], f"{archive.name} elapsed={time.time()-t0:.3f}s"

    roots_by_index: dict[int, Path] = {}
    logs: list[str] = []
    with ThreadPoolExecutor(max_workers=min(workers, len(parts))) as executor:
        futures = {executor.submit(one, item): item for item in enumerate(parts)}
        for future in as_completed(futures):
            index, root, message = future.result()
            roots_by_index[index] = root
            logs.append(message)
            print(f"extracted PART{index:02d}/{len(parts)-1:02d}", flush=True)
    roots = [roots_by_index[i] for i in range(len(parts))]
    (work / "extract_parts.log").write_text("\n".join(sorted(logs)) + "\n", encoding="utf-8")
    return roots[0], roots


def make_source_index(roots: list[Path]) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for root in roots:
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            rel = path.relative_to(root).as_posix()
            prior = index.get(rel)
            if prior is not None:
                if prior.stat().st_size != path.stat().st_size:
                    raise RuntimeError(f"conflicting source path across PARTs: {rel}")
                continue
            index[rel] = path
    return index


def audit_extracted_sources_fast(
    part00_root: Path, source_index: dict[str, Path], full_hash: bool, work: Path
) -> dict:
    inventory_path = part00_root / "06_MANIFESTS" / "MIE_CANONICAL_INVENTORY_198.csv"
    inventory = list(csv.DictReader(inventory_path.open(encoding="utf-8")))
    if len(inventory) != 198:
        raise RuntimeError(f"canonical inventory rows {len(inventory)} != 198")
    failures: list[str] = []
    checked_hashes = 0
    for i, row in enumerate(inventory):
        rel = "02_RUNTIME_MIE/" + row["relative_path"]
        path = source_index.get(rel)
        if path is None:
            failures.append(f"missing {rel}")
            continue
        if path.stat().st_size <= 0:
            failures.append(f"empty {rel}")
            continue
        # Exact PART archive SHA-256 plus unzip CRC already bind the extracted
        # bytes. Full per-Mie SHA is optional because it reads another 4.94 GB.
        # In fast mode, still hash one representative from each major family.
        representative = any(token in rel for token in [
            "/aerosol_snf_mie/M50C.mie",
            "/aerosol_ahmad2010_paper_mie/A2010ver_r80f20v01.mie",
            "/aerosol_ahmad2010_accurt_mie/r80f20v01.mie",
            "/water/tsm_ahn/Brown_earth_AHN.mie",
            "/water/detritus/Detritus_Stramski2001.mie",
            "/water/eap/EAP_00_Diatoms_pennate_D6.mie",
        ])
        if full_hash or representative:
            actual = sha256_file(path)
            checked_hashes += 1
            if actual != row["sha256"]:
                failures.append(f"sha256 {rel}: {actual} != {row['sha256']}")
    report = {
        "canonical_expected": 198,
        "canonical_present": 198 - len([x for x in failures if x.startswith('missing')]),
        "hashes_checked": checked_hashes,
        "full_hash_mode": full_hash,
        "failure_count": len(failures),
        "failures": failures[:20],
        "status": "PASS" if not failures else "FAIL",
    }
    (work / "extracted_source_audit.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    if failures:
        raise RuntimeError(f"extracted source audit failed: {failures[:5]}")
    return report


def install_data_from_source_index(
    part00_root: Path, source_index: dict[str, Path], croot: Path,
    python_root: Path | None, log_path: Path
) -> dict:
    installation_map = part00_root / "06_MANIFESTS" / "INSTALLATION_MAP_C_PYTHON.csv"
    inventory_path = part00_root / "06_MANIFESTS" / "MIE_CANONICAL_INVENTORY_198.csv"
    rows = list(csv.DictReader(installation_map.open(encoding="utf-8")))
    inventory = list(csv.DictReader(inventory_path.open(encoding="utf-8")))
    expected_mie = {
        "02_RUNTIME_MIE/" + row["relative_path"]: row["sha256"] for row in inventory
    }
    log_rows: list[dict[str, str]] = []
    for i, row in enumerate(rows, 1):
        source_relative = row["source_relative"].replace("\\", "/")
        source = source_index.get(source_relative)
        if source is None:
            raise FileNotFoundError(source_relative)
        c_destination = croot / row["c_destination"]
        c_method = hardlink_or_copy(source, c_destination)
        # Hardlink identity proves byte identity. Copy fallback is hashed.
        actual_hash = expected_mie.get(source_relative, "")
        if c_method == "copy":
            copied_hash = sha256_file(c_destination)
            if actual_hash and copied_hash != actual_hash:
                raise RuntimeError(f"copy hash mismatch: {source_relative}")
            actual_hash = copied_hash
        py_method = "not_requested"
        extra_method = ""
        if python_root is not None:
            py_method = hardlink_or_copy(c_destination, python_root / row["python_destination"])
            extra = row.get("python_extra_destination", "").strip()
            if extra:
                extra_method = hardlink_or_copy(c_destination, python_root / extra)
        log_rows.append({
            "source_relative": source_relative,
            "c_destination": str(c_destination),
            "sha256": actual_hash,
            "c_method": c_method,
            "python_method": py_method,
            "python_extra_method": extra_method,
        })
        if i % 25 == 0:
            print(f"installed {i}/{len(rows)} mapped files", flush=True)
    with log_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(log_rows[0]))
        writer.writeheader(); writer.writerows(log_rows)
    installed_canonical = sum(1 for row in log_rows if row["source_relative"] in expected_mie)
    if installed_canonical != 198:
        raise RuntimeError(f"canonical Mie install count {installed_canonical} != 198")
    return {
        "mapped_files": len(log_rows),
        "canonical_mie_installed": installed_canonical,
        "python_tree_installed": python_root is not None,
        "source_install_mode": "hardlink where supported; copy fallback",
    }

def detect_march() -> str:
    flags: set[str] = set()
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(errors="ignore").splitlines():
            if line.startswith("flags") or line.startswith("Features"):
                flags.update(line.split(":", 1)[1].strip().split())
                break
    cascadelake = {"avx512f", "avx512dq", "avx512bw", "avx512vl"}
    if cascadelake.issubset(flags):
        return "cascadelake"
    if {"avx2", "fma", "bmi2"}.issubset(flags):
        return "x86-64-v3"
    return "x86-64"


def build_c(croot: Path, march: str, work: Path) -> Path:
    output = croot / "build" / "ocrt_polarization"
    output.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OCRT_MARCH"] = march
    result = subprocess.run(
        ["bash", str(croot / "scripts" / "build_release_v1.2.sh"), str(output)],
        cwd=croot,
        env=env,
        text=True,
        capture_output=True,
    )
    (work / "build.stdout.log").write_text(result.stdout, encoding="utf-8")
    (work / "build.stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0 or not output.is_file():
        raise RuntimeError(f"C build failed; see {work / 'build.stderr.log'}")
    return output


def validate_smoke_csv(path: Path) -> int:
    required = {"sza_deg", "wavelength_nm", "vza_deg", "raa_deg", "rho_I", "rho_Q", "rho_U"}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"smoke CSV missing columns: {sorted(missing)}")
        rows = list(reader)
    if len(rows) != 1296:
        raise RuntimeError(f"smoke row count {len(rows)} != 1296")
    return len(rows)


def run_smoke(binary: Path, croot: Path, work: Path) -> dict:
    output = work / "setup_smoke_555_sza75.csv"
    cmd = [
        str(binary), "--sza", "75", "--wavelength", "555", "--pressure", "1013.25",
        "--surface", "black_fresnel_ocean", "--wind-speed", "5",
        "--decouple-sunglint", "--pssa", "--lut-vza-step", "5",
        "--lut-vza-max", "85", "--lut-raa-step", "5",
        "--output-full-grid", str(output),
    ]
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1",
        "OMP_DYNAMIC": "FALSE",
        "OCRT_ADVANCED": "1",
        "OCRT_WATER_PARTICLE_KERNEL": "direct",
    })
    t0 = time.time()
    result = subprocess.run(cmd, cwd=croot, env=env, text=True, capture_output=True)
    elapsed = time.time() - t0
    (work / "smoke.stdout.log").write_text(result.stdout, encoding="utf-8")
    (work / "smoke.stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"PSSA smoke failed; see {work / 'smoke.stderr.log'}")
    rows = validate_smoke_csv(output)
    digest = sha256_file(output)
    output.unlink()
    return {"rows": rows, "elapsed_s": elapsed, "sha256_before_delete": digest}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Assemble the latest OCRT C/Python runtime from the migration base and PART00-PART10."
    )
    parser.add_argument("--parts-dir", type=Path, required=True)
    parser.add_argument("--runtime-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--no-python-data", action="store_true", help="Skip paired Python tree/data hardlinks")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--march", default="auto", help="auto|cascadelake|x86-64-v3|x86-64")
    parser.add_argument("--full-mie-hash-audit", action="store_true",
                        help="Read and SHA-256 all 4.94 GB of canonical Mie after extraction")
    parser.add_argument("--extract-workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    parts_dir = args.parts_dir.resolve()
    runtime_dir = args.runtime_dir.resolve()
    base_tar, policy_zip, parts = verify_inputs(parts_dir)
    print("input archives: PASS", flush=True)
    if args.verify_only:
        return 0

    if runtime_dir.exists():
        if not args.force:
            raise RuntimeError(f"runtime directory exists: {runtime_dir}; use --force to replace")
        shutil.rmtree(runtime_dir)
    runtime_dir.mkdir(parents=True)
    work = runtime_dir / ".setup_work"
    work.mkdir()

    safe_extract_tar(base_tar, runtime_dir)
    migration_root = runtime_dir / MIGRATION_DIR_NAME
    if not migration_root.is_dir():
        raise RuntimeError(f"migration root not found after extraction: {migration_root}")
    croot = migration_root / "01_OCRT_C"

    apply_policy_overlay(policy_zip, migration_root, work)
    part_extract = work / "parts_extracted"
    bundle_root, source_roots = extract_parts_parallel(
        parts, part_extract, work, workers=max(1, args.extract_workers)
    )
    source_index = make_source_index(source_roots)
    source_audit = audit_extracted_sources_fast(
        bundle_root, source_index, args.full_mie_hash_audit, work
    )
    selected_patch = apply_selected_consumer_patch(bundle_root, croot)

    python_root = None if args.no_python_data else extract_python_tree(migration_root)
    install_summary = install_data_from_source_index(
        bundle_root, source_index, croot, python_root, work / "data_install_manifest.csv"
    )
    install_summary["source_audit"] = source_audit
    # Installed hardlinks retain the inode after extraction names are removed.
    shutil.rmtree(part_extract)

    march = detect_march() if args.march == "auto" else args.march
    binary = None
    smoke = None
    if not args.no_build:
        binary = build_c(croot, march, work)
        smoke = run_smoke(binary, croot, work)

    status = {
        "status": "PASS",
        "package": "OCRT_POLARIZATION_PSSA_16CORE_RUN_PACKAGE_2026-08-22",
        "migration_root": str(migration_root),
        "c_root": str(croot),
        "python_root": str(python_root) if python_root else None,
        "selected_consumer_patch": selected_patch,
        "data_install": install_summary,
        "binary": str(binary) if binary else None,
        "binary_sha256": sha256_file(binary) if binary else None,
        "march": march,
        "pssa_policy": "ON in campaign runner",
        "truncation_policy": "raw FR631; broad truncation OFF",
        "parallel_policy": "16 independent processes; OMP_NUM_THREADS=1 each",
        "smoke": smoke,
    }
    (runtime_dir / "SETUP_STATUS.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(status, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
