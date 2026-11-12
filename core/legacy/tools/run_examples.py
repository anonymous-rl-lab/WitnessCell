#!/usr/bin/env python3
"""Verify or rebuild saved RNA/Jurkat response displays without model fitting."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import stat
import subprocess
import sys
import time
import zipfile
from pathlib import Path, PurePosixPath

VERSION = "25.0.0"
EXPERIMENT = "experiments/28_rna_jurkat_response_examples"
ARCHIVE_NAME = "WitnessCell_RNA_Jurkat_Response_Examples_20261019.zip"
ARCHIVE_ROOT = "WitnessCell_RNA_Jurkat_Response_Examples_20261019"
ARCHIVE_BYTES = 7679888
ARCHIVE_SHA256 = "7f9f63bc812df611955a824b2573066d893c5c6d63b8bfc6294e550db562ec41"
REPORT_NAME = "RESPONSE_EXAMPLES_REPRODUCTION_CHECK.json"
FIGURE_NAMES = (
    "RNA_four_corner_responses", "RNA_all_contrasts_and_limitation",
    "Jurkat_specific_combinations", "Jurkat_all_pairs_and_limitation",
)


def digest_file(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def safe_member(name):
    """Require canonical relative POSIX paths; never accept archive links."""
    if not isinstance(name, str) or not name:
        raise ValueError("Empty/non-string ZIP or manifest path")
    stripped = name[:-1] if name.endswith("/") else name
    p = PurePosixPath(stripped)
    if (not stripped or "\\" in name or "\x00" in name or ":" in name
            or p.is_absolute() or ".." in p.parts or str(p) != stripped
            or str(p) == "."):
        raise ValueError("Unsafe/noncanonical archive path: " + repr(name))
    return p


def inside(path, parent):
    # Path.is_relative_to is unavailable on Python 3.8.
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def validate_archive(repo):
    archive = Path(repo) / EXPERIMENT / ARCHIVE_NAME
    if archive.is_symlink() or not archive.is_file():
        raise ValueError("Missing ordinary response-examples archive")
    if archive.stat().st_size != ARCHIVE_BYTES:
        raise ValueError("Response-examples archive byte count differs")
    if digest_file(archive) != ARCHIVE_SHA256:
        raise ValueError("Response-examples archive SHA-256 differs")

    with zipfile.ZipFile(archive) as z:
        entries = {}
        for info in z.infolist():
            p = safe_member(info.filename)
            if p.parts[0] != ARCHIVE_ROOT:
                raise ValueError("Unexpected ZIP root")
            key = str(p)
            if key in entries:
                raise ValueError("Duplicate ZIP member: " + key)
            kind = stat.S_IFMT(info.external_attr >> 16)
            if kind not in (0, stat.S_IFREG, stat.S_IFDIR):
                raise ValueError("ZIP links and special files are forbidden")
            if info.flag_bits & 1:
                raise ValueError("Encrypted ZIP members are unsupported")
            entries[key] = info
        manifest_name = ARCHIVE_ROOT + "/MANIFEST.json"
        manifest = json.loads(z.read(manifest_name).decode("utf-8"))
        if manifest.get("package") != ARCHIVE_ROOT:
            raise ValueError("Internal manifest package identity differs")
        records = manifest.get("files")
        if not isinstance(records, list):
            raise ValueError("Internal manifest has no file list")
        listed = set()
        for record in records:
            rel = safe_member(record["path"])
            name = ARCHIVE_ROOT + "/" + str(rel)
            if name in listed or name == manifest_name:
                raise ValueError("Duplicate/self-referential manifest entry")
            listed.add(name)
            info = entries.get(name)
            if info is None or info.is_dir():
                raise ValueError("Missing manifest file: " + name)
            if info.file_size != record["bytes"]:
                raise ValueError("Manifest byte count differs: " + name)
            h = hashlib.sha256()
            with z.open(info) as handle:
                for block in iter(lambda: handle.read(1 << 20), b""):
                    h.update(block)
            if h.hexdigest() != record["sha256"]:
                raise ValueError("Manifest SHA-256 differs: " + name)
        actual = {n for n, m in entries.items() if not m.is_dir()}
        if actual != listed | {manifest_name}:
            raise ValueError("Manifest does not exactly cover archive files")
    return {
        "archive_bytes": ARCHIVE_BYTES, "archive_sha256": ARCHIVE_SHA256,
        "manifest_files_verified": len(records), "archive_members": len(entries),
        "scope": "Saved-table numerical checks and figure reconstruction; no model training",
    }


def extract_checked(archive, destination):
    # validate_archive has already checked identity, every name, and every file.
    destination.mkdir(parents=True, exist_ok=False)
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            target = destination.joinpath(*safe_member(info.filename).parts)
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as source, target.open("xb") as dest:
                    shutil.copyfileobj(source, dest)


def run_script(script, payload, log, env, deadline, extra_args=None):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("Reproduction time budget exhausted")
    command = [sys.executable, "-B", str(script)] + list(extra_args or [])
    with log.open("w", encoding="utf-8") as output:
        subprocess.run(command, cwd=str(payload), env=env, stdout=output,
                       stderr=subprocess.STDOUT, check=True, timeout=remaining)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--mode", choices=("verify", "reproduce"), default="verify")
    parser.add_argument("--out", type=Path,
                        help="New absolute directory outside the repository; required for reproduce")
    parser.add_argument("--max-seconds", type=float, default=300.0)
    args = parser.parse_args()
    repo = args.repo.resolve()
    if not math.isfinite(args.max_seconds) or args.max_seconds <= 0:
        parser.error("--max-seconds must be positive and finite")
    if args.mode == "reproduce" and args.out is None:
        parser.error("--out is required for reproduce")
    out = None
    if args.out is not None:
        if not args.out.is_absolute():
            parser.error("--out must be absolute")
        if args.out.exists() or args.out.is_symlink():
            parser.error("--out must not already exist")
        out = args.out.resolve()
        if out.exists() or inside(out, repo):
            parser.error("--out must be a new directory outside the repository")

    started = time.monotonic()
    report = {
        "passed": False, "repository_version": VERSION, "mode": args.mode,
        "new_models_trained": 0, "installed_software": False,
        "external_downloads": 0,
    }
    created_output = False
    try:
        if (repo / "VERSION").read_text(encoding="utf-8").strip() != VERSION:
            raise ValueError("Expected grn_repo_v25 (VERSION 25.0.0)")
        archive = repo / EXPERIMENT / ARCHIVE_NAME
        report["archive_verification"] = validate_archive(repo)
        if out is not None:
            out.mkdir(parents=True, exist_ok=False)
            created_output = True
        if args.mode == "reproduce":
            extract_checked(archive, out / "payload")
            payload = out / "payload" / ARCHIVE_ROOT
            env = os.environ.copy()
            # Keep all generated files, including Matplotlib caches, outside repo.
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env.pop("PYTHONOPTIMIZE", None)  # Preserve assertions in table checks.
            env["MPLBACKEND"] = "Agg"
            env["MPLCONFIGDIR"] = str(out / "matplotlib_cache")
            env["XDG_CACHE_HOME"] = str(out / "runtime_cache")
            deadline = started + args.max_seconds
            run_script(payload / "code/check_display_tables.py", payload,
                       out / "check_display_tables.log", env, deadline)
            table_result = json.loads((out / "check_display_tables.log").read_text(encoding="utf-8"))
            if table_result.get("status") != "PASS":
                raise ValueError("Display source-table checks did not pass")
            report["table_verification"] = table_result
            # Remove only extracted copies before rebuilding, so stale bundled
            # figures cannot be mistaken for successful newly generated output.
            for stem in FIGURE_NAMES:
                for ext in ("pdf", "svg", "png"):
                    (payload / "figures" / (stem + "." + ext)).unlink()
            run_script(repo / "tools/rebuild_response_figures.py", payload,
                       out / "build_display_figures.log", env, deadline,
                       ["--payload", str(payload)])
            report["renderer_compatibility"] = json.loads(
                (out / "build_display_figures.log").read_text(encoding="utf-8"))
            if report["renderer_compatibility"].get("passed") is not True:
                raise ValueError("Figure metadata compatibility check failed")
            figures = []
            for stem in FIGURE_NAMES:
                for ext in ("pdf", "svg", "png"):
                    path = payload / "figures" / (stem + "." + ext)
                    if not path.is_file() or path.stat().st_size == 0:
                        raise ValueError("Missing rebuilt figure: " + path.name)
                    figures.append({"path": str(path.relative_to(out)),
                                    "bytes": path.stat().st_size,
                                    "sha256": digest_file(path)})
            report["rebuilt_figures"] = figures
            report["figure_identity_note"] = (
                "Numerical tables are checked against archived values. Figure bytes "
                "can vary with Matplotlib versions, fonts and PDF/SVG metadata. "
                "The extracted MANIFEST describes original bundled files before rebuilding."
            )
        report["passed"] = True
    except Exception as exc:
        report["error_type"] = type(exc).__name__
        report["error"] = str(exc)
    finally:
        report["elapsed_seconds"] = time.monotonic() - started
        rendered = json.dumps(report, indent=2) + "\n"
        if created_output:
            (out / REPORT_NAME).write_text(rendered, encoding="utf-8")
        print(rendered, end="")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
