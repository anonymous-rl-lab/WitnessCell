#!/usr/bin/env python3
"""Audit GitHub readiness, release-manifest integrity and repository hygiene."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAX_FILE_BYTES = 50 * 1024 * 1024
FORBIDDEN_PARTS = {".git", ".mypy_cache", ".pytest_cache", ".ruff_cache", ".venv", "__pycache__", "build", "dist"}
FORBIDDEN_SUFFIXES = {".pyc", ".pyo"}
FORBIDDEN_NAMES = {"grn_repo_v18.zip", "WitnessCell_preregistered_replication_20260810T105615Z.tar.gz"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def files() -> list[Path]:
    return sorted(path for path in ROOT.rglob("*") if path.is_file())


def verify_manifest() -> list[str]:
    manifest = ROOT / "RELEASE_MANIFEST.sha256"
    if not manifest.exists():
        return ["RELEASE_MANIFEST.sha256 is missing"]
    errors: list[str] = []
    for line_number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        parts = value.split(maxsplit=1)
        if len(parts) != 2:
            errors.append(f"invalid manifest line {line_number}")
            continue
        expected, relative = parts
        path = ROOT / relative
        if not path.is_file():
            errors.append(f"manifest file missing: {relative}")
        elif sha256(path) != expected:
            errors.append(f"manifest mismatch: {relative}")
    return errors


def main() -> int:
    errors: list[str] = []
    scanned = files()
    for path in scanned:
        relative = path.relative_to(ROOT)
        if any(part in FORBIDDEN_PARTS for part in relative.parts):
            errors.append(f"forbidden generated path: {relative}")
        if any(part.endswith(".egg-info") for part in relative.parts):
            errors.append(f"forbidden package metadata path: {relative}")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES or path.name in FORBIDDEN_NAMES:
            errors.append(f"forbidden generated/archive file: {relative}")
        if path.stat().st_size > MAX_FILE_BYTES:
            errors.append(f"file exceeds 50 MiB Git gate: {relative}")
    errors.extend(verify_manifest())

    anonymity = subprocess.run(
        [sys.executable, str(ROOT / "scripts/audit_anonymity.py"), str(ROOT)],
        check=False,
    )
    if anonymity.returncode:
        errors.append("anonymity audit failed")

    if errors:
        print("REPOSITORY AUDIT FAILED")
        for error in errors:
            print(f"- {error}")
        return 1
    largest = max(scanned, key=lambda path: path.stat().st_size)
    print(
        "REPOSITORY AUDIT PASSED: "
        f"{len(scanned)} files; largest={largest.relative_to(ROOT)} "
        f"({largest.stat().st_size} bytes)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
