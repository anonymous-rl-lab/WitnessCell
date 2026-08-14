#!/usr/bin/env python3
"""Fail when release files contain common double-blind identity leaks."""

from __future__ import annotations

import argparse
import re
import tarfile
from collections.abc import Iterable
from pathlib import Path, PurePosixPath
from zipfile import BadZipFile, ZipFile

MAX_TEXT_BYTES = 8 * 1024 * 1024
SKIP_PARTS = {".git", ".mypy_cache", ".pytest_cache", ".ruff_cache", ".venv", "__pycache__"}
TEXT_SUFFIXES = {
    "", ".cff", ".cfg", ".csv", ".ini", ".json", ".md", ".py", ".rst",
    ".toml", ".txt", ".yaml", ".yml",
}
PATTERNS = {
    "email address": re.compile(r"(?<![\w.-])[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}(?![\w.-])"),
    "ORCID": re.compile(r"\b0000-000[0-3]-\d{4}-\d{3}[\dX]\b", re.IGNORECASE),
    "local Unix path": re.compile(r"/(?:home|Users|workspace|root)/[^\s\"']+"),
    "local Windows path": re.compile(r"\b[A-Z]:\\(?:Users|Documents and Settings)\\[^\s\"']+", re.IGNORECASE),
    "preregistration account URL": re.compile(r"https?://(?:www\.)?osf\.io/[a-z0-9]{5}/?", re.IGNORECASE),
    "grant identifier": re.compile(r"\b(?:grant|award)\s*(?:no\.?|number|#)?\s*[A-Z]{1,6}[- ]?\d{4,}\b", re.IGNORECASE),
}


def _scan_text(label: str, data: bytes) -> list[str]:
    if len(data) > MAX_TEXT_BYTES:
        return []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return []
    failures = []
    for description, pattern in PATTERNS.items():
        match = pattern.search(text)
        if match:
            line = text.count("\n", 0, match.start()) + 1
            failures.append(f"{label}:{line}: {description}")
    return failures


def _safe_member(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts


def _scan_zip(path: Path) -> list[str]:
    failures: list[str] = []
    try:
        with ZipFile(path) as archive:
            for info in archive.infolist():
                if not _safe_member(info.filename):
                    failures.append(f"{path}!{info.filename}: unsafe archive path")
                    continue
                if Path(info.filename).suffix.lower() in TEXT_SUFFIXES and info.file_size <= MAX_TEXT_BYTES:
                    failures.extend(_scan_text(f"{path}!{info.filename}", archive.read(info)))
    except BadZipFile:
        failures.append(f"{path}: invalid ZIP-compatible archive")
    return failures


def _scan_tar(path: Path) -> list[str]:
    failures: list[str] = []
    with tarfile.open(path, "r:*") as archive:
        for member in archive.getmembers():
            if not _safe_member(member.name):
                failures.append(f"{path}!{member.name}: unsafe archive path")
                continue
            if not member.isfile() or member.size > MAX_TEXT_BYTES:
                continue
            if Path(member.name).suffix.lower() not in TEXT_SUFFIXES:
                continue
            stream = archive.extractfile(member)
            if stream is not None:
                failures.extend(_scan_text(f"{path}!{member.name}", stream.read()))
    return failures


def _files(targets: Iterable[Path]) -> Iterable[Path]:
    for target in targets:
        if target.is_file():
            yield target
            continue
        for path in target.rglob("*"):
            if path.is_file() and not any(part in SKIP_PARTS for part in path.parts):
                yield path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="+", type=Path)
    arguments = parser.parse_args()
    failures: list[str] = []
    scanned = 0
    for path in _files(arguments.targets):
        scanned += 1
        lower = path.name.lower()
        if lower.endswith((".whl", ".zip", ".wcell")):
            failures.extend(_scan_zip(path))
        elif lower.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tar.xz")):
            failures.extend(_scan_tar(path))
        elif path.suffix.lower() in TEXT_SUFFIXES and path.stat().st_size <= MAX_TEXT_BYTES:
            failures.extend(_scan_text(str(path), path.read_bytes()))
    if failures:
        print("ANONYMITY AUDIT FAILED")
        for failure in failures:
            print(f"- {failure}")
        return 1
    print(f"ANONYMITY AUDIT PASSED: {scanned} files checked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
