#!/usr/bin/env python3
"""Dependency-free structural checks for built wheel and source distribution."""

from __future__ import annotations

import argparse
import base64
import tarfile
from email.parser import BytesParser
from email.policy import default
from hashlib import sha256
from pathlib import Path, PurePosixPath
from zipfile import ZipFile

REQUIRED_WHEEL = {
    "witnesscell/__init__.py",
    "witnesscell/_version.py",
    "witnesscell/cli.py",
    "witnesscell/core.py",
    "witnesscell/data.py",
    "witnesscell/model.py",
    "witnesscell/py.typed",
    "witnesscell/risk.py",
    "witnesscell/selective.py",
    "witnesscell/serialization.py",
}


def _safe(names: list[str]) -> None:
    for name in names:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"unsafe archive path: {name}")


def verify_wheel(path: Path) -> None:
    with ZipFile(path) as archive:
        names = archive.namelist()
        _safe(names)
        if len(names) != len(set(names)):
            raise ValueError("wheel contains duplicate entries")
        missing = REQUIRED_WHEEL - set(names)
        if missing:
            raise ValueError(f"wheel is missing package files: {sorted(missing)}")
        metadata_name = next(name for name in names if name.endswith(".dist-info/METADATA"))
        record_name = next(name for name in names if name.endswith(".dist-info/RECORD"))
        metadata = BytesParser(policy=default).parsebytes(archive.read(metadata_name))
        if metadata["Name"] != "witnesscell" or metadata["Version"] != "0.1.0":
            raise ValueError("wheel name/version metadata mismatch")
        if metadata["License-Expression"] != "Apache-2.0":
            raise ValueError("wheel license expression mismatch")
        if metadata.get("Author") or metadata.get("Author-email") or metadata.get("Home-page"):
            raise ValueError("wheel exposes identity-bearing metadata")
        record = archive.read(record_name).decode("utf-8").splitlines()
        for line in record:
            name, encoded_hash, encoded_size = line.rsplit(",", 2)
            if name == record_name:
                if encoded_hash or encoded_size:
                    raise ValueError("RECORD must leave its own digest empty")
                continue
            payload = archive.read(name)
            algorithm, digest = encoded_hash.split("=", 1)
            if algorithm != "sha256":
                raise ValueError(f"unexpected RECORD digest algorithm: {algorithm}")
            actual = base64.urlsafe_b64encode(sha256(payload).digest()).rstrip(b"=").decode("ascii")
            if digest != actual or int(encoded_size) != len(payload):
                raise ValueError(f"RECORD integrity mismatch: {name}")


def verify_sdist(path: Path) -> None:
    with tarfile.open(path, "r:*") as archive:
        names = archive.getnames()
        _safe(names)
        if len(names) != len(set(names)):
            raise ValueError("sdist contains duplicate entries")
        prefix = "witnesscell-0.1.0/"
        required = {
            prefix + "pyproject.toml",
            prefix + "LICENSE",
            prefix + "README.md",
            prefix + "src/witnesscell/core.py",
            prefix + "tests/test_model.py",
            prefix + "scripts/audit_anonymity.py",
            prefix + "scripts/verify_reference_parity.py",
        }
        missing = required - set(names)
        if missing:
            raise ValueError(f"sdist is missing release files: {sorted(missing)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    arguments = parser.parse_args()
    verify_wheel(arguments.wheel)
    verify_sdist(arguments.sdist)
    print("DISTRIBUTION STRUCTURE AND RECORD INTEGRITY PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
