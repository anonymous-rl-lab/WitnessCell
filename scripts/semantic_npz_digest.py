#!/usr/bin/env python3
"""Compute or verify metadata-independent semantic digests for NPZ arrays."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np

DOMAIN = b"WITNESSCELL-NPZ-SEMANTIC-v1\0"


def semantic_digest(path: Path) -> str:
    digest = hashlib.sha256(DOMAIN)
    with np.load(path, allow_pickle=False) as archive:
        for name in sorted(archive.files):
            array = np.asarray(archive[name])
            if array.dtype == object:
                raise ValueError(f"object array is forbidden: {path}:{name}")
            contiguous = np.ascontiguousarray(array)
            fields = (
                name.encode("utf-8"),
                contiguous.dtype.str.encode("ascii"),
                ",".join(str(value) for value in contiguous.shape).encode("ascii"),
                contiguous.tobytes(order="C"),
            )
            for field in fields:
                digest.update(len(field).to_bytes(8, "big"))
                digest.update(field)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        value = line.strip()
        if not value or value.startswith("#"):
            continue
        parts = value.split(maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64:
            raise ValueError(f"invalid manifest line {path}:{line_number}")
        rows.append((parts[0], parts[1]))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="skip manifest entries absent from a deliberately partial reproduction",
    )
    arguments = parser.parse_args()

    failures = 0
    if arguments.manifest:
        for expected, relative in read_manifest(arguments.manifest):
            path = arguments.root / relative
            if not path.is_file():
                if arguments.allow_missing:
                    print(f"SKIP {relative}")
                else:
                    print(f"MISSING {relative}")
                    failures += 1
                continue
            observed = semantic_digest(path)
            if observed != expected:
                print(f"MISMATCH {relative} expected={expected} observed={observed}")
                failures += 1
            else:
                print(f"OK {relative}")

    for path in arguments.paths:
        print(f"{semantic_digest(path)}  {path}")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
