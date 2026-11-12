#!/usr/bin/env python3
"""Verify the compact experiment payload's manifest without any large inputs."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "SUPPLEMENT_MANIFEST.json").read_text())
    checked = []
    for row in manifest["files"]:
        path = root / row["path"]
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != row["size_bytes"]:
            raise ValueError("File size mismatch: " + row["path"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise ValueError("File hash mismatch: " + row["path"])
        checked.append(row["path"])
    files = [p for p in root.rglob("*") if p.is_file() and "__pycache__" not in p.parts]
    actual = {p.relative_to(root).as_posix() for p in files}
    if actual != set(checked) | {"SUPPLEMENT_MANIFEST.json"}:
        raise ValueError("Unexpected or missing payload files")
    size = sum(p.stat().st_size for p in files)
    if size >= 5_000_000:
        raise ValueError("Uncompressed supplement exceeds its 5,000,000-byte limit")
    print(json.dumps({"status": "PASS", "verified_files": len(checked),
                      "payload_files": len(files), "uncompressed_bytes": size,
                      "limit_bytes": 5_000_000,
                      "large_canonical_assets_verified": False}))


if __name__ == "__main__":
    main()
