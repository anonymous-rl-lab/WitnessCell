#!/usr/bin/env python3
"""Verify immutable upstream source locks and write a machine-readable audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True).strip()


def audit_repo(spec: dict, repo: Path) -> dict:
    expected_commit = spec["commit"]
    expected_tree = spec["tree"]
    actual_commit = git_value(repo, "rev-parse", "HEAD")
    actual_tree = git_value(repo, "rev-parse", "HEAD^{tree}")
    files = []
    for entry in spec["files"]:
        path = repo / entry["path"]
        actual = sha256(path) if path.is_file() else None
        files.append(
            {
                "path": entry["path"],
                "expected_sha256": entry["sha256"],
                "actual_sha256": actual,
                "pass": actual == entry["sha256"],
            }
        )
    return {
        "repository_path": str(repo.resolve()),
        "expected_commit": expected_commit,
        "actual_commit": actual_commit,
        "expected_tree": expected_tree,
        "actual_tree": actual_tree,
        "commit_pass": actual_commit == expected_commit,
        "tree_pass": actual_tree == expected_tree,
        "files": files,
        "pass": actual_commit == expected_commit
        and actual_tree == expected_tree
        and all(item["pass"] for item in files),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument("--normative-repo", type=Path, required=True)
    parser.add_argument("--provenance-repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    source_lock = json.loads((args.experiment_root / "SOURCE_LOCK.json").read_text())
    report = {
        "protocol_id": source_lock["protocol_id"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "normative": audit_repo(source_lock["normative_source"], args.normative_repo),
        "provenance_only": audit_repo(
            source_lock["provenance_only_source"], args.provenance_repo
        ),
    }
    report["pass"] = report["normative"]["pass"] and report["provenance_only"]["pass"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print("SOURCE_AUDIT_PASS" if report["pass"] else "SOURCE_AUDIT_FAIL")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

