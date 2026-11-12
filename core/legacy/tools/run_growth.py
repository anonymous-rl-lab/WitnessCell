#!/usr/bin/env python3
"""Verify saved Horlbeck results or reproduce the frozen two-stage experiment.

Uses the current Python interpreter and its installed NumPy/SciPy. Installs
nothing, downloads nothing, and never writes into the release. Verification
runs the original verify.py in a temporary extraction. Reproduction copies the
unchanged payload to a fresh external directory and invokes its reproduce.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import zipfile


ARCHIVE_RELATIVE = "experiments/26_horlbeck_two_stage/WitnessCell_Horlbeck_TwoStage_Repro_20261018.zip"
ARCHIVE_SHA256 = "b60e01230f3f43f156271697fe1fdc52b6c4ce5f68b522c3664c6e4aa6ff4082"
PAYLOAD_NAME = "WitnessCell_Horlbeck_TwoStage_Repro_20261018"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def safe_relative(name: str) -> PurePosixPath:
    p = PurePosixPath(name)
    if (not name or "\\" in name or "\x00" in name or ":" in name
            or p.is_absolute() or ".." in p.parts or str(p) != name):
        raise ValueError(f"Unsafe archive path: {name!r}")
    return p


def validate_archive(archive: Path) -> dict:
    """Validate archive identity, safe members, exact inventory and all hashes."""
    if not archive.is_file() or archive.is_symlink():
        raise ValueError("Expected an ordinary frozen Horlbeck ZIP file")
    actual_sha = sha(archive)
    if actual_sha != ARCHIVE_SHA256:
        raise ValueError(f"Frozen Horlbeck archive SHA-256 mismatch: {actual_sha}")
    with zipfile.ZipFile(archive) as z:
        members = z.infolist()
        names = [m.filename for m in members]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate archive entries")
        files = {}
        for member in members:
            path = safe_relative(member.filename.rstrip("/"))
            if path.parts[0] != PAYLOAD_NAME:
                raise ValueError("Unexpected archive root")
            file_type = stat.S_IFMT(member.external_attr >> 16)
            if file_type not in (0, stat.S_IFREG, stat.S_IFDIR):
                raise ValueError("Archive links and special files are forbidden")
            if not member.is_dir():
                files[PurePosixPath(*path.parts[1:]).as_posix()] = member
        manifest = json.loads(z.read(PAYLOAD_NAME + "/PACKAGE_MANIFEST.json"))
        entries = manifest["files"]
        if set(files) != set(entries) | {"PACKAGE_MANIFEST.json"}:
            raise ValueError("Internal manifest and archive file inventory disagree")
        for name, expected in entries.items():
            safe_relative(name)
            data = z.read(files[name])
            if len(data) != expected["bytes"] or hashlib.sha256(data).hexdigest() != expected["sha256"]:
                raise ValueError(f"Internal manifest mismatch: {name}")
    return {"archive_sha256": actual_sha, "archive_bytes": archive.stat().st_size,
            "internal_manifest_files_verified": len(entries),
            "archive_members": len(members)}


def extract_checked(archive: Path, destination: Path) -> Path:
    # validate_archive has already rejected special, duplicate and unsafe paths.
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            relative = safe_relative(member.filename.rstrip("/"))
            target = destination / relative
            if not target.resolve().is_relative_to(destination.resolve()):
                raise ValueError("Archive extraction would escape destination")
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(member) as source, target.open("xb") as output:
                    shutil.copyfileobj(source, output)
    return destination / PAYLOAD_NAME


def run_bounded(command: list[str], cwd: Path, log: Path, seconds: float) -> None:
    if seconds <= 0:
        raise TimeoutError("Horlbeck wall-clock limit reached")
    env = os.environ.copy()
    for name in ("PYTHONDONTWRITEBYTECODE", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                 "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[name] = "1"
    env["MPLBACKEND"] = "Agg"
    env["PYTHONOPTIMIZE"] = "0"
    env.pop("PYTHONPATH", None)
    with log.open("w", encoding="utf-8") as stream:
        process = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=os.name == "posix")
        try:
            code = process.wait(timeout=seconds)
        except BaseException:
            # The original reproduce.py starts children; stop the complete group.
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            raise
    if code:
        raise RuntimeError(f"Frozen entry point exited {code}:\n{log.read_text()[-5000:]}")


def main() -> None:
    # Let an enclosing quick-check timeout stop our child's separate group too.
    def interrupted(signum, frame):
        raise InterruptedError(f"Interrupted by signal {signum}")
    if os.name == "posix":
        signal.signal(signal.SIGTERM, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--mode", choices=("verify", "reproduce"), default="verify")
    parser.add_argument("--stage", choices=("all", "K562", "Jurkat"), default="all",
                        help="Stage to refit; verify always checks the complete two-stage payload.")
    parser.add_argument("--out", type=Path,
                        help="Fresh external directory for receipt/logs; required for reproduce.")
    parser.add_argument("--max-seconds", type=float, default=300.0)
    args = parser.parse_args()
    repo = args.repo.resolve()
    if (repo / "VERSION").read_text().strip() != "25.0.0":
        parser.error("Expected grn_repo_v25 VERSION 25.0.0")
    if not 0 < args.max_seconds < float("inf"):
        parser.error("--max-seconds must be finite and positive")
    if args.mode == "reproduce" and args.out is None:
        parser.error("--out is required for reproduction")
    out = args.out.resolve() if args.out else None
    if out is not None:
        if out.is_relative_to(repo):
            parser.error("--out must be outside the release")
        if out.exists():
            parser.error("--out must be a new directory")
    archive = repo / ARCHIVE_RELATIVE
    start = time.monotonic()
    integrity = validate_archive(archive)
    if out is not None:
        out.mkdir(parents=True)
    report = {"status": "RUNNING", "repository_version": "25.0.0", "mode": args.mode,
              "requested_stage": args.stage, "max_seconds": args.max_seconds,
              "python": sys.executable, **integrity,
              "scope": "saved evidence verification" if args.mode == "verify" else
                       "numerical reproduction of existing results; not a new confirmation",
              "software_installed": False, "gpu_requested": False, "external_downloads": 0}
    error = None
    try:
        with tempfile.TemporaryDirectory(prefix="witnesscell-v25-horlbeck-") as temporary:
            work = Path(temporary)
            destination = out / "payload" if args.mode == "reproduce" else work / "payload"
            payload = extract_checked(archive, destination)
            log = (out if out is not None else work) / "horlbeck.log"
            if args.mode == "verify":
                receipt = (out if out is not None else work) / "saved_results_check.json"
                command = [sys.executable, "-B", str(payload / "verify.py"), "--out", str(receipt)]
            else:
                receipt = out / "reproduced" / "REPRODUCTION_CHECK.json"
                command = [sys.executable, "-B", str(payload / "reproduce.py"),
                           "--stage", args.stage, "--out", str(out / "reproduced")]
            run_bounded(command, payload, log, args.max_seconds - (time.monotonic() - start))
            result = json.loads(receipt.read_text())
            if result.get("status") != "PASS":
                raise ValueError("Original Horlbeck verifier/reproducer did not return PASS")
            report["original_result"] = result
            report["status"] = "PASS"
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = f"{type(exc).__name__}: {exc}"
        error = exc
    finally:
        report["archive_sha256_after"] = sha(archive)
        report["archive_unchanged"] = report["archive_sha256_after"] == ARCHIVE_SHA256
        if not report["archive_unchanged"]:
            report["status"] = "FAIL"
        report["elapsed_seconds"] = time.monotonic() - start
        if out is not None:
            (out / "horlbeck_receipt.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    if error is not None or report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
