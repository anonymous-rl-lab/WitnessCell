#!/usr/bin/env python3
"""Replay the frozen confirmation offline and verify named arrays exactly.

This delivery wrapper was written after scoring. It never changes the frozen
scientific implementation, cohort, features, hyperparameter grids or verdict.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_record(value):
    value = np.ascontiguousarray(value)
    return {"dtype": value.dtype.str, "shape": list(value.shape),
            "sha256": hashlib.sha256(value.tobytes(order="C")).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, help="New output directory; otherwise a temporary directory")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    expected = json.loads((root / "EXPECTED_SCIENTIFIC_OUTPUTS.json").read_text())
    out = args.out.resolve() if args.out else Path(tempfile.mkdtemp(prefix="witnesscell_confirmation_replay_"))
    if (out / "EXPRESSION_ACCESS_START.json").exists():
        raise SystemExit("Refusing to overwrite an existing run: " + str(out))
    subprocess.run([sys.executable, str(root / "model/run_confirmation.py"), "run", "--out", str(out)],
                   check=True, cwd=root)
    checks = {}
    array_count = 0
    for name, values in expected["named_arrays"].items():
        with np.load(out / name, allow_pickle=False) as saved:
            checks[name + ":keys"] = sorted(saved.files) == sorted(values)
            for key, record in values.items():
                checks[name + ":" + key] = array_record(saved[key]) == record
                array_count += 1
    for name, digest in expected["csv_file_sha256"].items():
        checks[name] = sha(out / name) == digest
    for name in expected["exact_json_files"]:
        checks[name] = json.loads((out / name).read_text()) == json.loads((root / "run1" / name).read_text())
    completion = json.loads((out / "COMPLETE.json").read_text())
    for key, value in expected["completion_scientific_fields"].items():
        checks["COMPLETE.json:" + key] = completion[key] == value
    report = {"status": "PASS" if all(checks.values()) else "FAIL",
              "output_directory": str(out), "named_arrays_checked": array_count,
              "checks": checks, "decision": completion["decision"],
              "verification": "Exact dtype/shape/byte hashes for every named array, exact CSV bytes and scientific JSON fields.",
              "limits": "Exact binary agreement is expected on the recorded NumPy/BLAS environment. Other numerical environments can differ at floating-point precision; a mismatch is reported, not silently tolerated.",
              "historical_record": "The original run1 UTC and prediction-freeze records remain unchanged. This new run is an offline reproduction after outcomes are known, not a second confirmation experiment."}
    (out / "REPRODUCTION_VERIFICATION.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "decision", "named_arrays_checked", "output_directory")}, indent=2))
    if report["status"] != "PASS":
        raise SystemExit("Reproduction mismatch: inspect REPRODUCTION_VERIFICATION.json")


if __name__ == "__main__":
    main()
