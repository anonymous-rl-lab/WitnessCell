#!/usr/bin/env python3
"""Record the deterministic CPU environment used by Experiment 22."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_info


PACKAGES = (
    "numpy",
    "scipy",
    "pandas",
    "scikit-learn",
    "h5py",
    "anndata",
    "scanpy",
    "tqdm",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    checks = {
        "python_at_least_3_12": sys.version_info >= (3, 12),
        "float64_ieee": np.finfo(np.float64).bits == 64,
    }
    report = {
        "status": "PASS_ENVIRONMENT_AUDIT" if all(checks.values()) else "FAIL_ENVIRONMENT_AUDIT",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": {name: version(name) for name in PACKAGES},
        "threadpools": threadpool_info(),
        "gpu_required": False,
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
