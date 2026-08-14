#!/usr/bin/env python3
"""scPerturBench entry point for the external WitnessCell adapter.

The benchmark repository uses one script per method.  WitnessCell keeps its
audited implementation in a standalone package and this entry point delegates
without installing or importing GEARS/CPA/PyTorch.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adapter", type=Path, default=os.environ.get("WITNESSCELL_ADAPTER_ROOT"))
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--seed", type=int, required=True, choices=(1, 2, 3))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--gears-assets",
        type=Path,
        default=os.environ.get("WITNESSCELL_GEARS_ASSETS"),
    )
    args = parser.parse_args()
    if args.adapter is None:
        raise SystemExit("pass --adapter or set WITNESSCELL_ADAPTER_ROOT")
    if args.gears_assets is None:
        args.gears_assets = args.adapter / "data" / "gears_assets"
    runner = args.adapter / "run_witnesscell.py"
    exporter = args.adapter / "export_result_h5ad.py"
    prediction_dir = args.out / "prediction"
    result = args.out / "result.h5ad"
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(args.adapter) + (
        os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else ""
    )
    subprocess.run([
        sys.executable, str(runner), "--data", str(args.data),
        "--dataset", args.dataset, "--seed", str(args.seed),
        "--gears-assets", str(args.gears_assets),
        "--out", str(prediction_dir),
    ], check=True, env=environment)
    subprocess.run([
        sys.executable, str(exporter), "--data", str(args.data),
        "--prediction", str(prediction_dir / "deploy_predictions.npz"),
        "--out", str(result),
    ], check=True, env=environment)
    print(f"WITNESSCELL_SCPERTURBENCH_ENTRY_PASS {result}")


if __name__ == "__main__":
    main()
