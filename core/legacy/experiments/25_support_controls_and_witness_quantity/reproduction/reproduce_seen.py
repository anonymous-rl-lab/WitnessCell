#!/usr/bin/env python3
"""Reproduce the delivered seen statistics using only retained package inputs.

Requires Python 3.9+ and NumPy. The frozen analysis script is copied verbatim;
its two commands receive an explicit, temporary --source path. No original
sandbox path, raw cell matrix, mean-predictor fit, or new risk method is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


SOURCE_FILENAMES = (
    "training_oof_records.csv",
    "risk_models.json",
    "frozen_query_outputs.csv",
    "evaluation_metadata.json",
    "split_manifest.json",
    "frozen_selections.csv",
    "scored_queries.csv",
)
RESULT_KEY = ("dataset", "scope", "coverage", "loss", "comparison")
EFFECT_KEY = ("dataset", "scope", "comparison")
ATOL = 1e-12
RTOL = 1e-12


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def compare_values(expected, actual, path, audit):
    """Compare every statistical field; counts and text must match exactly."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(expected) != set(actual):
            raise ValueError("Object fields differ at " + path)
        for key in expected:
            compare_values(expected[key], actual[key], path + "." + key, audit)
        return
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise ValueError("Array length differs at " + path)
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare_values(left, right, f"{path}[{index}]", audit)
        return
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError("Non-numeric value differs at " + path)
        return
    if not isinstance(actual, (int, float)) or isinstance(actual, bool):
        raise ValueError("Expected a numeric value at " + path)
    if isinstance(expected, int):
        if not isinstance(actual, int) or actual != expected:
            raise ValueError("Integer count differs at " + path)
    else:
        if not (math.isfinite(expected) and math.isfinite(actual)):
            raise ValueError("Nonfinite reported statistic at " + path)
        error = abs(expected - actual)
        audit["maximum_absolute_numeric_difference"] = max(
            audit["maximum_absolute_numeric_difference"], error
        )
        if not math.isclose(expected, actual, rel_tol=RTOL, abs_tol=ATOL):
            raise ValueError(
                f"Numeric mismatch at {path}: expected {expected!r}, got {actual!r}"
            )
    audit["numeric_values_compared"] += 1


def indexed(rows, fields, name):
    result = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        if key in result:
            raise ValueError(f"Duplicate {name} identity: {key!r}")
        result[key] = row
    return result


def compare_statistics(expected, actual):
    audit = {
        "status": "PASS",
        "comparison_rule": "all statistical fields; exact text/counts and finite floats within atol/rtol",
        "absolute_tolerance": ATOL,
        "relative_tolerance": RTOL,
        "numeric_values_compared": 0,
        "maximum_absolute_numeric_difference": 0.0,
        "excluded_top_level_fields": ["completed_utc", "elapsed_seconds", "freeze_sha256"],
    }
    for key in ("protocol_id", "design", "bootstrap"):
        compare_values(expected[key], actual[key], key, audit)
    for name, fields, required_count in (
        ("results", RESULT_KEY, 480),
        ("effect_quartile_control", EFFECT_KEY, 2),
    ):
        left = indexed(expected[name], fields, name)
        right = indexed(actual[name], fields, name)
        if len(left) != required_count or set(left) != set(right):
            raise ValueError(f"Expected the same {required_count} {name} identities")
        for key in sorted(left):
            compare_values(left[key], right[key], name + repr(key), audit)
        audit[name + "_groups_compared"] = len(left)
    audit["confidence_intervals_compared"] = sum(
        sum(key.endswith("_ci95") for key in row)
        for name in ("results", "effect_quartile_control")
        for row in expected[name]
    )
    return audit


def reproduce(bundle, output):
    seen = bundle / "seen"
    plan_path = bundle / "analysis_plan.json"
    script_path = seen / "run_seen.py"
    expected_path = seen / "evaluation" / "results.json"
    plan = read_json(plan_path)
    freeze = read_json(seen / "all_new_decisions_frozen.json")
    if digest(script_path) != freeze["script_sha256"]:
        raise ValueError("The delivered frozen analysis script has changed")
    if digest(plan_path) != freeze["plan_sha256"]:
        raise ValueError("The delivered frozen analysis plan has changed")
    if not expected_path.is_file():
        raise FileNotFoundError(expected_path)
    cfg = plan["seen_comparison"]
    for dataset in cfg["datasets"]:
        for seed in cfg["seeds"]:
            source_batch = seen / "results" / f"{dataset}_seed{seed}"
            for name in SOURCE_FILENAMES:
                path = source_batch / ("original_" + name)
                if not path.is_file():
                    raise FileNotFoundError(path)
    output.mkdir(parents=True, exist_ok=False)
    copied_seen = output / "seen"
    copied_seen.mkdir()
    shutil.copyfile(plan_path, output / "analysis_plan.json")
    shutil.copyfile(script_path, copied_seen / "run_seen.py")
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="witnesscell_seen_source_") as temporary:
        source_root = Path(temporary) / "source"
        for dataset in cfg["datasets"]:
            for seed in cfg["seeds"]:
                batch_name = f"{dataset}_seed{seed}"
                restored = source_root / "results" / batch_name
                restored.mkdir(parents=True)
                for name in SOURCE_FILENAMES:
                    shutil.copyfile(
                        seen / "results" / batch_name / ("original_" + name),
                        restored / name,
                    )
        for stage in ("prepare", "evaluate"):
            command = [sys.executable, str(copied_seen / "run_seen.py"), stage,
                       "--source", str(source_root)]
            with (output / (stage + ".log")).open("w", encoding="utf-8") as log:
                subprocess.run(command, cwd=output, stdout=log,
                               stderr=subprocess.STDOUT, check=True)
    reproduced_path = copied_seen / "evaluation" / "results.json"
    audit = compare_statistics(read_json(expected_path), read_json(reproduced_path))
    audit.update({
        "reproduction": "delivered seen analysis, unchanged protocol and risk methods",
        "reference_results_sha256": digest(expected_path),
        "reproduced_results_sha256": digest(reproduced_path),
        "frozen_script_sha256": digest(script_path),
        "frozen_plan_sha256": digest(plan_path),
        "python_version": sys.version,
        "elapsed_seconds": time.monotonic() - start,
        "temporary_source_removed": True,
        "output_directory": str(output),
        "reproduced_results": str(reproduced_path),
    })
    write_json(output / "reproduction_status.json", audit)
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path,
                        default=Path(__file__).resolve().parents[1],
                        help="Directory containing analysis_plan.json and seen/ (default: script parent package)")
    parser.add_argument("--output", required=True, type=Path,
                        help="A new output directory; existing directories are never overwritten")
    args = parser.parse_args()
    bundle, output = args.package_root.resolve(), args.output.resolve()
    if output == bundle or (bundle / "seen") in output.parents:
        parser.error("The output must be outside the frozen seen directory")
    if output.exists():
        parser.error("The output directory already exists; choose a new directory")
    try:
        result = reproduce(bundle, output)
    except Exception as exc:
        result = {"status": "FAIL", "error_type": type(exc).__name__,
                  "error": str(exc), "output_directory": str(output)}
        if output.is_dir():
            write_json(output / "reproduction_status.json", result)
        print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
        return 1
    print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
