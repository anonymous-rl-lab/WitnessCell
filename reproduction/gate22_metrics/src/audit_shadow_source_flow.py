#!/usr/bin/env python3
"""Static taint-boundary audit for the training-only shadow endpoint gate."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    tree = ast.parse(args.source.read_text())
    forbidden_split_reads = []
    train_reads = 0
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "split":
            if node.attr == "train":
                train_reads += 1
            elif node.attr in {"validation", "test", "test_subgroup"}:
                forbidden_split_reads.append({"attribute": node.attr, "line": node.lineno})
    text = args.source.read_text()
    checks = {
        "official_training_split_explicitly_read": train_reads >= 1,
        "validation_or_test_split_never_read": not forbidden_split_reads,
        "formal_h5ad_not_read": "cell_level" not in text and ".h5ad" not in text,
        "weights_built_inside_outer_loo": "scanpy_overestim_scores_from_moments" in text
        and "source_weight_transform" in text,
        "counterfactual_does_not_write_model_archives": "deploy_predictions" not in text
        and "np.savez" not in text,
    }
    report = {
        "status": "PASS_SHADOW_TRAINING_ONLY_SOURCE_FLOW" if all(checks.values()) else "FAIL_SHADOW_SOURCE_FLOW",
        "source": str(args.source),
        "source_sha256": sha256(args.source),
        "train_attribute_reads": train_reads,
        "forbidden_split_reads": forbidden_split_reads,
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
