#!/usr/bin/env python3
"""Audit comparator source, local availability and prespecified scope."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--source-repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    lock_path = args.repo / "experiments/22_metric_calibration_stress/COMPARATOR_SOURCE_LOCK.json"
    scope_path = args.repo / "experiments/22_metric_calibration_stress/COMPARATOR_SCOPE.json"
    lock = json.loads(lock_path.read_text())
    scope = json.loads(scope_path.read_text())
    commit = subprocess.check_output(
        ["git", "-C", str(args.source_repo), "rev-parse", "HEAD"], text=True
    ).strip()
    tree = subprocess.check_output(
        ["git", "-C", str(args.source_repo), "rev-parse", "HEAD^{tree}"], text=True
    ).strip()
    file_checks = {}
    for relative, expected in lock["files"].items():
        actual = sha256(args.source_repo / relative)
        file_checks[relative] = {"expected": expected, "actual": actual, "pass": actual == expected}
    raw_predictions = list(
        (args.repo / "experiments/15_scperturbench_sota/module").rglob("pred.tsv")
    ) + list((args.repo / "experiments/15_scperturbench_sota/module").rglob("result.h5ad"))
    mean_source = (args.source_repo / "Perturbation_generalization/Genetic/linearModel_mean.py").read_text()
    checks = {
        "commit_exact": commit == lock["commit"],
        "tree_exact": tree == lock["tree"],
        "source_files_exact": all(value["pass"] for value in file_checks.values()),
        "trainmean_combo_four_branches_present": all(
            token in mean_source
            for token in (
                "geneA in adata_train_single_pert and geneB not in adata_train_single_pert",
                "geneA not in adata_train_single_pert and geneB in adata_train_single_pert",
                "geneA not in adata_train_single_pert and geneB not in adata_train_single_pert",
                "geneA in adata_train_single_pert and geneB in adata_train_single_pert",
            )
        ),
        "basecontrol_source_present": "senario == 'controlMean'" in mean_source,
        "comparator_scope_fixed": scope["fixed_before_prediction_scoring"],
        "linear_not_silently_imputed": scope["methods"]["linearModel"].startswith("NOT_ADJUDICATED"),
        "published_pool_not_silently_imputed": scope["published_method_pool"].startswith("NOT_ADJUDICATED"),
    }
    report = {
        "status": "PASS_COMPARATOR_SOURCE_AND_AVAILABILITY_AUDIT" if all(checks.values()) else "FAIL_COMPARATOR_AUDIT",
        "commit": commit,
        "tree": tree,
        "files": file_checks,
        "raw_published_prediction_assets_found": len(raw_predictions),
        "rscript_available": shutil.which("Rscript") is not None,
        "linearModel_formal_status": "NOT_ADJUDICATED",
        "trainMean_formal_status": "AVAILABLE_EXACT_CONDITION_MEAN_REGENERATION",
        "baseControl_formal_status": "AVAILABLE_EXACT_CONDITION_MEAN_REGENERATION",
        "published_metric_parity_boundary": "released generateExp samples unseeded Gaussian cells; published rounded distribution-derived rows cannot be a bitwise oracle for the deterministic underlying condition mean",
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
