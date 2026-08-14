#!/usr/bin/env python3
"""Aggregate and seal candidate-blind Phase M outputs without binary DRF gates."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


DATASETS = ("Norman", "Wessels", "Schmidt", "Replogle_exp6")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-m-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.phase_m_dir
    summary_rows = []
    nir_rows = []
    eligibility_rows = []
    for dataset in DATASETS:
        validity = pd.read_csv(root / f"{dataset}.metric_validity.csv")
        nir_controls = pd.read_csv(root / f"{dataset}.nir_controls.csv")
        for row in validity.itertuples(index=False):
            eligibility_rows.append(
                {
                    "dataset": dataset,
                    "seed": int(row.seed),
                    "condition": str(row.condition),
                    "evaluable": bool(row.evaluable),
                }
            )
        for column in [value for value in validity.columns if value.startswith("drf__")]:
            _, negative, positive, metric = column.split("__", maxsplit=3)
            values = validity.loc[validity["evaluable"].astype(bool), column].to_numpy(float)
            finite = values[np.isfinite(values)]
            summary_rows.append(
                {
                    "dataset": dataset,
                    "metric": metric,
                    "negative_control": negative,
                    "positive_control": positive,
                    "valid_n": int(finite.size),
                    "mean_drf": float(np.mean(finite)) if finite.size else float("nan"),
                    "median_drf": float(np.median(finite)) if finite.size else float("nan"),
                    "fraction_drf_le_zero": float(np.mean(finite <= 0)) if finite.size else float("nan"),
                    "binary_validity_threshold_applied": False,
                }
            )
        for (seed, control), local in nir_controls.groupby(["seed", "control"]):
            nir_rows.append(
                {
                    "dataset": dataset,
                    "seed": int(seed),
                    "control": control,
                    "conditions": len(local),
                    "mean_nir": float(local["nir"].mean()),
                    "median_nir": float(local["nir"].median()),
                    "drf_certified": False,
                }
            )
    drf_summary = pd.DataFrame(summary_rows)
    nir_summary = pd.DataFrame(nir_rows)
    eligibility = pd.DataFrame(eligibility_rows)
    drf_path = root / "drf_summary.csv"
    nir_path = root / "nir_control_summary.csv"
    eligibility_path = root / "eligible_operations.csv"
    drf_summary.to_csv(drf_path, index=False)
    nir_summary.to_csv(nir_path, index=False)
    eligibility.to_csv(eligibility_path, index=False)

    artifact_paths = sorted(
        path
        for path in root.iterdir()
        if path.is_file() and path.name not in {"PHASE_M_MANIFEST.sha256", "metric_validity_verdict.json"}
    )
    manifest_path = root / "PHASE_M_MANIFEST.sha256"
    manifest_path.write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in artifact_paths)
    )
    result = {
        "status": "PASS_PHASE_M_EXECUTION_AND_SEAL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_predictions_read": False,
        "binary_drf_cutoff_used": False,
        "fixed_metrics": ["wmse", "weighted_r2_deltapert", "nir"],
        "formal_rows": int(len(eligibility)),
        "formal_evaluable_rows": int(eligibility["evaluable"].sum()),
        "phase_m_manifest_sha256": sha256(manifest_path),
        "interpretation": "continuous candidate-independent control/DRF characterization; NIR reported but not DRF-certified at the locked commit",
    }
    verdict_path = root / "metric_validity_verdict.json"
    verdict_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
