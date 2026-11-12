#!/usr/bin/env python3
"""Cross-audit and seal the four candidate-blind Phase M contracts."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


DATASETS = ("Norman", "Wessels", "Schmidt", "Replogle_exp6")
REQUIRED_KEYS = {
    "dataset",
    "source_sha256",
    "conditions",
    "genes",
    "full_counts",
    "full_variances",
    "first_counts",
    "second_counts",
    "full_means",
    "first_means",
    "second_means",
    "weights",
    "weight_evaluable",
    "second_pvalues_adjusted",
    "control",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--contracts", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    rows = []
    all_checks = []
    for dataset in DATASETS:
        contract_path = args.contracts / f"{dataset}.phase_m_contract.npz"
        ledger_path = args.contracts / f"{dataset}.split_ledger.json"
        report_path = args.contracts / f"{dataset}.phase_m_report.json"
        ledger = json.loads(ledger_path.read_text())
        upstream_report = json.loads(report_path.read_text())
        with np.load(contract_path, allow_pickle=False) as contract:
            keys = set(contract.files)
            conditions = np.asarray(contract["conditions"], dtype=str)
            genes = np.asarray(contract["genes"], dtype=str)
            weights = np.asarray(contract["weights"], dtype=np.float64)
            checks = {
                "required_schema_exact": keys == REQUIRED_KEYS,
                "no_object_dtype": all(contract[key].dtype.kind != "O" for key in contract.files),
                "all_numeric_finite": all(
                    bool(np.all(np.isfinite(contract[key])))
                    for key in contract.files
                    if contract[key].dtype.kind in "fci"
                ),
                "condition_identity_exact": conditions.tolist() == ledger["condition_order"],
                "genes_unique": len(set(genes.tolist())) == genes.size,
                "weights_shape": weights.shape == (conditions.size, genes.size),
                "weight_rows_nonnegative": bool(np.all(weights >= 0)),
                "counts_partition": bool(
                    np.array_equal(
                        contract["first_counts"] + contract["second_counts"],
                        contract["full_counts"],
                    )
                ),
                "upstream_report_pass": upstream_report["status"]
                == "PASS_PREFREEZE_CANDIDATE_BLIND_CONTRACT",
                "upstream_contract_hash_exact": upstream_report["contract_sha256"]
                == sha256(contract_path),
                "upstream_ledger_hash_exact": upstream_report["split_ledger_sha256"]
                == sha256(ledger_path),
            }
            for seed in (1, 2, 3):
                prediction_path = (
                    repo
                    / "experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo"
                    / dataset
                    / f"seed{seed}/deploy_predictions.npz"
                )
                with np.load(prediction_path, allow_pickle=False) as prediction:
                    split = ledger["split_contract"][str(seed)]
                    checks[f"seed{seed}_test_order_exact"] = (
                        prediction["conditions"].astype(str).tolist() == split["test"]
                    )
                    checks[f"seed{seed}_gene_order_exact"] = np.array_equal(
                        prediction["genes"].astype(str), genes
                    )
                    condition_index = {value: index for index, value in enumerate(conditions)}
                    expected_counts = np.asarray(
                        [contract["full_counts"][condition_index[value]] for value in split["test"]]
                    )
                    checks[f"seed{seed}_test_counts_exact"] = np.array_equal(
                        prediction["test_cell_counts"], expected_counts
                    )
            all_checks.extend(checks.values())
            rows.append(
                {
                    "dataset": dataset,
                    "contract_sha256": sha256(contract_path),
                    "ledger_sha256": sha256(ledger_path),
                    "report_sha256": sha256(report_path),
                    "conditions": int(conditions.size),
                    "genes": int(genes.size),
                    "evaluable_weight_conditions": int(contract["weight_evaluable"].sum()),
                    "checks": checks,
                    "pass": all(checks.values()),
                }
            )

    payload = {
        "status": "PASS_SEALED_PHASE_M_CONTRACTS" if all(all_checks) else "FAIL_PHASE_M_CONTRACTS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_predictions_used_to_construct_splits_or_weights": False,
        "datasets": rows,
        "contracts": {
            row["dataset"]: {
                "contract_sha256": row["contract_sha256"],
                "split_ledger_sha256": row["ledger_sha256"],
                "report_sha256": row["report_sha256"],
            }
            for row in rows
        },
        "pass": all(all_checks),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
