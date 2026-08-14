#!/usr/bin/env python3
"""Audit immutable WitnessCell, Gate 19 and Gate 21 assets for Experiment 22."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_ARCHIVE_KEYS = {
    "dataset",
    "seed",
    "conditions",
    "subgroups",
    "genes",
    "prediction",
    "factorized_prediction",
    "baseline_prediction",
    "baseline_factorized_prediction",
    "control",
    "control_variance",
    "control_count",
    "test_cell_counts",
}
FORBIDDEN_TRUTH_TOKENS = ("truth", "target", "observed", "heldout")
FALLBACK_SPLITS = {("Schmidt", 3), ("Wessels", 1), ("Wessels", 2), ("Wessels", 3)}
GATE21_THRESHOLD = 0.0923227147328771


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def identity_hash(values: list[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def archive_row(path: Path) -> tuple[dict, list[str]]:
    with np.load(path, allow_pickle=False) as archive:
        keys = set(archive.files)
        object_free = all(archive[key].dtype.kind != "O" for key in archive.files)
        finite = all(
            bool(np.all(np.isfinite(archive[key])))
            for key in archive.files
            if archive[key].dtype.kind in "fc"
        )
        dataset = str(archive["dataset"].item())
        seed = int(archive["seed"].item())
        conditions = archive["conditions"].astype(str).tolist()
        unit_ids = [f"{dataset}::{seed}::{condition}" for condition in conditions]
        return (
            {
                "path": str(path),
                "sha256": sha256(path),
                "dataset": dataset,
                "seed": seed,
                "rows": len(conditions),
                "genes": int(archive["genes"].size),
                "schema_pass": keys == REQUIRED_ARCHIVE_KEYS,
                "object_free": object_free,
                "finite": finite,
                "truth_key_free": not any(
                    token in key.lower() for key in keys for token in FORBIDDEN_TRUTH_TOKENS
                ),
            },
            unit_ids,
        )


def compare_archives(v14: Path, v13: Path) -> dict:
    with np.load(v14, allow_pickle=False) as left, np.load(v13, allow_pickle=False) as right:
        common = sorted(set(left.files) & set(right.files))
        exact = {key: bool(np.array_equal(left[key], right[key])) for key in common}
        return {
            "v14": str(v14),
            "v13": str(v13),
            "common_keys": common,
            "key_exact": exact,
            "all_common_exact": all(exact.values()),
            "schemas_equal": set(left.files) == set(right.files),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root
    v14_root = root / "experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo"
    v13_root = root / "experiments/18_dual_head_evidence_gate/predictions/formal_combo"

    rows: list[dict] = []
    unit_ids: list[str] = []
    for path in sorted(v14_root.glob("*/seed*/deploy_predictions.npz")):
        row, ids = archive_row(path)
        rows.append(row)
        unit_ids.extend(ids)

    paired = pd.read_csv(
        root
        / "experiments/19_v14_incremental_amplitude_gate/results/21_formal_all12_v14_vs_v13/paired_by_condition.csv"
    )
    paired_ids = [
        f"{row.dataset}::{int(row.seed)}::{row.condition}" for row in paired.itertuples(index=False)
    ]
    split_audit = pd.read_csv(
        root
        / "experiments/19_v14_incremental_amplitude_gate/results/21_formal_all12_v14_vs_v13/split_fallback_audit.csv"
    )
    observed_inactive = {
        (str(row.dataset), int(row.seed))
        for row in split_audit.itertuples(index=False)
        if not bool(row.correction_active)
    }

    fallbacks = []
    for dataset, seed in sorted(FALLBACK_SPLITS):
        fallbacks.append(
            compare_archives(
                v14_root / dataset / f"seed{seed}" / "deploy_predictions.npz",
                v13_root / dataset / f"seed{seed}" / "deploy_predictions.npz",
            )
        )

    gate21_rows = pd.read_csv(
        root
        / "experiments/21_frozen_selective_prediction/results/formal_reveal/revealed_target_rows.csv"
    )
    gate21_verdict = json.loads(
        (
            root
            / "experiments/21_frozen_selective_prediction/results/formal_reveal/FORMAL_VERDICT.json"
        ).read_text()
    )
    recalculated_accept = gate21_rows["estimated_witness_risk"].to_numpy() <= GATE21_THRESHOLD
    accepted_ids = gate21_rows.loc[recalculated_accept, "query_id"].astype(str).tolist()
    stored_accepted_ids = gate21_rows.loc[gate21_rows["accepted_primary"], "query_id"].astype(str).tolist()

    checks = {
        "archive_count_12": len(rows) == 12,
        "archive_schema_object_finite_truth_free": all(
            row["schema_pass"] and row["object_free"] and row["finite"] and row["truth_key_free"]
            for row in rows
        ),
        "unit_count_654": len(unit_ids) == 654 and len(set(unit_ids)) == 654,
        "unit_identity_matches_gate19": set(unit_ids) == set(paired_ids),
        "unit_identity_order_hash_matches": identity_hash(sorted(unit_ids))
        == identity_hash(sorted(paired_ids)),
        "inactive_split_identity_8_4": observed_inactive == FALLBACK_SPLITS,
        "inactive_archives_exact": all(
            item["schemas_equal"] and item["all_common_exact"] for item in fallbacks
        ),
        "gate21_threshold_bitwise": np.float64(
            gate21_verdict["primary"]["risk_threshold"]
        ).tobytes()
        == np.float64(GATE21_THRESHOLD).tobytes(),
        "gate21_query_count_213": len(gate21_rows) == 213,
        "gate21_acceptance_identity_exact": accepted_ids == stored_accepted_ids,
    }
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "archives": rows,
        "formal_unit_identity_sha256": identity_hash(sorted(unit_ids)),
        "inactive_expected": [f"{dataset}::seed{seed}" for dataset, seed in sorted(FALLBACK_SPLITS)],
        "inactive_archive_comparisons": fallbacks,
        "gate21": {
            "threshold_repr": repr(GATE21_THRESHOLD),
            "query_count": len(gate21_rows),
            "accepted_count": len(accepted_ids),
            "accepted_query_identity_sha256": identity_hash(sorted(accepted_ids)),
        },
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print("FROZEN_ASSET_AUDIT_PASS" if report["pass"] else "FROZEN_ASSET_AUDIT_FAIL")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
