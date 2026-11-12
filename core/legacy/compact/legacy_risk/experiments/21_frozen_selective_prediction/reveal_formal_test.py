#!/usr/bin/env python3
"""Once-only formal reveal under the frozen selective-prediction protocol."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from selective_core import (
    cluster_bootstrap_ratio,
    evaluate_threshold,
    sha256,
    within_seed_permutation,
)


ROOT = Path(__file__).resolve().parent
REPLICATES = 20000
PRIMARY = 0.5


def verify_pre_reveal() -> None:
    expected = (ROOT / "FROZEN_PROTOCOL.sha256").read_text().split()[0]
    assert sha256(ROOT / "FROZEN_PROTOCOL.json") == expected
    for line in (ROOT / "PRE_REVEAL_MANIFEST.sha256").read_text().splitlines():
        digest, relative = line.split(maxsplit=1)
        assert sha256(ROOT / relative) == digest, relative
    protocol = json.loads((ROOT / "FROZEN_PROTOCOL.json").read_text())
    assert protocol["freeze_status"] == "FROZEN_BEFORE_FORMAL_OUTCOME_REVEAL"
    assert protocol["formal_outcomes_read_at_freeze"] is False


def main() -> None:
    verify_pre_reveal()
    outcome_path = ROOT / "sealed/formal_outcomes.csv"
    asset_audit = json.loads((ROOT / "assets/ASSET_AUDIT.json").read_text())
    assert sha256(outcome_path) == asset_audit["sealed_outcomes_sha256"]
    os.chmod(outcome_path, 0o400)
    risk = pd.read_csv(ROOT / "assets/formal_risk_only.csv")
    outcomes = pd.read_csv(outcome_path)
    assert risk.query_id.is_unique and outcomes.query_id.is_unique
    frame = risk.merge(outcomes, on="query_id", validate="one_to_one")
    assert len(frame) == asset_audit["formal_eligible_rows"]
    assert frame.pair.nunique() == asset_audit["formal_eligible_pairs"]

    witness_thresholds = pd.read_csv(
        ROOT / "results/frozen_calibration/witness_thresholds.csv"
    )
    geometry_thresholds = pd.read_csv(
        ROOT / "results/frozen_calibration/geometry_thresholds.csv"
    )
    rows = []
    for row in witness_thresholds.itertuples(index=False):
        result = evaluate_threshold(
            frame, "estimated_witness_risk", "estimated_realized_mse",
            float(row.risk_threshold),
        )
        result.update({
            "score": "estimated_witness_risk",
            "target_coverage": float(row.target_coverage),
            "risk_threshold": float(row.risk_threshold),
        })
        rows.append(result)
    witness_results = pd.DataFrame(rows)

    geometry_rows = []
    for row in geometry_thresholds.itertuples(index=False):
        result = evaluate_threshold(
            frame, "geometry_risk", "estimated_realized_mse",
            float(row.risk_threshold),
        )
        result.update({
            "score": "geometry_risk_control",
            "target_coverage": float(row.target_coverage),
            "risk_threshold": float(row.risk_threshold),
        })
        geometry_rows.append(result)
    geometry_results = pd.DataFrame(geometry_rows)

    primary = witness_results.loc[
        np.isclose(witness_results.target_coverage, PRIMARY)
    ].iloc[0].to_dict()
    threshold = float(primary["risk_threshold"])
    bootstrap = cluster_bootstrap_ratio(
        frame, "estimated_witness_risk", "estimated_realized_mse", threshold,
        replicates=REPLICATES, seed=20260811,
    )
    permuted = within_seed_permutation(
        frame, "estimated_witness_risk", "estimated_realized_mse", threshold,
        replicates=REPLICATES, seed=20260812,
    )
    observed = float(primary["accepted_mse"])
    permutation_p = float((1 + np.sum(permuted <= observed)) / (len(permuted) + 1))
    curve = pd.concat([
        witness_results[["pair_balanced_coverage", "accepted_mse"]],
        pd.DataFrame([{
            "pair_balanced_coverage": 1.0,
            "accepted_mse": float(primary["all_mse"]),
        }]),
    ], ignore_index=True).sort_values("pair_balanced_coverage")
    curve_spearman = float(spearmanr(
        curve.pair_balanced_coverage, curve.accepted_mse
    ).statistic)

    accepted = frame.estimated_witness_risk.to_numpy(float) <= threshold
    fallback_error = np.where(
        accepted,
        frame.estimated_realized_mse.to_numpy(float),
        frame.geometry_realized_mse.to_numpy(float),
    )
    baseline_error = frame.estimated_realized_mse.to_numpy(float)
    fallback_ratio = float(np.mean(fallback_error) / np.mean(baseline_error))

    criteria = {
        "coverage_non_degenerate": 0.35 <= primary["pair_balanced_coverage"] <= 0.65,
        "practical_effect": primary["accepted_over_all_mse"] <= 0.80,
        "cluster_uncertainty": float(np.quantile(bootstrap, 0.975)) < 1.0,
        "random_selection_control": permutation_p < 0.05,
        "rejected_separation": primary["accepted_mse"] < primary["rejected_mse"],
    }
    result = {
        "status": (
            "PASS_FROZEN_SELECTIVE_PREDICTION_GATE"
            if all(criteria.values()) else "FAIL_FROZEN_SELECTIVE_PREDICTION_GATE"
        ),
        "protocol_sha256": sha256(ROOT / "FROZEN_PROTOCOL.json"),
        "test_rows": int(len(frame)),
        "test_pairs": int(frame.pair.nunique()),
        "pair_identity_overlap_with_calibration": 0,
        "primary_target_coverage": PRIMARY,
        "primary": primary,
        "bootstrap_accepted_over_all_ci95": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "within_seed_random_selection_p_one_sided": permutation_p,
        "risk_coverage_spearman": curve_spearman,
        "secondary_risk_coverage_support": curve_spearman >= 0.8,
        "secondary_fallback_over_always_witness_mse": fallback_ratio,
        "criteria": criteria,
        "scope": "Norman development-panel pseudobulk; calibration-unseen target identities; retrospective frozen reanalysis",
    }
    out = ROOT / "results/formal_reveal"
    out.mkdir(parents=True, exist_ok=True)
    witness_results.to_csv(out / "witness_risk_coverage.csv", index=False)
    geometry_results.to_csv(out / "geometry_control_risk_coverage.csv", index=False)
    curve.to_csv(out / "capability_curve.csv", index=False)
    frame.assign(accepted_primary=accepted).to_csv(out / "revealed_target_rows.csv", index=False)
    np.savez_compressed(
        out / "inference_arrays.npz",
        bootstrap_accepted_over_all_ratio=bootstrap,
        random_selection_accepted_mse=permuted,
    )
    (out / "FORMAL_VERDICT.json").write_text(json.dumps(result, indent=2) + "\n")
    os.chmod(outcome_path, 0o400)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
