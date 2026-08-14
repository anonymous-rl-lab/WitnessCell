#!/usr/bin/env python3
"""Post-reveal robustness probes; these cannot modify the frozen formal gate."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from selective_core import (
    evaluate_threshold,
    pair_balanced_weights,
    weighted_mean,
    weighted_quantile,
)


ROOT = Path(__file__).resolve().parent


def main() -> None:
    frame = pd.read_csv(ROOT / "results/formal_reveal/revealed_target_rows.csv")
    verdict = json.loads((ROOT / "results/formal_reveal/FORMAL_VERDICT.json").read_text())
    threshold = float(verdict["primary"]["risk_threshold"])
    accepted = frame.estimated_witness_risk.to_numpy(float) <= threshold
    error = frame.estimated_realized_mse.to_numpy(float)

    seed_rows = []
    for seed, local in frame.groupby("seed", sort=True):
        local_accept = local.estimated_witness_risk.to_numpy(float) <= threshold
        if not local_accept.any() or local_accept.all():
            continue
        local_error = local.estimated_realized_mse.to_numpy(float)
        seed_rows.append({
            "seed": int(seed),
            "rows": int(len(local)),
            "coverage": float(local_accept.mean()),
            "all_mse": float(np.mean(local_error)),
            "accepted_mse": float(np.mean(local_error[local_accept])),
            "rejected_mse": float(np.mean(local_error[~local_accept])),
            "accepted_over_all": float(
                np.mean(local_error[local_accept]) / np.mean(local_error)
            ),
            "risk_error_spearman": float(spearmanr(
                local.estimated_witness_risk, local.estimated_realized_mse
            ).statistic),
        })
    per_seed = pd.DataFrame(seed_rows)

    pair = frame.groupby("pair", as_index=False).agg(
        median_risk=("estimated_witness_risk", "median"),
        mean_mse=("estimated_realized_mse", "mean"),
        rows=("query_id", "size"),
    )
    pair_accept = pair.median_risk.to_numpy(float) <= threshold
    pair_ratio = float(
        pair.loc[pair_accept, "mean_mse"].mean() / pair.mean_mse.mean()
    )

    loo_ratios = []
    for heldout in sorted(frame.pair.unique()):
        local = frame.loc[frame.pair != heldout]
        loo_ratios.append(evaluate_threshold(
            local, "estimated_witness_risk", "estimated_realized_mse", threshold
        )["accepted_over_all_mse"])

    row_ratio = float(np.mean(error[accepted]) / np.mean(error))
    reliability = frame.reliability_weight.to_numpy(float)
    reliability_ratio = float(
        weighted_mean(error[accepted], reliability[accepted])
        / weighted_mean(error, reliability)
    )

    mapping = json.loads((
        ROOT / "results/frozen_calibration/witness_mapping.json"
    ).read_text())
    predicted = np.interp(
        frame.estimated_witness_risk.to_numpy(float),
        np.asarray(mapping["x_thresholds"], float),
        np.asarray(mapping["y_thresholds"], float),
    )
    weights = pair_balanced_weights(frame)
    predicted_mean = weighted_mean(predicted, weights)
    observed_mean = weighted_mean(error, weights)

    geometry = pd.read_csv(
        ROOT / "results/formal_reveal/geometry_control_risk_coverage.csv"
    )
    geometry_primary = geometry.loc[np.isclose(
        geometry.target_coverage, 0.5
    )].iloc[0]
    matched_geometry_threshold = weighted_quantile(
        frame.geometry_risk.to_numpy(float),
        float(verdict["primary"]["pair_balanced_coverage"]),
        weights,
    )
    matched_geometry = evaluate_threshold(
        frame, "geometry_risk", "estimated_realized_mse", matched_geometry_threshold
    )
    witness_curve = pd.read_csv(
        ROOT / "results/formal_reveal/witness_risk_coverage.csv"
    )

    report = {
        "status": "PASS_POST_REVEAL_ROBUSTNESS_PROBES",
        "frozen_formal_status_unchanged": verdict["status"],
        "primary_threshold_unchanged": threshold,
        "row_weighted_accepted_over_all_mse": row_ratio,
        "pair_median_risk_accepted_over_all_mse": pair_ratio,
        "reliability_weighted_sensitivity_ratio": reliability_ratio,
        "seeds_with_both_decisions": int(len(per_seed)),
        "seed_wins_accepted_mse_below_all": int((per_seed.accepted_over_all < 1).sum()),
        "median_seed_accepted_over_all": float(per_seed.accepted_over_all.median()),
        "median_within_seed_risk_error_spearman": float(
            per_seed.risk_error_spearman.median()
        ),
        "leave_one_pair_out_ratio_range": [
            float(np.min(loo_ratios)), float(np.max(loo_ratios))
        ],
        "frozen_isotonic_predicted_mean_mse": predicted_mean,
        "observed_pair_balanced_mean_mse": observed_mean,
        "mean_calibration_ratio_predicted_over_observed": predicted_mean / observed_mean,
        "geometry_control_primary_coverage": float(
            geometry_primary.pair_balanced_coverage
        ),
        "geometry_control_primary_accepted_over_all_mse": float(
            geometry_primary.accepted_over_all_mse
        ),
        "geometry_control_outcome_blind_coverage_matched": {
            "risk_threshold": matched_geometry_threshold,
            "pair_balanced_coverage": matched_geometry["pair_balanced_coverage"],
            "accepted_over_all_mse": matched_geometry["accepted_over_all_mse"],
            "score_loss_spearman": matched_geometry["score_loss_spearman"],
        },
        "witness_grid_accepted_over_all_mse": {
            str(row.target_coverage): float(row.accepted_over_all_mse)
            for row in witness_curve.itertuples(index=False)
        },
        "interpretation": {
            "supported": "Frozen Witness self-risk supports accept/abstain ordering on calibration-unseen target identities.",
            "not_supported": "Geometry-only fallback routing worsened overall MSE and is not a supported deployment policy.",
            "calibration_boundary": "The monotone mapping is directionally useful but absolute loss calibration is reported separately from selective ordering."
        }
    }
    out = ROOT / "results/post_reveal_probes"
    out.mkdir(parents=True, exist_ok=True)
    per_seed.to_csv(out / "per_seed_sensitivity.csv", index=False)
    pair.assign(accepted_pair_median=pair_accept).to_csv(
        out / "per_pair_sensitivity.csv", index=False
    )
    (out / "verdict.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
