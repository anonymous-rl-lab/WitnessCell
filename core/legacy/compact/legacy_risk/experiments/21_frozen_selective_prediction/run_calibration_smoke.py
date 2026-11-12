#!/usr/bin/env python3
"""Engineering smoke using calibration rows only; formal outcomes stay sealed."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from selective_core import (
    cluster_bootstrap_ratio,
    evaluate_threshold,
    fit_calibration,
    stable_bucket,
)


ROOT = Path(__file__).resolve().parent
GRID = [0.3, 0.5, 0.7, 0.9]
PRIMARY = 0.5


def main() -> None:
    frame = pd.read_csv(ROOT / "assets/calibration_rows.csv")
    frame["smoke_role"] = frame.pair.map(
        lambda pair: "fit" if stable_bucket(str(pair), "gate21-smoke-v1") < 7 else "check"
    )
    fit = frame.loc[frame.smoke_role == "fit"].copy()
    check = frame.loc[frame.smoke_role == "check"].copy()
    assert set(fit.pair).isdisjoint(set(check.pair))
    assert fit.pair.nunique() >= 45 and check.pair.nunique() >= 15

    mapping, thresholds = fit_calibration(
        fit, "estimated_witness_risk", "estimated_realized_mse", GRID
    )
    rows = []
    for row in thresholds.itertuples(index=False):
        result = evaluate_threshold(
            check,
            "estimated_witness_risk",
            "estimated_realized_mse",
            float(row.risk_threshold),
        )
        result.update({
            "target_coverage": float(row.target_coverage),
            "risk_threshold": float(row.risk_threshold),
        })
        rows.append(result)
    results = pd.DataFrame(rows)
    primary = results.loc[np.isclose(results.target_coverage, PRIMARY)].iloc[0].to_dict()
    threshold = float(primary["risk_threshold"])
    bootstrap = cluster_bootstrap_ratio(
        check, "estimated_witness_risk", "estimated_realized_mse",
        threshold, replicates=2000, seed=20260811,
    )
    smoke = {
        "status": "PASS_CALIBRATION_ONLY_ENGINEERING_SMOKE",
        "formal_outcomes_read": False,
        "fit_pairs": int(fit.pair.nunique()),
        "check_pairs": int(check.pair.nunique()),
        "mapping_knots": len(mapping["x_thresholds"]),
        "primary_target_coverage": PRIMARY,
        "primary_result": primary,
        "primary_bootstrap_ratio_ci95": [
            float(np.quantile(bootstrap, 0.025)),
            float(np.quantile(bootstrap, 0.975)),
        ],
        "scientific_gate_used_for_selection": False,
        "note": "This smoke checks direction, nesting and code on calibration identities only; it does not tune the frozen formal score, grid or pass criteria.",
    }
    out = ROOT / "results/calibration_smoke"
    out.mkdir(parents=True, exist_ok=True)
    thresholds.to_csv(out / "fit_thresholds.csv", index=False)
    results.to_csv(out / "heldout_results.csv", index=False)
    (out / "mapping.json").write_text(json.dumps(mapping, indent=2) + "\n")
    (out / "verdict.json").write_text(json.dumps(smoke, indent=2) + "\n")
    print(json.dumps(smoke, indent=2))


if __name__ == "__main__":
    main()
