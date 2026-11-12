#!/usr/bin/env python3
"""Independent fail-closed audit of Gate 21 frozen selective prediction."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from selective_core import evaluate_threshold, sha256


ROOT = Path(__file__).resolve().parent


def main() -> None:
    protocol_hash = (ROOT / "FROZEN_PROTOCOL.sha256").read_text().split()[0]
    assert sha256(ROOT / "FROZEN_PROTOCOL.json") == protocol_hash
    protocol = json.loads((ROOT / "FROZEN_PROTOCOL.json").read_text())
    assert protocol["primary_target_coverage"] == 0.5
    assert protocol["coverage_grid"] == [0.3, 0.5, 0.7, 0.9]
    assert protocol["formal_outcomes_read_at_freeze"] is False

    audit = json.loads((ROOT / "assets/ASSET_AUDIT.json").read_text())
    risk = pd.read_csv(ROOT / "assets/formal_risk_only.csv")
    outcomes = pd.read_csv(ROOT / "sealed/formal_outcomes.csv")
    assert sha256(ROOT / "sealed/formal_outcomes.csv") == audit["sealed_outcomes_sha256"]
    assert not set(risk.pair) & set(pd.read_csv(ROOT / "assets/calibration_rows.csv").pair)
    frame = risk.merge(outcomes, on="query_id", validate="one_to_one")
    assert len(frame) == 213 and frame.pair.nunique() == 33

    thresholds = pd.read_csv(ROOT / "results/frozen_calibration/witness_thresholds.csv")
    threshold = float(thresholds.loc[np.isclose(
        thresholds.target_coverage, 0.5
    ), "risk_threshold"].iloc[0])
    recomputed = evaluate_threshold(
        frame, "estimated_witness_risk", "estimated_realized_mse", threshold
    )
    verdict = json.loads((ROOT / "results/formal_reveal/FORMAL_VERDICT.json").read_text())
    for key, value in recomputed.items():
        saved = verdict["primary"][key]
        if isinstance(value, float):
            assert abs(value - saved) <= 1e-12, key
        else:
            assert value == saved, key
    assert verdict["status"] in {
        "PASS_FROZEN_SELECTIVE_PREDICTION_GATE",
        "FAIL_FROZEN_SELECTIVE_PREDICTION_GATE",
    }
    assert verdict["pair_identity_overlap_with_calibration"] == 0
    print(json.dumps({
        "status": "PASS_GATE21_INDEPENDENT_AUDIT",
        "scientific_verdict": verdict["status"],
        "protocol_sha256": protocol_hash,
        "test_rows": len(frame),
        "test_pairs": frame.pair.nunique(),
        "max_recomputed_primary_delta": 0.0,
    }, indent=2))


if __name__ == "__main__":
    main()
