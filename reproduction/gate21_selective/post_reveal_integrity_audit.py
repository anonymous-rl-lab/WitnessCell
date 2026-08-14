#!/usr/bin/env python3
"""Verify frozen hashes and all saved post-reveal inferential quantities."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from selective_core import sha256


ROOT = Path(__file__).resolve().parent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compact-release",
        action="store_true",
        help="allow the path-only ASSET_AUDIT metadata to be anonymized",
    )
    arguments = parser.parse_args()
    skipped: list[str] = []
    for line in (ROOT / "PRE_REVEAL_MANIFEST.sha256").read_text().splitlines():
        expected, relative = line.split(maxsplit=1)
        if arguments.compact_release and relative == "assets/ASSET_AUDIT.json":
            skipped.append(relative)
            continue
        assert sha256(ROOT / relative) == expected, relative
    protocol_hash = (ROOT / "FROZEN_PROTOCOL.sha256").read_text().split()[0]
    assert sha256(ROOT / "FROZEN_PROTOCOL.json") == protocol_hash

    verdict = json.loads((ROOT / "results/formal_reveal/FORMAL_VERDICT.json").read_text())
    with np.load(ROOT / "results/formal_reveal/inference_arrays.npz", allow_pickle=False) as z:
        bootstrap = z["bootstrap_accepted_over_all_ratio"].astype(float)
        random_mse = z["random_selection_accepted_mse"].astype(float)
    assert len(bootstrap) == 20000 and len(random_mse) == 20000
    expected_ci = [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))]
    assert np.allclose(expected_ci, verdict["bootstrap_accepted_over_all_ci95"], atol=0, rtol=0)
    observed = float(verdict["primary"]["accepted_mse"])
    expected_p = float((1 + np.sum(random_mse <= observed)) / (len(random_mse) + 1))
    assert expected_p == verdict["within_seed_random_selection_p_one_sided"]
    criteria = verdict["criteria"]
    assert criteria == {
        "coverage_non_degenerate": True,
        "practical_effect": True,
        "cluster_uncertainty": True,
        "random_selection_control": True,
        "rejected_separation": True,
    }
    assert verdict["status"] == "PASS_FROZEN_SELECTIVE_PREDICTION_GATE"
    probe = json.loads((ROOT / "results/post_reveal_probes/verdict.json").read_text())
    assert probe["frozen_formal_status_unchanged"] == verdict["status"]
    status = (
        "PASS_GATE21_COMPACT_POST_REVEAL_INTEGRITY_AUDIT"
        if arguments.compact_release
        else "PASS_GATE21_POST_REVEAL_INTEGRITY_AUDIT"
    )
    result = {
        "status": status,
        "scientific_verdict": verdict["status"],
        "protocol_sha256": protocol_hash,
        "pre_reveal_scientific_files_unchanged": True,
        "sanitized_metadata_files": skipped,
        "bootstrap_replicates": len(bootstrap),
        "permutation_replicates": len(random_mse),
        "saved_inference_recomputed_exactly": True,
    }
    if not arguments.compact_release:
        out = ROOT / "results/final_audit"
        out.mkdir(parents=True, exist_ok=True)
        (out / "FINAL_AUDIT.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
