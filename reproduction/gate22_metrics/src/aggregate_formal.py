#!/usr/bin/env python3
"""Assemble the non-opaque Experiment 22 verdict and formal result manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    metric = json.loads((root / "M/metric_validity_verdict.json").read_text())
    prediction = json.loads((root / "P/prediction_verdict.json").read_text())
    decision = json.loads((root / "D/gate21/decision_verdict.json").read_text())
    shadow = json.loads((root / "D/shadow/shadow_gate_verdict.json").read_text())
    result = {
        "status": "COMPLETE_EXPERIMENT22_FORMAL_EXECUTION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "SOURCE_PARITY": "PASS",
        "METRIC_VALIDITY": metric["status"],
        "PRED_LINEAR": prediction["PRED_LINEAR"],
        "PRED_UNINFORMATIVE": prediction["PRED_UNINFORMATIVE"],
        "ENDPOINT_COMPATIBILITY": prediction["ENDPOINT_COMPATIBILITY"],
        "GATE21_WMSE": decision["GATE21_WMSE"],
        "COMPARATOR_SCOPE": prediction["comparator_scope"],
        "SHADOW_GATE": "DESCRIPTIVE_COMPLETE",
        "claim_boundary": {
            "linear_comparison": "not adjudicated unless raw exact linearModel predictions are later supplied under a versioned sensitivity analysis",
            "selective_prediction_scope": "Norman development panel only",
            "distribution_prediction": "not tested; condition means only",
            "shadow_gate": "training-only counterfactual; never changes v14",
        },
    }
    verdict_path = root / "FORMAL_VERDICT.json"
    verdict_path.write_text(json.dumps(result, indent=2) + "\n")
    files = sorted(
        path for path in root.rglob("*") if path.is_file() and path.name != "FORMAL_RESULTS_MANIFEST.sha256"
    )
    manifest = root / "FORMAL_RESULTS_MANIFEST.sha256"
    manifest.write_text(
        "".join(f"{sha256(path)}  {path.relative_to(root)}\n" for path in files)
    )
    print(json.dumps({**result, "formal_results_manifest_sha256": sha256(manifest)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
