#!/usr/bin/env python3
"""Aggregate the frozen 3x3 six-way gate without treating fits as biology."""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


COMPARISONS = (
    ("GEARS+Witness", "GEARS"),
    ("CPA+Witness", "CPA"),
    ("WitnessCell", "GEARS"),
    ("WitnessCell", "CPA"),
)
EXPECTED_STRATEGIES = {
    "Saturation", "GEARS", "CPA", "GEARS+Witness", "CPA+Witness", "WitnessCell"
}
EPS = 1e-12


def exact_signflip_p(differences: np.ndarray) -> float:
    observed = float(np.mean(differences))
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(differences))))
    statistics = (signs * differences[None, :]).mean(axis=1)
    return float(np.mean(statistics >= observed - 1e-15))


def comparison(
    targets: pd.DataFrame,
    candidate: str,
    comparator: str,
    seed: int,
    gate: dict,
) -> dict:
    pair = targets.groupby(["pair", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean"),
        full_effect_cosine=("full_effect_cosine", "mean"),
    )
    mse = pair.pivot(index="pair", columns="strategy", values="residual_mse")
    cosine = pair.pivot(index="pair", columns="strategy", values="full_effect_cosine")
    base = mse[comparator].to_numpy(float)
    trial = mse[candidate].to_numpy(float)
    differences = base - trial
    relative = float(differences.mean() / max(base.mean(), EPS))
    rng = np.random.default_rng(seed)
    draw = rng.integers(0, len(base), size=(20000, len(base)))
    sampled = (base[draw].mean(axis=1) - trial[draw].mean(axis=1)) / np.maximum(
        base[draw].mean(axis=1), EPS
    )

    run = targets.groupby(["data_seed", "model_seed", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean")
    ).pivot(index=["data_seed", "model_seed"], columns="strategy", values="residual_mse")
    by_data = targets.groupby(["data_seed", "strategy"], as_index=False).agg(
        residual_mse=("residual_mse", "mean")
    ).pivot(index="data_seed", columns="strategy", values="residual_mse")
    data_improvement = (by_data[comparator] - by_data[candidate]) / by_data[comparator]
    cosine_delta = float((cosine[candidate] - cosine[comparator]).mean())
    conditions = {
        "relative_mse": relative >= float(gate["mean_relative_residual_mse_improvement_min"]),
        "pair_bootstrap_ci": float(np.quantile(sampled, 0.025)) > float(gate["pair_bootstrap_ci95_low_gt"]),
        "exact_signflip": exact_signflip_p(differences) < float(gate["exact_pair_signflip_p_lt"]),
        "all_data_seeds_positive": bool((data_improvement > 0).all()),
        "run_win_rate": float(np.mean(run[candidate] < run[comparator])) >= float(gate["run_win_rate_min"]),
        "full_effect_cosine_safety": cosine_delta >= float(gate["full_effect_cosine_delta_min"]),
    }
    return {
        "candidate": candidate,
        "comparator": comparator,
        "pairs": int(len(base)),
        "relative_residual_mse_improvement": relative,
        "pair_bootstrap_ci95": [float(np.quantile(sampled, .025)), float(np.quantile(sampled, .975))],
        "exact_pair_signflip_p_one_sided": exact_signflip_p(differences),
        "pair_win_rate": float(np.mean(trial < base)),
        "run_win_rate": float(np.mean(run[candidate] < run[comparator])),
        "data_seed_relative_improvements": {str(k): float(v) for k, v in data_improvement.items()},
        "full_effect_cosine_delta": cosine_delta,
        "conditions": conditions,
        "pass": bool(all(conditions.values())),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=Path(__file__).resolve().parent / "FROZEN_SCIENTIFIC_PROTOCOL.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    rows: list[pd.DataFrame] = []
    seen: set[tuple[int, int]] = set()
    for path in sorted(args.root.glob("d*_m*/per_target.csv")):
        local = pd.read_csv(path)
        key = (int(local.data_seed.iloc[0]), int(local.model_seed.iloc[0]))
        if key in seen:
            raise ValueError(f"duplicate run {key}")
        seen.add(key)
        if set(local.strategy) != EXPECTED_STRATEGIES or local.pair.nunique() != 14 or len(local) != 84:
            raise ValueError(f"invalid six-way target table: {path}")
        rows.append(local)
    expected = {
        (d, m) for d in protocol["factorial_design"]["data_seeds"]
        for m in protocol["factorial_design"]["model_seeds"]
    }
    if seen != expected:
        raise ValueError(f"formal matrix incomplete; got {sorted(seen)}, expected {sorted(expected)}")
    targets = pd.concat(rows, ignore_index=True)
    results = [
        comparison(targets, candidate, comparator, 20260820 + index,
                   protocol["success_rule_per_comparison"])
        for index, (candidate, comparator) in enumerate(COMPARISONS)
    ]
    lookup = {(row["candidate"], row["comparator"]): row for row in results}
    augmentation_gears = lookup[("GEARS+Witness", "GEARS")]["pass"]
    augmentation_cpa = lookup[("CPA+Witness", "CPA")]["pass"]
    independent = (
        lookup[("WitnessCell", "GEARS")]["pass"]
        and lookup[("WitnessCell", "CPA")]["pass"]
    )
    if augmentation_gears and augmentation_cpa and independent:
        status = "PASS_UNIVERSAL_AUGMENTATION_AND_INDEPENDENT_AIVC"
    elif independent:
        status = "PASS_INDEPENDENT_AIVC_HETEROGENEOUS_AUGMENTATION"
    elif augmentation_gears or augmentation_cpa:
        status = "PARTIAL_MODEL_DEPENDENT_AUGMENTATION"
    else:
        status = "STOP_EXTERNAL_NAMED_AIVC_SUPERIORITY"
    verdict = {
        "status": status,
        "runs": 9,
        "named_gpu_fits": 18,
        "target_pairs": 14,
        "comparisons": results,
        "gears_augmentation_pass": augmentation_gears,
        "cpa_augmentation_pass": augmentation_cpa,
        "independent_witnesscell_vs_both_named_pass": independent,
        "claim_boundary": "All targets share FDPS. Pair statistics describe this pathway panel; data/model seeds and cells are not biological replicates."
    }
    targets.to_csv(args.out / "all_target_rows.csv", index=False)
    pd.DataFrame(results).drop(columns=["conditions"]).to_json(
        args.out / "comparison_table.json", orient="records", indent=2
    )
    (args.out / "formal_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()

