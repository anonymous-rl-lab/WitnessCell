#!/usr/bin/env python3
"""Training-only signal-weighted shadow of the frozen v14 endpoint gate."""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from metric_core import WeightVector, source_weight_transform, weighted_delta_r2, wmse


DATASETS = ("Norman", "Wessels", "Schmidt", "Replogle_exp6")


def scanpy_overestim_scores_from_moments(
    held: str,
    remaining: list[str],
    means: dict[str, np.ndarray],
    variances: dict[str, np.ndarray],
    counts: dict[str, int],
) -> np.ndarray:
    """Exact summary-moment equivalent of Scanpy target-vs-rest t scores."""

    n_held = counts[held]
    remaining_counts = np.asarray([counts[node] for node in remaining], dtype=np.float64)
    remaining_means = np.stack([means[node] for node in remaining])
    total = float(remaining_counts.sum())
    rest_mean = np.average(remaining_means, axis=0, weights=remaining_counts)
    rest_population_var = np.sum(
        remaining_counts[:, None]
        * (
            np.stack([variances[node] for node in remaining])
            + np.square(remaining_means - rest_mean)
        ),
        axis=0,
    ) / total
    held_unbiased = variances[held] * n_held / max(n_held - 1, 1)
    rest_unbiased = rest_population_var * total / max(total - 1, 1)
    scores, _ = stats.ttest_ind_from_stats(
        mean1=means[held],
        std1=np.sqrt(held_unbiased),
        nobs1=n_held,
        mean2=rest_mean,
        std2=np.sqrt(rest_unbiased),
        # Scanpy's t-test_overestim_var deliberately substitutes target n.
        nobs2=n_held,
        equal_var=False,
    )
    return np.nan_to_num(scores, nan=0.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    v14_src = repo / "experiments/19_v14_incremental_amplitude_gate/src"
    official_src = repo / "experiments/15_scperturbench_sota/module"
    sys.path.insert(0, str(v14_src))
    sys.path.insert(0, str(official_src))
    from identity_head import one_sided_lower  # type: ignore
    from multihead_v14 import (  # type: ignore
        RECIPE,
        fit_correction,
        fit_frozen_v13,
        top100,
    )
    from official_split import (  # type: ignore
        filter_gears_supported_conditions,
        load_gears_supported_genes,
        make_official_split,
    )
    from probe_multhead_training import build_bank, sparse_features  # type: ignore

    gears_assets = official_src / "data/gears_assets"
    supported_genes = load_gears_supported_genes(gears_assets)
    with (gears_assets / "gene2go_all.pkl").open("rb") as handle:
        gene2go = pickle.load(handle)
    rows = []
    split_rows = []
    for dataset in DATASETS:
        cache_path = (
            repo
            / "experiments/18_dual_head_evidence_gate/cache"
            / f"{dataset}.condition_moments.npz"
        )
        data = np.load(cache_path, allow_pickle=False)
        conditions = data["conditions"].astype(str).tolist()
        genes = data["genes"].astype(str).tolist()
        matrix = data["means"].astype(np.float64)
        variance_matrix = data["variances"].astype(np.float64)
        count_array = data["counts"].astype(np.int64)
        means = dict(zip(conditions, matrix, strict=True))
        variances = dict(zip(conditions, variance_matrix, strict=True))
        counts = {condition: int(value) for condition, value in zip(conditions, count_array, strict=True)}
        control = means["control"]
        official_conditions = filter_gears_supported_conditions(conditions, supported_genes)
        gene_index = {gene: index for index, gene in enumerate(genes)}

        for seed in (1, 2, 3):
            split = make_official_split(official_conditions, seed)
            # Non-negotiable taint boundary: only training singles enter every
            # outcome, weight, fit and lower-bound calculation below.
            nodes = sorted(
                condition
                for condition in split.train
                if condition != "control" and "+" not in condition and condition in gene_index
            )
            effects = {node: means[node] - control for node in nodes}
            top_indices = {
                node: top100(
                    node,
                    means,
                    variances,
                    counts,
                    control,
                    variances["control"],
                    counts["control"],
                )
                for node in nodes
            }
            local_rows = []
            for outer in nodes:
                remaining_nodes = [node for node in nodes if node != outer]
                remaining_effects = {node: effects[node] for node in remaining_nodes}
                outer_baseline = fit_frozen_v13(
                    remaining_effects,
                    means,
                    variances,
                    counts,
                    control,
                    variances["control"],
                    counts["control"],
                    genes,
                    gene2go,
                )
                correction_weight, _ = fit_correction(
                    remaining_nodes,
                    remaining_effects,
                    outer_baseline,
                    genes,
                    gene2go,
                    top_indices,
                )
                base = outer_baseline.predict(outer, remaining_effects, gene_index, gene2go)
                outer_bank = build_bank(
                    outer,
                    remaining_effects,
                    genes,
                    gene2go,
                    10,
                    (5, 10, 20, 40),
                    (3, 5, 10),
                )
                correction = sum(
                    weight * feature
                    for weight, feature in zip(
                        correction_weight,
                        sparse_features(outer_bank, RECIPE),
                        strict=True,
                    )
                )
                final = base + correction
                truth = effects[outer]
                scores = scanpy_overestim_scores_from_moments(
                    outer, remaining_nodes, means, variances, counts
                )
                weight = source_weight_transform(scores, genes, genes)
                baseline_direction = np.mean(
                    np.stack([effects[node] for node in remaining_nodes]), axis=0
                )
                base_wmse = wmse(base, truth, weight)
                final_wmse = wmse(final, truth, weight)
                base_r2 = weighted_delta_r2(base, truth, baseline_direction, weight)
                final_r2 = weighted_delta_r2(final, truth, baseline_direction, weight)
                local_rows.append(
                    {
                        "dataset": dataset,
                        "seed": seed,
                        "condition": outer,
                        "base_wmse": base_wmse,
                        "final_wmse": final_wmse,
                        "wmse_gain": (base_wmse - final_wmse) / max(base_wmse, 1e-12),
                        "base_weighted_r2": base_r2,
                        "final_weighted_r2": final_r2,
                        "weighted_r2_delta": final_r2 - base_r2,
                        "weight_sum": float(weight.values.sum()),
                    }
                )
            local = pd.DataFrame(local_rows)
            rows.extend(local_rows)
            wmse_lower = one_sided_lower(local["wmse_gain"].to_numpy(), 0.95)
            r2_lower = one_sided_lower(local["weighted_r2_delta"].to_numpy(), 0.95)
            shadow_active = bool(wmse_lower > 0 and r2_lower > 0)
            manifest = json.loads(
                (
                    repo
                    / "experiments/19_v14_incremental_amplitude_gate/predictions/formal_combo"
                    / dataset
                    / f"seed{seed}/manifest.json"
                ).read_text()
            )
            original_active = bool(manifest["identity_head"]["correction_active"])
            split_rows.append(
                {
                    "dataset": dataset,
                    "seed": seed,
                    "training_single_loo_rows": len(local),
                    "wmse_gain_lower95_one_sided": float(wmse_lower),
                    "weighted_r2_delta_lower95_one_sided": float(r2_lower),
                    "shadow_active": shadow_active,
                    "original_active": original_active,
                    "concordant": shadow_active == original_active,
                }
            )

    loo = pd.DataFrame(rows)
    splits = pd.DataFrame(split_rows)
    loo.to_csv(args.out / "shadow_training_loo_rows.csv", index=False)
    splits.to_csv(args.out / "shadow_gate_by_split.csv", index=False)
    confusion = pd.crosstab(
        splits["original_active"], splits["shadow_active"], dropna=False
    ).to_dict()
    result = {
        "status": "PASS_TRAINING_ONLY_SHADOW_GATE_EXECUTION",
        "training_only": True,
        "test_or_validation_outcomes_read": False,
        "splits": len(splits),
        "exact_concordance": float(splits["concordant"].mean()),
        "active_to_inactive_flips": int(
            np.sum(splits["original_active"] & ~splits["shadow_active"])
        ),
        "inactive_to_active_flips": int(
            np.sum(~splits["original_active"] & splits["shadow_active"])
        ),
        "confusion": confusion,
        "interpretation": "descriptive counterfactual only; does not alter Gate 19",
    }
    (args.out / "shadow_gate_verdict.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
