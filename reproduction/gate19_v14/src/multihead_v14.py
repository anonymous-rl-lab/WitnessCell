#!/usr/bin/env python3
"""Incremental amplitude-residual head on top of frozen WitnessCell v13.

The v13 endpoint head is never replaced.  A three-scalar, one-coordinate
correction is learned from training singles only:

* generic on-target self response;
* sign-matched tail of the target gene response fingerprint;
* sign-matched RMS of that fingerprint.

The correction is active only when an honest outer LOO shows positive
one-sided 95% lower bounds for both all-gene MSE gain and top-100 PCC delta
relative to a v13 head refit without the outer endpoint.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from dual_head import DualHeadFit, fit_dual_head
from identity_head import EPS, one_sided_lower, pearson, welch_scores
from probe_multhead_training import (
    Bank,
    Recipe,
    build_bank,
    solve,
    sparse_features,
)


RECIPE = Recipe("inc_amp_fp_r02", "inc_amp_fp", ridge=0.2)


def fit_frozen_v13(
    effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    gene2go: Mapping,
) -> DualHeadFit:
    return fit_dual_head(
        single_effects=effects,
        single_means={node: single_means[node] for node in effects},
        single_variances={node: single_variances[node] for node in effects},
        single_counts={node: single_counts[node] for node in effects},
        control_mean=control_mean,
        control_variance=control_variance,
        control_count=control_count,
        genes=genes,
        gene2go=gene2go,
        dense_mode="mean_go",
        sparse_mode="joint_self",
        support_mode="go_residual",
        support_k=1,
        go_top_k=10,
        nested_gate=False,
        direction_gate_metric="top100_pcc",
        require_incremental_direction=False,
    )


def top100(
    node: str,
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
) -> np.ndarray:
    score = welch_scores(
        single_means[node], single_variances[node], single_counts[node],
        control_mean, control_variance, control_count,
    )
    return np.argsort(-np.abs(score), kind="stable")[: min(100, len(score))]


def fit_correction(
    nodes: Sequence[str],
    effects: Mapping[str, np.ndarray],
    baseline: DualHeadFit,
    genes: Sequence[str],
    gene2go: Mapping,
    top_indices: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, list[Bank]]:
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    banks = []
    residuals = []
    tops = []
    for held in nodes:
        remaining = {node: effects[node] for node in nodes if node != held}
        bank = build_bank(
            held, remaining, genes, gene2go, 10, (5, 10, 20, 40), (3, 5, 10)
        )
        banks.append(bank)
        residuals.append(
            effects[held] - baseline.predict(held, remaining, gene_index, gene2go)
        )
        tops.append(top_indices[held])
    weight = solve(
        [sparse_features(bank, RECIPE) for bank in banks],
        residuals,
        indices=tops,
        lower=-1.5,
        upper=1.5,
        ridge=1e-4,
        relative_ridge=RECIPE.ridge,
    )
    return weight, banks


@dataclass(frozen=True)
class IncrementalAmplitudeHead:
    baseline_head: DualHeadFit
    correction_active: bool
    correction_weights: tuple[float, ...]
    all_gene_upgrade_mean: float
    all_gene_upgrade_lower: float
    top100_mse_upgrade_mean: float
    top100_mse_upgrade_lower: float
    top100_pcc_upgrade_mean: float
    top100_pcc_upgrade_lower: float
    known_single_count: int
    records: tuple[dict, ...]
    ridge: float = 0.2
    selected_recipe: str = "inc_amp_fp_r02"

    @property
    def active(self) -> bool:
        return self.baseline_head.active or self.correction_active

    @property
    def dense_active(self) -> bool:
        return self.baseline_head.dense_active

    @property
    def sparse_active(self) -> bool:
        return self.correction_active

    def predict(
        self,
        node: str,
        effects: Mapping[str, np.ndarray],
        gene_index: Mapping[str, int],
        gene2go: Mapping,
    ) -> np.ndarray:
        baseline = self.baseline_head.predict(node, effects, gene_index, gene2go)
        if not self.correction_active:
            return baseline
        genes = [None] * len(gene_index)
        for gene, index in gene_index.items():
            genes[index] = gene
        bank = build_bank(
            node, effects, genes, gene2go, 10, (5, 10, 20, 40), (3, 5, 10)
        )
        correction = sum(
            weight * feature
            for weight, feature in zip(
                self.correction_weights, sparse_features(bank, RECIPE), strict=True
            )
        )
        return baseline + correction


def fit_incremental_amplitude_head(
    single_effects: Mapping[str, np.ndarray],
    single_means: Mapping[str, np.ndarray],
    single_variances: Mapping[str, np.ndarray],
    single_counts: Mapping[str, int],
    control_mean: np.ndarray,
    control_variance: np.ndarray,
    control_count: int,
    genes: Sequence[str],
    gene2go: Mapping,
    gate_confidence: float = 0.95,
) -> IncrementalAmplitudeHead:
    gene_index = {str(gene): index for index, gene in enumerate(genes)}
    nodes = sorted(node for node in single_effects if node in gene_index)
    baseline = fit_frozen_v13(
        single_effects, single_means, single_variances, single_counts,
        control_mean, control_variance, control_count, genes, gene2go,
    )
    if len(nodes) < 6:
        return IncrementalAmplitudeHead(
            baseline, False, (0.0, 0.0, 0.0), 0.0, float("-inf"),
            0.0, float("-inf"), 0.0, float("-inf"), len(nodes), (),
        )
    top_indices = {
        node: top100(
            node, single_means, single_variances, single_counts,
            control_mean, control_variance, control_count,
        )
        for node in nodes
    }

    # Honest outer LOO evidence relative to a v13 model refit after removing
    # the outer endpoint.  This is slower than the final fit but prevents the
    # correction gate from sharing the outer endpoint through v13 weights.
    records = []
    for outer in nodes:
        remaining_nodes = [node for node in nodes if node != outer]
        remaining_effects = {node: single_effects[node] for node in remaining_nodes}
        outer_baseline = fit_frozen_v13(
            remaining_effects, single_means, single_variances, single_counts,
            control_mean, control_variance, control_count, genes, gene2go,
        )
        correction_weight, _ = fit_correction(
            remaining_nodes, remaining_effects, outer_baseline, genes, gene2go,
            top_indices,
        )
        outer_bank = build_bank(
            outer, remaining_effects, genes, gene2go, 10,
            (5, 10, 20, 40), (3, 5, 10),
        )
        base = outer_baseline.predict(
            outer, remaining_effects, gene_index, gene2go
        )
        correction = sum(
            weight * feature
            for weight, feature in zip(
                correction_weight, sparse_features(outer_bank, RECIPE), strict=True
            )
        )
        final = base + correction
        truth = single_effects[outer]
        top = top_indices[outer]
        base_all = float(np.mean(np.square(base - truth)))
        final_all = float(np.mean(np.square(final - truth)))
        base_top = float(np.mean(np.square(base[top] - truth[top])))
        final_top = float(np.mean(np.square(final[top] - truth[top])))
        records.append({
            "condition": outer,
            "v13_all_gene_mse": base_all,
            "v14_all_gene_mse": final_all,
            "upgrade_all_gene_gain": (base_all - final_all) / max(base_all, EPS),
            "v13_top100_mse": base_top,
            "v14_top100_mse": final_top,
            "upgrade_top100_mse_gain": (base_top - final_top) / max(base_top, EPS),
            "v13_top100_pcc": pearson(base[top], truth[top]),
            "v14_top100_pcc": pearson(final[top], truth[top]),
            "upgrade_top100_pcc_delta": pearson(final[top], truth[top]) - pearson(base[top], truth[top]),
        })

    all_gain = np.asarray([row["upgrade_all_gene_gain"] for row in records])
    top_gain = np.asarray([row["upgrade_top100_mse_gain"] for row in records])
    pcc_delta = np.asarray([row["upgrade_top100_pcc_delta"] for row in records])
    all_lower = one_sided_lower(all_gain, gate_confidence)
    top_lower = one_sided_lower(top_gain, gate_confidence)
    pcc_lower = one_sided_lower(pcc_delta, gate_confidence)
    active = bool(all_lower > 0.0 and pcc_lower > 0.0)

    final_weight, _ = fit_correction(
        nodes, single_effects, baseline, genes, gene2go, top_indices
    )
    return IncrementalAmplitudeHead(
        baseline_head=baseline,
        correction_active=active,
        correction_weights=tuple(map(float, final_weight)),
        all_gene_upgrade_mean=float(all_gain.mean()),
        all_gene_upgrade_lower=float(all_lower),
        top100_mse_upgrade_mean=float(top_gain.mean()),
        top100_mse_upgrade_lower=float(top_lower),
        top100_pcc_upgrade_mean=float(pcc_delta.mean()),
        top100_pcc_upgrade_lower=float(pcc_lower),
        known_single_count=len(nodes),
        records=tuple(records),
    )
