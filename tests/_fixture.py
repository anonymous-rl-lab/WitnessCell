"""Deterministic synthetic fixtures shared by package tests."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from witnesscell import ConditionMoments, SplitSpec


def synthetic_problem() -> tuple[ConditionMoments, SplitSpec, dict[str, list[str]], tuple[str, ...]]:
    rng = np.random.default_rng(20260814)
    genes = tuple(f"G{index}" for index in range(8))
    control = rng.normal(0.5, 0.05, len(genes))
    single_effects: dict[str, np.ndarray] = {}
    for index, gene in enumerate(genes):
        effect = rng.normal(0.0, 0.22, len(genes))
        effect[index] -= 0.9 + 0.03 * index
        single_effects[gene] = effect
    means = {"control": control}
    means.update({gene: control + effect for gene, effect in single_effects.items()})
    pairs = tuple(f"{left}+{right}" for left, right in combinations(genes, 2))
    for pair_index, pair in enumerate(pairs[:12]):
        left, right = pair.split("+")
        additive = single_effects[left] + single_effects[right]
        saturation = additive / (1.0 + 0.1 * np.abs(additive))
        interaction = 0.04 * (pair_index + 1) * np.mean([single_effects[left], single_effects[right]], axis=0)
        means[pair] = control + saturation + interaction
    variances = {
        condition: np.full(len(genes), 0.18 + 0.001 * index)
        for index, condition in enumerate(means)
    }
    counts = {condition: 90 + index for index, condition in enumerate(means)}
    moments = ConditionMoments.from_mappings(
        genes=genes, means=means, variances=variances, counts=counts
    )
    split = SplitSpec.create(
        ("control", *genes, *pairs[:7]),
        pairs[7:10],
    )
    gene2go = {
        gene: [f"GO:{index % 3:07d}", f"GO:{(index + 1) % 4:07d}"]
        for index, gene in enumerate(genes)
    }
    return moments, split, gene2go, pairs
