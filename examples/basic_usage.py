"""Minimal, self-contained WitnessCell example."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from witnesscell import ConditionMoments, SplitSpec, WitnessCell


def main() -> None:
    rng = np.random.default_rng(17)
    genes = tuple(f"G{index}" for index in range(8))
    control = rng.normal(0.0, 0.1, len(genes))
    singles = {gene: rng.normal(0.0, 0.4, len(genes)) for gene in genes}
    means = {"control": control}
    means.update({gene: control + effect for gene, effect in singles.items()})
    pairs = []
    for left, right in list(combinations(genes, 2))[:8]:
        name = f"{left}+{right}"
        pairs.append(name)
        means[name] = control + singles[left] + singles[right]
    variances = {condition: np.full(len(genes), 0.25) for condition in means}
    counts = {condition: 80 for condition in means}
    moments = ConditionMoments.from_mappings(
        genes=genes, means=means, variances=variances, counts=counts
    )
    split = SplitSpec.create(("control", *genes, *pairs[:6]), pairs[6:])
    model = WitnessCell().fit(
        moments,
        split,
        gene2go={gene: [f"GO:{index % 3:07d}"] for index, gene in enumerate(genes)},
    )
    prediction = model.predict(["G2+G7"])
    print(prediction.means)
    print(model.diagnostics())


if __name__ == "__main__":
    main()
