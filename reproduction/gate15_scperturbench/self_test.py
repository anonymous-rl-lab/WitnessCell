#!/usr/bin/env python3
"""Fast positive, negative, and target-isolation controls."""
from __future__ import annotations

import numpy as np

from witnesscell_core import fit_predict


def synthetic_means(seed: int = 7) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    genes = 64
    control = rng.normal(1.0, 0.1, genes)
    node_effect = {node: rng.normal(0, 0.2, genes) for node in "ABCDE"}
    anchor = {node: rng.normal(0, 0.08, genes) for node in "ABCDE"}
    means = {"control": control}
    for node in "ABCDE":
        means[node] = control + node_effect[node]
    for left, right in ("AB", "AC", "BD", "CD", "AE", "BE"):
        means[f"{left}+{right}"] = (
            control + node_effect[left] + node_effect[right]
            + anchor[left] + anchor[right]
        )
    return means


def main() -> None:
    means = synthetic_means()
    train = ["control", "A", "B", "C", "D", "A+B", "A+C", "B+D", "C+D"]
    validation = ["A+E"]
    target = ["B+E"]
    restricted = {key: means[key] for key in [*train, *validation]}
    fitted = fit_predict(restricted, train, validation, target)
    truth = means["B+E"]
    base = fitted.factorized_predictions["B+E"]
    witness = fitted.predictions["B+E"]
    base_mse = float(np.mean((base - truth) ** 2))
    witness_mse = float(np.mean((witness - truth) ** 2))
    assert witness_mse < base_mse, (witness_mse, base_mse)

    # Target isolation: the fit API never receives B+E, hence replacing its
    # outcome cannot change a prediction.
    changed = dict(means)
    changed["B+E"] = changed["B+E"] + 1000.0
    fitted_again = fit_predict(restricted, train, validation, target)
    assert np.array_equal(
        fitted.predictions["B+E"], fitted_again.predictions["B+E"]
    )

    # Zero-witness degradation: an official split with no training double may
    # not fail or infer a residual from test targets.  It must reduce exactly
    # to the factorized backbone with gamma=0.
    no_witness_fit = fit_predict(
        {key: means[key] for key in ("control", "A", "B", "C")},
        ["control", "A", "B", "C"],
        [],
        ["A+B"],
    )
    assert no_witness_fit.gamma == 0.0
    assert no_witness_fit.training_doubles == ()
    assert np.array_equal(
        no_witness_fit.predictions["A+B"],
        no_witness_fit.factorized_predictions["A+B"],
    )
    print({
        "status": "PASS_SELF_TEST",
        "positive_control_relative_mse_gain": (base_mse - witness_mse) / base_mse,
        "target_isolation": True,
        "zero_witness_degradation": True,
    })


if __name__ == "__main__":
    main()
