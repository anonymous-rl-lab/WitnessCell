from __future__ import annotations

from itertools import combinations

import numpy as np
import pytest

from witnesscell import ConditionMoments, SplitSpec, ValidationError, WitnessCell
from witnesscell.core import (
    closed_form_gamma,
    incidence,
    one_sided_lower,
    parse_condition,
    saturate,
)


def _small_problem(
    *, train_doubles: bool, validation_doubles: bool
) -> tuple[ConditionMoments, SplitSpec, dict[str, tuple[str, ...]]]:
    genes = ("A", "B", "C")
    control = np.array([1.0, 1.1, 0.9])
    effects = {
        "A": np.array([-0.7, 0.1, 0.0]),
        "B": np.array([0.0, -0.6, 0.1]),
        "C": np.array([0.1, 0.0, -0.8]),
    }
    means = {"control": control}
    means.update({gene: control + effect for gene, effect in effects.items()})
    pairs = tuple(f"{left}+{right}" for left, right in combinations(genes, 2))
    for pair in pairs:
        left, right = pair.split("+")
        means[pair] = control + saturate(effects[left] + effects[right], 0.1)
    variances = {condition: np.full(3, 0.2) for condition in means}
    counts = {condition: 50 for condition in means}
    moments = ConditionMoments.from_mappings(
        genes=genes,
        means=means,
        variances=variances,
        counts=counts,
    )
    train = ["control", *genes]
    validation: list[str] = []
    if train_doubles:
        train.append(pairs[0])
    if validation_doubles:
        validation.append(pairs[1])
    return moments, SplitSpec(tuple(train), tuple(validation)), {gene: () for gene in genes}


@pytest.mark.parametrize(
    "condition,message",
    [
        ("", "non-empty"),
        ("A+", "invalid"),
        ("A+A", "distinct"),
        ("A+B+C", "two-endpoint"),
    ],
)
def test_condition_parser_rejects_invalid_labels(condition: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_condition(condition)


def test_small_core_fallback_paths() -> None:
    moments, split, gene2go = _small_problem(
        train_doubles=False,
        validation_doubles=False,
    )
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    prediction = model.predict(["control", "A", "A+B", "A+D"])
    assert prediction.means.shape == (4, 3)
    assert not model.state.endpoint_head.baseline_head.active
    assert np.isnan(model.state.validation_mse_witness)
    np.testing.assert_array_equal(prediction.means[0], model.state.control_mean)


def test_validation_only_and_training_only_interaction_paths() -> None:
    moments, validation_split, gene2go = _small_problem(
        train_doubles=False,
        validation_doubles=True,
    )
    validation_model = WitnessCell().fit(moments, validation_split, gene2go=gene2go)
    assert validation_model.state.noise_ratio == 0.0
    assert np.isfinite(validation_model.state.validation_mse_factorized)

    moments, training_split, gene2go = _small_problem(
        train_doubles=True,
        validation_doubles=False,
    )
    training_model = WitnessCell().fit(moments, training_split, gene2go=gene2go)
    assert training_model.state.noise_ratio == 0.1
    assert training_model.state.gamma == 1.0
    assert np.all(np.isfinite(training_model.predict(["B+C"]).means))


def test_core_scalar_helpers_and_incidence() -> None:
    effect = np.array([-2.0, 0.0, 2.0])
    np.testing.assert_allclose(saturate(effect, 0.5), [-1.0, 0.0, 1.0])
    assert closed_form_gamma(np.zeros(3), np.ones(3)) == 0.0
    assert closed_form_gamma(np.ones(3), np.full(3, 5.0)) == 1.0
    assert one_sided_lower(np.array([0.3])) == float("-inf")
    np.testing.assert_array_equal(
        incidence(["A", "A+B"], ["A", "B"]),
        [[1.0, 0.0], [1.0, 1.0]],
    )
    with pytest.raises(KeyError):
        incidence(["A+C"], ["A", "B"])


def test_model_rejects_empty_prediction() -> None:
    moments, split, gene2go = _small_problem(
        train_doubles=False,
        validation_doubles=False,
    )
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    with pytest.raises(ValidationError, match="must not be empty"):
        model.predict([])
