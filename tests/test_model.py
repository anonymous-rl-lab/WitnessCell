from __future__ import annotations

import numpy as np
import pytest

from witnesscell import NotFittedError, ValidationError, WitnessCell

from ._fixture import synthetic_problem


def test_fit_predict_shapes_and_diagnostics() -> None:
    moments, split, gene2go, pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    prediction = model.predict((pairs[10], pairs[11], "G0"))
    assert prediction.means.shape == (3, len(moments.genes))
    assert prediction.factorized_means.shape == prediction.means.shape
    assert np.all(np.isfinite(prediction.means))
    diagnostics = model.diagnostics()
    assert diagnostics["package_contract"] == "witnesscell-v14-frozen"
    assert diagnostics["known_single_count"] == 8
    assert diagnostics["training_double_count"] == 7
    assert diagnostics["validation_double_count"] == 3


def test_final_target_outcome_is_not_consumed() -> None:
    moments, split, gene2go, pairs = synthetic_problem()
    target = pairs[10]
    first = WitnessCell().fit(moments, split, gene2go=gene2go)
    changed_means = dict(moments.means)
    changed_means[target] = changed_means[target] + 1_000_000.0
    changed = type(moments).from_mappings(
        genes=moments.genes,
        means=changed_means,
        variances=moments.variances,
        counts=moments.counts,
    )
    second = WitnessCell().fit(changed, split, gene2go=gene2go)
    np.testing.assert_array_equal(first.predict([target]).means, second.predict([target]).means)
    assert first.diagnostics() == second.diagnostics()


def test_predict_before_fit_and_invalid_condition() -> None:
    with pytest.raises(NotFittedError):
        WitnessCell().predict(["G0+G1"])
    moments, split, gene2go, _pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    with pytest.raises(ValueError, match="two-endpoint"):
        model.predict(["G0+G1+G2"])


def test_split_overlap_is_rejected() -> None:
    from witnesscell import SplitSpec

    with pytest.raises(ValidationError, match="overlap"):
        SplitSpec.create(["control", "G0"], ["G0"])
