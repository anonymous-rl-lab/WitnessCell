from __future__ import annotations

from dataclasses import FrozenInstanceError
from itertools import combinations

import numpy as np
import pytest

from witnesscell import (
    ConditionMoments,
    PredictionBatch,
    SelectivePolicy,
    SplitSpec,
    ValidationError,
    WitnessCell,
    WitnessCellConfig,
    WitnessRiskConfig,
    WitnessRiskEstimator,
)
from witnesscell.selective import weighted_quantile

from ._fixture import synthetic_problem


def _moment_arguments() -> dict[str, object]:
    return {
        "genes": ("A", "B"),
        "means": {"control": np.array([1.0, 2.0])},
        "variances": {"control": np.array([0.1, 0.2])},
        "counts": {"control": 3},
    }


@pytest.mark.parametrize(
    "update,message",
    [
        ({"genes": ()}, "must not be empty"),
        ({"genes": ("A", "A")}, "unique"),
        ({"genes": ("A", " ")}, "non-empty names"),
        ({"means": {}}, "at least one"),
        ({"variances": {"other": np.zeros(2)}}, "identical conditions"),
        ({"counts": {"other": 2}}, "identical conditions"),
        ({"means": {"A+B+C": np.zeros(2)}, "variances": {"A+B+C": np.zeros(2)}, "counts": {"A+B+C": 2}}, "two-endpoint"),
        ({"means": {"control": np.zeros(3)}}, "shape"),
        ({"means": {"control": np.array([np.nan, 0.0])}}, "NaN"),
        ({"variances": {"control": np.array([-0.1, 0.2])}}, "negatives"),
        ({"counts": {"control": 3.7}}, "positive integer"),
        ({"counts": {"control": True}}, "positive integer"),
    ],
)
def test_direct_moment_constructor_is_strict(update: dict[str, object], message: str) -> None:
    arguments = _moment_arguments()
    arguments.update(update)
    with pytest.raises(ValidationError, match=message):
        ConditionMoments(**arguments)  # type: ignore[arg-type]


def test_moment_and_prediction_objects_are_defensively_immutable() -> None:
    moments = ConditionMoments(**_moment_arguments())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        moments.counts["control"] = 4  # type: ignore[index]
    with pytest.raises(ValueError):
        moments.means["control"][0] = 9.0

    matrix = np.zeros((1, 2))
    prediction = PredictionBatch(
        ("A",), ("A", "B"), matrix, matrix, matrix, matrix
    ).with_decisions([0.1], ["accept"])
    matrix[:] = 99.0
    assert prediction.means[0, 0] == 0.0
    with pytest.raises(ValueError):
        prediction.means[0, 0] = 1.0
    with pytest.raises(FrozenInstanceError):
        prediction.conditions = ("B",)  # type: ignore[misc]


def test_direct_split_constructor_rejects_leakage_and_bad_labels() -> None:
    with pytest.raises(ValidationError, match="duplicates"):
        SplitSpec(("control", "A", "A"), ())
    with pytest.raises(ValidationError, match="overlap"):
        SplitSpec(("control", "A+B"), ("A+B",))
    with pytest.raises(ValidationError, match="two-endpoint"):
        SplitSpec(("control", "A+B+C"), ())


@pytest.mark.parametrize(
    "kwargs",
    [
        {"control_label": ""},
        {"control_label": "A+B"},
        {"alpha_grid": ()},
        {"alpha_grid": (np.nan,)},
        {"alpha_grid": (0.1, 0.1)},
        {"noise_grid": (0.0,)},
        {"noise_grid": (np.inf,)},
        {"noise_grid": (0.1, 0.1)},
        {"gate_confidence": np.nan},
        {"gate_confidence": 0.5},
    ],
)
def test_model_config_rejects_nonfinite_or_ambiguous_values(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        WitnessCellConfig(**kwargs)  # type: ignore[arg-type]


def test_model_fit_revalidates_public_boundaries() -> None:
    moments, split, gene2go, _pairs = synthetic_problem()
    with pytest.raises(ValidationError, match="moments"):
        WitnessCell().fit(object(), split, gene2go=gene2go)  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="split"):
        WitnessCell().fit(moments, object(), gene2go=gene2go)  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="gene2go"):
        WitnessCell().fit(moments, split, gene2go=[])  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="keys"):
        WitnessCell().fit(moments, split, gene2go={"": ()})
    with pytest.raises(ValidationError, match="iterable"):
        WitnessCell().fit(moments, split, gene2go={"G0": "GO:1"})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"pca_components": 0},
        {"inner_splits": 0},
        {"inner_fraction": 0.0},
        {"random_state": -1},
        {"length_grid": (np.nan,)},
        {"rho_grid": (1.1,)},
        {"noise_grid": (0.0,)},
        {"noise_grid": (0.1, 0.1)},
    ],
)
def test_risk_config_rejects_invalid_values(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        WitnessRiskConfig(**kwargs)  # type: ignore[arg-type]


def test_risk_estimator_rejects_nan_reliability_and_bad_predictions() -> None:
    rng = np.random.default_rng(3)
    nodes = tuple(f"N{index}" for index in range(5))
    pairs = tuple(f"{left}+{right}" for left, right in combinations(nodes, 2))
    arguments = {
        "node_names": nodes,
        "single_profiles": rng.normal(size=(5, 4)),
        "train_pairs": pairs,
        "residual_response": rng.normal(size=(10, 3)),
    }
    with pytest.raises(ValidationError, match="reliability"):
        WitnessRiskEstimator().fit(**arguments, reliability=np.full(10, np.nan))
    estimator = WitnessRiskEstimator(
        WitnessRiskConfig(
            inner_splits=1,
            length_grid=(1.0,),
            rho_grid=(0.0,),
            noise_grid=(0.1,),
        )
    ).fit(**arguments)
    with pytest.raises(ValidationError, match="must not be empty"):
        estimator.predict([])
    with pytest.raises(ValidationError, match="absent"):
        estimator.predict(["N0+missing"])


@pytest.mark.parametrize(
    "values,quantile,weights",
    [
        ([], 0.5, []),
        ([0.1, np.nan], 0.5, [1.0, 1.0]),
        ([0.1], np.nan, [1.0]),
        ([0.1], 0.5, [-1.0]),
    ],
)
def test_weighted_quantile_rejects_invalid_inputs(
    values: list[float], quantile: float, weights: list[float]
) -> None:
    with pytest.raises(ValidationError):
        weighted_quantile(values, quantile, weights)


def test_selective_calibration_rejects_nan_and_zero_coverage() -> None:
    with pytest.raises(ValidationError):
        SelectivePolicy.calibrate([0.1, np.nan])
    with pytest.raises(ValidationError):
        SelectivePolicy.calibrate([0.1], coverage=0.0)
    with pytest.raises(ValidationError):
        SelectivePolicy(0.1, "")
