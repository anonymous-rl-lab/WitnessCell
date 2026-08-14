from __future__ import annotations

from itertools import combinations

import numpy as np

from witnesscell import WitnessRiskConfig, WitnessRiskEstimator, exact_witness_risk


def test_exact_witness_risk_decomposition() -> None:
    design = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    target = np.array([1.0, 0.5])
    covariance = np.array([[1.4, 0.1, 0.2], [0.1, 1.3, 0.15], [0.2, 0.15, 1.5]])
    cross = np.array([0.4, 0.25, 0.35])
    result = exact_witness_risk(design, target, covariance, cross, 1.2)
    np.testing.assert_allclose(design.T @ result.weights, target, atol=1e-12)
    assert np.isclose(result.risk, result.direct_quadratic_risk)
    assert np.isclose(result.risk, result.adequacy_term + result.geometry_term)


def test_estimated_risk_is_finite_and_target_conditioned() -> None:
    rng = np.random.default_rng(9)
    nodes = tuple(f"N{index}" for index in range(5))
    pairs = tuple(f"{left}+{right}" for left, right in combinations(nodes, 2))
    estimator = WitnessRiskEstimator(
        WitnessRiskConfig(
            pca_components=3,
            inner_splits=2,
            inner_fraction=0.2,
            length_grid=(0.5, 1.0),
            rho_grid=(0.0, 0.5),
            noise_grid=(0.1, 1.0),
        )
    ).fit(
        node_names=nodes,
        single_profiles=rng.normal(size=(5, 7)),
        train_pairs=pairs,
        residual_response=rng.normal(scale=0.2, size=(len(pairs), 6)),
    )
    prediction = estimator.predict([pairs[0], pairs[-1]])
    assert prediction.residual_means.shape == (2, 6)
    assert prediction.risks.shape == (2,)
    assert np.all(np.isfinite(prediction.risks))
    assert np.all(prediction.risks > 0)
    assert estimator.diagnostics()["candidate_count"] == 6
