from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.metric_core import (  # noqa: E402
    MetricInputError,
    WeightVector,
    assert_gene_order,
    drf,
    nir,
    source_weight_transform,
    weighted_delta_r2,
    wmse,
)


class MetricCoreTests(unittest.TestCase):
    def test_source_weight_transform_and_duplicate_max(self) -> None:
        weights = source_weight_transform(
            scores=[-1.0, 3.0, 2.0, 5.0],
            score_genes=["g1", "g2", "g2", "unused"],
            evaluation_genes=["g2", "g1", "g3"],
        )
        # abs scores min/max 1/5 -> squares: 0, .25, .0625, 1;
        # duplicate g2 keeps .25 before alignment.
        np.testing.assert_allclose(weights.values, [0.25, 0.0, 0.0])

    def test_perfect_prediction(self) -> None:
        weights = WeightVector(np.array(["a", "b", "c"]), np.array([1.0, 2.0, 3.0]))
        truth = np.array([2.0, -1.0, 4.0])
        baseline = np.array([0.0, 0.0, 0.0])
        self.assertEqual(wmse(truth, truth, weights), 0.0)
        self.assertAlmostEqual(weighted_delta_r2(truth, truth, baseline, weights), 1.0)

    def test_zero_sum_weights_fail_closed(self) -> None:
        with self.assertRaises(MetricInputError):
            WeightVector(np.array(["a", "b"]), np.array([0.0, 0.0]))
        with self.assertRaises(MetricInputError):
            source_weight_transform([1.0, 1.0], ["a", "b"], ["a", "b"])

    def test_nan_fails_closed(self) -> None:
        weights = WeightVector(np.array(["a", "b"]), np.array([1.0, 1.0]))
        with self.assertRaises(MetricInputError):
            wmse([1.0, np.nan], [1.0, 2.0], weights)
        with self.assertRaises(MetricInputError):
            source_weight_transform([np.nan, np.nan], ["a", "b"], ["a", "b"])

    def test_gene_reordering_and_duplicates_fail_closed(self) -> None:
        with self.assertRaises(MetricInputError):
            assert_gene_order(["a", "b"], ["b", "a"])
        with self.assertRaises(MetricInputError):
            assert_gene_order(["a", "b"], ["a", "a"])

    def test_nir_strict_comparison_and_ties(self) -> None:
        truth = np.array([[0.0], [2.0], [4.0]])
        perfect = nir(truth.copy(), truth, ["c_a", "c_b", "c_c"])
        self.assertEqual(perfect, {"c_a": 1.0, "c_b": 1.0, "c_c": 1.0})
        tied = nir(np.array([[1.0], [2.0], [4.0]]), truth, ["c_a", "c_b", "c_c"])
        # For c_a, truth a and b are tied; strict comparison makes b a loss.
        self.assertEqual(tied["c_a"], 0.5)

    def test_nir_omits_singleton_covariate(self) -> None:
        output = nir(
            np.array([[0.0], [1.0], [9.0]]),
            np.array([[0.0], [1.0], [9.0]]),
            ["a_x", "a_y", "b_z"],
        )
        self.assertEqual(set(output), {"a_x", "a_y"})

    def test_drf_locked_branches_and_clipping(self) -> None:
        self.assertAlmostEqual(drf(4.0, 1.0, higher_better=False), 3.0 / 4.000001)
        self.assertAlmostEqual(drf(0.2, 0.8, higher_better=True), 0.6 / 0.800001)
        self.assertEqual(drf(1.0, 3.0, higher_better=False), -1.0)
        self.assertEqual(drf(0.5, -2.0, higher_better=True), -1.0)


if __name__ == "__main__":
    unittest.main()

