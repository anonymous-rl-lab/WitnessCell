from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from phase_m_metric_validity import technical_split  # noqa: E402


class PhaseMPrimitiveTests(unittest.TestCase):
    def test_technical_split_matches_locked_rng_sequence(self) -> None:
        labels = np.asarray(["b", "a", "b", "a", "b", "c", "c", "a"], dtype=str)
        cell_ids = np.asarray([f"cell{i}" for i in range(labels.size)], dtype=str)
        observed, _ = technical_split(labels, cell_ids)

        np.random.seed(42)
        expected = np.full(labels.size, -1, dtype=np.int8)
        for condition in ["b", "a", "c"]:
            indices = np.flatnonzero(labels == condition)
            permutation = np.random.permutation(indices.size)
            split = indices.size // 2
            expected[indices[permutation[:split]]] = 0
            expected[indices[permutation[split:]]] = 1
        np.testing.assert_array_equal(observed, expected)

    def test_split_is_balanced_and_complete(self) -> None:
        labels = np.repeat(np.asarray(["a", "b", "c"]), [10, 11, 12])
        cell_ids = np.asarray([f"x{i}" for i in range(labels.size)])
        assignment, ledger = technical_split(labels, cell_ids)
        self.assertTrue(np.all(assignment >= 0))
        for row in ledger:
            self.assertLessEqual(abs(row["first_count"] - row["second_count"]), 1)
            self.assertEqual(row["first_count"] + row["second_count"], row["cells"])


if __name__ == "__main__":
    unittest.main()
