from __future__ import annotations

import unittest

import numpy as np

from comparator_core import regenerate_mean_baselines


class ComparatorCoreTests(unittest.TestCase):
    def test_released_trainmean_branch_semantics(self) -> None:
        names = ["control", "A", "B", "C+D"]
        means = np.asarray([[1.0, 2.0], [3.0, 4.0], [7.0, 8.0], [9.0, 10.0]])
        counts = np.asarray([10, 1, 3, 5])
        targets = ["A+B", "A+Z", "Y+B", "Y+Z", "Z"]
        train_mean, base_control = regenerate_mean_baselines(
            targets, ["control", "A", "B", "C+D"], names, means, counts, means[0]
        )
        equal = (means[1] + means[2]) / 2.0
        cell_weighted = (means[1] + 3 * means[2]) / 4.0
        expected = np.stack(
            [
                means[1] + means[2] - means[0],
                means[1] + equal - means[0],
                means[2] + equal - means[0],
                2 * equal - means[0],
                cell_weighted,
            ]
        )
        np.testing.assert_allclose(train_mean, expected)
        np.testing.assert_array_equal(base_control, np.repeat(means[[0]], len(targets), axis=0))


if __name__ == "__main__":
    unittest.main()
