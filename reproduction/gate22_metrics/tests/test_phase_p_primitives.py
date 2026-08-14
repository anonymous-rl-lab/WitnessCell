from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from phase_p_prediction_stress import cluster_bootstrap, comparison_point


class PhasePPrimitiveTests(unittest.TestCase):
    def test_condition_cluster_bootstrap_preserves_exact_constant_contrast(self) -> None:
        rows = []
        for dataset in ("Norman", "Wessels", "Schmidt", "Replogle_exp6"):
            for seed in (1, 2, 3):
                for condition in ("A", "B", "C"):
                    rows.extend(
                        [
                            {
                                "dataset": dataset,
                                "seed": seed,
                                "condition": condition,
                                "method": "WitnessCell_v14",
                                "wmse": 1.0,
                                "weighted_r2_deltapert": 0.75,
                            },
                            {
                                "dataset": dataset,
                                "seed": seed,
                                "condition": condition,
                                "method": "baseControl",
                                "wmse": 2.0,
                                "weighted_r2_deltapert": 0.25,
                            },
                        ]
                    )
        frame = pd.DataFrame(rows)
        point = comparison_point(frame, "baseControl")
        self.assertAlmostEqual(point["dataset_equal_wmse_ratio"], 0.5)
        self.assertAlmostEqual(point["dataset_equal_weighted_r2_difference"], 0.5)
        ratio, difference = cluster_bootstrap(frame, "baseControl")
        np.testing.assert_allclose(ratio, 0.5)
        np.testing.assert_allclose(difference, 0.5)


if __name__ == "__main__":
    unittest.main()
