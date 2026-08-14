from __future__ import annotations

import numpy as np

from witnesscell import WitnessCell

from ._fixture import synthetic_problem


def test_frozen_v14_golden_prediction() -> None:
    moments, split, gene2go, pairs = synthetic_problem()
    model = WitnessCell().fit(moments, split, gene2go=gene2go)
    prediction = model.predict(pairs[10:12])
    expected = np.array(
        [
            [
                0.6758489911041444,
                0.19122250020563591,
                0.16859062828736932,
                0.39568983196822666,
                0.1860371344775461,
                -0.35340761850787106,
                0.23099523255178053,
                0.10713849635570694,
            ],
            [
                0.7424210336017792,
                -0.035090584529302915,
                0.2566420224594557,
                0.22897391713487192,
                0.44959476389385794,
                0.5734561748305553,
                -0.15240046823877057,
                0.3035825624955683,
            ],
        ]
    )
    np.testing.assert_allclose(prediction.means, expected, rtol=0.0, atol=2e-12)
    diagnostics = model.diagnostics()
    assert diagnostics["package_contract"] == "witnesscell-v14-frozen"
    assert diagnostics["dense_active"] is True
    assert diagnostics["sparse_active"] is True
    assert diagnostics["amplitude_active"] is False
    assert diagnostics["alpha"] == 0.0
    assert diagnostics["noise_ratio"] == 0.01
    assert diagnostics["gamma"] == 0.0
