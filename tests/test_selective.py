from __future__ import annotations

import numpy as np
import pytest

from witnesscell import FROZEN_NORMAN_THRESHOLD, SelectivePolicy, ValidationError


def test_calibrated_policy_and_frozen_policy() -> None:
    policy = SelectivePolicy.calibrate([0.1, 0.2, 0.3], coverage=0.5, pair_weights=[1, 2, 1])
    assert policy.threshold == 0.2
    assert policy.decide([0.1, 0.2, 0.21]) == ("accept", "accept", "abstain")
    assert SelectivePolicy.frozen_norman().threshold == FROZEN_NORMAN_THRESHOLD


def test_policy_rejects_invalid_risk() -> None:
    with pytest.raises(ValidationError):
        SelectivePolicy(0.2).decide(np.array([np.nan]))
