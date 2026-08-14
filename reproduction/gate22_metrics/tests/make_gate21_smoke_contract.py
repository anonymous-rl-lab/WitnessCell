#!/usr/bin/env python3
"""Create a synthetic Gate 21 contract for Phase D engineering smoke only."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


THRESHOLD = 0.0923227147328771


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20260811)
    pairs = [f"G{i:02d}+H{i:02d}" for i in range(33)]
    pair_rows = np.repeat(pairs, 6)
    n = len(pair_rows)
    genes = np.asarray([f"gene_{index:03d}" for index in range(120)], dtype=str)
    risk = np.linspace(0.03, 0.16, n)
    accepted = risk <= THRESHOLD
    truth = rng.normal(size=(n, len(genes)))
    noise = rng.normal(size=truth.shape)
    scale = np.where(accepted, 0.15, 0.55)[:, None]
    estimated = truth + scale * noise
    geometry = truth + 0.40 * rng.normal(size=truth.shape)
    weights = np.square(rng.uniform(0.05, 1.0, size=truth.shape))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        query_ids=np.asarray([f"{100 + i % 6:03d}::{pair}" for i, pair in enumerate(pair_rows)], dtype=str),
        seeds=np.asarray([100 + i % 6 for i in range(n)], dtype=np.int64),
        pairs=np.asarray(pair_rows, dtype=str),
        genes=genes,
        estimated_prediction=estimated,
        geometry_prediction=geometry,
        truth=truth,
        canonical_weights=weights,
        scoring_evaluable=np.ones(n, dtype=bool),
        estimated_witness_risk=risk,
        geometry_risk=np.linspace(0.02, 0.18, n)[::-1],
        accepted_primary=accepted,
        threshold=np.asarray(THRESHOLD, dtype=np.float64),
    )


if __name__ == "__main__":
    main()
