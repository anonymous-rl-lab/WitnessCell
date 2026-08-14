#!/usr/bin/env python3
"""Create a small deterministic non-formal h5ad for Experiment 22 smoke tests."""

from __future__ import annotations

import argparse
import pickle
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with (args.gears_assets / "gene2go_all.pkl").open("rb") as handle:
        gene2go = pickle.load(handle)
    with (args.gears_assets / "essential_all_data_pert_genes.pkl").open("rb") as handle:
        essential = pickle.load(handle)
    targets = sorted(map(str, set(gene2go).intersection(essential)))[:10]
    extra = [f"smoke_gene_{index:03d}" for index in range(190)]
    genes = targets + extra
    conditions = ["control", *targets, *(f"{left}+{right}" for left, right in combinations(targets, 2))]
    rng = np.random.RandomState(20260811)
    cells_per_condition = 10
    matrix = []
    labels = []
    for condition in conditions:
        values = rng.normal(0.0, 0.25, size=(cells_per_condition, len(genes))).astype(np.float32)
        if condition != "control":
            members = condition.split("+")
            for member in members:
                values[:, genes.index(member)] += 1.5
            interaction = sum(genes.index(member) for member in members) % len(extra)
            values[:, len(targets) + interaction] += 0.5 * len(members)
        matrix.append(values)
        labels.extend([condition] * cells_per_condition)
    obs = pd.DataFrame(
        {"perturbation": pd.Categorical(labels)},
        index=[f"smoke_cell_{index:05d}" for index in range(len(labels))],
    )
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=np.vstack(matrix), obs=obs, var=var)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"SMOKE_ASSET_READY {args.out} shape={adata.shape} conditions={len(conditions)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
