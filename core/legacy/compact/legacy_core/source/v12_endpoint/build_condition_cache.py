#!/usr/bin/env python3
"""Compress an official scPerturBench h5ad into condition moments.

The cache contains no split assignment and is reusable across all three
official seeds. It preserves population means/variances, condition counts and
gene identity, which are sufficient for WitnessCell fitting and directional
scoring. The original h5ad remains required for distributional metrics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def condition_moments(adata: ad.AnnData, chunk_size: int = 512):
    labels = adata.obs["perturbation"].astype(str).to_numpy()
    conditions = list(pd.unique(labels))
    condition_id = {condition: index for index, condition in enumerate(conditions)}
    codes = np.asarray([condition_id[value] for value in labels], dtype=np.int64)
    sums = np.zeros((len(conditions), adata.n_vars), dtype=np.float64)
    sumsq = np.zeros_like(sums)
    counts = np.bincount(codes, minlength=len(conditions)).astype(np.int64)
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        block = adata.X[start:stop]
        if hasattr(block, "toarray"):
            block = block.toarray()
        block = np.asarray(block, dtype=np.float64)
        local = codes[start:stop]
        for code in np.unique(local):
            selected = block[local == code]
            sums[code] += selected.sum(axis=0)
            sumsq[code] += np.square(selected).sum(axis=0)
    means = sums / counts[:, None]
    variances = np.maximum(sumsq / counts[:, None] - np.square(means), 0.0)
    return conditions, counts, means, variances


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.data, backed="r")
    conditions, counts, means, variances = condition_moments(adata, args.chunk_size)
    np.savez_compressed(
        args.out,
        dataset=np.asarray(args.dataset),
        conditions=np.asarray(conditions, dtype=str),
        genes=np.asarray(adata.var_names.astype(str), dtype=str),
        counts=counts,
        means=means.astype(np.float32),
        variances=variances.astype(np.float32),
    )
    report = {
        "status": "PASS_CONDITION_MOMENT_CACHE",
        "dataset": args.dataset,
        "source": str(args.data.resolve()),
        "source_sha256": sha256(args.data),
        "source_shape": [int(adata.n_obs), int(adata.n_vars)],
        "conditions": len(conditions),
        "cache": str(args.out.resolve()),
        "cache_sha256": sha256(args.out),
        "contract": "population moments by perturbation; float32 cache generated from float64 accumulation",
    }
    args.out.with_suffix(".manifest.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
