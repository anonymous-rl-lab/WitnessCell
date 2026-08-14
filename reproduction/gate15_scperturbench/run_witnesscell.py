#!/usr/bin/env python3
"""Run WitnessCell on one official scPerturBench genetic dataset."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)
from witnesscell_core import fit_predict


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
    codes = np.asarray([condition_id[value] for value in labels], dtype=int)
    sums = np.zeros((len(conditions), adata.n_vars), dtype=np.float64)
    sumsq = np.zeros_like(sums)
    counts = np.bincount(codes, minlength=len(conditions)).astype(int)
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
    parser.add_argument("--seed", type=int, required=True, choices=(1, 2, 3))
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.data, backed="r")
    conditions, counts, matrix, variances = condition_moments(adata)
    means = {condition: matrix[index] for index, condition in enumerate(conditions)}
    supported_genes = load_gears_supported_genes(args.gears_assets)
    official_conditions = filter_gears_supported_conditions(conditions, supported_genes)
    split = make_official_split(official_conditions, args.seed)
    fit_condition_names = {"control", *split.train, *split.validation}
    fit_means = {
        condition: means[condition]
        for condition in fit_condition_names
    }
    result = fit_predict(
        fit_means,
        list(split.train),
        list(split.validation),
        list(split.test),
    )
    prediction = np.stack([result.predictions[c] for c in split.test]).astype(np.float32)
    factorized = np.stack([result.factorized_predictions[c] for c in split.test]).astype(np.float32)
    truth = np.stack([means[c] for c in split.test]).astype(np.float32)
    truth_variance = np.stack([
        variances[conditions.index(c)] for c in split.test
    ]).astype(np.float32)
    control = means["control"].astype(np.float32)
    gene_names = np.asarray(adata.var_names.astype(str), dtype=str)
    test_counts = np.asarray([counts[conditions.index(c)] for c in split.test], dtype=np.int64)
    subgroup = {}
    for name, values in split.test_subgroup.items():
        subgroup.update({condition: name for condition in values})
    subgroup_array = np.asarray([subgroup[c] for c in split.test], dtype=str)
    np.savez_compressed(
        args.out / "predictions.npz",
        dataset=np.asarray(args.dataset),
        seed=np.asarray(args.seed, dtype=np.int64),
        conditions=np.asarray(split.test, dtype=str),
        subgroups=subgroup_array,
        genes=gene_names,
        prediction=prediction,
        factorized_prediction=factorized,
        truth=truth,
        truth_variance=truth_variance,
        control=control,
        control_variance=variances[conditions.index("control")].astype(np.float32),
        control_count=np.asarray(counts[conditions.index("control")], dtype=np.int64),
        test_cell_counts=test_counts,
    )
    # Submission-side archive: deliberately contains no target expression.
    np.savez_compressed(
        args.out / "deploy_predictions.npz",
        dataset=np.asarray(args.dataset),
        seed=np.asarray(args.seed, dtype=np.int64),
        conditions=np.asarray(split.test, dtype=str),
        subgroups=subgroup_array,
        genes=gene_names,
        prediction=prediction,
        factorized_prediction=factorized,
        control=control,
        control_variance=variances[conditions.index("control")].astype(np.float32),
        control_count=np.asarray(counts[conditions.index("control")], dtype=np.int64),
        test_cell_counts=test_counts,
    )
    manifest = {
        "status": "PASS_WITNESSCELL_MEAN_PREDICTION",
        "method": "WitnessCell",
        "dataset": args.dataset,
        "seed": args.seed,
        "data": str(args.data.resolve()),
        "data_sha256": sha256(args.data),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "official_gears_filter": {
            "raw_conditions": len(conditions),
            "supported_conditions": len(official_conditions),
            "removed_conditions": sorted(set(conditions) - set(official_conditions)),
        },
        "split_counts": {
            "train": len(split.train),
            "validation": len(split.validation),
            "test": len(split.test),
        },
        "test_subgroup_counts": {
            key: len(value) for key, value in split.test_subgroup.items()
        },
        "selected": {
            "alpha": result.alpha,
            "noise_ratio": result.noise_ratio,
            "gamma": result.gamma,
            "validation_mse_factorized": result.validation_mse_factorized,
            "validation_mse_witness": result.validation_mse_witness,
        },
        "training_doubles": len(result.training_doubles),
        "validation_doubles": len(result.validation_doubles),
        "known_training_singles": len(result.known_single_genes),
        "leakage_contract": "witnesscell_core.fit_predict receives a dictionary restricted to control, official train, and official validation conditions; deploy_predictions.npz contains no target expression",
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
