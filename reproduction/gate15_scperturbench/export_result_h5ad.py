#!/usr/bin/env python3
"""Convert a target-free WitnessCell archive to scPerturBench result.h5ad."""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def dense(matrix) -> np.ndarray:
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--compression", default="gzip")
    args = parser.parse_args()
    archive = np.load(args.prediction, allow_pickle=False)
    forbidden = {"truth", "truth_variance"}.intersection(archive.files)
    if forbidden:
        raise ValueError(
            f"refusing non-deployment archive containing target arrays: {sorted(forbidden)}"
        )
    source = ad.read_h5ad(args.data, backed="r")
    conditions = archive["conditions"].astype(str)
    prediction = archive["prediction"].astype(float)
    genes = archive["genes"].astype(str)
    if not np.array_equal(genes, source.var_names.astype(str).to_numpy()):
        raise ValueError("prediction genes do not exactly match source var_names")
    labels = source.obs["perturbation"].astype(str).to_numpy()
    control_index = np.flatnonzero(labels == "control")
    stimulated_index = np.flatnonzero(np.isin(labels, conditions))
    stimulated = dense(source.X[stimulated_index]).astype(np.float32)
    control = dense(source.X[control_index]).astype(np.float32)

    rng = np.random.default_rng(1000003 + int(archive["seed"]))
    std = np.sqrt(np.maximum(archive["control_variance"].astype(float), 0.0))
    condition_id = {condition: index for index, condition in enumerate(conditions)}
    imputed = np.empty_like(stimulated)
    stimulated_labels = labels[stimulated_index]
    for condition in conditions:
        selected = np.flatnonzero(stimulated_labels == condition)
        mean = prediction[condition_id[condition]]
        imputed[selected] = rng.normal(
            loc=mean,
            scale=std,
            size=(len(selected), len(mean)),
        ).astype(np.float32)

    keep_columns = ["perturbation"]
    stimulated_obs = source.obs.iloc[stimulated_index][keep_columns].copy()
    imputed_obs = stimulated_obs.copy()
    control_obs = source.obs.iloc[control_index][keep_columns].copy()
    imputed_obs.index = [f"WitnessCell_imputed_{i}" for i in range(len(imputed_obs))]
    stimulated_obs.index = [f"WitnessCell_stimulated_{i}" for i in range(len(stimulated_obs))]
    control_obs.index = [f"WitnessCell_control_{i}" for i in range(len(control_obs))]
    imputed_obs["Expcategory"] = "imputed"
    stimulated_obs["Expcategory"] = "stimulated"
    control_obs["Expcategory"] = "control"
    output = ad.AnnData(
        X=np.vstack([imputed, stimulated, control]),
        obs=pd.concat([imputed_obs, stimulated_obs, control_obs], axis=0),
        var=source.var.copy(),
    )
    output.uns["WitnessCell"] = {
        "method": "WitnessCell",
        "seed": int(archive["seed"]),
        "target_free_prediction_archive": str(args.prediction),
        "synthetic_noise": "independent Gaussian with official-control per-gene standard deviation",
        "synthetic_rng_seed": 1000003 + int(archive["seed"]),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    output.write_h5ad(args.out, compression=args.compression)
    print({
        "status": "PASS_SCPERTURBENCH_RESULT_H5AD",
        "out": str(args.out),
        "shape": list(output.shape),
        "role_counts": output.obs.Expcategory.value_counts().to_dict(),
        "conditions": len(conditions),
    })


if __name__ == "__main__":
    main()
