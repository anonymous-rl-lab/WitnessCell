#!/usr/bin/env python3
"""Synthetic truth test for the target-conditioned witness law."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge


def incidence(edges: np.ndarray, n: int) -> np.ndarray:
    b = np.zeros((len(edges), n))
    b[np.arange(len(edges)), edges[:, 0]] = 1.0
    b[np.arange(len(edges)), edges[:, 1]] = 1.0
    return b


def graph_information(edges: np.ndarray, n: int) -> dict:
    b = incidence(edges, n)
    q = b.T @ b
    eig = np.linalg.eigvalsh(q)
    return {
        "q_nullity": int(np.sum(eig <= 1e-9)),
        "q_min_eig": float(eig.min()),
        "q_trace_inv_1e3": float(np.trace(np.linalg.inv(q + 1e-3 * np.eye(n)))),
    }


def fit_anchor(edges: np.ndarray, y: np.ndarray, n: int) -> np.ndarray:
    b = incidence(edges, n)
    return Ridge(alpha=1e-5, fit_intercept=False).fit(b, y).coef_.T


def one_seed(seed: int) -> list[dict]:
    rng = np.random.default_rng(20260807 + seed)
    n, d = 40, 48
    left = np.arange(n // 2)
    right = np.arange(n // 2, n)
    cross = np.asarray([(i, j) for i in left for j in right], dtype=int)
    within = np.asarray(
        [(i, j) for group in (left, right) for z, i in enumerate(group) for j in group[z + 1:]],
        dtype=int,
    )
    base = rng.normal(scale=0.5, size=(n, d))
    contrast = rng.normal(scale=1.5, size=d)
    sign = np.r_[np.ones(len(left)), -np.ones(len(right))]
    truth_a = base + sign[:, None] * contrast
    all_e = np.vstack([cross, within])
    true_all = truth_a[all_e[:, 0]] + truth_a[all_e[:, 1]]
    sigma = 0.10 * math.sqrt(float(np.mean(true_all * true_all)))
    designs = {
        "dense_bipartite_400": cross,
        "odd_witness_1": np.vstack([cross, within[[0]]]),
        "odd_witness_4": np.vstack([cross, within[[0, 40, 190, 230]]]),
        "odd_witness_16": np.vstack([cross, within[rng.choice(len(within), 16, replace=False)]]),
        "random_120": all_e[rng.choice(len(all_e), 120, replace=False)],
    }
    rows = []
    for name, train_e in designs.items():
        train_y = truth_a[train_e[:, 0]] + truth_a[train_e[:, 1]]
        train_y += rng.normal(scale=sigma, size=train_y.shape)
        ahat = fit_anchor(train_e, train_y, n)
        pred = ahat[all_e[:, 0]] + ahat[all_e[:, 1]]
        rows.append({
            "seed": seed,
            "design": name,
            "edges": len(train_e),
            "corr": float(np.corrcoef(pred.ravel(), true_all.ravel())[0, 1]),
            "nmse": float(np.mean((pred - true_all) ** 2) / np.mean(true_all ** 2)),
            **graph_information(train_e, n),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--out", type=Path, default=Path("results"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([row for s in range(args.seeds) for row in one_seed(s)])
    summary = frame.groupby("design").agg(
        corr=("corr", "mean"), nmse=("nmse", "mean"),
        q_nullity=("q_nullity", "mean"), edges=("edges", "mean"), n=("seed", "size"),
    ).reset_index()
    frame.to_csv(args.out / "topology_rows.csv", index=False)
    summary.to_csv(args.out / "topology_summary.csv", index=False)
    verdict = {
        "law": "Exposure is not identifiability: target recovery is governed by the unsigned incidence row space and signless-Laplacian conditioning.",
        "dense_bipartite_corr": float(summary.loc[summary.design == "dense_bipartite_400", "corr"].iloc[0]),
        "odd_witness_1_corr": float(summary.loc[summary.design == "odd_witness_1", "corr"].iloc[0]),
    }
    (args.out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
