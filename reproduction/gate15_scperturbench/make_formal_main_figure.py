#!/usr/bin/env python3
"""Create the frozen three-panel formal scPerturBench result figure."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--aggregate", type=Path, default=Path("results/formal_score/aggregate"))
    parser.add_argument("--out", type=Path, default=Path("figures"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    board = pd.read_csv(args.aggregate / "formal_top100_leaderboard.csv")
    metric = pd.read_csv(args.aggregate / "official_top100_per_metric_leaderboard.csv")
    cells = pd.read_csv(args.aggregate / "dataset_panel_leaderboard.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 12,
        "axes.labelsize": 9,
        "figure.dpi": 150,
    })
    fig = plt.figure(figsize=(14, 7.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(1.25, 1.0), height_ratios=(1, 1))
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])

    ordered = board.sort_values("mean_operation_rank", ascending=False)
    colors = ["#7B2CBF" if name == "WitnessCell" else "#B8C2CC" for name in ordered.method]
    ax_a.barh(ordered.method, ordered.mean_operation_rank, color=colors, edgecolor="white")
    for y, value in enumerate(ordered.mean_operation_rank):
        ax_a.text(value + 0.12, y, f"{value:.2f}", va="center", fontsize=8)
    ax_a.set_xlim(0, max(ordered.mean_operation_rank) + 1.25)
    ax_a.set_xlabel("Mean per-operation rank (lower is better)")
    ax_a.set_title("A  Official top-100 six-metric leaderboard", loc="left", fontweight="bold")
    ax_a.grid(axis="x", alpha=0.2)
    ax_a.spines[["top", "right", "left"]].set_visible(False)
    ax_a.tick_params(axis="y", length=0)

    chosen = ["WitnessCell", "scouter", "linearModel"]
    metric_order = ["PCC-delta", "MSE", "E-distance", "symmetric KL", "Wasserstein", "Common-DEGs"]
    pivot = metric[metric.method.isin(chosen)].pivot(index="metric", columns="method", values="position").loc[metric_order]
    x = np.arange(len(metric_order))
    width = 0.24
    palette = {"WitnessCell": "#7B2CBF", "scouter": "#3A506B", "linearModel": "#A9B4C2"}
    for index, method in enumerate(chosen):
        ax_b.bar(x + (index - 1) * width, pivot[method], width, label=method, color=palette[method])
    ax_b.set_xticks(x, ["PCC", "MSE", "E-dist", "sym-KL", "Wass.", "DEG"], rotation=0)
    ax_b.set_ylabel("Metric position (lower is better)")
    ax_b.set_ylim(0, 16)
    ax_b.invert_yaxis()
    ax_b.set_title("B  Different routes to the composite rank", loc="left", fontweight="bold")
    ax_b.legend(frameon=False, ncol=3, fontsize=8, loc="lower center")
    ax_b.grid(axis="y", alpha=0.2)
    ax_b.spines[["top", "right"]].set_visible(False)

    wc = cells[cells.method == "WitnessCell"].copy()
    panel_names = ["Norman", "Wessels", "Replogle_exp6", "Schmidt"]
    labels = ["Norman", "Wessels", "Replogle", "Schmidt"]
    top = [float(wc[(wc.DataSet == name) & (wc.DEG == 100)].dataset_panel_position.iloc[0]) for name in panel_names]
    allg = [float(wc[(wc.DataSet == name) & (wc.DEG == 5000)].dataset_panel_position.iloc[0]) for name in panel_names]
    x = np.arange(len(labels))
    ax_c.bar(x - 0.17, top, 0.34, label="Top 100 / six metrics", color="#7B2CBF")
    ax_c.bar(x + 0.17, allg, 0.34, label="Up to 5,000 / four metrics", color="#C8A2E8")
    ax_c.set_xticks(x, labels)
    ax_c.set_ylabel("Dataset-panel position")
    ax_c.set_ylim(0, 10)
    ax_c.invert_yaxis()
    ax_c.set_title("C  Cross-dataset sensitivity", loc="left", fontweight="bold")
    ax_c.legend(frameon=False, fontsize=8, loc="lower left")
    ax_c.grid(axis="y", alpha=0.2)
    ax_c.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "WitnessCell leads the official scPerturBench genetic-combination point ranking",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        -0.01,
        "15 published eligible methods + WitnessCell; four datasets; three split seeds; n=654 operation-seed units. "
        "Panel C shows that the aggregate lead is not uniform across environments.",
        ha="center",
        fontsize=8,
        color="#4B5563",
    )
    fig.savefig(args.out / "figure_main_scperturbench.png", bbox_inches="tight", dpi=300)
    fig.savefig(args.out / "figure_main_scperturbench.pdf", bbox_inches="tight")
    print(args.out / "figure_main_scperturbench.png")


if __name__ == "__main__":
    main()
