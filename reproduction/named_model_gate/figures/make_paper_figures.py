#!/usr/bin/env python3
"""Generate the frozen paper figure and its diagnostic companion."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


ORDER = ["Saturation", "GEARS", "CPA", "GEARS+Witness", "CPA+Witness", "WitnessCell"]
COLORS = {
    "Saturation": "#8A8A8A",
    "GEARS": "#2C73B9",
    "CPA": "#D88920",
    "GEARS+Witness": "#16857B",
    "CPA+Witness": "#7651A6",
    "WitnessCell": "#C33C54",
}


def style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
    })


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.12, 1.08, label, transform=axis.transAxes, fontsize=11, fontweight="bold", va="top")


def box(axis: plt.Axes, xy: tuple[float, float], width: float, height: float,
        text: str, face: str, edge: str = "#374151", fontsize: float = 6.4) -> None:
    patch = FancyBboxPatch(
        xy, width, height, boxstyle="round,pad=0.015,rounding_size=0.025",
        linewidth=0.8, edgecolor=edge, facecolor=face,
    )
    axis.add_patch(patch)
    axis.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=fontsize)


def arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=8,
                                   linewidth=0.8, color="#4B5563"))


def design_panel(axis: plt.Axes) -> None:
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    box(axis, (0.02, 0.63), 0.24, 0.21, "Source half A\n12 train pairs\n5 validation pairs", "#E8F1FA", fontsize=5.8)
    box(axis, (0.02, 0.20), 0.24, 0.21, "Destination\nhalf B\n5 calibration | 14 test", "#F8E8EC", fontsize=5.6)
    box(axis, (0.37, 0.65), 0.24, 0.17, "Named AIVCs\nGEARS | CPA", "#F3F4F6", fontsize=6.0)
    box(axis, (0.37, 0.20), 0.24, 0.21, "Witness layer\ngeometry\n+ calibration", "#E7F5F2", fontsize=5.8)
    box(axis, (0.72, 0.37), 0.26, 0.29, "Frozen score\nsix strategies\n14 strict-seen1\ntarget pairs\n3 splits x 3 seeds", "#FFF5D8", fontsize=5.5)
    arrow(axis, (0.26, 0.735), (0.37, 0.735))
    arrow(axis, (0.26, 0.305), (0.37, 0.305))
    arrow(axis, (0.61, 0.735), (0.72, 0.57))
    arrow(axis, (0.61, 0.305), (0.72, 0.46))
    axis.text(0.02, 0.95, "Target outcomes are scoring-only", color="#9B1C31", fontsize=7.2, fontweight="bold")
    axis.text(0.50, 0.04, "18 named-model fits; inference over 14 biological pairs", ha="center", fontsize=6.8, color="#4B5563")


def residual_panel(axis: plt.Axes, strategy: pd.DataFrame) -> None:
    local = strategy.set_index("strategy").loc[ORDER].reset_index()
    y = np.arange(len(local))[::-1]
    for position, row in zip(y, local.itertuples(index=False)):
        low = row.residual_mse_ci95_low
        high = row.residual_mse_ci95_high
        axis.errorbar(row.residual_mse, position,
                      xerr=[[row.residual_mse - low], [high - row.residual_mse]],
                      fmt="o", color=COLORS[row.strategy], ecolor=COLORS[row.strategy],
                      markersize=5.2, elinewidth=1.3, capsize=2.5, zorder=3)
        axis.text(row.residual_mse * 1.12, position, f"{row.residual_mse:.3f}",
                  va="center", fontsize=6.5, color="#252525")
    axis.set_yticks(y, local.strategy)
    axis.set_xscale("log")
    axis.set_xlim(0.008, 0.27)
    axis.set_xlabel("Interaction-residual MSE (log scale; lower is better)")
    axis.grid(axis="x", which="major", color="#E5E7EB", linewidth=0.7)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    axis.set_title("Absolute predictive error", loc="left", fontweight="bold")


def paired_panel(axis: plt.Axes, targets: pd.DataFrame) -> None:
    pair = targets.groupby(["pair", "strategy"], as_index=False).residual_mse.mean()
    wide = pair.pivot(index="pair", columns="strategy", values="residual_mse")
    x = wide["Saturation"].to_numpy(float)
    y = wide["WitnessCell"].to_numpy(float)
    lo = min(x.min(), y.min()) * 0.92
    hi = max(x.max(), y.max()) * 1.06
    axis.plot([lo, hi], [lo, hi], linestyle="--", color="#9CA3AF", linewidth=1)
    axis.scatter(x, y, s=31, color=COLORS["WitnessCell"], edgecolor="white", linewidth=0.55, zorder=3)
    axis.set_xlim(lo, hi)
    axis.set_ylim(lo, hi)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Saturation residual MSE")
    axis.set_ylabel("WitnessCell residual MSE")
    axis.grid(color="#EEF0F2", linewidth=0.6)
    axis.spines[["top", "right"]].set_visible(False)
    axis.set_title("Every held-out target improves", loc="left", fontweight="bold")
    axis.text(0.04, 0.94, "14/14 pairs\n42.49% mean reduction\np = 6.1 x 10$^{-5}$",
              transform=axis.transAxes, va="top", fontsize=7,
              bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#D1D5DB"})


def forest_panel(axis: plt.Axes, frozen: pd.DataFrame, conservative: pd.DataFrame) -> None:
    rows = pd.concat([frozen, conservative], ignore_index=True)
    keys = [
        ("WitnessCell", "Saturation", "WitnessCell vs saturation"),
        ("WitnessCell", "GEARS", "WitnessCell vs GEARS"),
        ("WitnessCell", "CPA", "WitnessCell vs CPA"),
        ("GEARS+Witness", "GEARS", "GEARS+Witness vs GEARS"),
        ("CPA+Witness", "CPA", "CPA+Witness vs CPA"),
        ("GEARS+Witness", "WitnessCell", "GEARS+Witness vs WitnessCell"),
        ("CPA+Witness", "WitnessCell", "CPA+Witness vs WitnessCell"),
    ]
    selected = []
    for candidate, comparator, label in keys:
        local = rows[(rows.candidate == candidate) & (rows.comparator == comparator)]
        assert len(local) == 1
        row = local.iloc[0].to_dict()
        row["label"] = label
        selected.append(row)
    local = pd.DataFrame(selected)
    y = np.arange(len(local))[::-1]
    estimates = 100 * local.relative_residual_mse_improvement.to_numpy(float)
    low = 100 * local.pair_bootstrap_ci95_low.to_numpy(float)
    high = 100 * local.pair_bootstrap_ci95_high.to_numpy(float)
    colors = [COLORS[row.candidate] for row in local.itertuples(index=False)]
    for index, (yy, estimate, lower, upper, color) in enumerate(zip(y, estimates, low, high, colors)):
        filled = not (lower <= 0 <= upper)
        axis.errorbar(estimate, yy, xerr=[[estimate - lower], [upper - estimate]],
                      fmt="o", markersize=5, color=color, ecolor=color,
                      markerfacecolor=color if filled else "white", markeredgewidth=1.1,
                      elinewidth=1.2, capsize=2.3)
        axis.text(min(97.5, upper + 2.2), yy, f"{estimate:.1f}%", va="center", fontsize=6.4)
    axis.axvline(0, linestyle="--", linewidth=0.8, color="#6B7280")
    axis.axhline(3.5, color="#D1D5DB", linewidth=0.7)
    axis.axhline(1.5, color="#D1D5DB", linewidth=0.7)
    axis.set_yticks(y, local.label)
    axis.set_xlim(-4, 101)
    axis.set_xlabel("Residual-MSE reduction (%, pair bootstrap 95% CI)")
    axis.grid(axis="x", color="#EEF0F2", linewidth=0.6)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="y", length=0)
    axis.set_title("Effect sizes and model complementarity", loc="left", fontweight="bold")
    axis.text(0.985, 0.02, "Open point: CI crosses zero", transform=axis.transAxes,
              ha="right", va="bottom", fontsize=6.2, color="#6B7280")


def make_main(audit: Path, output: Path) -> None:
    strategy = pd.read_csv(audit / "strategy_metrics.csv")
    targets = pd.read_csv(audit / "independent_target_rows.csv")
    frozen = pd.read_csv(audit / "frozen_comparisons.csv")
    conservative = pd.read_csv(audit / "conservative_comparisons.csv")
    fig = plt.figure(figsize=(7.2, 7.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.92, 1.08], width_ratios=[1.03, 0.97])
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    design_panel(axes[0])
    residual_panel(axes[1], strategy)
    paired_panel(axes[2], targets)
    forest_panel(axes[3], frozen, conservative)
    for axis, label in zip(axes, "ABCD"):
        panel_label(axis, label)
    fig.suptitle("Target-conditioned witness geometry improves strict-seen1 combinatorial prediction",
                 fontsize=10.5, fontweight="bold")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output.with_suffix(f".{suffix}"), bbox_inches="tight")
    plt.close(fig)


def make_diagnostic(audit: Path, output: Path) -> None:
    strategy = pd.read_csv(audit / "strategy_metrics.csv").set_index("strategy").loc[ORDER].reset_index()
    metadata = pd.read_csv(audit / "run_metadata.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    y = np.arange(len(strategy))[::-1]
    for position, row in zip(y, strategy.itertuples(index=False)):
        axes[0].errorbar(row.full_effect_cosine, position,
                         xerr=[[row.full_effect_cosine - row.full_effect_cosine_ci95_low],
                               [row.full_effect_cosine_ci95_high - row.full_effect_cosine]],
                         fmt="o", color=COLORS[row.strategy], markersize=5, capsize=2.5)
    axes[0].set_yticks(y, strategy.strategy)
    axes[0].set_xlabel("Full-effect cosine (pair bootstrap 95% CI)")
    axes[0].set_title("Expression-effect fidelity", loc="left", fontweight="bold")
    axes[0].grid(axis="x", color="#EEF0F2")
    axes[0].spines[["top", "right", "left"]].set_visible(False)
    axes[0].tick_params(axis="y", length=0)

    values = [metadata.lambda_gears, metadata.lambda_cpa]
    labels = ["GEARS fusion", "CPA fusion"]
    colors = [COLORS["GEARS+Witness"], COLORS["CPA+Witness"]]
    for index, (series, label, color) in enumerate(zip(values, labels, colors)):
        axes[1].scatter(np.full(len(series), index) + np.linspace(-0.09, 0.09, len(series)), series,
                        s=22, color=color, alpha=0.88)
        axes[1].plot([index - 0.16, index + 0.16], [series.mean(), series.mean()], color="#111827", linewidth=1.3)
    axes[1].set_xticks([0, 1], labels)
    axes[1].set_ylabel("Witness fusion weight $\\lambda$")
    axes[1].set_ylim(0.88, 1.005)
    axes[1].set_title("Calibration chooses mostly Witness", loc="left", fontweight="bold")
    axes[1].grid(axis="y", color="#EEF0F2")
    axes[1].spines[["top", "right"]].set_visible(False)
    for axis, label in zip(axes, "AB"):
        panel_label(axis, label)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output.with_suffix(f".{suffix}"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, default=Path(__file__).resolve().parents[1] / "independent_audit")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    style()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    make_main(args.audit, args.out_dir / "Fig_main_named_AIVC_gate")
    make_diagnostic(args.audit, args.out_dir / "Fig_extended_named_AIVC_diagnostics")
    print(f"paper figures written to {args.out_dir}")


if __name__ == "__main__":
    main()
