#!/usr/bin/env python3
"""Render the frozen confirmation figure from saved results only.

No model fitting, RNA access, resampling or additional statistical analysis.
Input: run1/{RESULTS.csv,CONDITIONAL_STABILITY.json,RUN_TARGET_RESULTS.csv,COMPLETE.json}
Output: PNG + PDF. Paths can be overridden for a portable results package.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "witnesscell_matplotlib_cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parent = Path(__file__).resolve().parents[1]
    parser.add_argument("--results-dir", type=Path, default=parent / "run1")
    parser.add_argument("--output-dir", type=Path, default=parent / "deliverables")
    args = parser.parse_args()
    src = args.results_dir
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(src / "RESULTS.csv").set_index("method")
    runs = pd.read_csv(src / "RUN_TARGET_RESULTS.csv")
    stability = json.loads((src / "CONDITIONAL_STABILITY.json").read_text())
    verdict = json.loads((src / "COMPLETE.json").read_text())
    assert verdict["decision"] == "NOT_CONFIRMED"
    assert int(summary.loc["system_fixed_reference", "n_pairs"]) == 374

    methods = ["additive", "scalar_saturation", "endpoint_control",
               "fixed_reference", "system_fixed_reference"]
    labels = ["Additive", "Scalar saturation", "Endpoint control",
              "Saturation + G", "Endpoint control + G"]
    ratios = summary.loc[methods, "full_sse"].to_numpy() / summary.loc["endpoint_control", "full_sse"]
    observed = np.array([
        summary.loc["system_fixed_reference", "full_gain_vs_endpoint"],
        summary.loc["system_fixed_reference", "P_gain_vs_endpoint"],
    ]) * 100
    ranges = np.array([[stability["results"][key]["q025"],
                        stability["results"][key]["q975"]]
                       for key in ["full_vs_endpoint_control", "P_vs_endpoint_control"]]) * 100
    run_rows = runs[runs.method == "system_fixed_reference"].set_index("run").loc[["rep1", "rep2", "rep3"]]
    run_gain = run_rows.gain_vs_endpoint.to_numpy() * 100

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 10,
        "axes.edgecolor": "#B1BAC2",
        "axes.linewidth": 0.8,
        "text.color": "#202D36",
        "axes.labelcolor": "#202D36",
        "xtick.color": "#52616C",
        "ytick.color": "#202D36",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    ink = "#202D36"
    secondary = "#5A6A76"
    accent = "#1C6474"
    pale = "#8299A4"
    neutral = "#8A979F"
    fig = plt.figure(figsize=(14.8, 6.0), facecolor="white")
    gs = fig.add_gridspec(1, 3, left=0.148, right=0.972, top=0.73,
                          bottom=0.30, width_ratios=[1.04, 1.06, 0.96], wspace=0.58)
    axa, axb, axc = [fig.add_subplot(gs[0, i]) for i in range(3)]

    fig.text(0.04, 0.925, "Frozen RNA program confirmation", size=18, weight="bold", color=ink)
    fig.text(0.04, 0.873,
             "Combi-seq  |  374 drug pairs  |  36 drugs  |  3 experimental runs  |  51 G2M genes",
             size=10.5, color=secondary)
    fig.text(0.967, 0.935, "NOT CONFIRMED", ha="right", va="center", size=12,
             weight="bold", color="#4A555D",
             bbox=dict(boxstyle="round,pad=0.5,rounding_size=0.2", facecolor="#EDF0F2", edgecolor="#CDD5DA"))

    def style(ax, letter, title):
        for side in ["top", "right", "left"]:
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=8)
        ax.tick_params(axis="x", length=3.5)
        ax.set_title(title, loc="left", weight="bold", pad=23)
        ax.text(-0.09, 1.087, letter, transform=ax.transAxes, size=15,
                weight="bold", va="center", ha="right")
        ax.grid(axis="x", color="#E5E9EC", linewidth=0.65, zorder=0)
        ax.set_axisbelow(True)

    style(axa, "A", "Full program error")
    ypos = np.arange(5)[::-1]
    axa.axvline(1.0, color="#697781", linewidth=1, linestyle=(0, (3, 3)), zorder=1)
    for i, (y, ratio) in enumerate(zip(ypos, ratios)):
        col = accent if methods[i] == "system_fixed_reference" else neutral
        if methods[i] == "endpoint_control":
            col = "#465D6B"
        axa.plot([1, ratio], [y, y], color=col, linewidth=2.6, solid_capstyle="round", zorder=2)
        axa.scatter(ratio, y, s=52 if methods[i] == "system_fixed_reference" else 37,
                    color=col, edgecolor="white", linewidth=0.6, zorder=3)
        axa.annotate(f"{ratio:.5f}", xy=(ratio, y), xytext=(7, 0),
                     textcoords="offset points", va="center", ha="left", size=9.4, color=col)
    axa.set_yticks(ypos, labels)
    axa.get_yticklabels()[-1].set_weight("bold")
    axa.set_ylim(-0.62, 4.62)
    axa.set_xlim(0.998, 1.041)
    axa.set_xticks([1.00, 1.01, 1.02, 1.03, 1.04])
    axa.set_xlabel("Relative SSE (endpoint control = 1)", labelpad=12)
    axa.text(0, -0.26, "Lower is better", transform=axa.transAxes, color=secondary, size=9)

    style(axb, "B", "Increment beyond endpoint control")
    by = np.array([1, 0])
    axb.axvline(0, color="#697781", linewidth=1, linestyle=(0, (3, 3)), zorder=1)
    for i, y in enumerate(by):
        lo, hi = ranges[i]
        axb.plot([lo, hi], [y, y], color=pale, linewidth=2.2, zorder=2)
        axb.plot([lo, lo], [y-0.045, y+0.045], color=pale, linewidth=1.3, zorder=2)
        axb.plot([hi, hi], [y-0.045, y+0.045], color=pale, linewidth=1.3, zorder=2)
        axb.scatter(observed[i], y, s=58, color=accent, edgecolor="white", linewidth=0.7, zorder=3)
        axb.text(observed[i], y+0.16, f"{observed[i]:+.4f}%", ha="center", color=accent, weight="bold", size=10)
        axb.text((lo+hi)/2, y-0.18, f"[{lo:+.3f}, {hi:+.3f}]%", ha="center", color=secondary, size=9)
    axb.set_yticks(by, ["Full program", "Endpoint-\northogonal (P)"])
    axb.set_ylim(-0.65, 1.65)
    axb.set_xlim(-1.22, 0.83)
    axb.set_xticks([-1.0, -0.5, 0.0, 0.5])
    axb.set_xlabel("Gain vs endpoint control (%)", labelpad=12)
    axb.text(0, -0.26, "Positive values favor the partner increment", transform=axb.transAxes,
             color=secondary, size=9)

    style(axc, "C", "Across experimental runs")
    cy = np.array([2, 1, 0])
    axc.axvline(0, color="#697781", linewidth=1, linestyle=(0, (3, 3)), zorder=1)
    for y, value in zip(cy, run_gain):
        axc.plot([0, value], [y, y], color=pale, linewidth=3, solid_capstyle="round", zorder=2)
        axc.scatter(value, y, s=53, color=accent, edgecolor="white", linewidth=0.7, zorder=3)
        axc.annotate(f"{value:+.3f}%", xy=(value, y), xytext=(0, 15),
                     textcoords="offset points", ha="center", va="center", size=10, color=accent)
    axc.set_yticks(cy, ["Run 1", "Run 2", "Run 3"])
    axc.set_ylim(-0.57, 2.6)
    axc.set_xlim(-0.36, 0.57)
    axc.set_xticks([-0.25, 0, 0.25, 0.5])
    axc.set_xlabel("Gain vs endpoint control (%)", labelpad=12)
    axc.text(0, -0.26, "Full program; run-specific targets", transform=axc.transAxes,
             color=secondary, size=9)

    fig.text(0.04, 0.116,
             "B: dots are observed gains. Lines are 2.5–97.5% conditional endpoint-multiplier percentile ranges (4,000 draws).",
             fontsize=9.7, color=secondary)
    fig.text(0.04, 0.079,
             "Predictions remain fixed in these draws; the ranges are not unconditional biological confidence intervals. "
             "P denotes endpoint-orthogonal error.",
             fontsize=9.7, color=secondary)
    fig.text(0.04, 0.034,
             "Independent source; drug-combination RNA extension. Prespecified confirmation criteria were not met.",
             fontsize=9.7, color=ink)

    base = out / "WitnessCell_Frozen_Program_Confirmation_20261020_Figure"
    fig.savefig(base.with_suffix(".png"), dpi=180, facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), facecolor="white",
                metadata={"Title": "WitnessCell frozen RNA program confirmation",
                          "Subject": "Exact rendering of saved confirmation results; no new analysis",
                          "Author": "WitnessCell analysis"})
    plt.close(fig)
    manifest = {name: hashlib.sha256((src / name).read_bytes()).hexdigest()
                for name in ["RESULTS.csv", "CONDITIONAL_STABILITY.json", "RUN_TARGET_RESULTS.csv", "COMPLETE.json"]}
    print(json.dumps({"input_sha256": manifest,
                      "outputs": {str(base.with_suffix(ext)): base.with_suffix(ext).stat().st_size
                                  for ext in [".png", ".pdf"]}}, indent=2))


if __name__ == "__main__":
    main()
