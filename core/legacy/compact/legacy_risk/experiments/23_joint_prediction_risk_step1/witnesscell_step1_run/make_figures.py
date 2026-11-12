#!/usr/bin/env python3
"""Plot completed Step 1 outputs; no fitting, evaluation, or resampling.

Run after evaluate.py and, when available, summarize_baselines.py:
  python make_figures.py --evaluation evaluation --results results --out figures

Reads evaluation.json, selection_manifest.csv, baseline_summary.json (optional),
and results/*/scored_queries.csv. All intervals come from the existing 5,000-draw
paired-cluster summaries. No expression matrices, model files, evaluation-only
arrays, or bootstrap generators are opened/imported. Missing baseline results
are labelled pending. Scouter is never assigned a numerical result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


PROTOCOL_ID = "WitnessCell-step1-v1-20261017"
DATASETS = ("Norman", "Wessels")
SCORERS = ("S", "B", "effect", "variance")
SCORE_COLORS = {"S": "#0072B2", "B": "#D55E00", "effect": "#009E73", "variance": "#CC79A7"}
DATASET_COLORS = {"Norman": "#0072B2", "Wessels": "#D55E00"}
COMPARATORS = ("trainMean", "shared_endpoint_additive", "saturation_only", "factorized_prediction")
BASELINE_LABELS = ("Completed additive", "Shared-endpoint additive", "Saturation only", "Factorized prediction")


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite(value):
    return value is not None and math.isfinite(float(value))


def unique_row(rows, **criteria):
    found = [r for r in rows if all(r.get(k) == v for k, v in criteria.items())]
    if len(found) > 1:
        raise ValueError(f"Duplicate summary rows: {criteria}")
    return found[0] if found else None


def result(report, dataset, loss="loss_mse", comparison="S/B", coverage=.5):
    return unique_row(report["datasets"][dataset]["results"], loss=loss,
                      comparison=comparison, coverage=coverage)


def load_inputs(evaluation, results):
    report = read_json(evaluation / "evaluation.json")
    if report.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Unexpected evaluation protocol")
    # Require complete formal summaries BEFORE opening scored outcome tables.
    for dataset in DATASETS:
        ds = report.get("datasets", {}).get(dataset, {})
        if sorted(ds.get("seeds", [])) != [1, 2, 3]:
            raise ValueError(f"Wait for all three evaluated seeds: {dataset}")
        if ds.get("bootstrap_replicates_requested") != 5000:
            raise ValueError(f"Expected existing 5,000-draw evaluation: {dataset}")
    inventory = {(r["dataset"], int(r["seed"])): r for r in report["input_inventory"]}
    scored = {}
    for dataset in DATASETS:
        for seed in (1, 2, 3):
            path = results / f"{dataset}_seed{seed}" / "scored_queries.csv"
            reference = inventory.get((dataset, seed))
            if reference is None or sha256(path) != reference["sha256"]:
                raise ValueError(f"Scored CSV does not match completed evaluation: {path}")
            rows = read_csv(path)
            for row in rows:
                key = (row["dataset"], int(row["seed"]), row["pair"])
                if key[:2] != (dataset, seed) or key in scored:
                    raise ValueError(f"Invalid or duplicate query identity: {key}")
                effect = float(row["true_effect"])
                if not math.isfinite(effect) or effect < 0:
                    raise ValueError(f"Invalid true effect: {key}")
                scored[key] = effect
    selected = {}
    for row in read_csv(evaluation / "selection_manifest.csv"):
        if float(row["coverage_requested"]) != .5 or row["scorer"] not in ("S", "B"):
            continue
        key = (row["dataset"], int(row["seed"]), row["pair"], row["scorer"])
        if key in selected or row["accepted"] not in ("0", "1"):
            raise ValueError(f"Invalid selection row: {key}")
        selected[key] = bool(int(row["accepted"]))
    expected = {(*key, scorer) for key in scored for scorer in ("S", "B")}
    if set(selected) != expected:
        raise ValueError("Selection manifest and scored query identities differ")
    for dataset in DATASETS:
        for seed in (1, 2, 3):
            for scorer in ("S", "B"):
                flags = [v for (d, s, _, m), v in selected.items()
                         if (d, s, m) == (dataset, seed, scorer)]
                if sum(flags) != len(flags) // 2:
                    raise ValueError("Selection manifest violates the frozen half quota")
    path = evaluation / "baseline_summary.json"
    baseline = read_json(path) if path.is_file() else None
    if baseline is not None and baseline.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Unexpected baseline protocol")
    return report, scored, selected, baseline


def clean_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=3, width=.7)
    ax.grid(axis="y", color="#E5E5E5", linewidth=.6)
    ax.set_axisbelow(True)


def panel_title(fig, spec, letter, title):
    bounds = spec.get_position(fig)
    fig.text(bounds.x0, bounds.y1 + .018, letter, weight="bold", fontsize=14)
    fig.text(bounds.x0 + .022, bounds.y1 + .018, title, fontsize=11, weight="bold")


def pending(ax, text):
    ax.text(.5, .5, text, transform=ax.transAxes, ha="center", va="center",
            color="#666666", fontsize=10, wrap=True)


def draw_interval(ax, y, row, color, kind="gain", scale=1., marker="o"):
    """Draw stored point/CI directly; never recompute a bootstrap interval."""
    if row is None or not finite(row.get(kind)):
        return False
    point = float(row[kind]) * scale
    ci = row.get(kind + "_ci95")
    if ci is not None and len(ci) == 2 and all(finite(v) for v in ci):
        lo, hi = (float(v) * scale for v in ci)
        if lo > hi:
            raise ValueError("Reversed confidence interval")
        # CI need not contain the point estimate; do not force asymmetric xerr.
        ax.plot([lo, hi], [y, y], color=color, linewidth=1.8, solid_capstyle="round")
        ax.plot([lo, lo], [y-.045, y+.045], color=color, linewidth=1.)
        ax.plot([hi, hi], [y-.045, y+.045], color=color, linewidth=1.)
    ax.plot(point, y, marker=marker, color=color, markersize=5.8, zorder=4)
    if row.get("valid_bootstrap_replicates", 5000) != 5000:
        ax.annotate("*", (point, y), xytext=(5, 4), textcoords="offset points",
                    fontsize=10, color=color)
    return True


def draw_a(fig, spec, report):
    axes = [fig.add_subplot(s) for s in spec.subgridspec(1, 2, wspace=.27)]
    for ax, dataset in zip(axes, DATASETS):
        ds = report["datasets"][dataset]
        for scorer in SCORERS:
            xs, ys = [], []
            for coverage in (.3, .5, .7, .9):
                row = result(report, dataset, comparison=scorer + "/random_expectation", coverage=coverage)
                seeds = [r for r in ds["seed_results"]
                         if r["loss"] == "loss_mse" and r["coverage_requested"] == coverage]
                if len(seeds) != 3 or {r["seed"] for r in seeds} != {1, 2, 3}:
                    raise ValueError(f"Incomplete risk-coverage summary: {dataset}, {coverage}")
                if row is None or not finite(row.get("ratio")):
                    raise ValueError(f"Undefined normalized MSE: {dataset}, {scorer}, {coverage}")
                xs.append(100 * np.mean([r["coverage_actual"] for r in seeds]))
                ys.append(row["ratio"])
            xs.append(100.); ys.append(1.)
            ax.plot(xs, ys, "o-", color=SCORE_COLORS[scorer], label=scorer,
                    linewidth=1.5, markersize=3.7)
        ax.axhline(1, color="#777777", linestyle="--", linewidth=.8)
        ax.set(title=dataset, xlabel="Accepted queries (%)", xlim=(25, 102))
        ax.set_xticks([30, 50, 70, 90, 100])
        ax.set_ylim(bottom=0)
        clean_axis(ax)
    axes[0].set_ylabel("Non-target MSE / full-query MSE")
    axes[0].legend(ncol=2, frameon=False, fontsize=8, loc="best")
    panel_title(fig, spec, "A", "Risk at fixed coverage")


def draw_b(fig, spec, report):
    ax = fig.add_subplot(spec)
    for y, dataset in zip([1, 0], DATASETS):
        row = result(report, dataset)
        if not draw_interval(ax, y, row, DATASET_COLORS[dataset], scale=100):
            ax.text(0, y, "Not estimable", va="center", color="#666666")
    ax.axvline(0, color="#555555", linewidth=.8)
    ax.axvline(10, color="#999999", linewidth=.8, linestyle="--")
    ax.set(yticks=[1, 0], yticklabels=DATASETS, ylim=(-.5, 1.5),
           xlabel="Non-target MSE reduction: 100 × (1 − S/B) (%)")
    ax.text(.99, .96, "50% requested coverage\nDashed: prespecified 10% threshold",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#555555")
    clean_axis(ax)
    panel_title(fig, spec, "B", "Increment over B")


def weighted_ecdf(scored, selected, dataset, scorer, accepted):
    values, weights = [], []
    for seed in (1, 2, 3):
        x = [value for key, value in scored.items()
             if key[:2] == (dataset, seed) and selected[(*key, scorer)] == accepted]
        if not x:
            raise ValueError("Empty distribution stratum; cannot give seeds equal weight")
        values.extend(x)
        weights.extend([1. / (3 * len(x))] * len(x))
    order = np.argsort(values, kind="stable")
    x = np.asarray(values)[order]
    y = np.cumsum(np.asarray(weights)[order])
    # Collapse tied x values so a mass is one ECDF step, rather than draw-order dependent.
    ends = np.r_[x[1:] != x[:-1], True]
    return x[ends], y[ends]


def draw_c(fig, spec, report, scored, selected):
    grid = spec.subgridspec(2, 2, height_ratios=[1.05, .9], hspace=.85, wspace=.27)
    for i, dataset in enumerate(DATASETS):
        ax = fig.add_subplot(grid[0, i])
        for scorer in ("S", "B"):
            for accept, style in ((True, "-"), (False, "--")):
                x, y = weighted_ecdf(scored, selected, dataset, scorer, accept)
                ax.step(np.r_[0., x], np.r_[0., y], where="post", linestyle=style,
                        color=SCORE_COLORS[scorer], linewidth=1.3,
                        label=f"{scorer} {'accepted' if accept else 'rejected'}")
        ax.set(title=dataset, xlabel="True non-target effect (mean square)", ylim=(0, 1.03), xlim=(0, None))
        ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        ax.set_yticks([0, .5, 1])
        clean_axis(ax)
        if i == 0:
            ax.set_ylabel("Seed-balanced ECDF")
            ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")
    ax = fig.add_subplot(grid[1, :])
    for i, dataset in enumerate(DATASETS):
        offset = .12 if i == 0 else -.12
        rows = [result(report, dataset, loss="loss_relative"),
                report["datasets"][dataset].get("within_effect_quartile_control")]
        for y, row in zip([1, 0], rows):
            if not draw_interval(ax, y+offset, row, DATASET_COLORS[dataset], scale=100):
                ax.text(.98, (y+offset+.5)/2, f"{dataset}: not estimable", transform=ax.transAxes,
                        ha="right", va="center", color=DATASET_COLORS[dataset], fontsize=7)
    ax.axvline(0, color="#777777", linewidth=.8)
    ax.set(yticks=[1, 0], yticklabels=["Relative error", "Effect strata"],
           ylim=(-.5, 1.5), xlabel="S/B reduction (%) · stored 95% CI")
    clean_axis(ax)
    panel_title(fig, spec, "C", "Effect-size sensitivity (descriptive)")


def draw_d(fig, spec, baseline, pool):
    ax = fig.add_subplot(spec)
    any_plotted = False
    for j, dataset in enumerate(DATASETS):
        ds = (baseline or {}).get("datasets", {}).get(dataset, {})
        ready = (sorted(ds.get("seeds", [])) == [1, 2, 3]
                 and ds.get("bootstrap_replicates") == 5000
                 and not (baseline or {}).get("diagnostic_only", False))
        for i, comparator in enumerate(COMPARATORS):
            row = unique_row(ds.get("results", []), pool=pool, metric="loss_wmse", comparator=comparator) if ready else None
            y = 3-i + (.12 if j == 0 else -.12)
            drawn = draw_interval(ax, y, row, DATASET_COLORS[dataset], kind="ratio")
            any_plotted |= drawn
            if not drawn:
                ax.text(.99, (y+.55)/4.1, f"{dataset}: pending" if not ready or row is None else f"{dataset}: not estimable",
                        transform=ax.transAxes, ha="right", va="center", fontsize=7,
                        color=DATASET_COLORS[dataset])
    ax.axvline(1, color="#777777", linewidth=.8, linestyle="--")
    ax.set(yticks=[3, 2, 1, 0], yticklabels=BASELINE_LABELS, ylim=(-.55, 3.55),
           xlabel="WitnessCell / comparator WMSE (lower is better)")
    if not any_plotted:
        ax.set_xlim(0, 2)
    else:
        ax.set_xlim(left=0)
    handles = [Line2D([], [], marker="o", color=DATASET_COLORS[d], linewidth=1.4, label=d) for d in DATASETS]
    ax.legend(handles=handles, frameon=False, ncol=2, fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.19))
    label = "All queries" if pool == "all_queries" else "S-accepted half; common queries"
    ax.text(.99, 1.04, label, transform=ax.transAxes, ha="right", fontsize=8, color="#555555")
    ax.text(0, -.38, "Scouter: pending — no locked prediction archive", transform=ax.transAxes,
            fontsize=8, color="#666666")
    clean_axis(ax)
    panel_title(fig, spec, "D", "Prediction controls")


def make_figure(report, scored, selected, baseline, out, pool):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.titlesize": 10, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.linewidth": .7})
    fig = plt.figure(figsize=(14.4, 10.9), facecolor="white")
    grid = fig.add_gridspec(2, 2, left=.065, right=.975, bottom=.15, top=.9,
                           hspace=.43, wspace=.40, height_ratios=[.9, 1.25])
    draw_a(fig, grid[0, 0], report)
    draw_b(fig, grid[0, 1], report)
    draw_c(fig, grid[1, 0], report, scored, selected)
    draw_d(fig, grid[1, 1], baseline, pool)
    fig.suptitle("WitnessCell Step 1 | Retrospective risk screen", x=.065, y=.975,
                 ha="left", fontsize=17, weight="bold")
    fig.text(.065, .944, "Frozen scores and predictions · Norman + Wessels · Three seeds per dataset",
             fontsize=10, color="#555555")
    fig.text(.065, .038,
             "A: geometric mean of seed-level selected/full-query MSE ratios; points use mean actual coverage. "
             "B–D: existing 95% paired cluster-bootstrap intervals (5,000 draws).\n"
             "C: accepted/rejected ECDFs weight seeds equally; true effect is descriptive. Relative error uses the frozen OOF floor; "
             "quartile control uses TRAIN-OOF predicted-effect cutpoints.\n"
             "Ratios weight seeds equally. * Fewer than 5,000 valid draws; inspect the stored summary. "
             "No new resampling, prospective validation, or complete SOTA claim.",
             fontsize=7.7, color="#555555", linespacing=1.5)
    out.mkdir(parents=True, exist_ok=True)
    stem = out / "witnesscell_step1_four_panel"
    for suffix in (".png", ".pdf"):
        path = stem.with_suffix(suffix)
        fig.savefig(path, dpi=300, facecolor="white", metadata={"Title": "WitnessCell Step 1: retrospective risk screen"})
        print(path.resolve())
    plt.close(fig)


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, default=here / "evaluation")
    parser.add_argument("--results", type=Path, default=here / "results")
    parser.add_argument("--out", type=Path, default=here / "figures")
    parser.add_argument("--baseline-pool", choices=("all_queries", "S_accept_50pct"), default="all_queries",
                        help="Fixed comparison pool for panel D (default: all queries)")
    args = parser.parse_args()
    report, scored, selected, baseline = load_inputs(args.evaluation, args.results)
    make_figure(report, scored, selected, baseline, args.out, args.baseline_pool)


if __name__ == "__main__":
    main()
