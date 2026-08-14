#!/usr/bin/env python3
"""Audit locally reproduced test conditions against published result rows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import pandas as pd

from official_split import (
    filter_gears_supported_conditions,
    load_gears_supported_genes,
    make_official_split,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--published", type=Path, required=True)
    parser.add_argument("--gears-assets", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.data, backed="r")
    raw_conditions = list(pd.unique(adata.obs["perturbation"].astype(str)))
    supported_genes = load_gears_supported_genes(args.gears_assets)
    conditions = filter_gears_supported_conditions(raw_conditions, supported_genes)
    published = pd.read_csv(args.published)
    rows = []
    for seed in (1, 2, 3):
        split = make_official_split(conditions, seed)
        expected = set(
            published[
                (published.DataSet == args.dataset)
                & (published.op.str.endswith(f"_{seed}"))
            ].op.str.rsplit("_", n=1).str[0]
        )
        observed = set(split.test)
        rows.append({
            "seed": seed,
            "train": len(split.train),
            "validation": len(split.validation),
            "test": len(split.test),
            "expected_test": len(expected),
            "exact_match": observed == expected,
            "missing": sorted(expected - observed),
            "unexpected": sorted(observed - expected),
        })
    report = {
        "status": "PASS" if all(row["exact_match"] for row in rows) else "FAIL",
        "dataset": args.dataset,
        "raw_conditions": len(raw_conditions),
        "official_supported_conditions": len(conditions),
        "gears_filtered_conditions": sorted(set(raw_conditions) - set(conditions)),
        "rows": rows,
    }
    (args.out / "split_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
