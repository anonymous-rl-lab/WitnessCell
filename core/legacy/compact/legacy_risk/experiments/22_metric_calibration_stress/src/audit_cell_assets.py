#!/usr/bin/env python3
"""Checksum and schema audit for the four formal cell-level datasets."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad


EXPECTED = {
    "Norman": {
        "compressed_bytes": 494401663,
        "compressed_md5": "01d79ecd504186b574e124e7819e67bb",
        "h5ad_sha256": "72748687e175983d7456d825334f1a553650800c2510dc6dc36e20941d0b8b63",
        "shape": [96994, 5025],
        "conditions": 227,
    },
    "Wessels": {
        "compressed_bytes": 97345484,
        "compressed_md5": "8429ad50fafa82dc9d6d2de482e8f0b2",
        "h5ad_sha256": "02f35f93c250af2b718374c4dfb2fa77f3559b1eefbb7066b26cf8e6a9e02fa5",
        "shape": [19775, 5020],
        "conditions": 128,
    },
    "Schmidt": {
        "compressed_bytes": 254709272,
        "compressed_md5": "a0740b30e263a7a353d1da91906cb709",
        "h5ad_sha256": "7797165c9982bf40dc2d7d6b530d29e7a1f2796b8c96a5367baaaee93433cf95",
        "shape": [60572, 5005],
        "conditions": 148,
    },
    "Replogle_exp6": {
        "compressed_bytes": 121717611,
        "compressed_md5": "902c42c1b08d50fa8ee2ece79fd8f24c",
        "h5ad_sha256": "87762c9e34bc36f5ab9e6e31d89d162fc35a528fb86aa1d051f0a55479bf3949",
        "shape": [27104, 5019],
        "conditions": 70,
    },
}


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for dataset, expected in EXPECTED.items():
        compressed = args.data_dir / f"{dataset}.h5ad.gz"
        raw = args.data_dir / f"{dataset}.h5ad"
        row = {"dataset": dataset, "compressed": str(compressed), "h5ad": str(raw)}
        row["compressed_bytes"] = compressed.stat().st_size if compressed.is_file() else None
        row["compressed_md5"] = digest(compressed, "md5") if compressed.is_file() else None
        row["h5ad_sha256"] = digest(raw, "sha256") if raw.is_file() else None
        if raw.is_file():
            adata = ad.read_h5ad(raw, backed="r")
            row["shape"] = [int(adata.n_obs), int(adata.n_vars)]
            row["conditions"] = int(adata.obs["perturbation"].astype(str).nunique())
            row["condition_key"] = "perturbation"
            row["genes_unique"] = bool(adata.var_names.is_unique)
            adata.file.close()
        else:
            row.update({"shape": None, "conditions": None, "condition_key": None, "genes_unique": False})
        row["checks"] = {
            "compressed_bytes": row["compressed_bytes"] == expected["compressed_bytes"],
            "compressed_md5": row["compressed_md5"] == expected["compressed_md5"],
            "h5ad_sha256": row["h5ad_sha256"] == expected["h5ad_sha256"],
            "shape": row["shape"] == expected["shape"],
            "conditions": row["conditions"] == expected["conditions"],
            "genes_unique": row["genes_unique"],
        }
        row["pass"] = all(row["checks"].values())
        rows.append(row)
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": "scPerturBench Zenodo record 14638780 via the frozen Gate 15 data manifest",
        "datasets": rows,
        "pass": all(row["pass"] for row in rows),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print("CELL_ASSET_AUDIT_PASS" if report["pass"] else "CELL_ASSET_AUDIT_FAIL")
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
