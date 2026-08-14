# Data and asset policy

## Included derived assets

The repository includes four condition-moment archives under
`reproduction/assets/condition_moments/`. They contain condition labels, gene
order, means, variances and cell counts required for CPU regeneration of the
v14 mean-response predictions. They do not contain individual cells.

The required GEARS gene universe and GO mapping are under
`reproduction/gate15_scperturbench/data/gears_assets/`. Frozen result tables,
not raw model checkpoints, are retained for the benchmark and named-model
audits.

All included asset hashes are recorded in `RELEASE_MANIFEST.sha256` and the
specialized manifests under `reproduction/assets/checksums/`.

## Public cell-level data excluded from Git

| Dataset | Compressed bytes | MD5 | Source |
|---|---:|---|---|
| Norman | 494,401,663 | `01d79ecd504186b574e124e7819e67bb` | Zenodo record 14638780 |
| Wessels | 97,345,484 | `8429ad50fafa82dc9d6d2de482e8f0b2` | Zenodo record 14638780 |
| Schmidt | 254,709,272 | `a0740b30e263a7a353d1da91906cb709` | Zenodo record 14638780 |
| Replogle_exp6 | 121,717,611 | `902c42c1b08d50fa8ee2ece79fd8f24c` | Zenodo record 14638780 |

`scripts/download_full_data.sh` delegates to the frozen downloader, resumes
partial transfers, verifies every MD5 and validates the decompressed HDF5
structure.

## Deliberately excluded

- raw `.h5ad` and compressed cell matrices;
- historical v12/v13/v14 duplicate deployment archives;
- prediction archives containing held-out truth;
- nested source and result archives;
- model checkpoints, package caches and compiled bytecode;
- machine-specific logs and absolute local paths;
- preregistration account metadata that could break double-blind review.

The original full v18 archive remains scientifically frozen by its SHA-256 in
`reproduction/assets/checksums/ORIGINAL_V18.sha256`.
