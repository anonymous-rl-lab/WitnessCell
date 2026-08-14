# Experiment 22 public-data asset policy

The compact `grn_repo_v17` release does not duplicate the four public
scPerturBench cell matrices or their compressed copies. Together they occupy
approximately 7.9 GB and are already identified by the frozen public record.

The following release assets are retained instead:

- the source record (`scPerturBench` Zenodo record 14638780);
- compressed byte counts and MD5 values;
- decompressed h5ad SHA-256 values;
- matrix shapes, condition counts and condition keys;
- the four candidate-blind Phase-M contracts and split ledgers derived before
  prediction access;
- all formal per-operation tables, inference arrays and verdicts.

The authoritative public-data ledger is
`audit/cell_asset_ledger.json`. The scientific `FROZEN_MANIFEST.sha256` seals
that ledger and the derived contracts, but deliberately does not list the
multi-gigabyte matrices themselves. A full from-cell replay must place the
eight files at the exact paths recorded in the ledger, rerun the frozen
cell-asset audit, and obtain both a cell-ledger pass and a complete
`sha256sum -c FROZEN_MANIFEST.sha256` pass before execution. The compact-release
audit allows only the eight ledger-identified public files to remain external;
no protocol, code, contract, result or amendment may be absent.

This policy reduces transport size without replacing public data by synthetic
or summary-only evidence. The preserved contracts are sufficient to audit and
re-score the completed experiment, whereas a from-cell replay requires the
original public matrices.
