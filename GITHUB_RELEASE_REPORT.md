# GitHub compact-release report

## Outcome

Status: **GO for anonymous GitHub publication**, subject to the hosted Actions
run succeeding after upload.

Repository version: `18.0.0-github.1`; embedded Python package version:
`witnesscell 0.1.0`.

The compact repository is derived from the authoritative `grn_repo_v18.zip`
whose SHA-256 is
`fbac03b0101992a99093c802fa59d35a22750d89af46720b945ee6d52bbed64c`.
The full archive is not modified and is not represented as the Git repository.

## Engineering changes

- integrated the hardened `witnesscell 0.1.0` API, CLI and package tests;
- added a two-environment dependency lock for v14 prediction and metric stress;
- added reviewer smoke, frozen-result, anonymity and Git-size audits;
- added a one-command four-dataset × three-seed v14 reproduction runner;
- added canonical semantic NPZ digests independent of ZIP timestamps;
- added pinned, SHA-addressed GitHub Actions workflows;
- added English/Chinese navigation, data policy and upload checklist;
- removed generated caches, bytecode, historical nested archives, raw public
  matrices, duplicate predictions and machine-specific paths.

## Scientific preservation

- Gate 19: 12 official splits, 8 active corrections, 4 exact fallbacks and 654
  condition-seed units;
- Gate 15: 6,540 official metric rows and the exact primary aggregation;
- Gate 21: 213 test rows over 33 calibration-unseen target pairs;
- Gate 22: formal result tables and the original
  `NOT_ADJUDICATED` boundaries;
- exact and estimated Witness Risk verification code and compact results.

## Publication boundary

The compact repository is a new release identity, not a byte-preserving copy of
the full evidence archive. Excluded large assets are represented by hashes and
regeneration instructions. No author, affiliation, personal repository, email,
ORCID, project DOI or preregistration-account identifier is included.
