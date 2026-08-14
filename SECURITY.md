# Security policy

## Supported version

Security fixes are provided for the latest released minor version. Version
0.1.x is the current supported line.

## Reporting

Do not disclose suspected vulnerabilities publicly. During double-blind review,
use the submission system's confidential correspondence channel. After review,
use the private reporting channel published on the package project page.

Model bundles are pickle-free and integrity checked. The loader also enforces
outer-entry, nested-NPZ, expanded-byte and array-element limits and rejects
non-finite or semantically inconsistent fitted state. A successful integrity
check establishes internal bundle consistency, not publisher identity; verify
release checksums or PyPI attestations before trusting provenance.
