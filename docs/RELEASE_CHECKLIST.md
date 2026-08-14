# Release checklist

## Contract and tests

- [ ] Confirm version in `pyproject.toml`, `_version.py`, changelog, and CFF.
- [ ] Run unit, integration, CLI, serialization, tamper, and leakage tests.
- [ ] Run the 12-split reference-parity audit against the trusted frozen repo.
- [ ] Confirm gate activations, selected hyperparameters, factorized predictions,
      and WitnessCell predictions match the reference tolerances.
- [ ] Run lint, type checks, coverage, and dependency audit.

## Double-blind anonymity

- [ ] Run `python scripts/audit_anonymity.py .`.
- [ ] Run the same audit over unpacked wheel, sdist, and source ZIP.
- [ ] Inspect METADATA, RECORD, PKG-INFO, CFF, README, license notice and ZIP paths.
- [ ] Confirm no name, affiliation, email, repository user, DOI, ORCID, grant,
      machine path, VCS remote, or original research/manuscript asset is present.
- [ ] Confirm neutral attribution is `WitnessCell Authors` only.

## Distribution

- [ ] Build twice with a fixed `SOURCE_DATE_EPOCH`, normalize both sdists, and
      require byte-identical wheel and sdist artifacts.
- [ ] Run `twine check --strict dist/*`.
- [ ] Install the wheel into a fresh virtual environment and run smoke tests.
- [ ] Install the sdist into another fresh environment and run smoke tests.
- [ ] Check wheel contents and verify that tests/docs do not leak target data.
- [ ] Generate `SHA256SUMS` and provenance/release report.

## PyPI

- [ ] Reserve or verify the `witnesscell` project name without exposing identity.
- [ ] Configure PyPI Trusted Publishing for the release workflow and protected
      `pypi` environment; require manual approval.
- [ ] Publish to TestPyPI first and verify installation on all supported Python
      versions.
- [ ] Create immutable `v0.1.0` tag; verify it with
      `scripts/verify_release_tag.py`, then publish to PyPI.
- [ ] Verify project classifiers, rendered README, files, hashes and attestations.

The included release workflow intentionally stops at environment approval and
OIDC publishing. It contains no long-lived PyPI token.
