# WitnessCell 0.1.0 release report

## Release identity

- Distribution name: `witnesscell`
- Version: `0.1.0`
- License: Apache-2.0
- Python: 3.10–3.14
- Algorithm contract: `witnesscell-v14-frozen`
- Review identity: fully anonymous; neutral attribution `WitnessCell Authors`
- Publication state: locally release-qualified; not uploaded to TestPyPI/PyPI

## Frozen-reference parity

The product implementation was evaluated on the complete Gate 19 matrix:
Norman, Replogle_exp6, Schmidt and Wessels, each at seeds 1, 2 and 3.

- Splits passed: 12 / 12
- Gate states matched: 12 / 12
- Selected `alpha`, `noise_ratio`, and `gamma` matched: 12 / 12
- Maximum absolute WitnessCell prediction delta: `2.384185791015625e-07`
- Maximum absolute factorized prediction delta: `2.384185791015625e-07`
- Acceptance tolerance: `2e-5`

The residual delta is consistent with comparison against the frozen float32
deployment archives. The audit also reran the original reference class to read
the v13 baseline sparse state directly, avoiding a legacy manifest property
whose `sparse_active` label was mapped to the v14 amplitude gate.

Machine-readable evidence: `reports/reference_parity.json` in the complete
source bundle.

## Locally verified controls

- Python compilation over package, scripts and tests.
- Sixty-seven unit, integration, golden-contract, adversarial-input and security
  tests passed under real `pytest`.
- Branch-aware coverage: `89.91%` against a required `80%` threshold.
- `ruff` and strict `mypy` completed with zero findings.
- `pip-audit` found no known vulnerabilities in the resolved runtime dependency
  set; `twine check --strict` passed for wheel and sdist.
- Direct-constructor, NaN/Inf, fractional-count, split-overlap, malformed-risk,
  semantically forged model and resource-limit regressions are locked by tests.
- Wheel install and import from a clean virtual environment.
- Source-distribution build/install and prediction smoke test from a second
  clean virtual environment.
- Wheel METADATA, required package contents and every RECORD digest/size.
- Source-distribution path safety and required source/test/audit contents.
- Pickle-free model and prediction paths.
- Anonymous-release scan over source, wheel and sdist.
- Plain-text scan for known identity and common author/contact/path identifiers.
- Two independent builds produced byte-identical wheel and normalized sdist
  artifacts under the same `SOURCE_DATE_EPOCH`.

## Security and operations

- `.wcell` models contain JSON and non-object NumPy arrays only.
- Bundle entry set, nested expansion limits, format version, JSON schema,
  semantic state/configuration consistency, finite arrays, fitted ranges,
  declared shapes and SHA-256 values are validated before state reconstruction.
- Saves are atomic and deterministic; bundle timestamps are fixed.
- GitHub Actions are pinned to full commit SHAs. Publishing uses PyPI Trusted
  Publishing/OIDC and a protected `pypi` environment, with no long-lived token.
- Production publication is tag-only, requires an exact tag/version match, and
  is gated by Python 3.10–3.14 tests, coverage, lint, typing, dependency audit,
  anonymity, reproducibility, distribution integrity and clean-install smoke tests.
- A separate manual TestPyPI rehearsal workflow is included.
- The package performs no telemetry or network access.

## Remaining publisher-owned actions

1. Verify that the `witnesscell` project name is available or controlled on
   TestPyPI and PyPI under an anonymous review account.
2. Configure the Trusted Publisher and manual approval protection for the
   `pypi` environment.
3. Run the included CI on a clean hosted runner and retain the green run URL.
4. Publish to TestPyPI, install and inspect the rendered project page.
5. Create the immutable `v0.1.0` tag and approve the PyPI production job.

No author name, affiliation, email, repository account, DOI, ORCID, grant,
machine-specific path, research manuscript, dataset or pretrained asset is
included in the distribution artifacts.
