# Contributing

Thank you for improving WitnessCell. During double-blind review, do not add
names, affiliations, email addresses, repository usernames, DOI values, grant
identifiers, or machine-specific paths to issues, commits, source files,
documentation, metadata, fixtures, or generated artifacts.

Create a focused branch, add tests for behavioral changes, and run:

```bash
python -m pip install -e ".[dev]"
ruff check .
pytest
python scripts/audit_anonymity.py .
python -m build
twine check dist/*
```

Changes to the numerical contract require an explicit contract version, a
reference-parity report, and release notes. Never silently change an evidence
gate, fallback, target role, or risk interpretation.
