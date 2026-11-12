#!/usr/bin/env python3
from importlib.metadata import PackageNotFoundError, version
import sys

EXPECTED = {
    "numpy": "1.26.4",
    "pandas": "2.2.3",
    "scipy": "1.14.1",
    "scikit-learn": "1.5.2",
    "joblib": "1.5.3",
    "threadpoolctl": "3.6.0",
    "python-dateutil": "2.9.0.post0",
    "pytz": "2026.3.post1",
    "tzdata": "2026.3",
    "six": "1.17.0",
}

if sys.version_info[:2] != (3, 12):
    raise SystemExit(f"FAIL_PYTHON_VERSION expected=3.12 actual={sys.version.split()[0]}")

bad = []
for package, expected in EXPECTED.items():
    try:
        actual = version(package)
    except PackageNotFoundError:
        bad.append((package, expected, "MISSING"))
        continue
    if actual != expected:
        bad.append((package, expected, actual))

if bad:
    details = "; ".join(f"{p}: expected {e}, actual {a}" for p, e, a in bad)
    raise SystemExit("FAIL_ENVIRONMENT " + details)

print("PASS_FROZEN_ENVIRONMENT")
