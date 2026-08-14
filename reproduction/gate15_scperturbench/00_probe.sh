#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"

echo "[1/3] Python dependency probe; no installation is performed"
"$PYTHON_BIN" - <<'PY'
import importlib.metadata as md
import numpy, pandas, scipy
import anndata, h5py
print("python dependencies: PASS")
for package in ("numpy", "pandas", "scipy", "anndata", "h5py"):
    print(package, md.version(package))
PY

echo "[2/3] Syntax and synthetic estimator controls"
PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
  "$PYTHON_BIN" -m py_compile "$ROOT"/*.py
PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
  "$PYTHON_BIN" "$ROOT/self_test.py"

echo "[3/3] Protocol identity"
"$PYTHON_BIN" - "$ROOT/protocol.json" <<'PY'
import json, sys
payload=json.load(open(sys.argv[1]))
assert payload["official_seeds"] == [1,2,3]
assert len(payload["formal_combo_datasets"]) == 4
print("PROBE_PASS", payload["method_name"], payload["primary_sota_rule"])
PY
