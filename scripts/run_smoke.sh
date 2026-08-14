#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${PYTHON_BIN:-python}"
export PYTHONDONTWRITEBYTECODE=1

cd "$repo_root"

"$python_bin" -m pytest -q -p no:cacheprovider
"$python_bin" scripts/audit_compact_results.py
"$python_bin" reproduction/gate15_scperturbench/experiments/02_official_full_metric_sota/audit.py
"$python_bin" reproduction/gate21_selective/audit.py
"$python_bin" reproduction/gate21_selective/post_reveal_integrity_audit.py --compact-release

PYTHONPATH="$repo_root/reproduction/gate22_metrics/src" \
  "$python_bin" -m unittest discover \
  -s reproduction/gate22_metrics/tests -p 'test_metric_core.py'
PYTHONPATH="$repo_root/reproduction/gate22_metrics/src" \
  "$python_bin" -m unittest discover \
  -s reproduction/gate22_metrics/tests -p 'test_comparator_core.py'

"$python_bin" reproduction/theory/06_exact_risk_optimality/verify_exact_risk.py \
  --replicates 5000 --random-competitors 50 --out /tmp/witnesscell-exact-risk-smoke
"$python_bin" reproduction/theory/07_estimated_witness_norman/validate_estimator.py \
  --seeds 6 --out /tmp/witnesscell-estimator-smoke

"$python_bin" scripts/audit_anonymity.py .

echo "WITNESSCELL_GITHUB_SMOKE_PASS"
