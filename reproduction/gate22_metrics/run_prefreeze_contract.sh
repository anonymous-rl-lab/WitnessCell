#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
E22="$ROOT/experiments/22_metric_calibration_stress"
PY="$E22/.venv/bin/python"
DATASET="${1:?usage: run_prefreeze_contract.sh DATASET}"

case "$DATASET" in
  Norman|Wessels|Schmidt) DATA_FILE="$DATASET.h5ad" ;;
  Replogle_exp6) DATA_FILE="Replogle_exp6.h5ad" ;;
  *) echo "unsupported dataset: $DATASET" >&2; exit 2 ;;
esac

mkdir -p "$E22/.cache/matplotlib" "$E22/audit/firewall"
AUDIT_LOG="$E22/audit/firewall/${DATASET}.candidate_read_denials.jsonl"
if [[ -e "$AUDIT_LOG" ]]; then
  echo "refusing to overwrite existing firewall audit: $AUDIT_LOG" >&2
  exit 3
fi

export MPLCONFIGDIR="$E22/.cache/matplotlib"
export E22_PHASE="M_PREFREEZE_CONTRACT"
export E22_FORBIDDEN_READ_ROOTS
E22_FORBIDDEN_READ_ROOTS="$(printf '["%s","%s","%s"]' \
  "$ROOT/experiments/18_dual_head_evidence_gate" \
  "$ROOT/experiments/19_v14_incremental_amplitude_gate" \
  "$ROOT/experiments/21_frozen_selective_prediction/results")"
export E22_FIREWALL_AUDIT_LOG="$AUDIT_LOG"
export PYTHONPATH="$E22/firewall:$E22/src"

"$PY" "$E22/src/phase_m_metric_validity.py" \
  --data "$E22/assets/cell_level/$DATA_FILE" \
  --dataset "$DATASET" \
  --official-module "$ROOT/experiments/15_scperturbench_sota/module" \
  --gears-assets "$ROOT/experiments/15_scperturbench_sota/module/data/gears_assets" \
  --out "$E22/assets/formal_contracts" \
  --prepare-only

if [[ -e "$AUDIT_LOG" ]]; then
  echo "candidate-read denial occurred; contract is invalid: $AUDIT_LOG" >&2
  exit 4
fi
echo "PASS_${DATASET}_PREFREEZE_CONTRACT_FIREWALL_CLEAN"
