#!/usr/bin/env bash
set -euo pipefail

E22="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$E22/../.." && pwd)"
PY="$E22/.venv/bin/python"
RESULTS="$E22/results/formal_e22"
CONTRACTS="$E22/assets/formal_contracts"
SEAL="$E22/assets/CONTRACT_SEAL.json"

if [[ ! -f "$E22/FROZEN_MANIFEST.sha256" ]]; then
  echo "Experiment 22 is not frozen" >&2
  exit 2
fi
(cd "$E22" && sha256sum -c FROZEN_MANIFEST.sha256)
if [[ -e "$RESULTS" ]]; then
  echo "formal output already exists; once-only runner refuses overwrite: $RESULTS" >&2
  exit 3
fi
mkdir -p "$RESULTS/M" "$E22/.cache/matplotlib" "$E22/audit/firewall"

export MPLCONFIGDIR="$E22/.cache/matplotlib"
export PYTHONPATH="$E22/firewall:$E22/src"
export E22_PHASE="M_FORMAL"
export E22_FORBIDDEN_READ_ROOTS
E22_FORBIDDEN_READ_ROOTS="$(printf '["%s","%s","%s"]' \
  "$ROOT/experiments/18_dual_head_evidence_gate" \
  "$ROOT/experiments/19_v14_incremental_amplitude_gate" \
  "$ROOT/experiments/21_frozen_selective_prediction/results")"

for DATASET in Norman Wessels Schmidt Replogle_exp6; do
  AUDIT_LOG="$E22/audit/firewall/formal_${DATASET}.candidate_read_denials.jsonl"
  if [[ -e "$AUDIT_LOG" ]]; then
    echo "stale formal firewall log exists: $AUDIT_LOG" >&2
    exit 4
  fi
  export E22_FIREWALL_AUDIT_LOG="$AUDIT_LOG"
  "$PY" "$E22/src/phase_m_score_controls.py" \
    --contract "$CONTRACTS/$DATASET.phase_m_contract.npz" \
    --split-ledger "$CONTRACTS/$DATASET.split_ledger.json" \
    --contract-seal "$SEAL" \
    --out "$RESULTS/M"
  if [[ -e "$AUDIT_LOG" ]]; then
    echo "formal Phase M attempted a candidate read: $AUDIT_LOG" >&2
    exit 5
  fi
done
"$PY" "$E22/src/aggregate_phase_m.py" --phase-m-dir "$RESULTS/M"

unset E22_PHASE E22_FORBIDDEN_READ_ROOTS E22_FIREWALL_AUDIT_LOG
export PYTHONPATH="$E22/src"
"$PY" "$E22/src/phase_p_prediction_stress.py" \
  --repo "$ROOT" \
  --contracts "$CONTRACTS" \
  --contract-seal "$SEAL" \
  --phase-m-manifest "$RESULTS/M/PHASE_M_MANIFEST.sha256" \
  --scope "$E22/COMPARATOR_SCOPE.json" \
  --out "$RESULTS/P"

"$PY" "$E22/src/phase_d_decision_stress.py" \
  --repo "$ROOT" \
  --gene-contract "$E22/assets/gate21_gene_contract.npz" \
  --out "$RESULTS/D/gate21"

"$PY" "$E22/src/shadow_endpoint_gate.py" \
  --repo "$ROOT" \
  --out "$RESULTS/D/shadow"

"$PY" "$E22/src/aggregate_formal.py" --root "$RESULTS"
echo "PASS_EXPERIMENT22_FORMAL_COMPLETE"
