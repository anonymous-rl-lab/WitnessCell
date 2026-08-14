#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${1:?usage: $0 DATA_ROOT [PREDICTION_ROOT] [SCORE_ROOT]}"
PREDICTION_ROOT="${2:-$ROOT/results/formal_combo}"
SCORE_ROOT="${3:-$ROOT/results/formal_score/full}"
PYTHON_BIN="${PYTHON_BIN:-python}"
THREADS="${WITNESSCELL_SCORE_THREADS:-4}"

for DATASET in Norman Wessels Schmidt Replogle_exp6; do
  echo "[$DATASET] official full-condition scoring; existing seed checkpoints are reused"
  OPENBLAS_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
  JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false MPLCONFIGDIR=/tmp/mpl \
  PYTHONPATH="$ROOT${SCORER_PYTHONPATH:+:$SCORER_PYTHONPATH}${PYTHONPATH:+:$PYTHONPATH}" \
    "$PYTHON_BIN" "$ROOT/score_official_full.py" \
      --data "$DATA_ROOT/$DATASET.h5ad" \
      --prediction "$PREDICTION_ROOT/$DATASET/seed1/deploy_predictions.npz" \
      --prediction "$PREDICTION_ROOT/$DATASET/seed2/deploy_predictions.npz" \
      --prediction "$PREDICTION_ROOT/$DATASET/seed3/deploy_predictions.npz" \
      --out "$SCORE_ROOT/$DATASET"
done

echo "OFFICIAL_FULL_SCORING_COMPLETE $SCORE_ROOT"
