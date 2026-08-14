#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${1:?usage: $0 DATA_ROOT PREDICTION_ROOT [OUTPUT_ROOT]}"
PREDICTION_ROOT="${2:?usage: $0 DATA_ROOT PREDICTION_ROOT [OUTPUT_ROOT]}"
OUTPUT_ROOT="${3:-$ROOT/results/official_h5ad}"
PYTHON_BIN="${PYTHON_BIN:-python}"

for DATASET in Norman Wessels Schmidt Replogle_exp6; do
  for SEED in 1 2 3; do
    echo "[$DATASET seed $SEED/3] export official result.h5ad"
    TARGET="$OUTPUT_ROOT/DataSet2/$DATASET/hvg5000/WitnessCell/savedModels$SEED/result.h5ad"
    PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
      "$PYTHON_BIN" "$ROOT/export_result_h5ad.py" \
      --data "$DATA_ROOT/$DATASET.h5ad" \
      --prediction "$PREDICTION_ROOT/$DATASET/seed$SEED/deploy_predictions.npz" \
      --out "$TARGET"
  done
done

echo "OFFICIAL_H5AD_EXPORT_COMPLETE $OUTPUT_ROOT"
