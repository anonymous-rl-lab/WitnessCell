#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${1:?usage: $0 DATA_ROOT [OUTPUT_ROOT] [GEARS_ASSETS]}"
OUTPUT_ROOT="${2:-$ROOT/results/formal_combo}"
GEARS_ASSETS="${3:-$ROOT/data/gears_assets}"
PYTHON_BIN="${PYTHON_BIN:-python}"
mkdir -p "$OUTPUT_ROOT"

bash "$ROOT/00_probe.sh"
for DATASET in Norman Wessels Schmidt Replogle_exp6; do
  DATA="$DATA_ROOT/$DATASET.h5ad"
  [[ -f "$DATA" ]] || { echo "missing $DATA" >&2; exit 2; }
  for SEED in 1 2 3; do
    echo "[$DATASET seed $SEED/3] WitnessCell CPU fit + prediction"
    PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
      "$PYTHON_BIN" "$ROOT/run_witnesscell.py" \
      --data "$DATA" --dataset "$DATASET" --seed "$SEED" \
      --gears-assets "$GEARS_ASSETS" \
      --out "$OUTPUT_ROOT/$DATASET/seed$SEED"
  done
done

echo "PREDICTIONS_COMPLETE $OUTPUT_ROOT"
