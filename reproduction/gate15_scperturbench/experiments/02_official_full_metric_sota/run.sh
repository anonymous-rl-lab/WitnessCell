#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DATA_ROOT="${1:?usage: $0 DATA_ROOT [PREDICTION_ROOT] [SCORE_ROOT]}"
PREDICTION_ROOT="${2:-$ROOT/results/formal_combo}"
SCORE_ROOT="${3:-$ROOT/results/formal_score/full}"

bash "$ROOT/04_score_official_full.sh" "$DATA_ROOT" "$PREDICTION_ROOT" "$SCORE_ROOT"
bash "$ROOT/05_aggregate_official_full.sh" "$SCORE_ROOT" "$ROOT/results/formal_score/aggregate"
python "$HERE/audit.py" --package-root "$ROOT"
