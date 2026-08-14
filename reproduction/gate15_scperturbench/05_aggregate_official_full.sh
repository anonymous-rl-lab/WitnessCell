#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCORE_ROOT="${1:-$ROOT/results/formal_score/full}"
OUT="${2:-$ROOT/results/formal_score/aggregate}"
PYTHON_BIN="${PYTHON_BIN:-python}"

PYTHONPATH="$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
  "$PYTHON_BIN" "$ROOT/aggregate_official_full.py" \
    --published-raw "$ROOT/reference/all_combo_genetic_published_raw.csv.gz" \
    --published-top100 "$ROOT/reference/genetic_combo_performance_top100.csv" \
    --published-top5000 "$ROOT/reference/genetic_combo_performance_top5000.csv" \
    --witness-raw "$SCORE_ROOT/Norman/official_six_metric_raw.csv" \
    --witness-raw "$SCORE_ROOT/Wessels/official_six_metric_raw.csv" \
    --witness-raw "$SCORE_ROOT/Schmidt/official_six_metric_raw.csv" \
    --witness-raw "$SCORE_ROOT/Replogle_exp6/official_six_metric_raw.csv" \
    --out "$OUT"

"$PYTHON_BIN" "$ROOT/make_formal_main_figure.py" --aggregate "$OUT" --out "$ROOT/figures"
echo "FORMAL_AGGREGATION_COMPLETE $OUT"
