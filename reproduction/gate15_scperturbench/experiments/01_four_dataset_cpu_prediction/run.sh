#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DATA_ROOT="${1:?usage: $0 DATA_ROOT [OUTPUT_ROOT] [GEARS_ASSETS]}"
OUTPUT_ROOT="${2:-$ROOT/results/formal_combo}"
GEARS_ASSETS="${3:-$ROOT/data/gears_assets}"

bash "$ROOT/02_run_combo_predictions.sh" "$DATA_ROOT" "$OUTPUT_ROOT" "$GEARS_ASSETS"
python "$HERE/audit.py" --package-root "$ROOT"
