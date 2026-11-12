#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: bash run_preregistered_replication.sh /absolute/path/grn_repo_v14.zip [output_parent]" >&2
  exit 2
fi

SOURCE_ZIP="$1"
OUTPUT_PARENT="${2:-$PWD}"
EXPECTED_SHA="955d446a0cfe65a6f4f191d74c7a7d16a6083488b9beb40a5e155ab94d528bd7"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -f "$SOURCE_ZIP" ]; then
  echo "Source archive not found: $SOURCE_ZIP" >&2
  exit 3
fi

ACTUAL_SHA="$(sha256sum "$SOURCE_ZIP" | awk '{print $1}')"
if [ "$ACTUAL_SHA" != "$EXPECTED_SHA" ]; then
  echo "FAIL_SOURCE_SHA256 expected=$EXPECTED_SHA actual=$ACTUAL_SHA" >&2
  exit 4
fi

PYTHON_BIN="${PYTHON_BIN:-python3.12}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Python 3.12 not found. Set PYTHON_BIN to the Python 3.12 executable." >&2
  exit 5
fi

"$PYTHON_BIN" "$SCRIPT_DIR/src/check_environment.py"

RUN_TAG="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_ROOT="$OUTPUT_PARENT/WitnessCell_preregistered_replication_$RUN_TAG"
mkdir -p "$RUN_ROOT/source" "$RUN_ROOT/gate19_regenerated" "$RUN_ROOT/results"
date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_ROOT/run_started_utc.txt"

cp "$SCRIPT_DIR/FROZEN_PROTOCOL.json" "$RUN_ROOT/FROZEN_PROTOCOL.json"
sha256sum "$SOURCE_ZIP" > "$RUN_ROOT/source_archive.sha256"
"$PYTHON_BIN" --version > "$RUN_ROOT/environment.txt" 2>&1
"$PYTHON_BIN" -m pip freeze >> "$RUN_ROOT/environment.txt"

unzip -q "$SOURCE_ZIP" -d "$RUN_ROOT/source"
REPO="$RUN_ROOT/source/grn_repo_v14"
GATE19="$REPO/experiments/19_v14_incremental_amplitude_gate"

for dataset in Norman Replogle_exp6 Schmidt Wessels; do
  for seed in 1 2 3; do
    echo "[regenerate] dataset=$dataset seed=$seed"
    PYTHONPATH="$GATE19/src" "$PYTHON_BIN" "$GATE19/src/run_v14_from_cache.py" \
      --cache "$REPO/experiments/18_dual_head_evidence_gate/cache/${dataset}.condition_moments.npz" \
      --dataset "$dataset" \
      --seed "$seed" \
      --gears-assets "$REPO/experiments/15_scperturbench_sota/module/data/gears_assets" \
      --out "$RUN_ROOT/gate19_regenerated/${dataset}/seed${seed}" \
      | tee "$RUN_ROOT/results/regeneration_${dataset}_seed${seed}.log"
  done
done

PYTHONPATH="$GATE19/src" "$PYTHON_BIN" "$GATE19/src/compare_v14_v13.py" \
  --v14-root "$RUN_ROOT/gate19_regenerated" \
  --v13-root "$REPO/experiments/18_dual_head_evidence_gate/predictions/formal_combo" \
  --out "$RUN_ROOT/results/v14_vs_v13" \
  > "$RUN_ROOT/results/v14_vs_v13.log" 2>&1

"$PYTHON_BIN" "$SCRIPT_DIR/src/confirmatory_audit.py" \
  --repo "$REPO" \
  --regenerated "$RUN_ROOT/gate19_regenerated" \
  --paired "$RUN_ROOT/results/v14_vs_v13/verdict.json" \
  --out "$RUN_ROOT/CONFIRMATORY_VERDICT.json" \
  | tee "$RUN_ROOT/confirmatory_audit.log"

date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_ROOT/run_completed_utc.txt"
RESULT_ARCHIVE="$OUTPUT_PARENT/WitnessCell_preregistered_replication_${RUN_TAG}.tar.gz"
RUN_BASENAME="$(basename "$RUN_ROOT")"
tar \
  --exclude="$RUN_BASENAME/source" \
  --exclude="$RUN_BASENAME/gate19_regenerated/*/seed*/predictions.npz" \
  -czf "$RESULT_ARCHIVE" -C "$OUTPUT_PARENT" "$RUN_BASENAME"
sha256sum "$RESULT_ARCHIVE" | tee "$RESULT_ARCHIVE.sha256"
echo "RESULT_ARCHIVE=$RESULT_ARCHIVE"
