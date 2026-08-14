#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${PYTHON_BIN:-python}"
datasets="${DATASETS:-Norman Replogle_exp6 Schmidt Wessels}"
seeds="${SEEDS:-1 2 3}"
output_dir="${OUTPUT_DIR:-$repo_root/artifacts/generated/v14}"
keep_internal_truth="${KEEP_INTERNAL_TRUTH:-0}"
export PYTHONDONTWRITEBYTECODE=1

for dataset in $datasets; do
  case "$dataset" in
    Norman|Replogle_exp6|Schmidt|Wessels) ;;
    *) echo "Unsupported DATASETS entry: $dataset" >&2; exit 2 ;;
  esac
  cache="$repo_root/reproduction/assets/condition_moments/${dataset}.condition_moments.npz"
  [[ -f "$cache" ]] || { echo "Missing cache: $cache" >&2; exit 1; }
  for seed in $seeds; do
    case "$seed" in 1|2|3) ;; *) echo "Unsupported SEEDS entry: $seed" >&2; exit 2 ;; esac
    run_dir="$output_dir/formal_combo/$dataset/seed$seed"
    mkdir -p "$run_dir"
    PYTHONPATH="$repo_root/reproduction/gate19_v14/src" \
      "$python_bin" "$repo_root/reproduction/gate19_v14/src/run_v14_from_cache.py" \
      --cache "$cache" --dataset "$dataset" --seed "$seed" \
      --gears-assets "$repo_root/reproduction/gate15_scperturbench/data/gears_assets" \
      --out "$run_dir"
    [[ -f "$run_dir/deploy_predictions.npz" ]] || {
      echo "Missing generated deployment archive: $run_dir/deploy_predictions.npz" >&2
      exit 1
    }
    if [[ "$keep_internal_truth" != "1" ]]; then
      rm -f "$run_dir/predictions.npz"
    fi
  done
done

manifest="$repo_root/reproduction/assets/checksums/V14_DEPLOY_SEMANTIC.sha256"
"$python_bin" "$repo_root/scripts/semantic_npz_digest.py" \
  --manifest "$manifest" --root "$output_dir/formal_combo" --allow-missing

echo "WITNESSCELL_V14_REPRODUCTION_PASS output=$output_dir"
