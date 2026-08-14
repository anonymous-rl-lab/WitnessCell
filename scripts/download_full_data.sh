#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_dir="${1:-$repo_root/data/raw}"

bash "$repo_root/reproduction/gate15_scperturbench/01_download_official_combo_data.sh" \
  "$output_dir"
