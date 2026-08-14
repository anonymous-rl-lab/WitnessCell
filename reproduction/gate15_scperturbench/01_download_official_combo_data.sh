#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$ROOT/data}"
PYTHON_BIN="${PYTHON_BIN:-python}"
mkdir -p "$OUT"

while IFS=$'\t' read -r DATASET BYTES MD5 URL; do
  [[ "$DATASET" == "dataset" ]] && continue
  ARCHIVE="$OUT/$DATASET.h5ad.gz"
  DATA="$OUT/$DATASET.h5ad"
  if [[ -f "$ARCHIVE" ]] && [[ "$(md5sum "$ARCHIVE" | awk '{print $1}')" == "$MD5" ]]; then
    echo "REUSE verified $ARCHIVE"
  else
    echo "DOWNLOAD $DATASET ($BYTES bytes)"
    curl -L --fail --retry 5 --continue-at - --output "$ARCHIVE" "$URL"
    echo "$MD5  $ARCHIVE" | md5sum -c -
  fi
  H5_OK=0
  if [[ -f "$DATA" ]]; then
    if "$PYTHON_BIN" - "$DATA" <<'PY' >/dev/null 2>&1
import h5py, sys
with h5py.File(sys.argv[1], "r") as handle:
    assert "X" in handle and "obs" in handle and "var" in handle
PY
    then
      H5_OK=1
      echo "REUSE validated $DATA"
    else
      echo "INVALID derived h5ad; rebuild $DATA"
    fi
  fi
  if [[ "$H5_OK" -eq 0 ]]; then
    echo "DECOMPRESS $ARCHIVE"
    TMP="$DATA.partial"
    gzip -dc "$ARCHIVE" > "$TMP"
    "$PYTHON_BIN" - "$TMP" <<'PY'
import h5py, sys
with h5py.File(sys.argv[1], "r") as handle:
    assert "X" in handle and "obs" in handle and "var" in handle
print("H5AD_VALID", sys.argv[1])
PY
    mv "$TMP" "$DATA"
  fi
done < "$ROOT/data_manifest.tsv"

echo "OFFICIAL_COMBO_DATA_READY $OUT"
