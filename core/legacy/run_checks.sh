#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 -B tools/verify.py
if [ "$#" -gt 0 ]; then
 python3 -B tools/reproduce.py --out "$1"
fi
