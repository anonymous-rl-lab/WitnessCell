#!/usr/bin/env python3
"""One entry for the unchanged frozen Combi-seq offline refit."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    if out == ROOT or ROOT in out.parents:
        parser.error('Use an output directory outside the immutable experiment module.')
    if out.exists() and any(out.iterdir()):
        parser.error('Use a new or empty output directory.')
    subprocess.run([sys.executable, str(ROOT / 'REPRODUCE.py'), '--out', str(out)], check=True)


if __name__ == '__main__':
    main()
