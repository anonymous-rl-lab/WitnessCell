#!/usr/bin/env python3
"""Bounded calibration replay with results in a new external directory."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
from replay import replay

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    if out == HERE or HERE in out.parents:
        raise RuntimeError('--out must be external to the archived experiment directory')
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('--out must be a new or empty directory')
    manifest = json.loads((HERE / 'MODULE_MANIFEST.json').read_text())
    for row in manifest['files']:
        path = HERE / row['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
            raise RuntimeError('Module checksum mismatch: ' + row['path'])
    report = replay()
    report['module_files_verified'] = len(manifest['files'])
    out.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, allow_nan=False) + '\n'
    (out / 'REPLAY_REPORT.json').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
