#!/usr/bin/env python3
"""Verify release bytes, unchanged manuscript baselines and retained Han files."""
from pathlib import Path
import argparse, hashlib, json

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    rows = json.loads((ROOT/'MANIFEST.json').read_text())['files']
    expected = {row['path'] for row in rows}
    assert len(expected) == len(rows)
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    assert actual == expected | {'MANIFEST.json'}, sorted(actual ^ (expected | {'MANIFEST.json'}))
    for row in rows:
        path = ROOT/row['path']
        assert path.resolve().is_relative_to(ROOT) and not path.is_symlink()
        assert path.stat().st_size == row['bytes'] and digest(path) == row['sha256'], row['path']
    baseline = json.loads((ROOT/'PAPER_BASELINE.json').read_text())
    count = 0
    for row in baseline['manuscript_files']:
        path = ROOT/'paper'/row['path']
        if path.exists():
            assert digest(path) == row['sha256'], row['path']
            count += 1
    assert count == 4, 'Expected unchanged main/SI Markdown and both figure legends'
    han = ROOT/'experiments/han2017'
    ledger = json.loads((han/'COMPACTION_LEDGER.json').read_text())
    for row in ledger['retained_files']:
        assert digest(han/'snapshot'/row['path']) == row['sha256'], row['path']
    assert len(json.loads((han/'RECIPES.json').read_text())['derived_files']) == 180
    report = dict(status='PASS', release='28.0.0', paper='v36',
                  files_verified=len(rows), unchanged_manuscript_files=count,
                  retained_han_files=len(ledger['retained_files']),
                  scope='Stored byte integrity; numerical reproduction is a separate command.')
    if args.out:
        path = args.out.resolve()
        if path.is_relative_to(ROOT): parser.error('Write reports outside the release.')
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
