#!/usr/bin/env python3
"""Read-only verification of every file covered by the release manifest."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    manifest=json.loads((ROOT/'MANIFEST_SHA256.json').read_text())
    failed=[]
    for item in manifest['files']:
        p=ROOT/item['path']
        if not p.is_file() or p.stat().st_size!=item['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:failed.append(item['path'])
    print(json.dumps({'checked':len(manifest['files']),'failed':failed,'all_passed':not failed}))
    if failed:raise SystemExit(1)
if __name__=='__main__':main()
