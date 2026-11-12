#!/usr/bin/env python3
"""Verify compact release, protected experiments, dependencies and archive CRCs."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile

ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def safe(root,s):
    q=PurePosixPath(s)
    if q.is_absolute() or '..' in q.parts or str(q)!=s:raise ValueError('Unsafe member '+s)
    p=root/s
    if p.is_symlink() or not p.is_file() or not p.resolve().is_relative_to(root):raise ValueError('Missing/unsafe '+s)
    return p

def check(root):
    if (root/'VERSION').read_text().strip()!='27.0.0':raise ValueError('Expected version27.0.0')
    m=json.loads((root/'PACKAGE_MANIFEST.json').read_text());rows=m['files'];paths={x['path'] for x in rows}
    if len(paths)!=len(rows):raise ValueError('Duplicate manifest path')
    for x in rows:
        p=safe(root,x['path'])
        if p.stat().st_size!=x['bytes'] or sha(p)!=x['sha256']:raise ValueError('Release digest mismatch '+x['path'])
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.suffix!='.pyc' and not any(t.startswith('.') or t=='__pycache__' for t in p.relative_to(root).parts)}
    if actual != paths|{'PACKAGE_MANIFEST.json'}:raise ValueError('Unexpected or absent files '+str(sorted(actual^paths^{'PACKAGE_MANIFEST.json'})[:10]))
    protected=json.loads((root/'provenance/PROTECTED_EXPERIMENTS.json').read_text())['files']
    for x in protected:
        p=safe(root,x['path'])
        if p.stat().st_size!=x['bytes'] or sha(p)!=x['sha256']:raise ValueError('Already-compact experiment changed '+x['path'])
    mapping=json.loads((root/'experiments/25_support_controls_and_witness_quantity/CANONICAL_INPUT_MAP.json').read_text())
    count=0
    for x in mapping['files']:
        base=root if x['scope']=='repository' else root/'experiments/25_support_controls_and_witness_quantity'
        p=safe(base,x['path'])
        if sha(p)!=x['sha256']:raise ValueError('E25 dependency changed '+x['path'])
        count+=1
    zipped=[]
    for s in paths:
        if s.endswith('.zip'):
            with zipfile.ZipFile(root/s) as z:
                if z.testzip() is not None or len(z.namelist())!=len(set(z.namelist())):raise ValueError('Archive CRC/members failure '+s)
            zipped.append(s)
    if (root/'reproduction_store').exists():raise ValueError('Historical recovery store must not be included')
    return {'status':'PASS','release':'27.0.0','manifest_files':len(rows),'protected_experiment_files':len(protected),'E25_hash_locked_inputs':count,'embedded_zip_crc_verified':len(zipped),'scope':'Exact bytes and dependency closure; numerical reproduction is a separate command'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path);a=p.parse_args();r=check(ROOT)
    if a.out:
        out=a.out.resolve()
        if out.is_relative_to(ROOT):p.error('Output must be outside release')
        out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))

if __name__=='__main__':main()
