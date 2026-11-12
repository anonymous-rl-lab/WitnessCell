#!/usr/bin/env python3
"""Atomic packaging with a strict100,000,000-byte ceiling and complete member verification."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
LIMIT=100_000_000

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=ROOT.parent/'grn_repo_v25.zip');p.add_argument('--manifest-only',action='store_true');a=p.parse_args();out=a.out.resolve()
    if out.is_relative_to(ROOT):p.error('Output must be outside release')
    files=[]
    for x in sorted(ROOT.rglob('*')):
        rel=x.relative_to(ROOT)
        if any(t.startswith('.') or t=='__pycache__' for t in rel.parts) or x.suffix=='.pyc':continue
        if x.is_symlink():raise ValueError('Symlinks are not packaged')
        if x.is_file() and rel.as_posix()!='PACKAGE_MANIFEST.json':files.append(x)
    rows=[{'path':x.relative_to(ROOT).as_posix(),'bytes':x.stat().st_size,'sha256':sha(x)} for x in files]
    manifest={'schema':'witnesscell-compact-v25-v1','version':'25.0.0','paper':'22','archive_limit_bytes_exclusive':LIMIT,'files':rows}
    (ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if a.manifest_only:print(json.dumps({'status':'MANIFEST_CREATED','files':len(rows)}));return
    subprocess.run([sys.executable,'-B',str(ROOT/'tools/verify.py')],check=True)
    files.append(ROOT/'PACKAGE_MANIFEST.json');out.parent.mkdir(parents=True,exist_ok=True);tmp=out.with_suffix('.complete.tmp')
    with zipfile.ZipFile(tmp,'w',allowZip64=True) as z:
        for x in files:
            info=zipfile.ZipInfo('grn_repo_v25/'+x.relative_to(ROOT).as_posix(),(2026,10,26,0,0,0));info.external_attr=(0o100755 if x.suffix=='.sh' else 0o100644)<<16
            info.compress_type=zipfile.ZIP_STORED if x.suffix=='.zip' else zipfile.ZIP_DEFLATED
            z.writestr(info,x.read_bytes(),compresslevel=9)
    if tmp.stat().st_size>=LIMIT:raise ValueError('Budget exceeded: '+str(tmp.stat().st_size))
    with zipfile.ZipFile(tmp) as z:
        if z.testzip() is not None or len(z.namelist())!=len(files):raise ValueError('ZIP CRC/member count failure')
        for x in rows:
            b=z.read('grn_repo_v25/'+x['path'])
            if len(b)!=x['bytes'] or hashlib.sha256(b).hexdigest()!=x['sha256']:raise ValueError('Archived bytes differ '+x['path'])
    with tmp.open('rb') as f:os.fsync(f.fileno())
    os.replace(tmp,out);rec={'status':'PASS','file':out.name,'bytes':out.stat().st_size,'sha256':sha(out),'members':len(files),'verified_files':len(rows),'limit_exclusive':LIMIT,'under_limit':True,'atomic_install':True}
    out.with_name(out.name+'.sha256').write_text(rec['sha256']+'  '+out.name+'\n');out.with_name('grn_repo_v25_build_receipt.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))

if __name__=='__main__':main()
