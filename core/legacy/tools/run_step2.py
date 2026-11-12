#!/usr/bin/env python3
"""Verify or rescore Step 2 without installing packages or retraining models."""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, stat, subprocess, sys, tempfile, time, zipfile, signal
from pathlib import Path, PurePosixPath
from run_growth import run_bounded

EXPERIMENT='experiments/27_step2_contrast_headtohead_compact'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()

def safe_member(name):
    p=PurePosixPath(name.rstrip('/'))
    if not name or '\\' in name or '\x00' in name or ':' in name or p.is_absolute() or '..' in p.parts:
        raise ValueError('Unsafe ZIP member: '+repr(name))
    return p

def inspect_zip(path,root=None):
    with zipfile.ZipFile(path) as z:
        seen=set()
        for m in z.infolist():
            p=safe_member(m.filename)
            if str(p) in seen:raise ValueError('Duplicate ZIP member')
            seen.add(str(p))
            if root and p.parts[0]!=root:raise ValueError('Unexpected archive root')
            if stat.S_IFMT(m.external_attr>>16) not in (0,stat.S_IFREG,stat.S_IFDIR):
                raise ValueError('Special/archive link is forbidden')
        if z.testzip() is not None:raise ValueError('ZIP CRC failed')
    return len(seen)

def validate_archive(repo):
    experiment=repo/EXPERIMENT
    spec=json.loads((experiment/'COMPACT_ARCHIVE.json').read_text())
    archive=experiment/spec['filename']
    if not archive.is_file() or archive.is_symlink():raise ValueError('Missing ordinary compact archive')
    if archive.stat().st_size>=10_000_000:raise ValueError('Step 2 compact ZIP must be below 10,000,000 bytes')
    if archive.stat().st_size!=spec['bytes'] or sha(archive)!=spec['sha256']:raise ValueError('Compact archive identity differs')
    count=inspect_zip(archive,spec['root'])
    return {'archive_bytes':spec['bytes'],'archive_sha256':spec['sha256'],'members':count,'scope':'Exact result rescore from lossless float64 sufficient statistics; original per-gene prediction archives remain separate external source artifacts'}

def extract_checked(path,destination):
    inspect_zip(path)
    destination.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(path) as z:
        for m in z.infolist():
            target=destination/safe_member(m.filename)
            if m.is_dir():target.mkdir(parents=True,exist_ok=True)
            else:
                target.parent.mkdir(parents=True,exist_ok=True)
                with z.open(m) as a,target.open('xb') as b:shutil.copyfileobj(a,b)

def main():
    # Relay an enclosing runner timeout into run_bounded so its separate
    # child process group is terminated as well. No scientific code changes.
    def interrupted(signum, frame):
        raise InterruptedError(f"Interrupted by signal {signum}")
    if os.name == "posix":
        signal.signal(signal.SIGTERM, interrupted)
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--mode',choices=['verify','reproduce'],default='verify')
    p.add_argument('--out',type=Path)
    p.add_argument('--max-seconds',type=float,default=300.)
    a=p.parse_args();repo=a.repo.resolve()
    if (repo/'VERSION').read_text().strip()!='25.0.0':p.error('Expected grn_repo_v25')
    if not 0<a.max_seconds<float('inf'):p.error('Positive finite time limit required')
    if a.mode!='verify' and a.out is None:p.error('--out is required for rescore modes')
    out=a.out.resolve() if a.out else None
    if out and (out.exists() or out.is_relative_to(repo)):p.error('--out must be a new external directory')
    start=time.monotonic();integrity=validate_archive(repo);exp=repo/EXPERIMENT
    spec=json.loads((exp/'COMPACT_ARCHIVE.json').read_text())
    if out:out.mkdir(parents=True)
    report={'passed':False,'repository_version':'25.0.0','mode':a.mode,'compact_archive':integrity,'new_models_trained':0,'installed_software':False,'external_downloads':0}
    try:
        with tempfile.TemporaryDirectory(prefix='witnesscell-v25-step2-') as temp:
            work=out if out else Path(temp)
            extract_checked(exp/spec['filename'],work/'payload')
            payload=work/'payload'/spec['root']
            def run(command,log):run_bounded(command,payload,work/log,a.max_seconds-(time.monotonic()-start))
            if a.mode=='verify':
                run([sys.executable,'-B',str(payload/'verify.py')],'verify.log')
                text=(work/'verify.log').read_text();result=json.loads(text)
                if result.get('passed') is not True:raise ValueError('Compact verifier did not pass')
                report['compact_verification']=result
            elif a.mode=='reproduce':
                run([sys.executable,'-B',str(payload/'reproduce.py'),'--output',str(work/'reproduced')],'reproduce.log')
                result=json.loads((work/'reproduced/REPRODUCTION_CHECK.json').read_text())
                if result.get('passed') is not True:raise ValueError('Compact rescore did not pass')
                report['reproduction']=result
            report['passed']=True
    finally:
        report['elapsed_seconds']=time.monotonic()-start
        if out:(out/'STEP2_REPRODUCTION_CHECK.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2))
    if not report['passed']:raise SystemExit(1)

if __name__=='__main__':main()
