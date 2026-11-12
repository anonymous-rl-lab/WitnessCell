#!/usr/bin/env python3
"""Restore 64 historical Jurkat caches using unchanged, hash-checked source.

python rebuild_jurkat_cache.py --bridge /restored/WitnessCell_NM_NCS_Bridge_20261031 --restore
Only absent cache files are created; original JSON, source, inputs and results
already present are never overwritten. No query-answer file is opened.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import contextlib
import hashlib
import importlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import time

PREFIX='growth_assets/extracted/WitnessCell_Horlbeck_Jurkat_Confirmation_20261018/'
HISTORICAL_MISSING='protocol/v1_before_availability_amendment.zip'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bridge',required=True,type=Path)
    p.add_argument('--restore',action='store_true')
    p.add_argument('--report',type=Path)
    a=p.parse_args();bridge=a.bridge.resolve();source=bridge/PREFIX
    assert bridge.is_dir() and source.is_dir()
    if a.report:
        assert not a.report.resolve().is_relative_to(bridge),'Keep report outside frozen bridge'
    manifest=json.loads((bridge/'RELEASE_MANIFEST.json').read_text())['files']
    def verify_source(rel):
        path=source/rel
        assert PREFIX+rel in manifest,rel
        assert path.is_file() and sha(path)==manifest[PREFIX+rel],rel
        return path
    freeze=json.loads(verify_source('FREEZE.json').read_text())
    lock=json.loads(verify_source('PREDICTIONS_LOCKED.json').read_text())
    provenance=json.loads(verify_source('INPUT_PROVENANCE.json').read_text())
    plan=json.loads(verify_source('PLAN.json').read_text())
    required=['FREEZE.json','INPUT_PROVENANCE.json','PLAN.json']
    unavailable=[]
    for rel,expected in freeze['files'].items():
        path=source/rel
        if not path.exists():
            unavailable.append(rel)
            assert rel==HISTORICAL_MISSING,'Required scientific file absent: '+rel
            assert PREFIX+rel not in manifest,'Historical file unexpectedly missing from restored archive'
            continue
        verify_source(rel)
        assert sha(path)==expected,rel
        required.append(rel)
    for rel in ('inputs/training.npz','inputs/query.npz'):
        verify_source(rel)
        assert sha(source/rel)==provenance['file_sha256'][rel],rel
        required.append(rel)
    targets={n[len(PREFIX):]:h for n,h in manifest.items() if n.startswith(PREFIX+'results/high/') and
             (('/inner_models/' in n and n.endswith('.npz')) or n.endswith('_inner_predictions.npz'))}
    assert len(targets)==64, len(targets)
    assert sum('/inner_models/' in n for n in targets)==48
    for rel,expected in targets.items():
        assert lock['files'][rel]==expected,rel
        if (source/rel).exists(): assert sha(source/rel)==expected,'Existing cache differs: '+rel
    missing=[rel for rel in targets if not (source/rel).exists()]
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='wc_jurkat_cache_') as tmp:
        work=Path(tmp)
        for rel in sorted(set(required)):
            dst=work/rel;dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source/rel,dst)
        (work/'results/high/inner_models').mkdir(parents=True)
        sys.path.insert(0,str(work/'core'))
        import numpy as np
        frozen=importlib.import_module('run_smoke')
        assert Path(frozen.__file__).resolve()==(work/'core/run_smoke.py').resolve()
        with np.load(work/'inputs/training.npz',allow_pickle=False) as z:
            data={k:z[k] for k in z.files}
        with np.load(work/'inputs/query.npz',allow_pickle=False) as z:
            assert 'y' not in z.files
            edges=z['edges'];ids=z['global_ids']
        assert not set(ids)&set(data['global_ids'])
        frozen.ROOT=work
        with contextlib.redirect_stdout(io.StringIO()):
            frozen.run_budget('high',data,edges,ids,plan,time.monotonic()+plan['real_fit_budget_seconds'])
        for rel,expected in targets.items(): assert sha(work/rel)==expected,'Rebuilt SHA mismatch: '+rel
        if a.restore:
            for rel in missing:
                target=source/rel;target.parent.mkdir(parents=True,exist_ok=True)
                assert not target.exists(),'Concurrent cache write: '+rel
                shutil.copy2(work/rel,target)
                assert sha(target)==targets[rel]
    result={'status':'PASS','cache_files_exact':len(targets),'restored_files':len(missing) if a.restore else 0,
            'mode':'restore_missing_only' if a.restore else 'rebuild_and_verify_only',
            'seconds':time.monotonic()-started,'scientific_source_modified':False,'query_answers_read':False,
            'source_archive_existing_boundary':{'unavailable_historical_protocol_files':unavailable,
               'note':'Historical availability snapshot ZIP was already absent from original bridge archive; all present frozen scientific source, plan and permitted inputs are SHA-checked.'},
            'expected_sha256':targets}
    if a.report:
        a.report.parent.mkdir(parents=True,exist_ok=True)
        a.report.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='expected_sha256'},indent=2))

if __name__=='__main__':main()
