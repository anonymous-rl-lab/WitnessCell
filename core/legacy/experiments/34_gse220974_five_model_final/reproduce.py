#!/usr/bin/env python3
"""Offline exact-input replay of the completed GSE220974 five-model experiment.

Default: all 420 choices, five-model RNA/partner metrics, six paired intervals.
Optional --reference-zip: additionally replay all 15 historical/seed result rows
from the lossless full prediction archive. No model training or network access.
"""
import argparse
import csv
import hashlib
import io
import json
import lzma
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import time
import zipfile
import numpy as np

ROOT = Path(__file__).resolve().parent
MAIN = ['R_G', 'GEARS_Ensemble', 'Scouter_Ensemble', 'StandardRidge_RNA', 'SingleB_additive']

def digest_bytes(b):
    return hashlib.sha256(b).hexdigest()

def digest(p):
    return digest_bytes(p.read_bytes())

def write_json(p, x):
    p.write_text(json.dumps(x, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')

def verify_package():
    m = json.loads((ROOT / 'MODULE_MANIFEST.json').read_text())
    expected = {x['path'] for x in m['files']}
    if len(expected) != len(m['files']):
        raise ValueError('Duplicate manifest member')
    for row in m['files']:
        q = PurePosixPath(row['path'])
        if q.is_absolute() or '..' in q.parts:
            raise ValueError('Unsafe manifest member')
        p = ROOT / q
        if p.is_symlink() or not p.is_file() or not p.resolve().is_relative_to(ROOT):
            raise ValueError('Missing or unsafe file: ' + str(q))
        if p.stat().st_size != row['bytes'] or digest(p) != row['sha256']:
            raise ValueError('Package hash mismatch: ' + str(q))
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    if actual != expected | {'MODULE_MANIFEST.json'}:
        raise ValueError('Unexpected or missing package files')
    return len(expected)

def load_compact():
    meta = json.loads((ROOT / 'data/ARRAYS.json').read_text())
    raw = lzma.decompress((ROOT / 'data/arrays.xz').read_bytes())
    result = {}; offset = 0
    for row in meta:
        dtype = np.dtype(row['dtype']); size = row['bytes']
        shuffled = np.frombuffer(raw[offset:offset + size], dtype=np.uint8)
        restored = shuffled.reshape(dtype.itemsize, -1).T.copy().tobytes()
        if digest_bytes(restored) != row['sha256']:
            raise ValueError('Restored array SHA256 mismatch: ' + row['key'])
        result[row['key']] = np.frombuffer(restored, dtype=dtype).reshape(row['shape'])
        offset += size
    if offset != len(raw):
        raise ValueError('Unexpected decoded bytes')
    return result

def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def compare_csv(actual, expected, keys, columns=None):
    def read(p):
        with p.open(newline='', encoding='utf-8') as f:
            return {tuple(r[k] for k in keys): r for r in csv.DictReader(f)}
    aa, ee = read(actual), read(expected)
    if set(aa) != set(ee):
        raise ValueError('Different rows: ' + actual.name)
    maximum = 0.; numeric = 0; texts = 0
    for key, row in aa.items():
        other = ee[key]
        for k in (columns or row.keys()):
            a, e = row[k], other[k]
            if a == e:
                texts += 1; continue
            try:
                x, y = float(a), float(e)
            except ValueError:
                raise ValueError('Different text: ' + str((key, k, a, e)))
            if not np.isfinite(x) or not np.isfinite(y):
                raise ValueError('Non-finite comparison')
            np.testing.assert_allclose(x, y, rtol=1e-12, atol=1e-14)
            maximum = max(maximum, abs(x-y)); numeric += 1
    return {'rows':len(aa), 'maximum_absolute_difference':maximum,
            'byte_identical_csv':digest(actual)==digest(expected),
            'equal_formatted_fields':texts, 'numerically_checked_different_fields':numeric}

def score(d, out):
    names=d['targets'].astype(str); pairs=d['conditions'].astype(str)
    x=np.array([[n in p.split('+') for n in names] for p in pairs],float)
    geometry=[]; groups=[]
    for fold in range(5):
        ids=np.flatnonzero(d['fold']==fold); rank=np.linalg.matrix_rank(x[ids])
        left,_,_=np.linalg.svd(x[ids],full_matrices=True); basis=left[:,rank:]
        projection=np.eye(len(ids))-x[ids]@np.linalg.pinv(x[ids])
        np.testing.assert_allclose(basis@basis.T,projection,atol=1e-12)
        geometry.append((ids,basis,projection))
        for a,anchor in enumerate(names):
            pool=ids[x[ids,a]>0]
            if len(pool)>=2:
                pool=np.array(sorted(pool,key=lambda i:pairs[i]))
                groups.append((fold,a,anchor,pool))
    models=sorted(k[3:] for k in d if k.startswith('U__'))
    if len(models)!=15 or len(groups)!=28:
        raise ValueError('Wrong frozen cohort')
    records=[]; summaries={}; degrees=sum(b.shape[1] for _,b,_ in geometry)
    for model in models:
        rows=[]; pred=d['U__'+model]
        for f,a,anchor,pool in groups:
            utilities=pred[pool,a]
            chosen=pool[np.flatnonzero(utilities>=utilities.max()-1e-12)[0]]
            actual=d['truth_U'][chosen,a]; best=d['truth_U'][pool,a].max()
            row=dict(model=model,group=f'f{f}:{anchor}',fold=f,anchor=anchor,
                     physical_anchor=anchor.rsplit('-',1)[0],candidate_count=len(pool),pair=pairs[chosen],
                     actual_U=float(actual),best_U=float(best),regret=float(best-actual),
                     hit=int(actual>=best-1e-12),half0_U=float(d['half_U'][0,chosen,a]),
                     half1_U=float(d['half_U'][1,chosen,a]))
            rows.append(row);records.append(row)
        if model not in MAIN:
            continue
        large=[r for r in rows if r['candidate_count']>=4]
        err=d['ERR__'+model]
        svd_score=sum(np.sum((b.T@err[ids])**2) for ids,b,_ in geometry)/(degrees*err.shape[1])
        projection_score=sum(np.sum((p@err[ids])**2) for ids,_,p in geometry)/(degrees*err.shape[1])
        np.testing.assert_allclose(svd_score,projection_score,rtol=1e-12,atol=1e-14)
        summaries[model]=dict(model=model,tasks=len(rows),mean_U=float(np.mean([r['actual_U'] for r in rows])),
            mean_regret=float(np.mean([r['regret'] for r in rows])),hits=sum(r['hit'] for r in rows),
            RNA_MSE=float(np.mean(err**2)),partner_MSE=float(svd_score),partner_degrees=degrees,
            large_tasks=len(large),large_regret=float(np.mean([r['regret'] for r in large])),
            half0_U=float(np.mean([r['half0_U'] for r in rows])),half1_U=float(np.mean([r['half1_U'] for r in rows])))
    write_csv(out/'MAIN_FIVE_MODELS.csv',[summaries[n] for n in MAIN])
    write_csv(out/'CHOICES.csv',records)
    groups_sorted=sorted({r['group'] for r in records})
    by={(r['model'],r['group']):r for r in records}
    anchors=np.array([by[('R_G',g)]['physical_anchor'] for g in groups_sorted])
    keys=sorted(set(anchors));rng=np.random.default_rng(20261102)
    draws=rng.integers(0,len(keys),size=(10000,len(keys)));comparisons=[]
    for ref in ['PRESAGE_Combo_seed24','StandardRidge_RNA','StandardRidge_Choice','SingleB_additive','GEARS_Ensemble','Scouter_Ensemble']:
        delta=np.array([by[('R_G',g)]['actual_U']-by[(ref,g)]['actual_U'] for g in groups_sorted])
        sums=np.array([delta[anchors==k].sum() for k in keys]);counts=np.array([np.sum(anchors==k) for k in keys])
        boot=sums[draws].sum(1)/counts[draws].sum(1)
        comparisons.append(dict(model='R_G',reference=ref,mean_paired_U_difference=float(delta.mean()),
            cluster_bootstrap_low=float(np.quantile(boot,.025)),cluster_bootstrap_high=float(np.quantile(boot,.975)),
            physical_anchor_clusters=len(keys)))
    write_csv(out/'PAIRED_COMPARISONS.csv',comparisons)
    checks={n:compare_csv(out/n,ROOT/'reference'/n,k) for n,k in [
        ('MAIN_FIVE_MODELS.csv',['model']),('CHOICES.csv',['model','group']),
        ('PAIRED_COMPARISONS.csv',['model','reference'])]}
    return checks

def full_replay(reference, out, compact):
    source=json.loads((ROOT/'EXTERNAL_ASSETS.json').read_text())['full_comparison']
    if reference.stat().st_size!=source['bytes'] or digest(reference)!=source['sha256']:
        raise ValueError('Full prediction archive identity mismatch')
    folder=out/'full_prediction_replay';folder.mkdir()
    with zipfile.ZipFile(reference) as z:
        for name in ['score.py','scoring_inputs.npz']:
            b=z.read(name)
            if digest_bytes(b)!=source['members'][name]:
                raise ValueError('Full archive member hash mismatch')
            (folder/name).write_bytes(b)
    with np.load(folder/'scoring_inputs.npz',allow_pickle=False) as z:
        original={k:z[k] for k in z.files}
    checks=[]
    for k,a in compact.items():
        if k.startswith('ERR__'):
            b=original['RNA__'+k[5:]]-original['truth']
        elif k=='truth_U':
            b=original['truth']@original['objectives'].T
        else:
            b=original[k]
        equal=a.dtype==b.dtype and a.shape==b.shape and a.tobytes()==b.tobytes()
        difference = float(np.max(np.abs(a-b))) if a.dtype.kind=='f' else 0.0
        if not equal:
            if k != 'truth_U':
                raise ValueError('Compact/source array mismatch: '+k)
            # Only this derived matrix product can vary at its last bits with
            # BLAS thread count. Stored arrays themselves remain SHA256 exact.
            np.testing.assert_allclose(a,b,rtol=1e-12,atol=1e-14)
        checks.append({'array':k,'bitwise_identical':equal,'maximum_absolute_difference':difference,
                       'sha256':digest_bytes(a.tobytes())})
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    with (folder/'score.log').open('w') as log:
        subprocess.run([sys.executable,'-B',str(folder/'score.py')],cwd=folder,env=env,check=True,stdout=log,stderr=subprocess.STDOUT)
    result={n:compare_csv(folder/n,ROOT/'reference'/n,k) for n,k in [
        ('SUMMARY.csv',['model']),('MAIN_FIVE_MODELS.csv',['model']),('CHOICES.csv',['model','group'])]}
    return {'status':'PASS','compact_source_arrays':checks,'full_15_model_replay':result}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--reference-zip',type=Path);a=p.parse_args();out=a.out.resolve()
    if out.exists() or out.is_relative_to(ROOT):
        p.error('Output directory must be new and outside the release')
    count=verify_package();d=load_compact();out.mkdir(parents=True)
    start=time.monotonic();checks=score(d,out)
    report={'status':'PASS','scope':'Lossless residual/utility replay, not training',
        'package_files_verified':count,'arrays_sha256_exact':len(d),'main_models':5,'choice_records':420,
        'paired_comparisons':6,'projection_cross_checks':5,'numeric_checks':checks,
        'numpy':np.__version__,'python':sys.version,'full_prediction_replay':None}
    if a.reference_zip:
        report['full_prediction_replay']=full_replay(a.reference_zip.resolve(),out,d)
    report['seconds']=time.monotonic()-start
    write_json(out/'REPRODUCTION.json',report)
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
