#!/usr/bin/env python3
"""Redraw 7 main and 10 ED figures from preserved publication intermediates."""
from pathlib import Path
import argparse
import csv
import math
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parent

def check_original_extension():
    """Verify ED9/S11/S12 arithmetic from original six-split summary tables."""
    folder=ROOT/'source_data/original_system_extension'
    with (folder/'SPLIT_RESULTS.csv').open(newline='') as f: splits=list(csv.DictReader(f))
    with (folder/'DATASET_RESULTS.csv').open(newline='') as f: datasets=list(csv.DictReader(f))
    baseline={(r['dataset'],r['seed']):r for r in splits if r['method']=='zero'}
    discrepancies=[]
    for r in splits:
        base=baseline[r['dataset'],r['seed']]
        for loss,gain in [('non_target_mse','non_target_gain_vs_original'),('endpoint_removed_mse','endpoint_removed_gain_vs_original')]:
            calculated=1-float(r[loss])/float(base[loss])
            discrepancies.append(abs(calculated-float(r[gain])))
    for r in datasets:
        subset=[x for x in splits if x['dataset']==r['dataset'] and x['method']==r['method']]
        if len(subset)!=3:raise ValueError('Expected exactly three historical seeds')
        if sum(int(x['test_pairs']) for x in subset)!=int(r['pair_seed_records']):raise ValueError('Pair-seed count mismatch')
        for loss,gain in [('non_target_mse','non_target_geometric_gain'),('endpoint_removed_mse','endpoint_removed_geometric_gain')]:
            ratios=[float(x[loss])/float(baseline[x['dataset'],x['seed']][loss]) for x in subset]
            calculated=1-math.exp(sum(map(math.log,ratios))/len(ratios))
            discrepancies.append(abs(calculated-float(r[gain])))
    maximum=max(discrepancies)
    if maximum>1e-12:raise ValueError('Original extension aggregate discrepancy '+str(maximum))
    return {'status':'PASS','split_method_rows':len(splits),'dataset_method_rows':len(datasets),'maximum_absolute_difference':maximum,'scope':'Arithmetic reconstruction from archived scalar losses; not prediction-vector or model-fit replay'}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,type=Path)
    args=ap.parse_args();out=args.out.resolve()
    if out.exists() or out.is_relative_to(ROOT):ap.error('--out must be new and outside this module')
    identity=json.loads((ROOT/'SOURCE_IDENTITIES.json').read_text())
    for name,row in identity['source_files'].items():
        p=ROOT/name
        if p.stat().st_size!=row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('Source changed: '+name)
    out.mkdir(parents=True)
    for name in ('source_data','figure_code'):shutil.copytree(ROOT/name,out/name)
    env=os.environ.copy();env.pop('PYTHONPATH',None)
    env.update(PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    subprocess.run([sys.executable,'-B',str(out/'figure_code/build_figures.py')],cwd=out,env=env,check=True)
    names=['Fig'+str(i) for i in range(1,8)]+['ED'+str(i) for i in range(1,11)]
    for name in names:
        for ext in ('pdf','png','svg'):
            p=out/'figures'/(name+'.'+ext)
            if not p.is_file() or p.stat().st_size==0:raise ValueError('Missing figure '+str(p))
    table_ids=re.findall(r'^### Supplementary Table (S\d+)\s', (ROOT/'PUBLISHED_TABLES_v22.md').read_text(),flags=re.M)
    if table_ids!=['S'+str(i) for i in range(1,20)]:raise ValueError('Table reference incomplete')
    validation=json.loads((out/'figure_code/FIGURE_VALIDATION.json').read_text())
    result={'status':'PASS','figures':names,'rendered_files':51,'published_table_references':table_ids,'figure_builder_source_members_verified':len(identity['source_files']),'original_builders_byte_identical':True,'source_checks':validation,'original_extension_aggregate_check':check_original_extension(),'scope':'Exact-source numerical redraw and published-table reference; tables are not represented as new raw-cell fits or independent validation'}
    (out/'PAPER_REPRODUCTION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','rendered_files','published_table_references','scope']},indent=2))

if __name__=='__main__':main()
