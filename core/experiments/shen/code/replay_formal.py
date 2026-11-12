#!/usr/bin/env python3
"""Refit all three contexts with the locked model, then compare saved arrays."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,subprocess,sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    out=a.out.resolve()
    if out.exists():p.error('Fresh replay directory required')
    out.mkdir(parents=True)
    def run(script,*args):subprocess.run([sys.executable,str(ROOT/'code'/script),*map(str,args)],check=True)
    run('prepare_shen.py','--counts',ROOT/'measurement/canonical_counts.csv','--provenance',ROOT/'measurement/SOURCE_PROVENANCE.json','--out',out/'prepared')
    run('fit_prepared.py','--prepared',out/'prepared','--out',out/'fits','--stage','formal','--max-seconds-per-cell',1800)
    run('score_prepared.py','--prepared',out/'prepared','--fits',out/'fits','--out',out/'scores')
    checks=[]
    for cell in ['HeLa','A549','293T']:
        with np.load(ROOT/f'runs/formal/{cell}/results/high/predictions.npz') as old,np.load(out/f'fits/{cell}/results/high/predictions.npz') as new:
            for model in ['M0','rank1_selected','spectral_selected']:
                diff=float(np.max(np.abs(old[model]-new[model])));checks.append(dict(cell_line=cell,model=model,maximum_difference=diff,pass_tolerance=diff<1e-12))
    result={'checks':checks,'all_passed':all(x['pass_tolerance'] for x in checks),'scope':'Numerical reproducibility, not independent biological evidence'}
    (out/'REPLAY_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['all_passed']:raise RuntimeError('Replay differs beyond tolerance')
    print(json.dumps(result))
if __name__=='__main__':main()
