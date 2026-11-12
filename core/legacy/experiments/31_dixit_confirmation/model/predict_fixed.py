#!/usr/bin/env python3
"""Identity-only fixed-reference inference from one saved fit.

Example:
 python predict_fixed.py --state run_compact/FITTED_STATE_0.npz --lambda 0.1 \
     --pairs GENE_A+GENE_B GENE_A+GENE_C

Use the lambda recorded under choices.global in the corresponding FOLD JSON.
This command receives no observed query response or query-batch reference.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,sys
from pathlib import Path
import numpy as np

def predict_from_state(path,pairs,lam):
    with np.load(path,allow_pickle=False) as archive:s={k:archive[k] for k in archive.files}
    ni={str(g):i for i,g in enumerate(s['nodes'])};edges=np.array([[ni[g] for g in pair.split('+')] for pair in pairs],int)
    if np.any(edges[:,0]==edges[:,1]):raise ValueError('Distinct measured-single endpoints required')
    x=np.zeros((len(pairs),len(ni)));x[np.arange(len(pairs)),edges[:,0]]=1.;x[np.arange(len(pairs)),edges[:,1]]=1.
    additive=np.mean((x@s['e'])[:,s['pi']]/s['sigma'][s['pi']],axis=1,keepdims=True)
    base=additive/(1+float(s['alpha'])*abs(additive));ga,gb=edges.T;ta,tb=s['et'].T;gram=s['gram']
    kq=.5*(gram[np.ix_(ga,ta)]*gram[np.ix_(gb,tb)]+gram[np.ix_(ga,tb)]*gram[np.ix_(gb,ta)])
    raw=((kq-x@s['plus']@s['ktt'])@s['pt']/float(s['scale']))@s[f'c_{float(lam)}'] if float(s['scale']) else np.zeros_like(base)
    correction=raw-x@s[f'beta_global_{float(lam)}']
    return {'prediction':(base+correction).ravel(),'baseline':base.ravel(),'partner_correction':correction.ravel()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--state',type=Path,required=True);p.add_argument('--lambda',dest='lam',type=float,choices=(.1,1.,10.),required=True)
    p.add_argument('--pairs',nargs='+',required=True);a=p.parse_args();result=predict_from_state(a.state,a.pairs,a.lam)
    writer=csv.DictWriter(sys.stdout,fieldnames=['pair',*result]);writer.writeheader()
    for i,pair in enumerate(a.pairs):writer.writerow({'pair':pair,**{k:float(v[i]) for k,v in result.items()}})

if __name__=='__main__':main()
