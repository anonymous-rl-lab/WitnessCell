#!/usr/bin/env python3
"""Compact Dixit G2M fixed-reference inference for one or more pair identities.

These states are frozen outer-fold evaluation models, not confirmed production
models. The endpoint domain, program scale and reference are fixed in the state.
No labels, full expression matrices, or training routines are read at inference.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,sys
from pathlib import Path
import numpy as np


def predict(state_path,pairs):
    """Return predictions on a fixed fitted state; pairs may be singleton lists."""
    with np.load(state_path,allow_pickle=False) as archive:
        state={key:archive[key] for key in archive.files}
    if int(state['format_version'])!=1:raise ValueError('Unsupported compact state version')
    nodes={str(node):i for i,node in enumerate(state['nodes'])}
    edges=[]
    for pair in pairs:
        genes=pair.split('+')
        if len(genes)!=2 or genes[0]==genes[1]:raise ValueError('Use two distinct endpoints, e.g. CREB1+YY1')
        if any(gene not in nodes for gene in genes):raise ValueError('Endpoint outside this fitted state: '+pair)
        edges.append([nodes[gene] for gene in genes])
    if not edges:return {name:np.empty(0) for name in ('prediction','baseline','partner_correction')}
    edges=np.asarray(edges,int);x=np.zeros((len(edges),len(nodes)))
    x[np.arange(len(edges)),edges[:,0]]=1.;x[np.arange(len(edges)),edges[:,1]]=1.
    score=np.mean((x@state['program_endpoint_effects'])/state['program_sigma'],axis=1,keepdims=True)
    baseline=score/(1.+float(state['alpha'])*abs(score))
    a,b=edges.T;u,v=state['training_edges'].T;gram=state['endpoint_gram']
    kernel=.5*(gram[np.ix_(a,u)]*gram[np.ix_(b,v)]+gram[np.ix_(a,v)]*gram[np.ix_(b,u)])
    correction=kernel@state['kernel_coefficients']-x@state['endpoint_offset']
    return {'prediction':(baseline+correction).ravel(),'baseline':baseline.ravel(),'partner_correction':correction.ravel()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state',type=Path,required=True)
    parser.add_argument('--pairs',nargs='+',required=True)
    args=parser.parse_args();result=predict(args.state,args.pairs)
    writer=csv.DictWriter(sys.stdout,fieldnames=['pair',*result]);writer.writeheader()
    for i,pair in enumerate(args.pairs):writer.writerow({'pair':pair,**{key:float(value[i]) for key,value in result.items()}})


if __name__=='__main__':main()
