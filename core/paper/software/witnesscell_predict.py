#!/usr/bin/env python3
"""New-input adapter for the unchanged released WitnessCell core.

Condition means and variances must share the provided normalized gene axis.
This adapter does not implement the RNA R extension, PairPick, or growth model.
"""
from pathlib import Path
import argparse, csv, hashlib, json, sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor'))
from witnesscell_v4_core import fit_predict

def canonical(s):
    if not isinstance(s,str) or not s or any(c.isspace() for c in s):
        raise ValueError('Condition identities must be nonempty strings without whitespace')
    parts=s.split('+')
    if len(parts)>2 or len(set(parts))!=len(parts):
        raise ValueError('Only distinct single or double gene interventions are supported')
    return '+'.join(sorted(parts))

def validate(data):
    genes=data.get('genes',[])
    if len(genes)<2 or len(set(genes))!=len(genes) or any(not isinstance(g,str) or '+' in g for g in genes):
        raise ValueError('genes must be a unique ordered gene axis')
    observed=data.get('observations',{})
    train=data.get('train_conditions',[]); val=data.get('validation_conditions',[]);query=data.get('query_conditions',[])
    if not train or not val or not query: raise ValueError('Explicit nonempty train, validation, and query lists are required')
    if any(not isinstance(q,str) for q in query):raise ValueError('Queries contain identities only; query outcomes are not accepted')
    seen={}
    for group,items in [('train',train),('validation',val),('query',query)]:
        for c in items:
            key=canonical(c)
            if key in seen:raise ValueError('Duplicate or leaking intervention identity: '+key)
            seen[key]=group
            if c=='control' or any(n not in genes for n in c.split('+')):raise ValueError('Interventions must be genes on the gene axis; control is implicit')
    if set(observed)!=set(train+val+['control']):raise ValueError('Observations must contain exactly control, train and validation; never query outcomes')
    if not any('+' in c for c in train) or not any('+' in c for c in val):raise ValueError('This interface requires training and validation doubles')
    if not any('+' not in c for c in train):raise ValueError('At least one measured training single is required')
    for c,row in observed.items():
        for key in ['mean','variance']:
            a=np.asarray(row.get(key),dtype=float)
            if a.shape!=(len(genes),) or not np.isfinite(a).all():raise ValueError(c+': invalid '+key)
            if key=='variance' and (a<0).any():raise ValueError(c+': negative variance')
        if not isinstance(row.get('n'),int) or isinstance(row['n'],bool) or row['n']<2:raise ValueError(c+': n must be an integer >= 2')
    go=data.get('gene2go')
    if not isinstance(go,dict) or any(not isinstance(v,list) or any(not isinstance(t,str) for t in v) for v in go.values()):
        raise ValueError('gene2go must map gene identities to lists of annotation strings')
    return data

def predict(data):
    validate(data); o=data['observations']
    fit=fit_predict(means={k:np.asarray(v['mean'],float) for k,v in o.items()},
        variances={k:np.asarray(v['variance'],float) for k,v in o.items()},counts={k:v['n'] for k,v in o.items()},
        genes=data['genes'],train_conditions=data['train_conditions'],validation_conditions=data['validation_conditions'],
        target_conditions=data['query_conditions'],identity_mode='multihead_v14',gene2go=data['gene2go'])
    metadata={'configuration':'released_core/multihead_v14','alpha':fit.alpha,'noise_ratio':fit.noise_ratio,'gamma':fit.gamma,
        'known_single_genes':list(fit.known_single_genes),'training_doubles':list(fit.training_doubles),
        'validation_doubles':list(fit.validation_doubles),'input_sha256':hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest(),
        'normalization':'Inherited from user input; no raw-count normalization performed',
        'output_scope':'Predicted normalized condition means and nonlinear baseline; no calibrated uncertainty or partner-extension score'}
    return fit,metadata

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    try: fit,meta=predict(json.loads(a.input.read_text()))
    except (ValueError,KeyError,TypeError) as e:ap.error(str(e))
    data=json.loads(a.input.read_text());a.output.mkdir(parents=True,exist_ok=True)
    with (a.output/'predictions.csv').open('w',newline='') as h:
        w=csv.writer(h);w.writerow(['condition','gene','prediction','nonlinear_baseline'])
        for c in data['query_conditions']:
            for j,g in enumerate(data['genes']):w.writerow([c,g,float(fit.predictions[c][j]),float(fit.factorized_predictions[c][j])])
    (a.output/'fit_metadata.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'queries':len(data['query_conditions']),'readouts':len(data['genes']),'output':str(a.output)}))
if __name__=='__main__':main()
