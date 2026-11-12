#!/usr/bin/env python3
"""Fit the archived functional procedure to user-supplied, outcome-free queries.

The unchanged native engine fits separate repeat baselines and selects rank-one
and spectral residual corrections by three training-only folds. No query truth
is accepted. This adapter adds validation and file I/O, not a new estimator.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '1'
from pathlib import Path
import argparse, csv, hashlib, json, sys, time
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'growth_vendor'))
import run_smoke as native

PLAN = dict(tau_fractions=[.3, .1, .03, .01], nu=[.1, 1.],
            optimizer_maxiter=3000, optimizer_tolerance=1e-6)

def validate(data):
    allowed = {'schema', 'genes', 'singles', 'training', 'queries'}
    if set(data) != allowed or data.get('schema') != 'witnesscell-growth-input-v1':
        raise ValueError('Expected schema, genes, singles, training and queries only')
    genes = data['genes']
    if (not isinstance(genes, list) or len(genes) < 3 or
            any(not isinstance(g, str) or not g or '|' in g for g in genes) or
            len(set(genes)) != len(genes)):
        raise ValueError('genes must be unique nonempty identities without |')
    singles = np.asarray(data['singles'], float)
    if singles.shape != (len(genes), 2) or not np.isfinite(singles).all():
        raise ValueError('singles must be finite, genes-by-two-repeat measurements')
    lookup = {g:i for i,g in enumerate(genes)}
    used = set()
    def edge(pair):
        if not isinstance(pair, list) or len(pair) != 2 or pair[0] == pair[1]:
            raise ValueError('Each pair must contain two different gene identities')
        if any(g not in lookup for g in pair):
            raise ValueError('Unknown endpoint')
        key = tuple(sorted(pair))
        if key in used:
            raise ValueError('Duplicate or training-query overlapping physical pair')
        used.add(key)
        return [lookup[g] for g in key], '|'.join(key)
    train_edges, train_ids, values = [], [], []
    if not isinstance(data['training'], list) or not data['training']:
        raise ValueError('training must be a nonempty list')
    for item in data['training']:
        if not isinstance(item, dict) or set(item) != {'pair', 'response'}:
            raise ValueError('A training record contains pair and response only')
        ed, ident = edge(item['pair'])
        y = np.asarray(item['response'], float)
        if y.shape != (2,) or not np.isfinite(y).all():
            raise ValueError('Each training response must have two finite repeats')
        train_edges.append(ed); train_ids.append(ident); values.append(y)
    if not isinstance(data['queries'], list) or not data['queries']:
        raise ValueError('queries must be a nonempty list of identity pairs')
    query_edges, query_ids = [], []
    for pair in data['queries']:
        ed, ident = edge(pair)
        query_edges.append(ed); query_ids.append(ident)
    folds = np.array([int.from_bytes(hashlib.sha256(
        ('wc_spectral_smoke_v1|' + s).encode()).digest()[:8], 'big') % 3
        for s in train_ids], dtype=np.int64)
    if set(folds.tolist()) != {0,1,2}:
        raise ValueError('Training pairs must populate all three identity-hash folds')
    training = dict(genes=np.array(genes), singles=singles,
        edges=np.array(train_edges,dtype=np.int64), y=np.array(values),
        global_ids=np.array(train_ids), inner_fold=folds,
        low_positions=np.arange(len(train_edges)))
    return training, np.array(query_edges,dtype=np.int64), np.array(query_ids)

def predict(data, output, max_seconds=600):
    training, edges, ids = validate(data)
    if not np.isfinite(max_seconds) or max_seconds <= 0:
        raise ValueError('max-seconds must be finite and positive')
    output = Path(output).resolve()
    if output.exists():
        raise ValueError('Use a fresh output directory; existing results are preserved')
    output.mkdir(parents=True)
    native.ROOT = output
    begun = time.monotonic()
    native.run_budget('high', training, edges, ids, PLAN, begun + max_seconds)
    with np.load(output/'results/high/predictions.npz', allow_pickle=False) as z:
        with (output/'predictions.csv').open('w',newline='') as stream:
            writer=csv.writer(stream)
            writer.writerow(['pair','repeat','baseline','rank1_prediction',
                             'rank1_increment','spectral_prediction'])
            for i, ident in enumerate(ids):
                for j in range(2):
                    writer.writerow([ident,j+1,float(z['M0'][i,j]),
                        float(z['rank1_selected'][i,j]),
                        float(z['rank1_selected'][i,j]-z['M0'][i,j]),
                        float(z['spectral_selected'][i,j])])
    metadata=dict(configuration='functional/frozen_rank1_primary',
        source='unchanged Jurkat frozen engine from grn_repo_v27',
        input_sha256=hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest(),
        query_outcomes_supplied=False, training_pairs=len(training['edges']),
        query_pairs=len(ids), readouts=1, repeats=2,
        units='Inherited original functional response units',
        primary='rank1_selected', secondary='spectral_selected',
        seconds=time.monotonic()-begun,
        scope='New-input fitting interface; supplied real-data replay is not new validation',
        selections=json.loads((output/'results/high/selection.json').read_text()))
    (output/'fit_metadata.json').write_text(json.dumps(native.clean(metadata),indent=2)+'\n')
    return metadata

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--max-seconds',type=float,default=600)
    args=ap.parse_args()
    try:
        metadata=predict(json.loads(args.input.read_text()),args.output,args.max_seconds)
    except (ValueError,KeyError,TypeError) as error:
        ap.error(str(error))
    print(json.dumps(metadata))
if __name__=='__main__':
    main()
