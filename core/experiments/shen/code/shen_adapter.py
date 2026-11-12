"""Leakage-resistant adapter for a canonical Shen construct/time-series table.

This schema is a normalized interchange contract, NOT a claim about the unseen
publisher workbook headers. Missing/low counts are omitted, never pseudo-counted.
The primary outcome is mean construct-level growth slope, not published GI score.
"""
from __future__ import annotations
import csv, hashlib, itertools, json
from pathlib import Path
from collections import defaultdict
import numpy as np

SCHEMA = ['cell_line','replicate','day','construct_id','gene_a','guide_a','gene_b','guide_b','count']
NTC = 'NTC'
CELLS = ('HeLa','A549','293T')
DAYS = (3.,14.,21.,28.)
SPLIT_SALT = 'WitnessCell-Shen2017-confirmatory-pairs-v1-20261107'
SMOKE_SALT = 'WitnessCell-Shen2017-training-only-smoke-v1-20261107'


def keyed_hash(text: str, salt: str) -> str:
    return hashlib.sha256((salt+'|'+text).encode()).hexdigest()


def read_counts(path):
    """Read canonical CSV, fail rather than guessing workbook columns."""
    with open(path, newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        if list(reader.fieldnames or []) != SCHEMA:
            raise ValueError('Expected exactly the canonical schema: '+','.join(SCHEMA))
        rows=[]
        for row in reader:
            row['replicate']=int(row['replicate']); row['day']=float(row['day'])
            row['count']=float(row['count'])
            rows.append(row)
    return rows


def validate_rows(rows):
    seen=set(); annotations={}; out=[]
    for row in rows:
        if set(row)!=set(SCHEMA): raise ValueError('Unexpected or missing columns')
        cell=row['cell_line']; rep=row['replicate']; day=float(row['day'])
        count=float(row['count']); cid=row['construct_id']
        if cell not in CELLS or isinstance(rep,bool) or rep not in (1,2):
            raise ValueError('Unexpected cell line or replicate')
        if day not in DAYS or not np.isfinite(count) or count<0:
            raise ValueError('Unexpected day or invalid count')
        for key in ('construct_id','gene_a','guide_a','gene_b','guide_b'):
            if not isinstance(row[key],str) or not row[key] or '|' in row[key]:
                raise ValueError('Invalid identity')
        key=(cell,rep,day,cid)
        if key in seen:raise ValueError('Duplicated construct/time record')
        seen.add(key)
        anno=tuple(row[k] for k in ('gene_a','guide_a','gene_b','guide_b'))
        if cid in annotations and annotations[cid]!=anno:
            raise ValueError('Conflicting construct annotation')
        annotations[cid]=anno; out.append(dict(row,day=day,count=count))
    return out


def construct_slopes(rows, min_reads=20, min_points=3, min_control_constructs=3):
    """Use log2(count) minus per-time mean log2 dual-NTC count, then OLS.

    Thresholds are analyst-defined preflight defaults, NOT Shen's unpublished
    per-time empirical thresholds. Source-matched thresholds may replace them
    only in a logged measurement-only amendment before source outcome fitting.
    """
    rows=validate_rows(rows)
    controls=defaultdict(list); groups=defaultdict(list)
    for r in rows:
        key=(r['cell_line'],r['replicate'],r['day'])
        if r['gene_a']==NTC and r['gene_b']==NTC and r['count']>=min_reads:
            controls[key].append(np.log2(r['count']))
        groups[(r['cell_line'],r['replicate'],r['construct_id'])].append(r)
    refs={k:float(np.mean(v)) for k,v in controls.items() if len(v)>=min_control_constructs}
    needed={(r['cell_line'],r['replicate'],r['day']) for r in rows}
    if needed-set(refs):
        raise ValueError('Insufficient dual-NTC reference in one or more time samples')
    slopes=[]; rejected=[]
    for (cell,rep,cid), rs in sorted(groups.items()):
        kept=[r for r in rs if r['count']>=min_reads]
        if len(kept)<min_points:
            rejected.append({'cell_line':cell,'replicate':rep,'construct_id':cid,'reason':'low_timepoint_coverage'});continue
        t=np.array([r['day'] for r in kept]);t=t-t.mean()
        if np.dot(t,t)==0:raise ValueError('No time span')
        z=np.array([np.log2(r['count'])-refs[(cell,rep,r['day'])] for r in kept])
        slope=float(np.dot(t,z-z.mean())/np.dot(t,t))
        residual=z-z.mean()-slope*t
        r=rs[0]
        slopes.append({**{k:r[k] for k in ('cell_line','replicate','construct_id','gene_a','guide_a','gene_b','guide_b')},
                       'fitness':slope,'timepoints':len(kept),'fit_sse':float(np.dot(residual,residual))})
    return slopes,dict(accepted_construct_repeats=len(slopes),rejected=rejected,
                      reference_samples=len(refs),units='control-relative log2 growth slope per day')


def aggregate(slopes, min_constructs=3):
    """Equal physical-construct weighting, with matched available constructs in repeats.

    No rank weighting learned from double outcomes. Same-gene dual targeting is
    excluded. No author whole-map single or interaction estimate is imported.
    """
    tmp=defaultdict(dict)
    for r in slopes:
        a,b=r['gene_a'],r['gene_b']
        if a==NTC and b==NTC:continue
        if a==b:continue
        ga,gb=sorted((a,b)); key=(r['cell_line'],ga,gb)
        k=(r['construct_id'],r['replicate'])
        if k in tmp[key]:raise ValueError('Duplicate slope')
        tmp[key][k]=r['fitness']
    tables={}; exclusions=[]
    for cell in CELLS:
        singles={};pairs={};support={}
        for (c,a,b),values in sorted(tmp.items()):
            if c!=cell:continue
            ids=sorted({cid for cid,rep in values if (cid,1) in values and (cid,2) in values})
            if len(ids)<min_constructs:
                exclusions.append({'cell_line':c,'pair':[a,b],'common_constructs':len(ids)});continue
            v=np.array([[values[(cid,rep)] for rep in (1,2)] for cid in ids])
            mean=v.mean(axis=0)
            if NTC in (a,b):singles[b if a==NTC else a]=mean
            else:pairs[(a,b)]=mean;support[a+'|'+b]=len(ids)
        genes=sorted(singles)
        if not genes:continue
        pairs={p:y for p,y in pairs.items() if p[0] in singles and p[1] in singles}
        tables[cell]={'genes':genes,'singles':np.array([singles[g] for g in genes]),'pairs':pairs,'support':support}
    return tables,exclusions


def split_identities(pairs, test_fraction=.2):
    """Make global physical-pair assignment before response-dependent quality filters."""
    pair_list=sorted(set(tuple(sorted(x)) for x in pairs))
    ordered=sorted(pair_list,key=lambda p:keyed_hash('|'.join(p),SPLIT_SALT))
    ntest=int(np.ceil(len(ordered)*test_fraction))
    query=set(ordered[:ntest]); train=set(ordered[ntest:])
    smoke_order=sorted(train,key=lambda p:keyed_hash('|'.join(p),SMOKE_SALT))
    smoke_val=set(smoke_order[:int(np.ceil(.1*len(smoke_order)))])
    return dict(train=sorted(train),query=sorted(query),smoke_train=sorted(train-smoke_val),smoke_validation=sorted(smoke_val))


def make_input(table, train_ids, query_ids):
    available=set(table['pairs']);tr=sorted(set(map(tuple,train_ids))&available)
    qr=sorted(set(map(tuple,query_ids))&available)
    if set(tr)&set(qr):raise ValueError('Physical-pair leakage')
    return dict(schema='witnesscell-growth-input-v1',genes=table['genes'],
                singles=table['singles'].tolist(),
                training=[{'pair':list(p),'response':table['pairs'][p].tolist()} for p in tr],
                queries=[list(p) for p in qr])


def projection_info(genes, query_pairs):
    lookup={g:i for i,g in enumerate(genes)}
    X=np.zeros((len(query_pairs),len(genes)))
    for k,(a,b) in enumerate(query_pairs):X[k,lookup[a]]=1;X[k,lookup[b]]=1
    U,s,Vt=np.linalg.svd(X,full_matrices=False)
    rank=int(np.sum(s>max(X.shape)*np.finfo(float).eps*(s[0] if len(s) else 1)))
    return X,U[:,:rank],rank


def score(genes,pairs,y,base,pred):
    y,base,pred=map(lambda x:np.asarray(x,float),(y,base,pred))
    if y.shape!=base.shape or y.shape!=pred.shape or y.shape!=(len(pairs),2):
        raise ValueError('Evaluation shapes must agree')
    if not all(np.isfinite(x).all() for x in (y,base,pred)):raise ValueError('Nonfinite scores')
    X,U,rank=projection_info(genes,pairs)
    e=y-base; d=pred-base
    pe=e-U@(U.T@e);pd=d-U@(U.T@d)
    E=float(np.sum(pe**2));A=float(np.sum(pe*pd));B=float(np.sum(pd**2))
    def gain(a,b):return None if a<=1e-24 else 1-b/a
    return dict(query_pairs=len(pairs),endpoint_rank=rank,partner_dimensions=len(pairs)-rank,
                full_sse_M0=float(np.sum(e**2)),full_sse_model=float(np.sum((y-pred)**2)),
                full_gain=gain(float(np.sum(e**2)),float(np.sum((y-pred)**2))),
                partner_sse_M0=E,partner_sse_model=float(np.sum((pe-pd)**2)),
                partner_gain=gain(E,float(np.sum((pe-pd)**2))),
                partner_gain_each_repeat=[gain(float(np.sum(pe[:,r]**2)),float(np.sum((pe[:,r]-pd[:,r])**2))) for r in range(2)],
                full_gain_each_repeat=[gain(float(np.sum(e[:,r]**2)),float(np.sum((y[:,r]-pred[:,r])**2))) for r in range(2)],
                cosine=None if E*B<=1e-24 else A/np.sqrt(E*B),
                oracle_beta=None if B<=1e-24 else A/B,
                E=E,A=A,B=B,identity_error=abs((E-np.sum((pe-pd)**2))-(2*A-B)),
                posthoc_alignment_only=True)
