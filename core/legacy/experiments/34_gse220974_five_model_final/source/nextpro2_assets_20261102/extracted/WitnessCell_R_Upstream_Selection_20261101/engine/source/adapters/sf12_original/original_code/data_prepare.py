"""Freeze equal-budget support exchange and opposite-half evaluation tasks.

Only the 22 already-used FIT combinations supply responses.  Historical test
and the former eight validation responses are never opened.  Training files
contain one half only; opposite-half query responses live in separate files.
The common singles, control reference and four ABC measurements are permitted
shared inputs.  Accordingly this is a development experiment conditional on
those measurements, not independent biological replication.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(key,'2')
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np

R=Path(__file__).resolve().parent
ROOT=R.parent
PACKAGE=R
COMPACT_SOURCE_OVERRIDE=None
sys.path.insert(0,str(ROOT/'dynamic_factorial_run'))
from baselines_metrics import _m0_tune, _m0_state, predict_state, endpoint_basis

TASKS=tuple(f'{layout}_h{half}' for layout in ('low','high','seen') for half in (0,1))
QUERY=('CREB1+EGR1','CREB1+GABPA','EGR1+ELK1','ELK1+GABPA')
REMOVED={'low':('CREB1+E2F4','CREB1+ETS1'),
         'high':('ELF1+YY1','IRF1+NR2C2')}
FIT_SOURCE_KEYS=('fit_pair_ids','fit_edges','fit_conditions','fit_counts',
                 'fit_half_counts','fit_half0_effects','fit_half1_effects',
                 'fit_full_effects','reference_full')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def save_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n');temp.replace(path)
def save_npz(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.atomic_tmp.npz')
    with temp.open('wb') as handle:
        np.savez_compressed(handle,**arrays)
        handle.flush();os.fsync(handle.fileno())
    with np.load(temp,allow_pickle=False) as z:
        assert set(z.files)==set(arrays)
        for k,a in arrays.items():
            if not np.array_equal(z[k],np.asarray(a)):raise RuntimeError('roundtrip '+k)
    temp.replace(path)

def resolve_source(compact_source=None):
    """An explicit source wins; packaged compact source is the default."""
    if compact_source is not None:return Path(compact_source).resolve()
    if COMPACT_SOURCE_OVERRIDE is not None:return Path(COMPACT_SOURCE_OVERRIDE).resolve()
    compact=PACKAGE/'source/compact_source.npz'
    return compact if compact.exists() else ROOT/'density_diagnosis_run/probe_inputs.npz'

def source_provenance(path=None):
    path=resolve_source(path)
    with np.load(path,allow_pickle=False) as z:
        original_sha=str(z['original_source_sha256'].item()) if 'original_source_sha256' in z.files else sha(path)
        original_path=str(z['original_source_path'].item()) if 'original_source_path' in z.files else 'density_diagnosis_run/probe_inputs.npz'
    return dict(used_path=str(path),used_sha256=sha(path),original_source_path=original_path,
                original_source_sha256=original_sha,compact_source=path.name=='compact_source.npz')

def export_compact_source():
    """Copy allowed22FIT sufficient statistics; no validation/test responses."""
    original=ROOT/'density_diagnosis_run/probe_inputs.npz'
    with np.load(original,allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in FIT_SOURCE_KEYS}
        # Capture data are metadata only, restricted to FIT identities+control.
        names=z['capture_names'].copy();conditions=z['capture_conditions']
        capture=z['capture_half_counts']
        wanted=arrays['fit_conditions'].tolist()+['control_2guide']
        idx=[conditions.tolist().index(c) for c in wanted]
        arrays.update(capture_names=names,capture_conditions=conditions[idx],capture_half_counts=capture[idx])
    arrays.update(original_source_path=np.array('density_diagnosis_run/probe_inputs.npz'),
                  original_source_sha256=np.array(sha(original)))
    target=PACKAGE/'source/compact_source.npz'
    save_npz(target,**arrays)
    with np.load(original,allow_pickle=False) as z,np.load(target,allow_pickle=False) as compact:
        assert all(np.array_equal(z[k],compact[k]) for k in FIT_SOURCE_KEYS)
        assert not any(k.startswith('val_') or k.startswith('pair_raw_') for k in compact.files)
    provenance=source_provenance(target)
    provenance.update(bytes=target.stat().st_size,keys=list(arrays),fit_fields_exact=True,
         response_scope='Only22originalFIT combinations; no former validation or historical test responses.',
         capture_metadata_scope='22FIT identities and control2 only.')
    save_json(PACKAGE/'source/COMPACT_SOURCE_PROVENANCE.json',provenance)
    return provenance

def state_from_task(common,train):
    return dict(single_effects=common['single_effects'],coef=train['m0_coef'],
                alpha=float(train['m0_alpha']),gamma=float(train['m0_gamma']))

def baseline(common,train,u):
    u=np.asarray(u,dtype=np.float64)
    pred=predict_state(state_from_task(common,train),u)
    one=np.isclose(u.sum(1),1);zero=np.isclose(u.sum(1),0)
    pred[one]=u[one]@common['single_effects'];pred[zero]=0
    return pred

def load_task(name,include_evaluation=False):
    """Default returns (common, train); evaluation must be requested explicitly.

    common has only allowed single/control/ABC information and identities.
    train contains y, u, pair_ids, half_counts, M0 state and loss normalizers.
    evaluation has y/u for the fixed four queries from the opposite half.
    """
    if name not in TASKS:raise ValueError(name)
    common=read_npz(R/'data/common.npz')
    train=read_npz(R/f'data/{name}_train.npz')
    if include_evaluation:
        return common,train,read_npz(R/f'data/{name}_eval.npz')
    return common,train

def capture_audit():
    """Metadata-only comparison of capture allocation in the matched layouts."""
    protocol=json.loads((R/'DATA_PROTOCOL.json').read_text())
    with np.load(resolve_source(),allow_pickle=False) as z:
        names=z['capture_names'].tolist();conditions=z['capture_conditions'].tolist()
        counts=z['capture_half_counts']
    sums={}
    for layout in ('low','high','seen'):
        idx=[conditions.index(c) for c in protocol['layouts'][layout]['pair_conditions']]
        sums[layout]=counts[idx].sum(0)
    result=dict(capture_names=names,layout_capture_half_counts={k:v.tolist() for k,v in sums.items()},
        half_direction_total_variation_low_high=[],
        interpretation='Same total cells and direct AB edge count; capture composition is not exactly matched. Training-edge replacement changes biology and measurement distribution jointly.')
    for half in (0,1):
        a=sums['low'][:,half]/sums['low'][:,half].sum()
        b=sums['high'][:,half]/sums['high'][:,half].sum()
        result['half_direction_total_variation_low_high'].append(float(np.abs(a-b).sum()/2))
    save_json(R/'CAPTURE_MATCH_AUDIT.json',result)
    return result

def add_seen_all_evaluation():
    """Add a secondary same-condition test across all22FIT combinations.

    Primary four-query files and every training file remain byte-identical.
    This helper opens only named FIT half arrays and the fixed shared control.
    """
    common=read_npz(R/'data/common.npz')
    frozen={p.name:sha(p) for p in sorted((R/'data').glob('*_train.npz'))}
    frozen.update({p.name:sha(p) for p in sorted((R/'data').glob('*_eval.npz'))})
    source=resolve_source()
    with np.load(source,allow_pickle=False) as z:
        allowed=('fit_half0_effects','fit_half1_effects','fit_half_counts',
                 'fit_conditions','fit_pair_ids','reference_full')
        h={k:z[k].copy() for k in allowed}
    if not np.array_equal(h['fit_pair_ids'],common['fit_pair_ids']):
        raise RuntimeError('All-seen identity mismatch')
    if not np.array_equal(h['fit_conditions'],common['fit_conditions']):
        raise RuntimeError('All-seen condition mismatch')
    reference_shift=h['reference_full']-common['reference']
    if np.max(np.abs(reference_shift))>1e-12:raise RuntimeError('All-seen reference mismatch')
    x=common['fit_u'].astype(float);contrast=endpoint_basis(x)
    for half in (0,1):
        evaluation_half=1-half
        save_npz(R/f'data/seen_h{half}_eval_all.npz',
            y=h[f'fit_half{evaluation_half}_effects']+reference_shift,
            u=x,counts=h['fit_half_counts'][:,evaluation_half],
            conditions=h['fit_conditions'],pair_ids=h['fit_pair_ids'],
            half=np.array(evaluation_half),contrast=contrast)
    if any(sha(R/'data'/name)!=digest for name,digest in frozen.items()):
        raise RuntimeError('Existing primary or training data changed')
    protocol=json.loads((R/'DATA_PROTOCOL.json').read_text())
    protocol['secondary_seen_all_evaluation']={
        'purpose':'Same22combinations cross-half reproducibility diagnostic; not novel-combination prediction.',
        'files':['seen_h0_eval_all.npz','seen_h1_eval_all.npz'],
        'condition_count':int(len(x)),
        'partner_contrast_dimension':int(len(contrast)),
        'direction':'Each frozen seen_h model evaluates all22FIT combinations in opposite half.',
        'primary_query':'The original four-query endpoint stays unchanged.',
        'data_changes':'No existing training or primary-evaluation bytes changed.',
        'limitations':'Technical cell split within one experiment; shared singles/control/ABC; not independent biology.',
    }
    save_json(R/'DATA_PROTOCOL.json',protocol)
    manifest={p.name:dict(sha256=sha(p),bytes=p.stat().st_size)
              for p in sorted((R/'data').iterdir()) if p.is_file()}
    save_json(R/'DATA_MANIFEST.json',manifest)
    audit=dict(secondary=protocol['secondary_seen_all_evaluation'],
        preserved_file_sha256=frozen,source_sha256=sha(source),
        control_cancellation_max_abs=float(np.max(np.abs(contrast@np.ones(len(x))))),
        new_sha256={f'seen_h{half}_eval_all.npz':sha(R/f'data/seen_h{half}_eval_all.npz') for half in (0,1)})
    save_json(R/'SEEN_ALL_DATA_AUDIT.json',audit)
    return audit

def prepare():
    common_source=ROOT/'dynamic_factorial_run/data/common.npz'
    triple_source=ROOT/'dynamic_factorial_run/data/triples.npz'
    half_source=resolve_source()
    a=read_npz(common_source);triple=read_npz(triple_source)
    # Explicit key allowlist: do not materialize val_* or mixed 30-row arrays.
    allowed=FIT_SOURCE_KEYS
    with np.load(half_source,allow_pickle=False) as z:
        h={k:z[k].copy() for k in allowed}
    ids=h['fit_pair_ids'];labels=h['fit_conditions'].astype(str)
    if len(ids)!=22 or not np.array_equal(a['pairs'][ids],labels):raise RuntimeError('fit identity mismatch')
    if not np.array_equal(a['edges'][ids],h['fit_edges']):raise RuntimeError('edge mismatch')
    reference_error=float(np.max(np.abs(a['reference']-h['reference_full'])))
    if reference_error>1e-12:raise RuntimeError('reference mismatch')
    e=a['single_effects'];floor=max(float(np.quantile(np.sqrt(np.mean(e**2,axis=0)),.5)),.01)
    scale=np.sqrt(np.mean(e**2,axis=0)+floor**2)
    desc=np.column_stack([a['reference']/(np.std(a['reference'])+1e-8),(e/scale).T]).astype('float32')
    scale_error=float(np.max(np.abs(scale-a['scale'])))
    desc_error=float(np.max(np.abs(desc-a['descriptor'])))
    if scale_error!=0 or desc_error!=0:raise RuntimeError('single-only descriptors do not replay')
    # Two archive accumulation orders differ by <= floating point roundoff.
    # Put every half in the exact common.npz reference coordinate explicitly.
    reference_shift=h['reference_full']-a['reference']
    yhalves=[h['fit_half0_effects']+reference_shift,h['fit_half1_effects']+reference_shift]
    counts=h['fit_counts'];hn=h['fit_half_counts']
    recon=(yhalves[0]*hn[:,0,None]+yhalves[1]*hn[:,1,None])/counts[:,None]
    recon_error=float(np.max(np.abs(recon-h['fit_full_effects'])))
    if recon_error>1e-12:raise RuntimeError('half reconstruction mismatch')
    q=np.array([np.flatnonzero(labels==label).item() for label in QUERY],dtype=np.int64)
    masks={layout:np.array([label not in set(QUERY)|set(REMOVED[layout]) for label in labels])
           for layout in ('low','high')}
    masks['seen']=np.ones(len(ids),dtype=bool)
    qx=a['incidence'][ids[q]].astype(float);cq=endpoint_basis(qx)
    if len(cq)!=1:raise RuntimeError('Expected one query contrast')
    shared_keys=('reference','single_effects','scale','descriptor','targets','non_target',
                 'nodes','genes','gene_symbols')
    common={k:a[k] for k in shared_keys}
    common.update(triple_conditions=triple['conditions'],triple_u=triple['u'],
                  triple_y=triple['effects'],triple_counts=triple['counts'],
                  query_conditions=labels[q],query_pair_ids=ids[q],query_u=qx,
                  fit_conditions=labels,fit_pair_ids=ids,fit_u=a['incidence'][ids])
    save_npz(R/'data/common.npz',**common)
    metadata={}
    for layout,mask in masks.items():
        x=a['incidence'][ids[mask]].astype(float)
        deg=x.sum(0).astype(int)
        metadata[layout]=dict(pair_conditions=labels[mask].tolist(),pair_ids=ids[mask].tolist(),
           pair_count=int(mask.sum()),full_cell_count=int(counts[mask].sum()),
           half_cell_counts=hn[mask].sum(0).tolist(),
           endpoint_rank=int(np.linalg.matrix_rank(x)),partner_contrast_dimension=int(len(endpoint_basis(x))),
           degrees=dict(zip(a['nodes'].tolist(),deg.tolist())),
           query_overlap=labels[np.flatnonzero(mask&np.isin(np.arange(22),q))].tolist())
    assert metadata['low']['pair_count']==metadata['high']['pair_count']==16
    assert metadata['low']['full_cell_count']==metadata['high']['full_cell_count']==2325
    assert metadata['low']['half_cell_counts']==metadata['high']['half_cell_counts']==[1163,1162]
    assert metadata['low']['endpoint_rank']==metadata['high']['endpoint_rank']==10
    triple_pairs=[]
    for name,u in zip(triple['conditions'],triple['u']):
        members=set(a['nodes'][u>0].tolist())
        triple_pairs.append(dict(condition=str(name),query_pairs=[p for p in QUERY if set(p.split('+'))<=members],
           exchange_pairs=[p for p in REMOVED['low']+REMOVED['high'] if set(p.split('+'))<=members]))
    protocol=dict(protocol_version=1,seed=20261027,
      source='Dixit2016 K562 TF Cas9 knockout; existing development fit pool',
      query_conditions=labels[q].tolist(),query_counts=counts[q].tolist(),query_contrast_dimension=1,
      tasks=list(TASKS),layouts=metadata,exchange_rules=REMOVED,
      task_rules={'low_high':'Exclude every query AB response from training; use only selected other AB half.',
                  'seen':'All22 FIT AB half visible, including query AB; score their opposite half. Diagnostic, not novel-combination prediction.',
                  'half_direction':'Train half h; evaluate half 1-h. Never select checkpoints from evaluation scores.'},
      shared_information={'single':'All10 measured singles, permitted known-single information, same in every task.',
            'control':'One fixed full control2 reference shared by all. Conditional cross-half test; shared-control uncertainty cancels from C Y since C1=0.',
            'descriptor':'Exact original11dim descriptor = normalized control2 +10 scaled measured single effects; no combination labels.',
            'scale':'Per-RNA scale from singles alone, unchanged; all combination-dependent denominators recomputed from task training half.',
            'anchor':'Original200epoch single-only anchor may be reused; all combination-finetuned weights forbidden as starts.',
            'triples':'All4 original measured ABC (49cells) shared unchanged at auxiliary weight0.1; they contain related query/exchange subpairs.'},
      triple_subpair_disclosure=triple_pairs,
      baseline='Original M0 saturation and endpoint shrinkage, refit/tuned on each task training half only with identity-based internal folds.',
      subset_isolation='Exact query AB measurements hidden in low/high. Related higher-order ABC remain allowed, so this tests value of additional direct AB evidence conditional on fixed ABC.',
      interpretation='Existing development panel; one query contrast; two cell halves are not biological replicates. Equal global cells do not imply equal library/state composition.',
      forbidden_data=['dynamic_factorial_run/data/test.npz','dynamic_factorial_run/data/validation.npz','probe_inputs val_* response keys','mixed pair_raw_* response arrays'],
      source_sha256={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):sha(p)
                     for p in (common_source,triple_source,half_source)},
      half_source_provenance=source_provenance(half_source),
      audit=dict(single_scale_reconstruction_max_abs=scale_error,descriptor_reconstruction_max_abs=desc_error,
                 source_reference_roundoff_max_abs=reference_error,
                 half_weighted_reconstruction_max_abs=recon_error,
                 common_control_partner_cancellation_max_abs=float(np.max(np.abs(cq@np.ones(4))))))
    save_json(R/'DATA_PROTOCOL.json',protocol)
    for layout in ('low','high','seen'):
        mask=masks[layout];x=a['incidence'][ids[mask]].astype(float)
        local_ids=np.arange(mask.sum(),dtype=np.int64);cfit=endpoint_basis(x)
        for half in (0,1):
            name=f'{layout}_h{half}';y=yhalves[half][mask].copy()
            pars=_m0_tune(e,y,x,local_ids,a['non_target'])
            state=_m0_state(e,y,x,local_ids,pars)
            den=np.array([max(float(np.mean(e**2)),1e-8),max(float(np.mean(y**2)),1e-8),
                          max(float(np.mean(triple['effects']**2)),1e-8)])
            pden=np.array(max(float(np.mean((cfit@y)**2)),1e-8))
            save_npz(R/f'data/{name}_train.npz',y=y,u=x,pair_ids=ids[mask],conditions=labels[mask],
                     counts=hn[mask,half],full_counts=counts[mask],half=np.array(half),
                     contrast=cfit,denominators=den,partner_denominator=pden,
                     m0_coef=state['coef'],m0_alpha=np.array(state['alpha']),m0_gamma=np.array(state['gamma']),
                     m0_noise=np.array(pars['noise']))
            save_npz(R/f'data/{name}_eval.npz',y=yhalves[1-half][q],u=qx,pair_ids=ids[q],
                     conditions=labels[q],counts=hn[q,1-half],half=np.array(1-half),contrast=cq)
            save_json(R/f'data/{name}_M0_SELECTION.json',pars)
            print(json.dumps(dict(task=name,pairs=int(mask.sum()),train_half_cells=int(hn[mask,half].sum()),
                 m0_alpha=pars['alpha'],m0_noise=pars['noise'],m0_gamma=pars['gamma']),ensure_ascii=False),flush=True)
    manifest={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted((R/'data').iterdir()) if p.is_file()}
    save_json(R/'DATA_MANIFEST.json',manifest)
    capture_audit()
    print(json.dumps(dict(total_bytes=sum(v['bytes'] for v in manifest.values()),tasks=TASKS,
                         audits=protocol['audit']),ensure_ascii=False),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','capture-audit','seen-all-eval','export-compact'],default='prepare',nargs='?')
    ap.add_argument('--compact-source',type=Path,help='Explicit compact FIT source; otherwise packaged source is preferred.')
    ap.add_argument('--output-dir',type=Path,help='Optional isolated replay destination; source paths stay package-relative.')
    args=ap.parse_args()
    COMPACT_SOURCE_OVERRIDE=args.compact_source
    if args.output_dir is not None:R=args.output_dir.resolve()
    if args.action=='prepare':prepare()
    elif args.action=='capture-audit':print(json.dumps(capture_audit(),ensure_ascii=False))
    elif args.action=='seen-all-eval':print(json.dumps(add_seen_all_evaluation(),ensure_ascii=False))
    else:print(json.dumps(export_compact_source(),ensure_ascii=False))
