"""Input-only SF12 and neural interaction kernels for the unchanged R backbone.

No function here reads a combination response. Each feature scale that uses
combination identities is recomputed from the current fitting edges. Every
kernel uses the original R endpoint projection and train-trace normalization.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
import hashlib,json,sys
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source/adapters/sf12_original/original_code'))
import legacy_shared_core as legacy
import global_pair_model as gm

def silu(x):
    """Stable float64 SiLU, without large positive exponential arguments."""
    x=np.asarray(x,dtype=np.float64)
    e=np.exp(-np.abs(x))
    sigmoid=np.where(x>=0,1./(1.+e),e/(1.+e))
    return x*sigmoid

def neural_pair_features(a,b):
    """Symmetric joint-minus-single tanh/SiLU features; zero if either is zero."""
    a,b=np.asarray(a,dtype=np.float64),np.asarray(b,dtype=np.float64)
    if a.shape!=b.shape:raise ValueError('Neural endpoint feature shapes differ')
    return np.concatenate([np.tanh(a+b)-(np.tanh(a)+np.tanh(b)),
                           silu(a+b)-(silu(a)+silu(b))],axis=-1)

def _ids(ids,total,name):
    ids=np.asarray(ids)
    if ids.ndim!=1 or ids.dtype.kind not in 'iu' or len(set(ids.tolist()))!=len(ids):
        raise ValueError(name+' must contain unique integer indices')
    if np.any(ids<0) or np.any(ids>=total):raise ValueError(name+' outside pair axis')
    return ids.astype(np.int64,copy=False)

def _inputs(common,single_effects,score_dataset,train_ids,query_ids):
    e=np.asarray(single_effects,dtype=np.float64)
    score=np.asarray(common['score_gene_mask'],dtype=bool)
    nt=np.asarray(common['non_target'],dtype=bool)
    names=np.asarray(common['target_names']).astype(str)
    knames=np.asarray(common['knowledge_names']).astype(str)
    if not np.array_equal(names,knames):raise ValueError('Knowledge target identity mismatch')
    if not np.array_equal(names,np.asarray(score_dataset.nodes).astype(str)):
        raise ValueError('Dataset target identity mismatch')
    knowledge=np.asarray(common['knowledge'],dtype=np.float64)
    if e.ndim!=2 or e.shape[0]!=len(names) or e.shape[1]!=len(score) or nt.shape!=score.shape:
        raise ValueError('Single-response or RNA mask shape mismatch')
    if knowledge.ndim!=2 or knowledge.shape[0]!=len(names):raise ValueError('Knowledge shape mismatch')
    if not np.isfinite(e).all() or not np.isfinite(knowledge).all():raise ValueError('Nonfinite input feature')
    x=np.asarray(score_dataset.x,dtype=np.float64)
    if x.shape!=(len(score_dataset.pairs),len(names)) or not np.all((x==0)|(x==1)) or not np.all(x.sum(1)==2):
        raise ValueError('Expected distinct-double pair incidence')
    edges=np.array([np.flatnonzero(row) for row in x],dtype=np.int64)
    tr=_ids(train_ids,len(x),'train_ids');q=_ids(query_ids,len(x),'query_ids')
    if len(tr)<2:raise ValueError('Need at least two fitting pairs')
    kr=np.asarray(score_dataset.pair_kernel,dtype=np.float64)
    if kr.shape!=(len(x),len(x)) or not np.isfinite(kr).all():raise ValueError('Invalid original R kernel')
    return e,score,nt,knowledge,x,edges,tr,q,kr

def build_kernels(common,single_effects,score_dataset,train_ids,query_ids):
    """Return R/SF12/NEURAL projected float64 ktt/kqt and input-only audit.

    ``score_dataset`` needs only pairs, nodes, x and pair_kernel attributes.
    It deliberately need not expose a y attribute. Query IDs can also contain
    fitting IDs when a caller requests training predictions.
    """
    e,mask,nt,knowledge,x,edges,tr,q,kr=_inputs(common,single_effects,score_dataset,train_ids,query_ids)
    # Preserve all archived SF12 per-RNA preprocessing, then select only the
    # intersection with the official scoring axis. One single-only intervention
    # port is retained by scoring but protected by SF12, yielding 4994 features.
    centered=e[:,nt]-e[:,nt].mean(0,keepdims=True)
    _,full_sf=legacy.features(centered,knowledge,edges[tr],edges)
    f12=full_sf[:,mask[nt],:].reshape(len(edges),-1)
    if f12.shape[1]==0:raise ValueError('No SF12 scoring features')
    k12=(f12@f12.T)/f12.shape[1]
    # Neural-style nonlinear interaction geometry. Its inputs are uncentered
    # measured single effects; no paired-cell or response labels are consumed.
    es=e[:,mask]
    global_rms=float(np.sqrt(np.mean(es*es)))
    scale=np.maximum(np.sqrt(np.mean(es*es,axis=0)),global_rms)
    z=np.divide(es,scale,out=np.zeros_like(es),where=scale>0)
    fn=neural_pair_features(z[edges[:,0]],z[edges[:,1]])
    kn=(fn@fn.T)/fn.shape[1]
    out={};audits={}
    for name,k in [('R',kr),('SF12',k12),('NEURAL',kn)]:
        tt,qt,qq,audit=gm.projected_blocks(x[tr],x[q],k[np.ix_(tr,tr)],k[np.ix_(q,tr)],k[np.ix_(q,q)])
        out[name]={'ktt':np.asarray(tt,dtype=np.float64),'kqt':np.asarray(qt,dtype=np.float64),
                   'kqq':np.asarray(qq,dtype=np.float64),'raw_all_pairs':k}
        audits[name]={**audit,'minimum_train_eigenvalue':float(np.linalg.eigvalsh(tt).min()),
                      'normalized_train_trace':float(np.trace(tt)),
                      'query_training_endpoint_annihilation':float(np.max(abs(qt@x[tr]),initial=0.))}
    out['audit']={'no_response_attribute_access':True,'fit_pair_count':len(tr),'query_pair_count':len(q),
                  'train_ids':tr.tolist(),'query_ids':q.tolist(),'R_single_columns':int(mask.sum()),
                  'SF12_single_columns':int(np.sum(mask&nt)),'SF12_flat_dimension':f12.shape[1],
                  'neural_flat_dimension':fn.shape[1],'neural_single_global_rms':global_rms,
                  'knowledge_shape':list(knowledge.shape),'knowledge_identity_checked':True,
                  'neural_features':'concat(tanh(a+b)-tanh(a)-tanh(b),silu(a+b)-silu(a)-silu(b)); no multiplication by RNA scale',
                  'kernels':audits}
    return out

def mix_kernels(bundle,family,eta):
    """Return (ktt,kqt), using native R exactly when eta is zero."""
    if family not in ('native','sf12','neural','joint'):raise ValueError('Unknown kernel family')
    if isinstance(eta,(bool,np.bool_)) or not np.isfinite(eta) or eta<0:raise ValueError('eta must be nonnegative')
    if eta==0:return bundle['R']['ktt'].copy(),bundle['R']['kqt'].copy()
    if family=='native':raise ValueError('native kernel requires eta=0')
    def added(key):
        if family=='sf12':return bundle['SF12'][key]
        if family=='neural':return bundle['NEURAL'][key]
        return .5*(bundle['SF12'][key]+bundle['NEURAL'][key])
    return tuple((bundle['R'][key]+float(eta)*added(key))/(1.+float(eta)) for key in ('ktt','kqt'))

def smoke():
    """Identity/feature-only smoke: do not load any combination response array."""
    from types import SimpleNamespace
    inp=ROOT/'source/adapters/sf12_original/inputs'
    with np.load(inp/'common.npz',allow_pickle=False) as f:
        common={k:f[k].copy() for k in ('genes','target_names','knowledge_names','knowledge','non_target','score_gene_mask','all_pair_names')}
    with np.load(inp/'fold_0.npz',allow_pickle=False) as f:
        singles=f['single_effects'].copy();train_names=f['train_names'].astype(str);query_names=f['query_names'].astype(str)
    pairs=common['all_pair_names'].astype(str).tolist();names=common['target_names'].astype(str).tolist()
    edges=np.array([[names.index(a) for a in pair.split('+')] for pair in pairs])
    es=singles[:,common['score_gene_mask']]
    # Namespace has no response attribute, so any accidental y access fails.
    d=SimpleNamespace(pairs=pairs,nodes=names,x=gm.incidence(pairs,names),
                      pair_kernel=gm.symmetric_kernel(edges,es@es.T/es.shape[1]))
    tr=np.array([pairs.index(n) for n in train_names]);q=np.array([pairs.index(n) for n in query_names])
    b=build_kernels(common,singles,d,tr,q)
    checks=[]
    def check(name,value):
        checks.append({'name':name,'passed':bool(value)})
        if not value:raise AssertionError(name)
    for name in ('R','SF12','NEURAL'):
        k=b[name];ev=np.linalg.eigvalsh(k['ktt'])
        check(name+'_PSD',float(ev.min())>=-1e-10*max(1.,float(ev.max())))
        check(name+'_training_endpoint_annihilation',np.max(abs(d.x[tr].T@k['ktt']))<1e-9)
        check(name+'_query_endpoint_annihilation',np.max(abs(k['kqt']@d.x[tr]))<1e-9)
        check(name+'_trace',abs(np.trace(k['ktt'])-len(tr))<1e-10)
    tt,qt,_,_=gm.projected_blocks(d.x[tr],d.x[q],d.pair_kernel[np.ix_(tr,tr)],d.pair_kernel[np.ix_(q,tr)],d.pair_kernel[np.ix_(q,q)])
    rng=np.random.default_rng(20261101);resid=rng.normal(size=(len(tr),5))
    for family in ('native','sf12','neural','joint'):
        a,c=mix_kernels(b,family,0.)
        check(family+'_eta0_exact',np.array_equal(a,tt) and np.array_equal(c,qt))
        for lam in (.1,1.,10.):
            v,u=np.linalg.eigh(a);op=(c@u/(v+lam))@u.T
            v0,u0=np.linalg.eigh(tt);op0=(qt@u0/(v0+lam))@u0.T
            check(family+'_eta0_operator_'+str(lam),np.array_equal(op@resid,op0@resid))
    aa=rng.normal(size=(9,es.shape[1]));bb=rng.normal(size=aa.shape)
    check('neural_AB_symmetry',np.array_equal(neural_pair_features(aa,bb),neural_pair_features(bb,aa)))
    check('neural_single_zero',np.array_equal(neural_pair_features(aa,np.zeros_like(aa)),np.zeros((len(aa),2*aa.shape[1]))))
    nt=common['non_target'];cent=singles[:,nt]-singles[:,nt].mean(0,keepdims=True)
    ft,fq=legacy.features(cent,common['knowledge'],edges[tr],edges[q])
    fs,qs=legacy.features(cent,common['knowledge'],edges[tr][:,::-1],edges[q][:,::-1])
    check('SF12_AB_symmetry',np.array_equal(ft,fs) and np.array_equal(fq,qs))
    wrong={**common,'knowledge_names':common['knowledge_names'][::-1]}
    try:build_kernels(wrong,singles,d,tr,q)
    except ValueError:rejected=True
    else:rejected=False
    check('wrong_knowledge_identity_rejected',rejected)
    for family in ('sf12','neural','joint'):
        for eta in (.25,1.):
            a,c=mix_kernels(b,family,eta)
            check(family+'_'+str(eta)+'_PSD',np.linalg.eigvalsh(a).min()>-1e-9)
            check(family+'_'+str(eta)+'_endpoint',np.max(abs(c@d.x[tr]))<1e-9)
    result={'status':'PASS','checks':checks,'check_count':len(checks),'input_only_audit':b['audit'],
            'combination_response_arrays_loaded':False,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'KERNEL_FEATURE_SMOKE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':len(checks),'dimensions':{'SF12':b['audit']['SF12_flat_dimension'],'NEURAL':b['audit']['neural_flat_dimension']}}))

if __name__=='__main__':smoke()
