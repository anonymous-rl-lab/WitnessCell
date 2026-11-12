"""Per-edge reproducibility masks for effective intervention-response relations.

This score is NOT a probability that an edge is causal, a p value, an FDR
certificate, or a guarantee of predictive improvement.  It combines agreement
between existing complementary cell halves with conditional sensitivity to
resampling the permitted training double conditions.  The original cells and
the available intervention identities are not expanded by this resampling.
"""
from __future__ import annotations
from typing import Any
import numpy as np
from relation_model import fit_relation

DEFAULT_BOOTSTRAPS = 32
THRESHOLDS = (0.0, 0.25, 0.5, 0.75)


def _matrix(value: Any, name: str) -> np.ndarray:
    a=np.asarray(value,dtype=np.float64)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or not np.isfinite(a).all():
        raise ValueError(name+' must be a finite square matrix')
    return a


def stability_scores(W: Any, W_first: Any, W_second: Any,
                     bootstrap_W: Any) -> dict[str,np.ndarray]:
    """Calculate sign-aware, scale-invariant reproducibility for every edge.

    q_half = max(0, 2*w1*w2/(w1**2+w2**2)), additionally requiring both half
    signs to equal the full-fit sign.  s_boot=max(0,2*p_same_sign-1).
    Final r=q_half*s_boot.  Negative but stable coefficients are retained.
    Diagonal and exact-zero full-fit coefficients always receive score zero.
    """
    w=_matrix(W,'W');a=_matrix(W_first,'W_first');b=_matrix(W_second,'W_second')
    boot=np.asarray(bootstrap_W,dtype=np.float64)
    if a.shape!=w.shape or b.shape!=w.shape or boot.ndim!=3 or boot.shape[1:]!=w.shape or not len(boot) or not np.isfinite(boot).all():
        raise ValueError('All coefficient matrices must share the same shape; bootstrap_W must be nonempty and finite')
    if np.any(np.diag(w)!=0) or np.any(np.diag(a)!=0) or np.any(np.diag(b)!=0) or np.any(np.diagonal(boot,axis1=1,axis2=2)!=0):
        raise ValueError('Relation coefficients must have exactly zero diagonals')
    scale=np.maximum(np.abs(a),np.abs(b))
    an=np.divide(a,scale,out=np.zeros_like(a),where=scale>0)
    bn=np.divide(b,scale,out=np.zeros_like(b),where=scale>0)
    denominator=an*an+bn*bn
    agreement=np.divide(2*an*bn,denominator,out=np.zeros_like(an),where=denominator>0).clip(0,1)
    full_sign=np.sign(w)
    consistent=(full_sign!=0)&(np.sign(a)==full_sign)&(np.sign(b)==full_sign)
    agreement*=consistent
    fraction=np.mean(np.sign(boot)==full_sign[None,:,:],axis=0)
    boot_support=np.maximum(0.,2*fraction-1.)
    boot_support[full_sign==0]=0.
    score=(agreement*boot_support).clip(0,1)
    np.fill_diagonal(score,0.);np.fill_diagonal(fraction,0.);np.fill_diagonal(agreement,0.);np.fill_diagonal(boot_support,0.)
    return {'edge_stability':score,'half_agreement':agreement,
            'bootstrap_sign_fraction':fraction,'bootstrap_sign_support':boot_support}


def fit_edge_stability(U: Any, Z: Any, Z_first: Any, Z_second: Any,
                       scales: Any, ridge: float, observation_weights: Any=None,
                       n_single: int=14, seed: int=20261101,
                       n_bootstrap: int=DEFAULT_BOOTSTRAPS,
                       mode: str='masked') -> dict[str,Any]:
    """Fit full/half/condition-bootstrap graphs using only supplied fit rows.

    Rows [0:n_single] are the complete allowed single-intervention reference,
    kept once in every bootstrap.  Remaining rows must be allowed doubles;
    each bootstrap samples exactly that many double rows with replacement.
    The same supplied full-fit scales and weights are used in both halves.
    This is conditional condition-resampling, not a new biological replicate.
    """
    u=np.asarray(U);z=np.asarray(Z,dtype=np.float64)
    za=np.asarray(Z_first,dtype=np.float64);zb=np.asarray(Z_second,dtype=np.float64)
    if z.ndim!=2 or u.shape!=z.shape or za.shape!=z.shape or zb.shape!=z.shape:
        raise ValueError('U, Z, Z_first, Z_second must have the same two-dimensional shape')
    if not np.isfinite(z).all() or not np.isfinite(za).all() or not np.isfinite(zb).all():
        raise ValueError('All fitting responses must be finite')
    for value,name,minimum in ((n_single,'n_single',1),(n_bootstrap,'n_bootstrap',1),(seed,'seed',0)):
        if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or value<minimum:
            raise ValueError(name+' must be a valid integer')
    n_single=int(n_single);n_bootstrap=int(n_bootstrap);seed=int(seed)
    if n_single>len(z) or not np.all(u[:n_single].sum(axis=1)==1) or not np.all(u[n_single:].sum(axis=1)==2):
        raise ValueError('First rows must be singles; remaining fitting rows must be distinct-target doubles')
    weights=None if observation_weights is None else np.asarray(observation_weights,dtype=np.float64)
    if weights is not None and weights.shape not in ((len(z),),z.shape):
        raise ValueError('observation_weights has the wrong shape')
    full=fit_relation(u,z,scales,ridge,mode=mode,observation_weights=weights)
    first=fit_relation(u,za,scales,ridge,mode=mode,observation_weights=weights)
    second=fit_relation(u,zb,scales,ridge,mode=mode,observation_weights=weights)
    rng=np.random.default_rng(seed);single_ids=np.arange(n_single);double_ids=np.arange(n_single,len(z))
    draws=[];graphs=[]
    for _ in range(n_bootstrap):
        sampled=rng.choice(double_ids,size=len(double_ids),replace=True) if len(double_ids) else np.empty(0,dtype=int)
        ids=np.concatenate((single_ids,sampled));draws.append(sampled)
        graph=fit_relation(u[ids],z[ids],scales,ridge,mode=mode,
                           observation_weights=None if weights is None else weights[ids])
        graphs.append(graph['W'])
    bootstrap=np.stack(graphs)
    score=stability_scores(full['W'],first['W'],second['W'],bootstrap)
    r=score['edge_stability'];off=~np.eye(r.shape[0],dtype=bool)
    audit={'definition':'r = half coefficient agreement * max(0,2*bootstrap sign agreement-1)',
        'orientation':'W[response,source]; preserves signed coefficients when masking',
        'interpretation':'reproducibility score; not causal probability, FDR or independent replication',
        'bootstrap_unit':'allowed training double-condition row; singles retained once',
        'bootstrap_count':n_bootstrap,'seed':seed,'single_rows':n_single,'double_rows':len(double_ids),
        'same_full_fit_scales_and_weights_in_halves':True,'mode':mode,
        'positive_score_edge_count':int(np.sum(r[off]>0)),
        'score_quantiles':np.quantile(r[off],[0,.25,.5,.75,1]).tolist(),
        'full_relation':full['audit'],'first_relation':first['audit'],'second_relation':second['audit']}
    return {'W':full['W'],'W_first':first['W'],'W_second':second['W'],
            'scales':full['scales'],'bootstrap_W':bootstrap,
            'bootstrap_double_indices':np.stack(draws),**score,'audit':audit}


def threshold_relation(W: Any, edge_stability: Any, threshold: float) -> dict[str,Any]:
    """Hard edge mask; threshold zero is the exact unfiltered comparator.

    The retained coefficient is W itself, not W*r.  This separates threshold
    selection from an additional continuous coefficient-shrinkage hypothesis.
    """
    w=_matrix(W,'W');r=np.asarray(edge_stability,dtype=np.float64)
    if r.shape!=w.shape or not np.isfinite(r).all() or np.any(r<0) or np.any(r>1):
        raise ValueError('edge_stability must lie in [0,1] and match W')
    if np.any(np.diag(w)!=0):raise ValueError('W must have an exactly zero diagonal')
    if isinstance(threshold,(bool,np.bool_)) or not np.isscalar(threshold) or not np.isfinite(threshold) or not 0<=float(threshold)<=1:
        raise ValueError('threshold must be a finite scalar in [0,1]')
    tau=float(threshold)
    mask=(w!=0) if tau==0 else ((r>=tau)&(w!=0))
    np.fill_diagonal(mask,False)
    masked=np.where(mask,w,0.)
    return {'W':masked,'mask':mask,'audit':{'threshold':tau,
        'threshold_zero_means_exact_unfiltered':tau==0,'retained_edges':int(mask.sum()),
        'retained_positive_edges':int(np.sum(masked>0)),'retained_negative_edges':int(np.sum(masked<0)),
        'soft_coefficient_weighting_used':False}}
