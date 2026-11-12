"""Training-only intervention relation features for an unchanged R-G predictor.

The small structural model is a hypothesis about effective RNA responses, not
an identified direct gene regulatory network. Targeted equations are excluded
from relation fitting. All other RNA output labels remain in the original R
regression. Outer-query target RNA is never an input to this module.
"""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'source/adapters/sf12_original/original_code'))
import global_pair_model as gm
from relation_model import fit_relation, predict_pair_features

MODALITIES = ('SINGLE','DUAL','UNMASKED','SHUFFLED0','SHUFFLED1','SHUFFLED2')
SOLVE_RIDGE = 1e-3
SHUFFLE_SEEDS = (20261101,20261102,20261103)

def load_moments():
    with np.load(ROOT/'source/adapters/original_moments.npz',allow_pickle=False) as z:
        return {k:z[k].copy() for k in z.files}

def _ids(ids,n):
    ids=np.asarray(ids)
    if ids.ndim!=1 or ids.dtype.kind not in 'iu' or len(np.unique(ids))!=len(ids):
        raise ValueError('unique integer condition identities required')
    if np.any(ids<0) or np.any(ids>=n):raise ValueError('condition identity out of bounds')
    return ids.astype(int)

def _weights(var):
    # Stabilized inverse mean-estimation variance. These are optimization
    # weights, not independent-observation likelihoods: controls are shared.
    med=np.median(var,axis=0)
    floor=np.maximum(med,np.finfo(float).eps)
    w=1./(var+floor)
    return w/w.mean(axis=0,keepdims=True)

def relation_inputs(common,data,full,train,moments):
    names=common['target_names'].astype(str)
    ports=common['target_indices'].astype(int)
    if not np.array_equal(common['genes'][ports],names):raise ValueError('target mapping mismatch')
    if not np.array_equal(common['genes'],moments['genes']):raise ValueError('moment gene mapping mismatch')
    lookup={str(x):i for i,x in enumerate(moments['conditions'])}
    ci=lookup['control']; si=np.array([lookup[n] for n in names])
    ti=np.array([lookup[full.pairs[i]] for i in train])
    E=np.asarray(data['single_effects'][:,ports],float)
    control=common['control'][ports]
    if not np.allclose(E,moments['means'][si][:,ports]-control,rtol=0,atol=1e-12):
        raise ValueError('single response mismatch')
    V=moments['variances'];N=moments['counts']
    cv=V[ci,ports]/N[ci]
    sv=V[si][:,ports]/N[si,None]+cv
    tv=V[ti][:,ports]/N[ti,None]+cv
    # Only allowed fitting double responses are sliced. Validation and outer
    # query target readouts are unavailable to the graph learner.
    Zt=np.asarray(full.y[train][:,ports],float)
    if not np.isfinite(Zt).all():raise ValueError('missing allowed double response')
    if not np.allclose(Zt,moments['means'][ti][:,ports]-control,rtol=0,atol=1e-12):
        raise ValueError('double response mismatch')
    scales=np.sqrt(np.mean(E*E+sv,axis=0))
    scales=np.maximum(scales,np.finfo(float).eps)
    a=moments['first_means'][si][:,ports]-moments['first_means'][ci,ports]
    b=moments['second_means'][si][:,ports]-moments['second_means'][ci,ports]
    aa=np.diag(a);bb=np.diag(b);den=.5*(aa*aa+bb*bb)
    conf=np.divide(aa*bb,den,out=np.zeros_like(aa),where=den>0).clip(0,1)
    return dict(single_ports=E,scales=scales,single_var=sv,double_var=tv,
                U_single=np.eye(len(names)),Z_single=E,U_double=full.x[train],Z_double=Zt,
                target_confidence=conf,ports=ports,names=names)

def build_relation_kernels(common,data,full,score,train,query,moments,ridge):
    train=_ids(train,len(full.pairs));query=_ids(query,len(full.pairs))
    inp=relation_inputs(common,data,full,train,moments)
    Es=inp['single_ports']/inp['scales']
    all_edges=np.array([np.flatnonzero(u) for u in full.x],int)
    if all_edges.shape!=(len(full.pairs),2):raise ValueError('distinct double conditions required')
    confidence=np.ones(len(all_edges),dtype=np.float64)  # Own-target response repeatability is not an edge-confidence gate.
    us=inp['U_single'];zs=inp['Z_single'];scales=inp['scales']
    ud=np.concatenate((us,inp['U_double']));zd=np.concatenate((zs,inp['Z_double']))
    vs=inp['single_var'];vd=np.concatenate((vs,inp['double_var']))
    graphs={
        'SINGLE':fit_relation(us,zs,scales,ridge,mode='masked',observation_weights=_weights(vs)),
        'DUAL':fit_relation(ud,zd,scales,ridge,mode='masked',observation_weights=_weights(vd)),
        'UNMASKED':fit_relation(ud,zd,scales,ridge,mode='unmasked',observation_weights=_weights(vd)),
    }
    permutations={}
    for i,seed in enumerate(SHUFFLE_SEEDS):
        p=np.random.default_rng(seed).permutation(len(scales));permutations[f'SHUFFLED{i}']=p
        graphs[f'SHUFFLED{i}']={'W':graphs['DUAL']['W'][np.ix_(p,p)],'scales':scales,
            'audit':{'seed':seed,'permutation':p.tolist(),
                     'meaning':'graph relabelled while intervention and single-response identities stay fixed; not node-wise degree preserving'}}
    out={};audit={};state={'scales':scales,'target_confidence':inp['target_confidence'],
                          'pair_confidence':confidence,'train_ids':train,'query_ids':query}
    kr=score.pair_kernel
    tt,qt,qq,a=gm.projected_blocks(score.x[train],score.x[query],kr[np.ix_(train,train)],kr[np.ix_(query,train)],kr[np.ix_(query,query)])
    out['R']={'ktt':tt,'kqt':qt,'kqq':qq};audit['R']=a
    zero_pred=predict_pair_features(np.zeros((len(scales),len(scales))),Es,all_edges,solve_ridge=SOLVE_RIDGE)
    state['zero_graph_port_prediction']=zero_pred['prediction']*scales
    for name,graph in graphs.items():
        pred=predict_pair_features(graph['W'],Es,all_edges,solve_ridge=SOLVE_RIDGE)
        # Subtract the SAME clamp/solver with W=0 before applying the gate.
        # Otherwise a pair-specific gate can turn clamp bookkeeping into a
        # nonadditive kernel even in the complete absence of relation evidence.
        f=confidence[:,None]*(pred['prediction']-zero_pred['prediction'])
        k=f@f.T/f.shape[1]
        tt,qt,qq,ka=gm.projected_blocks(score.x[train],score.x[query],k[np.ix_(train,train)],k[np.ix_(query,train)],k[np.ix_(query,query)])
        out[name]={'ktt':tt,'kqt':qt,'kqq':qq}
        audit[name]={'relation':graph['audit'],'query_solver':pred['audit'],'kernel':ka,
                     'feature_rms':float(np.sqrt(np.mean(f*f))),
                     'minimum_kernel_eigenvalue':float(np.linalg.eigvalsh(tt).min()),
                     'effective_training_rank':int(np.linalg.matrix_rank(tt))}
        state['W__'+name]=graph['W'];state['features__'+name]=f
        state['port_prediction__'+name]=(pred['additive']+confidence[:,None]*pred['features'])*scales
        state['raw_port_prediction__'+name]=pred['prediction']*scales
    state['port_additive']=np.array([inp['single_ports'][a]+inp['single_ports'][b] for a,b in all_edges])
    out['audit']={'fit_pair_ids':train.tolist(),'query_pair_ids':query.tolist(),'ridge':float(ridge),
        'solve_ridge':SOLVE_RIDGE,'target_names':inp['names'].tolist(),
        'single_rows':len(us),'double_rows':len(train),'only_training_double_target_RNA_used':True,
        'target_self_half_agreement':inp['target_confidence'].tolist(),
        'target_self_half_agreement_used_for_features':False,
        'pair_feature_confidence':'one for every pair; own-target response agreement retained only as a diagnostic',
        'scale_source':'single RMS plus single/control mean-estimation variance only',
        'relation_feature':'port_prediction_W - port_prediction_W_zero; no own-target-confidence gate',
        'confidence_scope':'all pair weights are one; own-target half agreement is a diagnostic only',
        'interpretation':'effective intervention-masked RNA response relations; not certified direct gene regulation',
        'modalities':audit}
    out['state']=state
    return out
