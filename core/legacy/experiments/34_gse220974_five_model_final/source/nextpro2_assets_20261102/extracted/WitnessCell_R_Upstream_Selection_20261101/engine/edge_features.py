"""Per-edge reproducibility masks for the unchanged R-G response predictor.

All directed relations, scores and masks are rebuilt on the current training
subset. Reproducibility scores are not probabilities of biological causation.
"""
from pathlib import Path
import numpy as np
from bridge_features import relation_inputs,load_moments,_weights,_ids,gm,SOLVE_RIDGE
from relation_model import predict_pair_features
from edge_stability import fit_edge_stability,threshold_relation

ROOT=Path(__file__).resolve().parent
THRESHOLDS=(0.,.25,.5,.75)
BOOTSTRAP_SEED=20261101
BOOTSTRAP_COUNT=32
SHUFFLE_SEEDS=(20261101,20261102,20261103)
FAMILIES=('STABLE','SHUFFLED0','SHUFFLED1','SHUFFLED2')

def fitting_evidence(common,data,full,train,moments):
    inp=relation_inputs(common,data,full,train,moments)
    ports=inp['ports'];names=inp['names']
    lookup={str(n):i for i,n in enumerate(moments['conditions'])};ci=lookup['control']
    ids=[lookup[str(n)] for n in names]+[lookup[full.pairs[i]] for i in train]
    evidence=dict(U=np.concatenate((inp['U_single'],inp['U_double'])),
                  Z=np.concatenate((inp['Z_single'],inp['Z_double'])),
                  scales=inp['scales'],
                  observation_weights=_weights(np.concatenate((inp['single_var'],inp['double_var']))),
                  n_single=len(names))
    for arg,key in [('Z_first','first_means'),('Z_second','second_means')]:
        evidence[arg]=moments[key][ids][:,ports]-moments[key][ci,ports]
    return inp,evidence

def build_edge_kernels(common,data,full,score,train,query,moments,ridge):
    train=_ids(train,len(full.pairs));query=_ids(query,len(full.pairs))
    # Canonical row order makes bootstrap identities invariant to the ordering
    # of a caller's training-index list. Kernel rows still follow that caller.
    canonical_train=np.sort(train)
    inp,evidence=fitting_evidence(common,data,full,canonical_train,moments)
    stable=fit_edge_stability(**evidence,ridge=ridge,seed=BOOTSTRAP_SEED,
                             n_bootstrap=BOOTSTRAP_COUNT,mode='masked')
    E=inp['single_ports']/inp['scales']
    pairs=np.array([np.flatnonzero(u) for u in full.x],int)
    zero=predict_pair_features(np.zeros_like(stable['W']),E,pairs,solve_ridge=SOLVE_RIDGE)
    kr=score.pair_kernel;x=score.x
    tt,qt,qq,a=gm.projected_blocks(x[train],x[query],kr[np.ix_(train,train)],
                                  kr[np.ix_(query,train)],kr[np.ix_(query,query)])
    out={'R':dict(ktt=tt,kqt=qt,kqq=qq)}
    state={k:v for k,v in stable.items() if isinstance(v,np.ndarray)}
    state.update(scales=inp['scales'],train_ids=train,canonical_train_ids=canonical_train,
                 query_ids=query,single_ports=inp['single_ports'])
    graphs={family:{} for family in FAMILIES}
    permutations={f'SHUFFLED{i}':np.random.default_rng(seed).permutation(len(E))
                  for i,seed in enumerate(SHUFFLE_SEEDS)}
    for family in FAMILIES:out[family]={}
    for tau in THRESHOLDS:
        masked=threshold_relation(stable['W'],stable['edge_stability'],tau)
        W0=masked['W']
        for family in FAMILIES:
            W=W0 if family=='STABLE' else W0[np.ix_(permutations[family],permutations[family])]
            pred=predict_pair_features(W,E,pairs,solve_ridge=SOLVE_RIDGE)
            features=pred['prediction']-zero['prediction']
            k=features@features.T/features.shape[1]
            tt,qt,qq,ka=gm.projected_blocks(x[train],x[query],k[np.ix_(train,train)],
                                          k[np.ix_(query,train)],k[np.ix_(query,query)])
            out[family][float(tau)]=dict(ktt=tt,kqt=qt,kqq=qq)
            key=f'{family}__t{tau:g}'
            state['masked_W__'+key]=W
            state['features__'+key]=features
            state['port_prediction__'+key]=pred['prediction']*inp['scales']
            count=int(np.count_nonzero(W));pos=int(np.count_nonzero(W>0));neg=int(np.count_nonzero(W<0))
            assert count==np.count_nonzero(W0)
            graphs[family][str(float(tau))]={
                'threshold':float(tau),'retained_edge_count':count,'positive_edges':pos,'negative_edges':neg,
                'threshold_meaning':'exact unfiltered reference' if tau==0 else 'edge_stability >= threshold',
                'identity_permutation':None if family=='STABLE' else permutations[family].tolist(),
                'kernel':ka,'solver':pred['audit'],
                'feature_rms':float(np.sqrt(np.mean(features*features))),
                'min_train_eigenvalue':float(np.linalg.eigvalsh(tt).min())}
    out['UNFILTERED']=out['STABLE'][0.]
    out['audit']={'fit_pair_ids':canonical_train.tolist(),'query_pair_ids':query.tolist(),
        'graph_ridge':float(ridge),'target_names':inp['names'].tolist(),'thresholds':list(THRESHOLDS),
        'single_rows':len(E),'double_rows':len(train),'bootstrap_seed':BOOTSTRAP_SEED,
        'bootstrap_count':BOOTSTRAP_COUNT,'edge_score':stable['audit'],
        'graphs':graphs,'native_kernel':a,
        'permutation_scope':'after thresholding; exact retained-edge and signed-weight multiset match per candidate; not node-wise degree preserving',
        'selection_scope':'inner-training conditions only; neither validation nor query RNA informs relation masks',
        'interpretation':'empirical reproducibility of effective RNA relations; not causal probability or FDR guarantee'}
    out['state']=state
    return out
