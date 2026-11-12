"""Portable intervention identities and read-only, explicitly masked fit views."""
from pathlib import Path
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import hashlib,importlib.util,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
ASSETS=ROOT.parent/'nextpro2_assets_20261102/extracted'
UP=ASSETS/'WitnessCell_R_Upstream_Selection_20261101'
PK=ASSETS/'WitnessCell_PairPick_External_PK_20261101/external_pairpick_pk_20261101'
PRE=ASSETS/'WitnessCell_TargetRank_External_Preflight_20261101'
sys.path.insert(0,str(UP/'engine'))
sys.path.insert(0,str(UP/'engine/source/adapters/sf12_original/original_code'))
import global_pair_model as gm
import baseline_core as bc
import relation_model as rm
import edge_stability as es
import bridge_features as br
def module(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
rank=module('archived_ranker',UP/'study/train_ranker.py');rank.configure(UP/'engine')
portable=module('portable_targetrank',PRE/'code/portable_targetrank.py')
linear=module('archived_linear',PK/'train_linear.py')
def load(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
def writej(p,o):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(o,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def family(p):return '+'.join(sorted(x.rsplit('-',1)[0] for x in p.split('+')))
def folds(allowed,pairs,blocks=None,k=3,tag='aux'):
 buckets={}
 for i in allowed:buckets.setdefault(family(pairs[i]),[]).append(int(i))
 keys=sorted(buckets,key=lambda s:hashlib.sha256(('GSE220974_TargetRank_20261102|'+tag+'|'+s.replace('+','|')).encode()).hexdigest())
 fs=[[] for _ in range(k)]
 for key in keys:fs[min(range(k),key=lambda j:(len(fs[j]),j))]+=buckets[key]
 return [np.array(sorted(x),int) for x in fs if x]
gm.preserved_folds=folds
def rectangles(pairs,nodes):
 lookup={frozenset(p.split('+')):p for p in pairs};rows=set()
 from itertools import combinations
 for a,b in combinations(nodes,2):
  shared=[n for n in nodes if n not in (a,b) and frozenset((a,n)) in lookup and frozenset((b,n)) in lookup]
  for c,d in combinations(shared,2):rows.add(tuple(sorted(lookup[frozenset(e)] for e in [(a,c),(a,d),(b,c),(b,d)])))
 return [list(x) for x in sorted(rows)]
def setup(design='detected49'):
 m=load(ROOT/'data/moments.npz');m['conditions']=np.char.replace(m['conditions'],'|','+')
 singles=m['target_names'].astype(str).tolist();pairs=[p for p in m['conditions'].astype(str) if '+' in p]
 if design!='detected49':
  whitelist=json.loads((ROOT/'data/DESIGN_WHITELIST.json').read_text())['regular_pairs'];pairs=[p for p in pairs if p in whitelist]
 nodes=singles;lookup={p:i for i,p in enumerate(m['conditions'])};ids=np.array([lookup[p] for p in pairs]);mask=m['score_gene_mask']
 common={k:m[k] for k in ['genes','target_names','target_indices','control','score_gene_mask']};common['all_pair_names']=np.array(pairs)
 blocks=rectangles(pairs,nodes);full=gm.make_dataset(pairs,dict(zip(nodes,m['single_effects'])),np.full((len(pairs),len(m['genes'])),np.nan),blocks)
 score=gm.make_dataset(pairs,dict(zip(nodes,m['single_effects'][:,mask])),np.full((len(pairs),int(mask.sum())),np.nan),blocks)
 return dict(m=m,pair_moment_ids=ids,common=common,full=full,score=score,design=design)
def masked_view(base,train):
 # Only allowed double rows enter fit datasets. Evaluation arrays stay in moments.
 full=base['full'];score=base['score'];m=base['m'];train=np.asarray(train,int)
 y=np.full_like(full.y,np.nan);y[train]=m['means'][base['pair_moment_ids'][train]]
 f=gm.Dataset(full.pairs,full.nodes,y,full.additive,full.x,full.pair_kernel,full.blocks)
 s=gm.Dataset(score.pairs,score.nodes,y[:,m['score_gene_mask']],score.additive,score.x,score.pair_kernel,score.blocks)
 moments={k:m[k].copy() for k in ['genes','conditions','first_means','second_means']}
 allowed={full.pairs[i] for i in train};forbidden=np.array(['+' in n and n not in allowed for n in moments['conditions']])
 for k in ('first_means','second_means'):moments[k][forbidden]=np.nan
 return f,s,moments
def graph_features(base,train):
 """19 physical RNA ports, 35 distinct interventions; no duplicate gene equations."""
 m=base['m'];full=base['full'];train=np.sort(train);names=m['target_names'].astype(str)
 genes=m['genes'].astype(str);physical=m['base_genes'].astype(str);ports=np.array([np.flatnonzero(genes==g)[0] for g in physical]);pmap={g:i for i,g in enumerate(physical)}
 ip=np.array([pmap[n.rsplit('-',1)[0]] for n in names]);U0=np.eye(len(physical))[ip]
 edges=np.array([np.flatnonzero(r) for r in full.x]);ep=ip[edges]
 UD=np.zeros((len(full.pairs),len(physical)))
 for j,(a,b) in enumerate(ep):UD[j,[a,b]]=1
 si=np.arange(1,1+len(names));di=base['pair_moment_ids'][train];ids=np.r_[si,di]
 E=m['single_effects'][:,ports]
 sv=m['mean_variances'][si][:,ports]+m['mean_variances'][0,ports]
 scales=np.maximum(np.sqrt(np.mean(E*E+sv,0)),np.finfo(float).eps)
 variance=m['mean_variances'][ids][:,ports]+m['mean_variances'][0,ports]
 stable=es.fit_edge_stability(np.r_[U0,UD[train]],m['means'][ids][:,ports],m['first_means'][ids][:,ports],m['second_means'][ids][:,ports],
  scales,ridge=1.,observation_weights=br._weights(variance),n_single=len(si),n_bootstrap=32,seed=20261101)
 W=es.threshold_relation(stable['W'],stable['edge_stability'],.5)['W'];En=E/scales
 def solve(W):
  pred=full.x@En
  for j,((a,b),fixed) in enumerate(zip(edges,ep)):
   free=np.setdiff1d(np.arange(len(physical)),fixed);clamp=En[[a,b],fixed];pred[j,fixed]=clamp
   eq=np.eye(len(free))-W[np.ix_(free,free)];rhs=W[np.ix_(free,fixed)]@clamp
   design=np.vstack([eq,np.sqrt(1e-3)*np.eye(len(free))]);response=np.r_[rhs-eq@pred[j,free],np.zeros(len(free))]
   pred[j,free]+=np.linalg.lstsq(design,response,rcond=None)[0]
  return pred
 raw_gene=solve(W)-solve(np.zeros_like(W))
 # Ranking head coordinates remain intervention-indexed, as in the portable API.
 # Expansion happens AFTER the physical graph solve; no a/i alias is treated as a distinct RNA gene.
 raw=raw_gene[:,ip];coef=np.linalg.pinv(full.x[train])@raw[train];F=raw-full.x@coef
 norm=float(np.sqrt(np.mean(np.sum(F[train]**2,axis=1))));rnorm=float(np.sqrt(np.mean(np.sum(raw[train]**2,axis=1))))
 if norm*norm<=1e-12*max(rnorm*rnorm,1e-300):F=np.zeros_like(F);norm=0.
 else:F/=norm
 assert np.isfinite(F).all() and np.max(abs(full.x[train].T@F[train]))<1e-7
 return F,raw,dict(W=W,unmasked_W=stable['W'],stability=stable['edge_stability'],fit_ids=train,physical_ports=physical,intervention_to_port=ip,normalization=np.array(norm))
def fit_upstream(base,train,query,out):
 out=Path(out);train=np.sort(train);query=np.asarray(query,int)
 key=hashlib.sha256((sha(__file__)+'|'+sha(ROOT/'data/moments.npz')+'|'+base['design']+'|'+','.join(map(str,train))).encode()).hexdigest()[:16]
 folder=out/'cache'/key;folder.mkdir(parents=True,exist_ok=True)
 required=np.setdiff1d(np.arange(len(base['full'].pairs)),train)
 if (folder/'predictions.npz').exists():
  z=load(folder/'predictions.npz');assert np.array_equal(z['train'],train) and np.array_equal(z['query'],required)
 else:
  full,score,mom=masked_view(base,train);inner=folds(train,full.pairs)
  pred,meta=rank.fit_native_rg(base['common'],dict(single_effects=base['m']['single_effects']),full,score,train,required,inner,mom)
  F,raw,state=graph_features(base,train)
  z=dict(train=train,query=required,prediction=pred,graph=F,raw=raw)
  np.savez_compressed(folder/'predictions.npz',**z);np.savez_compressed(folder/'graph_state.npz',**state);writej(folder/'selection.json',meta)
 loc={int(i):j for j,i in enumerate(required)};take=[loc[int(i)] for i in query]
 return z['prediction'][take],z['graph'][query],z['raw'][query]
