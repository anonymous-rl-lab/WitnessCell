"""Fresh nested R-G, TargetRank, Stable and ridge fits. Never scores outer labels."""
import argparse,time
from common import *
def run(f,design,smoke=False):
 t=time.monotonic();base=setup(design);full=base['full'];m=base['m'];P=len(full.pairs);T=len(full.nodes);V=m['objectives'][:,m['score_gene_mask']]
 out=ROOT/('smoke_pipeline' if smoke else 'results')/design/'classical';folder=out/f'fold_{f}';folder.mkdir(parents=True,exist_ok=True)
 allids=np.arange(P);query=folds(allids,full.pairs,k=5,tag='outer')[f];train=np.setdiff1d(allids,query);inner=folds(train,full.pairs,tag='inner'+str(f))
 ff,ss,mom=masked_view(base,train);Y=ss.y@V.T;assert np.isnan(Y[query]).all()
 d=dict(X=full.x,single_effects=m['single_effects'][:,m['score_gene_mask']],objectives=V,single_utility=m['single_effects'][:,m['score_gene_mask']]@V.T,
  RG=np.full((P,T),np.nan),Y=Y,target_names=m['target_names'],pair_names=np.array(full.pairs),train=train,query=query)
 cases=[];oof_graph=np.zeros((P,T));oof_raw=np.zeros((P,T));groups=np.full(P,-1,int)
 for j,val in enumerate(inner):
  fit=np.setdiff1d(train,val);pr,gr,raw=fit_upstream(base,fit,val,out);d['RG'][val]=pr@V.T;oof_graph[val]=gr;oof_raw[val]=raw;groups[val]=j
  c=dict(train=fit,query=val,train_blocks=np.full(len(fit),-1,int),query_RG=pr@V.T,query_graph=gr,query_raw=raw,
   train_RG=np.zeros((len(fit),T)),train_graph=np.zeros((len(fit),T)),train_raw=np.zeros((len(fit),T)))
  sub=folds(fit,full.pairs,tag='aux'+str(f)+'_'+str(j));loc={int(i):k for k,i in enumerate(fit)}
  for k,hold in enumerate(sub):
   tt=np.setdiff1d(fit,hold);pred,graph,raw=fit_upstream(base,tt,hold,out);ix=np.array([loc[int(i)] for i in hold])
   c['train_RG'][ix]=pred@V.T;c['train_graph'][ix]=graph;c['train_raw'][ix]=raw;c['train_blocks'][ix]=k
  cases.append(c)
  print(json.dumps(dict(fold=f,stage='nested_features',inner=j,seconds=time.monotonic()-t)),flush=True)
 pr,gr,raw=fit_upstream(base,train,query,out);d['RG'][query]=pr@V.T
 outer=dict(train=train,query=query,train_blocks=groups[train],train_RG=d['RG'][train],query_RG=d['RG'][query],train_graph=oof_graph[train],query_graph=gr,train_raw=oof_raw[train],query_raw=raw)
 updated,meta=portable.train_ranker(d,cases,outer);stable,state=portable.old_branch(d,outer)
 predictions=dict(R_G=pr,SingleB_additive=full.additive[query][:,m['score_gene_mask']])
 utilities=dict(PairPick_TargetRank=updated,PairPick_Stable=d['RG'][query]+stable,R_G=d['RG'][query],SingleB_additive=predictions['SingleB_additive']@V.T)
 records=[]
 for lam in linear.RIDGES:
  sse=0;count=0;choices=[]
  for j,val in enumerate(inner):
   fit=np.setdiff1d(train,val);coef=linear.fit_subset(ss,m['single_effects'][:,m['score_gene_mask']],fit,lam);p=linear.predict(full.x[val],coef)
   err=p-ss.y[val];sse+=float(np.sum(err*err));count+=err.size
   choices+=linear.decision_rows(full.nodes,[full.pairs[i] for i in val],p,ss.y[val],V,j)
  assert choices
  records.append(dict(ridge=lam,rna_mse=sse/count,choice_utility=float(np.mean([r['actual_utility'] for r in choices]))))
 selected={'StandardRidge_RNA':linear.select_ridge(records,'rna_mse'),'StandardRidge_Choice':linear.select_ridge(records,'choice_utility',True)}
 for name,sel in selected.items():
  coef=linear.fit_subset(ss,m['single_effects'][:,m['score_gene_mask']],train,sel['ridge']);p=linear.predict(full.x[query],coef);predictions[name]=p;utilities[name]=p@V.T
  np.save(folder/(name+'_coefficients.npy'),coef)
 np.savez_compressed(folder/'predictions.npz',query=query,train=train,conditions=np.array(full.pairs)[query],genes=m['genes'][m['score_gene_mask']],target_names=m['target_names'],
  **{'pred_effect_'+k:v for k,v in predictions.items()},**{'utility_'+k:v for k,v in utilities.items()})
 np.savez_compressed(folder/'ranker_inputs.npz',**{k:v for k,v in d.items() if k not in ('pair_names','target_names')},**{'outer__'+k:v for k,v in outer.items()},
  **{f'inner{j}__{k}':v for j,c in enumerate(cases) for k,v in c.items()})
 writej(folder/'selection.json',dict(rank=meta,stable=state,ridge_selected=selected,ridge_grid=records))
 writej(folder/'FREEZE.json',dict(status='SMOKE_ONLY' if smoke else 'PREDICTIONS_FROZEN_UNSCORED',fold=f,design=design,train=train.tolist(),query=query.tolist(),
  source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'common.py',ROOT/'data/moments.npz']},seconds=time.monotonic()-t,
  output_sha256=sha(folder/'predictions.npz'),outer_scores_read=False))
 print(json.dumps(dict(fold=f,status='SMOKE_COMPLETE' if smoke else 'FROZEN',seconds=time.monotonic()-t)),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--fold',type=int);ap.add_argument('--design',default='detected49');ap.add_argument('--smoke',action='store_true');a=ap.parse_args()
 for f in ([a.fold] if a.fold is not None else range(5)):run(f,a.design,a.smoke)
