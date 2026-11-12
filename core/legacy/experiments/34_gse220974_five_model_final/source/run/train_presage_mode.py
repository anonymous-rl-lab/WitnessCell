"""Published ComboPRESAGE with an explicit gene-by-mode direct-sum input adapter."""
from common import *
import contextlib,os,time,types,shutil
import torch,pandas as pd,anndata
tp=module('original_presage_runner',PK/'train_presage.py')
PR=ROOT/'priors/presage';SOURCE=PK/'source/presage_probe'
tp.SRC=SOURCE
def features(base,fit,out,tag):
 out.mkdir(parents=True,exist_ok=True);m=base['m'];full=base['full'];Y=np.r_[m['single_effects'],m['means'][base['pair_moment_ids'][fit]]]
 labels=list(m['target_names'].astype(str))+[full.pairs[i].replace('+','_') for i in fit]
 ad=anndata.AnnData(X=Y,obs=pd.DataFrame({'perturbation':labels},index=labels),var=pd.DataFrame(index=m['genes'].astype(str)))
 dm=types.SimpleNamespace(batch_size=32,train_dataset=types.SimpleNamespace(adata=ad,X=Y),perturb_field='perturbation',split_path=tag+'.json')
 cfg=tp.CONFIG.copy();cfg['pathway_files']=str((PR/'knowledge_5006.txt').resolve());np.random.seed(24)
 here=Path.cwd()
 try:
  os.chdir(PR)
  with (out/'features.log').open('w') as log,contextlib.redirect_stdout(log):emb=tp.official_module().PrepareInputs(dm,cfg)._prep_inputs()
 finally:os.chdir(here)
 assert emb.shape==(len(m['genes']),128,42) and np.isfinite(emb).all()
 norms=np.median(np.linalg.norm(emb,axis=1,keepdims=True),axis=(0,1),keepdims=True);norms[norms==0]=1
 normalized=(emb/norms).astype(np.float32);small=normalized[m['target_indices']]
 a=np.array([str(n).endswith('-a') for n in m['target_names']])[:,None,None]
 expanded=np.concatenate([small*a,small*(~a)],axis=1)
 np.savez_compressed(out/'features.npz',normalized_targets=expanded,normalizers=norms,fit_ids=fit,fit_names=np.array(full.pairs)[fit],target_names=m['target_names'])
 writej(out/'FEATURES.json',dict(sources=42,physical_feature_dimension=128,mode_feature_dimension=256,fit_ids=fit.tolist(),fit_names=np.array(full.pairs)[fit].tolist(),query_outcomes_used=False))
 return expanded
def run(f,smoke=False):
 torch.set_num_threads(2);b=setup();m=b['m'];full=b['full'];P=len(full.pairs)
 q=folds(range(P),full.pairs,k=5,tag='outer')[f];tr=np.setdiff1d(range(P),q);inner=folds(tr,full.pairs,tag='inner'+str(f));val=inner[0];fit=np.setdiff1d(tr,val)
 out=ROOT/('smoke_neural' if smoke else 'results/detected49')/'presage'/f'fold_{f}';out.mkdir(parents=True,exist_ok=True)
 if (out/'FREEZE.json').exists():print('already frozen',f,flush=True);return
 # The archived fitter consumes local outer-training indices; query labels absent.
 d=dict(single_effects=m['single_effects'],train_y=m['means'][b['pair_moment_ids'][tr]],train_names=np.array(full.pairs)[tr],train_u=full.x[tr])
 c={k:m[k] for k in ('genes','target_names','target_indices','control')};order=np.arange(len(full.nodes));local={int(i):j for j,i in enumerate(tr)}
 fit_local=np.array([local[int(i)] for i in fit]);val_local=np.array([local[int(i)] for i in val])
 selection=out/'selection';small=features(b,fit,selection,f'new_f{f}_inner0')
 model,info=tp.fit_network(c,d,fit_local,order,small,selection,1000,val_ids=val_local,smoke_limit=1 if smoke else None)
 if smoke:
  model.eval();u=torch.tensor(full.x[val],dtype=torch.float32)
  with torch.no_grad():pred=model(u)[0]
  assert torch.isfinite(pred).all();writej(out/'SMOKE.json',dict(status='PASS',parameters=info['parameters'],pred_shape=list(pred.shape),outer_scores_read=False));return
 del model,small
 refit=out/'refit';small=features(b,tr,refit,f'new_f{f}_refit')
 model,refit_info=tp.fit_network(c,d,np.arange(len(tr)),order,small,refit,info['best_epoch'])
 model.eval()
 with torch.no_grad():pred=model(torch.tensor(full.x[q],dtype=torch.float32))[0].numpy()
 mask=m['score_gene_mask'];assert np.isfinite(pred).all()
 np.savez_compressed(out/'predictions.npz',query=q,conditions=np.array(full.pairs)[q],genes=m['genes'],pred_effect=pred,target_names=m['target_names'],utility=pred@m['objectives'].T)
 writej(out/'FREEZE.json',dict(status='FROZEN_UNSCORED',model='PRESAGE_Combo_mode_adapter_seed24',fold=f,selected_epoch=info['best_epoch'],
  source_model_sha256=sha(SOURCE/'upstream/src/presage.py'),source_wrapper_sha256=sha(SOURCE/'presage_wrapper.py'),adapter_sha256=sha(__file__),
  predictions_sha256=sha(out/'predictions.npz'),training_ids=tr.tolist(),query_ids=q.tolist(),outer_scores_read=False))
 print(json.dumps(dict(model='PRESAGE',fold=f,status='FROZEN')),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--fold',type=int);p.add_argument('--smoke',action='store_true');a=p.parse_args()
 for f in ([a.fold] if a.fold is not None else range(5)):run(f,a.smoke)
