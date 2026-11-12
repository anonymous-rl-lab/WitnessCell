"""Faithful CellSimBench ComboPRESAGE adapter, fixed split and epoch selection."""
from pathlib import Path
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,ast,contextlib,copy,hashlib,json,sys,time,types
import numpy as np,pandas as pd,torch,anndata
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'source/presage_probe'
INPUTS=ROOT.parent/'four_model_pk_20261101/adapters/sf12_original/inputs'
SEED=24
CONFIG=dict(item_hidden_size=512,item_nlayers=0,pathway_item_hidden_size=128,pathway_item_nlayers=2,
 pathway_pool_type='sum',pathway_weight_type='gat',pool_nlayers=2,softmax_temperature=.15,gat_weight=.85,
 n_nmf_embedding=128,pca_dim=None,learnable_gene_embedding=False,batch_norm=True,
 pathway_embedding_dimension=5,pathway_embedding_files='./cache/pathway_embeddings/',
 node2vec_walk_length=100,node2vec_context_size=5,node2vec_walks_per_node=3,node2vec_num_negative_samples=3,
 node2vec_p=1.,node2vec_q=2.,node2vec_batchsize=32,dim_red_alg='Node2Vec',n_neigh_prune='5',
 dataset='cellsimbench',source='cellsimbench',input_preparation='prep_gene_embeddings',embedding_files='None',
 mse_loss_scale=1.,cosine_loss_scale=1.,vector_norm_loss_scale=1.,num_heads=4,
 lr=.0024,weight_decay=2.17e-13,momentum=.9,optimizer='Adam')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writej(p,x):Path(p).write_text(json.dumps(x,indent=2,sort_keys=True,allow_nan=False)+'\n')
def inputs(fold):
 with np.load(INPUTS/'common.npz',allow_pickle=False) as z:c={k:z[k].copy() for k in z.files}
 with np.load(INPUTS/f'fold_{fold}.npz',allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
 assert not set(d['query_names'])&set(d['train_names'])
 return c,d

def official_module():
 sys.path.insert(0,str(SRC/'upstream/src'))
 import presage
 return presage

def combo_class(prep):
 upstream=official_module();source=(SRC/'presage_wrapper.py').read_text()
 node=next(n for n in ast.parse(source).body if isinstance(n,ast.ClassDef) and n.name=='ComboPRESAGE')
 ns=dict(torch=torch,nn=torch.nn,np=np,PrepareInputs=prep,GeneEmbeddingTransformation=upstream.GeneEmbeddingTransformation,Pool=upstream.Pool,ItemNet=upstream.ItemNet)
 exec(compile(ast.Module(body=[node],type_ignores=[]),str(SRC/'presage_wrapper.py'),'exec'),ns)
 return ns['ComboPRESAGE']

def prepare_features(c,d,fit_ids,tag,out):
 """Only fit response rows plus singles reach original PrepareInputs."""
 fit_ids=np.asarray(fit_ids,int);Y=np.concatenate([d['single_effects'],d['train_y'][fit_ids]])
 labels=list(map(str,c['target_names']))+[str(x).replace('+','_') for x in d['train_names'][fit_ids]]
 assert np.isfinite(Y).all()
 ad=anndata.AnnData(X=Y.copy(),obs=pd.DataFrame({'perturbation':labels},index=labels),var=pd.DataFrame(index=c['genes'].astype(str)))
 dm=types.SimpleNamespace(batch_size=32,train_dataset=types.SimpleNamespace(adata=ad,X=Y),perturb_field='perturbation',split_path=tag+'.json')
 genehash=hashlib.sha256('\n'.join(c['genes'].astype(str)).encode()).hexdigest()
 prior=SRC/f'required_knowledge_5020_{genehash[:12]}.txt'
 sources=(SRC/'upstream/sample_files/prior_files/sample.knowledge_experimental.txt').read_text().splitlines()
 subset_manifest=json.loads((SRC/'subset_5020/MANIFEST.json').read_text())
 assert subset_manifest['status']=='PASS' and subset_manifest['genes_sha256']==genehash
 prior.write_text('\n'.join(str((SRC/'subset_5020'/Path(x).name).resolve()) for x in sources if x.strip())+'\n')
 cfg=CONFIG.copy();cfg['pathway_files']=str(prior.resolve())
 np.random.seed(SEED);oldcwd=Path.cwd();start=time.monotonic()
 try:
  os.chdir(SRC)
  with (out/'features.log').open('w') as log,contextlib.redirect_stdout(log):emb=official_module().PrepareInputs(dm,cfg)._prep_inputs()
 finally:os.chdir(oldcwd)
 assert emb.shape==(len(c['genes']),128,42) and np.isfinite(emb).all(),emb.shape
 norms=np.linalg.norm(emb,axis=1,keepdims=True);norms=np.median(norms,axis=(0,1),keepdims=True);norms[norms==0]=1
 normalized=(emb/norms).astype('float32')
 order=np.argsort(c['target_indices']);rna_indices=c['target_indices'][order]
 target_normalized=normalized[rna_indices].copy()
 np.savez_compressed(out/'features.npz',normalized_targets=target_normalized,normalizers=norms,genes=c['genes'],target_names=c['target_names'][order],target_rna_indices=rna_indices,target_order=order,fit_ids=fit_ids,fit_names=d['train_names'][fit_ids],single_names=c['target_names'])
 writej(out/'FEATURES.json',dict(status='PASS',fit_names=list(map(str,d['train_names'][fit_ids])),single_names=list(map(str,c['target_names'])),outer_query_outcomes_used=False,source_count=42,source_sha256=sha(SRC/'upstream/sample_files/prior_files/sample.knowledge_experimental.txt'),seconds=time.monotonic()-start))
 return target_normalized,order,normalized

def make_model(normalized,c):
 class FixedInputs:
  def __init__(self,*args):pass
  def _prep_inputs(self):return normalized
 ad=types.SimpleNamespace(var=pd.DataFrame(index=c['genes'].astype(str)))
 dm=types.SimpleNamespace(batch_size=32,train_dataset=types.SimpleNamespace(adata=ad))
 model=combo_class(FixedInputs)(CONFIG.copy(),dm,len(c['genes']),len(c['genes']))
 # Preserve full-5020 normalization already applied by original initialization.
 model.gene_embeddings=torch.as_tensor(normalized,dtype=torch.float32).clone()
 model.mask=model.gene_embeddings.sum(dim=1,keepdim=True)!=0
 return model

def equivalence(c,d,small,order,full,out):
 torch.manual_seed(SEED);a=make_model(full,c);b=make_model(small,c);b.load_state_dict(a.state_dict())
 u=d['train_u'][:5];full_u=np.zeros((len(u),len(c['genes'])),dtype='float32');full_u[:,c['target_indices']]=u
 compact_u=u[:,order].astype('float32')
 checks=[]
 for mode in ['eval','train']:
  a.train(mode=='train');b.train(mode=='train');a.zero_grad(set_to_none=True);b.zero_grad(set_to_none=True)
  pa=a(torch.from_numpy(full_u))[0];pb=b(torch.from_numpy(compact_u))[0]
  np.testing.assert_allclose(pa.detach().numpy(),pb.detach().numpy(),rtol=1e-6,atol=1e-7)
  checks.append(dict(mode=mode,output_max_abs=float((pa-pb).abs().max().detach())))
  if mode=='train':
   target=torch.tensor(d['train_y'][:len(u)],dtype=torch.float32)
   ((pa-target)**2).mean().backward();((pb-target)**2).mean().backward()
   maxdiff=0.
   for (name,p),(name2,q) in zip(a.named_parameters(),b.named_parameters()):
    assert name==name2 and (p.grad is None)==(q.grad is None)
    if p.grad is None:continue
    torch.testing.assert_close(p.grad,q.grad,rtol=1e-5,atol=1e-7);maxdiff=max(maxdiff,float((p.grad-q.grad).abs().max()))
   checks[-1]['gradient_max_abs']=maxdiff
 writej(out/'INDEX_EQUIVALENCE.json',dict(status='PASS',checks=checks,query_outcomes_used=False))
 del a,b

def train_matrix(c,d,fit_ids,order):
 y=np.concatenate([d['single_effects'],d['train_y'][fit_ids]]).astype('float32')
 u=np.concatenate([np.eye(len(order)),d['train_u'][fit_ids]])[:,order].astype('float32')
 return torch.from_numpy(u),torch.from_numpy(y)

def fit_network(c,d,fit_ids,order,features,out,epochs,val_ids=None,smoke_limit=None):
 torch.manual_seed(SEED);np.random.seed(SEED)
 model=make_model(features,c);u,y=train_matrix(c,d,fit_ids,order)
 ds=torch.utils.data.TensorDataset(u,y)
 gen=torch.Generator().manual_seed(SEED)
 loader=torch.utils.data.DataLoader(ds,batch_size=32,shuffle=True,generator=gen,drop_last=False,num_workers=0)
 opt=torch.optim.Adam(model.parameters(),lr=CONFIG['lr'],weight_decay=CONFIG['weight_decay']);scheduler=torch.optim.lr_scheduler.StepLR(opt,step_size=200,gamma=1.)
 best=float('inf');best_epoch=0;wait=0;curve=[];start=time.monotonic();epoch_times=[]
 v_u=None if val_ids is None else torch.tensor(d['train_u'][val_ids][:,order],dtype=torch.float32)
 v_y=None if val_ids is None else torch.tensor(d['train_y'][val_ids],dtype=torch.float32)
 for epoch in range(1,epochs+1):
  t=time.monotonic();model.train();total=0.;count=0
  for ub,yb in loader:
   opt.zero_grad(set_to_none=True);pred=model(ub)[0];loss=model.compute_loss(pred,yb,None)
   if not torch.isfinite(loss):raise RuntimeError('nonfinite training loss')
   loss.backward();opt.step();total+=float(loss.detach())*len(ub);count+=len(ub)
  scheduler.step();value=None
  if val_ids is not None:
   model.eval()
   with torch.no_grad():value=float(torch.nn.functional.mse_loss(model(v_u)[0],v_y))
   if not np.isfinite(value):raise RuntimeError('nonfinite validation loss')
   if value<best:
    best=value;best_epoch=epoch;wait=0;torch.save(model.state_dict(),out/'best_state.pt')
   else:wait+=1
  elapsed=time.monotonic()-t;epoch_times.append(elapsed)
  row=dict(epoch=epoch,train_loss=total/count,validation_loss=value,epoch_seconds=elapsed,best_epoch=best_epoch if val_ids is not None else None,wait=wait)
  curve.append(row);writej(out/'curve.json',curve)
  if epoch<=3 or epoch%10==0 or wait>=10:print(json.dumps(dict(stage=out.name,**row)),flush=True)
  if smoke_limit is not None and epoch>=smoke_limit:break
  if val_ids is not None and wait>=10:break
 torch.save(model.state_dict(),out/'last_state.pt')
 report=dict(status='SMOKE_ONLY' if smoke_limit is not None else 'COMPLETE',epochs_completed=len(curve),best_epoch=best_epoch if val_ids is not None else len(curve),best_validation_loss=best if val_ids is not None else None,stop='smoke' if smoke_limit is not None else 'validation_patience' if val_ids is not None and wait>=10 else 'budget',seconds=time.monotonic()-start,parameters=sum(p.numel() for p in model.parameters()),outer_query_outcomes_used=False)
 writej(out/'RUN.json',report)
 return model,report

def run_fold(fold,outroot,smoke=False,prepare_only=False):
 c,d=inputs(fold);torch.set_num_threads(1);out=outroot.resolve()/f'fold_{fold}';out.mkdir(parents=True,exist_ok=True)
 if (out/'FREEZE.json').exists():raise FileExistsError('already frozen')
 selected=np.flatnonzero(d['inner_fold_id']==0);fit=np.flatnonzero(d['inner_fold_id']!=0)
 assert not set(selected)&set(fit) and len(selected)>0
 selection=out/'selection';selection.mkdir(exist_ok=True)
 small,order,full=prepare_features(c,d,fit,f'fold{fold}_inner0',selection)
 if fold==0:equivalence(c,d,small,order,full,selection)
 del full
 if prepare_only:return
 model,info=fit_network(c,d,fit,order,small,selection,1000,val_ids=selected,smoke_limit=1 if smoke else None)
 del model,small
 if smoke:return
 refit=out/'refit';refit.mkdir(exist_ok=True)
 small,order,full=prepare_features(c,d,np.arange(len(d['train_names'])),f'fold{fold}_refit',refit);del full
 model,refit_info=fit_network(c,d,np.arange(len(d['train_names'])),order,small,refit,info['best_epoch'])
 model.eval()
 with torch.no_grad():pred=model(torch.tensor(d['query_u'][:,order],dtype=torch.float32))[0].numpy()
 assert pred.shape==(len(d['query_names']),len(c['genes'])) and np.isfinite(pred).all()
 np.savez_compressed(out/'predictions.npz',conditions=d['query_names'],genes=c['genes'],control=c['control'],pred_effect=pred,pred_mean=pred+c['control'],fold=fold,seed=SEED)
 writej(out/'MODEL.json',dict(model='PRESAGE_CellSimBench_Combo',seed=SEED,fold=fold,best_epoch=info['best_epoch'],config=CONFIG,source_provenance=json.loads((SRC/'SOURCE_PROVENANCE.json').read_text()),query_names=list(map(str,d['query_names'])),train_names=list(map(str,d['train_names'])),selection_fit_names=list(map(str,d['train_names'][fit])),selection_validation_names=list(map(str,d['train_names'][selected])),outer_query_outcomes_used=False,control_added_once=True))
 names=['predictions.npz','MODEL.json','selection/features.npz','selection/RUN.json','selection/best_state.pt','refit/features.npz','refit/RUN.json','refit/last_state.pt']
 writej(out/'FREEZE.json',dict(status='FROZEN_BEFORE_EXTERNAL_SCORING',code_sha256=sha(__file__),protocol_sha256=sha(ROOT/'PRESAGE_TRAINING_PROTOCOL.txt'),sources={str(p.relative_to(SRC)):sha(p) for p in [SRC/'presage_wrapper.py',SRC/'upstream/src/presage.py']},files={n:sha(out/n) for n in names},outer_query_outcomes_used=False))
 print(json.dumps(dict(fold=fold,status='FROZEN',best_epoch=info['best_epoch'])),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('--fold',type=int,required=True);p.add_argument('--output',type=Path,default=ROOT/'presage');p.add_argument('--smoke',action='store_true');p.add_argument('--prepare-only',action='store_true');a=p.parse_args();run_fold(a.fold,a.output,a.smoke,a.prepare_only)
if __name__=='__main__':main()
