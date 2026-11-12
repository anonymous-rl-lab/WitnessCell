"""Original GEARS / Scouter networks and losses on fresh mode-aware cell inputs."""
from common import *
import argparse,random,time,types,pickle
import h5py,torch,pandas as pd
from torch.utils.data import TensorDataset,DataLoader
ORIG=PK/'assets/original_run'
GM=module('original_gears_model',ORIG/'gears/vendor/source/gears/model.py')
GU=module('original_gears_utils',ORIG/'gears/vendor/source/gears/utils.py')
SM=module('original_scouter_model',ORIG/'scouter/upstream/scouter/_model.py')
SU=module('original_scouter_loss',ORIG/'scouter/upstream/scouter/_utils.py')
def prepare_graph():
 out=ROOT/'priors/gears';out.mkdir(exist_ok=True);path=out/'graphs.pt'
 if path.exists():return torch.load(path,weights_only=False)
 start=time.monotonic();m=load(ROOT/'data/moments.npz');assets=ORIG/'gears/assets'
 gene2go=pickle.load(open(assets/'gene2go_all.pkl','rb'));essential=pickle.load(open(assets/'essential_all_data_pert_genes.pkl','rb'))
 perts=np.unique([p for p in essential if p in gene2go]);lookup={g:i for i,g in enumerate(perts)}
 required={str(t).rsplit('-',1)[0] for t in m['target_names']};assert required<=set(perts),required-set(perts)
 go=pd.read_csv(assets/'go_essential_all/go_essential_all_top21.csv');net=GU.GeneSimNetwork(go,perts,lookup)
 G=torch.cat([net.edge_index,net.edge_index+len(perts)],dim=1);w=torch.cat([net.edge_weight,net.edge_weight])
 with h5py.File(ROOT/'data/cells.h5') as f:
  labels=f['condition_names'].asstr()[:];ci=f['condition_index'][:];ids=np.flatnonzero(np.array(['|' not in s for s in labels])[ci]);X=f['X'][ids].astype(np.float64)
 corr=np.abs(GU.np_pearson_cor(X,X));del X;corr[~np.isfinite(corr)]=0
 indices=np.argsort(corr,axis=1)[:,-21:];rows=np.repeat(np.arange(len(m['genes'])),21);cols=indices.ravel();values=corr[rows,cols];keep=values>.4
 frame=pd.DataFrame(dict(source=m['genes'][cols[keep]],target=m['genes'][rows[keep]],importance=values[keep]));co=GU.GeneSimNetwork(frame,m['genes'],{str(g):i for i,g in enumerate(m['genes'])})
 typed=np.r_[np.char.add(perts.astype(str),'-a'),np.char.add(perts.astype(str),'-i')]
 payload=dict(perts=typed,G_go=G,G_go_weight=w,G_coexpress=co.edge_index,G_coexpress_weight=co.edge_weight,genes=m['genes'])
 torch.save(payload,path);writej(out/'GRAPH.json',dict(status='PASS',mode_rule='two disjoint copies of original GO graph',physical_genes=len(m['genes']),perturbation_identities=len(typed),
  go_edges=int(G.shape[1]),coexpression_edges=int(co.edge_index.shape[1]),controls_and_singles_only=True,seconds=time.monotonic()-start))
 return payload
def config(graph):
 return dict(hidden_size=64,num_go_gnn_layers=1,num_gene_gnn_layers=1,decoder_hidden_size=16,uncertainty=False,uncertainty_reg=1,direction_lambda=.1,
  G_go=graph['G_go'],G_go_weight=graph['G_go_weight'],G_coexpress=graph['G_coexpress'],G_coexpress_weight=graph['G_coexpress_weight'],device='cpu',num_genes=len(graph['genes']),num_perts=len(graph['perts']),no_perturb=False)
class Data:
 def __init__(self,fold,kind):
  self.base=setup();b=self.base;self.m=b['m'];self.full=b['full'];self.kind=kind
  P=len(self.full.pairs);self.q=folds(range(P),self.full.pairs,k=5,tag='outer')[fold];self.tr=np.setdiff1d(range(P),self.q)
  self.val=folds(self.tr,self.full.pairs,tag='inner'+str(fold))[0];self.fit=np.setdiff1d(self.tr,self.val)
  self.allowed=['control']+self.m['target_names'].astype(str).tolist()+np.array(self.full.pairs)[self.tr].tolist()
  with h5py.File(ROOT/'data/cells.h5') as h:
   names=np.char.replace(h['condition_names'].asstr()[:].astype(str),'|','+');ci=h['condition_index'][:];ii=np.flatnonzero(np.isin(names[ci],self.allowed))
   self.X=h['X'][ii];self.labels=names[ci[ii]]
  self.ctrl=self.X[self.labels=='control'];self.ctrlmean=self.ctrl.mean(0,dtype=np.float64)
  perm=np.random.default_rng(20261018).permutation(len(self.ctrl));self.valctrl=self.ctrl[perm[:64]];self.testctrl=self.ctrl[perm[:300]]
  self.emb=None;self.graph=None
  if kind=='gears':self.graph=prepare_graph();self.mapping={n:i for i,n in enumerate(self.graph['perts'])}
  else:
   z=load(ROOT/'priors/genept/mode_embeddings.npz');self.emb=torch.tensor(z['embedding']);self.mapping={n:i for i,n in enumerate(z['names'])}
 def names(self,ids):return np.array(self.full.pairs)[ids].tolist()
 def newmodel(self):
  if self.kind=='gears':return GM.GEARS_Model(config(self.graph))
  return SM.ScouterModel(len(self.m['genes']),self.emb,n_encoder=(2048,512),n_out_encoder=64,n_decoder=(2048,),use_batch_norm=True,use_layer_norm=False,dropout_rate=0.)
 def indices(self,condition):
  if condition=='control':return [-1] if self.kind=='gears' else [self.mapping['control'],self.mapping['control']]
  p=[self.mapping[g] for g in condition.split('+')]
  if self.kind=='scouter' and len(p)==1:p.append(self.mapping['control'])
  return p
 def train_arrays(self,ids,seed):
  names=self.m['target_names'].astype(str).tolist()+self.names(ids)
  if self.kind=='gears':names+=['control']
  ii=np.flatnonzero(np.isin(self.labels,names));Y=self.X[ii];labels=self.labels[ii];N=len(ii)
  if self.kind=='gears':
   rng=np.random.default_rng(seed);ctrl=self.ctrl[rng.integers(0,len(self.ctrl),size=N)].copy();sel=labels=='control';ctrl[sel]=Y[sel]
  else:
   # Same fixed BalancedDataset mapping, without materializing a dense CSR copy.
   rng=np.random.RandomState(24);pp=np.arange(N);cc=np.resize(np.arange(len(self.ctrl)),N);rng.shuffle(pp);rng.shuffle(cc)
   Y=Y[pp];labels=labels[pp];ctrl=self.ctrl[cc]
  pi=[self.indices(n) for n in labels];nz={n:np.flatnonzero(Y[labels==n].mean(0)!=0) for n in names if n!='control'}
  return Y,labels,ctrl,pi,nz
# Memory layout and graph operations only; original model/loss/fit rules retained.
from lazy_cell_data import LazyData
import fast_gears_ops
Data=LazyData
_original_newmodel=Data.newmodel
def _runtime_newmodel(self):
 model=_original_newmodel(self)
 if self.kind=='gears':fast_gears_ops.enable(model)
 return model
Data.newmodel=_runtime_newmodel

def runtime_fingerprint():
 return {p:sha(ROOT/p) for p in ['train_cell_models.py','train_cell_models_runtime.py','lazy_cell_data.py','fast_gears_ops.py']}

def batch(X,pi):
 n,g=X.shape
 return types.SimpleNamespace(x=torch.as_tensor(X).reshape(-1,1),pert_idx=pi,batch=torch.arange(n).repeat_interleave(g))
@torch.no_grad()
def predict(model,d,ids,validation=False):
 model.eval();ctrl=d.valctrl if validation else d.testctrl;output=[]
 for name in d.names(ids):
  pi=d.indices(name)
  if d.kind=='gears':
   pp=model(batch(ctrl[:2],[pi]*2)).numpy().astype(float);delta=pp-ctrl[:2];n=len(ctrl)
   p=ctrl.mean(0,dtype=np.float64)+delta[0]/n+delta[1]*(n-1)/n
  else:
   vals=[]
   for start in range(0,len(ctrl),64):
    c=torch.from_numpy(ctrl[start:start+64]);vals.append(model(torch.tensor([pi]*len(c)),c).numpy())
   p=np.concatenate(vals).mean(0,dtype=np.float64)
  output.append(p-d.ctrlmean)
 return np.array(output)
def validation_key(pred,d,epoch):
 ids=d.val;truth=d.m['means'][d.base['pair_moment_ids'][ids]];mask=d.m['score_gene_mask'];err=pred[:,mask]-truth[:,mask]
 lookup={int(i):j for j,i in enumerate(ids)};loss=[]
 for block in gm.available_blocks(d.full,ids):
  _,sign=gm.block_contrast(np.zeros((len(d.full.pairs),1)),block,d.full.x);loss.append(float(np.mean((sign@err[[lookup[int(i)] for i in block]])**2)))
 assert loss,'No inner rectangle; do not substitute outer metric'
 return (float(np.mean(loss)),float(np.mean(err*err)),epoch)
def fit(d,ids,seed,out,epochs,selection=False,schedule=5,smoke=False):
 out.mkdir(parents=True,exist_ok=True);random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
 model=d.newmodel();opt=torch.optim.Adam(model.parameters(),lr=.001,weight_decay=.0005 if d.kind=='gears' else 0.)
 sched=torch.optim.lr_scheduler.StepLR(opt,step_size=schedule,gamma=.5) if d.kind=='gears' else torch.optim.lr_scheduler.ExponentialLR(opt,gamma=.9)
 Y,labels,ctrl,pi,nz=d.train_arrays(ids,seed);N=len(Y);bs=32 if d.kind=='gears' else 256
 rng=np.random.default_rng(seed);history=[];best=(float('inf'),float('inf'),0);best_epoch=0;stale=0;start_epoch=0;cap=epochs;extended=False
 ckpath=out/'resume.pt';partial=None
 if (out/'RUN.json').exists():
  info=json.loads((out/'RUN.json').read_text())
  if info['status']=='COMPLETE':
   ck=torch.load(ckpath,weights_only=False)
   assert ck['runtime_fingerprint']==runtime_fingerprint()
   model.load_state_dict(torch.load(out/'last.pt',weights_only=False));return model,info

 if ckpath.exists():
  ck=torch.load(ckpath,weights_only=False);assert ck['runtime_fingerprint']==runtime_fingerprint(),'Changed source cannot resume silently'
  model.load_state_dict(ck['model']);opt.load_state_dict(ck['optimizer']);sched.load_state_dict(ck['scheduler']);torch.set_rng_state(ck['torch_rng']);rng.bit_generator.state=ck['numpy_rng']
  history=ck['history'];start_epoch=ck['epoch'];best=tuple(ck['best']);best_epoch=ck['best_epoch'];stale=ck['stale'];cap=ck['cap'];extended=ck['extended'];partial=ck.get('partial')
 t=time.monotonic();epoch=start_epoch
 def checkpoint(completed_epoch,partial_state=None):
  tmp=out/'resume.pt.tmp'
  torch.save(dict(model=model.state_dict(),optimizer=opt.state_dict(),scheduler=sched.state_dict(),torch_rng=torch.get_rng_state(),numpy_rng=rng.bit_generator.state,history=history,
   epoch=completed_epoch,best=best,best_epoch=best_epoch,stale=stale,cap=cap,extended=extended,adapter_sha256=sha(__file__),runtime_fingerprint=runtime_fingerprint(),partial=partial_state),tmp)
  tmp.replace(ckpath)

 while epoch<cap:
  epoch+=1;model.train();tic=time.monotonic();losses=[]
  if partial is None:
   order=rng.permutation(N) if d.kind=='gears' else torch.randperm(N).numpy();start_batch=0
  else:
   order=partial['order'];start_batch=partial['next_start'];losses=partial['losses'];partial=None
  for start in range(start_batch,N-bs+1,bs):
   ix=order[start:start+bs];opt.zero_grad(set_to_none=True);x=torch.from_numpy(ctrl[ix]);y=torch.from_numpy(Y[ix])
   if d.kind=='gears':
    pred=model(batch(x,[pi[i] for i in ix]));lab=np.where(labels[ix]=='control','ctrl',labels[ix]);loss=GU.loss_fct(pred,y,lab,torch.tensor(d.ctrlmean,dtype=torch.float32),.1,nz)
   else:
    pred=model(torch.tensor([pi[i] for i in ix]),x);loss=SU.gears_loss(pred,y,x,labels[ix],nz,gamma=0.,lambda_=.05)
   if not torch.isfinite(loss):raise FloatingPointError('nonfinite loss')
   loss.backward()
   if not all(torch.isfinite(p.grad).all() for p in model.parameters() if p.grad is not None):raise FloatingPointError('nonfinite gradient')
   if d.kind=='gears':torch.nn.utils.clip_grad_value_(model.parameters(),1.)
   else:torch.nn.utils.clip_grad_norm_(model.parameters(),1.)
   opt.step();losses.append(float(loss.detach()))
   if not smoke and len(losses) % (100 if d.kind=='gears' else 25)==0:
    checkpoint(epoch-1,dict(order=order,next_start=start+bs,losses=losses))
   if not smoke and (len(losses)==1 or len(losses)%50==0):
    print(json.dumps(dict(model=d.kind,run=str(out.relative_to(ROOT)),epoch=epoch,batch=len(losses),batches_per_epoch=N//bs,elapsed_seconds=time.monotonic()-tic,loss=losses[-1])),flush=True)
   if smoke and len(losses)>=3:break
  sched.step();key=None
  if selection:
   key=validation_key(predict(model,d,d.val,True),d,epoch)
   if key<best:best=key;best_epoch=epoch;stale=0;torch.save(model.state_dict(),out/'best.pt')
   else:stale+=1
  row=dict(epoch=epoch,loss=float(np.mean(losses)),steps=len(losses),seconds=time.monotonic()-tic,validation_C=None if key is None else key[0],validation_RNA=None if key is None else key[1],best_epoch=best_epoch,lr=opt.param_groups[0]['lr'])
  history.append(row);writej(out/'curve.json',history);print(json.dumps(dict(model=d.kind,run=str(out.relative_to(ROOT)),**row)),flush=True)
  if selection and d.kind=='scouter' and epoch==40 and len(history)>=40:
   early=min(x['validation_C'] for x in history[:-5]);late=min(x['validation_C'] for x in history[-5:])
   if best_epoch>35 and late<.99*early:cap=80;extended=True
  checkpoint(epoch)
  if smoke or (selection and d.kind=='gears' and stale>=10):break
 torch.save(model.state_dict(),out/'last.pt')
 info=dict(status='SMOKE' if smoke else 'COMPLETE',epochs=epoch,selected_epoch=best_epoch if selection else epoch,best=None if not selection else list(best),seconds=time.monotonic()-t,parameters=sum(p.numel() for p in model.parameters()),extended=extended,training_cells=N)
 writej(out/'RUN.json',info);return model,info
def run(f,kind,seed,smoke=False,fixed40=False):
 torch.set_num_threads(2);d=Data(f,kind);name=kind+('_fixed40' if fixed40 else '')+'_seed'+str(seed)
 out=ROOT/('smoke_neural' if smoke else 'results/detected49')/name/f'fold_{f}';out.mkdir(parents=True,exist_ok=True)
 if (out/'FREEZE.json').exists():print('already frozen',name,f,flush=True);return
 if smoke:
  model,info=fit(d,d.fit,seed,out/'selection',1,selection=True,smoke=True)
  if kind=='gears':
   model.eval();cc=d.valctrl;pi=d.indices(d.names(d.val)[0])
   with torch.no_grad():direct=model(batch(cc,[pi]*len(cc))).numpy().mean(0,dtype=np.float64)-d.ctrlmean
   fast=predict(model,d,d.val[:1],True)[0];assert np.max(abs(direct-fast))<5e-5
  writej(out/'SMOKE.json',{**info,'status':'PASS','finite_loss_and_gradients':True,'outer_scores_read':False});return
 if fixed40:epochs=40;schedule=5
 elif kind=='gears':
  selected=[]
  for schedule,cap in ((1,20),(5,60)):
   model,info=fit(d,d.fit,seed,out/f'select_schedule{schedule}',cap,selection=True,schedule=schedule);selected.append(dict(schedule=schedule,**info));del model
  chosen=min(selected,key=lambda x:(x['best'][0],x['best'][1],x['schedule'],x['selected_epoch']));epochs=chosen['selected_epoch'];schedule=chosen['schedule'];writej(out/'SELECTION.json',dict(selected=chosen,candidates=selected))
 else:
  model,info=fit(d,d.fit,seed,out/'selection',40,selection=True);epochs=info['selected_epoch'];schedule=5;del model;writej(out/'SELECTION.json',info)
 model,info=fit(d,d.tr,seed,out/'refit',epochs,schedule=schedule);pred=predict(model,d,d.q)
 assert np.isfinite(pred).all()
 np.savez_compressed(out/'predictions.npz',query=d.q,conditions=np.array(d.full.pairs)[d.q],genes=d.m['genes'],pred_effect=pred,utility=pred@d.m['objectives'].T,target_names=d.m['target_names'])
 writej(out/'FREEZE.json',dict(status='FROZEN_UNSCORED',model=name+'_mode_adapter',fold=f,seed=seed,epochs=epochs,schedule=schedule,adapter_sha256=sha(__file__),
  runtime_fingerprint=runtime_fingerprint(),model_source_sha256=sha(ORIG/('gears/vendor/source/gears/model.py' if kind=='gears' else 'scouter/upstream/scouter/_model.py')),predictions_sha256=sha(out/'predictions.npz'),train=d.tr.tolist(),query=d.q.tolist(),outer_scores_read=False))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['gears','scouter'],required=True);ap.add_argument('--fold',type=int,default=0);ap.add_argument('--seed',type=int,default=0);ap.add_argument('--smoke',action='store_true');ap.add_argument('--fixed40',action='store_true');a=ap.parse_args();run(a.fold,a.kind,a.seed,a.smoke,a.fixed40)
