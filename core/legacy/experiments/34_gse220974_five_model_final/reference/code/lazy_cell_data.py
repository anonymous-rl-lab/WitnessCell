"""Index-only training views. Values and cell pairings match native Data exactly."""
from train_cell_models import *

class Rows:
 def __init__(self,source,indices):self.source=source;self.indices=np.asarray(indices)
 def __len__(self):return len(self.indices)
 def __getitem__(self,key):return np.asarray(self.source[self.indices[key]])

class LazyData(Data):
 def __init__(self,fold,kind):
  self.base=setup();b=self.base;self.m=b['m'];self.full=b['full'];self.kind=kind
  P=len(self.full.pairs);self.q=folds(range(P),self.full.pairs,k=5,tag='outer')[fold];self.tr=np.setdiff1d(range(P),self.q)
  self.val=folds(self.tr,self.full.pairs,tag='inner'+str(fold))[0];self.fit=np.setdiff1d(self.tr,self.val)
  self.allowed=['control']+self.m['target_names'].astype(str).tolist()+np.array(self.full.pairs)[self.tr].tolist()
  self.source=np.load(ROOT/'data/cells_X.npy',mmap_mode='r')
  index=load(ROOT/'data/cells_index.npz');names=np.char.replace(index['condition_names'].astype(str),'|','+');ci=index['condition_index']
  self.cell_ids=np.flatnonzero(np.isin(names[ci],self.allowed));self.labels=names[ci[self.cell_ids]]
  self.X=Rows(self.source,self.cell_ids);self.ctrl_ids=self.cell_ids[self.labels=='control'];self.ctrl=self.source[self.ctrl_ids]
  self.ctrlmean=self.ctrl.mean(0,dtype=np.float64)
  perm=np.random.default_rng(20261018).permutation(len(self.ctrl));self.valctrl=self.ctrl[perm[:64]];self.testctrl=self.ctrl[perm[:300]]
  self.emb=None;self.graph=None
  if kind=='gears':self.graph=prepare_graph();self.mapping={n:i for i,n in enumerate(self.graph['perts'])}
  else:
   z=load(ROOT/'priors/genept/mode_embeddings.npz');self.emb=torch.tensor(z['embedding']);self.mapping={n:i for i,n in enumerate(z['names'])}
 def train_arrays(self,ids,seed):
  names=self.m['target_names'].astype(str).tolist()+self.names(ids)
  if self.kind=='gears':names+=['control']
  ii=np.flatnonzero(np.isin(self.labels,names));labels=self.labels[ii];yi=self.cell_ids[ii];N=len(ii)
  if self.kind=='gears':
   rng=np.random.default_rng(seed);cc=self.ctrl_ids[rng.integers(0,len(self.ctrl_ids),size=N)].copy();sel=labels=='control';cc[sel]=yi[sel]
  else:
   rng=np.random.RandomState(24);pp=np.arange(N);ci=np.resize(np.arange(len(self.ctrl)),N);rng.shuffle(pp);rng.shuffle(ci)
   yi=yi[pp];labels=labels[pp];cc=self.ctrl_ids[ci]
  Y=Rows(self.source,yi);ctrl=Rows(self.source,cc)
  pi=[self.indices(n) for n in labels]
  nz={n:np.flatnonzero(Y[labels==n].mean(0)!=0) for n in names if n!='control'}
  return Y,labels,ctrl,pi,nz
