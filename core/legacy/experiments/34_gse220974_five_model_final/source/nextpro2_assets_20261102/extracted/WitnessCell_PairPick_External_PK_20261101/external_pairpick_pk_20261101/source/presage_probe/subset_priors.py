"""Memory-map official NumPy/Pandas pickle payloads before original row filter."""
from pathlib import Path
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import gc,hashlib,json,mmap,pickle,resource,time
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
class MappedReader:
 def __init__(self,path):
  self.file=open(path,'rb');self.mm=mmap.mmap(self.file.fileno(),0,access=mmap.ACCESS_READ);self.pos=0
 def read(self,n=-1):
  if n<0:n=len(self.mm)-self.pos
  start=self.pos;self.pos+=n
  return memoryview(self.mm)[start:self.pos] if n>1024*1024 else self.mm[start:self.pos]
 def readline(self):
  end=self.mm.find(b'\n',self.pos);end=len(self.mm) if end<0 else end+1
  out=self.mm[self.pos:end];self.pos=end;return out
 def close(self):self.mm.close();self.file.close()
class ZeroCopyUnpickler(pickle._Unpickler):
 dispatch=pickle._Unpickler.dispatch.copy()
 def load_build(self):
  state=self.stack[-1];inst=self.stack[-2]
  if isinstance(inst,np.ndarray) and isinstance(state,tuple) and len(state)==5 and isinstance(state[-1],memoryview):
   _,shape,dtype,fortran,buf=state
   if dtype.hasobject:raise ValueError('large object arrays cannot be zero copied')
   array=np.ndarray(shape=shape,dtype=dtype,buffer=buf,order='F' if fortran else 'C')
   self.stack.pop();self.stack[-1]=array
   for key,value in list(self.memo.items()):
    if value is inst:self.memo[key]=array
  else:pickle._Unpickler.load_build(self)
 dispatch[pickle.BUILD[0]]=load_build

def mapped_load(p):
 reader=MappedReader(p);obj=ZeroCopyUnpickler(reader).load();return obj,reader

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def values_hash(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def main():
 out=ROOT/'subset_5020';out.mkdir(exist_ok=True)
 # Regression against the standard loader, including huge-byte opcode branch.
 for protocol in [4,5]:
  test=pd.DataFrame(np.arange(600000,dtype='float64').reshape(2000,300),index=[f'g{i}' for i in range(2000)])
  path=out/f'loader_test_{protocol}.pkl';test.to_pickle(path,protocol=protocol)
  obs,reader=mapped_load(path);pd.testing.assert_frame_equal(obs,pd.read_pickle(path));del obs;gc.collect();reader.close();path.unlink()
 common=ROOT.parents[2]/'four_model_pk_20261101/adapters/sf12_original/inputs/common.npz'
 with np.load(common,allow_pickle=False) as z:genes=z['genes'].astype(str)
 srcfile=ROOT/'upstream/sample_files/prior_files/sample.knowledge_experimental.txt'
 sources=[x.strip() for x in srcfile.read_text().splitlines() if x.strip()]
 manifest=[];start=time.monotonic()
 for source in sources:
  p=ROOT/source.replace('../','',1);dest=out/p.name;t=time.monotonic()
  original,reader=mapped_load(p)
  shared=np.intersect1d(original.index,genes)
  sub=original.loc[shared,:].copy()
  expected_hash=values_hash(sub.to_numpy());idx=list(map(str,sub.index));shape=list(original.shape)
  sub.to_pickle(dest,protocol=4)
  restored=pd.read_pickle(dest);pd.testing.assert_frame_equal(sub,restored,check_exact=True)
  assert values_hash(restored.to_numpy())==expected_hash
  entry=dict(source=source,source_shape=shape,subset_shape=list(sub.shape),subset=str(dest.relative_to(ROOT)),subset_sha256=sha(dest),values_sha256=expected_hash,index_sha256=hashlib.sha256('\n'.join(idx).encode()).hexdigest(),exact_reload=True,seconds=time.monotonic()-t,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
  manifest.append(entry);print(json.dumps(entry),flush=True)
  del original,sub,restored;gc.collect();reader.close()
  (out/'MANIFEST.json').write_text(json.dumps(dict(status='RUNNING',files=manifest),indent=2)+'\n')
 (out/'MANIFEST.json').write_text(json.dumps(dict(status='PASS',loader_regression_protocols=[4,5],genes_sha256=hashlib.sha256('\n'.join(genes).encode()).hexdigest(),files=manifest,seconds=time.monotonic()-start,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),indent=2)+'\n')
if __name__=='__main__':main()
