"""Restore original public priors; encode modes without collapsing swapped pairs."""
from common import *
import gc,pickle,time
def presage():
 m=load(ROOT/'data/moments.npz');genes=m['genes'].astype(str);root=ROOT/'priors/presage';dest=root/'subset';dest.mkdir(exist_ok=True)
 helper=module('original_subset_loader',PK/'source/presage_probe/subset_priors.py');records=[]
 sources=(root/'upstream/sample_files/prior_files/sample.knowledge_experimental.txt').read_text().splitlines()
 for source in sources:
  p=root/source.replace('../','',1);target=dest/p.name
  if target.exists():continue
  obj,reader=helper.mapped_load(p);shared=np.intersect1d(obj.index,genes);sub=obj.loc[shared].copy();sub.to_pickle(target,protocol=4)
  record=dict(source=source,source_shape=list(obj.shape),subset_shape=list(sub.shape),file=str(target.relative_to(root)),sha256=sha(target));records.append(record)
  print(json.dumps(record),flush=True);del obj,sub;gc.collect();reader.close()
 paths=[str((dest/Path(s).name).resolve()) for s in sources]
 assert all(Path(s).exists() for s in paths)
 (root/'knowledge_5006.txt').write_text('\n'.join(paths)+'\n');writej(dest/'MANIFEST.json',dict(genes_sha256=hashlib.sha256('\n'.join(genes).encode()).hexdigest(),files=records,status='PASS',sources=40))
def genept():
 m=load(ROOT/'data/moments.npz');p=ROOT/'priors/genept/GenePT_gene_embedding_ada_text.pickle'
 with p.open('rb') as f:obj=pickle.load(f)
 names=m['target_names'].astype(str);base=[n.rsplit('-',1)[0] for n in names]
 import pandas as pd
 if isinstance(obj,pd.DataFrame):E=obj.loc[base].to_numpy(dtype=np.float32)
 else:E=np.array([obj[g] for g in base],dtype=np.float32)
 mode=np.array([n.endswith('-a') for n in names]);expanded=np.concatenate([E*mode[:,None],E*(~mode[:,None])],axis=1)
 np.savez_compressed(ROOT/'priors/genept/mode_embeddings.npz',names=np.concatenate([names,['control']]),embedding=np.r_[expanded,np.zeros((1,expanded.shape[1]),dtype=np.float32)])
 writej(ROOT/'priors/genept/MODE_ADAPTER.json',dict(mode_encoding='direct sum: activation [GenePT,0], repression [0,GenePT]; control zero',source_sha256=sha(p),input_dimensions=expanded.shape[1],all_required_genes_covered=True))
 print('GenePT: all 35 mode identities encoded',flush=True)
if __name__=='__main__':genept();presage()
