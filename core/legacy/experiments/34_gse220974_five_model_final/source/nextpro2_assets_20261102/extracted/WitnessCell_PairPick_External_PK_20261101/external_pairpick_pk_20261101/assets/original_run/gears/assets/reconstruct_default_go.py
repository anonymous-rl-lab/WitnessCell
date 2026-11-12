#!/usr/bin/env python3
"""Reconstruct default essential GO network by official GEARS make_GO rule.
Preserves essential-file iteration order and exact set Jaccard > 0.1.
This is a mathematical reconstruction, NOT claimed byte-identical downloaded CSV.
"""
from pathlib import Path
import csv,hashlib,json,pickle,time
import numpy as np
from scipy import sparse
P=Path(__file__).resolve().parent
start=time.monotonic()
gene2go=pickle.load((P/'gene2go_all.pkl').open('rb'))
essential=pickle.load((P/'essential_all_data_pert_genes.pkl').open('rb'))
# Official PertData filters the default set to entries in the GO dictionary.
genes=list(dict.fromkeys(str(g) for g in essential if g in gene2go))
terms=sorted(set().union(*(gene2go[g] for g in genes)));tidx={t:i for i,t in enumerate(terms)}
rows=[];cols=[]
for i,g in enumerate(genes):
 for t in gene2go[g]:rows.append(i);cols.append(tidx[t])
A=sparse.csr_matrix((np.ones(len(rows),dtype=np.int32),(rows,cols)),shape=(len(genes),len(terms)))
counts=np.asarray(A.sum(axis=1)).ravel();assert counts.min()>0
out=P/'go_essential_all'/'go_essential_all.csv';out.parent.mkdir(exist_ok=True)
count=0;verified=0
check=set(np.random.default_rng(20261018).choice(len(genes),100,replace=False))
with out.with_suffix('.tmp').open('w',newline='') as f:
 w=csv.writer(f,lineterminator='\n');w.writerow(['source','target','importance'])
 for lo in range(0,len(genes),128):
  inter=(A[lo:lo+128]@A.T).toarray()
  union=counts[lo:lo+128,None]+counts[None,:]-inter
  scores=inter/union
  for j in range(len(inter)):
   i=lo+j;idx=np.flatnonzero(scores[j]>.1)
   if i in check:
    literal=[(k,len(gene2go[genes[i]]&gene2go[g])/len(gene2go[genes[i]]|gene2go[g])) for k,g in enumerate(genes)]
    literal=[(k,s) for k,s in literal if s>.1]
    assert len(literal)==len(idx)
    assert all(k==idx[z] and s==scores[j,k] for z,(k,s) in enumerate(literal));verified+=1
   w.writerows((genes[i],genes[k],float(scores[j,k])) for k in idx);count+=len(idx)
  if lo%1280==0:print(f'rows {min(lo+128,len(genes))}/{len(genes)}, edges {count}, seconds {time.monotonic()-start:.2f}',flush=True)
out.with_suffix('.tmp').replace(out)
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
receipt={'algorithm':'Official GEARS make_GO exact Jaccard intersection/union; retain score > 0.1','upstream_rule_file':'../vendor/source/gears/utils.py:get_GO_edge_list and make_GO','essential_gene_order':'file order after GO support filtering and stable deduplication','genes':len(genes),'terms':len(terms),'edges':count,'literal_set_rows_verified':verified,'input_sha256':{p.name:h(p) for p in (P/'gene2go_all.pkl',P/'essential_all_data_pert_genes.pkl')},'output_sha256':h(out),'output_bytes':out.stat().st_size,'wall_seconds':time.monotonic()-start,'provenance_limitation':'Reconstructed using official custom-graph algorithm and the official default gene set. Original Dataverse default graph could not be downloaded (HTTP 403); byte identity to that archive is not established.'}
(P/'go_graph_reconstruction.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
