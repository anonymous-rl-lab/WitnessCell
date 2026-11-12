"""GSE220974 counts -> capture-corrected, mode-aware common data.

Gene selection and objective directions use only true controls and singles.
No outcome-dependent condition selection or examination of query winners.
The design whitelist is applied later, independently of these moments.
"""
from pathlib import Path
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
import hashlib,json,time
import h5py,numpy as np,pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parent
AUDIT=ROOT.parent/'nextpro2_assets_20261102/extracted/WitnessCell_GSE220974_Intake_Audit_20261102/audit'
OUT=ROOT/'data'; OUT.mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 t=time.monotonic();obs=pd.read_csv(AUDIT/'obs_audited.csv.gz',index_col=0)
 with h5py.File(ROOT/'raw/GSE220974.h5ad','r') as f:
  idx=f['obs/_index'].asstr()[:].astype(str);genes=f['var/_index'].asstr()[:].astype(str)
  assert np.array_equal(obs.index.to_numpy(),idx)
  X=sparse.csr_matrix((f['X/data'][:].astype(np.float64),f['X/indices'][:],f['X/indptr'][:]),shape=(len(idx),len(genes)))
 sums=np.asarray(X.sum(1)).ravel(); assert np.all(sums>0)
 X=sparse.diags(1e4/sums)@X;X.data=np.log1p(X.data)
 allowed=obs.audit_class.isin(['NTC','single_gene','distinct_gene_double']).to_numpy()
 single=(obs.audit_class=='single_gene').to_numpy()
 assert single.sum()==13226
 base=(obs.audit_class=='NTC').to_numpy(); assert base.sum()==3944
 permitted=single|base
 ref=X[permitted];n=permitted.sum()
 mu=np.asarray(ref.mean(0)).ravel();var=np.asarray(ref.power(2).mean(0)).ravel()-mu*mu
 targets=sorted(obs.loc[single,'audit_identity'].unique())
 base_genes=sorted({n.rsplit('-',1)[0] for n in targets});gi={g:i for i,g in enumerate(genes)}
 assert set(base_genes)<=set(gi)
 chosen=set(np.argsort(-var,kind='stable')[:5000]);chosen.update(gi[g] for g in base_genes)
 keep=np.array(sorted(chosen));genes=genes[keep];X=X[:,keep].toarray().astype(np.float32)
 del ref
 captures=obs['orig.ident'].astype(str).to_numpy()
 control={c:X[base&(captures==c)].mean(0,dtype=np.float64) for c in sorted(set(captures))}
 common_control=np.mean([control[c] for c in ('D1','D2','D3')],axis=0)
 # Cell and mean representations agree exactly up to float32 export rounding.
 corrected=X.astype(np.float64)
 for c,v in control.items():corrected[captures==c]+=common_control-v
 names=np.where(base,'control',obs.audit_identity.astype(str)).astype(str)
 conditions=['control']+targets+sorted(set(names[allowed])-{'control'}-set(targets))
 lookup={c:i for i,c in enumerate(conditions)};ci=np.array([lookup.get(c,-1) for c in names])
 halves=np.full(len(obs),-1,int)
 rng=np.random.default_rng(20261102)
 means=[];first=[];second=[];variances=[];counts=[];cap_effects=[]
 for c in conditions:
  rows=np.flatnonzero(allowed&(names==c));caplist=sorted(set(captures[rows]));full=[];h1=[];h2=[];v=[];ne=[]
  cs=np.full((3,len(genes)),np.nan)
  for cap in caplist:
   ids=rows[captures[rows]==cap];perm=rng.permutation(ids)
   halves[perm[::2]]=0;halves[perm[1::2]]=1
   a=corrected[perm[::2]];b=corrected[perm[1::2]]
   # Subtract separately estimated control halves below, after all labels exist.
   full.append(corrected[ids].mean(0));h1.append(a.sum(0));h2.append(b.sum(0))
   v.append(corrected[ids].var(0,ddof=1)/len(ids) if len(ids)>1 else np.zeros(len(genes)));ne.append(len(ids))
   if cap in ('D1','D2','D3'):cs[('D1','D2','D3').index(cap)]=full[-1]-common_control
  means.append(np.average(full,axis=0,weights=ne)-common_control)
  first.append(np.sum(h1,axis=0)/np.sum(halves[rows]==0)-common_control)
  second.append(np.sum(h2,axis=0)/np.sum(halves[rows]==1)-common_control)
  variances.append(np.sum(np.asarray(v)*np.asarray(ne)[:,None]**2,axis=0)/sum(ne)**2);counts.append(sum(ne));cap_effects.append(cs)
 # Shared full-control subtraction would induce artificial half agreement.
 for h,arr in enumerate((first,second)):
  shifts={cap:corrected[base&(captures==cap)&(halves==h)].mean(0)-common_control for cap in control}
  for j,c in enumerate(conditions):
   caps=sorted(set(captures[allowed&(names==c)]))
   nh=np.array([np.sum(allowed&(names==c)&(captures==cap)&(halves==h)) for cap in caps])
   arr[j]=arr[j]-np.average([shifts[cap] for cap in caps],axis=0,weights=nh)
 means=np.asarray(means); first=np.asarray(first);second=np.asarray(second)
 # Controls are defined as zero effects; raw control cells retain heterogeneity.
 means[0]=0;first[0]=0;second[0]=0
 gene_lut={g:i for i,g in enumerate(genes)};target_indices=np.array([gene_lut[t.rsplit('-',1)[0]] for t in targets])
 mask=~np.isin(genes,base_genes)
 E=means[1:1+len(targets)];a=first[1:1+len(targets)];b=second[1:1+len(targets)]
 agreement=np.divide(2*a*b,a*a+b*b,out=np.zeros_like(a),where=a*a+b*b>0).clip(0,1)
 V=E*agreement;V[:,~mask]=0;norm=np.abs(V).sum(1)
 if np.any(norm==0):raise ValueError('An objective has no reproducible single-effect direction')
 V/=norm[:,None]
 np.savez_compressed(OUT/'moments.npz',genes=genes,conditions=np.array(conditions),means=means,first_means=first,second_means=second,
  mean_variances=np.array(variances),counts=np.array(counts),capture_effects=np.array(cap_effects),capture_names=np.array(['D1','D2','D3']),
  target_names=np.array(targets),target_indices=target_indices,base_genes=np.array(base_genes),single_effects=E,control=common_control,
  score_gene_mask=mask,objectives=V,objective_agreement=agreement,source_gene_indices=keep)
 with h5py.File(OUT/'cells.h5','w') as f:
  ids=np.flatnonzero(allowed)
  f.create_dataset('X',data=corrected[ids].astype(np.float32),compression='gzip',compression_opts=1,chunks=(64,len(genes)))
  for name,value in [('genes',genes),('condition_names',np.array(conditions)),('cell_ids',idx[ids]),('capture',captures[ids])]:f.create_dataset(name,data=value.astype(object),dtype=h5py.string_dtype())
  f.create_dataset('condition_index',data=ci[ids]);f.create_dataset('half',data=halves[ids]);f.create_dataset('score_gene_mask',data=mask)
  # Equal condition mass; all observed cells within a condition equally weighted.
  weights=np.zeros(len(ids))
  for c in conditions:
   sel=names[ids]==c;weights[sel]=1/sel.sum()
  f.create_dataset('sampling_weight',data=weights)
  f.attrs['expression']='log1p(10000*counts/library_size), capture NTC centered and rebased to equal-D1/D2/D3 control; may be negative'
 obs.loc[allowed,['audit_identity','audit_class','orig.ident','guide_group']].to_csv(OUT/'cell_metadata.csv.gz')
 record=dict(status='PREPARED_NOT_SCORED',source_sha256='ca6f2f11401b691cafea16b3682a3ef24aed61ba769a4d5796c2672ece9499be',
  cells=int(allowed.sum()),genes=len(genes),score_genes=int(mask.sum()),single_identities=len(targets),base_genes=len(base_genes),
  conditions=len(conditions),normalization='log1p library-size 10000; cell-weighted condition effects with capture-matched NTC; rare one-cell single/capture strata kept once only',
  gene_selection='5000 largest variances in true NTC and single cells plus all target genes; no double cells',
  split='condition-by-capture hash-independent RNG seed20261102; controls separately split; no shared control noise between halves',
  objective='existing TargetRank L1-normalized single effect times half agreement, all physical target genes excluded',
  double_design_status='not yet whitelisted; no double outcomes used for model or condition selection',
  seconds=time.monotonic()-t)
 (OUT/'PREPROCESS.json').write_text(json.dumps(record,indent=2));print(json.dumps(record),flush=True)
if __name__=='__main__':main()
