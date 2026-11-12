"""Fixed-model sensitivity; no outcome cutoff or hyperparameter search."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np
from common import weights, split, fit, metric, configs
R=Path(__file__).resolve().parents[1];O=R/'results'
d=pd.read_csv(O/'ENRICHED_INPUT.csv');d['zero_to_reference']=d.zero_to_reference.fillna(False).astype(bool)
subsets={
 'all':d,
 'no_terminal_checkpoints':d[~d.model.str.contains('_terminal_')],
 'no_Joung':d[d.panel!='Joung'],
 'RNA_only_no_Joung':d[(d.modality!='growth')&(d.panel!='Joung')],
 'growth_selected_rank1':d[(d.modality=='growth')&d.model.isin(['rank1','rank1_selected'])]}
rows=[];held=[]
for label,z in subsets.items():
 z=z.reset_index(drop=True)
 for name,(cols,cats) in configs.items():
  for scheme,group in [('fit',None),('hold_query','cohort'),('hold_panel','panel'),('hold_study','study')]:
   if group and z[group].nunique()<2:continue
   p=np.empty(len(z))
   if group:
    for val in z[group].unique():
     tr,te=split(z,group,val)
     if not tr.any():p[te]=np.nan
     else:p[te]=fit(z[tr],z[te],cols,cats)[0]
   else:p=fit(z,z,cols,cats)[0]
   p[z.zero_to_reference]=0
   for scope,ix in [('all',np.ones(len(z),bool)),('active',~z.zero_to_reference.to_numpy())]:
    ix&=np.isfinite(p)
    rows.append(dict(subset=label,regression=name,scheme=scheme,scope=scope,n=int(ix.sum()),panels=z[ix].panel.nunique(),studies=z[ix].study.nunique(),**metric(z[ix],p[ix])))
   if scheme=='hold_study':
    for panel in z.panel.unique():
     ix=(z.panel==panel).to_numpy()&np.isfinite(p)
     held.append(dict(subset=label,regression=name,panel=panel,**metric(z[ix],p[ix])))
out=pd.DataFrame(rows);out.to_csv(O/'REGRESSION_SENSITIVITY.csv',index=False)
pd.DataFrame(held).to_csv(O/'REGRESSION_STUDY_HELD_BY_PANEL.csv',index=False)
print(out.query("scope=='all' and scheme=='hold_study'")[['subset','regression','n','rmse_pp','R2']].to_string(index=False))
