import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--compatible',action='store_true');args=ap.parse_args()
R=args.root.resolve()
sys.path.insert(0,str(R/'code'));sys.path.insert(0,str(R/'verification'))
from common import loadnpz,innerfold,savez,digest
from compact_replay import matrix
import baselines
O=R;rows=[];expected={x['path']:x for x in json.loads(args.manifest.read_text())}
for t in json.loads((R/'TASKS.json').read_text()):
 d=loadnpz(R/'inputs'/f"{t['panel']}_inputs.npz");ids=json.loads((R/'inputs/splits'/t['panel']/str(t['seed'])/str(t['fold'])/'identities.json').read_text());tr,qu=np.array(ids['train']),np.array(ids['query']);lut={g:i for i,g in enumerate(d['genes'])};ed=np.array([[lut[a],lut[b]] for a,b in d['pairs']]);n=len(d['genes']);p=Path('fitted')/t['panel']/t['baseline']/str(t['seed'])/str(t['fold']);src=R/p;sel=json.loads((src/'SELECTION.json').read_text());nr=d['singles'].shape[1]
 b=baselines.predict(json.loads((src/'baseline.json').read_text()),d['singles'],ed[qu]);folds=innerfold(d['pairs'][tr]);ib=np.zeros((len(tr),nr));q={'query_ids':qu,'B':b};inn={'B':ib,'folds':folds}
 for k in range(3):
  va=folds==k;ib[va]=baselines.predict(json.loads((src/f'inner{k}_baseline.json').read_text()),d['singles'],ed[tr[va]])
 for fam in ('rank1','spectral'):
  G=matrix(src/f'{fam}_final.npz');dd=np.repeat(G[ed[qu,0],ed[qu,1]][:,None],nr,axis=1)
  if sel[fam]['selected']=='G0':dd=np.zeros_like(b)
  q[fam+'_D']=dd;di=np.zeros_like(ib)
  for k in range(3):
   va=folds==k;G=matrix(src/f'{fam}_inner{k}.npz');z=np.repeat(G[ed[tr[va],0],ed[tr[va],1]][:,None],nr,axis=1)
   if sel[fam]['selected']=='G0':z*=0
   di[va]=(ib[va]+z)-ib[va]
  inn[fam+'_D']=di
 for name,vals in [('QUERY.npz',q),('INNER.npz',inn)]:
  target=O/p/name;target.parent.mkdir(parents=True,exist_ok=True);savez(target,**vals);ref=expected[str(p/name)]
  same=digest(target)==ref['sha256']
  if not args.compatible and not same:raise RuntimeError('Cache byte mismatch: '+str(p/name)+'; use the recorded environment or --compatible followed by numerical weight replay.')
  rows.append({'path':str(p/name),'byte_identical':same,'reference_sha256':ref['sha256'],'reconstructed_sha256':digest(target)})
report={'status':'PASS','files':len(rows),'byte_identical':sum(x['byte_identical'] for x in rows),'query_outcomes_read':False,'predictor_weights_refitted':False,'precision':'compatible; numerical acceptance requires weight replay' if args.compatible else 'exact','checks':rows}
args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}))
