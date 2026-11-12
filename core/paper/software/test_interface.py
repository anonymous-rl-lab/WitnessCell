#!/usr/bin/env python3
"""Meaningful interface checks: native parity, leakage, shape, input dependence."""
import copy,json,hashlib
from pathlib import Path
import numpy as np
from witnesscell_predict import predict,validate,fit_predict
p=Path(__file__).resolve().parent;d=json.loads((p/'example_input.json').read_text());checks={}
f,meta=predict(d);o=d['observations']
direct=fit_predict(means={k:np.array(v['mean']) for k,v in o.items()},variances={k:np.array(v['variance']) for k,v in o.items()},counts={k:v['n'] for k,v in o.items()},genes=d['genes'],train_conditions=d['train_conditions'],validation_conditions=d['validation_conditions'],target_conditions=d['query_conditions'],identity_mode='multihead_v14',gene2go=d['gene2go'])
checks['native_prediction_exact_parity']=all(np.array_equal(f.predictions[k],direct.predictions[k]) for k in f.predictions)
checks['finite_output']=all(np.isfinite(v).all() for v in f.predictions.values())
checks['query_truth_absent']=not set(d['query_conditions'])&set(o)
checks['unseen_single_supported']='G7' not in f.known_single_genes and any('G7' in q for q in d['query_conditions'])
for kind in ['query_leak','reversed_pair_overlap','negative_variance','wrong_gene_axis']:
 bad=copy.deepcopy(d)
 if kind=='query_leak':bad['observations'][bad['query_conditions'][0]]=copy.deepcopy(o['control'])
 elif kind=='reversed_pair_overlap':bad['query_conditions'][0]='+'.join(reversed(next(c for c in d['train_conditions'] if '+' in c).split('+')))
 elif kind=='negative_variance':bad['observations']['control']['variance'][0]=-1
 else:bad['observations']['control']['mean'].pop()
 try:validate(bad);checks[kind+'_rejected']=False
 except ValueError:checks[kind+'_rejected']=True
changed=copy.deepcopy(d)
for row in changed['observations'].values():row['mean']=[v+.3 for v in row['mean']]
f2,_=predict(changed)
checks['common_shift_equivariance']=all(np.allclose(f2.predictions[k],f.predictions[k]+.3,rtol=1e-10,atol=1e-10) for k in f.predictions)
report={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (p/'vendor').glob('*.py')}}
(p/'INTERFACE_TESTS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));assert all(checks.values())
