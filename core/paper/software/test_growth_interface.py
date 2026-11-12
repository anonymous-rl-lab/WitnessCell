#!/usr/bin/env python3
from pathlib import Path
import copy,json,hashlib
import numpy as np
from witnesscell_growth import validate
R=Path(__file__).resolve().parent
x=json.loads((R/'growth_example_input.json').read_text());checks=[]
def check(n,p,value=None):
 checks.append(dict(check=n,passed=bool(p),value=value));assert p,n
train,q,ids=validate(x)
check('Real input has training outcomes but identity-only queries',len(train['edges'])==2464 and len(q)==617)
with np.load(R/'growth_example_expected.npz',allow_pickle=False) as a,np.load(R/'growth_example_output/results/high/predictions.npz',allow_pickle=False) as b:
 check('Query identity order matches native archive',np.array_equal(a['global_ids'],b['global_ids']))
 for k in ['M0','rank1_selected','spectral_selected']:
  delta=float(np.max(np.abs(a[k]-b[k])));check('Native functional prediction parity '+k,delta==0,delta)
for name,mutate in [('Reverse-pair train/query leakage',lambda y:y['queries'].__setitem__(0,list(reversed(y['training'][0]['pair'])))),('Query outcomes',lambda y:y['queries'].__setitem__(0,{'pair':y['queries'][0],'response':[0,0]})),('Nonfinite training response',lambda y:y['training'][0].__setitem__('response',[float('nan'),0]))]:
 y=copy.deepcopy(x);mutate(y)
 try:validate(y);rejected=False
 except (ValueError,TypeError):rejected=True
 check(name+' rejected',rejected)
p=json.loads((R/'GROWTH_INPUT_PROVENANCE.json').read_text())
check('Three native functional modules unmodified',all(hashlib.sha256((R/'growth_vendor'/n).read_bytes()).hexdigest()==h for n,h in p['vendor'].items()))
result=dict(passed=all(c['passed'] for c in checks),checks=checks,scope='Interface and exact archived-prediction parity; not new biological confirmation')
(R/'GROWTH_INTERFACE_TESTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
