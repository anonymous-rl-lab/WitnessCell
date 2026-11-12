#!/usr/bin/env python3
"""Generate synthetic input for API demonstration, not biological validation."""
from pathlib import Path
import itertools,json,numpy as np
rng=np.random.default_rng(25); genes=['G'+str(i) for i in range(24)];endpoints=genes[:8]
control=np.full(24,2.);effects=rng.normal(0,.15,(8,24))
for i in range(8):effects[i,i]-=.5
pairs=list(itertools.combinations(range(8),2));train=pairs[:18];val=pairs[18:23];query=pairs[23:]
def label(p):return '+'.join(genes[i] for i in p)
def response(p):
 e=effects[list(p)].sum(0);return control+e/(1+.1*np.abs(e))
obs={'control':{'mean':control.tolist(),'variance':[.1]*24,'n':100}}
for i in range(7):obs[genes[i]]={'mean':response((i,)).tolist(),'variance':[.1]*24,'n':80}
for p in train+val:obs[label(p)]={'mean':response(p).tolist(),'variance':[.1]*24,'n':80}
d={'description':'Synthetic interface example. No query outcomes. G7 has no measured training single.',
'genes':genes,'observations':obs,'train_conditions':genes[:7]+[label(p) for p in train],
'validation_conditions':[label(p) for p in val],'query_conditions':[label(p) for p in query],
'gene2go':{g:['SYNTHETIC_GROUP_'+str(i%3)] for i,g in enumerate(genes)}}
out=Path(__file__).resolve().parent/'example_input.json';out.write_text(json.dumps(d,indent=2)+'\n');print(out)
