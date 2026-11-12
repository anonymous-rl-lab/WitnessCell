#!/usr/bin/env python3
"""Descriptive contrasts from existing Jurkat predictions; no fitting or outcome selection."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
R=Path(__file__).resolve().parents[1]
p=json.loads((R/'revision/FOUR_CORNER_SELECTION_V28.json').read_text())
source=R/p['input']
assert hashlib.sha256(source.read_bytes()).hexdigest()==p['source_sha256']
with source.open() as f: data=list(csv.DictReader(f))
lookup={(frozenset([r['gene_a'],r['gene_b']]),int(r['repeat'])):r for r in data}
assert len(lookup)==1234
genes=sorted({r[k] for r in data for k in ['gene_a','gene_b']})
A,B,D=[p['fixed_'+k] for k in ['A','B','D']]
eligible=[C for C in genes if C not in [A,B,D] and all((frozenset(pair),rep) in lookup for pair in [(A,B),(A,D),(C,B),(C,D)] for rep in [1,2])]
assert eligible==p['eligible_C'] and eligible[0]==p['selected_C']
raw=[];contrasts=[];components=['observed','M0','rank1_selected']
for C in eligible:
 for rep in [1,2]:
  vals={}
  for label,pair,sign in [('AB',(A,B),1),('AD',(A,D),-1),('CB',(C,B),-1),('CD',(C,D),1)]:
   r=lookup[frozenset(pair),rep]
   assert r['stage']=='Jurkat' and r['budget']=='high' and int(r['fold'])==0
   out=dict(background=C,selected=C==p['selected_C'],repeat=rep,corner=label,gene_left=pair[0],gene_right=pair[1],sign=sign,pair_id=r['pair_id'])
   for k in components:out[k]=float(r[k])
   raw.append(out);vals[label]=out
  out=dict(background=C,selected=C==p['selected_C'],repeat=str(rep))
  for k in components:out[k]=sum(vals[c]['sign']*vals[c][k] for c in vals)
  out['baseline_squared_error']=(out['M0']-out['observed'])**2
  out['rank1_squared_error']=(out['rank1_selected']-out['observed'])**2
  contrasts.append(out)
 rr=[r for r in contrasts if r['background']==C]
 out=dict(background=C,selected=C==p['selected_C'],repeat='mean')
 for k in components:out[k]=float(np.mean([r[k] for r in rr]))
 out['baseline_squared_error']=(out['M0']-out['observed'])**2
 out['rank1_squared_error']=(out['rank1_selected']-out['observed'])**2
 contrasts.append(out)
for name,rows in [('JURKAT_FOUR_CORNERS_V28.csv',raw),('JURKAT_FOUR_CORNER_CONTRASTS_V28.csv',contrasts)]:
 with (R/'revision'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps([r for r in contrasts if r['repeat']=='mean'],indent=2))
