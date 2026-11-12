#!/usr/bin/env python3
"""All-task application summaries; uses every frozen task, no outcome filtering."""
from pathlib import Path
import csv,json
import numpy as np
R=Path(__file__).resolve().parents[1];D=R/'revision';G=R/'evidence/source_data/gse220974'
def read(p):
 with p.open(newline='') as h:return list(csv.DictReader(h))
def write(name,rows):
 with (D/name).open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
rows=read(G/'CHOICES.csv');groups=json.loads((G/'GROUPS.json').read_text());idx={(r['group'],r['model']):r for r in rows}
assert len(rows)==420 and len(groups)==28
models=['R_G','StandardRidge_RNA','Scouter_Ensemble','GEARS_Ensemble','SingleB_additive','PRESAGE_Combo_seed24']
ordered=sorted(groups,key=lambda r:(r['anchor'].rsplit('-',1)[0],r['anchor'],r['fold']))
out=[];summary=[]
for g in ordered:
 row={'group':g['group'],'physical_anchor':g['anchor'].rsplit('-',1)[0],'anchor':g['anchor'],'fold':g['fold'],'candidate_count':g['n']}
 for m in models:
  r=idx[(g['group'],m)];row[m+'_regret']=float(r['regret']);row[m+'_utility']=float(r['actual_U']);row[m+'_choice']=r['pair'];row[m+'_hit']=int(r['hit'])
 out.append(row)
for m in models:
 rr=[idx[(g['group'],m)] for g in ordered];diff=np.array([float(idx[(g['group'],'R_G')]['actual_U'])-float(r['actual_U']) for g,r in zip(ordered,rr)])
 summary.append({'reference':m,'tasks':28,'physical_anchors':8,'mean_regret':np.mean([float(r['regret']) for r in rr]),'median_regret':np.median([float(r['regret']) for r in rr]),'hits':sum(int(r['hit']) for r in rr),'RG_wins':int((diff>1e-12).sum()),'RG_ties':int((abs(diff)<=1e-12).sum()),'RG_losses':int((diff<-1e-12).sum()),'RG_same_choice':sum(idx[(g['group'],'R_G')]['pair']==r['pair'] for g,r in zip(ordered,rr))})
write('PAIRPICK_ALL_28_TASKS.csv',out);write('PAIRPICK_ALL_TASK_SUMMARY.csv',summary)
assert len({r['physical_anchor'] for r in out})==8
(D/'APPLICATION_ANALYSIS.json').write_text(json.dumps({'status':'PASS','all_tasks_included':28,'source_choice_records':420,'physical_anchors':8,'ordering':'physical anchor, mode-specific anchor, fold; outcomes not used','tie_tolerance':1e-12},indent=2)+'\n')
print(json.dumps(summary,indent=2))
