"""Recreate current Tables S44-S53 from frozen full-precision source CSVs."""
from pathlib import Path
import argparse
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=R/'table_replay');a=ap.parse_args();O=a.out;O.mkdir(parents=True,exist_ok=True)
N=R/'inputs/native_results';pool=pd.read_csv(N/'POOLED.csv');part=pd.read_csv(N/'BY_PARTITION.csv');rep=pd.read_csv(N/'BY_REPEAT.csv');sel=pd.read_csv(N/'SELECTIONS.csv');rank=pd.read_csv(N/'SPECTRAL_RANKS.csv')
headers={i:list(pd.read_csv(R/f'tables/Table_S{i}.csv').columns) for i in range(44,54)}
pn={'full_raw':'Reconstructed','full_published':'Published'};bn={'native':'Native','saturation':'Saturation'};fn={'rank1':'Rank one','spectral':'Spectral'}
configs=[(p,b,f) for p in pn for b in bn for f in fn]
def fmt(v,n=4):return f'{0 if abs(float(v))<5e-9 else float(v):.{n}f}'
def row(p,b,m):return pool[(pool.panel==p)&(pool.baseline==b)&(pool.method==m)].iloc[0]
def labels(p,b,f):return [pn[p],bn[b],fn[f]]
def write(i,rows):
 d=pd.DataFrame(rows,columns=headers[i]);d.to_csv(O/f'Table_S{i}.csv',index=False)
 text='| '+' | '.join(headers[i])+' |\n|'+'|'.join('---' for _ in headers[i])+'|\n'
 text+=''.join('| '+' | '.join(str(x) for x in r)+' |\n' for r in rows)
 (O/f'Table_S{i}.md').write_text(text)
 ref=pd.read_csv(R/f'tables/Table_S{i}.csv',dtype=str).fillna('')
 got=pd.read_csv(O/f'Table_S{i}.csv',dtype=str).fillna('')
 assert got.equals(ref),f'Table S{i} differs'
for i in [44,51]:write(i,pd.read_csv(R/f'inputs/Table_S{i}_definition.csv',dtype=str).fillna('').values.tolist())
for num in [45,46,47,48,49]:
 rows=[]
 for p,b,f in configs:
  r=labels(p,b,f);x=row(p,b,f+'_original');c=row(p,b,f+'_controlled')
  if num==45:
   r+=[fmt(x.full_gain_vs_own_baseline_pct),fmt(x.partner_gain_vs_own_baseline_pct),fmt(c.full_gain_vs_own_baseline_pct),fmt(c.partner_gain_vs_own_baseline_pct)]
  elif num==46:
   for pol in ['original','controlled']:
    v=part[(part.panel==p)&(part.baseline==b)&(part.method==f+'_'+pol)].partner_gain_vs_own_baseline_pct
    r += [f'{v.min():.4f} to {v.max():.4f}',f'{(v>1e-9).sum()}/5']
   z=rep[(rep.panel==p)&(rep.baseline==b)&(rep.method==f+'_original')].sort_values('repeat')
   r+=[' / '.join(fmt(v) for v in z.partner_gain_vs_own_baseline_pct) if p=='full_raw' else 'One aggregate']
  elif num==47:
   z=sel[(sel.panel==p)&(sel.baseline==b)&(sel.family==f)];rr=rank[(rank.panel==p)&(rank.baseline==b)]['rank']
   r += [f'{z.selected.eq("G0").sum()}/15',f'{z.alpha.min():.6f}–{z.alpha.max():.6f}',f'{z.gamma.min():.6f}–{z.gamma.max():.6f}',f'{rr.min()}–{rr.max()}' if f=='spectral' else 'At most 1']
  elif num==48:
   for pol in ['original','off']:
    z=row(p,b,f+'_'+pol)
    r += [fmt(100*(1-c[m+'_sse']/z[m+'_sse'])) for m in ['full','partner']]
  else:
   z=row(p,b,f+'_whole_same');w=row(p,b,f+'_whole_full')
   r += [fmt(100*(1-c.full_sse/z.full_sse)),fmt(100*(1-c.full_sse/w.full_sse)),fmt(100*(1-c.partner_sse/w.partner_sse))]
  rows.append(r)
 write(num,rows)
rows=[]
for p in pn:
 for b in bn:
  x=row(p,b,'baseline');z=row(p,b,'additive')
  rows.append([pn[p],bn[b],fmt(x.full_sse,6),fmt(x.partner_sse,6),fmt(100*(1-x.full_sse/z.full_sse)),fmt(100*(1-x.partner_sse/z.partner_sse))])
write(50,rows)
q=pd.read_csv(R/'inputs/diagnostics/FUNCTIONAL_REPEATABILITY.csv')
rows=[[x.feature_id.replace('_',' '),str(x.training_pairs),str(x.query_pairs)]+[f'{v:.4f}' for v in [x.training_repeat_pearson,x.query_repeat_pearson,x.original_partner_gain_pct,x.controlled_partner_gain_pct]] for x in q.itertuples()]
write(52,rows)
q=pd.read_csv(R/'inputs/diagnostics/FUNCTIONAL_ORACLES.csv').set_index('feature_id')
rows=[]
for name in ['HeLa','A549','293T','Jurkat']:
 x=q.loc[name+'_high'];rows.append([name]+[f'{x[k]:.4f}' for k in ['actual_alpha','unconstrained_alpha','bounded_alpha','actual_gain_pct','bounded_oracle_gain_pct','nonnegative_oracle_gain_pct']])
write(53,rows)
print('Tables S44-S53 reconstructed; all displayed values match.')
