"""Independent table-level verification for v41; no fitting or network."""
from pathlib import Path
import re,json,argparse,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--paper',type=Path,default=R.parent/'paper');args=ap.parse_args();P=args.paper
N=R/'inputs/native_results';out=R/'verification';out.mkdir(exist_ok=True)
checks=[]
def check(name,ok,detail=None):
 checks.append(dict(check=name,passed=bool(ok),detail=detail))
 if not ok:print('FAILED',name,detail)
fold=pd.read_csv(N/'FOLD_REPEAT_SCORES.csv');pool=pd.read_csv(N/'POOLED.csv')
cols=['full_sse','partner_sse','baseline_full_sse','baseline_partner_sse'];keys=['panel','baseline','method']
a=fold.groupby(keys)[cols].sum().sort_index();b=pool.set_index(keys)[cols].sort_index()
check('48 pooled rows reconstructed from 1,080 fold-repeat rows',a.index.equals(b.index) and len(a)==48 and len(fold)==1080)
check('Pooled SSE matches fold sums',np.allclose(a,b,rtol=1e-12,atol=1e-10),float(np.max(np.abs(a.to_numpy()-b.to_numpy()))))
for metric in ['full','partner']:
 v=100*(1-a[metric+'_sse']/a['baseline_'+metric+'_sse'])
 z=pool.set_index(keys)[metric+'_gain_vs_own_baseline_pct'].sort_index()
 check(metric+' pooled gain recalculation',np.allclose(v,z,rtol=1e-12,atol=1e-11),float(np.max(abs(v-z))))
panels=['full_raw','full_published'];bases=['native','saturation'];fams=['rank1','spectral']
sel=pd.read_csv(N/'SELECTIONS.csv');part=pd.read_csv(N/'BY_PARTITION.csv');rep=pd.read_csv(N/'BY_REPEAT.csv');ranks=pd.read_csv(N/'SPECTRAL_RANKS.csv')
t45=pd.read_csv(R/'tables/Table_S45.csv');t46=pd.read_csv(R/'tables/Table_S46.csv');t47=pd.read_csv(R/'tables/Table_S47.csv')
t48=pd.read_csv(R/'tables/Table_S48.csv');t49=pd.read_csv(R/'tables/Table_S49.csv');t50=pd.read_csv(R/'tables/Table_S50.csv')
def fmt(v,n=4):return f'{0 if abs(float(v))<5e-9 else float(v):.{n}f}'
for i,(pn,base,fam) in enumerate((p,b,f) for p in panels for b in bases for f in fams):
 original=a.loc[(pn,base,fam+'_original')];controlled=a.loc[(pn,base,fam+'_controlled')]
 vals=[]
 for x in [original,controlled]:
  vals.extend([100*(1-x.full_sse/x.baseline_full_sse),100*(1-x.partner_sse/x.baseline_partner_sse)])
 check(f'S45 row {i+1}',np.allclose(t45.iloc[i,3:].to_numpy(float),np.array(vals).round(4),atol=1e-9))
 z=sel[(sel.panel==pn)&(sel.baseline==base)&(sel.family==fam)]
 check(f'S47 row {i+1} zero/strength',t47.iloc[i,3]==f'{z.selected.eq("G0").sum()}/15' and t47.iloc[i,4]==f'{z.alpha.min():.6f}–{z.alpha.max():.6f}' and t47.iloc[i,5]==f'{z.gamma.min():.6f}–{z.gamma.max():.6f}')
 for pol,j in [('original',3),('controlled',5)]:
  v=part[(part.panel==pn)&(part.baseline==base)&(part.method==fam+'_'+pol)].partner_gain_vs_own_baseline_pct
  check(f'S46 row {i+1} {pol}',t46.iloc[i,j]==f'{v.min():.4f} to {v.max():.4f}' and t46.iloc[i,j+1]==f'{(v>1e-9).sum()}/5')
 for j,ref in enumerate(['original','off']):
  rr=a.loc[(pn,base,fam+'_'+ref)]
  vals=[100*(1-controlled[m+'_sse']/rr[m+'_sse']) for m in ['full','partner']]
  check(f'S48 row {i+1} vs {ref}',np.allclose(t48.iloc[i,3+j*2:5+j*2].to_numpy(float),np.round(vals,4),atol=1e-9))
 same=a.loc[(pn,base,fam+'_whole_same')];whole=a.loc[(pn,base,fam+'_whole_full')]
 vals=[100*(1-controlled.full_sse/same.full_sse),100*(1-controlled.full_sse/whole.full_sse),100*(1-controlled.partner_sse/whole.partner_sse)]
 check(f'S49 row {i+1}',np.allclose(t49.iloc[i,3:].to_numpy(float),np.round(vals,4),atol=1e-9))
 # Exact structural identity is checked on original sums before converting to percentages.
 check(f'Preservation identity {pn}/{base}/{fam}',abs((original.full_sse-controlled.full_sse)-(original.partner_sse-controlled.partner_sse))<1e-8)
 check(f'Same-alpha partner equality {pn}/{base}/{fam}',abs(same.partner_sse-controlled.partner_sse)<1e-8)
for i,(pn,base) in enumerate((p,b) for p in panels for b in bases):
 x=a.loc[(pn,base,'baseline')];y=a.loc[(pn,base,'additive')]
 vals=[x.full_sse,x.partner_sse,100*(1-x.full_sse/y.full_sse),100*(1-x.partner_sse/y.partner_sse)]
 check(f'S50 row {i+1}',np.allclose(t50.iloc[i,2:].to_numpy(float),[round(vals[0],6),round(vals[1],6),round(vals[2],4),round(vals[3],4)],rtol=0,atol=1e-9))
main=(P/'WitnessCell_v41_Main_Methods.md').read_text();si=(P/'WitnessCell_v41_Supplementary_Information.md').read_text()
check('57 sequential SI tables',[int(x) for x in re.findall(r'^### Supplementary Table S(\d+) \|',si,re.M)]==list(range(1,58)))
check('24 sequential SI notes',[int(x) for x in re.findall(r'^## Supplementary Note (\d+) \|',si,re.M)]==list(range(1,25)))
check('No main-text source-function transfer results',all(x not in main for x in ['SF12','V4','0.2617','0.6513','12.839']))
check('Historical material results retained',all(x in si for x in ['−0.9961%','+0.0350%','0.175 / 0.193','−0.0081% / −0.0222%']))
check('Native primary and stronger-reference results present',all(x in main for x in ['9.09%','1.91%','4.82%','7.04%','2.44%','2.77%']))
check('Scientific core retained',all(x in main for x in ['92.84%','0.0625','3.515','6.43%','6.18%','15.96%','6.94%']))
for fn in ['WitnessCell_v41_Main_Methods.md','WitnessCell_v41_Supplementary_Information.md','Figure_Legends.md','Supplementary_Figure_Legends.md']:
 text=(P/fn).read_text()
 check('No control characters '+fn,not any(ord(c)<32 and c not in '\n\t' for c in text))
 check('Even display-math delimiters '+fn,text.count('$$')%2==0)
 check('No placeholder/editorial TODO '+fn,not re.search(r'\b(?:TODO|FIXME|INSERT HERE|TBD)\b',text))
leg=(P/'Supplementary_Figure_Legends.md').read_text()
check('Seven supplement figure legends',[int(x) for x in re.findall(r'^## Supplementary Figure S(\d+)',leg,re.M)]==list(range(1,8)))
check('Supplement legends synchronized (whitespace-normalized)',leg.split('# Supplementary figure legends',1)[1].split()==si.split('## Supplementary figure legends',1)[1].split('## Supplementary reference correspondence',1)[0].split())
# Paragraphs/headings, abstract and refs are counted without skipping the first introduction paragraph.
abstract,body=main.split('## Abstract\n',1)[1].strip().split('\n\n',1)
body=body.split('## Methods\n',1)[0];body='\n'.join(x for x in body.splitlines() if not x.startswith('#'))
methods=main.split('## Methods\n',1)[1].split('## Data availability',1)[0]
words={'version':'v41','rule':'Whitespace-delimited; body includes Introduction, Results, Discussion; excludes headings, abstract, Methods, availability, references and legends. Math tokens are counted under the same rule.','abstract':len(abstract.split()),'main_body':len(body.split()),'methods_excluding_headings':len('\n'.join(x for x in methods.splitlines() if not x.startswith('#')).split()),'references':len(re.findall(r'^\d+\.',main.split('## References\n',1)[1],re.M)),'main_figures':5,'supplementary_figures':7,'supplementary_notes':24,'supplementary_tables':57}
check('Abstract <=150 and body <=5000',words['abstract']<=150 and words['main_body']<=5000)
check('46 reference entries',words['references']==46)
(P/'WORD_COUNTS.json').write_text(json.dumps(words,indent=2)+'\n')
result={'scope':'Table/report verification from existing saved scores, selections and sources; no model retraining or raw sequencing reconstruction.','checks':checks,'passed':sum(x['passed'] for x in checks),'total':len(checks),'all_passed':all(x['passed'] for x in checks),'word_counts':words}
(out/'V41_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['passed','total','all_passed','word_counts']},ensure_ascii=False))
assert result['all_passed']
