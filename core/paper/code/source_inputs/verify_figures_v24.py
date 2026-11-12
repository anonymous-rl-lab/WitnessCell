#!/usr/bin/env python3
"""Check frozen numerical inputs and exported figure integrity; no model fitting."""
from pathlib import Path
import csv,hashlib,json,re
import numpy as np
from PIL import Image
from pypdf import PdfReader

R=Path(__file__).resolve().parents[1]
P=R/'evidence/source_data/publication'
F=R/'evidence/source_data/formal_bridge'
G=R/'evidence/source_data/gse220974'
checks=[]
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def check(name,ok,value=None):
    checks.append({'check':name,'passed':bool(ok),'value':value})
    if not ok:raise AssertionError(name)
def near(a,b):return np.isclose(a,b,rtol=1e-10,atol=1e-12)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

manifest=json.loads((R/'evidence/SOURCE_MANIFEST.json').read_text())['files']
check('All archived source bytes match SI manifest',all(sha(R/'evidence'/r['path'])==r['sha256'] for r in manifest),len(manifest))
check('SI text preserved',sha(R/'WitnessCell_v24_Supplementary_Information.md')=='26b979b09c1e35f0186e41ce662c57769d88aa292d2c6b3b28ab4c21456e33a8')
mi=json.loads((R/'MANUSCRIPT_FIGURE_INTEGRATION.json').read_text());t=(R/'WitnessCell_v24_Main_Methods.md').read_text()
for e in reversed(mi['edits']):
    assert t.count(e['after'])==1
    t=t.replace(e['after'],e['before'])
check('Main text differs only by recorded figure callouts',hashlib.sha256(t.encode()).hexdigest()==mi['original_main_sha256'],len(mi['edits']))

ss={r['method']:r for r in rows(F/'scores/summary.csv')}
gain=100*(1-float(ss['OLS_high']['four_world_rectangle_MSE'])/float(ss['OLS_low']['four_world_rectangle_MSE']))
check('Figure 2 information gain rounds to 92.84%',round(gain,2)==92.84,float(gain))
check('Figure 2 comparison-panel lower bound',float(ss['OLS_low']['low_information_floor'])==1/16)
lb=rows(P/'released_benchmark/formal_top100_leaderboard.csv')
check('Figure 3 complete 16-entry leaderboard',len(lb)==16 and lb[0]['method']=='WitnessCell' and float(lb[0]['official_position'])==1)
check('Figure 3 headline rank and population',round(float(lb[0]['mean_operation_rank']),3)==3.803 and int(lb[0]['operation_seed_units'])==654)
metric=[int(r['position']) for r in rows(P/'released_benchmark/official_top100_per_metric_leaderboard.csv') if r['method']=='WitnessCell']
check('Figure 3 metric profile preserved',sorted(metric)==[1,1,1,1,8,9],metric)
matched={r['model']:r for r in rows(P/'Step2_Fig5_Matched55.csv')}; seeds=rows(P/'Step2_Fig5_SeedMetrics.csv')
for model in ['GEARS','scouter']:
    sub=[r for r in seeds if r['model']==model and r['variant']=='raw']
    check('Figure 4 '+model+' uses exactly three raw-scale seeds',len(sub)==3)
    for key in ['response_MSE','C_MSE']:
        mean=np.mean([float(r[key]) for r in sub])
        check('Figure 4 '+model+' '+key+' seed-loss aggregation',near(mean,float(matched[model][key])),float(mean))
a=np.load(P/'RNA_signed_rectangle_predictions.npz');truth=a['truth_effect'];e0=(a['M0']-truth)**2;er=(a['R']-truth)**2
check('Figure 4 array dimensions',truth.shape==(11,4995),list(truth.shape))
check('Figure 4 M0 rectangle loss from saved predictions',near(e0.mean(),float(matched['M0']['C_MSE'])),float(e0.mean()))
check('Figure 4 R rectangle loss from saved predictions',near(er.mean(),float(matched['R']['C_MSE'])),float(er.mean()))
check('Figure 4 improved rectangle count',int((er.mean(1)<e0.mean(1)).sum())==9)
positive=np.sort(np.maximum((e0-er).sum(0),0))[::-1];share=positive[:10].sum()/positive.sum()
check('Figure 4 ten-readout positive contribution',round(share*100,2)==92.13,float(share))
perms=rows(P/'RNA_weight_shuffles.csv');cg=1-er.mean()/e0.mean()
check('Figure 4 200 conditional weight controls',len(perms)==200 and max(float(r['contrast_gain']) for r in perms)<cg)

hr=rows(P/'Horlbeck_figure_metrics.csv')
def growth(stage,budget,model,key='sse',repeat=None):
    sub=[r for r in hr if r['stage']==stage and r['budget']==budget and r['fold']!='pooled' and r['repeat']!='pooled' and (repeat is None or r['repeat']==str(repeat))]
    agg=lambda m:sum(float(r[key]) for r in sub if r['model']==m)
    return 100*(1-agg(model)/agg('M0'))
for model,key,expected in [('rank1_selected','sse',6.43),('rank1_selected','projected_sse',6.18),('spectral_selected','sse',8.16),('spectral_selected','projected_sse',5.90)]:
    val=growth('Jurkat','high',model,key)
    check('Figure 5 Jurkat '+model+' '+key,round(val,2)==expected,float(val))
for row in rows(R/'figure_source_data/Fig5_growth_gains.csv'):
    rep={'Repeat 1':1,'Repeat 2':2,'Partner contrast':None}[row['metric']]
    key='projected_sse' if rep is None else 'sse'
    check('Figure 5 exported '+row['model']+' '+row['metric'],near(float(row['gain_percent']),growth('Jurkat','high',row['model'],key,rep)))
gr={r['model']:r for r in rows(G/'SUMMARY.csv')}; primary=['R_G','StandardRidge_RNA','Scouter_Ensemble','GEARS_Ensemble','SingleB_additive']
for key,expected in [('RNA_MSE','StandardRidge_RNA'),('partner_MSE','GEARS_Ensemble'),('mean_regret','R_G')]:
    best=min(primary,key=lambda m:float(gr[m][key]));check('Figure 6 best '+key,best==expected,best)
check('Figure 6 PairPick 19 of 28 hits',int(gr['R_G']['hits'])==19)
ps={r['reference']:r for r in rows(G/'PAIRED_COMPARISONS.csv')}
check('Figure 6 PRESAGE tied recommendation utility',float(ps['PRESAGE_Combo_seed24']['mean_paired_U_difference'])==0)
pool=rows(R/'figure_source_data/PairPick_largest_pool.csv');pick=max(pool,key=lambda r:float(r['R_G']));observed=max(pool,key=lambda r:float(r['observed']))
check('Figure 6 illustrative pool and recommendation',len(pool)==6 and pool[0]['task']=='f1:KLF1-a' and pick['partner']=='GATA1-i' and observed['partner']=='JUN-i',{'candidates':len(pool),'selected':pick['partner'],'observed_best':observed['partner']})

audit=json.loads((R/'FIGURE_AUDIT.json').read_text())
check('Six figures with all text within page',len(audit)==6 and all(not a['text_outside_page'] for a in audit))
check('Combined PDF contains six pages',len(PdfReader(R/'WitnessCell_v24_Main_Figures.pdf').pages)==6)
for a in audit:
    stem=R/'figures'/a['stem'];pdf=PdfReader(stem.with_suffix('.pdf'))
    check(a['stem']+' vector PDF page',len(pdf.pages)==1 and len(pdf.pages[0].extract_text())>50)
    svg=stem.with_suffix('.svg').read_text();check(a['stem']+' SVG retains editable text','<text' in svg)
    im=Image.open(stem.with_suffix('.png'));expected=[round(a['width_inches']*600),round(a['height_inches']*600)]
    check(a['stem']+' 600 dpi PNG',all(abs(x-y)<=1 for x,y in zip(im.size,expected)) and abs(im.info['dpi'][0]-600)<.1,list(im.size))
result={'status':'PASS','checks_passed':len(checks),'source_files_verified':len(manifest),'scope':'Frozen result verification and figure export checks; no retraining, no new biological experiment.','checks':checks}
(R/'FIGURE_NUMERIC_VERIFICATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','checks_passed','source_files_verified']}))
