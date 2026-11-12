"""Bridge regenerated full-precision tables to all displayed v41 added table rows."""
from pathlib import Path
import argparse,json
import pandas as pd
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
si=(R/'paper/manuscript_sources/WitnessCell_v41_Supplementary_Information.md').read_text();made=[]
fmt=lambda v:f'{0. if abs(v)<5e-9 else v:.4f}'
def row(values):
 t='| '+' | '.join(values)+' |';assert t in si,t;made.append(t)
x=pd.read_csv(a.results/'COMMON_REFERENCE.csv').query("group=='Jurkat_high'").set_index('method')
for name,label in [('native_rank1','Original baseline + original rank one'),('global_B','Expanded global baseline'),('global_rank1','Expanded baseline + refitted rank one'),('global_spectral','Expanded baseline + refitted spectral')]:
 z=x.loc[name];row([label,fmt(z.full_gain_original_pct),fmt(z.partner_gain_original_pct)])
x=pd.read_csv(a.results/'SATURATION_POOLED.csv')
for z in x.itertuples():
 label={'Jurkat_high':'Jurkat','K562_high':'K562 high'}[z.group]
 method={'baseline':'Reference only','rank1':'Rank one','rank1_control':'Controlled rank one','spectral':'Spectral','spectral_control':'Controlled spectral'}[z.method]
 row([label,method,fmt(z.full_gain_pct),fmt(z.partner_gain_pct),fmt(z.partner_gain_original_pct)])
 if z.method=='baseline':row([label,fmt(z.full_sse),fmt(z.partner_sse),fmt(z.full_gain_original_pct),fmt(z.partner_gain_original_pct)])
x=pd.read_csv(a.results/'HAN_ATTRIBUTION_POOLED.csv').query("template=='repeat_specific'")
for z in x.itertuples():
 row([{'full_raw':'Reconstructed','full_published':'Published'}[z.panel],z.baseline,{'rank1':'Rank one','spectral':'Spectral'}[z.family],fmt(z.total_gain_pct),fmt(z.template_gain_pct),fmt(z.orthogonal_gain_pct)])
assert len(made)==24
report=dict(status='PASS',displayed_rows_checked=24,sections=['S54D','S54E','S54F','S57B'],rows=made)
a.out.write_text(json.dumps(report,indent=2)+'\n');print('PASS 24 current display rows')
