"""Figures for the executed v20 extension, drawing only package-local tables."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import build_historical_panels as H
ROOT=Path(__file__).resolve().parents[1]

def rows(p):
 return list(csv.DictReader((ROOT/'source_data'/p).open()))

def ed9():
 d=rows('original_system_extension/DATASET_RESULTS.csv')
 s=rows('original_system_extension/SPLIT_RESULTS.csv')
 z=rows('original_system_extension/SEEN_RESULTS.csv')
 f,aa=plt.subplots(2,2,figsize=(7.7,7.8))
 f.subplots_adjust(left=.19,right=.96,top=.90,bottom=.13,wspace=.67,hspace=.67)
 names=['unit','uniform','amplitude_only','R_amplitude_noise','selected']
 labs=['Unit','Uniform','Amplitude-only','Amplitude-noise','Selected']
 colors={'Norman':'#376FA0','Wessels':'#C98930'}
 for ax,key,letter,title in [(aa[0,0],'non_target_geometric_gain','a','Complete non-target prediction'),(aa[0,1],'endpoint_removed_geometric_gain','b','Endpoint-eliminated prediction')]:
  H.panel(ax,letter,title,'x'); y=np.arange(len(names))
  for i,ds in enumerate(['Norman','Wessels']):
   vals=[100*float(next(r for r in d if r['dataset']==ds and r['method']==m)[key]) for m in names]
   ax.barh(y+(i-.5)*.29,vals,height=.26,color=colors[ds],label=ds)
  ax.set_yticks(y,labs);ax.invert_yaxis();ax.axvline(0,color=H.INK,lw=.8)
  ax.set_xlabel('Error reduction vs original v14 (%)')
  if letter=='a':ax.set_xlim(-3.8,1.15);ax.legend(frameon=False,fontsize=7.6,loc='lower left')
  else:ax.set_xlim(-.08,3.05)
 ax=aa[1,0];H.panel(ax,'c','Training-selected policy: all six splits','x')
 ss=sorted([r for r in s if r['method']=='selected'],key=lambda r:(r['dataset'],int(r['seed'])))
 yy=np.arange(6)
 for k,col,off,lab in [('non_target_gain_vs_original','#376FA0',-.15,'Complete'),('endpoint_removed_gain_vs_original','#8B6CA5',.15,'Projected')]:
  vals=[100*float(r[k]) for r in ss];ax.barh(yy+off,vals,height=.27,color=col,label=lab)
 ax.set_yticks(yy,[r['dataset']+' / '+r['seed'] for r in ss]);ax.invert_yaxis();ax.set_xlim(-.2,3.1);ax.axvline(0,color=H.INK,lw=.8);ax.set_xlabel('Error reduction (%)');ax.legend(frameon=False,fontsize=7.5,loc='lower right')
 ax=aa[1,1];H.panel(ax,'d','Selected policy by measured endpoints','x')
 zs=sorted([r for r in z if r['method']=='selected'],key=lambda r:(r['dataset'],int(r['seen'])))
 yy=np.arange(len(zs));vv=[100*float(r['pooled_gain_vs_original']) for r in zs]
 ax.barh(yy,vv,height=.53,color=[colors[r['dataset']] for r in zs])
 ax.set_yticks(yy,[r['dataset']+' / seen'+r['seen'] for r in zs]);ax.invert_yaxis();ax.axvline(0,color=H.INK,lw=.8);ax.set_xlabel('Pooled complete-error reduction (%)')
 f.text(.13,.035,'Original mean replay: maximum difference 0 in every split. 389 pair-seed records; 224 dataset-specific unique pairs.\nPanels a-c: three-seed or split descriptions. Panel d: pooled descriptive strata. No independent-replicate inference.',fontsize=8,linespacing=1.6)
 H.save(f,'ED9')

if __name__=='__main__':ed9()

def bdata():
 import json
 return json.loads((ROOT/'source_data/independent_RNA/EVALUATION.json').read_text())

def fig6():
 j=bdata();m=j['methods'];old=H.save
 def finish(f,n):
  f.set_size_inches(7.7,11.3)
  for ax in f.axes:
   p=ax.get_position();ax.set_position([p.x0,.28+.72*p.y0,p.width,.72*p.height])
  for leg in f.legends:leg.set_bbox_to_anchor((.55,.315),transform=f.transFigure)
  a=f.add_axes([.14,.065,.34,.17]);H.panel(a,'f','Independent-source RNA: additional gains','x')
  values=[100*m['selected']['full_gain_vs_M0'],100*m['selected']['projected_gain_vs_M0']]
  for k,v in enumerate(values):
   a.barh(k,v,height=.4,color=['#376FA0','#A26473'][k]);a.text(v+(.015 if v>=0 else -.015),k,f'{v:+.3f}%',ha='left' if v>=0 else 'right',va='center',fontsize=8.5)
  a.set_yticks([0,1],['Complete','Projected']);a.invert_yaxis();a.axvline(0,color=H.INK,lw=.8);a.set_xlim(-.19,.47);a.set_ylim(1.65,-.6);a.set_xlabel('Reduction vs RNA M0 (%)')
  a.text(0,-.34,'M0 vs additive: 23.43% lower complete error.\n45 pairs / 10 TFs; primary increment not confirmed.',transform=a.transAxes,fontsize=7.8,linespacing=1.5)
  b=f.add_axes([.66,.065,.29,.17]);H.panel(b,'g','Fixed G2M RNA program','x')
  keys=['additive','zero','selected'];yy=np.arange(3)
  for k,field in enumerate(['program_full_sse','program_projected_sse']):
   vals=[m[z][field]/m['zero'][field] for z in keys];b.barh(yy+(k-.5)*.28,vals,height=.25,color=['#376FA0','#8B6CA5'][k],label=['Complete','Projected'][k])
  b.axvline(1,color=H.INK,lw=.8,ls='--');b.set_yticks(yy,['Additive','M0','Selected']);b.invert_yaxis();b.set_xlim(.72,1.055);b.set_ylim(3.25,-.5);b.set_xlabel('Program error / M0');b.legend(frameon=False,fontsize=7.5,loc='lower left')
  b.text(0,-.34,'52 / 54 predefined genes.\nAdditive retains the lowest complete error.',transform=b.transAxes,fontsize=7.8,linespacing=1.5)
  old(f,'Fig6')
 H.save=finish
 try:H.fig5()
 finally:H.save=old

def ed10():
 j=bdata();m=j['methods'];f,aa=plt.subplots(3,2,figsize=(7.7,9.5));f.subplots_adjust(left=.18,right=.96,top=.91,bottom=.10,wspace=.63,hspace=.88)
 a,b,c,d,e,g=aa.ravel()
 H.panel(a,'a','Complete-response core and extension','x')
 keys=['additive','saturation_only','zero','selected'];labs=['Additive','Saturation-only','M0','Selected']
 vals=[m[z]['full_sse']/m['zero']['full_sse'] for z in keys];a.barh(range(4),vals,height=.53,color=['#B4BFC6','#7F929C','#596D78','#376FA0']);a.set_yticks(range(4),labs);a.invert_yaxis();a.axvline(1,color=H.INK,lw=.8,ls='--');a.set_xlim(.92,1.36);a.set_xlabel('Complete MSE / M0')
 H.panel(b,'b','Primary projected gain: every control','x')
 keys=['unit','uniform','amplitude_only','R_amplitude_noise','selected'];labs=['Unit','Uniform','Amplitude-only','Amplitude-noise','Selected'];vals=[100*m[z]['projected_gain_vs_M0'] for z in keys]
 b.barh(range(5),vals,height=.53,color='#A26473');b.set_yticks(range(5),labs);b.invert_yaxis();b.axvline(0,color=H.INK,lw=.8);b.set_xlim(-.15,.015);b.set_xlabel('Projected reduction vs M0 (%)')
 H.panel(c,'c','Selected outcomes in every fold','x')
 yy=np.arange(3)
 for k,field in enumerate(['selected_full_gain','selected_projected_gain']):c.barh(yy+(k-.5)*.28,[100*r[field] for r in j['fold_results']],height=.25,color=['#376FA0','#8B6CA5'][k],label=['Complete','Projected'][k])
 c.set_yticks(yy,['Fold 0','Fold 1','Fold 2']);c.invert_yaxis();c.axvline(0,color=H.INK,lw=.8);c.set_xlim(-.11,.47);c.set_xlabel('Reduction vs M0 (%)');c.legend(frameon=False,fontsize=7.5,loc='lower right')
 H.panel(d,'d','Prespecified G2M program','x')
 keys=['additive','zero','selected'];yy=np.arange(3)
 for k,field in enumerate(['program_full_sse','program_projected_sse']):d.barh(yy+(k-.5)*.28,[m[z][field]/m['zero'][field] for z in keys],height=.25,color=['#376FA0','#8B6CA5'][k],label=['Complete','Projected'][k])
 d.set_yticks(yy,['Additive','M0','Selected']);d.invert_yaxis();d.axvline(1,color=H.INK,lw=.8,ls='--');d.set_xlim(.72,1.055);d.set_ylim(3.25,-.5);d.set_xlabel('Program error / M0');d.legend(frameon=False,fontsize=7.5,loc='lower left')
 H.panel(e,'e','Delete one endpoint: primary gain','x')
 rr=j['leave_one_endpoint_out'];e.scatter([100*r['projected_gain'] for r in rr],range(10),s=19,color='#8B6CA5');e.set_yticks(range(10),[r['excluded_endpoint'] for r in rr]);e.tick_params(axis='y',labelsize=7.2);e.invert_yaxis();e.axvline(0,color=H.INK,lw=.8);e.set_xlabel('Projected reduction vs M0 (%)')
 H.panel(g,'f','Frozen evaluation identity',None);g.axis('off')
 lines=[('10:10:16','Source + analysis protocol'),('10:21:38','Multiplier distribution'),('10:26:22','Runner + dependencies'),('10:26:30','All outer predictions'),('10:26:36','Scoring completed')]
 for i,(t,txt) in enumerate(lines):g.text(0,.93-i*.17,t+'  '+txt,transform=g.transAxes,fontsize=7.8)
 g.text(0,-.06,'UTC, 19 October 2026\nPrimary resource gate: not passed\nAll results retained; no post-result refitting.',transform=g.transAxes,fontsize=7.8,linespacing=1.6)
 H.save(f,'ED10')
