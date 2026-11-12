from pathlib import Path
import sys,json,importlib.util
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pypdf import PdfWriter
import argparse
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,default=R/'results');ap.add_argument('--out',type=Path,default=R/'figure_replay');args=ap.parse_args()
old=R/'inputs/figures';P=args.out.resolve();O=args.results.resolve();P.mkdir(parents=True,exist_ok=True);(P/'previews').mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('old_figures',old/'code/make_figures.py');f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
audit=[]
def save(fig,name,directory):
 d=P/directory;d.mkdir(exist_ok=True);fig.savefig(d/(name+'.pdf'));fig.savefig(d/(name+'.svg'));fig.savefig(P/'previews'/(name+'.png'),dpi=150)
 fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[]
 for t in fig.findobj(matplotlib.text.Text):
  if not t.get_visible() or not t.get_text():continue
  b=t.get_window_extent(renderer);w,h=fig.bbox.width,fig.bbox.height
  if b.x0< -1 or b.y0< -1 or b.x1>w+1 or b.y1>h+1:bad.append(t.get_text())
 audit.append(dict(name=name,outside=bad));plt.close(fig)
def panel(fig,rect,letter,title,dy=.022):
 ax=fig.add_axes(rect);fig.text(rect[0]-.035,rect[1]+rect[3]+dy,letter,fontweight='bold',fontsize=10,va='bottom');fig.text(rect[0],rect[1]+rect[3]+dy,title,fontweight='bold',fontsize=8.1,va='bottom');return ax
f.panel=panel
# Preserve original main evidence in panels a-f.
s=(old/'provenance/v33/make_figure4_v33.py').read_text().split('f.save(fig,4,')[0];s=s[s.index('fig=plt.figure'):]
repl={'figsize=(7.15,9.6)':'figsize=(7.15,10.8)','[.13,.705,.33,.205]':'[.13,.79,.33,.15]','[.64,.705,.30,.205]':'[.64,.79,.30,.15]', '.64,.657':'.64,.743','[.13,.405,.33,.175]':'[.13,.55,.33,.14]','[.64,.405,.30,.175]':'[.64,.55,.30,.14]', '.13,.345':'.13,.495','(.10,.315)':'(.61,.477)','[.13,.085,.33,.185]':'[.13,.31,.33,.13]','[.64,.085,.30,.185]':'[.64,.31,.30,.13]','(.12,.015)':'(.12,.25)','Shen: retain complete-response benefit':'Shen: retain endpoint benefit','Shen: adjust partner correction':'Shen: regulate partner change'}
for a,b in repl.items():s=s.replace(a,b)
df=pd.read_csv(old/'partner_control/results/FOUR_ARM_BY_CONDITION.csv');colors={'original':f.GRAY,'whole_shrink':f.GOLD,'partner_only':f.TEAL};labels={'original':'Original','whole_shrink':'Whole attenuation','partner_only':'Component control'}
sourceR=R;R=old
exec(compile(s,'preserved_figure4_a_f','exec'))
R=sourceR
han=pd.read_csv(old/'revision/han_calibration/POOLED_SCORES.csv');models=['HeLa__sf12','293T__sf12','A549__sf12','HeLa__v4','293T__v4','A549__v4']
readouts=[('full_raw',f.TEAL,'D',-.13),('full_published',f.BLUE,'o',.13)]
ax=panel(fig,[.16,.075,.30,.13],'g','Han: net gain vs addition')
for i,m in enumerate(models):
 for pn,col,mk,off in readouts:
  z=han[(han.panel==pn)&(han.model==m)&(han['repeat']=='pooled')].set_index('method');g=100*(1-z.loc['target_control','partner_sse']/z.loc['partner_off','partner_sse']);ax.plot(g,i+off,mk,ms=4,color=col,mfc=col if mk=='D' else 'white')
ax.set_yticks(range(6),[m.replace('__',' ').replace('sf12','SF12').replace('v4','V4') for m in models]);ax.set_ylim(5.6,-.6);ax.set_xlim(-.05,.8);ax.axvline(0,color=f.GRAY,lw=.7);ax.set_xlabel('Partner error reduction (%)',fontsize=7);f.style(ax)
ax=panel(fig,[.64,.075,.30,.13],'h','Han: gain after target fitting')
new=pd.read_csv(O/'POOLED_SCORES.csv');sc={'HeLa':f.PURPLE,'293T':f.TEAL}
for j,source in enumerate(['HeLa','293T']):
 for pn,_,mk,off in readouts:
  vals=[new[(new.panel==pn)&(new.source==source)&(new.baseline==b)&(new.policy=='controlled')].partner_gain_pct.iloc[0] for b in ['additive','saturation','local4']]
  ax.plot(np.arange(3)+(j-.5)*.06+off*.15,vals,marker=mk,color=sc[source],mfc=sc[source] if mk=='D' else 'white',lw=.8,ls='-' if mk=='D' else ':',ms=4)
ax.axhline(0,color=f.GRAY,lw=.7);ax.set_xticks(range(3),['Addition','Saturation','Four terms']);ax.tick_params(axis='x',labelsize=6.6);ax.set_ylim(-.055,.77);ax.set_xlim(-.15,2.2);ax.set_ylabel('Additional partner gain (%)');f.style(ax,'y');ax.legend(handles=[Line2D([0],[0],color=sc[x],label=x+' SF12') for x in sc],fontsize=6.2,frameon=False,loc='upper right')
fig.legend(handles=[Line2D([0],[0],marker='D',ls='',color=f.GRAY,label='Reconstructed'),Line2D([0],[0],marker='o',ls='',color=f.GRAY,mfc='white',label='Published')],bbox_to_anchor=(.12,.012),loc='lower left',frameon=False,ncol=2,fontsize=7,handlelength=.8)
fig.text(.61,.022,'Han target: K562; labels identify source models',fontsize=6.7,color=f.GRAY)
save(fig,'Fig4_Functional_Growth','figures')

# Diagnostics are supplementary, preserving the main five-figure scope.
q=pd.read_csv(O/'FUNCTIONAL_REPEATABILITY.csv');orac=pd.read_csv(O/'FUNCTIONAL_ORACLES.csv').set_index('feature_id');ho=pd.read_csv(O/'HAN_ORACLE_SUMMARY.csv')
fig=plt.figure(figsize=(7.15,7.6));ax=panel(fig,[.24,.57,.66,.35],'a','Repeat agreement on matched training and query graphs')
for i,x in enumerate(q.itertuples()):
 ax.plot([x.training_repeat_pearson,x.query_repeat_pearson],[i,i],color='#CCD6DC',lw=.8);ax.plot(x.training_repeat_pearson,i,'o',ms=4,color=f.TEAL);ax.plot(x.query_repeat_pearson,i,'s',ms=4,color=f.BLUE,mfc='white')
ax.set_yticks(range(len(q)),[x.replace('_',' ').replace('high','high').replace('fold','f') for x in q.feature_id]);ax.set_ylim(len(q)-.5,-.5);ax.set_xlim(-.02,.57);ax.set_xlabel('Endpoint-eliminated repeat correlation');f.style(ax);ax.legend(handles=[Line2D([0],[0],marker='o',color=f.TEAL,ls='',label='Training'),Line2D([0],[0],marker='s',color=f.BLUE,mfc='white',ls='',label='Query')],loc='lower right',frameon=False,fontsize=7)
ax=panel(fig,[.15,.17,.31,.25],'b','Shen: scalar opportunity')
for i,name in enumerate(['HeLa','A549','293T']):
 x=orac.loc[name+'_high'];ax.plot(x.actual_gain_pct,i,'D',color=f.TEAL,ms=5);ax.plot(x.bounded_oracle_gain_pct,i,'o',color=f.GRAY,mfc='white',ms=5)
ax.axvline(0,color=f.GRAY,lw=.7);ax.set_yticks(range(3),['HeLa','A549','293T']);ax.set_ylim(2.5,-.5);ax.set_xlim(-.05,.18);ax.set_xlabel('Partner gain (%)');f.style(ax)
ax=panel(fig,[.65,.17,.29,.25],'c','Han: exact fold oracle')
for i,x in enumerate(ho.itertuples()):
 ax.plot(x.actual_gain_pct,i,'D',color=f.TEAL,ms=5);ax.plot(x.fold_oracle_gain_pct,i,'o',color=f.GRAY,mfc='white',ms=5)
ax.set_yticks(range(4),[x.source+' / '+('raw' if x.panel=='full_raw' else 'pub') for x in ho.itertuples()]);ax.set_ylim(3.6,-.5);ax.set_xlim(-.02,.8);ax.set_xlabel('Partner gain (%)');f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker='D',color=f.TEAL,ls='',label='Internal coefficient'),Line2D([0],[0],marker='o',color=f.GRAY,mfc='white',ls='',label='Outcome-informed [0,1] oracle')],bbox_to_anchor=(.12,.06),loc='lower left',ncol=2,frameon=False,fontsize=7)
fig.text(.15,.025,'Jurkat: actual = [0,1] oracle = 6.175%; unconstrained nonnegative oracle = 6.782%.',fontsize=7,color=f.GRAY)
save(fig,'SuppFig6_Direction_Repeatability','supplementary_figures')

rnd=pd.read_csv(O/'RANDOM_POOLED.csv');rs=pd.read_csv(O/'RANDOM_SUMMARY.csv');base=pd.read_csv(O/'BASELINE_COMPARISON.csv')
fig=plt.figure(figsize=(7.15,7.0))
for rect,letter,pn,title in [([.13,.60,.33,.29],'a','full_raw','Reconstructed outcomes'),([.63,.60,.31,.29],'b','full_published','Published phenotypes')]:
 ax=panel(fig,rect,letter,title)
 for source in ['HeLa','293T']:
  vals=rnd[(rnd.panel==pn)&(rnd.source==source)].gain_pct
  ax.hist(vals,bins=np.linspace(-.1,1.05,35),histtype='step',lw=1.3,color=sc[source],label=source)
  t=rs[(rs.panel==pn)&(rs.source==source)].iloc[0];ax.axvline(t.observed_gain_pct,color=sc[source],ls='--',lw=1.1)
 ax.set_xlabel('Partner gain (%)');ax.set_ylabel('Random directions');ax.set_xlim(-.1,1.05);f.style(ax,'y');ax.legend(frameon=False,fontsize=7)
fig.text(.13,.53,'Dashed lines: learned coefficients. Histograms: 999 energy-matched random directions.',fontsize=7,color=f.GRAY)
ax=panel(fig,[.13,.15,.33,.25],'c','Nonlinear baseline gain')
for j,(pn,col,mk,off) in enumerate(readouts):
 vals=[base[(base.panel==pn)&(base.baseline==b)].partner_gain_vs_additive.iloc[0] for b in ['saturation','local4']];ax.bar(np.arange(2)+(j-.5)*.3,vals,width=.28,color=col,label='Reconstructed' if pn=='full_raw' else 'Published')
ax.set_xticks([0,1],['Saturation','Four terms']);ax.set_ylabel('Partner gain vs addition (%)');f.style(ax,'y');ax.legend(frameon=False,fontsize=6.7)
ax=panel(fig,[.66,.15,.28,.25],'d','Increment beyond nonlinear baseline')
for j,source in enumerate(['HeLa','293T']):
 for pn,col,mk,off in readouts:
  vals=[new[(new.panel==pn)&(new.source==source)&(new.baseline==b)&(new.policy=='controlled')].partner_gain_pct.iloc[0] for b in ['saturation','local4']];ax.plot(np.arange(2)+(j-.5)*.12+off*.2,vals,marker=mk,color=sc[source],mfc=sc[source] if mk=='D' else 'white',ls='',ms=4.5)
ax.set_xticks([0,1],['Saturation','Four terms']);ax.axhline(0,color=f.GRAY,lw=.7);ax.set_ylim(-.026,.012);ax.set_xlim(-.3,1.3);ax.set_ylabel('Additional partner gain (%)');f.style(ax,'y')
fig.text(.13,.045,'Nonlinear references and amplitudes use target calibration pairs; source predictors stay frozen.',fontsize=7,color=f.GRAY)
save(fig,'SuppFig7_Han_Attribution','supplementary_figures')
(P/'FIGURE_CHECKS.json').write_text(json.dumps(audit,indent=2))
assert not any(x['outside'] for x in audit),audit
print(audit)
