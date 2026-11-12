"""Rebuild v38 Figures 1,4 and Supplementary Figures 5,6 from fixed source tables."""
from pathlib import Path
import argparse,json,importlib.util
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pypdf import PdfWriter
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'figure_replay');args=ap.parse_args()
P=args.out.resolve();(P/'previews').mkdir(parents=True,exist_ok=True)
old=ROOT/'inputs/figures'
spec=importlib.util.spec_from_file_location('frozen_style',old/'code/make_figures.py');f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
audit=[]
def save(fig,name,directory):
 d=P/directory;d.mkdir(exist_ok=True)
 fig.savefig(d/(name+'.pdf'));fig.savefig(d/(name+'.svg'));fig.savefig(P/'previews'/(name+'.png'),dpi=160)
 fig.canvas.draw();r=fig.canvas.get_renderer();outside=[]
 for t in fig.findobj(matplotlib.text.Text):
  if not t.get_visible() or not t.get_text():continue
  b=t.get_window_extent(r);w,h=fig.bbox.width,fig.bbox.height
  if b.x0< -1 or b.y0< -1 or b.x1>w+1 or b.y1>h+1:outside.append(t.get_text())
 audit.append({'name':name,'outside':outside});plt.close(fig)
def panel(fig,rect,letter,title,dy=.022):
 ax=fig.add_axes(rect);fig.text(rect[0]-.035,rect[1]+rect[3]+dy,letter,fontweight='bold',fontsize=10,va='bottom')
 fig.text(rect[0],rect[1]+rect[3]+dy,title,fontweight='bold',fontsize=8.1,va='bottom');return ax
f.panel=panel
# Fig.1: remove the old transfer-specific box; retain the complete method.
s=(old/'provenance/make_method_v36.py').read_text();s=s[s.index('fig=plt.figure'):s.index("f.save(fig,1,")]
s=s.replace(".03,.69,.94,.25", ".03,.56,.94,.34").replace(".03,.38,.94,.25", ".03,.08,.94,.34")
s="\n".join(x for x in s.splitlines() if 'Function: shared SF12 / V4' not in x)
exec(compile(s,'method_figure','exec'));save(fig,'Fig1_Method','figures')
# Fig.4: earlier scientific evidence remains in a-d.
s=(old/'provenance/v33/make_figure4_v33.py').read_text();s=s[s.index('fig=plt.figure'):s.index('for rect,letter,title,metric,limits')]
repl={'figsize=(7.15,9.6)':'figsize=(7.15,10.8)','[.13,.705,.33,.205]':'[.13,.79,.33,.15]','[.64,.705,.30,.205]':'[.64,.79,.30,.15]', '.64,.657':'.64,.743','[.13,.405,.33,.175]':'[.13,.55,.33,.14]','[.64,.405,.30,.175]':'[.64,.55,.30,.14]', '.13,.345':'.13,.495','(.10,.315)':'(.61,.477)'}
for a,b in repl.items():s=s.replace(a,b)
R=old
exec(compile(s,'preserved_functional_panels','exec'))
native=ROOT/'inputs/native_results';pool=pd.read_csv(native/'POOLED.csv');parts=pd.read_csv(native/'BY_PARTITION.csv')
def row(p,b,m):return pool[(pool.panel==p)&(pool.baseline==b)&(pool.method==m)].iloc[0]
readouts=[('full_raw',f.TEAL,'D',-.16),('full_published',f.BLUE,'o',.16)]
for rect,letter,base,title in [([.13,.31,.33,.13],'e','native','Han: native partner learning'),([.64,.31,.30,.13],'f','saturation','Han: beyond nonlinear baseline')]:
 ax=panel(fig,rect,letter,title)
 for i,fam in enumerate(['rank1','spectral']):
  for pn,col,mk,off in readouts:
   a=row(pn,base,fam+'_original');v=a.partner_gain_vs_own_baseline_pct
   zz=parts[(parts.panel==pn)&(parts.baseline==base)&(parts.method==fam+'_original')].partner_gain_vs_own_baseline_pct
   ax.plot([zz.min(),zz.max()],[i+off,i+off],color=col,lw=1.1,alpha=.65)
   ax.plot(v,i+off,mk,color=col,ms=4.9,mfc=col if mk=='D' else 'white')
   ax.text(min(8.75,max(v,zz.max())+.18),i+off,f'{v:.2f}',fontsize=6.8,va='center',color=col)
 ax.axvline(0,color=f.GRAY,lw=.7);ax.set_yticks([0,1],['Rank one','Spectral']);ax.set_ylim(1.55,-.55);ax.set_xlim(-.35,9.5);ax.set_xticks([0,3,6,9]);ax.set_xlabel('Partner error reduction (%)');f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker='D',ls='',color=f.TEAL,label='Reconstructed'),Line2D([0],[0],marker='o',ls='',color=f.BLUE,mfc='white',label='Published')],bbox_to_anchor=(.12,.247),loc='lower left',frameon=False,ncol=2,fontsize=7,handlelength=.8)
fig.text(.64,.258,'Points: pooled  |  Segments: five-partition range',fontsize=6.3,color=f.GRAY)
df=pd.read_csv(old/'partner_control/results/FOUR_ARM_BY_CONDITION.csv')
colors={'original':f.GRAY,'whole_shrink':f.GOLD,'partner_only':f.TEAL}
ax=panel(fig,[.13,.075,.33,.13],'g','Shen: preserve endpoint benefit')
for j,pol in enumerate(colors):
 vals=[float(df[(df.condition.eq(c+'_high'))&df.policy.eq(pol)].full_gain_pp.iloc[0]) for c in ['HeLa','A549','293T']]
 ax.bar(np.arange(3)+(j-1)*.24,vals,width=.22,color=colors[pol])
ax.set_xticks(range(3),['HeLa','A549','293T']);ax.set_ylim(0,7.5);ax.set_ylabel('Complete error reduction (%)');f.style(ax,'y')
ax=panel(fig,[.69,.075,.25,.13],'h','Han: control a learned correction')
policies=[('original','Original',f.GRAY),('whole_full',r'Whole ($\gamma$)',f.GOLD),('off','Partner off',f.PURPLE),('controlled','Controlled',f.TEAL)]
for i,(key,label,col) in enumerate(policies):
 a=row('full_raw','native','rank1_'+key)
 for metric,marker,fill,off in [('full','o',True,-.08),('partner','s',False,.08)]:
  v=a[metric+'_gain_vs_own_baseline_pct']
  ax.plot(v,i+off,marker,color=col,mfc=col if fill else 'white',ms=4.6)
ax.set_yticks(range(4),[x[1] for x in policies]);ax.tick_params(axis='y',labelsize=6.8);ax.set_ylim(3.6,-.5);ax.set_xlim(-.25,10.5);ax.set_xticks([0,5,10]);ax.set_xlabel('Error reduction vs native baseline (%)',fontsize=7);f.style(ax)
fig.legend(handles=[Rectangle((0,0),1,1,fc=colors[p],label=l) for p,l in [('original','Original'),('whole_shrink','Whole attenuation'),('partner_only','Component control')]],bbox_to_anchor=(.095,.008),loc='lower left',ncol=3,frameon=False,fontsize=6.3,columnspacing=.7,handlelength=.8)
fig.text(.66,.023,'Han h: reconstructed; rank one\nFilled circle: complete; open square: partner',fontsize=6.3,color=f.GRAY,linespacing=1.35)
save(fig,'Fig4_Functional_Growth','figures')
# Supplementary Fig.5: diagnostics only, no reinterpreted source-transfer results.
q=pd.read_csv(ROOT/'inputs/diagnostics/FUNCTIONAL_REPEATABILITY.csv')
orac=pd.read_csv(ROOT/'inputs/diagnostics/FUNCTIONAL_ORACLES.csv').set_index('feature_id')
fig=plt.figure(figsize=(7.15,7.3));ax=panel(fig,[.24,.51,.65,.38],'a','Repeat agreement on matched training and query graphs')
for i,x in enumerate(q.itertuples()):
 ax.plot([x.training_repeat_pearson,x.query_repeat_pearson],[i,i],color='#CCD6DC',lw=.8);ax.plot(x.training_repeat_pearson,i,'o',ms=4,color=f.TEAL);ax.plot(x.query_repeat_pearson,i,'s',ms=4,color=f.BLUE,mfc='white')
ax.set_yticks(range(len(q)),[x.replace('_',' ').replace('fold','f') for x in q.feature_id]);ax.set_ylim(len(q)-.5,-.5);ax.set_xlim(-.02,.57);ax.set_xlabel('Endpoint-eliminated repeat correlation');f.style(ax)
ax.legend(handles=[Line2D([0],[0],marker='o',color=f.TEAL,ls='',label='Training'),Line2D([0],[0],marker='s',color=f.BLUE,mfc='white',ls='',label='Query')],loc='lower right',frameon=False,fontsize=7)
ax=panel(fig,[.24,.17,.65,.22],'b','Shen: fixed-direction scalar opportunity')
for i,name in enumerate(['HeLa','A549','293T']):
 x=orac.loc[name+'_high'];ax.plot(x.actual_gain_pct,i,'D',color=f.TEAL,ms=5);ax.plot(x.bounded_oracle_gain_pct,i,'o',color=f.GRAY,mfc='white',ms=5)
ax.axvline(0,color=f.GRAY,lw=.7);ax.set_yticks(range(3),['HeLa','A549','293T']);ax.set_ylim(2.5,-.5);ax.set_xlim(-.05,.18);ax.set_xlabel('Partner gain (%)');f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker='D',color=f.TEAL,ls='',label='Internal coefficient'),Line2D([0],[0],marker='o',color=f.GRAY,mfc='white',ls='',label='Outcome-informed [0,1] optimum')],bbox_to_anchor=(.16,.058),loc='lower left',ncol=2,frameon=False,fontsize=7)
fig.text(.16,.025,'Jurkat: actual = [0,1] optimum = 6.175%; nonnegative optimum = 6.782%.',fontsize=7,color=f.GRAY)
save(fig,'SuppFig5_Direction_Repeatability','supplementary_figures')
# Supplementary Fig.6: every native learner/baseline/readout and both control states.
fig=plt.figure(figsize=(7.15,5.2))
for rect,letter,pn,title in [([.25,.22,.27,.62],'a','full_raw','Reconstructed'),([.66,.22,.27,.62],'b','full_published','Published')]:
 ax=panel(fig,rect,letter,title,dy=.027)
 for i,(base,fam) in enumerate([('native','rank1'),('native','spectral'),('saturation','rank1'),('saturation','spectral')]):
  col=f.TEAL if base=='native' else f.BLUE
  for pol,mk,fill,off in [('original','o',False,-.13),('controlled','D',True,.13)]:
   a=row(pn,base,fam+'_'+pol);zz=parts[(parts.panel==pn)&(parts.baseline==base)&(parts.method==fam+'_'+pol)].partner_gain_vs_own_baseline_pct
   ax.plot([zz.min(),zz.max()],[i+off,i+off],color=col,lw=1)
   ax.plot(a.partner_gain_vs_own_baseline_pct,i+off,mk,color=col,mfc=col if fill else 'white',ms=4.7)
 ax.axvline(0,color=f.GRAY,lw=.7);ax.axhline(1.5,color='#DEE6EA',lw=.7)
 ax.set_yticks(range(4),['Native / rank one','Native / spectral','Saturation / rank one','Saturation / spectral'] if letter=='a' else ['Rank one','Spectral','Rank one','Spectral'])
 ax.set_ylim(3.65,-.65);ax.set_xlim(-.35,9.5);ax.set_xticks([0,3,6,9]);ax.set_xlabel('Partner error reduction (%)');f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc='white',color=f.GRAY,label='Original correction'),Line2D([0],[0],marker='D',ls='',color=f.GRAY,label='Component control')],bbox_to_anchor=(.20,.075),loc='lower left',ncol=2,frameon=False,fontsize=8)
fig.text(.20,.035,'Points: pooled SSE gains. Segments: ranges across five overlapping partitions.',fontsize=7,color=f.GRAY)
save(fig,'SuppFig6_Han_Native_Controls','supplementary_figures')
(P/'FIGURE_CHECKS.json').write_text(json.dumps(audit,indent=2)+'\n')
assert not any(x['outside'] for x in audit),audit
print(audit)
