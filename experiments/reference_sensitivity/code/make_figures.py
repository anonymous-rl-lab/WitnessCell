"""Publication vectors from archived and v40 matched-reference results."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--paper',type=Path,required=True);args=ap.parse_args();O=args.paper
INK='#233747';TEAL='#087F8C';BLUE='#4477AA';GOLD='#D59435';GRAY='#8998A5';PURPLE='#8567A8'
plt.rcParams.update({'font.family':'Liberation Sans','font.size':8,'axes.labelsize':7.6,'xtick.labelsize':7,'ytick.labelsize':7,'axes.edgecolor':'#9AA8B3','axes.labelcolor':INK,'text.color':INK,'xtick.color':INK,'ytick.color':INK,'axes.linewidth':.65,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False,'figure.facecolor':'white'})
audit=[]
def panel(fig,rect,letter,title,dy=.022):
 ax=fig.add_axes(rect);fig.text(rect[0]-.045,rect[1]+rect[3]+dy,letter,fontweight='bold',fontsize=10,va='bottom');fig.text(rect[0],rect[1]+rect[3]+dy,title,fontweight='bold',fontsize=8,va='bottom');return ax
def style(ax,axis='x'):
 for s in ['top','right']:ax.spines[s].set_visible(False)
 ax.set_axisbelow(True);ax.grid(axis=axis,color='#E8EDF0',lw=.6)
def save(fig,stem,directory):
 d=O/directory;d.mkdir(exist_ok=True,parents=True);fig.savefig(d/(stem+'.pdf'));fig.savefig(d/(stem+'.svg'));p=R/'previews';p.mkdir(exist_ok=True);fig.savefig(p/(stem+'.png'),dpi=150)
 fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
 for t in fig.findobj(matplotlib.text.Text):
  if t.get_visible() and t.get_text():
   b=t.get_window_extent(renderer)
   if b.x0< -1 or b.y0< -1 or b.x1>fig.bbox.width+1 or b.y1>fig.bbox.height+1:outside.append(t.get_text())
 audit.append(dict(figure=stem,outside=outside));plt.close(fig)
F=R/'figure_inputs';growth=pd.read_csv(F/'Horlbeck_figure_metrics.csv');pool=pd.read_csv(F/'HAN_POOLED.csv');parts=pd.read_csv(F/'HAN_BY_PARTITION.csv');new=pd.read_csv(R/'results/MATCHED_REFERENCE_POOLED.csv');load=pd.read_csv(R/'results/GENE_LOADINGS.csv');md=pd.read_csv(R/'results/MODE_DIAGNOSTICS.csv');fold=pd.read_csv(R/'results/MATCHED_REFERENCE_SCORES.csv')
def gg(stage,budget,model,metric,rep='pooled'):
 a=growth[(growth.stage==stage)&(growth.budget==budget)&(growth.fold.astype(str)=='pooled')&(growth['repeat'].astype(str)==str(rep))]
 if a.empty:a=growth[(growth.stage==stage)&(growth.budget==budget)&(growth['repeat'].astype(str)==str(rep))]
 vals=a.groupby('model')[metric].sum();return 100*(1-vals[model]/vals['M0'])
def hr(p,b,m):return pool[(pool.panel==p)&(pool.baseline==b)&(pool.method==m)].iloc[0]
fig=plt.figure(figsize=(7.15,10.8))
ax=panel(fig,[.13,.79,.33,.15],'a','Learning: K562 support')
for model,col in [('rank1_selected',TEAL),('spectral_selected',BLUE)]:
 for metric,mk,fill in [('sse','o',True),('projected_sse','s',False)]:
  vv=[gg('K562',b,model,metric) for b in ['low','high']];shift=(-.045 if col==TEAL else .045)+(-.018 if fill else .018);ax.plot(np.arange(2)+shift,vv,ls='',marker=mk,ms=4.6,color=col,mfc=col if fill else 'white')
ax.axhline(0,color=GRAY,lw=.7);ax.set_xticks([0,1],['421-422','2,317-2,318']);ax.set_xlim(-.15,1.15);ax.set_ylim(-2,16);ax.set_xlabel('Training pairs per fold');ax.set_ylabel('Error reduction vs original B (%)');style(ax,'y')
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=TEAL,label='Rank one'),Line2D([],[],marker='o',ls='',color=BLUE,label='Spectral'),Line2D([],[],marker='o',ls='',color=GRAY,label='Complete'),Line2D([],[],marker='s',ls='',color=GRAY,mfc='white',label='Partner')],loc='upper left',ncol=2,frameon=False,fontsize=6.4,columnspacing=.6,handlelength=.8)
ax=panel(fig,[.64,.79,.30,.15],'b','Learning: Jurkat procedure')
for i,metric in enumerate(['sse','projected_sse']):
 for j in [1,2]:ax.plot(gg('Jurkat','high','rank1_selected',metric,j),i+(j-1.5)*.12,'o',mfc='white',mec=TEAL,ms=4)
 v=gg('Jurkat','high','rank1_selected',metric);ax.plot(v,i,'D',color=TEAL,ms=5);ax.text(v,i+.23,f'{v:.2f}%',ha='center',fontsize=7.5,color=TEAL)
ax.set_yticks([0,1],['Complete','Partner']);ax.set_ylim(1.65,-.5);ax.set_xlim(0,11);ax.set_xlabel('Error reduction vs original B (%)');style(ax)
fig.text(.64,.743,'617 queries; 2,464 training pairs\nFilled: pooled | Open: measured repeats',fontsize=6.7,color=GRAY)
ax=panel(fig,[.16,.55,.30,.14],'c','Single-product reference')
for i,g in enumerate(['K562_high','Jurkat_high']):
 a=new[new.group==g].set_index('method');v=[a.loc[k,'partner_gain_pct'] for k in ['native_rank1','global_rank1']];ax.plot(v,[i,i],color='#B4C2CA',lw=1.2);ax.plot(v[0],i,'o',mfc='white',mec=GRAY,ms=5);ax.plot(v[1],i,'D',color=TEAL,ms=5);ax.text(v[1]+.35,i+.20,f'{v[1]:.2f}',fontsize=7,color=TEAL)
ax.set_yticks([0,1],['K562 high','Jurkat']);ax.set_xlim(-.3,14);ax.set_ylim(1.6,-.6);ax.set_xlabel('Rank-one partner error reduction (%)');style(ax)
fig.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',color=GRAY,label='Original baseline'),Line2D([],[],marker='D',ls='',color=TEAL,label='Global cross-term baseline')],bbox_to_anchor=(.10,.482),loc='lower left',ncol=1,frameon=False,fontsize=6.5,handlelength=.8,labelspacing=.3)
ax=panel(fig,[.64,.55,.30,.14],'d','Attribution: K562 loadings')
q=load[load.configuration=='K562_high_fold0'];ax.scatter(q.single_mean,q.fitted_rank1_loading,s=13,c=np.where(q.mitochondrial==1,TEAL,'#C2CBD1'),edgecolors='none');ax.axhline(0,color=GRAY,lw=.6)
for gene,dx,dy in [('MRPL33',-.025,.018),('MRPL3',.005,.015),('PSMD12',-.13,-.035),('ACTL6A',.005,-.035)]:
 a=q[q.gene==gene].iloc[0];ax.annotate(gene,(a.single_mean,a.fitted_rank1_loading),xytext=(a.single_mean+dx,a.fitted_rank1_loading+dy),fontsize=6.2,color=TEAL if a.mitochondrial else INK)
ax.set_xlabel('Mean single effect');ax.set_ylabel('Fitted rank-one loading');ax.set_ylim(-.29,.30);ax.set_xlim(q.single_mean.min()-.04,q.single_mean.max()+.03);ax.set_xticks([-.6,-.4,-.2,0]);style(ax,'y');fig.text(.64,.495,'Teal: mitochondrial set | Fold 0, |r| = 0.002',fontsize=6.6,color=GRAY)
for rect,letter,base,title in [([.13,.31,.33,.13],'e','native','Transfer: Han native learning'),([.64,.31,.30,.13],'f','saturation','Transfer: nonlinear reference')]:
 ax=panel(fig,rect,letter,title)
 for i,fam in enumerate(['rank1','spectral']):
  for pn,col,mk,off in [('full_raw',TEAL,'D',-.16),('full_published',BLUE,'o',.16)]:
   a=hr(pn,base,fam+'_original');v=a.partner_gain_vs_own_baseline_pct;zz=parts[(parts.panel==pn)&(parts.baseline==base)&(parts.method==fam+'_original')].partner_gain_vs_own_baseline_pct
   ax.plot([zz.min(),zz.max()],[i+off,i+off],color=col,lw=1.1,alpha=.65);ax.plot(v,i+off,mk,color=col,ms=4.9,mfc=col if mk=='D' else 'white');ax.text(min(8.75,max(v,zz.max())+.18),i+off,f'{v:.2f}',fontsize=6.8,va='center',color=col)
 ax.axvline(0,color=GRAY,lw=.7);ax.set_yticks([0,1],['Rank one','Spectral']);ax.set_ylim(1.55,-.55);ax.set_xlim(-.35,9.5);ax.set_xticks([0,3,6,9]);ax.set_xlabel('Partner error reduction (%)');style(ax)
fig.legend(handles=[Line2D([],[],marker='D',ls='',color=TEAL,label='Reconstructed'),Line2D([],[],marker='o',ls='',color=BLUE,mfc='white',label='Published')],bbox_to_anchor=(.12,.247),loc='lower left',frameon=False,ncol=2,fontsize=7,handlelength=.8);fig.text(.64,.258,'Points: pooled | Segments: partition range',fontsize=6.3,color=GRAY)
df=pd.read_csv(F/'FOUR_ARM_BY_CONDITION.csv');colors={'original':GRAY,'whole_shrink':GOLD,'partner_only':TEAL}
ax=panel(fig,[.13,.075,.33,.13],'g','Control: preserve endpoints')
for j,pol in enumerate(colors):ax.bar(np.arange(3)+(j-1)*.24,[float(df[(df.condition==c+'_high')&(df.policy==pol)].full_gain_pp.iloc[0]) for c in ['HeLa','A549','293T']],width=.22,color=colors[pol])
ax.set_xticks(range(3),['HeLa','A549','293T']);ax.set_ylim(0,7.5);ax.set_ylabel('Complete error reduction (%)');style(ax,'y')
ax=panel(fig,[.69,.075,.25,.13],'h','Control: Han rank one')
pols=[('original','Original',GRAY),('whole_full',r'Whole ($\gamma$)',GOLD),('off','Partner off',PURPLE),('controlled','Controlled',TEAL)]
for i,(key,label,col) in enumerate(pols):
 a=hr('full_raw','native','rank1_'+key)
 for met,mk,fill,off in [('full','o',True,-.08),('partner','s',False,.08)]:ax.plot(a[met+'_gain_vs_own_baseline_pct'],i+off,mk,color=col,mfc=col if fill else 'white',ms=4.6)
ax.set_yticks(range(4),[a[1] for a in pols]);ax.tick_params(axis='y',labelsize=6.8);ax.set_ylim(3.6,-.5);ax.set_xlim(-.25,10.5);ax.set_xticks([0,5,10]);ax.set_xlabel('Error reduction vs native baseline (%)',fontsize=7);style(ax)
fig.legend(handles=[Rectangle((0,0),1,1,fc=colors[p],label=l) for p,l in [('original','Original'),('whole_shrink','Whole attenuation'),('partner_only','Component control')]],bbox_to_anchor=(.095,.008),loc='lower left',ncol=3,frameon=False,fontsize=6.3,columnspacing=.7,handlelength=.8);fig.text(.66,.023,'Han h: reconstructed; rank one\nFilled: complete; open: partner',fontsize=6.3,color=GRAY)
save(fig,'Fig4_Functional_Growth','figures')

# All-configuration sensitivity and the original descriptive Jurkat example.
fig=plt.figure(figsize=(7.15,10.5));labels=md.configuration.tolist();short=[x.replace('K562_','K562 ').replace('_fold',' f').replace('Jurkat_high','Jurkat') for x in labels]
for rect,letter,fam,title in [([.22,.70,.28,.22],'a','rank1','Rank-one sensitivity'),([.67,.70,.27,.22],'b','spectral','Spectral sensitivity')]:
 ax=panel(fig,rect,letter,title)
 for i,c in enumerate(labels):
  q=fold[(fold.configuration==c)&(fold.repeat=='pooled')].set_index('method');v=[q.loc[k,'partner_gain_pct'] for k in ['native_'+fam,'global_'+fam]];ax.plot(v,[i,i],color='#C4CED4',lw=.8);ax.plot(v[0],i,'o',mfc='white',mec=GRAY,ms=4);ax.plot(v[1],i,'D',color=TEAL,ms=4)
 ax.axvline(0,color=GRAY,lw=.6);ax.set_yticks(range(len(labels)),short if letter=='a' else []);ax.set_ylim(len(labels)-.5,-.5);ax.set_xlabel('Partner error reduction (%)');style(ax)
fig.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',color=GRAY,label='Original baseline'),Line2D([],[],marker='D',ls='',color=TEAL,label='Expanded global reference')],bbox_to_anchor=(.22,.638),loc='lower left',ncol=2,frameon=False,fontsize=7)
ax=panel(fig,[.22,.405,.28,.20],'c','Fitted and residual modes')
q=md[md.configuration.str.contains('high')]
for i,z in enumerate(q.itertuples()):ax.plot(z.fitted_mode_single_abs_r,i,'D',color=TEAL,ms=4.5);ax.plot(z.residual_mode_single_abs_r,i,'o',mfc='white',mec=BLUE,ms=4.5)
ax.set_yticks(range(4),['K562 f0','K562 f1','K562 f2','Jurkat']);ax.set_ylim(3.5,-.5);ax.set_xlim(-.03,1);ax.set_xlabel('|Pearson r| with mean singles');style(ax);ax.legend(handles=[Line2D([],[],marker='D',ls='',color=TEAL,label='Fitted'),Line2D([],[],marker='o',ls='',mfc='white',color=BLUE,label='Residual')],loc='upper right',frameon=False,fontsize=6.5)
ax=panel(fig,[.67,.405,.27,.20],'d','Mitochondrial-set coverage')
q=pd.read_csv(R/'results/MODULE_GAIN_ATTRIBUTION.csv');q=q[(q.pair_class=='any')]
for i,c in enumerate(['K562_high_fold0','K562_high_fold1','K562_high_fold2']):
 z=q[q.configuration==c].set_index('model');ax.bar(i-.23,100*z.loc['original','query_share'],width=.21,color=GRAY);ax.bar(i,100*z.loc['original','net_gain_share'],width=.21,color=BLUE);ax.bar(i+.23,100*z.loc['global_refit','net_gain_share'],width=.21,color=TEAL)
ax.set_xticks(range(3),['Fold 0','Fold 1','Fold 2']);ax.set_ylim(0,100);ax.set_ylabel('Share (%)');style(ax,'y');ax.legend(handles=[Rectangle((0,0),1,1,fc=c,label=l) for c,l in [(GRAY,'Queries'),(BLUE,'Original gain'),(TEAL,'Expanded gain')]],loc='upper left',frameon=False,fontsize=6)
styles=[('observed',INK,'D','Observed'),('M0',GRAY,'o','Original B'),('rank1_selected',TEAL,'o','Rank one')]
ex=pd.read_csv(F/'JURKAT_FOUR_CORNERS_V28.csv');ex=ex[ex.background=='ACTR10'];ax=panel(fig,[.22,.105,.28,.18],'e','Jurkat: four response corners')
for j,corner in enumerate(['AB','AD','CB','CD']):
 rr=ex[ex.corner==corner]
 for key,col,m,label in styles:ax.plot(j,float(rr[key].mean()),m,color=col,ms=4.7)
ax.set_xticks(range(4),['A+B','A+D','C+B','C+D']);ax.set_xlim(-.4,3.4);ax.set_ylabel('Mean growth response');style(ax,'y')
ax.legend(handles=[Line2D([],[],marker=m,ls='',color=col,label=label) for key,col,m,label in styles],loc='lower center',bbox_to_anchor=(.5,1.01),ncol=3,frameon=False,fontsize=5.7,handlelength=.6,columnspacing=.6)
cc=pd.read_csv(F/'JURKAT_FOUR_CORNER_CONTRASTS_V28.csv');cc=cc[cc.background=='ACTR10'];ax=panel(fig,[.72,.105,.22,.18],'f','Partner-exchange contrast')
for j,(key,col,m,label) in enumerate(styles):
 ax.plot(float(cc[cc['repeat']=='mean'][key].iloc[0]),j,m,color=col,ms=5)
 for rep,off in [('1',-.09),('2',.09)]:ax.plot(float(cc[cc['repeat']==rep][key].iloc[0]),j+off,'o',mfc='white',mec=col,ms=3)
ax.axvline(0,color=GRAY,lw=.6,ls='--');ax.set_yticks(range(3),['Observed','Original B','Rank one']);ax.set_ylim(2.6,-.6);ax.set_xlim(-.075,.008);ax.set_xticks([-.06,-.03,0]);ax.set_xlabel(r'$\Delta$ (growth units)');style(ax)
fig.text(.22,.042,'A = TTI2; C = ACTR10; B = PGK1; D = OGFOD1.\nFour-corner example uses the original baseline and selected correction.',fontsize=7,color=GRAY)
save(fig,'SuppFig7_Global_Reference_Diagnostics','supplementary_figures')

# Preserve S6 values and add the missing baseline labels to both panels.
fig=plt.figure(figsize=(7.15,5.2))
for rect,letter,pn,title in [([.24,.22,.27,.62],'a','full_raw','Reconstructed'),([.71,.22,.24,.62],'b','full_published','Published')]:
 ax=panel(fig,rect,letter,title,dy=.027)
 for i,(base,fam) in enumerate([('native','rank1'),('native','spectral'),('saturation','rank1'),('saturation','spectral')]):
  col=TEAL if base=='native' else BLUE
  for pol,mk,fill,off in [('original','o',False,-.13),('controlled','D',True,.13)]:
   a=hr(pn,base,fam+'_'+pol);z=parts[(parts.panel==pn)&(parts.baseline==base)&(parts.method==fam+'_'+pol)].partner_gain_vs_own_baseline_pct;ax.plot([z.min(),z.max()],[i+off,i+off],color=col,lw=1);ax.plot(a.partner_gain_vs_own_baseline_pct,i+off,mk,color=col,mfc=col if fill else 'white',ms=4.7)
 ax.axvline(0,color=GRAY,lw=.7);ax.axhline(1.5,color='#DEE6EA',lw=.7);ax.set_yticks(range(4),['Native / rank one','Native / spectral','Saturation / rank one','Saturation / spectral']);ax.tick_params(axis='y',labelsize=6.6);ax.set_ylim(3.65,-.65);ax.set_xlim(-.35,9.5);ax.set_xticks([0,3,6,9]);ax.set_xlabel('Partner error reduction (%)',fontsize=7);style(ax)
fig.legend(handles=[Line2D([],[],marker='o',ls='',mfc='white',color=GRAY,label='Original correction'),Line2D([],[],marker='D',ls='',color=GRAY,label='Component control')],bbox_to_anchor=(.20,.075),loc='lower left',ncol=2,frameon=False,fontsize=8);fig.text(.20,.035,'Points: pooled SSE gains. Segments: five-partition ranges.',fontsize=7,color=GRAY)
save(fig,'SuppFig6_Han_Native_Controls','supplementary_figures')
(R/'results/FIGURE_CHECKS.json').write_text(json.dumps(audit,indent=2));assert not any(v['outside'] for v in audit),audit
print(audit)
