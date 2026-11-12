#!/usr/bin/env python3
"""Rebuild all v18 figures from package-local source tables and saved contrasts.
No fitting, model selection, external files, network or additional experiments.
Python 3 + numpy + matplotlib. Run from any directory.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import csv,json,hashlib,textwrap
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];S=ROOT/'source_data';OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
INK='#20313E';GRID='#DDE4E9';PALE='#F2F5F7'
COL={'M0':'#7B8790','additive':'#ACB6BE','saturation_only':'#9BABB6','M1':'#326F9E','R':'#B47B27','scalar':'#587688','GEARS':'#32856F','scouter':'#8F6C9F','rank1_selected':'#326F9E','spectral_selected':'#B47B27','M0_repeat_mean':'#9CA6AE'}
NAMES={'additive':'Single additive','saturation_only':'Saturation only','M0':'M0','M1':'M1','R':'R','scalar':'Uniform shrinkage','GEARS':'GEARS','scouter':'scouter','rank1_selected':'Rank one','spectral_selected':'Spectral','M0_repeat_mean':'Repeat-mean M0'}
MODELS=['additive','M0','scouter','GEARS','M1','R'];PART=['rank1_selected','spectral_selected']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.titlesize':9.5,'axes.labelsize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':INK,'xtick.color':INK,'ytick.color':INK,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white'})
def read(n):return list(csv.DictReader((S/n).open()))
def js(n):return json.loads((S/n).read_text())
def panel(a,label,title,axis='y',pad=14):
 offset=min(.13,.04/a.get_position().width); a.text(-offset,1.07,label,transform=a.transAxes,fontsize=12,weight='bold',va='bottom'); title='\n'.join(textwrap.fill(line,width=max(24,int(a.get_position().width*a.figure.get_figwidth()*12))) for line in title.split('\n')); a.set_title(title,loc='left',pad=pad,weight='bold');a.set_axisbelow(True)
 if axis:a.grid(axis=axis,color=GRID,lw=.6)
 a.tick_params(length=3,pad=3)
def save(f,n):
 import io
 for ext in ['pdf','svg','png']:
  stream=io.BytesIO()
  f.savefig(stream,format=ext,dpi=300,bbox_inches='tight',pad_inches=.12)
  path=OUT/f'{n}.{ext}'
  temporary=path.with_suffix(path.suffix+'.tmp')
  temporary.write_bytes(stream.getvalue())
  temporary.replace(path)
 plt.close(f)
def box(a,x,y,w,h,t,color=PALE,edge=GRID,size=8.5):
 a.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.014,rounding_size=0.02',fc=color,ec=edge,lw=.8,transform=a.transAxes,clip_on=False));a.text(x+w/2,y+h/2,t,ha='center',va='center',transform=a.transAxes,fontsize=size,linespacing=1.45)
def arrow(a,x1,y1,x2,y2):a.annotate('',xy=(x2,y2),xytext=(x1,y1),xycoords=a.transAxes,arrowprops={'arrowstyle':'->','lw':1,'color':INK})
def legend_models(a,models,**kw):a.legend(handles=[Line2D([0],[0],color=COL[m],marker='o',ms=4,lw=1.3,label=NAMES[m]) for m in models],frameon=False,**kw)
M={r['model']:r for r in read('Step2_Fig5_Matched55.csv')};SEED=[r for r in read('Step2_Fig5_SeedMetrics.csv') if r['variant']=='raw'];HM=read('Horlbeck_figure_metrics.csv');CONTROL=js('RNA_controls.json');GENE=read('Step2_ED5_GeneGains.csv')
def hval(stage,budget,m,key='gain',fold='pooled',rep='pooled'):
 q=[r for r in HM if r['stage']==stage and r['budget']==budget and r['model']==m and str(r['fold'])==str(fold) and str(r['repeat'])==str(rep)];assert len(q)==1;return float(q[0][key])
def fig1():
 f=plt.figure(figsize=(7.7,8.1));a=f.add_axes([.07,.71,.88,.24]);a.axis('off');panel(a,'a','Learning additional combination responses',None)
 box(a,.01,.36,.21,.43,'Measured singles\n+ training combinations');arrow(a,.23,.58,.28,.69);arrow(a,.23,.58,.28,.28)
 box(a,.30,.57,.29,.25,'Baseline B\nEndpoint + nonlinear response',size=8.4)
 box(a,.30,.16,.29,.25,'Partner correction G\nAdditional pair information',color='#EAF0F6',edge=COL['M1'],size=8.4)
 arrow(a,.60,.69,.69,.69);arrow(a,.60,.28,.69,.49)
 box(a,.71,.38,.27,.44,'Unmeasured pair\nPrediction = B + wG',color='#FAF2E5',edge=COL['R'],size=9)
 a.text(.30,.03,'RNA: response-profile kernel + readout weights',transform=a.transAxes,color=COL['M1'],fontsize=8.1)
 a.text(.30,-.12,'Growth: shared, signed partner matrix across repeats',transform=a.transAxes,color=COL['M1'],fontsize=8.1)
 b=f.add_axes([.09,.37,.38,.21]);b.axis('off');panel(b,'b','Two different prediction questions',None)
 b.text(0,.80,'Complete response',weight='bold',transform=b.transAxes);b.text(0,.60,r'$\|\widehat Y-Y\|^2$',fontsize=15,transform=b.transAxes)
 b.text(0,.31,'After removing free endpoint effects',weight='bold',transform=b.transAxes);b.text(0,.10,r'$\|P(\widehat Y-Y)\|^2,\quad P=I-XX^+$',fontsize=12,transform=b.transAxes)
 b.text(0,-.18,r'Endpoint correction: $\Delta=XB\;\Rightarrow\;P\Delta=0$',fontsize=9,transform=b.transAxes,color=COL['M1'])
 c=f.add_axes([.68,.36,.27,.23]);panel(c,'c','Endpoint gains and\nprojected gains','x',pad=11)
 official=read('Step2_Fig5_Official12.csv');ys=np.arange(4)
 for k,r in enumerate(official):
  g=100*float(r['full_response_gain_vs_joint_backbone']);c.barh(k,g,color=COL['M1'],height=.50);c.scatter(100*float(r['projected_error_gain_vs_joint_backbone']),k,marker='D',color=COL['R'],s=25,zorder=4);c.text(g+.8,k,f'{g:.2f}',va='center',fontsize=8)
 c.set_yticks(ys,[r['dataset'].replace('Replogle_exp6','Replogle exp6') for r in official]);c.invert_yaxis();c.set_xlim(-1,39);c.set_xlabel('Error reduction versus backbone (%)')
 c.legend(handles=[Line2D([0],[0],color=COL['M1'],lw=6,label='Complete'),Line2D([0],[0],color=COL['R'],marker='D',lw=0,label='Projected')],frameon=False,ncol=2,fontsize=7.5,loc='upper center',bbox_to_anchor=(.45,-.29),handlelength=1)
 d=f.add_axes([.10,.08,.85,.18]);panel(d,'d','No additional projected prediction across 12 published splits',None)
 rows=read('Official12_MAIN_TABLE_all_genes_raw.csv');x=np.arange(12)
 d.axhline(1,ls='--',color=INK,lw=.7)
 d.plot(x,[float(r['v14_joint_saturation_backbone']) for r in rows],color=COL['M0'],marker='o',ms=6,mfc='white',lw=1,label='Joint backbone')
 d.scatter(x,[float(r['WitnessCell_v14']) for r in rows],color=COL['M1'],s=20,zorder=4,label='v14 with endpoint correction')
 d.set_xticks(x,[str(r['seed']) for r in rows]);d.set_ylim(.991,1.023);d.set_ylabel('Projected error / truth energy');d.set_xlim(-.5,11.5)
 for xx,lab in [(1,'Norman'),(4,'Wessels'),(7,'Schmidt'),(10,'Replogle exp6')]:d.text(xx,-.22,lab,transform=d.get_xaxis_transform(),ha='center',fontsize=8)
 for xx in [2.5,5.5,8.5]:d.axvline(xx,color=GRID,lw=.7)
 d.legend(frameon=False,ncol=2,fontsize=7.8,loc='upper right');save(f,'Fig1')
def fig2():
 f=plt.figure(figsize=(7.7,8.4));gs=f.add_gridspec(3,2,left=.12,right=.96,bottom=.08,top=.94,wspace=.62,hspace=.65,height_ratios=[.65,1,1.25]);a=f.add_subplot(gs[0,:]);a.axis('off');panel(a,'a','Matched prediction of held-out RNA combinations',None)
 box(a,.01,.30,.22,.48,'55 combinations\n4,995 non-target\nreadouts');arrow(a,.24,.54,.29,.54)
 box(a,.31,.30,.25,.48,'Five outer folds\nSame outer observations');arrow(a,.57,.54,.62,.54)
 box(a,.64,.30,.34,.48,'11 rectangle contrasts\n15 endpoint-free dimensions');a.text(.5,-.10,'M0 → explicit partner kernel M1 → readout-weighted R',ha='center',transform=a.transAxes,color=COL['M1'],fontsize=10)
 b=f.add_subplot(gs[1,0]);panel(b,'b','Combination-contrast prediction','x')
 for k,m in enumerate(MODELS):
  v=100*float(M[m]['C_gain_vs_M0']);b.barh(k,v,color=COL[m],height=.55)
  if m in ['GEARS','scouter']:
   vals=[100*float(r['C_gain_over_M0']) for r in SEED if r['model']==m];b.plot([min(vals),max(vals)],[k,k],color=INK,lw=.8);b.scatter(vals,k+np.linspace(-.1,.1,len(vals)),s=13,fc='white',ec=INK,lw=.5,zorder=5)
  b.text(v+.6,k,f'{v:.2f}',va='center',fontsize=8)
 b.set_yticks(range(6),[NAMES[m] for m in MODELS]);b.invert_yaxis();b.set_xlim(-.4,19);b.set_xlabel('Error reduction versus M0 (%)')
 c=f.add_subplot(gs[1,1]);panel(c,'c','Endpoint-projected prediction','x');c.axvline(1,ls='--',lw=.7,color=INK)
 for k,m in enumerate(MODELS):
  v=float(M[m]['projection_ratio']);lo=float(M[m]['projection_ratio_min']);hi=float(M[m]['projection_ratio_max']);c.plot([lo,hi],[k,k],color=COL[m],lw=1.3);c.scatter(v,k,color=COL[m],s=32,marker='D',zorder=4);c.text(1.035,k,f'{v:.3f}',va='center',fontsize=8)
 c.set_yticks(range(6),[NAMES[m] for m in MODELS]);c.invert_yaxis();c.set_xlim(.84,1.085);c.set_xticks([.85,.9,.95,1]);c.set_xlabel('Projected error / truth energy')
 d=f.add_subplot(gs[2,0]);panel(d,'d','All eleven held-out contrasts',None)
 rect=read('Step2_ED5_Rectangles.csv');keys=sorted({(int(r['fold']),int(r['rectangle'])) for r in rect});disp=['M1','R','GEARS','scouter']
 def rv(k,m):return np.mean([float(r['C_MSE']) for r in rect if (int(r['fold']),int(r['rectangle']))==k and r['model']==m])
 vals=np.array([[100*(1-rv(k,m)/rv(k,'M0')) for m in disp] for k in keys]);im=d.imshow(vals,cmap='RdBu',vmin=-40,vmax=40,aspect='auto')
 d.set_xticks(range(4),disp);d.set_yticks(range(11),[f'{k[0]+1}.{k[1]+1}' for k in keys]);d.set_ylabel('Outer fold.rectangle');d.tick_params(length=0)
 for i in range(11):
  for j in range(4):d.text(j,i,f'{vals[i,j]:.0f}',ha='center',va='center',color='white' if abs(vals[i,j])>25 else INK,fontsize=7.6)
 cb=f.colorbar(im,ax=d,fraction=.05,pad=.04);cb.set_label('Error reduction versus M0 (%)',fontsize=7.5);cb.ax.tick_params(labelsize=7)
 e=f.add_subplot(gs[2,1]);panel(e,'e','Partner gains across outer folds')
 for j,m in enumerate(['M1','R']):
  vals=[100*(1-r['C_MSE'][m]/r['C_MSE']['M0']) for r in CONTROL['folds']];e.bar(np.arange(5)+(j-.5)*.32,vals,width=.30,color=COL[m],label=m)
 e.axhline(0,color=INK,lw=.7);e.set_xticks(range(5),range(1,6));e.set_xlabel('Outer fold');e.set_ylabel('Contrast error reduction versus M0 (%)');e.set_ylim(-3,26);e.legend(frameon=False,ncol=2);save(f,'Fig2')
def rna_corners():
 rows=read('biological_examples/rna/selected_four_corners.csv')
 for case in dict.fromkeys(r['case'] for r in rows):
  rr=[r for r in rows if r['case']==case];assert len(rr)==4
  pp=[set(r['pair'].split('+')) for r in rr]
  ra=next(iter(pp[0]&pp[1]));rc=next(iter(pp[2]&pp[3]));cb=next(iter(pp[0]&pp[2]));cd=next(iter(pp[1]&pp[3]))
  for r,row,col in zip(rr,[ra,ra,rc,rc],[cb,cd,cb,cd]):
   assert set([row,col])==set(r['pair'].split('+'));r['row_endpoint']=row;r['column_endpoint']=col
 return rows
def fig3():
 f,axs=plt.subplots(3,2,figsize=(7.7,10.5));f.subplots_adjust(left=.12,right=.96,top=.94,bottom=.15,hspace=.79,wspace=.58);a,b,c,d,e,g=axs.ravel()
 panel(a,'a','Readout weights add to uniform shrinkage','x')
 for k,m in enumerate(['M1','scalar','R']):
  v=100*CONTROL['C_gain_over_M0'][m];a.barh(k,v,height=.53,color=COL[m]);a.text(v+.45,k,f'{v:.2f}',va='center')
 a.set_yticks(range(3),['M1','Uniform\nshrinkage','R']);a.invert_yaxis();a.set_xlim(0,20);a.set_xlabel('Contrast error reduction versus M0 (%)')
 panel(b,'b','Weights retain readout-specific information');v=[100*float(r['contrast_gain']) for r in read('RNA_weight_shuffles.csv')];b.hist(v,bins=18,color='#AABAC3',ec='white',lw=.5);r=100*CONTROL['C_gain_over_M0']['R'];b.axvline(r,color=COL['R'],lw=2);b.text(r-.4,b.get_ylim()[1]*.88,'R: 15.96%',ha='right',color=COL['R'],weight='bold');b.set_xlim(0,18);b.set_xlabel('Contrast error reduction versus M0 (%)');b.set_ylabel('Weight permutations (n = 200)')
 panel(c,'c','Gains concentrate in a few readouts')
 for m in ['M1','R']:
  vv=sorted([float(r['C_SSE_gain_vs_M0']) for r in GENE if r['model']==m and float(r['C_SSE_gain_vs_M0'])>0],reverse=True);yy=100*np.cumsum(vv)/sum(vv);c.plot(np.arange(1,len(vv)+1),yy,color=COL[m],lw=1.7,label=m);c.scatter(10,yy[9],color=COL[m],s=25);c.annotate(f'{yy[9]:.2f}%',(10,yy[9]),xytext=(7,-15 if m=='M1' else 6),textcoords='offset points',color=COL[m],weight='bold',fontsize=8)
 c.set_xscale('log');c.set_xlim(.9,3000);c.set_ylim(0,104);c.set_xticks([1,10,100,1000],['1','10','100','1,000']);c.axvline(10,color=INK,ls=':',lw=.7);c.set_xlabel('Readouts ranked by positive gain');c.set_ylabel('Cumulative positive gain (%)');c.legend(frameon=False,loc='lower right')
 panel(d,'d','Sensitivity to readout scale');scale=js('RNA_control_scale.json');keys=['condition_mean_variance_positive_quantile_floor_0.01','mean_within_cell_variance_positive_quantile_floor_0.01']
 for j,m in enumerate(['M1','R']):
  vals=[100*float(M[m]['C_gain_vs_M0'])]+[100*scale['scale_sensitivity'][k]['rectangle_gain_over_M0'][m] for k in keys]+[100*float(M[m]['drop10_C_gain_vs_M0'])];xx=np.arange(4)+(j-.5)*.32;d.bar(xx,vals,width=.3,color=COL[m],label=m)
  for x,v in zip(xx,vals):d.text(x,v+(1.1 if v>=0 else -1.1),f'{v:.2f}',ha='center',va='bottom' if v>=0 else 'top',fontsize=7,color=COL[m])
 d.axhline(0,color=INK,lw=.7);d.set_ylim(-15,21);d.set_ylabel('Contrast error reduction versus M0 (%)');d.set_xticks(range(4),['Original','Condition\nvariance','Cell\nvariance','Top ten\nremoved']);d.legend(frameon=False,ncol=2,loc='upper right',fontsize=8)
 corners=rna_corners();examples=js('biological_examples/rna/RNA_EXAMPLES.json')['example_records']
 for ax,letter,case in [(e,'e','positive_S100A9'),(g,'f','positive_LYZ')]:
  rows=[r for r in corners if r['case']==case];ex=next(r for r in examples if r['case']==case)
  panel(ax,letter,ex['readout']+': four held-out responses',None);ax.axhline(0,color=GRID,lw=.7)
  for m,marker,color in [('observed','x',INK),('M0','s',COL['M0']),('M1','^',COL['M1']),('R','D',COL['R'])]:
   vals=[float(r[m+'_effect']) for r in rows]
   for xx in [[0,1],[2.6,3.6]]:
    sl=slice(0,2) if xx[0]==0 else slice(2,4);ax.plot(xx,vals[sl],color=color,marker=marker,ms=5 if m=='observed' else 4,lw=1.7 if m=='observed' else 1.1,ls='--' if m=='M0' else '-',mfc='white' if m=='M0' else color,label=('Observed' if m=='observed' else m) if xx[0]==0 else None)
  ax.set_xticks([0,1,2.6,3.6],[r['column_endpoint'] for r in rows]);ax.tick_params(axis='x',labelsize=7.3);ax.set_xlim(-.35,3.95)
  for x,r in [(0.5,rows[0]),(3.1,rows[2])]:ax.text(x,-.23,r['row_endpoint']+' background',transform=ax.get_xaxis_transform(),ha='center',fontsize=7.6)
  ax.axvline(1.8,color=GRID,lw=.6,ls=':');ax.set_ylabel('RNA response relative to control')
  ax.text(.5,-.44,f"C: observed {ex['observed_contrast']:+.3f} | M0 {ex['M0_contrast']:+.3f}\nM1 {ex['M1_contrast']:+.3f} | R {ex['R_contrast']:+.3f}",transform=ax.transAxes,ha='center',va='top',fontsize=7.7,linespacing=1.4)
  if letter=='e':handles,labels=ax.get_legend_handles_labels()
 f.legend(handles,labels,frameon=False,ncol=4,loc='lower center',bbox_to_anchor=(.55,.005),fontsize=8);save(f,'Fig3')
def budgetpanel(a,key):
 for b,budget in enumerate(['low','high']):
  for j,m in enumerate(PART):
   x=b+(j-.5)*.30;vv=[100*hval('K562',budget,m,key,fold=k,rep=t) for k in range(3) for t in [1,2]];a.scatter(x+np.linspace(-.055,.055,6),vv,fc='white',ec=COL[m],s=22,lw=.8,zorder=3);v=100*hval('K562',budget,m,key);a.scatter(x,v,color=COL[m],s=50,marker='D',edgecolors='white',lw=.5,zorder=5);a.annotate(f'{v:.2f}',(x,v),xytext=(-8 if j==0 else 8,12 if budget=='high' else -17),textcoords='offset points',ha='center',color=COL[m],weight='bold',fontsize=8.5)
 a.axhline(0,color=INK,lw=.7);a.set_xticks([0,1],['421–422 pairs','2,317–2,318 pairs']);a.set_xlim(-.4,1.4);a.set_ylim(-7,31);a.set_ylabel('Error reduction versus matched M0 (%)');legend_models(a,PART,loc='upper left',fontsize=8)
def fig4():
 f,axs=plt.subplots(2,2,figsize=(7.7,6.9));f.subplots_adjust(left=.12,right=.96,bottom=.12,top=.90,wspace=.58,hspace=.73);a,b,c,d=axs.ravel();a.axis('off');panel(a,'a','Same prediction task, more observations',None)
 box(a,.02,.68,.95,.22,'80 genes · two growth repeats',size=9)
 box(a,.02,.25,.43,.32,'Low budget\n421–422\ntraining pairs');box(a,.55,.25,.43,.32,'High budget\n2,317–2,318\ntraining pairs',color='#EAF0F6',edge=COL['M1']);arrow(a,.46,.41,.53,.41)
 a.text(.5,.13,'Nested training sets; refit every model',ha='center',transform=a.transAxes,fontsize=8.3)
 box(a,.02,-.18,.95,.23,'Same 632 evaluation pairs\nacross three folds',size=8.7)
 panel(b,'b','Baseline estimation also improves')
 for fold in range(3):b.plot([0,1],[1000*hval('K562',bud,'M0','mse',fold=fold) for bud in ['low','high']],color='#ADB8C0',marker='o',mfc='white',ms=4,lw=1)
 vals=[1000*hval('K562',bud,'M0','mse') for bud in ['low','high']];b.plot([0,1],vals,color=COL['M0'],marker='D',ms=6,lw=2)
 for x,v in enumerate(vals):b.annotate(f'{v:.3f}',(x,v),xytext=(0,9),textcoords='offset points',ha='center',weight='bold')
 b.text(.96,.40,f'{100*(1-vals[1]/vals[0]):.2f}% lower',transform=b.transAxes,ha='right',weight='bold');b.set_xticks([0,1],['Low budget','High budget']);b.set_xlim(-.3,1.3);b.set_ylim(0,5.4);b.set_ylabel('Complete-response MSE (×10⁻³)')
 panel(c,'c','Additional complete-response prediction');budgetpanel(c,'gain')
 panel(d,'d','Additional endpoint-projected prediction');budgetpanel(d,'projected_gain');save(f,'Fig4')
def fig5():
 f=plt.figure(figsize=(7.7,8.6));gs=f.add_gridspec(3,2,left=.12,right=.96,bottom=.15,top=.93,wspace=.58,hspace=.74,height_ratios=[.7,1,1]);a=f.add_subplot(gs[0,:]);a.axis('off');panel(a,'a','Confirming a frozen learning workflow in Jurkat',None)
 for x,w,t,c in [(0,.19,'K562\nDevelopment',PALE),(.26,.20,'Freeze\nLearning rules','#EAF0F6'),(.53,.20,'Jurkat training\n2,464 pairs',PALE),(.80,.20,'Jurkat test\n617 pairs','#FAF2E5')]:box(a,x,.25,w,.56,t,color=c,size=8.5)
 for x1,x2 in [(.20,.25),(.47,.52),(.74,.79)]:arrow(a,x1,.53,x2,.53)
 a.text(.5,-.13,'79 genes · two independent growth screens · parameters fitted in Jurkat',ha='center',transform=a.transAxes,fontsize=8.3)
 b=f.add_subplot(gs[1,0]);panel(b,'b','Complete-response prediction')
 c=f.add_subplot(gs[1,1]);panel(c,'c','Endpoint-projected prediction')
 for ax,key in [(b,'gain'),(c,'projected_gain')]:
  for j,m in enumerate(PART):
   for k,rep in enumerate([1,2,'pooled']):
    v=100*hval('Jurkat','high',m,key,rep=rep);x=k+(j-.5)*.28;ax.scatter(x,v,color=COL[m] if rep=='pooled' else 'white',edgecolors=COL[m],marker='D' if rep=='pooled' else 'o',s=46 if rep=='pooled' else 33,lw=1.1);ax.annotate(f'{v:.2f}',(x,v),xytext=(0,9 if j==0 else -15),textcoords='offset points',ha='center',color=COL[m],fontsize=8,weight='bold' if rep=='pooled' else 'normal')
  ax.axhline(0,color=INK,lw=.7);ax.set_xticks(range(3),['Repeat 1','Repeat 2','Pooled']);ax.set_xlim(-.4,2.4);ax.set_ylim(0,12);ax.set_ylabel('Error reduction versus M0 (%)')
 d=f.add_subplot(gs[2,0]);panel(d,'d','TTI2: changing the partner')
 rows=read('biological_examples/jurkat/Jurkat_selected_pairs_two_repeats.csv')
 def pair(k):return sorted([r for r in rows if r['pair_id']==k],key=lambda r:int(r['repeat']))
 mods=[('observed',INK,'x'),('M0',COL['M0'],'s'),('rank1_selected',COL['M1'],'D')]
 for shift,key in [(0,'PGK1|TTI2'),(2.6,'OGFOD1|TTI2')]:
  for rep,r in enumerate(pair(key)):
   for j,(m,color,marker) in enumerate(mods):d.scatter(shift+rep+(j-1)*.17,float(r[m]),color=color,marker=marker,s=35 if m=='observed' else 25,facecolors='white' if m=='M0' else color,lw=1,zorder=5)
   d.hlines(float(r['observed']),shift+rep-.33,shift+rep+.33,color=INK,lw=.6,alpha=.55)
 d.set_xticks([0,1,2.6,3.6],['1','2','1','2']);d.set_xlim(-.55,4.15);d.set_ylabel('Growth phenotype, γ');d.set_xlabel('Growth-screen repeat')
 for x,t in [(.5,'TTI2 + PGK1'),(3.1,'TTI2 + OGFOD1')]:d.text(x,-.30,t,transform=d.get_xaxis_transform(),ha='center',fontsize=8)
 d.axvline(1.8,color=GRID,lw=.65,ls=':');d.set_ylim(-.39,.015)
 e=f.add_subplot(gs[2,1]);panel(e,'e','MED1 + RPL10: both repeats improve')
 for rep,r in enumerate(pair('MED1|RPL10')):
  for j,(m,color,marker) in enumerate(mods):e.scatter(rep+(j-1)*.17,float(r[m]),color=color,marker=marker,s=37 if m=='observed' else 27,facecolors='white' if m=='M0' else color,lw=1,zorder=5)
  e.hlines(float(r['observed']),rep-.34,rep+.34,color=INK,lw=.6,alpha=.55)
 e.set_xticks([0,1],['Repeat 1','Repeat 2']);e.set_xlim(-.5,1.5);e.set_ylim(-.345,-.293);e.set_ylabel('Growth phenotype, γ')
 b.legend(handles=[Line2D([0],[0],color=COL[m],marker='D',lw=1.2,label=NAMES[m]+(' (primary)' if i==0 else '')) for i,m in enumerate(PART)],frameon=False,fontsize=7.2,loc='upper left')
 f.legend(handles=[Line2D([0],[0],color=color,marker=marker,lw=0,label={'observed':'Observed','M0':'M0','rank1_selected':'M0 + rank-one partner'}[m],mfc='white' if m=='M0' else color) for m,color,marker in mods],frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(.55,.012),fontsize=8);save(f,'Fig5')
def ed1():
 f,axs=plt.subplots(2,2,figsize=(7.7,6.5));f.subplots_adjust(left=.12,right=.97,top=.9,bottom=.1,wspace=.65,hspace=.68)
 files=['Official12_MAIN_TABLE_all_genes_raw.csv','Official12_MAIN_TABLE_exclude_panel_targets.csv','Official12_MAIN_TABLE_training_variance_standardized.csv']
 for ax,letter,n,title in zip(axs.ravel()[:3],['a','b','c'],files,['Original readout axis','Perturbed readouts excluded','Training variance standardized']):
  panel(ax,letter,title,None);rows=read(n);x=np.arange(12);ax.axhline(1,color=INK,ls='--',lw=.7)
  for key,color,marker in [('v14_joint_saturation_backbone',COL['M0'],'o'),('WitnessCell_v14',COL['M1'],'D')]:ax.plot(x,[float(r[key]) for r in rows],color=color,marker=marker,ms=3.5,lw=.8,mfc='white' if marker=='o' else color)
  ax.set_xticks([1,4,7,10],['Norman','Wessels','Schmidt','Replogle'],rotation=20);ax.set_ylabel('Projected error / truth energy')
 d=axs.ravel()[3];panel(d,'d','Full-response gains in Step 1','x');rows=read('Fig2_prediction_ratios.csv');rows=[r for r in rows if r['pool']=='all_queries' and r['metric']=='loss_wmse'];yy=np.arange(len(rows));v=[100*float(r['gain']) for r in rows]
 d.scatter(v,yy,color=COL['M1'],s=30);d.axvline(0,color=INK,lw=.7);d.set_yticks(yy,[r['dataset']+' / '+('additive' if r['comparator']=='shared_endpoint_additive' else 'saturation') for r in rows]);d.invert_yaxis();d.set_xlabel('Weighted error reduction (%)');save(f,'ED1')
def ed2():
 f,axs=plt.subplots(2,2,figsize=(7.7,6.8));f.subplots_adjust(left=.12,right=.96,top=.9,bottom=.12,wspace=.62,hspace=.72);a,b,c,d=axs.ravel();panel(a,'a','Complete-response error','x')
 for k,m in enumerate(MODELS):
  v=1000*float(M[m]['response_MSE']);a.scatter(v,k,color=COL[m],s=33)
  if m in ['GEARS','scouter']:
   vv=[1000*float(r['response_MSE']) for r in SEED if r['model']==m];a.plot([min(vv),max(vv)],[k,k],color=COL[m],lw=1.4)
 a.set_xscale('log');a.set_yticks(range(6),[NAMES[m] for m in MODELS]);a.invert_yaxis();a.set_xlabel('Complete-response MSE (×10⁻³; log scale)')
 panel(b,'b','Complete-response gains across folds')
 folds=read('Step2_ED5_Folds.csv')
 for j,m in enumerate(['M1','R']):
  vals=[]
  for k in range(5):
   def ss(mm):return np.mean([float(r['response_SSE']) for r in folds if int(r['fold'])==k and r['model']==mm])
   vals.append(100*(1-ss(m)/ss('M0')))
  b.bar(np.arange(5)+(j-.5)*.32,vals,.30,color=COL[m],label=m)
 b.axhline(0,color=INK,lw=.7);b.set_xticks(range(5),range(1,6));b.set_xlabel('Outer fold');b.set_ylabel('Error reduction versus M0 (%)');b.legend(frameon=False,ncol=2)
 panel(c,'c','GEARS: fold contributions to total error');panel(d,'d','scouter: fold contributions to total error')
 for ax,m in [(c,'GEARS'),(d,'scouter')]:
  for j,seed in enumerate([0,1,2]):
   rr=[r for r in folds if r['model']==m and int(r['seed'])==seed];vv=[float(next(r['response_SSE'] for r in rr if int(r['fold'])==k)) for k in range(5)];vv=np.array(vv)/sum(vv)*100;ax.plot(range(1,6),vv,marker=['o','s','D'][j],lw=1.2,ms=4,label=f'Seed {seed}')
  ax.set_xticks(range(1,6));ax.set_xlabel('Outer fold');ax.set_ylabel('Share of complete-response SSE (%)');ax.set_ylim(0,101);ax.legend(frameon=False,fontsize=7.6)
 save(f,'ED2')
def ed3():
 f,axs=plt.subplots(2,2,figsize=(7.7,7));f.subplots_adjust(left=.12,right=.96,top=.9,bottom=.12,wspace=.62,hspace=.72);a,b,c,d=axs.ravel();panel(a,'a','Largest positive contributors','x');pos=sorted([r for r in GENE if r['model']=='R' and float(r['C_SSE_gain_vs_M0'])>0],key=lambda r:float(r['C_SSE_gain_vs_M0']),reverse=True);total=sum(float(r['C_SSE_gain_vs_M0']) for r in pos);v=[100*float(r['C_SSE_gain_vs_M0'])/total for r in pos[:10]];a.barh(range(10),v,color=COL['R']);a.set_yticks(range(10),[r['gene'] for r in pos[:10]]);a.invert_yaxis();a.set_xlabel('Share of all positive gain (%)')
 panel(b,'b','Positive and negative contributions')
 for j,m in enumerate(['M1','R']):
  v=np.array([float(r['C_SSE_gain_vs_M0']) for r in GENE if r['model']==m]);b.bar(j-.15,v[v>0].sum(),.28,color=COL[m]);b.bar(j+.15,v[v<0].sum(),.28,color=COL[m],alpha=.4);b.text(j,1.04*v[v>0].sum(),f'{sum(v>0):,} improved',ha='center',fontsize=8)
 b.axhline(0,color=INK,lw=.7);b.set_xticks([0,1],['M1','R']);b.set_ylabel('Sum of readout contributions to SSE gain');b.text(.03,.05,'Pale: increased squared error',transform=b.transAxes,fontsize=7.6)
 panel(c,'c','Repairing the correct increment','x');sup=js('RNA_supplement_controls.json');v=[100*CONTROL['C_gain_over_M0']['R'],100*sup['additive_control']['additive_weighted']['C_gain_over_M0'],100*sup['random_control']['random_weighted']['C_gain_over_M0_quantiles_0_025_50_975_100'][2]];c.barh(range(3),v,color=[COL['R'],COL['M0'],COL['additive']]);c.set_yticks(range(3),['Partner increment','M0 − additive increment','Random increment']);c.invert_yaxis();c.axvline(0,color=INK,lw=.7);c.set_xlabel('Contrast error reduction versus M0 (%)')
 panel(d,'d','Uniform shrinkage across outer folds');d.plot(range(1,6),CONTROL['scalar_c'],color=COL['scalar'],marker='o',lw=1.5);d.set_xticks(range(1,6));d.set_ylim(0,1.08);d.set_xlabel('Outer fold');d.set_ylabel('Training-selected coefficient c');save(f,'ED3')
def ed4():
 f,axs=plt.subplots(2,2,figsize=(7.7,6.8));f.subplots_adjust(left=.12,right=.96,top=.9,bottom=.12,wspace=.58,hspace=.72)
 for ax,letter,budget,key,title in zip(axs.ravel(),['a','b','c','d'],['low','high','low','high'],['gain','gain','projected_gain','projected_gain'],['Low budget: complete response','High budget: complete response','Low budget: endpoint-projected','High budget: endpoint-projected']):
  panel(ax,letter,title);xx=np.arange(6)
  for j,m in enumerate(['M0_repeat_mean']+PART):
   vals=[100*hval('K562',budget,m,key,fold=k,rep=t) for k in range(3) for t in [1,2]];ax.plot(xx,vals,marker=['s','o','D'][j],color=COL[m],lw=1.2,ms=4,label=NAMES[m])
  ax.axhline(0,color=INK,lw=.7);ax.set_xticks(xx,['1.1','1.2','2.1','2.2','3.1','3.2']);ax.set_xlabel('Outer fold.repeat');ax.set_ylabel('Error reduction versus M0 (%)');
 f.legend(handles=[Line2D([0],[0],color=COL[m],marker='o',lw=1.2,label=NAMES[m]) for m in ['M0_repeat_mean']+PART],frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(.55,.01),fontsize=8);save(f,'ED4')
def ed5():
 f,axs=plt.subplots(1,3,figsize=(7.7,3.4));f.subplots_adjust(left=.09,right=.97,top=.79,bottom=.30,wspace=.65);a,b,c=axs
 panel(a,'a','Gene-composition reweighting','x');interval=read('Horlbeck_Jurkat_intervals.csv')
 for j,m in enumerate(PART):
  rr=next(r for r in interval if r['model']==m);lo=100*float(rr['lower']);hi=100*float(rr['upper']);v=100*hval('Jurkat','high',m);a.plot([lo,hi],[j,j],color=COL[m],lw=2);a.scatter(v,j,color=COL[m],marker='D',s=36)
 a.axvline(0,color=INK,lw=.7);a.set_yticks([0,1],['Rank one','Spectral']);a.set_ylim(1.6,-.6);a.set_xlabel('Error reduction versus M0 (%)')
 panel(b,'b','Averaging the baseline alone')
 for j,m in enumerate(['M0_repeat_mean']+PART):
  vals=[100*hval('Jurkat','high',m,key) for key in ['gain','projected_gain']];b.plot([0,1],vals,color=COL[m],marker='o',lw=1.2,ms=4,label=NAMES[m])
 b.axhline(0,color=INK,lw=.7);b.set_xticks([0,1],['Complete','Projected']);b.set_ylabel('Error reduction versus M0 (%)');b.legend(frameon=False,fontsize=7.2,loc='upper center',bbox_to_anchor=(.5,-.30))
 panel(c,'c','All deletion sensitivity values','x');deletion=read('Horlbeck_Jurkat_gene_deletion.csv')
 for j,m in enumerate(PART):
  vv=np.array([100*float(r['gain']) for r in deletion if r['model']==m]);c.scatter(vv,j+np.linspace(-.15,.15,len(vv)),s=9,color=COL[m],alpha=.65);c.scatter([vv.min(),vv.max()],[j,j],marker='|',s=70,color=INK)
 c.axvline(0,color=INK,lw=.7);c.set_yticks([0,1],['Rank one','Spectral']);c.set_ylim(1.5,-.5);c.set_xlim(0,10);c.set_xlabel('Gain after gene deletion (%)');save(f,'ED5')
def ed6():
 f,axs=plt.subplots(2,2,figsize=(7.7,7.3));f.subplots_adjust(left=.12,right=.96,top=.91,bottom=.13,wspace=.58,hspace=.77)
 rows=read('biological_examples/rna/selected_readouts_all_11_rectangles.csv')
 for ax,letter,gene in [(axs[0,0],'a','S100A9'),(axs[0,1],'b','LYZ')]:
  panel(ax,letter,gene+': all eleven contrasts',None);rr=sorted([r for r in rows if r['gene']==gene],key=lambda r:int(r['rectangle']))
  for m,mark,col in [('observed','x',INK),('M0','s',COL['M0']),('M1','^',COL['M1']),('R','D',COL['R'])]:ax.plot(range(1,12),[float(r['C_'+m]) for r in rr],color=col,marker=mark,ms=3.5,lw=1.2,ls='--' if m=='M0' else '-',label='Observed' if m=='observed' else m)
  ax.axhline(0,color=GRID,lw=.7);ax.set_xticks([1,3,5,7,9,11]);ax.set_xlabel('Archived rectangle');ax.set_ylabel('Signed RNA contrast')
 a=axs[1,0];panel(a,'c','LYZ: a harmful correction',None);rr=[r for r in rna_corners() if r['case']=='failure_global_largest_harm']
 for m,mark,col in [('observed','x',INK),('M0','s',COL['M0']),('M1','^',COL['M1']),('R','D',COL['R'])]:
  for k,xx in enumerate([[0,1],[2.6,3.6]]):a.plot(xx,[float(r[m+'_effect']) for r in rr[k*2:k*2+2]],color=col,marker=mark,ms=4,lw=1.2,ls='--' if m=='M0' else '-')
 a.set_xticks([0,1,2.6,3.6],[r['column_endpoint'] for r in rr]);a.tick_params(axis='x',labelsize=7);a.set_xlim(-.35,3.95);a.set_ylabel('RNA response relative to control');a.axhline(0,color=GRID,lw=.7);a.axvline(1.8,color=GRID,lw=.7,ls=':')
 for x,r in [(.5,rr[0]),(3.1,rr[2])]:a.text(x,-.23,r['row_endpoint']+' background',transform=a.get_xaxis_transform(),ha='center',fontsize=7.4)
 a=axs[1,1];panel(a,'d','TYROBP: all eleven contrasts',None);rr=read('RNA_displayed_signed_examples.csv');rr=[r for r in rr if r['gene']=='TYROBP']
 for m,mark,col in [('truth_effect','x',INK),('M0','s',COL['M0']),('M1','^',COL['M1']),('R','D',COL['R'])]:a.plot(range(1,12),[float(r[m]) for r in rr],color=col,marker=mark,ms=3.5,lw=1.2,ls='--' if m=='M0' else '-')
 a.axhline(0,color=GRID,lw=.7);a.set_xticks([1,3,5,7,9,11]);a.set_xlabel('Archived rectangle');a.set_ylabel('Signed RNA contrast')
 handles,labels=axs[0,0].get_legend_handles_labels();f.legend(handles,labels,frameon=False,ncol=4,loc='lower center',bbox_to_anchor=(.54,.013),fontsize=8);save(f,'ED6')
def ed7():
 f,axs=plt.subplots(2,2,figsize=(7.7,7.5));f.subplots_adjust(left=.12,right=.96,top=.91,bottom=.14,wspace=.58,hspace=.75);a,b,c,d=axs.ravel();allj=read('biological_examples/jurkat/Jurkat_all_617_pair_scores.csv')
 panel(a,'a','All 617 held-out combinations',None);xx=np.array([np.sqrt(float(r['M0_SSE'])/2) for r in allj]);yy=np.array([np.sqrt(float(r['rank1_selected_SSE'])/2) for r in allj]);a.scatter(xx,yy,s=11,color='#9BA8B3',alpha=.6,edgecolors='none');lim=max(xx.max(),yy.max())*1.05;a.plot([0,lim],[0,lim],color=INK,ls='--',lw=.8);a.set_xlim(-.004,lim);a.set_ylim(-.004,lim);a.set_xlabel('M0: pair RMSE');a.set_ylabel('M0 + partner: pair RMSE');a.text(.04,.96,'349/617 pairs improve',transform=a.transAxes,va='top',fontsize=8)
 for pair,offset in [('MED1|RPL10',(10,-8)),('NIP7|POLN',(-59,11))]:
  r=next(r for r in allj if r['pair_id']==pair);x=np.sqrt(float(r['M0_SSE'])/2);y=np.sqrt(float(r['rank1_selected_SSE'])/2);a.scatter(x,y,s=31,color=COL['R'],edgecolor='white',lw=.5,zorder=5);a.annotate(pair.replace('|',' + '),(x,y),xytext=offset,textcoords='offset points',fontsize=6.8,arrowprops={'arrowstyle':'-','lw':.6,'color':'#647786'})
 panel(b,'b','Direction of per-pair error changes');both=sum(r['rank1_both_repeats_improved']=='True' for r in allj);bad=sum(r['rank1_both_repeats_worsened']=='True' for r in allj);counts=[both,len(allj)-both-bad,bad];b.bar(range(3),counts,color=[COL['M1'],COL['M0'],COL['R']])
 for x,v in enumerate(counts):b.text(x,v+5,str(v),ha='center',fontsize=8)
 b.set_xticks(range(3),['Both\nimprove','Mixed\nor tied','Both\nworsen']);b.set_ylabel('Held-out combinations');b.set_ylim(0,max(counts)*1.17)
 panel(c,'c','NIP7 + POLN: both repeats worsen');rows=sorted([r for r in read('biological_examples/jurkat/Jurkat_selected_pairs_two_repeats.csv') if r['pair_id']=='NIP7|POLN'],key=lambda r:int(r['repeat']));mods=[('observed',INK,'x'),('M0',COL['M0'],'s'),('rank1_selected',COL['M1'],'D')]
 for rep,r in enumerate(rows):
  for j,(m,color,marker) in enumerate(mods):c.scatter(rep+(j-1)*.17,float(r[m]),color=color,marker=marker,s=35,facecolors='white' if m=='M0' else color,lw=1,zorder=5)
  c.hlines(float(r['observed']),rep-.34,rep+.34,color=INK,lw=.6,alpha=.6)
 c.set_xticks([0,1],['Repeat 1','Repeat 2']);c.set_xlim(-.5,1.5);c.set_ylim(-.53,-.38);c.set_ylabel('Growth phenotype, γ')
 panel(d,'d','Performance relative to zero contrast','x');ratio=read('Horlbeck_zero_contrast_ratios.csv')
 for j,m in enumerate(['M0']+PART):
  v=next(float(r['ratio']) for r in ratio if r['stage']=='Jurkat' and r['model']==m);d.scatter(v,j,color=COL[m],s=42,marker='D');d.text(v+.015,j,f'{v:.3f}',va='center',fontsize=9)
 d.axvline(1,color=INK,ls='--',lw=.7);d.set_xlim(.70,1.04);d.set_ylim(2.6,-.6);d.set_yticks(range(3),['M0','Rank one','Spectral']);d.set_xlabel('Projected error / projected truth energy')
 f.legend(handles=[Line2D([0],[0],color=color,marker=marker,lw=0,label={'observed':'Observed','M0':'M0','rank1_selected':'M0 + rank-one partner'}[m],mfc='white' if m=='M0' else color) for m,color,marker in mods],frameon=False,ncol=3,loc='lower center',bbox_to_anchor=(.55,.012),fontsize=8);save(f,'ED7')
def checks():
 off=read('Step2_Fig5_Official12.csv');assert sum(int(r['condition_seed_units']) for r in off)==540;assert sum(int(r['contrast_dimensions']) for r in off)==212
 signed=np.load(S/'RNA_signed_rectangle_predictions.npz');computed={m:float(np.mean((signed[m]-signed['truth_effect'])**2)) for m in ['M0','M1','R','scalar']};assert max(abs(computed[m]-CONTROL['C_MSE'][m]) for m in computed)<1e-14
 rs=[r for r in GENE if r['model']=='R'];vv=sorted([float(r['C_SSE_gain_vs_M0']) for r in rs if float(r['C_SSE_gain_vs_M0'])>0],reverse=True);share=sum(vv[:10])/sum(vv);assert abs(share-.921253)<1e-5
 example_checks={}
 corners=rna_corners()
 for e in js('biological_examples/rna/RNA_EXAMPLES.json')['example_records']:
  rr=[r for r in corners if r['case']==e['case']];assert [int(r['contrast_sign']) for r in rr]==[1,-1,-1,1]
  err=max(abs(sum(float(r[m+'_effect'])*int(r['contrast_sign']) for r in rr)-e[m+'_contrast']) for m in ['observed','M0','M1','R']);assert err<1e-12;example_checks[e['case']]={'four_corners':len(rr),'max_signed_contrast_difference':err}
 jr=read('biological_examples/jurkat/Jurkat_all_617_pairs_two_repeats.csv');assert len(jr)==1234
 mse0=sum((float(r['observed'])-float(r['M0']))**2 for r in jr);mse1=sum((float(r['observed'])-float(r['rank1_selected']))**2 for r in jr);jgain=1-mse1/mse0;assert abs(jgain-hval('Jurkat','high','rank1_selected'))<1e-12
 result={'new_experiment_or_model_fits':0,'RNA_four_corner_checks':example_checks,'Jurkat_1234_rows_recomputed_gain':jgain,'figures':['Fig'+str(i) for i in range(1,6)]+['ED'+str(i) for i in range(1,8)],'published_split_units':540,'published_projection_dimensions':212,'RNA_saved_contrast_MSE_recomputed':computed,'RNA_top_ten_share_of_positive_gain':share,'K562_high_rankone_gain':hval('K562','high','rank1_selected'),'Jurkat_primary_gain':hval('Jurkat','high','rank1_selected'),'uncertainty':'Neural seeds are optimization variability; folds/rectangles share genes. Jurkat reweighting ranges are conditional sensitivities, not calibrated confidence intervals.','source_files':{str(p.relative_to(S)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(S.rglob('*')) if p.is_file()}}
 (ROOT/'figure_code/FIGURE_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_files'},indent=2))
if __name__=='__main__':
 for fn in [fig1,fig2,fig3,fig4,fig5,ed1,ed2,ed3,ed4,ed5,ed6,ed7]:fn()
 checks()
