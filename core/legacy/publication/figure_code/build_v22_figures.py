#!/usr/bin/env python3
"""Rebuild v22 Figures 1 and 7 from package-local, frozen summaries. No model fitting."""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'source_data/v22_synthesis'
OUT=ROOT/'figures'; OUT.mkdir(exist_ok=True)
INK='#243544'; BLUE='#2476A5'; TEAL='#178579'; ORANGE='#C68A26'; RED='#B54A50'; GREY='#85929E'; GRID='#DDE4E9'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.edgecolor':GRID,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})

def save(fig,name):
    for ext in ('pdf','png','svg'):
        fig.savefig(OUT/f'{name}.{ext}',dpi=220,bbox_inches=None)
    plt.close(fig)

def title(ax,letter,text):
    pos=ax.get_position()
    ax.figure.text(max(.02,pos.x0-.032),pos.y1+.011,letter,fontsize=12,fontweight='bold',va='bottom')
    ax.text(0,1.065,text,transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom')

def clean(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='x',color=GRID,lw=.5,zorder=0)
    ax.set_axisbelow(True)

def box(ax,x,y,w,h,heading,body,fill='#F2F5F7',edge=GRID):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.009,rounding_size=.015',facecolor=fill,edgecolor=edge,lw=.9))
    ax.text(x+w/2,y+h*.74,heading,ha='center',va='center',weight='bold',fontsize=8.6)
    ax.text(x+w/2,y+h*.34,body,ha='center',va='center',fontsize=7.6,linespacing=1.45)

def arrow(ax,a,b):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=10,color=INK,lw=1))

def fig1():
    f=plt.figure(figsize=(8.4,9.3))
    a=f.add_axes([.075,.635,.88,.29]);a.set(xlim=(0,1),ylim=(0,1));a.axis('off')
    title(a,'a','Separate measurement calibration, endpoint refinement and partner learning')
    box(a,.01,.55,.27,.35,'Allowed observations','Single perturbations and controls\nMeasured training combinations\nFixed functional annotations')
    box(a,.37,.64,.26,.26,'Measurement reference','Common target composition\nOptional training-derived offsets',fill='#EAF5F2',edge=TEAL)
    box(a,.72,.55,.27,.35,'Complete response core','Measured or inferred endpoints\nNonlinear composition\nRefinement from training combinations')
    box(a,.37,.12,.26,.28,'Explicit partner component','Task-specific representation\nTraining-derived shrinkage',fill='#EAF1F7',edge=BLUE)
    box(a,.72,.12,.27,.28,'Delivered response','Complete core + partner increment\nSpecified response scale',fill='#FBF3E8',edge=ORANGE)
    arrow(a,(.29,.76),(.36,.76));arrow(a,(.64,.76),(.71,.76))
    arrow(a,(.145,.54),(.36,.28));arrow(a,(.50,.63),(.50,.41))
    arrow(a,(.64,.26),(.71,.26));arrow(a,(.855,.54),(.855,.41))
    a.text(.01,-.015,'Every fitted transformation and selection uses the permitted training observations.',fontsize=7.8)

    b=f.add_axes([.075,.325,.88,.245]);b.axis('off');title(b,'b','Shared construction; distinct fitted models and evaluation tasks')
    rows=[
      ['Released RNA core','Endpoint completion + composition + correction','Four-study benchmark; no explicit partner term'],
      ['Matched RNA','Complete M0 + response-profile kernel','55 combinations; readout weighting tested'],
      ['Growth','Repeat-specific M0 + signed pair matrix','K562 budgets; frozen-rule Jurkat evaluation'],
      ['Program extension','Fixed program + query-invariant reference','RNA development; frozen drug-combination test'],
      ['H1 calibration','Common-reference input + complete M0','Historical development; 21 held-out conditions'],
      ['H1 shared features\n(SF12)','Calibrated M0 + 12 shared coefficients','Historical development; five overlapping tasks']]
    table=b.table(cellText=rows,colLabels=['Implementation','Construction','Evidence scope'],cellLoc='left',colLoc='left',bbox=[0,-.035,1,1.02],colWidths=[.20,.40,.40])
    table.auto_set_font_size(False);table.set_fontsize(7.1)
    for (i,j),cell in table.get_celld().items():
      cell.set_edgecolor(GRID);cell.set_linewidth(.5)
      if i==0:cell.set_facecolor('#EAF0F5');cell.set_text_props(weight='bold')
      else:cell.set_facecolor('white' if i%2 else '#FAFBFC')

    c=f.add_axes([.075,.07,.41,.18]);c.axis('off');title(c,'c','Two complementary prediction questions')
    c.text(0,.88,r'$\widehat{Y}=F_{\mathrm{core}}+\Delta_{\mathrm{partner}}$',fontsize=13)
    c.text(0,.58,'Delivered accuracy',weight='bold',fontsize=8.5)
    c.text(0,.40,r'Complete-response loss: $\|\widehat{Y}-Y\|^2$',fontsize=8.8)
    c.text(0,.14,'Additional partner prediction',weight='bold',fontsize=8.5)
    c.text(0,-.04,r'Endpoint-eliminated loss: $\|P_Q(\widehat{Y}-Y)\|^2$',fontsize=8.8)
    c.text(0,-.24,r'$P_Q=I-X_QX_Q^+$ is an evaluation operator.',fontsize=7.8)

    d=f.add_axes([.575,.07,.38,.18]);d.axis('off');title(d,'d','Experimental support is part of the model')
    yy=[.85,.55,.25,-.05]
    heads=['Measurement reference','Repeated observations','Combination coverage','Frozen evaluation']
    tails=['Comparability of response targets','Readout-specific adjustment','Constraints on new-partner predictions','Transfer and limits beyond development']
    for y,h,t in zip(yy,heads,tails):
      d.text(0,y,h,fontsize=8.2,weight='bold')
      d.text(0,y-.13,t,fontsize=7.6)
    save(f,'Fig1')

def fig7():
    scores=json.loads((SRC/'qualified_rectangles_scores.json').read_text())
    cal=list(csv.DictReader((SRC/'calibration_main_summary.csv').open()))
    f=plt.figure(figsize=(8.5,9.8))
    a=f.add_axes([.21,.735,.28,.18]); title(a,'a','Calibration improves the complete M0')
    labels=['Raw mixture','Direct 7-library mean','7-library calibration','9-library calibration']
    vals=[float(r['MSE'])*1000 for r in cal]
    gains=[float(r['gain_vs_raw_pct']) for r in cal]
    for i,(v,g) in enumerate(zip(vals,gains)):
      a.barh(i,v,height=.55,color=[GREY,'#A5C2D2',TEAL,BLUE][i],zorder=3)
      a.text(v+.025,i,'reference' if i==0 else f'{g:.2f}%',ha='left',va='center',fontsize=7.7)
    a.set(yticks=range(4),yticklabels=labels,xlim=(0,1.98),xticks=[0,.5,1,1.5],xlabel=r'RNA MSE ($\times 10^{-3}$)');a.invert_yaxis();clean(a)
    a.text(-.50,-.45,'One fixed seven-library target; 21 held-out conditions.\nMatched seven-library comparison: 2.73% further reduction.',transform=a.transAxes,fontsize=7.5,linespacing=1.45)

    b=f.add_axes([.66,.735,.29,.18]);title(b,'b','Partner gains vary across SF12 tasks')
    tasks=list(scores['tasks']) if 'tasks' in scores else list(scores['per_task'])
    s=scores.get('tasks',scores.get('per_task'))
    for i,q in enumerate(tasks):
      m=s[q]['methods']['old12'];v1=m['full_gain_vs_M0_pct'];v2=m['partner_gain_vs_M0_pct']
      b.plot([v1,v2],[i,i],color=GRID,lw=1.4)
      b.scatter(v1,i,color=BLUE,s=31,marker='o',zorder=4)
      b.scatter(v2,i,color=ORANGE,s=29,marker='s',zorder=4)
    b.axvline(0,color=GREY,lw=.8);b.set(yticks=range(5),yticklabels=tasks,xlim=(-7,8),xticks=[-5,0,5],xlabel='Error reduction relative to M0 (%)');b.invert_yaxis();clean(b)
    b.legend(handles=[Line2D([],[],marker='o',color=BLUE,lw=0,label='Complete RNA'),Line2D([],[],marker='s',color=ORANGE,lw=0,label='Partner contrast')],frameon=False,fontsize=7.3,loc='upper left',bbox_to_anchor=(-.02,-.26),ncol=1)

    c=f.add_axes([.13,.39,.80,.20]);title(c,'c','Five tasks share eight combinations and span three contrast directions')
    edges=sorted({p for q in tasks for p in s[q]['query_pairs']})
    mat=np.zeros((5,8));signs=[.5,-.5,-.5,.5]
    for i,q in enumerate(tasks):
      for j,p in enumerate(s[q]['query_pairs']):mat[i,edges.index(p)]=signs[j]
    c.imshow(mat,cmap=ListedColormap(['#B8D6E6','#FAFBFC','#E6C68E']),vmin=-.5,vmax=.5,aspect='auto')
    for i in range(5):
      for j in range(8):
        if mat[i,j]:c.text(j,i,'+' if mat[i,j]>0 else '−',ha='center',va='center',fontsize=12,color=INK)
    c.set(yticks=range(5),yticklabels=tasks,xticks=range(8),xticklabels=[p.replace('+','\n') for p in edges]);c.tick_params(axis='both',length=0,pad=6)
    c.set_xticks(np.arange(-.5,8,1),minor=True);c.set_yticks(np.arange(-.5,5,1),minor=True);c.grid(which='minor',color='white',lw=1.5);c.tick_params(which='minor',length=0)
    c.text(0,-.40,'Each nonzero coefficient has magnitude ½. Q4 was the previously examined positive block.',transform=c.transAxes,fontsize=7.6)

    d=f.add_axes([.13,.105,.33,.15]);title(d,'d','Partner gain is sensitive to task removal')
    lod=scores['leave_one_task_out_partner_gain_vs_M0_pct']
    for i,q in enumerate(tasks):
      v=lod[q]['old12'];d.scatter(v,i,color=ORANGE,s=29,zorder=3)
      d.text(v+.09,i,f'{v:.2f}',va='center',fontsize=7.5)
    d.axvline(0,color=GREY,lw=.8);d.set(yticks=range(5),yticklabels=['Omit '+q for q in tasks],xlim=(-.7,2.25),xticks=[0,1,2],xlabel='Pooled partner error reduction (%)');d.invert_yaxis();clean(d)

    e=f.add_axes([.58,.095,.37,.16]);e.axis('off');title(e,'e','Aggregate outcomes retain model identity')
    a5=scores['aggregate']['all5']
    rows=[]
    for name,key in [('SF12','old12'),('Case shrinkage','case_R')]:
      r=a5[key];rows.append([name,f"{r['edge_balanced_full_gain_vs_M0_pct']:+.2f}%",f"{r['partner_gain_vs_M0_pct']:+.2f}%"])
    tab=e.table(cellText=rows,colLabels=['Model','Complete RNA','Partner'],cellLoc='center',colLoc='center',bbox=[-.02,.30,1.03,.70],colWidths=[.36,.36,.28]);tab.auto_set_font_size(False);tab.set_fontsize(7.4)
    for (i,j),cell in tab.get_celld().items():
      cell.set_edgecolor(GRID);cell.set_linewidth(.6)
      if i==0:cell.set_facecolor('#EAF0F5');cell.set_text_props(weight='bold')
      if i==2 and j==1:cell.set_text_props(color=RED)
    e.text(-.02,.16,'Complete RNA: edge-balanced losses.\nPartner: pooled four-corner squared errors.\nHistorical development; no independent replication.',fontsize=7.4,va='top',linespacing=1.5)
    save(f,'Fig7')

if __name__=='__main__':
    fig1();fig7()
