#!/usr/bin/env python3
"""Update Figures 1 and 4; add one complete-control supplementary figure."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pypdf import PdfReader,PdfWriter
R=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(R/'code'))
import make_figures as f
O=R/'partner_control/results';df=pd.read_csv(O/'FOUR_ARM_BY_CONDITION.csv')
colors={'original':f.GRAY,'whole_shrink':f.GOLD,'partner_only':f.TEAL}
labels={'original':'Original','whole_shrink':'Whole attenuation','partner_only':'Partner only'}

fig=plt.figure(figsize=(7.15,9.6))
ax=f.panel(fig,[.13,.705,.33,.205],'a','K562: expanded support',dy=.022)
for model,col in [('rank1_selected',f.TEAL),('spectral_selected',f.BLUE)]:
    for metric,marker,fill in [('sse','o',True),('projected_sse','s',False)]:
        vv=[f.growth_gain('K562',budget,model,metric) for budget in ['low','high']]
        shift=(-.045 if model=='rank1_selected' else .045)+(-.018 if fill else .018)
        ax.plot(np.arange(2)+shift,vv,ls='',marker=marker,ms=4.6,color=col,mfc=col if fill else 'white')
ax.axhline(0,color=f.GRAY,lw=.7);ax.set_xticks([0,1],['421-422','2,317-2,318']);ax.set_xlim(-.15,1.15);ax.set_ylim(-2,16);ax.set_xlabel('Training pairs per fold');ax.set_ylabel('Error reduction vs M0 (%)');f.style(ax,'y')
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=f.TEAL,label='Rank one'),Line2D([0],[0],marker='o',ls='',color=f.BLUE,label='Spectral'),Line2D([0],[0],marker='o',ls='',color=f.GRAY,label='Complete'),Line2D([0],[0],marker='s',ls='',color=f.GRAY,mfc='white',label='Partner')],loc='upper left',ncol=2,frameon=False,fontsize=6.8,columnspacing=.7,handlelength=.8)
ax=f.panel(fig,[.64,.705,.30,.205],'b','Jurkat: frozen-procedure gain',dy=.022)
for i,(metric,label) in enumerate([('sse','Complete'),('projected_sse','Partner')]):
    for j in [1,2]:ax.plot(f.growth_gain('Jurkat','high','rank1_selected',metric,j),i+(j-1.5)*.12,'o',mfc='white',mec=f.TEAL,ms=4)
    val=f.growth_gain('Jurkat','high','rank1_selected',metric)
    ax.plot(val,i,'D',color=f.TEAL,ms=5.5);ax.text(val,i+.20,f'{val:.2f}%',ha='center',fontsize=8,color=f.TEAL)
ax.set_yticks([0,1],['Complete','Partner']);ax.set_ylim(1.65,-.5);ax.set_xlim(0,11);ax.set_xlabel('Error reduction vs M0 (%)');f.style(ax)
ax.text(.02,.98,'617 held-out pairs\n2,464 training pairs',transform=ax.transAxes,va='top',fontsize=7.5)
fig.text(.64,.657,'Filled: pooled  |  Open: measured repeats',fontsize=6.8,color=f.GRAY)
styles=[('observed',f.INK,'D','Observed'),('M0',f.GRAY,'o','M0'),('rank1_selected',f.TEAL,'o','Rank one')]
ex=pd.read_csv(R/'revision/JURKAT_FOUR_CORNERS_V28.csv');ex=ex[ex.background.eq('ACTR10')]
ax=f.panel(fig,[.13,.405,.33,.175],'c','Four measured response corners',dy=.023)
for j,corner in enumerate(['AB','AD','CB','CD']):
    rr=ex[ex.corner.eq(corner)]
    for key,col,m,label in styles:ax.plot(j,float(rr[key].mean()),m,color=col,ms=4.7)
ax.set_xticks(range(4),['A+B','A+D','C+B','C+D']);ax.set_xlim(-.4,3.4);ax.set_ylabel('Mean growth response');f.style(ax,'y')
fig.text(.13,.345,'A = TTI2; C = ACTR10\nB = PGK1; D = OGFOD1',fontsize=7,color=f.GRAY,linespacing=1.4)
ax=f.panel(fig,[.64,.405,.30,.175],'d','Partner-exchange contrast',dy=.023)
cc=pd.read_csv(R/'revision/JURKAT_FOUR_CORNER_CONTRASTS_V28.csv');cc=cc[cc.background.eq('ACTR10')]
for j,(key,col,m,label) in enumerate(styles):
    val=float(cc[cc['repeat'].eq('mean')][key].iloc[0]);ax.plot(val,j,m,color=col,ms=5.4)
    for rep,offset in [('1',-.09),('2',.09)]:ax.plot(float(cc[cc['repeat'].eq(rep)][key].iloc[0]),j+offset,'o',mfc='white',mec=col,ms=3)
ax.axvline(0,color=f.GRAY,lw=.7,ls='--');ax.set_yticks(range(3),['Observed','M0','Rank one']);ax.set_ylim(2.6,-.6);ax.set_xlim(-.075,.008);ax.set_xticks([-.06,-.03,0]);ax.set_xlabel(r'$\Delta$ (growth units)');f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker=m,color=c,ls='',ms=4,label=l) for _,c,m,l in styles],bbox_to_anchor=(.10,.315),loc='lower left',ncol=3,frameon=False,fontsize=7,handlelength=.9,columnspacing=1)
for rect,letter,title,metric,limits in [([.13,.085,.33,.185],'e','Shen: retain complete-response benefit','full_gain_pp',(0,7.5)),([.64,.085,.30,.185],'f','Shen: adjust partner correction','partner_gain_pp',(-1.3,.3))]:
    ax=f.panel(fig,rect,letter,title,dy=.023)
    for j,pol in enumerate(['original','whole_shrink','partner_only']):
        vals=[float(df[(df.condition.eq(cell+'_high')) & df.policy.eq(pol)][metric].iloc[0]) for cell in ['HeLa','A549','293T']]
        ax.bar(np.arange(3)+(j-1)*.24,vals,width=.22,color=colors[pol],label=labels[pol])
    ax.axhline(0,color=f.GRAY,lw=.7);ax.set_xticks(range(3),['HeLa','A549','293T']);ax.set_ylim(*limits);ax.set_ylabel('Error reduction vs M0 (%)');f.style(ax,'y')
fig.legend(handles=[Rectangle((0,0),1,1,fc=colors[p],label=labels[p]) for p in colors],bbox_to_anchor=(.12,.015),loc='lower left',ncol=3,frameon=False,fontsize=7.2,columnspacing=1.3,handlelength=1)
f.save(fig,4,'Functional_Growth',['Horlbeck_figure_metrics.csv','JURKAT_FOUR_CORNERS_V28.csv','JURKAT_FOUR_CORNER_CONTRASTS_V28.csv','partner_control/results/FOUR_ARM_BY_CONDITION.csv'])
df.to_csv(R/'figure_source_data/Fig4_control_comparison.csv',index=False)

old=json.loads((R/'provenance/v32/FIGURE_AUDIT.json').read_text());done={a['figure'] for a in f.AUDIT};audit=sorted([a for a in old if a['figure'] not in done]+f.AUDIT,key=lambda a:a['figure'])
(R/'FIGURE_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
for kind,entries in [('Main',[a for a in audit if a['figure']<=5]),('Supplementary',[a for a in audit if a['figure']>5])]:
    w=PdfWriter()
    for a in entries:w.append(str(R/a['directory']/(a['stem']+'.pdf')))
    w.add_metadata({'/Title':'WitnessCell v33 '+kind+' figures','/Author':'WitnessCell'})
    with (R/('WitnessCell_v33_'+kind+'_Figures.pdf')).open('wb') as h:w.write(h)
print(json.dumps({a['figure']:a['text_outside_page'] for a in audit}))
