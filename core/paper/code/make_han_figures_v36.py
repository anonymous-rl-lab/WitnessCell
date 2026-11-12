#!/usr/bin/env python3
"""Retain Figure 4 a-f; replace Han panels and Supplementary Figure S5."""
from pathlib import Path
import json,sys,zipfile,os
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pypdf import PdfWriter
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'code'));import make_figures as f
H=R/'revision/han_calibration'
pooled=pd.read_csv(H/'POOLED_SCORES.csv');seeds=pd.read_csv(H/'SEED_SCORES.csv')
keys=pd.read_csv(H/'FULL_HAN_KEY_RESULTS.csv')
models=['HeLa__sf12','293T__sf12','A549__sf12','HeLa__v4','293T__v4','A549__v4']
readouts=[('full_raw','Reconstructed',f.TEAL,'D',-.14),('full_published','Published',f.BLUE,'o',.14)]

def panel(fig,rect,letter,title,dy=.022):
    ax=fig.add_axes(rect)
    fig.text(rect[0]-.035,rect[1]+rect[3]+dy,letter,fontweight='bold',fontsize=10,va='bottom')
    fig.text(rect[0],rect[1]+rect[3]+dy,title,fontweight='bold',fontsize=8.1,va='bottom')
    return ax
f.panel=panel

# Reuse the original six panels without changing their scientific inputs.
s=(R/'provenance/v33/make_figure4_v33.py').read_text().split('f.save(fig,4,')[0]
s=s[s.index('fig=plt.figure'):]
replacements={'figsize=(7.15,9.6)':'figsize=(7.15,10.8)','[.13,.705,.33,.205]':'[.13,.79,.33,.15]','[.64,.705,.30,.205]':'[.64,.79,.30,.15]', '.64,.657':'.64,.743','[.13,.405,.33,.175]':'[.13,.55,.33,.14]','[.64,.405,.30,.175]':'[.64,.55,.30,.14]', '.13,.345':'.13,.495','(.10,.315)':'(.61,.477)','[.13,.085,.33,.185]':'[.13,.31,.33,.13]','[.64,.085,.30,.185]':'[.64,.31,.30,.13]','(.12,.015)':'(.12,.25)','Shen: retain complete-response benefit':'Shen: complete-response gain','Shen: adjust partner correction':'Shen: partner-response gain'}
for a,b in replacements.items():s=s.replace(a,b)
df=pd.read_csv(R/'partner_control/results/FOUR_ARM_BY_CONDITION.csv')
colors={'original':f.GRAY,'whole_shrink':f.GOLD,'partner_only':f.TEAL}
labels={'original':'Original','whole_shrink':'Whole attenuation','partner_only':'Component control'}
exec(compile(s,'unchanged_figure4_a_to_f','exec'))
for rect,letter,title,column,limit in [([.13,.075,.33,.13],'g','Han K562: vs unchanged','partner_gain_vs_original_pct',(-.8,17)),([.64,.075,.30,.13],'h','Han K562: vs removed','partner_gain_vs_off_pct',(-.04,.80))]:
    ax=f.panel(fig,rect,letter,title)
    for i,model in enumerate(models[:2]):
        for panel_name,label,col,marker,offset in readouts:
            val=float(keys[(keys.panel==panel_name)&(keys.model==model)][column].iloc[0])
            ax.plot(val,i+offset,marker,color=col,mfc=col if marker=='D' else 'white',ms=4.3)
            ax.annotate(f'{val:.3f}' if val<1 else f'{val:.2f}',(val,i+offset),xytext=(4,0),textcoords='offset points',fontsize=6.7,color=col,va='center')
    ax.axvline(0,color=f.GRAY,lw=.7)
    ax.set_yticks([0,1],['HeLa source','293T source']);ax.set_ylim(1.5,-.5);ax.set_xlim(*limit)
    ax.set_xlabel('Partner error reduction (%)',fontsize=7.2);f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker='D',ls='',color=f.TEAL,label='Reconstructed'),Line2D([0],[0],marker='o',ls='',color=f.BLUE,mfc='white',label='Published')],bbox_to_anchor=(.115,.012),loc='lower left',frameon=False,ncol=2,fontsize=7,handlelength=.8,columnspacing=1.1)
fig.text(.59,.022,'SF12 sources; target: Han K562',fontsize=7.0,color=f.GRAY)
f.save(fig,4,'Functional_Growth',['evidence/source_data/publication/Horlbeck_figure_metrics.csv','revision/JURKAT_FOUR_CORNERS_V28.csv','revision/JURKAT_FOUR_CORNER_CONTRASTS_V28.csv','partner_control/results/FOUR_ARM_BY_CONDITION.csv','revision/han_calibration/POOLED_SCORES.csv'])
keys.to_csv(R/'figure_source_data/Fig4_Han_target_calibration.csv',index=False)

def gain(frame,panel_name,model,ref,metric):
    t=frame[(frame.panel==panel_name)&(frame.model==model)&(frame['repeat']=='pooled')].set_index('method')
    return 100*(1-float(t.loc['target_control',metric+'_sse'])/float(t.loc[ref,metric+'_sse']))

source_rows=[]
fig=plt.figure(figsize=(7.15,8.3))
for rect,letter,title,ref,metric in [([.18,.60,.29,.29],'a','Partner: vs unchanged','original','partner'),([.68,.60,.27,.29],'b','Partner: vs removed','partner_off','partner'),([.18,.16,.29,.29],'c','Complete: vs addition','single_additive','full')]:
    ax=f.panel(fig,rect,letter,title,dy=.026)
    for i,model in enumerate(models):
        for panel_name,label,col,marker,offset in readouts:
            val=gain(pooled,panel_name,model,ref,metric)
            values=[gain(sub,panel_name,model,ref,metric) for _,sub in seeds.groupby('seed')]
            ax.plot([min(values),max(values)],[i+offset]*2,color=col,lw=1.1,alpha=.75)
            ax.plot(val,i+offset,marker,color=col,mfc=col if marker=='D' else 'white',ms=4.0)
            source_rows.append(dict(panel=panel_name,model=model,comparator=ref,metric=metric,pooled_gain=val,partition_min=min(values),partition_max=max(values)))
    ax.axvline(0,color=f.GRAY,lw=.7);ax.axhline(2.5,color='#E0E5E8',lw=.7)
    ax.set_yticks(range(6),[m.replace('__',' ').replace('sf12','SF12').replace('v4','V4') for m in models]);ax.set_ylim(5.6,-.6)
    ax.set_xlabel('Error reduction (%)');f.style(ax)
    if letter=='a':ax.set_xlim(-.8,16)
    elif letter=='b':ax.set_xlim(-.06,.80)
    else:ax.set_xlim(-.04,.60)
pd.DataFrame(source_rows).to_csv(R/'figure_source_data/SuppFig5_Han_estimates_ranges.csv',index=False)

params=[]
han_snapshot=Path(os.environ['WITNESSCELL_HAN_SNAPSHOT'])
for path in sorted(han_snapshot.rglob('PARAMETERS.json')):
    name='han2017/'+path.relative_to(han_snapshot).as_posix()
    if name.endswith('/PARAMETERS.json') and '/results/full_' in name:
        parts=name.split('/');idx=parts.index('results');pn,seed,fold=parts[idx+1:idx+4]
        state=json.loads(path.read_text());order=['293T__sf12','293T__v4','A549__sf12','A549__v4','HeLa__sf12','HeLa__v4']
        # Frozen array order is HeLa, A549, 293T, each SF12 then V4.
        order=['HeLa__sf12','HeLa__v4','A549__sf12','A549__v4','293T__sf12','293T__v4']
        for m,alpha in zip(order,state['alpha']):params.append(dict(panel=pn,seed=int(seed),fold=int(fold),model=m,alpha=float(alpha)))
pa=pd.DataFrame(params);assert len(pa)==180
pa.to_csv(R/'figure_source_data/Han_target_strengths.csv',index=False)
ax=f.panel(fig,[.68,.16,.27,.29],'d','Estimated partner strength',dy=.026)
for i,model in enumerate(models):
    for pn,label,col,marker,offset in readouts:
        vals=pa[(pa.panel==pn)&(pa.model==model)].alpha.to_numpy();assert len(vals)==15
        yy=i+offset+np.linspace(-.055,.055,len(vals))
        ax.scatter(vals,yy,s=8,color=col,alpha=.45,edgecolors='none')
        ax.plot(np.median(vals),i+offset,marker,color=col,mfc=col if marker=='D' else 'white',ms=4.5)
ax.set_yticks(range(6),[m.replace('__',' ').replace('sf12','SF12').replace('v4','V4') for m in models]);ax.set_ylim(5.6,-.6);ax.set_xlim(-.05,1.05);ax.set_xlabel(r'Target coefficient $\widehat{\alpha}$');ax.axhline(2.5,color='#E0E5E8',lw=.7);f.style(ax)
fig.legend(handles=[Line2D([0],[0],marker='D',ls='',color=f.TEAL,label='Reconstructed'),Line2D([0],[0],marker='o',ls='',color=f.BLUE,mfc='white',label='Published')],bbox_to_anchor=(.14,.045),loc='lower left',frameon=False,ncol=2,fontsize=7.2,handlelength=.8)
fig.text(.18,.025,'Target: Han K562. Labels identify source models; segments show five-partition ranges.',fontsize=6.8,color=f.GRAY)
f.save(fig,10,'Han_Target_Calibration',['revision/han_calibration/POOLED_SCORES.csv','revision/han_calibration/SEED_SCORES.csv','reproduction/WitnessCell_Han2017_Frozen_Reassessment_v3.zip'])

old=json.loads((R/'provenance/v35/FIGURE_AUDIT.json').read_text()); f.AUDIT=json.loads((R/'revision/FIG1_AUDIT_V36.json').read_text())+f.AUDIT;changed={a['figure'] for a in f.AUDIT}
audit=sorted([a for a in old if a['figure'] not in changed]+f.AUDIT,key=lambda a:a['figure'])
(R/'FIGURE_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
for kind,entries in [('Main',[a for a in audit if a['figure']<=5]),('Supplementary',[a for a in audit if a['figure']>5])]:
    w=PdfWriter()
    for a in entries:w.append(str(R/a['directory']/(a['stem']+'.pdf')))
    w.add_metadata({'/Title':'WitnessCell v36 '+kind+' figures','/Author':'WitnessCell'})
    with (R/('WitnessCell_v36_'+kind+'_Figures.pdf')).open('wb') as h:w.write(h)
assert all(not x['text_outside_page'] for x in audit),audit
print(json.dumps({a['figure']:a['text_outside_page'] for a in audit}))
