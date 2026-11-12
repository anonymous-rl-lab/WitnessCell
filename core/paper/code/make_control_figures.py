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

fig=plt.figure(figsize=(7.15,8.5))
ax=f.panel(fig,[.07,.68,.88,.24],'a','Predict responses from measured interventions',dy=.023);f.off(ax)
f.box(ax,.02,.72,.40,.20,'Controls + single responses\nMeasured or inferred endpoint effects',fc='#EAF2F8',size=8)
f.box(ax,.58,.72,.40,.20,'Measured combinations\nPermitted training support',fc='#E5F2F0',size=8)
f.box(ax,.02,.31,.40,.23,'Nonlinear composition\n'+r'$B=\mu_0+S(e_g+e_h)$',fc='#EAF2F8',size=9)
f.box(ax,.58,.31,.40,.23,'Endpoint residual transfer\n'+r'$\delta_Q=X_Q\Theta$',fc='#E5F2F0',size=9)
f.arrow(ax,(.22,.70),(.22,.56),f.BLUE);f.arrow(ax,(.78,.70),(.78,.56),f.TEAL)
f.arrow(ax,(.43,.43),(.56,.43));ax.text(.49,.52,'residuals',ha='center',fontsize=6.8,color=f.GRAY)
f.arrow(ax,(.22,.29),(.43,.13),f.BLUE);f.arrow(ax,(.78,.29),(.57,.13),f.TEAL)
ax.text(.50,.05,r'$\widehat Y_{\rm core}=B+\gamma\delta$'+'  |  Released benchmark configuration',ha='center',fontsize=9)
ax=f.panel(fig,[.07,.40,.40,.21],'b','Remove endpoint contributions',dy=.023);f.off(ax)
ax.text(.5,.89,'A+B    A+D    C+B    C+D',ha='center',fontsize=8.8)
ax.text(.5,.64,r'$\Delta=Y_{AB}-Y_{AD}-Y_{CB}+Y_{CD}$',ha='center',fontsize=9)
ax.text(.5,.39,r'$P_Q=I-X_QX_Q^+$',ha='center',fontsize=11,color=f.TEAL)
ax.text(.5,.10,'Partner contrasts remain after\nseparate endpoint contributions cancel',ha='center',fontsize=8,linespacing=1.4)
ax=f.panel(fig,[.55,.40,.40,.21],'c','Learn additional partner responses',dy=.023);f.off(ax)
f.box(ax,.03,.56,.94,.33,'RNA: projected pair kernel\nTraining-half readout weights',fc='#EAF2F8',size=8.5)
f.box(ax,.03,.08,.94,.33,'Function: shared residual matrix\nTraining-selected rank-one correction',fc='#E5F2F0',size=8.3)
ax=f.panel(fig,[.07,.225,.88,.105],'d','Test the response component changed by each branch',dy=.021);f.off(ax)
for x,s in [(.01,'Added branch'),(.41,'Complete response'),(.71,'Partner contrast')]:ax.text(x,.87,s,fontweight='bold',fontsize=8)
ax.axhline(.72,color='#D8E2E8',lw=.8)
for y,l,p in [(.48,'Endpoint residual',r'Unchanged: $P_Q\delta_Q=0$'),(.12,'Explicit partner learner','Can change; evaluate its gain')]:
    ax.text(.01,y,l,fontsize=8);ax.text(.41,y,'Can change',fontsize=8);ax.text(.71,y,p,fontsize=7.4,color=f.TEAL)
ax=f.panel(fig,[.07,.035,.88,.115],'e','Adjust the partner correction while preserving endpoints',dy=.023);f.off(ax)
ax.text(.5,.79,r'$\widehat Y_\alpha=M_0+(I-P_Q)D+\alpha P_QD$',ha='center',fontsize=12,color=f.INK)
ax.text(.5,.46,'Retain endpoint part  |  Adjust partner part',ha='center',fontsize=8.2,color=f.TEAL)
ax.text(.5,.23,'Internal held-combination predictions set '+r'$\alpha\in[0,1]$',ha='center',fontsize=8.0)
ax.text(.5,.04,'Specified query panel  |  Query identities used; query outcomes excluded',ha='center',fontsize=7.6,color=f.GRAY)
f.save(fig,1,'Method',['Main Methods and Supplementary Note 17','partner_control/METHOD_DERIVATION.md'])

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
for rect,letter,title,metric,limits in [([.13,.085,.33,.185],'e','Shen: preserve complete prediction','full_gain_pp',(0,7.5)),([.64,.085,.30,.185],'f','Shen: adjust partner correction','partner_gain_pp',(-1.3,.3))]:
    ax=f.panel(fig,rect,letter,title,dy=.023)
    for j,pol in enumerate(['original','whole_shrink','partner_only']):
        vals=[float(df[(df.condition.eq(cell+'_high')) & df.policy.eq(pol)][metric].iloc[0]) for cell in ['HeLa','A549','293T']]
        ax.bar(np.arange(3)+(j-1)*.24,vals,width=.22,color=colors[pol],label=labels[pol])
    ax.axhline(0,color=f.GRAY,lw=.7);ax.set_xticks(range(3),['HeLa','A549','293T']);ax.set_ylim(*limits);ax.set_ylabel('Error reduction vs M0 (%)');f.style(ax,'y')
fig.legend(handles=[Rectangle((0,0),1,1,fc=colors[p],label=labels[p]) for p in colors],bbox_to_anchor=(.12,.015),loc='lower left',ncol=3,frameon=False,fontsize=7.2,columnspacing=1.3,handlelength=1)
f.save(fig,4,'Functional_Growth',['Horlbeck_figure_metrics.csv','JURKAT_FOUR_CORNERS_V28.csv','JURKAT_FOUR_CORNER_CONTRASTS_V28.csv','partner_control/results/FOUR_ARM_BY_CONDITION.csv'])
df.to_csv(R/'figure_source_data/Fig4_control_comparison.csv',index=False)

# Supplementary figure: all nonzero cases, all repeat coefficients, and panel effects.
fig=plt.figure(figsize=(7.15,8.8));allscore=pd.read_csv(O/'FOUR_ARM_SCORES.csv');dec=pd.read_csv(O/'DECISIONS.csv');rep=pd.read_csv(O/'REPEAT_DECISIONS.csv');pan=pd.read_csv(O/'QUERY_PANEL_BY_CONDITION.csv')
wide=allscore.pivot(index='feature_id',columns='policy',values='partner_gain_pp');active=dec[~dec.selected_zero].feature_id
ax=f.panel(fig,[.13,.59,.34,.30],'a','All active functional corrections',dy=.025)
ax.plot([-12,22],[-12,22],ls='--',lw=.8,color=f.GRAY)
ax.scatter(wide.loc[active,'original'],wide.loc[active,'partner_only'],s=26,color=f.TEAL,edgecolor='white',linewidth=.5)
ax.set_xlim(-12,22);ax.set_ylim(-12,22);ax.set_xlabel('Original partner gain (%)');ax.set_ylabel('Controlled partner gain (%)');f.style(ax,'both')
ax.text(.03,.96,'27 nonzero corrections\n28 zero-selected fits unchanged',transform=ax.transAxes,va='top',fontsize=7)
ax=f.panel(fig,[.64,.59,.30,.30],'b','Training-repeat agreement',dy=.025)
rw=rep.pivot(index='feature_id',columns='training_repeat',values='alpha')
ax.plot([0,1],[0,1],ls='--',lw=.8,color=f.GRAY);ax.scatter(rw.loc[active,1],rw.loc[active,2],s=25,color=f.BLUE,edgecolor='white',lw=.5)
ax.set_xlim(-.05,1.05);ax.set_ylim(-.05,1.05);ax.set_xlabel('Coefficient from repeat 1');ax.set_ylabel('Coefficient from repeat 2');f.style(ax,'both')
ax=f.panel(fig,[.13,.12,.81,.30],'c','Query grouping changes which directions can be adjusted',dy=.025)
conds=['Jurkat_high','K562_high','K562_low','K562_layout','HeLa_high','A549_high','293T_high'];rules=['full','identity_halves','identity_thirds','singletons'];cols=[f.TEAL,f.BLUE,f.GOLD,f.GRAY]
for j,(rule,col) in enumerate(zip(rules,cols)):
    vals=[float(pan[pan.condition.eq(c)&pan.panel_rule.eq(rule)].partner_gain_pp.iloc[0]) for c in conds]
    ax.plot(np.arange(7)+(j-1.5)*.16,vals,ls='',marker='o',ms=4.3,color=col,label=rule.replace('identity_','').replace('full','Full panel').replace('singletons','Singletons'))
ax.axhline(0,color=f.GRAY,lw=.7);ax.set_xticks(range(7),['Jurkat','K562\nhigh','K562\nlow','K562\nlayouts','HeLa','A549','293T']);ax.set_ylabel('Partner-error reduction (%)');ax.set_ylim(-1.5,13);f.style(ax,'y');ax.legend(frameon=False,ncol=4,fontsize=7,loc='upper right')
f.save(fig,9,'Component_Control',['partner_control/results/FOUR_ARM_SCORES.csv','partner_control/results/REPEAT_DECISIONS.csv','partner_control/results/QUERY_PANEL_BY_CONDITION.csv'])

old=json.loads((R/'provenance/v31/FIGURE_AUDIT.json').read_text());done={a['figure'] for a in f.AUDIT};audit=sorted([a for a in old if a['figure'] not in done]+f.AUDIT,key=lambda a:a['figure'])
(R/'FIGURE_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
for kind,entries in [('Main',[a for a in audit if a['figure']<=5]),('Supplementary',[a for a in audit if a['figure']>5])]:
    w=PdfWriter()
    for a in entries:w.append(str(R/a['directory']/(a['stem']+'.pdf')))
    w.add_metadata({'/Title':'WitnessCell v32 '+kind+' figures','/Author':'WitnessCell'})
    with (R/('WitnessCell_v32_'+kind+'_Figures.pdf')).open('wb') as h:w.write(h)
print(json.dumps({a['figure']:a['text_outside_page'] for a in audit}))
