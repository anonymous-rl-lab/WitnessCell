#!/usr/bin/env python3
"""Build the v22 full-system figures from package-local, versioned tables.

Historical predictions retain their original configuration. New experiment
tables are included only when available; no predictions are fitted here.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[k]='1'
import csv,json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import build_historical_panels as H
import build_execution_panels as E

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'source_data'
OUT=ROOT/'figures'
COL=H.COL; INK=H.INK; GRID=H.GRID
save=H.save; panel=H.panel; box=H.box; arrow=H.arrow

def rows(name):
    return list(csv.DictReader((S/name).open()))

def fig1():
    f=plt.figure(figsize=(7.7,8.4))
    a=f.add_axes([.07,.64,.88,.29]);a.axis('off')
    panel(a,'a','WitnessCell: from allowed observations to a complete prediction',None)
    box(a,.01,.54,.23,.34,'Measured singles\nTraining combinations\nGene annotations',size=8.5)
    box(a,.31,.63,.29,.25,'Endpoint responses\nMeasured or training-inferred',size=8.3)
    box(a,.31,.15,.29,.31,'Partner component\nResponse kernel or signed matrix',color='#EAF0F6',edge=COL['M1'],size=8.3)
    box(a,.70,.58,.28,.30,'Core prediction\nComposition + endpoint correction',size=8.2)
    box(a,.70,.10,.28,.31,'Final response\nCore + adjusted partner term',color='#FAF2E5',edge=COL['R'],size=8.3)
    arrow(a,.25,.73,.30,.75);arrow(a,.25,.64,.30,.32)
    arrow(a,.61,.75,.69,.74);arrow(a,.61,.30,.69,.25);arrow(a,.84,.56,.84,.43)
    a.text(.02,.05,'All fitted quantities use only the permitted training information.',transform=a.transAxes,fontsize=8.1)
    b=f.add_axes([.08,.37,.86,.18]);b.axis('off');panel(b,'b','One system definition; task-specific fitted implementations',None)
    text=[['Configuration','Core','Partner / adjustment'],
          ['Original RNA','Released v14 core','Original v14: G = 0\nSix-split extension: G and adjustment'],
          ['Measured-endpoint RNA','Complete nonlinear M0','Matched 55-pair and new-source tasks\nProjected kernel + readout adjustment'],
          ['Growth','Repeat-specific nonlinear M0','Shared signed symmetric matrix\nUnit weight']]
    t=b.table(cellText=text[1:],colLabels=text[0],cellLoc='left',colLoc='left',bbox=[0,-.10,1,1.03],colWidths=[.23,.38,.39]);t.auto_set_font_size(False);t.set_fontsize(8)
    for (i,j),c in t.get_celld().items():
        c.set_edgecolor(GRID);c.set_linewidth(.5)
        if i==0:c.set_facecolor('#EAF0F6');c.set_text_props(weight='bold')
    c=f.add_axes([.09,.07,.38,.18]);c.axis('off');panel(c,'c','Complete response and information attribution',None)
    c.text(0,.68,r'$\widehat{Y}=F_{\rm core}+G\,\mathrm{diag}(w)$',fontsize=12,transform=c.transAxes)
    c.text(0,.38,'Predict the full response first.',transform=c.transAxes,fontsize=8.5)
    c.text(0,.13,r'Then test $P_Q(\widehat{Y}-Y)$, $P_Q=I-X_QX_Q^+$.',transform=c.transAxes,fontsize=9)
    c.text(0,-.14,'Projection is an evaluation diagnostic.',transform=c.transAxes,fontsize=8)
    d=f.add_axes([.60,.07,.35,.18]);d.axis('off');panel(d,'d','Evidence preserved with its model identity',None)
    d.text(0,.75,'Released core and exact six-split integration',transform=d.transAxes,fontsize=8.4)
    d.text(0,.47,'Matched RNA: additional partner prediction',transform=d.transAxes,fontsize=8.4)
    d.text(0,.19,'K562: effect of training observations',transform=d.transAxes,fontsize=8.4)
    d.text(0,-.09,'Frozen evaluations: Jurkat and Dixit RNA',transform=d.transAxes,fontsize=8.4)
    save(f,'Fig1')

def fig2():
    f=plt.figure(figsize=(7.7,8.7))
    a=f.add_axes([.16,.60,.35,.30]);panel(a,'a','Released benchmark: top-100 panel','x')
    rr=rows('released_benchmark/formal_top100_leaderboard.csv')[:8]
    for i,r in enumerate(rr):
        v=float(r['mean_operation_rank']);a.barh(i,v,color=COL['M1'] if r['method']=='WitnessCell' else '#ACB6BE',height=.57)
        a.text(v+.10,i,f'{v:.3f}',va='center',fontsize=8)
    a.set_yticks(range(len(rr)),[r['method'] for r in rr]);a.invert_yaxis();a.set_xlim(0,9);a.set_xlabel('Mean operation rank (lower is better)')
    b=f.add_axes([.72,.60,.23,.30]);panel(b,'b','Metric-specific positions','x')
    mm=[r for r in rows('released_benchmark/official_top100_per_metric_leaderboard.csv') if r['method']=='WitnessCell']
    labels={'PCC-delta':'PCC','MSE':'MSE','E-distance':'Energy distance','Wasserstein':'Wasserstein','KL':'Symmetric KL','common_DEG':'Common DEG'}
    for i,r in enumerate(mm):
        v=float(r['position']);b.scatter(v,i,s=30,color=COL['M1']);b.text(v+.35,i,str(int(v)),va='center',fontsize=8)
    b.set_yticks(range(len(mm)),[labels.get(r['metric'],r['metric']) for r in mm]);b.invert_yaxis();b.set_xlim(.3,10.5);b.set_xticks([1,5,9]);b.set_xlabel('Position among 16 entrants')
    c=f.add_axes([.16,.28,.36,.19]);panel(c,'c','Non-target full-response comparisons','x')
    rr=[r for r in rows('Fig2_prediction_ratios.csv') if r['metric']=='loss_mse']
    for i,r in enumerate(rr):
        v=100*float(r['gain']);lo=100*float(r['gain_ci_low']);hi=100*float(r['gain_ci_high'])
        c.plot([lo,hi],[i,i],color=COL['M1'],lw=1.3);c.scatter(v,i,color=COL['M1'],s=29)
    c.axvline(0,lw=.7,color=INK);c.set_yticks(range(len(rr)),[r['dataset']+' / '+('additive' if r['comparator']=='shared_endpoint_additive' else 'saturation') for r in rr]);c.invert_yaxis();c.set_xlim(-3,46);c.set_xlabel('MSE reduction (%)')
    d=f.add_axes([.72,.28,.23,.19]);panel(d,'d','Weighted non-target error','x')
    rr=[r for r in rows('Fig2_prediction_ratios.csv') if r['metric']=='loss_wmse']
    for i,r in enumerate(rr):
        v=100*float(r['gain']);lo=100*float(r['gain_ci_low']);hi=100*float(r['gain_ci_high'])
        d.plot([lo,hi],[i,i],color=COL['R'],lw=1.3);d.scatter(v,i,color=COL['R'],s=29)
    d.axvline(0,lw=.7,color=INK);d.set_yticks(range(len(rr)),['Norman / add.','Norman / sat.','Wessels / add.','Wessels / sat.']);d.invert_yaxis();d.set_xlim(-8,60);d.set_xlabel('WMSE reduction (%)')
    e=f.add_axes([.09,.06,.86,.09]);e.axis('off');panel(e,'e','Benchmark scope and combination-only result',None)
    e.text(0,.65,'654 condition-seed records = 540 combinations + 114 unseen singles',transform=e.transAxes,fontsize=9)
    e.text(0,.05,'Combination-only WMSE ratio: 0.7513 [0.7031, 0.8019] versus measured-single additivity with mean imputation.',transform=e.transAxes,fontsize=8.25)
    e.text(0,-.47,'This historical WMSE includes targeted readouts. Panels c-d exclude them and retain separate denominators.',transform=e.transAxes,fontsize=8)
    save(f,'Fig2')

def fig3():
    f=plt.figure(figsize=(7.7,7.6))
    a=f.add_axes([.09,.65,.86,.22]);a.axis('off');panel(a,'a','What combination-trained endpoint correction can represent',None)
    a.text(.02,.68,r'$\delta_Q=\gamma K_{QT}(K_{TT}+\eta I)^{-1}R_T=X_Q\Theta_T$',fontsize=15,transform=a.transAxes)
    a.text(.02,.30,r'$P_QX_Q=0\quad\Longrightarrow\quad P_Q\delta_Q=0$',fontsize=14,transform=a.transAxes)
    a.text(.02,-.02,'Endpoint estimation can improve full responses without changing endpoint-eliminated predictions.',fontsize=9,transform=a.transAxes)
    b=f.add_axes([.20,.37,.72,.18]);panel(b,'b','Complete and projected gains against the paired backbone','x')
    rr=rows('Step2_Fig5_Official12.csv')
    for i,r in enumerate(rr):
        v=100*float(r['full_response_gain_vs_joint_backbone']);b.barh(i,v,height=.48,color=COL['M1']);b.scatter(100*float(r['projected_error_gain_vs_joint_backbone']),i,color=COL['R'],marker='D',s=26,zorder=5);b.text(v+.6,i,f'{v:.2f}',va='center',fontsize=8)
    b.set_yticks(range(4),[r['dataset'].replace('_exp6',' exp6') for r in rr]);b.invert_yaxis();b.set_xlim(-1,39);b.set_xlabel('Error reduction (%)')
    b.legend(handles=[Line2D([0],[0],color=COL['M1'],lw=5,label='Complete'),Line2D([0],[0],color=COL['R'],marker='D',lw=0,label='Projected')],frameon=False,ncol=2,loc='upper right',fontsize=8)
    c=f.add_axes([.12,.09,.82,.15]);panel(c,'c','All twelve released splits retain the same projected predictions',None)
    rr=rows('Official12_MAIN_TABLE_all_genes_raw.csv');x=np.arange(12)
    c.plot(x,[float(r['v14_joint_saturation_backbone']) for r in rr],marker='o',mfc='white',ms=6,color=COL['M0'],lw=1,label='Paired backbone')
    c.scatter(x,[float(r['WitnessCell_v14']) for r in rr],marker='D',s=18,color=COL['M1'],label='Released WitnessCell',zorder=5)
    c.axhline(1,color=INK,ls='--',lw=.7);c.set_ylim(.991,1.025);c.set_ylabel('Projected error / truth energy');c.set_xticks(x,[r['seed'] for r in rr]);c.set_xlim(-.5,11.5)
    for xx,lab in [(1,'Norman'),(4,'Wessels'),(7,'Schmidt'),(10,'Replogle exp6')]:c.text(xx,-.25,lab,transform=c.get_xaxis_transform(),ha='center',fontsize=8)
    c.legend(frameon=False,ncol=2,fontsize=8,loc='upper right');save(f,'Fig3')

def fig5():
    f,aa=plt.subplots(3,2,figsize=(7.7,9.1));f.subplots_adjust(left=.13,right=.96,top=.92,bottom=.09,wspace=.60,hspace=.83)
    a,b,c,d,e,g=aa.ravel();panel(a,'a','Training-measurement adjustment','x')
    for i,m in enumerate(['M1','scalar','R']):
        v=100*H.CONTROL['C_gain_over_M0'][m];a.barh(i,v,color=COL[m],height=.55);a.text(v+.5,i,f'{v:.2f}',va='center',fontsize=8)
    a.set_yticks(range(3),['Unit weight','Uniform shrinkage','Readout weight R']);a.invert_yaxis();a.set_xlim(0,20);a.set_xlabel('RNA contrast error reduction (%)')
    b.axis('off');panel(b,'b','Exact amplitude-noise representation',None)
    b.text(0,.82,r'$m=(Z_1+Z_2)/2,\quad d=(Z_1-Z_2)/2$',fontsize=9.8,transform=b.transAxes)
    b.text(0,.45,r'$w_j=\left[\frac{\|m_j\|^2-\|d_j\|^2}{\|m_j\|^2+\|d_j\|^2}\right]_{[0,1]}$',fontsize=15,transform=b.transAxes)
    b.text(0,.04,'Same rule as the half-measurement formula.\nIt is not an independent information source.',fontsize=8.4,transform=b.transAxes,linespacing=1.6)
    panel(c,'c','K562: same test pairs, more training pairs')
    for m in H.PART:
        y=[100*H.hval('K562',q,m) for q in ['low','high']];c.plot([0,1],y,marker='o',color=COL[m],lw=1.7,label=H.NAMES[m])
        for x,v in enumerate(y):c.annotate(f'{v:.2f}',(x,v),xytext=(-6 if m=='rank1_selected' else 5,7 if m=='rank1_selected' else (-7 if x==0 else -14)),textcoords='offset points',ha='right' if m=='rank1_selected' else 'left',fontsize=8,color=COL[m])
    c.axhline(0,color=INK,lw=.7);c.set_xticks([0,1],['421-422','2,317-2,318']);c.set_xlim(-.4,1.4);c.set_ylim(-2,14);c.set_xlabel('Training pairs per outer fold');c.set_ylabel('Full-response reduction vs M0 (%)');c.legend(frameon=False,fontsize=7.5)
    panel(d,'d','Gains beyond endpoint-additive prediction')
    for m in H.PART:
        y=[100*H.hval('K562',q,m,'projected_gain') for q in ['low','high']];d.plot([0,1],y,marker='D',color=COL[m],lw=1.7,label=H.NAMES[m])
    d.axhline(0,color=INK,lw=.7);d.set_xticks([0,1],['Low budget','High budget']);d.set_xlim(-.25,1.25);d.set_ylabel('Projected error reduction vs M0 (%)')
    panel(e,'e','RNA improvement is concentrated')
    for m in ['M1','R']:
        v=sorted([float(r['C_SSE_gain_vs_M0']) for r in H.GENE if r['model']==m and float(r['C_SSE_gain_vs_M0'])>0],reverse=True);v=100*np.cumsum(v)/sum(v);e.plot(np.arange(1,len(v)+1),v,color=COL[m],lw=1.5,label=m);e.scatter(10,v[9],color=COL[m],s=25)
    e.set_xscale('log');e.set_ylim(0,105);e.set_xlabel('Readouts ranked by positive gain');e.set_ylabel('Cumulative positive reduction (%)');e.legend(frameon=False,fontsize=8)
    panel(g,'f','Measurement budget changes more than rank',None);g.axis('off')
    g.text(0,.80,'Same 632 K562 evaluation pairs\nTwo measurements held out together',transform=g.transAxes,fontsize=8.5,linespacing=1.5)
    g.text(0,.37,'Endpoint-incidence rank: 80 at both budgets\nMedian degree: about 10 to 58-59',transform=g.transAxes,fontsize=8.5,linespacing=1.5)
    g.text(0,-.07,'More observations improve M0 by 10.12%.\nPartner gains are measured against each M0.',transform=g.transAxes,fontsize=8.3,linespacing=1.5)
    save(f,'Fig5')

def remapped(fn,name):
    old=H.save
    H.save=lambda f,n:old(f,name)
    try:fn()
    finally:H.save=old

def main():
    fig1();fig2();fig3();remapped(H.fig2,'Fig4');fig5();E.fig6()
    for fn in [H.ed1,H.ed2,H.ed3,H.ed4,H.ed5,H.ed6,H.ed7]:fn()
    remapped(H.fig3,'ED8');E.ed9();E.ed10()
    H.checks()
    p=ROOT/'figure_code/FIGURE_VALIDATION.json';j=json.loads(p.read_text());j['figures']=['Fig'+str(i) for i in range(1,7)]+['ED'+str(i) for i in range(1,11)];j['version']='v20';j['figure_only_reconstruction']=True;j['validation_scope']='figure regeneration from stored results; executed A/B experiment audits are separate';j['new_v20_executed_experiments']={'A':'completed six-split extension','B':'completed independent-source evaluation; negative primary partner outcome'};j['benchmark_source']='released_benchmark/formal_top100_leaderboard.csv';p.write_text(json.dumps(j,indent=2)+'\n')

if __name__=='__main__':
    main()
    # v22 first-citation order: structure then released benchmark.
    for ext in ('pdf','png','svg'):
        a=OUT/('Fig2.'+ext);b=OUT/('Fig3.'+ext)
        aa=a.read_bytes();a.write_bytes(b.read_bytes());b.write_bytes(aa)
    import runpy
    runpy.run_path(str(ROOT/'figure_code/build_v22_figures.py'),run_name='__main__')
    p=ROOT/'figure_code/FIGURE_VALIDATION.json'
    j=json.loads(p.read_text());j['version']='v22';j['figures']=['Fig'+str(i) for i in range(1,8)]+['ED'+str(i) for i in range(1,11)]
    j['legacy_main_figure_map']={'Fig2':'Fig3','Fig3':'Fig2'}
    j['validation_scope']='Drawing from stored predictions/statistics, without model fitting.'
    p.write_text(json.dumps(j,indent=2)+'\n')

