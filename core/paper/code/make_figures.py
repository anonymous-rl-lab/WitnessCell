#!/usr/bin/env python3
"""Rebuild five v31 main and three supplementary figures from frozen source data. No training or network."""
from pathlib import Path
import csv,json,hashlib,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
from matplotlib.lines import Line2D
from pypdf import PdfReader,PdfWriter

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'evidence/source_data/publication'; F=ROOT/'evidence/source_data/formal_bridge'; G=ROOT/'evidence/source_data/gse220974'
OUT=ROOT/'figures'; D=ROOT/'figure_source_data'
OUT.mkdir(exist_ok=True);D.mkdir(exist_ok=True)
INK='#233747'; TEAL='#087F8C'; BLUE='#4477AA'; GOLD='#D59435'; PURPLE='#8567A8'; GRAY='#8998A5'; PALE='#EEF3F6'; RED='#B76058'; WHITE='#FFFFFF'
plt.rcParams.update({'font.family':'Liberation Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
 'xtick.labelsize':7.3,'ytick.labelsize':7.3,'axes.edgecolor':'#9AA8B3','axes.labelcolor':INK,'text.color':INK,
 'xtick.color':INK,'ytick.color':INK,'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65,
 'xtick.major.size':3,'ytick.major.size':3,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'mathtext.fontset':'dejavusans','savefig.facecolor':'white','figure.facecolor':'white','axes.unicode_minus':False})
AUDIT=[]

def read(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as h:return list(csv.DictReader(h))
def dump(name,rows):
    with (D/name).open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def panel(fig,rect,letter,title,dy=.036):
    ax=fig.add_axes(rect)
    fig.text(rect[0]-.035,rect[1]+rect[3]+dy,letter,fontweight='bold',fontsize=11,va='bottom')
    fig.text(rect[0],rect[1]+rect[3]+dy,title,fontweight='bold',fontsize=9,va='bottom')
    return ax
def style(ax,axis='x'):
    for side in ['top','right']:ax.spines[side].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis=axis,color='#E8EDF0',lw=.6)
def off(ax):ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
def box(ax,x,y,w,h,text,fc=PALE,ec='none',size=8.3,bold=False):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.025',facecolor=fc,edgecolor=ec,linewidth=.75)
    ax.add_patch(p);ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=size,fontweight='bold' if bold else 'normal',linespacing=1.35)
    return p
def arrow(ax,a,b,color=GRAY,style='-',lw=1.05):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,linewidth=lw,color=color,linestyle=style,connectionstyle='arc3,rad=0'))
def node(ax,x,y,label,color=TEAL,r=.034):
    ax.add_patch(Circle((x,y),r,fc='white',ec=color,lw=1.2));ax.text(x,y,label,ha='center',va='center',fontsize=8.3,color=color,fontweight='bold')
def save(fig,num,name,inputs):
    stem=(f'Fig{num}_{name}' if num<=5 else f'SuppFig{num-5}_{name}')
    directory=OUT if num<=5 else ROOT/'supplementary_figures'
    directory.mkdir(exist_ok=True)
    fig.savefig(directory/(stem+'.pdf'))
    fig.savefig(directory/(stem+'.svg'))
    fig.savefig(directory/(stem+'.png'),dpi=600)
    fig.savefig(directory/(stem+'_preview.png'),dpi=160)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer(); outside=[]
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():continue
        b=text.get_window_extent(renderer); fw,fh=fig.bbox.width,fig.bbox.height
        if b.x0 < -1 or b.y0 < -1 or b.x1 > fw+1 or b.y1 > fh+1:outside.append(text.get_text())
    AUDIT.append({'figure':num,'stem':stem,'directory':directory.name,'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight(),'png_dpi':600,
                  'sources':inputs,'text_outside_page':outside})
    plt.close(fig)

def fig1():
    fig=plt.figure(figsize=(7.15,8.0))
    ax=panel(fig,[.07,.60,.88,.32],'a','Build a response prediction from measured support');off(ax)
    box(ax,.00,.74,.20,.19,'Controls\nReference response',fc=PALE)
    box(ax,.27,.74,.23,.19,'Single perturbations\nObserved / inferred effects',fc='#EAF2F8')
    box(ax,.60,.74,.28,.19,'Measured combinations\nTraining responses',fc='#E5F2F0')
    arrow(ax,(.38,.72),(.30,.61),BLUE);arrow(ax,(.74,.72),(.74,.61),TEAL)
    arrow(ax,(.10,.72),(.17,.61),GRAY)
    box(ax,.02,.36,.48,.23,'Nonlinear baseline\n'+r'$B_{gh}=\mu_0+S_\alpha(e_g+e_h)$',fc='#EAF2F8',size=9)
    box(ax,.60,.36,.36,.23,'Endpoint residual transfer\n'+r'$\delta_Q=K^{\rm end}_{QT}(K^{\rm end}_{TT}+\eta I)^{-1}R_T$',fc='#E5F2F0',size=8.1)
    arrow(ax,(.51,.48),(.59,.48),GRAY)
    ax.text(.55,.54,'residuals',fontsize=6.9,ha='center',color=GRAY)
    arrow(ax,(.25,.34),(.42,.20),BLUE);arrow(ax,(.77,.34),(.59,.20),TEAL)
    box(ax,.23,.02,.56,.17,r'$\widehat{Y}_{\rm core}=B+\gamma\delta$'+'\nReleased predictor: formal benchmark',fc='#F0F5F7',size=9.2,bold=False)
    # Endpoint / contrast decomposition uses only exact incidence identities.
    ax=panel(fig,[.07,.285,.40,.24],'b','Separate response components',dy=.024);off(ax)
    ax.text(.0,.94,'Query pairs',fontsize=8,fontweight='bold');ax.text(.56,.94,'Endpoint incidence',fontsize=8,ha='center',fontweight='bold')
    labels=['AB','AD','CB','CD'];X=np.array([[1,1,0,0],[1,0,0,1],[0,1,1,0],[0,0,1,1]])
    for j,l in enumerate(['A','B','C','D']):ax.text(.40+j*.12,.80,l,ha='center',fontsize=8)
    for i,label in enumerate(labels):
        yy=.66-i*.12;ax.text(.03,yy,label,va='center',fontsize=8.5)
        for j in range(4):
            ax.add_patch(Rectangle((.35+j*.12,yy-.043),.10,.086,fc=TEAL if X[i,j] else '#F0F3F5',ec='white',lw=.8))
            if X[i,j]:ax.text(.40+j*.12,yy,'1',ha='center',va='center',color='white',fontsize=8)
    ax.text(.47,.13,r'$P_Q=I-X_QX_Q^+$',ha='center',fontsize=10)
    ax.text(.47,.005,r'$P_Q\delta_Q=0$',ha='center',fontsize=10,color=TEAL)
    ax=panel(fig,[.56,.285,.39,.24],'c','Learn an explicit partner increment',dy=.024);off(ax)
    box(ax,.04,.65,.91,.28,'RNA partner extension\nProjected response-similarity kernel\nTraining-half readout weights',fc='#EAF2F8',size=8.0)
    box(ax,.04,.24,.91,.28,'Functional extension\nShared symmetric residual matrix\nTraining-only rank / penalty selection',fc='#E5F2F0',size=8.0)
    ax.text(.49,.06,'Each extension uses its own fitted baseline',ha='center',fontsize=7.1,color=GRAY)
    ax=panel(fig,[.07,.05,.88,.16],'d','Match the added capacity to a held-out test',dy=.022);off(ax)
    for x,label in [(.02,'Added correction'),(.38,'Complete response'),(.68,'Partner contrast')]:
        ax.text(x,.91,label,fontsize=8,fontweight='bold',va='center')
    ax.axhline(.80,color='#D8E2E8',lw=.8)
    for y,label,partner in [(.64,'Endpoint residual',r'Unchanged: $P_Q\delta_Q=0$'),(.36,'Explicit partner','Can change; test held-out gain')]:
        ax.text(.02,y,label,fontsize=8,va='center')
        ax.text(.38,y,'Can change',fontsize=8,va='center')
        ax.text(.68,y,partner,fontsize=7.6,va='center',color=TEAL)
    rr=read(P/'Official12_MAIN_TABLE_all_genes_raw.csv')
    delta=max(abs(float(r['WitnessCell_v14'])-float(r['v14_joint_saturation_backbone'])) for r in rr)
    dims=sum(int(r['d']) for r in rr)
    assert len(rr)==12 and dims==212 and np.isclose(delta,6.33775e-8,rtol=1e-5)
    ax.text(.02,.06,'Released splits: projected-error ratios agree to '+r'$6.34\times10^{-8}$'+' | 12 query graphs',fontsize=7.2,color=GRAY,va='center')
    dump('Fig1_branch_capacity_check.csv',[{'query_graphs':len(rr),'algebraic_dimensions':dims,'maximum_absolute_projected_error_ratio_change':delta,'reference':'actual jointly selected nonlinear backbone'}])
    save(fig,1,'Method',['Main Methods: released core, endpoint elimination, RNA kernel','Official12_MAIN_TABLE_all_genes_raw.csv','SI Notes 1, 3-6'])

def fig2():
    fig=plt.figure(figsize=(7.15,6.65))
    ax=panel(fig,[.07,.57,.39,.34],'a','Hidden downstream mechanism');off(ax)
    # Known upstream mechanism, hidden downstream response coefficients / gate.
    node(ax,.48,.86,'X',INK,.05)
    pos={c:(.08+i*.165,.56) for i,c in enumerate('ABCDEF')}
    for c,(x,y) in pos.items():arrow(ax,(.48,.80),(x,y+.043),GRAY,lw=.75);node(ax,x,y,c,BLUE,.040)
    box(ax,.02,.18,.62,.18,'One candidate gate\n'+r'$h_A,\ h_B\ \mathrm{or}\ h_C$',fc='#E5F2F0',size=8.4)
    for c in 'ABC':arrow(ax,(pos[c][0],.51),(.26,.38),TEAL,lw=.8)
    node(ax,.84,.26,'Y',INK,.056);arrow(ax,(.66,.26),(.78,.26),TEAL)
    arrow(ax,(.53,.83),(.87,.33),GOLD,'--',1)
    ax.text(.88,.70,'Direct\nsource term',ha='center',va='center',fontsize=7,color=GOLD)
    ax.text(.02,.02,'Affine mediator terms in every candidate',fontsize=7.3,color=GRAY)
    ax=panel(fig,[.56,.57,.39,.34],'b','Equal budget; different joint information');off(ax)
    layouts=[('Comparison',[(0,3),(1,3),(1,4),(2,5)],GRAY,.02),('Informative',[(0,1),(1,2),(3,4),(3,5)],TEAL,.55)]
    for name,edges,col,start in layouts:
        pts=[(start+.035+j*.15,.66) for j in range(3)]+[(start+.035+j*.15,.38) for j in range(3)]
        ax.text(start+.18,.91,name,ha='center',fontsize=8.3,fontweight='bold',color=col)
        for a,b in edges:ax.plot([pts[a][0],pts[b][0]],[pts[a][1],pts[b][1]],color=col,lw=1.8,zorder=1)
        for c,(x,y) in zip('ABCDEF',pts):node(ax,x,y,c,col,.043)
    ax.text(.48,.16,'Same endpoint degrees  |  12 conditions',ha='center',fontsize=7.2)
    ax.text(.48,.035,'256 observations / condition',ha='center',fontsize=7.6,color=GRAY)
    ax=panel(fig,[.18,.15,.29,.29],'c','Target-relevant population contrasts',dy=.05)
    vals=np.array([[0,-.5,-.5,0],[0,0,-.5,-.5],[0,-.5,0,-.5]])
    cmap=colors.LinearSegmentedColormap.from_list('truth',['#B8DDDC','#FFFFFF'])
    ax.imshow(vals,cmap=cmap,vmin=-.5,vmax=0,aspect='auto')
    for i in range(3):
        for j in range(4):ax.text(j,i,'0' if vals[i,j]==0 else '-0.5',ha='center',va='center',fontsize=9)
    ax.set_xticks(range(4),['Direct',r'$h_A$',r'$h_B$',r'$h_C$']);ax.set_yticks(range(3),['AB contrast','BC contrast',r'Query $\Delta$'])
    ax.axhline(1.5,color=INK,lw=1);ax.tick_params(length=0,pad=6)
    for s in ax.spines.values():s.set_visible(False)
    ax.set_xlabel('Generating world',labelpad=9)
    fig.text(.09,.065,r'$\Delta=\mu_{AC}-\mu_{AF}-\mu_{CE}+\mu_{EF}$',fontsize=9)
    ax=panel(fig,[.65,.15,.29,.29],'d','Target prediction risk',dy=.05)
    rs={r['method']:r for r in read(F/'scores/summary.csv')};names=['OLS_low','OLS_high','contrast_GLS_high'];yy=[2,1,0]
    for m,y,col in zip(names,yy,[GRAY,TEAL,GOLD]):
        r=rs[m];val=float(r['four_world_rectangle_MSE']);lo=float(r['MC95_low']);hi=float(r['MC95_high'])
        ax.errorbar(val,y,xerr=[[val-lo],[hi-val]],fmt='o',color=col,ms=5,capsize=3,lw=1.2)
        ax.text(val+.0045,y,f'{val:.4f}',va='center',fontsize=7.4,color=col)
    ax.axvline(.0625,color=RED,ls='--',lw=.9)
    ax.set_xlim(0,.087);ax.set_ylim(-.65,2.65);ax.set_yticks(yy,['Comparison','Informative','GLS']);ax.set_xticks([0,.02,.04,.06,.08]);ax.set_xlabel('Mean squared target error');style(ax)
    ax.text(.02,1.07,'92.84% lower risk',transform=ax.transAxes,color=TEAL,fontweight='bold',fontsize=8.2)
    ax.text(.0625,2.40,'Risk floor',color=RED,fontsize=7,ha='center')
    save(fig,2,'Joint_Intervention_Information',['formal_bridge/scores/summary.csv','formal_bridge/theory/THREE_GATE_DESIGN_ENUMERATION.json'])


def fig3():
    fig=plt.figure(figsize=(7.15,7.4));rr=read(ROOT/'revision/BENCHMARK_DECOMPOSITION.csv');idx={(r['population'],r['method']):r for r in rr}
    methods=sorted({r['method'] for r in rr},key=lambda m:float(idx['all',m]['mean_operation_rank']));n=len(methods)
    ax=panel(fig,[.23,.43,.31,.48],'a','Full track and double-only rank')
    for i,m in enumerate(methods):
        x=float(idx['all',m]['mean_operation_rank']);y=float(idx['double',m]['mean_operation_rank'])
        ax.plot([x,y],[i,i],color=GRAY,lw=.8);ax.plot(x,i,'o',mfc='white',mec=GRAY,ms=4);ax.plot(y,i,'o',color=TEAL if m=='WitnessCell' else BLUE,ms=4)
        if m=='WitnessCell':ax.text(y+.35,i-.35,f'{y:.3f}',color=TEAL,fontweight='bold',fontsize=7)
    ax.set_yticks(range(n),[m.replace('scouter','Scouter') for m in methods]);ax.set_ylim(n-.5,-.8);ax.set_xlim(0,16.4);ax.set_xticks([0,4,8,12,16]);ax.set_xlabel('Mean operation rank');style(ax)
    ax.get_yticklabels()[0].set_color(TEAL);ax.get_yticklabels()[0].set_fontweight('bold')
    ax.legend(handles=[Line2D([0],[0],marker='o',ls='',mfc='white',mec=GRAY,label='All 654'),Line2D([0],[0],marker='o',ls='',color=BLUE,label='Doubles 540')],frameon=False,fontsize=6.8,loc='lower left',bbox_to_anchor=(-.03,-.18),ncol=2,handlelength=1,columnspacing=1)
    bx=panel(fig,[.65,.43,.29,.48],'b','Double-only exposure positions')
    states=['combo_seen0','combo_seen1','combo_seen2'];v=np.array([[int(idx[p,m]['position']) for p in states] for m in methods]);cmap=colors.LinearSegmentedColormap.from_list('ranks',['#C5E5E1','#F9FBFB','#F0DFCE'])
    bx.imshow(v,vmin=1,vmax=16,cmap=cmap,aspect='auto')
    for i in range(n):
        for j in range(3):bx.text(j,i,str(v[i,j]),ha='center',va='center',fontsize=7,fontweight='bold' if i==0 else 'normal')
    bx.set_xticks(range(3),['0 seen\nn = 43','1 seen\nn = 374','2 seen\nn = 123']);bx.set_yticks([]);bx.tick_params(length=0);[p.set_visible(False) for p in bx.spines.values()]
    ax=panel(fig,[.19,.09,.34,.20],'c','Double-only paired differences')
    cr={r['reference']:r for r in read(ROOT/'revision/BENCHMARK_PAIRED_INTERVALS.csv') if r['population']=='double'}
    for i,m in enumerate(['scouter','linearModel','GEARS']):
        r=cr[m];v=float(r['difference_WitnessCell_minus_reference']);lo=float(r['low']);hi=float(r['high']);ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color=TEAL,ms=4.5,capsize=2.5,lw=1.1)
    ax.axvline(0,color=GRAY,ls='--',lw=.8);ax.set_yticks(range(3),['Scouter','Linear','GEARS']);ax.set_ylim(2.5,-.5);ax.set_xlim(-5.1,.3);ax.set_xlabel('WitnessCell minus reference rank');style(ax)
    ax.text(.02,-.43,'95% intervals; 287 dataset-pair clusters',transform=ax.transAxes,fontsize=6.7,color=GRAY)
    ax=panel(fig,[.72,.09,.22,.20],'d','Native-scale MSE')
    dd={r['population']:r for r in read(ROOT/'revision/NATIVE_MSE_VS_SCOUTER.csv')}
    for i,(ds,label) in enumerate([('Norman','Norman'),('Replogle_exp6','Replogle'),('Wessels','Wessels'),('Schmidt','Schmidt')]):
        v=100*float(dd[ds+'_double']['relative_MSE_gain']);col=TEAL if v>0 else RED
        ax.barh(i,v,height=.48,color=col)
        ax.text(v+(2 if v>=0 else -2),i,f'{v:+.1f}',ha='left' if v>=0 else 'right',va='center',fontsize=6.9,color=col)
    ax.axvline(0,color=GRAY,lw=.8);ax.set_yticks(range(4),['Norman','Replogle','Wessels','Schmidt']);ax.set_ylim(3.6,-.6);ax.set_xlim(-52,66);ax.set_xticks([-40,0,40]);ax.set_xlabel('MSE reduction vs Scouter (%)',fontsize=7.3);style(ax)
    save(fig,3,'Prediction_Benchmark',['revision/BENCHMARK_DECOMPOSITION.csv','revision/BENCHMARK_PAIRED_INTERVALS.csv','revision/NATIVE_MSE_VS_SCOUTER.csv'])

def rna_figure():
    fig=plt.figure(figsize=(7.15,8.0));rr=read(P/'Step2_Fig5_Matched55.csv');order=['additive','saturation_only','M0','M1','R','GEARS','scouter'];idx={r['model']:r for r in rr}
    labels=['Additive','Saturation','M0','M1','R','GEARS','Scouter']; cols=[GRAY,GOLD,INK,BLUE,TEAL,PURPLE,PURPLE]
    seed=read(P/'Step2_Fig5_SeedMetrics.csv')
    for rect,letter,title,key,xlim,log in [([.20,.70,.27,.22],'a','Complete RNA response','response_MSE',(1,35),True),([.66,.70,.28,.22],'b','Partner contrasts','C_MSE',(3.35,4.55),False)]:
        ax=panel(fig,rect,letter,title)
        for i,(m,col) in enumerate(zip(order,cols)):
            v=float(idx[m][key])*1000;ax.plot(v,i,'o',ms=5,color=col,zorder=3)
            if m in ['GEARS','scouter']:
                ss=[r for r in seed if r['model']==m and r['variant']=='raw']
                assert len(ss)==3
                values=[float(r[key])*1000 for r in ss]
                ax.hlines(i,min(values),max(values),color=col,alpha=.45,lw=1.2)
                ax.scatter(values,[i]*len(values),s=12,facecolor='white',edgecolor=col,lw=.7,zorder=2)
        if log:ax.set_xscale('log');ax.set_xticks([1,2,5,10,30],['1','2','5','10','30']);ax.minorticks_off()
        ax.set_xlim(xlim);ax.set_ylim(6.6,-.6);ax.set_yticks(range(7),labels);ax.set_xlabel(r'MSE $\times 10^3$'+(' (log axis)' if log else ''));style(ax)
        ax.get_yticklabels()[4].set_color(TEAL);ax.get_yticklabels()[4].set_fontweight('bold')
        if letter=='b':ax.text(.99,.95,'15.96%\nlower than M0',ha='right',va='top',transform=ax.transAxes,color=TEAL,fontsize=7.4,fontweight='bold')
    dat=np.load(P/'RNA_signed_rectangle_predictions.npz');truth=dat['truth_effect'];e0=(dat['M0']-truth)**2;eR=(dat['R']-truth)**2
    ax=panel(fig,[.13,.39,.34,.21],'c','Eleven same-fit rectangles')
    x=e0.mean(1)*1000;y=eR.mean(1)*1000;lim=max(max(x),max(y))*1.12
    ax.plot([0,lim],[0,lim],color=GRAY,ls='--',lw=.8)
    ax.scatter(x,y,c=[TEAL if v<=u else RED for u,v in zip(x,y)],s=27,edgecolors='white',lw=.6,zorder=3)
    ax.set_xlim(0,lim);ax.set_ylim(0,lim);ax.set_xlabel(r'M0 rectangle MSE $\times 10^3$');ax.set_ylabel(r'R rectangle MSE $\times 10^3$');style(ax,'both')
    ax.text(.04,.94,'9 / 11 improved',transform=ax.transAxes,va='top',color=TEAL,fontsize=8)

    ax=panel(fig,[.65,.39,.29,.21],'d','Shrinkage reduces adverse corrections')
    budget={r['model']:r for r in read(ROOT/'revision/RNA_ERROR_BUDGET.csv') if r['readout_set']=='all'}
    for i,m in enumerate(['M1','scalar','R']):
        r=budget[m];p=float(r['positive_percent_of_all_M0']);a=float(r['adverse_percent_of_all_M0']);net=p-a
        ax.bar(i,p,width=.58,color=TEAL,alpha=.82);ax.bar(i,-a,width=.58,color=RED,alpha=.85);ax.plot(i,net,'D',color=INK,ms=4)
        ax.text(i,net-2.5,f'{net:.2f}',ha='center',va='top',fontsize=7,color=INK)
    ax.axhline(0,color=GRAY,lw=.8);ax.set_xticks(range(3),['M1','Uniform','R']);ax.set_ylabel('Share of total M0 SSE (%)');ax.set_ylim(-15,28);style(ax,'y')
    ax.legend(handles=[Rectangle((0,0),1,1,fc=TEAL,label='Positive'),Rectangle((0,0),1,1,fc=RED,label='Adverse'),Line2D([0],[0],marker='D',color=INK,ls='',label='Net')],loc='upper right',frameon=False,fontsize=6.1,ncol=3,columnspacing=.5,handlelength=.8,handletextpad=.3)
    ax=panel(fig,[.13,.075,.34,.21],'e','Gains concentrate in a few readouts')
    conc=[]
    for m,col in [('M1',BLUE),('R',TEAL)]:
        delta=(e0-(dat[m]-truth)**2).sum(0);positive=np.sort(np.maximum(delta,0))[::-1];cs=np.cumsum(positive)/positive.sum()*100
        ax.plot(np.arange(1,len(cs)+1),cs,color=col,lw=1.6,label=m);ax.plot(10,cs[9],'o',color=col,ms=4)
        conc.append({'model':m,'top10_positive_gain_share':cs[9]/100})
    ax.axvline(10,color=GRAY,ls=':',lw=.8);ax.set_xscale('log');ax.set_xlim(1,4995);ax.set_ylim(0,104);ax.set_xticks([1,10,100,1000,4995],['1','10','100','1,000','4,995']);ax.minorticks_off()
    ax.set_xlabel('Readouts ranked by positive gain');ax.set_ylabel('Cumulative positive gain (%)');style(ax,'y');ax.legend(frameon=False,loc='lower right',fontsize=7)
    ax.text(14,72,'R: 92.13% in top 10',color=TEAL,fontsize=7.3)
    dump('Fig5_positive_gain_concentration.csv',conc)
    ax=panel(fig,[.65,.075,.29,.21],'f','Readout-scale sensitivity')
    xx=np.arange(3);keys=['C_gain_vs_M0','standardized_C_gain_vs_M0','drop10_C_gain_vs_M0']
    for m,col,offset in [('M1',BLUE,-.15),('R',TEAL,.15)]:
        vals=[100*float(idx[m][k]) for k in keys];ax.bar(xx+offset,vals,width=.27,color=col,label=m,zorder=3)
    ax.axhline(0,color=GRAY,lw=.8);ax.set_xticks(xx,['Raw','Across-condition\nstandardized','Ten-readout\ndeletion']);ax.set_ylabel('Rectangle-error reduction (%)');ax.set_ylim(-14,20);ax.set_yticks([-10,0,10,20]);style(ax,'y');ax.legend(frameon=False,fontsize=7,loc='upper right')
    save(fig,5,'RNA_Partner_Response',['Step2_Fig5_Matched55.csv','Step2_Fig5_SeedMetrics.csv','RNA_signed_rectangle_predictions.npz','revision/RNA_ERROR_BUDGET.csv','RNA_controls.json'])

def growth_gain(stage,budget,model,key='sse',repeat=None):
    rr=[r for r in read(P/'Horlbeck_figure_metrics.csv') if r['stage']==stage and r['budget']==budget and r['fold']!='pooled' and r['repeat']!='pooled' and (repeat is None or r['repeat']==str(repeat))]
    s=lambda m:sum(float(r[key]) for r in rr if r['model']==m)
    return 100*(1-s(model)/s('M0'))

def growth_figure():
    fig=plt.figure(figsize=(7.15,9.1))
    ax=panel(fig,[.07,.80,.88,.12],'a','Complete a partial functional screen',dy=.025);off(ax)
    box(ax,.00,.14,.26,.69,'Measured functional screen\n2,464 pairs (~80%) + singles\nTwo source repeats',fc=PALE,size=7.9)
    box(ax,.39,.14,.23,.69,'Fit baseline + partner\nSelect from training data\nRank one is primary',fc='#EAF2F8',size=8.2)
    box(ax,.75,.14,.24,.69,'617 query pairs (~20%)\nBaseline + partner increment\nComplete response',fc='#E5F2F0',size=7.9)
    arrow(ax,(.28,.48),(.37,.48),BLUE);arrow(ax,(.64,.48),(.73,.48),TEAL)
    ax.text(.49,-.16,'Jurkat: 79 endpoints, 3,081 possible pairs  |  Procedure frozen before refitting',ha='center',fontsize=7.1,color=GRAY)
    inp=json.loads((ROOT/'software/growth_example_input.json').read_text())
    total=len(inp['genes'])*(len(inp['genes'])-1)//2
    dump('Fig4_support_split.csv',[{'partition':key,'pairs':len(inp[key]),'fraction':len(inp[key])/total,'endpoints':len(inp['genes']),'possible_pairs':total} for key in ['training','queries']])
    ax=panel(fig,[.13,.52,.34,.19],'b','K562: expanded combination support',dy=.023)
    for m,col,l in [('rank1_selected',TEAL,'Rank one'),('spectral_selected',BLUE,'Spectral')]:
        for key,ls,mk in [('sse','-','o'),('projected_sse','--','s')]:
            vals=[growth_gain('K562',b,m,key) for b in ['low','high']]
            offset=(-.05 if m=='rank1_selected' else .05)+(-.025 if key=='sse' else .025)
            ax.plot(np.array([0,1])+offset,vals,ls='none',marker=mk,color=col,ms=4.5,mfc=col if key=='sse' else 'white')
    ax.axhline(0,color=GRAY,lw=.7);ax.set_xticks([0,1],['421-422','2,317-2,318']);ax.set_xlim(-.12,1.12);ax.set_ylim(-2,15);ax.set_ylabel('Error reduction vs same-budget M0 (%)');ax.set_xlabel('Training pairs per fold');style(ax,'y')
    handles=[Line2D([0],[0],color=TEAL,marker='o',ls='',label='Rank one'),Line2D([0],[0],color=BLUE,marker='o',ls='',label='Spectral'),Line2D([0],[0],color=GRAY,marker='o',ls='',label='Complete'),Line2D([0],[0],color=GRAY,marker='s',mfc='white',ls='',label='Partner')]
    ax.legend(handles=handles,frameon=False,ncol=2,fontsize=6.9,loc='upper left',columnspacing=.8,handlelength=1.4)
    ax=panel(fig,[.65,.52,.29,.19],'c','Jurkat: held-out functional gain',dy=.023)
    ints={r['model']:r for r in read(P/'Horlbeck_Jurkat_intervals.csv')}
    for y,m,col in [(1,'rank1_selected',TEAL),(0,'spectral_selected',BLUE)]:
        v=growth_gain('Jurkat','high',m);lo=100*float(ints[m]['lower']);hi=100*float(ints[m]['upper'])
        ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='o',color=col,ms=5.5,lw=1.3,capsize=3)
        ax.text(v,y+.25,f'{v:.2f}%',ha='center',color=col,fontweight='bold',fontsize=8.5)
    ax.axvline(0,color=GRAY,ls='--',lw=.8);ax.set_yticks([1,0],['Rank one','Spectral']);ax.set_ylim(-.6,1.6);ax.set_xlim(-.5,15);ax.set_xlabel('Complete-error reduction (%)');style(ax)
    ax=panel(fig,[.13,.285,.34,.135],'d','Jurkat: repeat and partner outcomes',dy=.023)
    rows=[]
    for i,(key,label,repeat) in enumerate([('sse','Repeat 1',1),('sse','Repeat 2',2),('projected_sse','Partner contrast',None)]):
        for m,col,shift in [('rank1_selected',TEAL,-.09),('spectral_selected',BLUE,.09)]:
            v=growth_gain('Jurkat','high',m,key,repeat);ax.plot(v,i+shift,'o',ms=5,color=col);ax.text(v+.28,i+shift,f'{v:.2f}',va='center',color=col,fontsize=7);rows.append({'metric':label,'model':m,'gain_percent':v})
    ax.set_yticks(range(3),['Repeat 1','Repeat 2','Partner contrast']);ax.set_ylim(2.5,-.5);ax.set_xlim(0,11);ax.set_xlabel('Error reduction (%)');style(ax)
    dump('Fig4_growth_gains.csv',rows)
    # Four corners come from the same frozen Jurkat fit. The extra background
    # is selected by name, not by observed or predicted outcome.
    ex=[r for r in read(ROOT/'revision/JURKAT_FOUR_CORNERS_V28.csv') if r['background']=='ACTR10']
    ax=panel(fig,[.64,.275,.30,.145],'e','Four same-fit response corners',dy=.023)
    styles=[('observed',INK,'D','Observed'),('M0',GRAY,'o','M0'),('rank1_selected',TEAL,'o','Rank one')]
    plotrows=[]
    for i,corner in enumerate(['AB','AD','CB','CD']):
        rr=[r for r in ex if r['corner']==corner]
        for k,col,m,label in styles:
            v=np.mean([float(r[k]) for r in rr]);ax.plot(i,v,m,color=col,ms=4.5,label=label if i==0 else None)
            plotrows.append(dict(corner=corner,gene_left=rr[0]['gene_left'],gene_right=rr[0]['gene_right'],component=k,repeat_mean=float(v)))
    ax.set_xticks(range(4),['A+B','A+D','C+B','C+D']);ax.set_xlim(-.4,3.4)
    ax.set_ylabel('Mean growth response',fontsize=7.2);style(ax,'y')
    fig.text(.64,.232,'A = TTI2; C = ACTR10\nB = PGK1; D = OGFOD1',fontsize=6.8,color=GRAY,linespacing=1.4)
    fig.legend(handles=[Line2D([0],[0],marker=m,color=c,ls='',ms=4,label=label) for _,c,m,label in styles],loc='upper left',bbox_to_anchor=(.635,.212),frameon=False,ncol=3,fontsize=6.3,handletextpad=.25,columnspacing=.6,handlelength=.9)
    ax=panel(fig,[.13,.055,.81,.105],'f','Endpoint-cancelled partner contrast',dy=.023)
    cc=[r for r in read(ROOT/'revision/JURKAT_FOUR_CORNER_CONTRASTS_V28.csv') if r['background']=='ACTR10']
    for j,(k,col,m,label) in enumerate(styles):
        mean=float(next(r[k] for r in cc if r['repeat']=='mean'))
        for rep,offset in [('1',-.08),('2',.08)]:
            v=float(next(r[k] for r in cc if r['repeat']==rep));ax.plot(v,j+offset,'o',mfc='white',mec=col,ms=3.2,alpha=.85)
        ax.plot(mean,j,m,color=col,ms=6);ax.text(mean+.003,j+.21,f'{mean:.4f}',fontsize=7,color=col,va='center')
    ax.axvline(0,color=GRAY,ls='--',lw=.7);ax.set_xlim(-.073,.008);ax.set_ylim(2.6,-.6)
    ax.set_yticks([0,1,2],['Observed','M0','Rank one']);ax.set_xlabel(r'$\Delta=Y_{AB}-Y_{AD}-Y_{CB}+Y_{CD}$'+' (growth units)',fontsize=8);style(ax)
    fig.text(.58,.175,'Large: repeat mean   |   Small: source repeats',fontsize=6.8,color=GRAY)
    dump('Fig4_four_corner_means.csv',plotrows)
    save(fig,4,'Functional_Growth',['Horlbeck_figure_metrics.csv','Horlbeck_Jurkat_intervals.csv','revision/JURKAT_FOUR_CORNERS_V28.csv','revision/JURKAT_FOUR_CORNER_CONTRASTS_V28.csv'])


def applications_figure():
    fig=plt.figure(figsize=(7.15,8.1))
    ax=panel(fig,[.07,.79,.88,.13],'a','PairPick: from an objective to a candidate partner');off(ax)
    box(ax,.00,.18,.24,.63,'Anchor single response\nDefine objective  '+r'$v_a$',fc='#EAF2F8',size=8.3)
    box(ax,.36,.18,.27,.63,'Predict candidate responses\n'+r'$\widehat U(a,b)=v_a^\top\widehat y_{ab}$',fc='#E5F2F0',size=8.4)
    box(ax,.75,.18,.24,.63,'Select a partner\nEvaluate held-out utility',fc=PALE,size=8.2)
    arrow(ax,(.26,.49),(.34,.49),TEAL);arrow(ax,(.65,.49),(.73,.49),TEAL)
    ax=panel(fig,[.16,.51,.28,.17],'b','Separate risk-selection evaluation')
    risk=[{'loss':'MSE','ratio':.5343,'low':.3646,'high':.7591,'records':213,'coverage':49.48},{'loss':'Weighted MSE','ratio':.4775,'low':.3175,'high':.7071,'records':207,'coverage':48.95}]
    for y,r,col in zip([1,0],risk,[TEAL,BLUE]):
        v=r['ratio'];ax.errorbar(v,y,xerr=[[v-r['low']],[r['high']-v]],fmt='o',color=col,ms=5,capsize=3,lw=1.2);ax.text(v,y+.3,f'{v:.3f}',ha='center',fontsize=8,color=col)
    ax.axvline(1,color=GRAY,ls='--',lw=.8);ax.set_xlim(.2,1.08);ax.set_ylim(-.6,1.7);ax.set_yticks([1,0],['MSE','Weighted MSE']);ax.set_xticks([.25,.5,.75,1]);ax.set_xlabel('Accepted / all-query error');style(ax)
    ax.text(.0,-.43,'49.48% / 48.95% coverage\n33 / 32 calibration-disjoint pair identities',transform=ax.transAxes,fontsize=6.6,color=GRAY,linespacing=1.4)
    dump('SuppFig1_risk_selection.csv',risk)
    ax=panel(fig,[.60,.16,.35,.52],'c','All 28 detected-label tasks')
    rr=read(ROOT/'revision/PAIRPICK_ALL_28_TASKS.csv');ms=['R_G','StandardRidge_RNA','Scouter_Ensemble','GEARS_Ensemble','SingleB_additive'];vv=np.array([[1000*float(r[m+'_regret']) for m in ms] for r in rr]);cmap=colors.LinearSegmentedColormap.from_list('regret',['#F5F8F8','#C7DEDD','#6CAAA5','#166B6D'])
    im=ax.imshow(vv,aspect='auto',cmap=cmap,norm=colors.AsinhNorm(linear_width=5,vmin=0,vmax=vv.max()));ax.set_yticks(range(28),[r['anchor']+' / '+r['fold'] for r in rr],fontsize=6.1);ax.set_xticks(range(5),['PairPick','Linear','Scouter','GEARS','Additive'],rotation=40,ha='right',fontsize=7)
    for i in range(28):
        for j in range(5):
            if vv[i,j]<1e-9:ax.plot(j,i,'.',color='#A1B4B7',ms=1.7)
            else:ax.text(j,i,f'{vv[i,j]:.1f}',ha='center',va='center',fontsize=5.5,color='white' if vv[i,j]>40 else INK)
        if i and rr[i]['physical_anchor']!=rr[i-1]['physical_anchor']:ax.axhline(i-.5,color='white',lw=2)
    ax.tick_params(length=0);[p.set_visible(False) for p in ax.spines.values()]
    cax=fig.add_axes([.66,.065,.23,.012]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_label(r'Regret $\times 10^3$ (lower is better)',fontsize=7);cb.set_ticks([0,5,25,150]);cb.set_ticklabels(['0','5','25','150']);cb.ax.tick_params(labelsize=6.5)
    ax=panel(fig,[.20,.16,.24,.20],'d','Paired recommendation utility')
    prs={r['reference']:r for r in read(G/'PAIRED_COMPARISONS.csv')};refs=['PRESAGE_Combo_seed24','StandardRidge_RNA','Scouter_Ensemble','GEARS_Ensemble','SingleB_additive'];rl=['PRESAGE','Linear','Scouter','GEARS','Additive']
    for i,m in enumerate(refs):
        r=prs[m];v=float(r['mean_paired_U_difference'])*1000;lo=float(r['cluster_bootstrap_low'])*1000;hi=float(r['cluster_bootstrap_high'])*1000;ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color=TEAL,ms=4,capsize=2.5,lw=1)
    ax.axvline(0,color=GRAY,ls='--',lw=.8);ax.set_yticks(range(5),rl);ax.set_ylim(4.6,-.6);ax.set_xlim(-9,25);ax.set_xticks([-5,0,10,20]);ax.set_xlabel(r'Utility difference $\times 10^3$');style(ax)
    ax.text(-.02,-.40,'R-G minus reference; 95% intervals\nEight physical anchors',transform=ax.transAxes,fontsize=6.6,color=GRAY,linespacing=1.4)
    fig.text(.60,.762,'15 / 49 labels match designed gene-mode pairs\nRetrospective demonstration; see identity audit',fontsize=6.9,color=INK,linespacing=1.4)
    save(fig,6,'Research_Applications',['SI Table S16 (risk, displayed precision)','revision/PAIRPICK_ALL_28_TASKS.csv','gse220974/PAIRED_COMPARISONS.csv','revision/GSE220974_IDENTITY_AUDIT.csv'])


def all_intervals_figure():
    fig=plt.figure(figsize=(7.15,6.4))
    ax=panel(fig,[.24,.14,.68,.76],'a','Double-only comparisons: all fifteen reference entries')
    rr=sorted(read(ROOT/'revision/DOUBLE_ONLY_ALL_COMPARATOR_INTERVALS.csv'),key=lambda r:float(r['difference_WitnessCell_minus_reference']),reverse=True)
    for i,r in enumerate(rr):
        v=float(r['difference_WitnessCell_minus_reference']);lo=float(r['low']);hi=float(r['high'])
        ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color=TEAL,ms=4.3,capsize=2.5,lw=1.0)
    ax.axvline(0,color=GRAY,ls='--',lw=.9);ax.set_yticks(range(len(rr)),[r['reference'].replace('scouter','Scouter') for r in rr]);ax.set_ylim(14.6,-.6);ax.set_xlim(-12,.5);ax.set_xticks([-12,-10,-8,-6,-4,-2,0])
    ax.set_xlabel('WitnessCell minus reference mean operation rank');style(ax)
    fig.text(.24,.035,'Pointwise 95% cluster intervals; fixed scoring archive\n540 records | 287 dataset-pair clusters | 20,000 draws',fontsize=7.2,color=GRAY,linespacing=1.4)
    save(fig,7,'All_Comparator_Intervals',['revision/DOUBLE_ONLY_ALL_COMPARATOR_INTERVALS.csv'])

def support_figure():
    fig=plt.figure(figsize=(7.15,4.1))
    ax=panel(fig,[.14,.20,.75,.64],'a','Complete-response benefit from measured combination support')
    rr=[r for r in read(P/'Fig3_ED4_witness_gains.csv') if r['reference']=='factorized_prediction']
    for metric,col,label in [('loss_mse',TEAL,'MSE'),('loss_wmse',BLUE,'Weighted MSE')]:
        sub=sorted([r for r in rr if r['metric']==metric],key=lambda r:float(r['fraction']));x=np.array([float(r['fraction']) for r in sub])*100;y=np.array([float(r['gain']) for r in sub])*100
        lo=np.array([float(r['gain_ci95_low']) for r in sub])*100;hi=np.array([float(r['gain_ci95_high']) for r in sub])*100
        ax.plot(x,y,'o-',color=col,ms=5,lw=1.5,label=label);ax.fill_between(x,lo,hi,color=col,alpha=.10,lw=0)
    ax.set_xticks([25,50,100]);ax.set_xlim(20,105);ax.set_ylim(0,53);ax.set_ylabel('Error reduction versus fixed backbone (%)');ax.set_xlabel('Training-combination residual library retained (%)');style(ax,'y');ax.legend(frameon=False,fontsize=8,loc='upper left')
    fig.text(.14,.055,'Wessels | Same nonlinear backbone | Archived pair-bootstrap intervals',fontsize=7.3,color=GRAY)
    save(fig,8,'Combination_Support',['Fig3_ED4_witness_gains.csv'])

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--figures',default='1,2,3,4,5,6,7,8');args=ap.parse_args()
    functions={1:fig1,2:fig2,3:fig3,4:growth_figure,5:rna_figure,6:applications_figure,7:all_intervals_figure,8:support_figure}
    for i in map(int,args.figures.split(',')):
        functions[i]();print('Rendered figure',i,flush=True)
    auditpath=ROOT/'FIGURE_AUDIT.json'
    previous=json.loads(auditpath.read_text()) if auditpath.exists() else []
    done={a['figure'] for a in AUDIT};AUDIT=sorted([a for a in previous if a['figure'] not in done]+AUDIT,key=lambda a:a['figure'])
    auditpath.write_text(json.dumps(AUDIT,indent=2)+'\n')
    for kind,entries in [('Main',[a for a in AUDIT if a['figure']<=5]),('Supplementary',[a for a in AUDIT if a['figure']>5])]:
        writer=PdfWriter()
        for r in entries:
            reader=PdfReader(ROOT/r['directory']/(r['stem']+'.pdf'));writer.add_page(reader.pages[0])
        writer.add_metadata({'/Title':'WitnessCell v31 - '+kind+' Figures','/Author':'WitnessCell','/Subject':'Vector figures; legends supplied separately'})
        with (ROOT/('WitnessCell_v31_'+kind+'_Figures.pdf')).open('wb') as h:writer.write(h)
    from PIL import Image,ImageDraw
    thumbs=[]
    for r in AUDIT:
        if r['figure']>5:continue
        im=Image.open(ROOT/r['directory']/(r['stem']+'_preview.png')).convert('RGB');im.thumbnail((740,835));canvas=Image.new('RGB',(780,895),'white');canvas.paste(im,((780-im.width)//2,40));ImageDraw.Draw(canvas).text((25,12),'FIGURE '+str(r['figure']),fill=INK);thumbs.append(canvas)
    contact=Image.new('RGB',(1560,895*math.ceil(len(thumbs)/2)),'#E6EDF1')
    for i,im in enumerate(thumbs):contact.paste(im,((i%2)*780,(i//2)*895))
    contact.save(ROOT/'Figures_Overview.png')
    print(json.dumps({'outside_page':{r['figure']:r['text_outside_page'] for r in AUDIT},'figures':len(AUDIT)}))
