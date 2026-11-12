#!/usr/bin/env python3
"""Rebuild six v24 main figures from frozen source data. No training or network."""
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
    stem=f'Fig{num}_{name}'
    fig.savefig(OUT/(stem+'.pdf'))
    fig.savefig(OUT/(stem+'.svg'))
    fig.savefig(OUT/(stem+'.png'),dpi=600)
    fig.savefig(OUT/(stem+'_preview.png'),dpi=160)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer(); outside=[]
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():continue
        b=text.get_window_extent(renderer); fw,fh=fig.bbox.width,fig.bbox.height
        if b.x0 < -1 or b.y0 < -1 or b.x1 > fw+1 or b.y1 > fh+1:outside.append(text.get_text())
    AUDIT.append({'figure':num,'stem':stem,'width_inches':fig.get_figwidth(),'height_inches':fig.get_figheight(),'png_dpi':600,
                  'sources':inputs,'text_outside_page':outside})
    plt.close(fig)

def fig1():
    fig=plt.figure(figsize=(7.15,6.65))
    ax=panel(fig,[.07,.52,.88,.40],'a','Build a response prediction from measured support');off(ax)
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
    ax=panel(fig,[.07,.09,.40,.33],'b','Separate response components');off(ax)
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
    ax=panel(fig,[.56,.09,.39,.33],'c','Learn an explicit partner increment');off(ax)
    box(ax,.04,.70,.91,.19,'Symmetric single-response kernel\n'+r'$K_{gh,uv}=\frac{1}{2}(L_{gu}L_{hv}+L_{gv}L_{hu})$',fc='#EAF2F8',size=8.4)
    arrow(ax,(.50,.67),(.50,.58),BLUE)
    box(ax,.04,.39,.91,.17,'Remove training endpoint directions\nFit ridge correction  '+r'$G$',fc=PALE,size=8.1)
    arrow(ax,(.50,.36),(.50,.28),TEAL)
    box(ax,.04,.04,.91,.22,'Training-half readout shrinkage\n'+r'$\widehat{Y}_{R}=M_0+G\,\mathrm{diag}(w)$',fc='#E5F2F0',size=9)
    save(fig,1,'Method',['Main Methods: released core, endpoint elimination, RNA kernel','SI Notes 1, 3-6'])

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
    fig=plt.figure(figsize=(7.15,7.0));lb=read(P/'released_benchmark/formal_top100_leaderboard.csv');n=len(lb);ys=np.arange(n)
    ax=panel(fig,[.23,.43,.31,.48],'a','Formal aggregate rank')
    for i,r in enumerate(lb):
        val=float(r['mean_operation_rank']);col=TEAL if i==0 else GRAY
        if i==0:ax.axhspan(-.46,.46,color='#EAF5F3',zorder=0)
        ax.hlines(i,0,val,color=col,lw=1 if i else 2);ax.plot(val,i,'o',color=col,ms=4 if i else 5.5)
        ax.text(val+.35,i,f'{val:.3f}',va='center',fontsize=7,color=col)
    ax.set_yticks(ys,[r['method'].replace('scouter','Scouter') for r in lb]);ax.invert_yaxis();ax.set_ylim(n-.5,-.5);ax.set_xlim(0,16.2);ax.set_xlabel('Mean operation rank (lower is better)');ax.set_xticks([0,4,8,12,16]);style(ax)
    ax.get_yticklabels()[0].set_fontweight('bold');ax.get_yticklabels()[0].set_color(TEAL)
    bx=panel(fig,[.64,.43,.30,.48],'b','Dataset-specific positions')
    ds=['Norman','Wessels','Schmidt','Replogle_exp6'];rr=read(P/'released_benchmark/dataset_panel_leaderboard.csv');index={(r['method'],r['DataSet']):float(r['dataset_panel_position']) for r in rr if r['DEG']=='100'}
    vv=np.array([[index[r['method'],d] for d in ds] for r in lb]);cmap=colors.LinearSegmentedColormap.from_list('rank',['#C5E5E1','#F9FBFB','#F0DFCE'])
    bx.imshow(vv,vmin=1,vmax=16,aspect='auto',cmap=cmap)
    for i in range(n):
        for j in range(4):bx.text(j,i,str(int(vv[i,j])),ha='center',va='center',fontsize=7.0,fontweight='bold' if i==0 else 'normal')
    bx.set_xticks(range(4),['Norman','Wessels','Schmidt','Replogle'],rotation=35,ha='right');bx.set_yticks([]);bx.tick_params(length=0);[s.set_visible(False) for s in bx.spines.values()]
    ax=panel(fig,[.12,.09,.37,.22],'c','Metric-specific positions')
    mm=[r for r in read(P/'released_benchmark/official_top100_per_metric_leaderboard.csv') if r['method']=='WitnessCell'];metricorder=['PCC-delta','MSE','E-distance','SymKL','Wasserstein','Common-DEG']
    # Use source names directly; display aliases only.
    aliases={'PCC-delta':'PCC change','MSE':'MSE','E-distance':'Energy','SymKL':'Sym. KL','Wasserstein':'Wasserstein','Common-DEG':'Common DEG'}
    for i,r in enumerate(mm):
        name=r['metric'];v=float(r['position']);col=TEAL if v==1 else GOLD
        ax.hlines(i,0,v,color=col,lw=2);ax.plot(v,i,'o',color=col,ms=4);ax.text(v+.3,i,str(int(v)),va='center',fontsize=7.3)
    ax.set_yticks(range(6),[aliases.get(r['metric'],r['metric']) for r in mm]);ax.set_ylim(5.65,-.65);ax.set_xlim(0,11);ax.set_xticks([1,5,10]);ax.set_xlabel('Position among 16 entries');style(ax)
    ax=panel(fig,[.64,.09,.30,.22],'d','Combination support in Wessels')
    rr=[r for r in read(P/'Fig3_ED4_witness_gains.csv') if r['reference']=='factorized_prediction']
    for metric,col,label in [('loss_mse',TEAL,'MSE'),('loss_wmse',BLUE,'Weighted MSE')]:
        sub=sorted([r for r in rr if r['metric']==metric],key=lambda r:float(r['fraction']));x=np.array([float(r['fraction']) for r in sub])*100;y=np.array([float(r['gain']) for r in sub])*100
        lo=np.array([float(r['gain_ci95_low']) for r in sub])*100;hi=np.array([float(r['gain_ci95_high']) for r in sub])*100
        ax.plot(x,y,'o-',color=col,ms=4,lw=1.3,label=label);ax.fill_between(x,lo,hi,color=col,alpha=.10,lw=0)
    ax.set_xticks([25,50,100]);ax.set_xlim(20,105);ax.set_ylim(0,53);ax.set_yticks([0,20,40]);ax.set_xlabel('Residual library retained (%)');ax.set_ylabel('Error reduction (%)');style(ax,'y');ax.legend(loc='upper left',frameon=False,fontsize=7)
    save(fig,3,'Prediction_Benchmark',['released_benchmark/*.csv','Fig3_ED4_witness_gains.csv'])

def fig4():
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
    ax=panel(fig,[.65,.39,.29,.21],'d','Readout-specific shrinkage control')
    perms=np.array([float(r['contrast_gain'])*100 for r in read(P/'RNA_weight_shuffles.csv')]);ax.hist(perms,bins=20,color='#C6D6DC',edgecolor='white',linewidth=.4)
    uniform=100*json.loads((P/'RNA_controls.json').read_text())['C_gain_over_M0']['scalar'];gain=100*float(idx['R']['C_gain_vs_M0'])
    ax.axvline(uniform,color=BLUE,lw=1.3,ls='--');ax.axvline(gain,color=TEAL,lw=1.8)
    ax.text(uniform-.3,ax.get_ylim()[1]*.88,'Uniform',rotation=90,va='top',ha='right',color=BLUE,fontsize=7)
    ax.text(gain-.3,ax.get_ylim()[1]*.88,'R',ha='right',color=TEAL,fontweight='bold',fontsize=8)
    ax.set_xlabel('Rectangle-error reduction (%)');ax.set_ylabel('Weight permutations');ax.set_xlim(min(-2,perms.min()-1),18);style(ax,'y')
    ax=panel(fig,[.13,.075,.34,.21],'e','Gains concentrate in a few readouts')
    conc=[]
    for m,col in [('M1',BLUE),('R',TEAL)]:
        delta=(e0-(dat[m]-truth)**2).sum(0);positive=np.sort(np.maximum(delta,0))[::-1];cs=np.cumsum(positive)/positive.sum()*100
        ax.plot(np.arange(1,len(cs)+1),cs,color=col,lw=1.6,label=m);ax.plot(10,cs[9],'o',color=col,ms=4)
        conc.append({'model':m,'top10_positive_gain_share':cs[9]/100})
    ax.axvline(10,color=GRAY,ls=':',lw=.8);ax.set_xscale('log');ax.set_xlim(1,4995);ax.set_ylim(0,104);ax.set_xticks([1,10,100,1000,4995],['1','10','100','1,000','4,995']);ax.minorticks_off()
    ax.set_xlabel('Readouts ranked by positive gain');ax.set_ylabel('Cumulative positive gain (%)');style(ax,'y');ax.legend(frameon=False,loc='lower right',fontsize=7)
    ax.text(14,72,'R: 92.13% in top 10',color=TEAL,fontsize=7.3)
    dump('Fig4_positive_gain_concentration.csv',conc)
    ax=panel(fig,[.65,.075,.29,.21],'f','Readout-scale sensitivity')
    xx=np.arange(3);keys=['C_gain_vs_M0','standardized_C_gain_vs_M0','drop10_C_gain_vs_M0']
    for m,col,offset in [('M1',BLUE,-.15),('R',TEAL,.15)]:
        vals=[100*float(idx[m][k]) for k in keys];ax.bar(xx+offset,vals,width=.27,color=col,label=m,zorder=3)
    ax.axhline(0,color=GRAY,lw=.8);ax.set_xticks(xx,['Raw','Across-condition\nstandardized','Ten-readout\ndeletion']);ax.set_ylabel('Rectangle-error reduction (%)');ax.set_ylim(-14,20);ax.set_yticks([-10,0,10,20]);style(ax,'y');ax.legend(frameon=False,fontsize=7,loc='upper right')
    save(fig,4,'RNA_Partner_Response',['Step2_Fig5_Matched55.csv','Step2_Fig5_SeedMetrics.csv','RNA_signed_rectangle_predictions.npz','RNA_weight_shuffles.csv','RNA_controls.json'])

def growth_gain(stage,budget,model,key='sse',repeat=None):
    rr=[r for r in read(P/'Horlbeck_figure_metrics.csv') if r['stage']==stage and r['budget']==budget and r['fold']!='pooled' and r['repeat']!='pooled' and (repeat is None or r['repeat']==str(repeat))]
    s=lambda m:sum(float(r[key]) for r in rr if r['model']==m)
    return 100*(1-s(model)/s('M0'))

def fig5():
    fig=plt.figure(figsize=(7.15,7.35))
    ax=panel(fig,[.07,.76,.88,.16],'a','Functional prediction across support budgets and cell contexts');off(ax)
    box(ax,.00,.14,.26,.69,'K562 development\n80 genes | 632 held-out pairs\nLow / high support',fc=PALE,size=8.4)
    box(ax,.39,.14,.23,.69,'Freeze the procedure\nRefit coefficients\nin Jurkat',fc='#EAF2F8',size=8.2)
    box(ax,.75,.14,.24,.69,'Jurkat evaluation\n79 genes | 617 held-out pairs\nTwo measured repeats',fc='#E5F2F0',size=8.2)
    arrow(ax,(.28,.48),(.37,.48),BLUE);arrow(ax,(.64,.48),(.73,.48),TEAL)
    ax=panel(fig,[.13,.44,.34,.22],'b','K562: added partner structure')
    for m,col,l in [('rank1_selected',TEAL,'Rank one'),('spectral_selected',BLUE,'Spectral')]:
        for key,ls,mk in [('sse','-','o'),('projected_sse','--','s')]:
            vals=[growth_gain('K562',b,m,key) for b in ['low','high']]
            ax.plot([0,1],vals,ls=ls,marker=mk,color=col,ms=4,lw=1.3,mfc=col if key=='sse' else 'white')
    ax.axhline(0,color=GRAY,lw=.7);ax.set_xticks([0,1],['421-422','2,317-2,318']);ax.set_xlim(-.12,1.12);ax.set_ylim(-2,15);ax.set_ylabel('Error reduction vs same-budget M0 (%)');ax.set_xlabel('Training pairs per fold');style(ax,'y')
    handles=[Line2D([0],[0],color=TEAL,lw=1.3,label='Rank one'),Line2D([0],[0],color=BLUE,lw=1.3,label='Spectral'),Line2D([0],[0],color=GRAY,lw=1.3,label='Complete'),Line2D([0],[0],color=GRAY,lw=1.3,ls='--',label='Partner')]
    ax.legend(handles=handles,frameon=False,ncol=2,fontsize=6.9,loc='upper left',columnspacing=.8,handlelength=1.4)
    ax=panel(fig,[.65,.44,.29,.22],'c','Jurkat: held-out functional gain')
    ints={r['model']:r for r in read(P/'Horlbeck_Jurkat_intervals.csv')}
    for y,m,col in [(1,'rank1_selected',TEAL),(0,'spectral_selected',BLUE)]:
        v=growth_gain('Jurkat','high',m);lo=100*float(ints[m]['lower']);hi=100*float(ints[m]['upper'])
        ax.errorbar(v,y,xerr=[[v-lo],[hi-v]],fmt='o',color=col,ms=5.5,lw=1.3,capsize=3)
        ax.text(v,y+.25,f'{v:.2f}%',ha='center',color=col,fontweight='bold',fontsize=8.5)
    ax.axvline(0,color=GRAY,ls='--',lw=.8);ax.set_yticks([1,0],['Rank one','Spectral']);ax.set_ylim(-.6,1.6);ax.set_xlim(-.5,15);ax.set_xlabel('Complete-error reduction (%)');style(ax)
    ax=panel(fig,[.13,.09,.34,.23],'d','Jurkat: repeat and partner outcomes')
    rows=[]
    for i,(key,label,repeat) in enumerate([('sse','Repeat 1',1),('sse','Repeat 2',2),('projected_sse','Partner contrast',None)]):
        for m,col,shift in [('rank1_selected',TEAL,-.09),('spectral_selected',BLUE,.09)]:
            v=growth_gain('Jurkat','high',m,key,repeat);ax.plot(v,i+shift,'o',ms=5,color=col);ax.text(v+.28,i+shift,f'{v:.2f}',va='center',color=col,fontsize=7);rows.append({'metric':label,'model':m,'gain_percent':v})
    ax.set_yticks(range(3),['Repeat 1','Repeat 2','Partner contrast']);ax.set_ylim(2.5,-.5);ax.set_xlim(0,11);ax.set_xlabel('Error reduction (%)');style(ax)
    dump('Fig5_growth_gains.csv',rows)
    ax=panel(fig,[.66,.09,.28,.23],'e','Partner-dependent corrections')
    # Mean over the two source repeats; both pairs share one original-unit axis.
    ex=read(P/'biological_examples/jurkat/Jurkat_selected_pairs_two_repeats.csv')
    for j,pair in enumerate(['TTI2|PGK1','TTI2|OGFOD1']):
        rr=[r for r in ex if {r['gene_a'],r['gene_b']}==set(pair.split('|'))]
        if not rr:raise ValueError('Missing named example '+pair)
        y=1-j;ax.text(-.35,y+.28,pair.replace('|',' + '),fontsize=8.0,fontweight='bold')
        vals={k:np.mean([float(r[k]) for r in rr]) for k in ['observed','M0','rank1_selected']}
        ax.hlines(y,-.36,0,color='#E3EAEE',lw=.8)
        for key,col,marker in [('observed',INK,'D'),('M0',GRAY,'o'),('rank1_selected',TEAL,'o')]:
            ax.plot(vals[key],y,marker=marker,ms=5,color=col,label={'observed':'Observed','M0':'M0','rank1_selected':'Rank one'}[key] if j==0 else None)
        ax.annotate('',xy=(vals['rank1_selected'],y+.13),xytext=(vals['M0'],y+.13),arrowprops=dict(arrowstyle='->',color=TEAL,lw=1,mutation_scale=8))
    ax.set_xlim(-.36,.005);ax.set_ylim(-.58,1.55);ax.set_yticks([]);ax.set_xticks([-.3,-.2,-.1,0]);ax.set_xlabel('Mean growth response');style(ax,'x')
    ax.legend(frameon=False,fontsize=6.5,loc='lower left',ncol=3,handletextpad=.25,columnspacing=.5,borderaxespad=0,handlelength=.9)
    save(fig,5,'Functional_Growth',['Horlbeck_figure_metrics.csv','Horlbeck_Jurkat_intervals.csv','biological_examples/jurkat/Jurkat_selected_pairs_two_repeats.csv'])

def fig6():
    fig=plt.figure(figsize=(7.15,7.45))
    ax=panel(fig,[.07,.77,.88,.15],'a','PairPick: from a response objective to a candidate partner');off(ax)
    box(ax,.00,.18,.24,.63,'Anchor single response\nDefine objective  '+r'$v_a$',fc='#EAF2F8',size=8.3)
    box(ax,.36,.18,.27,.63,'Predict candidate responses\n'+r'$\widehat U(a,b)=v_a^\top\widehat y_{ab}$',fc='#E5F2F0',size=8.4)
    box(ax,.75,.18,.24,.63,'Select a partner\nEvaluate held-out utility',fc=PALE,size=8.2)
    arrow(ax,(.26,.49),(.34,.49),TEAL);arrow(ax,(.65,.49),(.73,.49),TEAL)
    ax=panel(fig,[.13,.47,.32,.19],'b','Separate risk-selection evaluation')
    risk=[{'loss':'MSE','ratio':.5343,'low':.3646,'high':.7591,'records':213,'coverage':49.48},
          {'loss':'Weighted MSE','ratio':.4775,'low':.3175,'high':.7071,'records':207,'coverage':48.95}]
    for y,r,col in zip([1,0],risk,[TEAL,BLUE]):
        v=r['ratio'];ax.errorbar(v,y,xerr=[[v-r['low']],[r['high']-v]],fmt='o',color=col,ms=5.3,capsize=3,lw=1.25)
        ax.text(v,y+.27,f'{v:.3f}',ha='center',fontsize=8,color=col,fontweight='bold')
    ax.axvline(1,color=GRAY,ls='--',lw=.8);ax.set_xlim(.2,1.08);ax.set_ylim(-.55,1.6);ax.set_yticks([1,0],['MSE','Weighted MSE']);ax.set_xticks([.25,.5,.75,1]);ax.set_xlabel('Accepted / all-query error');style(ax)
    ax.text(.02,-.40,'49.48% / 48.95% pair-balanced coverage',fontsize=6.8,color=GRAY,transform=ax.transAxes)
    dump('Fig6_risk_selection.csv',risk)
    ax=panel(fig,[.59,.44,.36,.24],'c','Prediction and decision endpoints');off(ax)
    rr={r['model']:r for r in read(G/'SUMMARY.csv')};order=['R_G','StandardRidge_RNA','Scouter_Ensemble','GEARS_Ensemble','SingleB_additive'];labels=['PairPick','Linear','Scouter','GEARS','Additive'];keys=['RNA_MSE','partner_MSE','mean_regret'];vv=np.array([[float(rr[m][k])*1000 for k in keys] for m in order])
    for j,label in enumerate(['RNA','Partner','Regret','Hits']):ax.text(.41+j*.165,.97,label,ha='center',fontsize=7.5,fontweight='bold')
    for i,(m,label) in enumerate(zip(order,labels)):
        yy=.79-i*.155;ax.text(-.02,yy,label,va='center',fontsize=7.7,fontweight='bold' if i==0 else 'normal',color=TEAL if i==0 else INK)
        for j in range(3):
            x=.41+j*.165;best=i==int(np.argmin(vv[:,j]));ax.add_patch(Rectangle((x-.075,yy-.063),.15,.125,fc='#C7E6E1' if best else '#F2F5F7',ec='white',lw=.8));ax.text(x,yy,f'{vv[i,j]:.3f}',ha='center',va='center',fontsize=7.4,fontweight='bold' if best else 'normal')
        ax.text(.905,yy,rr[m]['hits']+'/28',ha='center',va='center',fontsize=7.1)
    ax.text(.42,-.09,r'Losses $\times 10^3$; shaded cell = lowest loss',ha='center',fontsize=6.9,color=GRAY)
    ax=panel(fig,[.20,.09,.27,.22],'d','Paired recommendation utility')
    prs={r['reference']:r for r in read(G/'PAIRED_COMPARISONS.csv')};refs=['PRESAGE_Combo_seed24','StandardRidge_RNA','Scouter_Ensemble','GEARS_Ensemble','SingleB_additive'];rl=['PRESAGE','Linear','Scouter','GEARS','Additive']
    for i,m in enumerate(refs):
        r=prs[m];v=float(r['mean_paired_U_difference'])*1000;lo=float(r['cluster_bootstrap_low'])*1000;hi=float(r['cluster_bootstrap_high'])*1000
        ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],fmt='o',color=TEAL,ms=4.5,capsize=2.5,lw=1.1)
    ax.axvline(0,color=GRAY,ls='--',lw=.85);ax.set_yticks(range(5),rl);ax.set_ylim(4.6,-.6);ax.set_xlim(-9,25);ax.set_xticks([-5,0,10,20]);ax.set_xlabel(r'R-G minus reference utility $\times 10^3$');style(ax)
    ax.text(.02,-.39,'95% intervals; eight physical anchors',transform=ax.transAxes,fontsize=6.8,color=GRAY)
    ax=panel(fig,[.66,.09,.28,.22],'e','KLF1 activation: largest pool')
    rr=read(D/'PairPick_largest_pool.csv');x=np.arange(len(rr));obs=np.array([float(r['observed']) for r in rr])*1000;pred=np.array([float(r['R_G']) for r in rr])*1000
    ax.plot(x,obs,'D-',color=INK,ms=3.8,lw=.9,label='Observed');ax.plot(x,pred,'o--',color=TEAL,ms=4,lw=1,label='PairPick')
    best=int(np.argmax(pred));ax.scatter(best,pred[best],s=95,facecolors='none',edgecolors=TEAL,lw=1.1,zorder=4)
    ax.set_xticks(x,[r['partner'] for r in rr],rotation=55,ha='right');ax.set_ylabel(r'Objective utility $\times 10^3$');ax.set_xlim(-.35,len(rr)-.65);style(ax,'y');ax.legend(frameon=False,fontsize=6.8,loc='best')
    save(fig,6,'Research_Applications',['SI Table S16 (risk selection, displayed precision)','gse220974/SUMMARY.csv','gse220974/PAIRED_COMPARISONS.csv','PairPick_largest_pool.csv'])

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--figures',default='1,2,3,4,5,6');args=ap.parse_args()
    for i in map(int,args.figures.split(',')):
        globals()['fig'+str(i)]();print('Rendered figure',i,flush=True)
    auditpath=ROOT/'FIGURE_AUDIT.json'
    previous=json.loads(auditpath.read_text()) if auditpath.exists() else []
    done={a['figure'] for a in AUDIT};AUDIT=sorted([a for a in previous if a['figure'] not in done]+AUDIT,key=lambda a:a['figure'])
    auditpath.write_text(json.dumps(AUDIT,indent=2)+'\n')
    writer=PdfWriter()
    for r in AUDIT:
        reader=PdfReader(OUT/(r['stem']+'.pdf'));writer.add_page(reader.pages[0])
    writer.add_metadata({'/Title':'WitnessCell v24 - Main Figures','/Author':'WitnessCell','/Subject':'Six main figures; legends provided separately'})
    with (ROOT/'WitnessCell_v24_Main_Figures.pdf').open('wb') as h:writer.write(h)
    from PIL import Image,ImageOps,ImageDraw
    thumbs=[]
    for r in AUDIT:
        im=Image.open(OUT/(r['stem']+'_preview.png')).convert('RGB');im.thumbnail((740,835));canvas=Image.new('RGB',(780,895),'white');canvas.paste(im,((780-im.width)//2,40));ImageDraw.Draw(canvas).text((25,12),'FIGURE '+str(r['figure']),fill=INK);thumbs.append(canvas)
    contact=Image.new('RGB',(1560,895*math.ceil(len(thumbs)/2)),'#E6EDF1')
    for i,im in enumerate(thumbs):contact.paste(im,((i%2)*780,(i//2)*895))
    contact.save(ROOT/'Figures_Overview.png')
    print(json.dumps({'outside_page':{r['figure']:r['text_outside_page'] for r in AUDIT},'figures':len(AUDIT)}))
