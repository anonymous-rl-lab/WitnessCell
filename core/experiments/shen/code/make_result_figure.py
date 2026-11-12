#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    score=json.loads((ROOT/'results/formal/SCORES.json').read_text())
    cells=['HeLa','A549','293T']
    full=[100*score['models'][c]['rank1_primary']['full_gain'] for c in cells]
    part=[100*score['models'][c]['rank1_primary']['partner_gain'] for c in cells]
    intervals=[100*np.array(score['conditional_sensitivity']['per_context'][c]['central95']) for c in cells]
    part.append(100*score['mean_within_context_partner_gain'])
    intervals.append(100*np.array(score['conditional_sensitivity']['equal_context_mean']['central95']))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,(ax,bx)=plt.subplots(1,2,figsize=(9.2,3.8),gridspec_kw={'width_ratios':[1,1.25]})
    ax.bar(cells,full,color='#327A93',width=.58)
    for i,v in enumerate(full):ax.text(i,v+.16,f'{v:.2f}%',ha='center',fontsize=11)
    ax.set_ylim(0,7.5);ax.set_ylabel('Error reduction relative to M0 (%)')
    ax.set_title('a  Complete response',loc='left',fontweight='bold')
    ys=np.arange(4)[::-1]
    for i,(v,ci) in enumerate(zip(part,intervals)):
        color='#B45B42' if i<3 else '#323C46'
        bx.plot(ci,[ys[i]]*2,color=color,lw=2)
        bx.scatter([v],[ys[i]],c=color,s=42,zorder=4)
        bx.text(3.1,ys[i],f'{v:+.2f}%',ha='left',va='center',fontsize=10)
    bx.axvline(0,color='#7B858D',ls='--',lw=1)
    bx.set_yticks(ys,cells+['3-context mean']);bx.set_xlim(-4.7,4.3);bx.set_ylim(-.5,3.5)
    bx.set_xticks([-4,-2,0,2]);bx.set_xlabel('Error reduction relative to M0 (%)')
    bx.set_title('b  Endpoint-eliminated partner response',loc='left',fontweight='bold',fontsize=10.5)
    fig.subplots_adjust(left=.085,right=.98,bottom=.22,top=.85,wspace=.40)
    fig.text(.085,.055,'Primary signed rank-one model. Positive values indicate lower error.\nLines: central 95% of fixed-prediction endpoint reweighting; conditional sensitivity only.',fontsize=8.5,color='#4D5760')
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    fig.savefig(out/'Shen2017_Formal_Results.png',dpi=200)
    fig.savefig(out/'Shen2017_Formal_Results.svg')
    plt.close(fig)
if __name__=='__main__':main()
