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
f.box(ax,.03,.69,.94,.25,'RNA: projected pair kernel\nTraining-half readout weights',fc='#EAF2F8',size=8.0)
f.box(ax,.03,.38,.94,.25,'Function: indexed residual matrix\nRank one / spectral correction',fc='#E5F2F0',size=8.0)
f.box(ax,.03,.07,.94,.25,'Function: shared SF12 / V4\nSingle effects + gene annotations',fc='#F1ECF6',size=8.0)
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


import json
(R/'revision/FIG1_AUDIT_V36.json').write_text(json.dumps(f.AUDIT,indent=2)+'\n')
