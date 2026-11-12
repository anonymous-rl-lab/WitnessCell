from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
df=pd.read_csv(root/'evaluation/SCORES.csv')
support=json.loads((root/'diagnostics/PRIOR_TRAINING_DIAGNOSTICS.json').read_text())['support']
cells=['HeLa','A549','293T','Jurkat'];s={r['dataset']:r for r in support};x=np.arange(4)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axes=plt.subplots(1,3,figsize=(12,3.5),layout='constrained')
ax=axes[0];ax.bar(x-.17,[s[c]['full_repeat_correlation'] for c in cells],.32,color='#a7bac7',label='Full response')
ax.bar(x+.17,[s[c]['partner_repeat_correlation'] for c in cells],.32,color='#196b8b',label='Endpoint removed')
ax.set(xticks=x,xticklabels=cells,ylim=(0,1),ylabel='Training-repeat correlation',title='a  Similar coverage, different repeatability');ax.legend(frameon=False,fontsize=8)
ax=axes[1]
methods=[('original_rank1','Original rank one','#8f8f8f'),('sf12','Recovered SF12','#b57a28'),('network_v4','Trainable V4 adapter','#679c91'),('combined_selected','Training-selected update','#196b8b')]
for k,(method,label,color) in enumerate(methods):
 vals=[df[(df.cell==c)&(df.method==method)].partner_gain_pct.iloc[0] for c in cells]
 ax.bar(x+(k-1.5)*.17,vals,.16,label=label,color=color)
ax.axhline(0,color='black',lw=.7);ax.set(xticks=x,xticklabels=cells,ylabel='Partner MSE reduction (%)',title='b  Original query sets and reference M0');ax.legend(frameon=False,fontsize=7,loc='upper left')
ax=axes[2]
for k,(m,label,color) in enumerate([('network_without_calibration','V4 before calibration','#b77e79'),('network_v4','V4 after calibration','#679c91')]):
 vals=[df[(df.cell==c)&(df.method==m)].partner_gain_pct.iloc[0] for c in cells]
 ax.bar(x+(k-.5)*.35,vals,.32,label=label,color=color)
ax.axhline(0,color='black',lw=.7);ax.set(xticks=x,xticklabels=cells,ylabel='Partner MSE reduction (%)',title='c  Correction amplitude matters');ax.legend(frameon=False,fontsize=8,loc='upper left')
out=root/'figures';out.mkdir(exist_ok=True)
fig.savefig(out/'Partner_Recovery_Diagnostics.png',dpi=220);fig.savefig(out/'Partner_Recovery_Diagnostics.pdf')
