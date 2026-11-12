import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'results'
dec=pd.read_csv(O/'CALIBRATION_DECISIONS.csv').set_index('feature_id');g=pd.read_csv(O/'CALIBRATION_BY_CONDITION.csv');c=pd.read_csv(O/'PARTNER_COMPONENT_BY_CONDITION.csv').set_index('group')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(2,2,figsize=(11,8),constrained_layout=True);grey='#94a3b8';green='#0f766e';orange='#d97706';blue='#2563eb'
cells=['Jurkat','HeLa','A549','293T'];x=np.arange(4)
a=ax[0,0]
for off,col,label,color in [(-.18,'inner_full_gain','Full response',blue),(.18,'inner_partner_gain','Partner response',green)]:a.bar(x+off,[dec.loc[k+'_high',col] for k in cells],.35,color=color,label=label)
a.set_xticks(x,cells);a.axhline(0,color='.6',lw=.8);a.set(title='A  What did the inner validation measure?',ylabel='Inner gain (%)');a.legend(frameon=False,fontsize=8,loc='lower left')
a=ax[0,1];shen=['HeLa','A549','293T'];x=np.arange(3)
for off,policy,color,label in [(-.18,'original',grey,'Original rank-one'),(.18,'inner_shrink',green,'Training-calibrated')]:
 vals=[g[(g.group==k+'_high')&(g.policy==policy)].partner_gain_pp.iloc[0] for k in shen];b=a.bar(x+off,vals,.35,color=color,label=label)
 a.bar_label(b,labels=[f'{v:+.3f}' for v in vals],fontsize=8,padding=3)
a.axhline(0,color='.6',lw=.8);a.set_xticks(x,shen);a.set(title='B  Shen: partner correction after calibration',ylabel='Outer partner gain (%)',ylim=(-1.4,.45));a.legend(frameon=False,fontsize=8,loc='lower right')
a=ax[1,0];keys=['K562_low','K562_layout','K562_high','Jurkat_high'];x=np.arange(4)
for off,policy,color,label in [(-.18,'original',grey,'Original rank-one'),(.18,'inner_shrink',green,'Training-calibrated')]:
 vals=[g[(g.group==k)&(g.policy==policy)].partner_gain_pp.iloc[0] for k in keys];b=a.bar(x+off,vals,.35,color=color,label=label)
 a.bar_label(b,labels=[f'{v:.2f}' for v in vals],fontsize=8,padding=3)
a.axhline(0,color='.6',lw=.8);a.set_xticks(x,['K562 low','K562 layouts*','K562 high','Jurkat']);a.set(title='C  The same rule across functional settings',ylabel='Outer partner gain (%)',ylim=(-1,13));a.legend(frameon=False,fontsize=8,loc='upper left')
a=ax[1,1];x=np.arange(4)
for off,policy,color,label in [(-.25,'original',grey,'Original'),(0,'inner_shrink',orange,'Shrink whole correction'),(.25,'partner_only',green,'Shrink partner part only*')]:
 vals=[c.loc[k+'_high','full_gain_pp'] if policy=='partner_only' else g[(g.group==k+'_high')&(g.policy==policy)].full_gain_pp.iloc[0] for k in cells]
 a.bar(x+off,vals,.24,color=color,label=label)
a.axhline(0,color='.6',lw=.8);a.set_xticks(x,cells);a.set(title='D  Preserve the useful endpoint contribution',ylabel='Outer full-response gain (%)',ylim=(0,8.9));a.legend(frameon=False,fontsize=8,loc='upper right')
fig.savefig(O/'Training_Calibration.png',dpi=220);fig.savefig(O/'Training_Calibration.pdf');plt.close(fig)
(O/'FIGURE_NOTES.txt').write_text('Model_Matched_Gain: all 788 matched records, including terminal checkpoints; panel B is sensitive to their inclusion. See REGRESSION_SENSITIVITY.csv. Points are not independent biological replicates.\nTraining_Calibration: coefficients use training-side inner predictions only; all effects are retrospective. B shows pooled two-repeat gains, which do not imply positivity in each repeat. C: K562 layouts pool 45 overlapping fold/seed/layout configurations, not 45 independent experiments. D: partner-only shrinkage uses the identities of the fixed query batch to project its correction, not query outcomes. It is a secondary, query-batch-dependent transformation. The full-response gains of whole-correction shrinkage and partner-only shrinkage must not be conflated.\n')
