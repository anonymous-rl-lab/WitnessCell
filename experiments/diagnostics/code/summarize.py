"""Regenerate the two pooled tables used in SI S55 and S57."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument('--results',type=Path,required=True)
a=ap.parse_args(); O=a.results
rows=[]
for (panel,source),g in pd.read_csv(O/'HAN_ORACLE_FOLDS.csv').groupby(['panel','source']):
    E,A,B=g[['E','A','B']].sum()
    actual=100*g.actual_improvement.sum()/E
    oracle=100*g.bounded_improvement.sum()/E
    alpha=float(np.clip(A/B,0,1))
    rows.append(dict(panel=panel,source=source,actual_gain_pct=actual,
        fold_oracle_gain_pct=oracle,pooled_single_oracle_gain_pct=100*(2*alpha*A-alpha*alpha*B)/E,
        actual_fraction_of_fold_oracle_pct=100*actual/oracle))
pd.DataFrame(rows).to_csv(O/'HAN_ORACLE_SUMMARY.csv',index=False)
p=pd.read_csv(O/'POOLED_SCORES.csv'); rows=[]
for panel,g in p[p.policy.eq('baseline')].groupby('panel'):
    base=g.set_index('baseline'); additive=base.loc['additive']
    for name in ['saturation','local4']:
        v=base.loc[name]
        rows.append(dict(panel=panel,baseline=name,
            partner_gain_vs_additive=100*(1-v.partner_sse/additive.partner_sse),
            full_gain_vs_additive=100*(1-v.full_sse/additive.full_sse)))
pd.DataFrame(rows).to_csv(O/'BASELINE_COMPARISON.csv',index=False)
print('Pooled oracle and baseline summaries regenerated.')
