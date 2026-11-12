#!/usr/bin/env python3
"""Four policies from exact squared-error sufficient statistics; no retraining."""
from pathlib import Path
import pandas as pd,numpy as np,json
R=Path(__file__).resolve().parent
d=pd.read_csv(R/'inputs/RNA_SUFFICIENT_STATISTICS.csv');rows=[]
for r in d.itertuples():
    for policy,ae,ap in [('baseline',0.,0.),('original',1.,1.),('whole_shrink',r.alpha,r.alpha),('partner_only',1.,r.alpha)]:
        le=r.outer_endpoint_E-2*ae*r.outer_endpoint_A+ae*ae*r.outer_endpoint_B
        lp=r.outer_partner_E-2*ap*r.outer_partner_A+ap*ap*r.outer_partner_B
        rows.append(dict(record_id=r.record_id,feature_id=r.feature_id,model=r.model,pool_id=r.pool_id,actual_family=r.actual_family,chosen_eta=r.chosen_eta,alpha=r.alpha,policy=policy,endpoint_sse=le,partner_sse=lp,full_sse=le+lp,base_endpoint_sse=r.outer_endpoint_E,base_partner_sse=r.outer_partner_E,base_full_sse=r.outer_endpoint_E+r.outer_partner_E,
            full_gain_pp=100*(1-(le+lp)/(r.outer_endpoint_E+r.outer_partner_E)),partner_gain_pp=100*(1-lp/r.outer_partner_E)))
out=pd.DataFrame(rows);out.to_csv(R/'results/RNA_FOUR_ARM_SCORES.csv',index=False)
summary=dict(matched_records=len(d),score_rows=len(out),training_folds=d.feature_id.nunique(),source_snapshot_pools=d.pool_id.nunique(),calibrated_nontrivially=int((d.alpha<1-1e-12).sum()),
    inner_replay_max_error=float(d.inner_replay_error.max()),outer_replay_max_error=float(d.outer_replay_error.max()),
    scope='Retrospective saved-candidate analysis on global query projection. Archive snapshots and candidate aliases share data; not independent experiments. Original eleven-rectangle R result is unchanged.')
(R/'results/RNA_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
