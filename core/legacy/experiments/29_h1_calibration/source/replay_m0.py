"""Exactly replay selection and full-gene scores from compact sufficient stats.
No source expression, numpy only. Includes independent selector and Gram scoring.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import csv,json,numpy as np
HERE=Path(__file__).resolve().parent
def main():
    rows=json.loads((HERE/'SCORE_GRAM_ROWS.json').read_text());K=np.load(HERE/'SCORE_GRAM.npz')['gram']
    truth={(r['global_id'],r['half']):i for i,r in enumerate(rows) if r['type']=='reference'}
    table=list(csv.DictReader((HERE/'PER_CONDITION_HALF_ERRORS.csv').open()))
    saved={(int(r['fold']),r['arm'],r['condition'],int(r['half'])):float(r['mse']) for r in table}
    reproduced={};maxdiff=0.
    for i,r in enumerate(rows):
      if r['type']=='prediction':
        for h in (0,1):
          j=truth[r['global_id'],h];v=float(K[i,i]+K[j,j]-2*K[i,j])
          key=(r['fold'],r['arm'],r['condition'],h)
          maxdiff=max(maxdiff,abs(v-saved[key]));reproduced[key]=v
    stats=np.load(HERE/'SELECTION_QUADRATICS.npz');selection=json.loads((HERE/'SELECTION.json').read_text())
    proto=json.loads((HERE/'PROTOCOL_BEFORE_FIT.json').read_text());etas=proto['baseline']['noises'];alphas=proto['baseline']['alphas']
    maxgamma=0.;maxrisk=0.
    for s in selection:
        p=f"fold{s['fold']}__{s['arm']}";t=stats[p+'__truth2'];dot=stats[p+'__dot'];r2=stats[p+'__raw2']
        allc=[]
        for ai in range(len(alphas)):
          for ni,eta in enumerate(etas):
            g=min(1.,max(0.,float(dot[ai,ni]/r2[ai,ni]))) if r2[ai,ni]>1e-12 else 0.
            risk=float(t[ai]-2*g*dot[ai,ni]+g*g*r2[ai,ni]);allc.append((risk,-eta,ai,ni,g))
        win=min(allc)
        assert [win[2],win[3]]==s['index']
        maxgamma=max(maxgamma,abs(win[4]-s['selected']['gamma']))
        maxrisk=max(maxrisk,abs(win[0]/(20*37292)-s['oof_mse']))
    assert maxdiff<1e-12 and maxgamma<1e-12 and maxrisk<1e-12
    arms=proto['arms'];conditions=sorted({k[2] for k in reproduced})
    errors={arm:np.array([np.mean([v for k,v in reproduced.items() if k[1]==arm and k[2]==c]) for c in conditions]) for arm in arms}
    diagnostics={}
    for arm in arms[1:]:
        delta=errors['raw']-errors[arm];loo=[]
        for j in range(len(conditions)):
            idx=np.arange(len(conditions))!=j
            loo.append(100*(1-errors[arm][idx].sum()/errors['raw'][idx].sum()))
        order=np.argsort(-delta)
        diagnostics[arm]={'improved_conditions_vs_raw':int(np.sum(delta>0)),'total_conditions':len(conditions),
                         'leave_one_condition_out_gain_pct_range':[min(loo),max(loo)],
                         'largest_positive_contributors':[{'condition':conditions[j],'absolute_mse_reduction':float(delta[j])} for j in order[:3]],
                         'half_gains_pct':{str(h):100*(1-np.mean([v for k,v in reproduced.items() if k[1]==arm and k[3]==h])/np.mean([v for k,v in reproduced.items() if k[1]=='raw' and k[3]==h])) for h in(0,1)}}
    d=errors['B1_direct7']-errors['B2_cal9'];ds=errors['B1_direct7']-errors['B2_cal7']
    result={'score_max_abs_error':maxdiff,'selection_gamma_max_abs':maxgamma,'selection_mse_max_abs':maxrisk,
            'candidate_selections_verified':len(selection),'condition_half_scores_verified':len(saved),
            'cal9_conditions_better_than_B1':int(np.sum(d>0)),'cal7_conditions_better_than_B1':int(np.sum(ds>0)),
            'per_arm_diagnostics':diagnostics}
    (HERE/'REPLAY_AND_INFLUENCE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
