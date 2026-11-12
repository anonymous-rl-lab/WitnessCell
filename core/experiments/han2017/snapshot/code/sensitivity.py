"""Post-score, identity-defined robustness checks; no model or rule selection."""
from common import *
import argparse
import pandas as pd
import run

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--old-published-predictions',type=Path,required=True);args=ap.parse_args()
    old=load(args.old_published_predictions)['pairs'];smoke={tuple(sorted(p)) for p in old};dump(ROOT/'audit/PRIOR_SCORED_PAIRS.json',sorted(smoke))
    missing=set(load(ROOT/'data/target_knowledge.npz')['missing'].astype(str))
    dump(ROOT/'audit/SENSITIVITY_PLAN.json',dict(created_utc=now(),status='Post-score robustness checks, no model changes or primary-result replacement',filters=['remove_prior_scored_pairs_from_calibration_and_query','remove_missing_knowledge_endpoints_from_calibration_and_query'],prior_scored_pairs=len(smoke),missing_genes=sorted(missing),fixed_seeds=True,fixed_alpha_rule=True,base_model_fits=0))
    rows=[];states=[]
    for panel in ['full_raw','full_published']:
        data=load(ROOT/f'data/{panel}_frozen_predictions.npz');pairs=data['pairs']
        masks={'exclude_smoke':np.array([tuple(sorted(p)) not in smoke for p in pairs]),
               'exclude_missing_knowledge':np.array([not (set(p)&missing) for p in pairs])}
        for filt,keep in masks.items():
            for seed in range(20261108,20261113):
                for fold in range(3):
                    sp=ROOT/'data/splits'/panel/str(seed)/str(fold);tr=load(sp/'calibration.npz');qu=load(sp/'query_outcomes.npz')
                    tm=keep[tr['ids']];qm=keep[qu['ids']];tids=tr['ids'][tm];qids=qu['ids'][qm];y=qu['y'][qm]
                    state=run.calibrate(data,tids,tr['y'][tm],seed,fold);p,arms=run.arms(data,qids,state,panel)
                    states.append(dict(panel=panel,filter=filt,seed=seed,fold=fold,n_train=len(tids),n_query=len(qids),alpha=state['alpha']))
                    for arm,pred in arms.items():
                        e=y[:,:,None]-pred;pe=p.P(e)
                        for j,model in enumerate(data['models']):rows.append(dict(panel=panel,filter=filt,seed=seed,fold=fold,model=model,arm=arm,n=len(qids)*y.shape[1],full_sse=np.sum(e[:,:,j]**2),partner_sse=np.sum(pe[:,:,j]**2)))
            print(json.dumps(dict(stage='sensitivity_complete',panel=panel,filter=filt,remaining_pairs=int(keep.sum()))),flush=True)
    df=pd.DataFrame(rows);df.to_csv(ROOT/'audit/SENSITIVITY_FOLD_SCORES.csv',index=False);dump(ROOT/'audit/SENSITIVITY_COEFFICIENTS.json',states)
    g=df.groupby(['panel','filter','model','arm'],as_index=False)[['n','full_sse','partner_sse']].sum();out=[]
    for (panel,filt,model),d in g.groupby(['panel','filter','model']):
        d=d.set_index('arm')
        for ref in ['original','single_additive']:
            out.append(dict(panel=panel,filter=filt,model=model,reference=ref,full_gain_pct=100*(1-d.loc['target_control','full_sse']/d.loc[ref,'full_sse']),partner_gain_pct=100*(1-d.loc['target_control','partner_sse']/d.loc[ref,'partner_sse'])))
    pd.DataFrame(out).to_csv(ROOT/'audit/SENSITIVITY_SUMMARY.csv',index=False)
if __name__=='__main__':main()
