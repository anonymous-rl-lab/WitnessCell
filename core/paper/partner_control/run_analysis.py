#!/usr/bin/env python3
"""Reproduce four policies and two fixed checks from archived arrays."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'software'))
from witnesscell_partner_control import fit_amplitude,predict_panel,project


def loss(y,p,e,n):
    endpoint,partner,dim=project(e,n,y-p)
    return dict(full_sse=float(np.sum((y-p)**2)),endpoint_sse=float(np.sum(endpoint**2)),partner_sse=float(np.sum(partner**2)),contrast_dim=dim)


def gains(s,b):
    return {k+'_gain_pp':100*(1-s[k+'_sse']/b[k+'_sse']) if b[k+'_sse']>1e-14 else np.nan for k in ['full','endpoint','partner']}


def aggregate(df,groups):
    cols=['full_sse','endpoint_sse','partner_sse','base_full_sse','base_endpoint_sse','base_partner_sse']
    out=df.groupby(groups,sort=False)[cols].sum().reset_index()
    for k in ['full','endpoint','partner']:out[k+'_gain_pp']=100*(1-out[k+'_sse']/out['base_'+k+'_sse'])
    return out


def main():
    O=R/'results';O.mkdir(exist_ok=True)
    meta=pd.read_csv(R/'inputs/CONFIGURATIONS.csv').sort_values('feature_id')
    inner=np.load(R/'inputs/INNER_GROWTH_REPLAY.npz');outer=np.load(R/'inputs/CALIBRATION_OUTCOMES.npz')
    decisions=[];repeat_decisions=[];fitted={}
    # Complete every training-side decision before consulting outer truth.
    for m in meta.itertuples():
        key=m.feature_id
        def ii(s):return inner[key+'__'+s]
        base=ii('M0');pred=base if bool(ii('selected_zero')) else ii('regularized')
        args=(ii('y'),base,pred,ii('edges'),ii('folds'),int(ii('nodes')))
        pooled=fit_amplitude(*args);fitted[key]=pooled
        d={k:v for k,v in pooled.items() if k!='folds'}
        d.update(feature_id=key,cell=m.cell,study=m.study,selected_zero=bool(ii('selected_zero')),
                 condition=key.rsplit('_fold',1)[0] if not key.startswith('K562_layout') else 'K562_layout')
        decisions.append(d)
        for j in range(ii('y').shape[1]):
            v=fit_amplitude(ii('y')[:,j],base[:,j],pred[:,j],ii('edges'),ii('folds'),int(ii('nodes')))
            repeat_decisions.append(dict(feature_id=key,condition=d['condition'],training_repeat=j+1,alpha=v['alpha'],rho=v['rho'],inner_gain_pp=v['inner_gain_pp'],cross=v['cross'],correction_energy=v['correction_energy']))
    dd=pd.DataFrame(decisions);rd=pd.DataFrame(repeat_decisions)
    dd.to_csv(O/'DECISIONS.csv',index=False);rd.to_csv(O/'REPEAT_DECISIONS.csv',index=False)
    digest=hashlib.sha256((O/'DECISIONS.csv').read_bytes()).hexdigest()
    rows=[];repeats=[];panels=[];checks=[];pred_archive={}
    for m in meta.itertuples():
        key=m.feature_id
        def oo(s):return outer[key+'__'+s]
        y,b,p,e,n=oo('truth'),oo('M0'),oo('original'),oo('edges'),int(oo('nodes'))
        alpha=fitted[key]['alpha'];condition=dd.loc[dd.feature_id.eq(key),'condition'].iloc[0]
        result=predict_panel(b,p,e,n,alpha);base=loss(y,b,e,n)
        original=loss(y,p,e,n);common=dict(feature_id=key,condition=condition,study=m.study,cell=m.cell,alpha=alpha,query_pairs=len(e),**{'base_'+k:v for k,v in base.items() if k.endswith('sse')})
        score={}
        for policy in ['baseline','original','whole_shrink','partner_only']:
            ss=loss(y,result[policy],e,n);score[policy]=ss
            rows.append(dict(**common,policy=policy,**ss,**gains(ss,base)))
            pred_archive[key+'__'+policy]=result[policy]
        pred_archive[key+'__edges']=e
        de=result['endpoint_correction'];dp=result['partner_correction']
        controlled_de,_,_=project(e,n,result['partner_only']-b)
        checks.append(dict(feature_id=key,
            decomposition_error=max(abs(s['full_sse']-s['endpoint_sse']-s['partner_sse']) for s in score.values()),
            endpoint_preservation_error=float(np.max(abs(controlled_de-de))),
            partner_equivalence_error=abs(score['whole_shrink']['partner_sse']-score['partner_only']['partner_sse']),
            sse_improvement_identity_error=abs((original['full_sse']-score['partner_only']['full_sse'])-(original['partner_sse']-score['partner_only']['partner_sse'])),
            original_archive_match_error=float(np.max(abs(result['original']-p)))))
        alpha_variants={'pooled':alpha}
        alpha_variants.update({'repeat_'+str(int(t.training_repeat)):float(t.alpha) for t in rd[rd.feature_id.eq(key)].itertuples()})
        for rule,aa in alpha_variants.items():
            rr=predict_panel(b,p,e,n,aa)
            for j in range(y.shape[1]):
                bb=loss(y[:,j],b[:,j],e,n)
                for pol in ['original','whole_shrink','partner_only']:
                    ss=loss(y[:,j],rr[pol][:,j],e,n)
                    repeats.append(dict(feature_id=key,condition=condition,training_alpha=rule,alpha=aa,evaluation_repeat=j+1,policy=pol,**ss,**{'base_'+k:v for k,v in bb.items() if k.endswith('sse')},**gains(ss,bb)))
        canonical=np.sort(e,axis=1);order=np.lexsort((canonical[:,1],canonical[:,0]));rank=np.empty(len(e),int);rank[order]=np.arange(len(e))
        for parts,label in [(1,'full'),(2,'identity_halves'),(3,'identity_thirds'),(len(e),'singletons')]:
            pp=np.empty_like(p);dims=0
            if label=='singletons':pp=p.copy()
            else:
                for block in range(parts):
                    mask=rank%parts==block
                    out=predict_panel(b[mask],p[mask],e[mask],n,alpha)
                    pp[mask]=out['partner_only'];dims+=int(out['contrast_dim'])
            ss=loss(y,pp,e,n)
            panels.append(dict(**common,panel_rule=label,blocks=parts,controllable_dimensions=dims,
                rms_prediction_change_vs_full=float(np.sqrt(np.mean((pp-result['partner_only'])**2))),
                max_prediction_change_vs_full=float(np.max(abs(pp-result['partner_only']))),
                **ss,**gains(ss,base)))
        print('completed',key,flush=True)
    scores=pd.DataFrame(rows);rep=pd.DataFrame(repeats);panel=pd.DataFrame(panels);cc=pd.DataFrame(checks)
    scores.to_csv(O/'FOUR_ARM_SCORES.csv',index=False)
    aggregate(scores,['condition','policy']).to_csv(O/'FOUR_ARM_BY_CONDITION.csv',index=False)
    rep.to_csv(O/'REPEAT_SCORES.csv',index=False)
    aggregate(rep,['condition','training_alpha','evaluation_repeat','policy']).to_csv(O/'REPEAT_BY_CONDITION.csv',index=False)
    panel.to_csv(O/'QUERY_PANEL_SCORES.csv',index=False)
    aggregate(panel,['condition','panel_rule']).to_csv(O/'QUERY_PANEL_BY_CONDITION.csv',index=False)
    cc.to_csv(O/'ALGEBRA_CHECKS.csv',index=False)
    np.savez_compressed(O/'FOUR_ARM_PREDICTIONS.npz',**pred_archive)
    old=pd.read_csv(R/'inputs/ARCHIVED_CALIBRATION_DECISIONS.csv')
    match=dd.merge(old[['feature_id','alpha']],on='feature_id',suffixes=('_new','_old'))
    alpha_difference=float(max(abs(match.alpha_new-match.alpha_old)))
    maxerr={c:float(cc[c].max()) for c in cc.columns if c!='feature_id'}
    assert max(maxerr.values())<1e-8
    assert alpha_difference<1e-9
    assert hashlib.sha256((O/'DECISIONS.csv').read_bytes()).hexdigest()==digest
    metrics=scores.pivot(index='feature_id',columns='policy',values='partner_sse')
    diff=metrics['original']-metrics['partner_only'];nonzero=~dd.set_index('feature_id').selected_zero
    counts=dict(improved=int((diff>1e-12).sum()),unchanged=int((abs(diff)<=1e-12).sum()),worsened=int((diff< -1e-12).sum()))
    repeat_wide=rd.pivot(index='feature_id',columns='training_repeat',values='alpha')
    deltas=abs(repeat_wide[1]-repeat_wide[2])
    summary=dict(configurations=len(meta),selected_zero=int(dd.selected_zero.sum()),nonzero=int(nonzero.sum()),four_arm_rows=len(scores),repeat_score_rows=len(rep),panel_score_rows=len(panel),
        alpha_replay_max_difference=alpha_difference,algebra_max_errors=maxerr,decisions_sha256=digest,
        component_vs_original_configuration_counts=counts,repeat_alpha_median_absolute_difference=float(deltas.median()),repeat_alpha_max_absolute_difference=float(deltas.max()),
        interpretation='Configuration counts are descriptive; K562 layouts share queries. No new response-model training and no new independent confirmation.',
        original_outcomes_preserved=True)
    (O/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    # A runnable real-data example, with training and query inputs physically separated.
    key='HeLa_high';sel=bool(inner[key+'__selected_zero']);eb=inner[key+'__M0']
    ex=R/'example';ex.mkdir(exist_ok=True)
    np.savez_compressed(ex/'training_oof.npz',y=inner[key+'__y'],baseline=eb,prediction=eb if sel else inner[key+'__regularized'],edges=inner[key+'__edges'],folds=inner[key+'__folds'],nodes=inner[key+'__nodes'])
    np.savez_compressed(ex/'query_predictions.npz',baseline=outer[key+'__M0'],prediction=outer[key+'__original'],edges=outer[key+'__edges'],nodes=outer[key+'__nodes'])
    np.savez_compressed(ex/'evaluation_truth.npz',truth=outer[key+'__truth'])
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
