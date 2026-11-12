"""Recompute v40 tables from saved predictions. No predictor fitting here."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parents[1]

def projector(edges,n):
    x=np.zeros((len(edges),n));x[np.arange(len(edges))[:,None],edges]=1
    u,s,_=np.linalg.svd(x,full_matrices=False)
    q=u[:,s>np.finfo(float).eps*max(x.shape)*s[0]]
    return lambda a:a-q@(q.T@a)

def common_reference():
    d=pd.read_csv(R/'inputs/V39_MATCHED_REFERENCE_POOLED.csv');out=[]
    for group,a in d.groupby('group'):
        native=a[a.method=='native_B'].iloc[0]
        for row in a.itertuples():
            out.append(dict(group=group,method=row.method,full_sse=row.full_sse,
                partner_sse=row.partner_sse,original_full_sse=native.full_sse,
                original_partner_sse=native.partner_sse,
                full_gain_original_pct=100*(1-row.full_sse/native.full_sse),
                partner_gain_original_pct=100*(1-row.partner_sse/native.partner_sse)))
    pd.DataFrame(out).to_csv(R/'results/COMMON_REFERENCE.csv',index=False)

def saturation():
    labels=json.loads((R/'ANALYSIS_PLAN.json').read_text())['configurations']; rows=[];checks=[]
    # Check every prediction lock before any query outcome is opened by this pass.
    for label in labels:
        lock=json.loads((R/'fitted'/label/'LOCK.json').read_text())
        actual=hashlib.sha256((R/'fitted'/label/'PREDICTIONS.npz').read_bytes()).hexdigest()
        assert actual==lock['predictions_sha256']
    for label in labels:
        tr=np.load(R/'inputs'/label/'training.npz'); p=np.load(R/'fitted'/label/'PREDICTIONS.npz')
        old=np.load(R/'inputs'/label/'original.npz');out=np.load(R/'inputs'/label/'outcomes.npz')
        assert np.array_equal(p['query_edges'],out['edges'])
        y=out['y'];proj=projector(out['edges'],len(tr['genes']));e0=y-p['B'];old_e=y-old['M0']
        preds={'baseline':p['B'],'rank1':p['B']+p['rank1_D'],
            'spectral':p['B']+p['spectral_D'], 'rank1_control':p['rank1_controlled'],
            'spectral_control':p['spectral_controlled']}
        for name,pred in preds.items():
            error=y-pred;pe=proj(error)
            for rep in range(y.shape[1]):
                rows.append(dict(configuration=label,group=label.split('_fold')[0],method=name,
                    repeat=rep+1,full_sse=float(np.sum(error[:,rep]**2)),
                    partner_sse=float(np.sum(pe[:,rep]**2)),
                    reference_full_sse=float(np.sum(e0[:,rep]**2)),
                    reference_partner_sse=float(np.sum(proj(e0)[:,rep]**2)),
                    original_full_sse=float(np.sum(old_e[:,rep]**2)),
                    original_partner_sse=float(np.sum(proj(old_e)[:,rep]**2))))
        for fam in ['rank1','spectral']:
            change=p[fam+'_controlled']-(p['B']+p[fam+'_D'])
            endpoint=change-proj(change);a=y-p['B']-p[fam+'_D'];b=y-p[fam+'_controlled']
            checks.append(dict(configuration=label,family=fam,
                endpoint_change_max=float(abs(endpoint).max()),
                sse_identity_gap=float(abs(np.sum(a*a-b*b)-np.sum(proj(a)**2-proj(b)**2)))))
    d=pd.DataFrame(rows);d.to_csv(R/'results/SATURATION_SCORES.csv',index=False)
    fields=['full_sse','partner_sse','reference_full_sse','reference_partner_sse','original_full_sse','original_partner_sse']
    for name,keys in [('SATURATION_POOLED',['group','method']),('SATURATION_BY_CONFIGURATION',['configuration','group','method']),('SATURATION_BY_REPEAT',['group','method','repeat'])]:
        a=d.groupby(keys,as_index=False)[fields].sum()
        for met in ['full','partner']:
            a[met+'_gain_pct']=100*(1-a[met+'_sse']/a['reference_'+met+'_sse'])
            a[met+'_gain_original_pct']=100*(1-a[met+'_sse']/a['original_'+met+'_sse'])
        a.to_csv(R/'results'/f'{name}.csv',index=False)
    return checks

def han_attribution(base):
    h=base/'experiments/han_native';data={};rows=[];maxgap=0.;maxtemplate=0.
    for panel in ['full_raw','full_published']:
        a=np.load(h/'inputs'/f'{panel}_inputs.npz');o=np.load(h/'inputs'/f'{panel}_outcomes.npz')
        assert np.array_equal(a['pairs'],o['pairs']);idx={g:i for i,g in enumerate(a['genes'])}
        edges=np.array([[idx[g],idx[j]] for g,j in a['pairs']])
        data[panel]=(a['singles'],edges,o['outcomes'])
    query_files=sorted((h/'fitted').rglob('QUERY.npz'));assert len(query_files)==60
    for f in query_files:
        panel,baseline,seed,fold,_=f.relative_to(h/'fitted').parts
        p=np.load(f);singles,edges,yall=data[panel];ids=p['query_ids'];ed=edges[ids]
        y=yall[ids];proj=projector(ed,len(singles));e=proj(y-p['B']);e0=np.sum(e*e,axis=0)
        mean=singles.mean(axis=1)
        templates={'repeat_specific':singles[ed[:,0]]*singles[ed[:,1]],
            'repeat_mean':np.repeat((mean[ed[:,0]]*mean[ed[:,1]])[:,None],y.shape[1],axis=1)}
        for family in ['rank1','spectral']:
            d=proj(p[family+'_D']);total=np.sum(2*e*d-d*d,axis=0)
            for typ,t0 in templates.items():
                t=proj(t0);den=np.sum(t*t,axis=0)
                coef=np.divide(np.sum(t*d,axis=0),den,out=np.zeros_like(den),where=den>1e-14)
                dt=t*coef;dp=d-dt;g1=np.sum(2*e*dt-dt*dt,axis=0);g2=np.sum(2*e*dp-dp*dp,axis=0)
                gap=abs(total-g1-g2);maxgap=max(maxgap,float(gap.max()))
                maxtemplate=max(maxtemplate,float(abs(np.sum(dt*dp,axis=0)).max()))
                for rep in range(y.shape[1]):
                    rows.append(dict(panel=panel,baseline=baseline,seed=int(seed),fold=int(fold),
                        family=family,template=typ,repeat=rep+1,n_queries=len(ids),
                        baseline_partner_sse=float(e0[rep]),total_sse_reduction=float(total[rep]),
                        template_sse_reduction=float(g1[rep]),orthogonal_sse_reduction=float(g2[rep]),
                        correction_energy=float(np.sum(d[:,rep]**2)),template_energy=float(np.sum(dt[:,rep]**2)),
                        decomposition_gap=float(gap[rep])))
    d=pd.DataFrame(rows);d.to_csv(R/'results/HAN_DIRECTION_ATTRIBUTION.csv',index=False)
    cols=['baseline_partner_sse','total_sse_reduction','template_sse_reduction','orthogonal_sse_reduction','correction_energy','template_energy']
    for name,keys in [('HAN_ATTRIBUTION_POOLED',['panel','baseline','family','template']),
                      ('HAN_ATTRIBUTION_BY_PARTITION',['panel','baseline','family','template','seed']),
                      ('HAN_ATTRIBUTION_BY_REPEAT',['panel','baseline','family','template','repeat'])]:
        a=d.groupby(keys,as_index=False)[cols].sum()
        for component in ['total','template','orthogonal']:
            a[component+'_gain_pct']=100*a[component+'_sse_reduction']/a.baseline_partner_sse
        a.to_csv(R/'results'/f'{name}.csv',index=False)
    return dict(saved_batches=60,model_outputs=120,rows=len(rows),max_sse_decomposition_gap=maxgap,max_orthogonality_gap=maxtemplate)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True);a=ap.parse_args()
    (R/'results').mkdir(exist_ok=True)
    common_reference();checks=saturation();han=han_attribution(a.base.resolve())
    assert max(c['endpoint_change_max'] for c in checks)<1e-10
    assert max(c['sse_identity_gap'] for c in checks)<1e-10
    assert han['max_sse_decomposition_gap']<1e-8
    report=dict(status='PASS',saturation_configurations=4,control_checks=checks,han_attribution=han,
        scope='Retrospective reference sensitivity and saved-prediction attribution; no new sources or independent confirmation.')
    (R/'results/ANALYSIS_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(pd.read_csv(R/'results/SATURATION_POOLED.csv').to_string(index=False))
    print(pd.read_csv(R/'results/HAN_ATTRIBUTION_POOLED.csv').query("baseline=='saturation' and family=='spectral'").to_string(index=False))
    print(json.dumps(report,indent=2))
