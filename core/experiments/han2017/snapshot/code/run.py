"""One scalar control rule; each fold is saved before its query outcomes are read."""
from common import *
import argparse
import pandas as pd
PANELS=['common_raw','common_published','full_raw','full_published']

def calibrate(data,train_ids,train_y,seed,fold):
    ids=np.asarray(train_ids);position={i:j for j,i in enumerate(ids)};k=data['D'].shape[2]
    cross=np.zeros(k);energy=np.zeros(k);geometry=[]
    for block in split(data['pairs'],ids,2,seed+1000+37*fold):
        local=np.array([position[i] for i in block]);p=Projection(data['pairs'][block])
        r=p.P(train_y[local,:,None]-data['M0'][block]);d=p.P(data['D'][block])
        cross+=np.sum(r*d,axis=(0,1));energy+=np.sum(d*d,axis=(0,1))
        geometry.append(dict(n=len(block),genes=len(p.genes),contrasts=p.dim))
    alpha=np.divide(cross,energy,out=np.zeros(k),where=energy>1e-14);alpha=np.clip(alpha,0,1)
    d=data['D'][ids];r=train_y[:,:,None]-data['M0'][ids];ef=np.sum(d*d,axis=(0,1));cf=np.sum(r*d,axis=(0,1))
    gamma=np.clip(np.divide(cf,ef,out=np.zeros(k),where=ef>1e-14),0,1)
    return dict(alpha=alpha.tolist(),gamma=gamma.tolist(),cross=cross.tolist(),energy=energy.tolist(),calibration_blocks=geometry)

def arms(data,ids,state,panel):
    p=Projection(data['pairs'][ids]);m=data['M0'][ids];d=data['D'][ids];b=data['B'][ids,:,None]
    pd=p.P(d);alpha=np.asarray(state['alpha']);gamma=np.asarray(state['gamma']);src=data['source_alpha']
    if panel.startswith('full'):
        out={'original':m+pd,'partner_off':m,'source_control':m+pd*src,'target_control':m+pd*alpha,'single_additive':np.broadcast_to(b,m.shape).copy()}
    else:
        off=m+d-pd
        out={'M0':m,'single_additive':np.broadcast_to(b,m.shape).copy(),
             'original':m+d,'partner_off':off,'source_control':off+pd*src,'target_control':off+pd*alpha,
             'whole_same_alpha':m+d*alpha,'whole_full_gamma':m+d*gamma}
        hb=p.H(b);pm=p.P(m)
        for name,value in [('original',1.),('partner_off',0.),('source_control',src),('target_control',alpha)]:
            out['anchored_'+name]=hb+pm+pd*value
    return p,out

def predict_fold(panel,seed,fold,data):
    dest=ROOT/'results'/panel/str(seed)/str(fold);dest.mkdir(parents=True,exist_ok=True)
    if (dest/'PREDICTIONS_LOCKED.json').exists():return
    splitdir=ROOT/'data/splits'/panel/str(seed)/str(fold)
    ident=json.loads((splitdir/'identities.json').read_text());q=np.asarray(ident['query']);tr=load(splitdir/'calibration.npz')
    assert set(tr['ids'])==set(ident['train']) and not set(tr['ids'])&set(q)
    state=calibrate(data,tr['ids'],tr['y'],seed,fold);p,outputs=arms(data,q,state,panel)
    dump(dest/'PARAMETERS.json',state)
    # Save controlled query values; other locked arms exactly reconstruct from frozen arrays + coefficients.
    save={'ids':q,'target_control':outputs['target_control']}
    if 'anchored_target_control' in outputs:save['anchored_target_control']=outputs['anchored_target_control']
    np.savez_compressed(dest/'QUERY_PREDICTIONS.npz',**save)
    hashes={name:hashlib.sha256(np.ascontiguousarray(pred).tobytes()).hexdigest() for name,pred in outputs.items()}
    dump(dest/'PREDICTIONS_LOCKED.json',dict(utc=now(),panel=panel,seed=seed,fold=fold,n_query=len(q),n_calibration=len(tr['ids']),query_contrasts=p.dim,query_labels_read=False,base_model_fit_calls=0,control_coefficients_target_fitted=True,all_arm_array_sha256=hashes,files={name:sha(dest/name) for name in ['PARAMETERS.json','QUERY_PREDICTIONS.npz']},input_prediction_sha256=sha(ROOT/f'data/{panel}_frozen_predictions.npz')))

def score_fold(panel,seed,fold,data):
    dest=ROOT/'results'/panel/str(seed)/str(fold)
    if (dest/'COMPLETE.json').exists():return
    lock=json.loads((dest/'PREDICTIONS_LOCKED.json').read_text())
    for name,h in lock['files'].items():assert sha(dest/name)==h
    state=json.loads((dest/'PARAMETERS.json').read_text());saved=load(dest/'QUERY_PREDICTIONS.npz');q=saved['ids'];p,outputs=arms(data,q,state,panel)
    for name,pred in outputs.items():assert hashlib.sha256(np.ascontiguousarray(pred).tobytes()).hexdigest()==lock['all_arm_array_sha256'][name]
    assert np.array_equal(saved['target_control'],outputs['target_control'])
    # Query outcomes are accessed only after the above durable prediction lock.
    truth=load(ROOT/'data/splits'/panel/str(seed)/str(fold)/'query_outcomes.npz');assert np.array_equal(q,truth['ids']);y=truth['y']
    rows=[];checks=[];reps=[('pooled',slice(None))]
    if y.shape[1]==2:reps=[('repeat1',slice(0,1)),('repeat2',slice(1,2))]+reps
    for method,pred in outputs.items():
        e=y[:,:,None]-pred;pe=p.P(e);he=e-pe
        for label,col in reps:
            fs=np.sum(e[:,col,:]**2,axis=(0,1));ps=np.sum(pe[:,col,:]**2,axis=(0,1));es=np.sum(he[:,col,:]**2,axis=(0,1))
            for j,model in enumerate(data['models']):rows.append(dict(panel=panel,seed=seed,fold=fold,model=model,method=method,repeat=label,n=len(q)*e[:,col,:].shape[1],full_sse=fs[j],partner_sse=ps[j],endpoint_sse=es[j]))
    delta=outputs['target_control']-outputs['original'];endpoint=float(np.max(abs(p.H(delta))))
    e0=y[:,:,None]-outputs['original'];e1=y[:,:,None]-outputs['target_control']
    equality=float(np.max(abs(np.sum(e0*e0-e1*e1,axis=(0,1))-np.sum(p.P(e0)**2-p.P(e1)**2,axis=(0,1)))))
    anchored=0.
    if 'anchored_target_control' in outputs:anchored=float(np.max(abs(p.P(outputs['anchored_target_control'])-p.P(outputs['target_control']))))
    assert endpoint<1e-7 and equality<1e-7*max(1.,float(np.sum(e0*e0))) and anchored<1e-7
    pd.DataFrame(rows).to_csv(dest/'SCORES.csv',index=False)
    dump(dest/'COMPLETE.json',dict(utc=now(),endpoint_max_abs=endpoint,sse_identity_error=equality,anchored_partner_change=anchored,n_rows=len(rows),score_sha256=sha(dest/'SCORES.csv'),prediction_lock_sha256=sha(dest/'PREDICTIONS_LOCKED.json')))
    print(json.dumps(dict(stage='fold_complete',panel=panel,seed=seed,fold=fold,n_query=len(q))),flush=True)

def summarize():
    paths=sorted((ROOT/'results').glob('*/*/*/SCORES.csv'));df=pd.concat([pd.read_csv(p) for p in paths],ignore_index=True)
    g=df.groupby(['panel','seed','model','method','repeat'],as_index=False)[['n','full_sse','partner_sse','endpoint_sse']].sum();g.to_csv(ROOT/'results/SEED_SCORES.csv',index=False)
    total=g.groupby(['panel','model','method','repeat'],as_index=False)[['n','full_sse','partner_sse','endpoint_sse']].sum()
    for key in ['full','partner','endpoint']:total[key+'_mse']=total[key+'_sse']/total['n']
    total.to_csv(ROOT/'results/POOLED_SCORES.csv',index=False)
    comparisons=[]
    for (panel,model,rep),d in total.groupby(['panel','model','repeat']):
        d=d.set_index('method')
        for controlled in ['target_control','anchored_target_control']:
            if controlled not in d.index:continue
            for ref in d.index:
                if ref==controlled:continue
                row=dict(panel=panel,model=model,repeat=rep,controlled=controlled,reference=ref)
                for key in ['full','partner']:
                    a=d.loc[ref,key+'_sse'];b=d.loc[controlled,key+'_sse'];row[key+'_sse_reduction']=a-b;row[key+'_reduction_pct']=100*(1-b/a) if a>0 else np.nan
                comparisons.append(row)
    pd.DataFrame(comparisons).to_csv(ROOT/'results/COMPARISONS.csv',index=False)
    print('SCORING_SUMMARY_COMPLETE',len(paths),'folds',len(df),'score rows',flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panels',nargs='+',default=PANELS);ap.add_argument('--phase',choices=['all','predict','score'],default='all');args=ap.parse_args()
    for panel in args.panels:
        data=load(ROOT/f'data/{panel}_frozen_predictions.npz')
        for seed in range(20261108,20261113):
            for fold in range(3):
                if args.phase in ['all','predict']:predict_fold(panel,seed,fold,data)
                if args.phase in ['all','score']:score_fold(panel,seed,fold,data)
    if args.phase in ['all','score']:summarize()
if __name__=='__main__':main()
