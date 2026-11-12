#!/usr/bin/env python3
"""Retrospective analyses of frozen v24 predictions. No model training."""
from pathlib import Path
import csv,gzip,json,hashlib
import numpy as np
R=Path(__file__).resolve().parents[1]; D=R/'revision'; P=R/'evidence/source_data/publication'; G=R/'evidence/source_data/gse220974'
def read(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))
def write(name,rows):
    with (D/name).open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
src=read(D/'source/benchmark_per_operation.csv.gz'); src=[r for r in src if r['DEG']=='100']
meta=read(D/'source/benchmark_operation_identity.csv.gz');mid={(r['dataset'],r['condition']+'_'+r['seed']):r['subgroup'] for r in meta}
methods=sorted({r['method'] for r in src});wi=methods.index('WitnessCell');assert len(methods)==16
keys=sorted({(r['DataSet'],r['op']) for r in src});assert len(keys)==654
table={(r['DataSet'],r['op'],r['method']):r for r in src}
rankcols=['Rank_pcc','Rank_mse','Rank_edistance','Rank_sym','Rank_was','Rank_deg']
for key in keys:
    sub=[table[key+(m,)] for m in methods]
    av=np.array([np.mean([float(r[c]) for c in rankcols]) for r in sub])
    actual=1+(av[:,None]>av[None,:]).sum(1)
    assert np.array_equal(actual,[float(r['Rank']) for r in sub]),key
    assert key in mid,key
def bootstrap(kk,values,seed):
    # Sum/count aggregation preserves the original operation weighting.
    cl=sorted({(d,o.rsplit('_',1)[0]) for d,o in kk});ci={x:i for i,x in enumerate(cl)}
    sums=np.zeros((len(cl),values.shape[1]));counts=np.zeros(len(cl))
    for k,v in zip(kk,values):i=ci[(k[0],k[1].rsplit('_',1)[0])];sums[i]+=v;counts[i]+=1
    rng=np.random.default_rng(seed);out=[]
    for _ in range(20):
        ix=rng.integers(0,len(cl),(1000,len(cl)));out.append(sums[ix].sum(1)/counts[ix].sum(1)[:,None])
    return np.concatenate(out),len(cl)
contexts={'all':keys,'double':[k for k in keys if '+' in k[1]],'single':[k for k in keys if '+' not in k[1]]}
for state in sorted({mid[k] for k in keys}):contexts[state]=[k for k in keys if mid[k]==state]
for ds in sorted({k[0] for k in keys}):contexts[ds+'_double']=[k for k in contexts['double'] if k[0]==ds]
leaders=[];intervals=[];absolute=[]
for n,(label,kk) in enumerate(contexts.items()):
    v=np.array([[float(table[k+(m,)]['Rank']) for m in methods] for k in kk]);means=v.mean(0);pos=1+(means[:,None]>means[None,:]).sum(1)
    boot,nc=bootstrap(kk,v,20261107+n)
    for j,m in enumerate(methods):
        leaders.append({'population':label,'method':m,'operation_seed_units':len(kk),'clusters':nc,'mean_operation_rank':means[j],'position':int(pos[j])})
        if j!=wi:
            diff=boot[:,wi]-boot[:,j];lo,hi=np.quantile(diff,[.025,.975])
            intervals.append({'population':label,'reference':m,'difference_WitnessCell_minus_reference':means[wi]-means[j],'low':lo,'high':hi,'clusters':nc,'draws':20000})
    if label.endswith('_double'):
        # Keep dataset scales separate. No cross-dataset raw-MSE pooling.
        for metric in ['mse','cor']:
            mv=np.array([[float(table[k+(m,)][metric]) for m in methods] for k in kk]);mb,_=bootstrap(kk,mv,20261207+n)
            for j,m in enumerate(methods):
                if j==wi:continue
                dif=mb[:,wi]-mb[:,j];lo,hi=np.quantile(dif,[.025,.975]);ref=mv[:,j].mean()
                absolute.append({'population':label,'metric':metric,'reference':m,'WitnessCell_mean':mv[:,wi].mean(),'reference_mean':ref,'paired_difference':mv[:,wi].mean()-ref,'low':lo,'high':hi,'relative_MSE_gain':1-mv[:,wi].mean()/ref if metric=='mse' and ref>0 else ''})
write('BENCHMARK_DECOMPOSITION.csv',leaders);write('BENCHMARK_PAIRED_INTERVALS.csv',intervals);write('BENCHMARK_DATASET_ABSOLUTE.csv',absolute)
write('BENCHMARK_OPERATION_MAP.csv',[{'dataset':k[0],'operation':k[1],'population':'double' if '+' in k[1] else 'single','exposure_state':mid[k]} for k in keys])

a=np.load(P/'RNA_signed_rectangle_predictions.npz');base=(a['M0']-a['truth_effect'])**2;base_sse=base.sum();names=a['genes'].astype(str)
fixed={'S100A9','S100A8','S100A12','LYZ','VCAN','VIM','MS4A6A','MYL9','LGALS1','PPP1R27'};mask=np.array([g in fixed for g in names]);assert mask.sum()==10
budgets=[];genes=[];rect=[]
for model in ['M1','scalar','R']:
    err=(a[model]-a['truth_effect'])**2;delta=(base-err).sum(0)
    for label,sel in [('all',np.ones(len(names),bool)),('fixed_ten',mask),('remaining',~mask)]:
        d=delta[sel];den=base[:,sel].sum();positive=np.maximum(d,0).sum();negative=np.maximum(-d,0).sum()
        budgets.append({'model':model,'readout_set':label,'readouts':int(sel.sum()),'positive_SSE_reduction':positive,'adverse_SSE_increment':negative,'net_SSE_reduction':d.sum(),'M0_SSE_same_set':den,'positive_percent_of_all_M0':100*positive/base_sse,'adverse_percent_of_all_M0':100*negative/base_sse,'net_gain_percent_same_set':100*d.sum()/den})
    for g,d,b,e in zip(names,delta,base.sum(0),err.sum(0)):genes.append({'model':model,'gene':g,'fixed_ten':g in fixed,'M0_SSE':b,'model_SSE':e,'SSE_reduction':d})
    for i in range(11):rect.append({'model':model,'rectangle':i+1,'M0_MSE':base[i].mean(),'model_MSE':err[i].mean(),'gain_percent':100*(1-err[i].sum()/base[i].sum())})
write('RNA_ERROR_BUDGET.csv',budgets);write('RNA_ALL_READOUT_CONTRIBUTIONS.csv',genes);write('RNA_ALL_RECTANGLE_GAINS.csv',rect)
b={r['model']:r for r in budgets if r['readout_set']=='all'}
mechanism={'adverse_increment_reduction_M1_to_R':1-b['R']['adverse_SSE_increment']/b['M1']['adverse_SSE_increment'],'positive_reduction_retained_M1_to_R':b['R']['positive_SSE_reduction']/b['M1']['positive_SSE_reduction'],'interpretation':'Readout-aggregated descriptive decomposition on the same developmental panel; not independent confirmation.'}
(D/'RNA_ERROR_BUDGET_SUMMARY.json').write_text(json.dumps(mechanism,indent=2)+'\n')
choices=read(G/'CHOICES.csv');print('Choice schema:',list(choices[0]),flush=True)
summary={'status':'PASS','source_operation_rank_identity':True,'methods':len(methods),'operation_units':len(keys),'population_counts':{k:len(v) for k,v in contexts.items()},'RNA':mechanism,'no_training':True}
(D/'ANALYSIS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary),flush=True)
