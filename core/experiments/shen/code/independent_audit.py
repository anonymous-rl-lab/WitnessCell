#!/usr/bin/env python3
"""Independent measurement and loss reconstruction, without package scoring code.

This is post-freeze validation, not another model search or independent experiment.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,csv,itertools,datetime,argparse
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,ok,detail=None):
    checks.append({'check':name,'pass':bool(ok),'detail':detail})
    if not ok:raise AssertionError(name+': '+str(detail))
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'audit_recomputed')
    args=parser.parse_args()
    if args.out.exists():parser.error('Use a fresh audit output directory')
    fits=ROOT/'runs/formal';prep=ROOT/'prepared'
    freeze=json.loads((fits/'PREDICTIONS_FROZEN.json').read_text())
    official=json.loads((ROOT/'results/formal/SCORES.json').read_text())
    lock=json.loads((ROOT/'audit/SOURCE_AND_FIT_LOCK.json').read_text())
    check('source lock precedes prediction freeze precedes scoring',lock['utc']<freeze['utc']<official['utc'])
    check('all three predictions frozen',freeze['all_predictions_complete'] and set(freeze['models'])=={'HeLa','A549','293T'})
    check('protocol unchanged',sha(ROOT/'protocol/PROTOCOL.json')==freeze['protocol_sha256']==lock['protocol_sha256'])
    check('prepared input receipt unchanged',sha(prep/'PREPARATION_RECEIPT.json')==freeze['prepared_receipt_sha256'])
    source=json.loads((ROOT/'measurement/SOURCE_PROVENANCE.json').read_text())
    for item in source['source_files']:check('source SHA '+item['file'],sha(ROOT/item['file'])==item['sha256'])
    check('canonical count SHA',sha(ROOT/'measurement/canonical_counts.csv')==source['counts_sha256'])
    # Different implementation: sufficient-statistic OLS and a wide repeat table.
    raw=pd.read_csv(ROOT/'measurement/canonical_counts.csv')
    check('raw sample identity unique',not raw.duplicated(['cell_line','replicate','day','construct_id']).any())
    keys=['cell_line','replicate','day']
    controls=raw[(raw.gene_a=='NTC')&(raw.gene_b=='NTC')&(raw['count']>=20)].copy()
    controls['logc']=np.log2(controls['count'])
    refs=controls.groupby(keys).agg(ref=('logc','mean'),n_controls=('logc','size')).reset_index()
    check('all 24 samples have sufficient dual-NTC reference',len(refs)==24 and (refs.n_controls>=3).all(),int(refs.n_controls.min()))
    kept=raw[raw['count']>=20].merge(refs,on=keys,validate='many_to_one')
    kept['v']=np.log2(kept['count'])-kept.ref
    kept['t2']=kept.day**2;kept['tv']=kept.day*kept.v
    ck=['cell_line','replicate','construct_id','gene_a','gene_b']
    ag=kept.groupby(ck).agg(n=('v','size'),t=('day','sum'),v=('v','sum'),t2=('t2','sum'),tv=('tv','sum')).reset_index()
    ag=ag[ag.n>=3].copy();ag['slope']=(ag.tv-ag.t*ag.v/ag.n)/(ag.t2-ag.t**2/ag.n)
    common=ag.pivot(index=['cell_line','construct_id','gene_a','gene_b'],columns='replicate',values='slope').dropna().reset_index()
    common=common[common.gene_a!=common.gene_b].copy()
    common['a']=np.where(common.gene_a<common.gene_b,common.gene_a,common.gene_b)
    common['b']=np.where(common.gene_a<common.gene_b,common.gene_b,common.gene_a)
    means=common.groupby(['cell_line','a','b']).agg(n=(1,'size'),r1=(1,'mean'),r2=(2,'mean'))
    means=means[means.n>=3]
    raw_qc=json.loads((prep/'PREPROCESSING_QC.json').read_text())
    check('independent accepted construct-repeat count',len(ag)==raw_qc['accepted_construct_repeats'],len(ag))
    rows=[];diags=[];max_gain_diff=0.;max_measurement_diff=0.;weights_input={}
    split=json.loads((prep/'IDENTITY_SPLIT.json').read_text())
    salt='WitnessCell-Shen2017-confirmatory-pairs-v1-20261107'
    for cell in ['HeLa','A549','293T']:
        inp_path=prep/'training_inputs'/f'{cell}.json';inp=json.loads(inp_path.read_text());rec=freeze['models'][cell]
        check(cell+' input SHA',sha(inp_path)==rec['input_sha256'])
        check(cell+' prediction SHA',sha(fits/rec['prediction_file'])==rec['predictions_sha256'])
        check(cell+' metadata SHA',sha(fits/cell/'fit_metadata.json')==rec['metadata_sha256'])
        check(cell+' query interface has identities only',all(isinstance(q,list) and len(q)==2 and all(isinstance(x,str) for x in q) for q in inp['queries']))
        tr={tuple(x['pair']) for x in inp['training']};q={tuple(x) for x in inp['queries']}
        check(cell+' no physical-pair leakage',not tr&q)
        check(cell+' globally fixed split preserved',tr<=set(map(tuple,split['train'])) and q<=set(map(tuple,split['query'])))
        pairs=list(itertools.combinations(sorted(inp['genes']),2));ordered=sorted(pairs,key=lambda x:hashlib.sha256((salt+'|'+'|'.join(x)).encode()).hexdigest())
        check(cell+' nominal split independently reproduced',set(ordered[:526])==set(map(tuple,split['query'])))
        sub=means.loc[cell]
        def values(pair):return sub.loc[tuple(sorted(pair)),['r1','r2']].to_numpy(float)
        ds=np.max(np.abs(np.array([values([g,'NTC']) for g in inp['genes']])-np.array(inp['singles'])))
        dt=np.max(np.abs(np.array([values(x['pair']) for x in inp['training']])-np.array([x['response'] for x in inp['training']])))
        with np.load(prep/'sealed_evaluation'/f'{cell}.npz') as truth:
            y=truth['y'];query_pairs=truth['pairs'].tolist();check(cell+' truth identities',query_pairs==inp['queries'])
        dq=np.max(np.abs(np.array([values(x) for x in query_pairs])-y))
        err=float(max(ds,dt,dq));max_measurement_diff=max(max_measurement_diff,err)
        check(cell+' independent raw-count to outcome reconstruction',err<1e-12,err)
        with np.load(fits/rec['prediction_file']) as z:pred={k:z[k].copy() for k in z.files}
        genes=inp['genes'];lookup={g:i for i,g in enumerate(genes)};edges=np.array([[lookup[a],lookup[b]] for a,b in query_pairs]);g,h=edges.T
        X=np.zeros((len(edges),len(genes)));X[np.arange(len(edges)),g]=1;X[np.arange(len(edges)),h]=1
        def residual(v):return v-X@np.linalg.lstsq(X,v,rcond=None)[0]
        base=pred['M0'];models=json.loads((fits/cell/'results/high/M0_models.json').read_text());scales=np.std(inp['singles'],axis=0);scales[scales==0]=1
        # Direct formula from saved coefficients, no model prediction method.
        rebuilt=[]
        for k,m in enumerate(models):
            a=np.array(m['a']);e=np.array(m['e']);db=np.array(m['d_beta']);dc=np.array(m['d_c'])
            b=.5*(a[g]+a[h]+(m['beta']+db[g])*a[h]+(m['beta']+db[h])*a[g]+(m['c']+dc[g])*a[h]**2+(m['c']+dc[h])*a[g]**2)+m['e0']+e[g]+e[h]
            rebuilt.append(b*scales[k])
        check(cell+' saved coefficient baseline replay',np.max(np.abs(np.array(rebuilt).T-base))<1e-12)
        selections=json.loads((fits/cell/'results/high/selection.json').read_text())
        inner=json.loads((fits/cell/'results/high/inner_summary.json').read_text())
        check(cell+' all 24 objective-order checks',len(inner['objective_probes'])==24 and all(p['pass'] for p in inner['objective_probes']))
        for family in ['rank1','spectral']:
            check(cell+' '+family+' convergence',selections[family]['final_converged'])
            cand=[a for a in inner['grid'] if a['family']==family and a['eligible']]
            best=min(cand,key=lambda a:(a['total'],-a['fraction'],-a['nu']))
            total=sum(inner['M0_sse']);use=best['total']<total-1e-12*max(1,total)
            check(cell+' '+family+' training-only selection',selections[family]['selected']==(best['key'] if use else 'G0'))
            G=pred[family+'_G'];replay=base+G[g,h,None] if use else base
            check(cell+' '+family+' saved matrix prediction replay',np.max(np.abs(replay-pred[family+'_selected']))<1e-12)
        singles=np.array(inp['singles']);add=singles[g]+singles[h]
        for method,p in [('M0',base),('rank1_primary',pred['rank1_selected']),('spectral_secondary',pred['spectral_selected']),('single_additive',add)]:
            e=y-base;f=y-p;pe=residual(e);pf=residual(f)
            fullb=float((e*e).sum());fullm=float((f*f).sum());partb=float((pe*pe).sum());partm=float((pf*pf).sum())
            fullgain=1-fullm/fullb;partgain=1-partm/partb;reference=official['models'][cell][method]
            diff=max(abs(fullgain-reference['full_gain']),abs(partgain-reference['partner_gain']));max_gain_diff=max(max_gain_diff,diff)
            check(cell+' '+method+' least-squares score reconstruction',diff<1e-12,diff)
            rows.append(dict(cell_line=cell,method=method,train_pairs=len(tr),query_pairs=len(q),endpoint_rank=int(np.linalg.matrix_rank(X)),full_sse_M0=fullb,full_sse=fullm,partner_sse_M0=partb,partner_sse=partm,full_gain=fullgain,partner_gain=partgain))
        p=pred['rank1_selected'];e=y-base;f=y-p;d=p-base;pe=residual(e);pf=residual(f);partner_d=residual(d)
        fullimprovement=float((e*e).sum()-(f*f).sum());partnerimprovement=float((pe*pe).sum()-(pf*pf).sum());endpointimprovement=fullimprovement-partnerimprovement
        check(cell+' orthogonal error decomposition',abs(float((e*e).sum())-float((pe*pe).sum()+((e-pe)**2).sum()))<1e-12)
        diags.append(dict(cell_line=cell,endpoint_error_reduction=endpointimprovement,partner_error_reduction=partnerimprovement,total_error_reduction=fullimprovement,partner_fraction_of_correction_energy=float((partner_d*partner_d).sum()/(d*d).sum()),observed_partner_repeat_correlation=float(np.corrcoef(residual(y).T)[0,1]),partner_residual_repeat_correlation=float(np.corrcoef(pe.T)[0,1]),cosine=float((pe*partner_d).sum()/np.sqrt((pe*pe).sum()*(partner_d*partner_d).sum())),oracle_beta=float((pe*partner_d).sum()/(partner_d*partner_d).sum()),posthoc_only=True))
        diags[-1]['observed_full_repeat_correlation']=float(np.corrcoef(y.T)[0,1])
        weights_input[cell]=(edges,(pe*pe).sum(1),(pf*pf).sum(1))
    # Independent draw-by-draw implementation of fixed-projector endpoint sensitivity.
    rng=np.random.default_rng(20261107);draws=[];per={c:[] for c in weights_input}
    for _ in range(10000):
        w=rng.exponential(size=73);gains=[]
        for c,(edges,a,b) in weights_input.items():
            wp=w[edges[:,0]]*w[edges[:,1]];v=1-float(wp@b)/float(wp@a);per[c].append(v);gains.append(v)
        draws.append(np.mean(gains))
    intervals={c:np.quantile(v,[.025,.975]).tolist() for c,v in per.items()};intervals['equal_context_mean']=np.quantile(draws,[.025,.975]).tolist()
    for c in per:
        check(c+' independently reconstructed sensitivity interval',np.max(np.abs(np.array(intervals[c])-official['conditional_sensitivity']['per_context'][c]['central95']))<1e-12)
    check('equal-context sensitivity interval',np.max(np.abs(np.array(intervals['equal_context_mean'])-official['conditional_sensitivity']['equal_context_mean']['central95']))<1e-12)
    output=args.out;output.mkdir(parents=True)
    for name,data in [('INDEPENDENT_SCORES.csv',rows),('POSTFREEZE_COMPONENT_DIAGNOSTICS.csv',diags)]:
        with (output/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'passed':sum(x['pass'] for x in checks),'failed':sum(not x['pass'] for x in checks),'maximum_measurement_difference':max_measurement_diff,'maximum_gain_difference':max_gain_diff,'sensitivity':intervals,'scope':'Reconstruction from raw counts and frozen predictions; no model retuning; not a new independent biological validation.'}
    (output/'INDEPENDENT_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['passed','failed','maximum_measurement_difference','maximum_gain_difference']},indent=2))
if __name__=='__main__':main()
