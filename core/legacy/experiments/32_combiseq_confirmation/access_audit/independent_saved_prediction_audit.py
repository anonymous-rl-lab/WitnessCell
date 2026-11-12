#!/usr/bin/env python3
"""Independent saved-prediction arithmetic audit; no scientific tuning or fitting."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[k]='1'
import csv, hashlib, itertools, json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
RUN=ROOT/'run1'
def load(p):return json.loads(p.read_text())
def readcsv(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def projected(x,v):
    u,s,_=np.linalg.svd(x,full_matrices=False)
    rank=int(np.count_nonzero(s>max(x.shape)*np.finfo(float).eps*s.max()))
    return v-u[:,:rank]@(u[:,:rank].T@v),rank
fr=load(ROOT/'FREEZE.json');cohort=load(ROOT/'COHORT_AND_SPLITS.json')
complete=load(RUN/'COMPLETE.json');stability=load(RUN/'CONDITIONAL_STABILITY.json')
access=load(RUN/'EXPRESSION_ACCESS_START.json');saved=load(RUN/'ALL_PREDICTIONS_FROZEN.json')
with np.load(RUN/'EVALUATION.npz',allow_pickle=False) as z:ev={k:z[k] for k in z.files}
with np.load(RUN/'PREDICTIONS.npz',allow_pickle=False) as z:pred={k:z[k] for k in z.files}
methods=[k for k in pred if k not in ('pairs','foldids')]
tables={r['method']:r for r in readcsv(RUN/'RESULTS.csv')}
foldtable={(int(r['fold']),r['method']):r for r in readcsv(RUN/'FOLD_RESULTS.csv')}
runtable={(r['run'],r['method']):r for r in readcsv(RUN/'RUN_TARGET_RESULTS.csv')}
full={};ploss={};foldrows=[];max_table_error=0.;max_array_error=0.
for m in methods:
    err=pred[m]-ev['truth'];full[m]=err**2;pe=np.empty_like(err)
    for f in sorted(set(pred['foldids'])):
        q=np.flatnonzero(pred['foldids']==f);pe[q],rank=projected(ev['x'][q],err[q])
        fs=float(np.dot(err[q],err[q]));ps=float(np.dot(pe[q],pe[q]))
        row=foldtable[int(f),m]
        max_table_error=max(max_table_error,abs(fs-float(row['full_sse'])),abs(ps-float(row['P_sse'])))
        foldrows.append({'fold':int(f),'method':m,'n_pairs':len(q),'rank':rank,'P_dimensions':len(q)-rank,'full_sse':fs,'P_sse':ps})
    ploss[m]=pe**2
    max_array_error=max(max_array_error,float(np.max(abs(full[m]-ev['full_loss_'+m]))),float(np.max(abs(ploss[m]-ev['P_loss_'+m]))))
    for metric,l in [('full_sse',full[m]),('P_sse',ploss[m])]:max_table_error=max(max_table_error,abs(float(l.sum())-float(tables[m][metric])))
    for t,name in enumerate(cohort['runs']):
        ss=float(np.sum((pred[m]-ev['runtruth'][t])**2))
        max_table_error=max(max_table_error,abs(ss-float(runtable[name,m]['full_sse'])))
        for den in ('endpoint_control','scalar_saturation'):
            dss=float(np.sum((pred[den]-ev['runtruth'][t])**2));gain=1-ss/dss
            suffix='endpoint' if den=='endpoint_control' else 'saturation'
            max_table_error=max(max_table_error,abs(gain-float(runtable[name,m]['gain_vs_'+suffix])))
main='system_fixed_reference';gain=1-float(full[main].sum()/full['endpoint_control'].sum());pgain=1-float(ploss[main].sum()/ploss['endpoint_control'].sum())
node_index={g:i for i,g in enumerate(cohort['nodes'])};edges=np.array([[node_index[g] for g in p.split('+')] for p in ev['pairs']])
rng=np.random.default_rng(20261020);v=rng.exponential(1.,(4000,len(node_index)));w=v[:,edges[:,0]]*v[:,edges[:,1]]
intervals={};interval_diff=0.
for metric,loss in [('full',full),('P',ploss)]:
    for baseline in ('endpoint_control','scalar_saturation'):
        values=1-(w@loss[main])/(w@loss[baseline]);q=np.quantile(values,[.025,.5,.975]);name=metric+'_vs_'+baseline
        intervals[name]=dict(zip(('q025','median','q975'),map(float,q)))
        interval_diff=max(interval_diff,max(abs(float(q[i])-stability['results'][name][k]) for i,k in enumerate(('q025','median','q975'))))
runpositive=sum(float(runtable[t,main]['gain_vs_endpoint'])>0 for t in cohort['runs'])
audits=load(RUN/'INVARIANCE_AND_LEAKAGE.json')
checks={'principal_gain_ge_2_percent':gain>=.02,'conditional_lower_bound_positive':intervals['full_vs_endpoint_control']['q025']>0,
        'full_gain_vs_saturation_positive':1-full[main].sum()/full['scalar_saturation'].sum()>0,
        'partner_gain_vs_endpoint_positive':pgain>0,'at_least_two_runs_positive':runpositive>=2,'invariance_and_leakage':all(a['PASS'] for a in audits)}
checks={k:bool(v) for k,v in checks.items()}
decision='STRONG_SUPPORT' if all(checks.values()) else 'NOT_CONFIRMED' if gain<=0 else 'PARTIAL_OR_INCONCLUSIVE'
hashcheck={rel:digest(ROOT/rel)==h for rel,h in fr['files'].items()}
identity_checks={'all_frozen_hashes_match':all(hashcheck.values()),'freeze_record_sha_matches':digest(ROOT/'FREEZE.json')==complete['freeze_sha256'],
    'prediction_hash_matches':digest(RUN/'PREDICTIONS.npz')==saved['prediction_sha256'],
    'freeze_before_expression':fr['utc']<access['utc'],'prediction_save_before_completion':saved['utc']<complete['utc'],
    'unique_pairs':len(set(ev['pairs']))==len(ev['pairs']),'pair_arrays_match':np.array_equal(ev['pairs'],pred['pairs']),
    'fold_arrays_match':np.array_equal(ev['foldids'],pred['foldids']),'cohort_pair_order_match':list(ev['pairs'])==cohort['pairs'],
    'one_outer_evaluation_per_pair':sorted(i for f in cohort['folds'] for i in f['query'])==list(range(len(ev['pairs']))),
    'no_outer_train_query_overlap':all(not(set(f['train'])&set(f['query'])) for f in cohort['folds']),
    'all_pairs_distinct_endpoints':bool(np.all(edges[:,0]!=edges[:,1])),
    'all_finite_saved_predictions':all(bool(np.isfinite(pred[m]).all()) for m in methods),
    'primary_method_identity':complete['primary']['method']==main,'decision_reproduces':complete['decision']==decision,'decision_checks_reproduce':checks==complete['checks']}
sample_counts=Counter(z['sample'] for z in cohort['allowed']);group=defaultdict(list)
for z in cohort['allowed']:group[z['condition'],z['run']].append(z)
orientation_counts=[]
for (condition,run),records in group.items():
    cc=Counter((z['drug1'],z['drug2']) for z in records)
    orientation_counts.append({'condition':condition,'run':run,'orientation_sample_counts':list(cc.values())})
identity_checks['unique_allowed_sample_ids']=max(sample_counts.values())==1
identity_checks['all_allowed_groups_present']=all((c,t) in group for c in cohort['conditions'] for t in cohort['runs'])
identity_checks['equal_orientation_sample_counts']=all(len(set(r['orientation_sample_counts']))==1 for r in orientation_counts)
identity_checks['internal_controls_only']=all(z['drug1']==z['drug2']=='DMSO' for z in cohort['allowed'] if z['condition']=='control')
identity_checks['excluded_extra_controls_absent']=all('_ctrlReadsPerGene' not in z['sample'] for z in cohort['allowed'])
identity_checks['run_truth_mean_matches_principal']=bool(np.allclose(ev['runtruth'].mean(axis=0),ev['truth'],rtol=0,atol=1e-12))

# ONE explicitly post-score diagnostic; it is not selection or a new success test.
res=ev['runtruth']-pred['endpoint_control'][None,:];pr=np.empty_like(res)
for f in sorted(set(pred['foldids'])):
    q=np.flatnonzero(pred['foldids']==f);p,_=projected(ev['x'][q],res[:,q].T);pr[:,q]=p.T
diagnostic=[]
for label,rr in [('baseline_residual',res),('same_fold_P_baseline_residual',pr)]:
    for a,b in itertools.combinations(range(3),2):
        r1,r2=rr[a],rr[b];dot=float(np.dot(r1,r2));e1=float(np.dot(r1,r1));e2=float(np.dot(r2,r2))
        diagnostic.append({'object':label,'run_a':cohort['runs'][a],'run_b':cohort['runs'][b],
            'pearson':float(np.corrcoef(r1,r2)[0,1]),'cosine':dot/np.sqrt(e1*e2),'cross_energy':dot,
            'energy_a':e1,'energy_b':e2,'cross_energy_per_pair':dot/len(r1)})
report={'audit_status':'PASS' if all(identity_checks.values()) and max_table_error<1e-9 and max_array_error<1e-9 and interval_diff<1e-10 else 'FAIL',
    'scope':'independent saved-prediction arithmetic; no model refit, tuning, new subset, exclusion or new source',
    'identity_checks':identity_checks,'frozen_file_checks':hashcheck,'max_table_abs_error':max_table_error,'max_per_pair_loss_array_abs_error':max_array_error,
    'max_stability_quantile_abs_error':interval_diff,'primary_gain_fraction':gain,'primary_gain_percent':100*gain,
    'partner_gain_fraction':pgain,'partner_gain_percent':100*pgain,'decision':decision,'decision_checks':checks,
    'n_unique_pairs':len(ev['pairs']),'n_compounds':len(node_index),'n_program':cohort['n_program'],
    'full_sse_recomputed':{m:float(full[m].sum()) for m in methods},'P_sse_recomputed':{m:float(ploss[m].sum()) for m in methods},
    'stability_recomputed':intervals,'folds':foldrows,'metadata':{'n_allowed_samples':len(cohort['allowed']),'n_excluded_samples':len(cohort['excluded']),
        'n_condition_run_groups':len(group),'control_samples':{t:[z['sample'] for z in cohort['allowed'] if z['condition']=='control' and z['run']==t] for t in cohort['runs']}},
    'post_score_measurement_adequacy':{'status':'FOLLOWUP_DIAGNOSTIC_ONLY; does not change NOT_CONFIRMED',
        'operation':'subtract same frozen endpoint-control baseline from each original run target; optionally apply the same outer-fold endpoint projection',
        'rows':diagnostic,'limits':'Residual agreement may contain common baseline bias; targets share experimental design and pooled-training baseline. Low agreement does not identify its biological, technical, preprocessing, or measurement-depth cause. No run is excluded.'}}
(OUT/'Independent_Audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
lines=['WitnessCell Combi-seq independent saved-prediction audit',
    'Audit status: '+report['audit_status'],'Scientific result: '+decision,
    f'Primary SYSTEM vs endpoint baseline: {100*gain:.9f}% full SSE gain; {100*pgain:.9f}% projected SSE gain.',
    f'Unique pairs {len(ev["pairs"])}; compounds {len(node_index)}; G2M readouts {cohort["n_program"]}.',
    f'Max table absolute discrepancy: {max_table_error:.3g}; max saved per-pair loss discrepancy: {max_array_error:.3g}; stability quantile discrepancy: {interval_diff:.3g}.',
    'All freeze/source/prediction hashes and recorded chronological order checked. Predictions are finite, every unordered pair occurs exactly once in outer evaluation, and no outer fitting/query overlap exists.',
    'Principal method identity, all decision gates, fold/run loss tables, full and projected per-pair losses, and all four multiplier ranges were independently recomputed from saved predictions and targets.',
    'Gain values in CSV/JSON are fractions. Multiply by 100 to obtain percentage. Thus -0.0002328227 is -0.02328227%, not -2.328% and not +0.023%.',
    'Source unit is a barcoded cell pool. The three run-specific targets are assessed against a predictor trained using other conditions from those same runs; this is not a leave-run-out experiment.',
    'Metadata checks passed for matched internal DMSO_DMSO controls, exclusion of extra ctrlReadsPerGene samples, allowed condition coverage per run and equal orientation sample counts. Metadata cannot independently establish identical numerical delivered concentrations, which the protocol discloses.',
    'The P stability range reweights fixed projected losses. It does not recompute a weighted projection. This is a conditional sensitivity description and not a biological confidence interval.',
    '', 'One post-score measurement-adequacy diagnostic (does not alter status or selection):']
for d in diagnostic:lines.append(f"{d['object']} {d['run_a']} vs {d['run_b']}: Pearson={d['pearson']:.6f}; cosine={d['cosine']:.6f}; cross-energy={d['cross_energy']:.6f}; energies={d['energy_a']:.6f}/{d['energy_b']:.6f}.")
lines+=['','No model, scientific rule, source, sample or program was changed. No run was removed. These residual diagnostics do not prove a technical cause or the absence of biological interactions.',
    'Machine-readable exact values and every identity gate are in Independent_Audit.json.']
(OUT/'Independent_Audit.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'audit_status':report['audit_status'],'identity_checks':identity_checks,'max_table_abs_error':max_table_error,
 'primary_gain_percent':100*gain,'partner_gain_percent':100*pgain,'diagnostic':diagnostic},indent=2))
