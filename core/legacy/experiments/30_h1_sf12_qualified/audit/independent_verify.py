"""Independent direct-vector and sufficient-statistic audit; never fits a new model.

Requires the completed freeze and scores, plus original source measurement asset.
Does not import the experiment's fitting or scoring functions.
"""
from pathlib import Path
import json, hashlib, datetime
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'run'

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    spec = json.loads((OUT/'PROTOCOL.json').read_text())['spec']
    freeze = json.loads((OUT/'ALL_PREDICTIONS_FROZEN.json').read_text())
    reported = json.loads((OUT/'SCORES.json').read_text())
    assert digest(OUT/'PROTOCOL.json') == freeze['protocol_sha256']
    for name, value in freeze['files'].items():
        assert digest(OUT/name) == value, name
    for name, value in spec['code_hashes'].items():
        assert digest(ROOT/name) == value, name
    source = next(Path(n) for n in spec['input_hashes'] if n.endswith('TRAINING_BATCH_ONLY.npz'))
    assert digest(source) == spec['input_hashes'][str(source)]
    with np.load(source, allow_pickle=False) as z:
        means, obs, batches, pi, nodes = [z[n].copy() for n in ['batch_values','obs_condition','obs_batch','pi','nodes']]
    names = list(spec['methods'])
    c = np.array(spec['contrast'])
    P = np.outer(c,c)
    assert np.isclose(c@c,1)
    independent = {}
    case_geometry = {}
    edge_losses = {}
    max_coefficient = max_stat_sum = max_prediction = max_score = max_decomposition = max_quadratic = 0.
    with np.load(OUT/'EVALUATION_STATS.npz', allow_pickle=False) as z:
        evalstats = {n:z[n].copy() for n in z.files}
    for task in spec['tasks']:
        name = task['name']; folder = OUT/name
        qids = task['query_ids']
        y = []
        for query in qids:
            rows = np.flatnonzero(obs == len(nodes)+query)
            weights = pi[batches[rows]]
            assert abs(weights.sum()-1) < 1e-12
            y.append(sum((float(w)*means[row] for row,w in zip(rows,weights) if w>0), np.zeros(means.shape[1])))
        y = np.stack(y)
        with np.load(folder/'PREDICTIONS.npz', allow_pickle=False) as z:
            pred, h = z['predictions'].copy(), z['h_case_R'].copy()
            assert list(z['names']) == names and list(z['query_ids']) == qids
            assert np.array_equal(z['contrast'],c)
        with np.load(folder/'TRAIN_STATS.npz', allow_pickle=False) as z:
            stats={n:z[n].copy() for n in z.files}
        M = sum((stats[f'e{i}_M'] for i in range(len(task['episodes']))), np.zeros((1,1)))
        b = sum((stats[f'e{i}_b'] for i in range(len(task['episodes']))), np.zeros(1))
        max_stat_sum=max(max_stat_sum,float(np.max(abs(M-stats['M']))),float(np.max(abs(b-stats['b']))))
        m=float(M[0,0]); rhs=float(b[0]); scale=max(np.sqrt(max(m,0)),1e-15)
        theta=rhs/(m+spec['case_R']['rho']*scale*scale)
        max_coefficient=max(max_coefficient,abs(theta-float(stats['theta'][0])))
        rebuild=pred[1]+h[:,:,0]*theta
        max_prediction=max(max_prediction,float(np.max(abs(rebuild-pred[2]))))
        assert set(task['train_ids']) == set(range(21))-set(qids)
        for ep in task['episodes']:
            assert not set(qids)&set(ep['train_ids']+ep['val_ids'])
            assert set(ep['train_ids']) == set(task['train_ids'])-set(ep['val_ids'])
        errors=y[None]-pred
        full=np.array([float(np.linalg.norm(e)**2) for e in errors])
        partner=np.array([float(np.linalg.norm(c@e)**2) for e in errors])
        endpoint=np.array([float(np.linalg.norm(e-P@e)**2) for e in errors])
        max_decomposition=max(max_decomposition,float(np.max(abs(full-partner-endpoint))))
        edge_losses[name]=np.array([[float(row@row) for row in e] for e in errors])
        rows={}
        for j,method in enumerate(names):
            r=c@(y-pred[0]); delta=c@(pred[j]-pred[0]); denom=np.linalg.norm(r)*np.linalg.norm(delta)
            row=dict(full_SSE=float(full[j]),partner_SSE=float(partner[j]),
                full_gain_vs_M0_pct=float(100*(1-full[j]/full[0])),partner_gain_vs_M0_pct=float(100*(1-partner[j]/partner[0])),
                full_gain_vs_old12_pct=float(100*(1-full[j]/full[1])),partner_gain_vs_old12_pct=float(100*(1-partner[j]/partner[1])),
                partner_direction_cosine_vs_M0=float(r@delta/denom) if denom>1e-30 else None,
                per_query_full_SSE=edge_losses[name][j].tolist())
            expected=reported['tasks'][name]['methods'][method]
            for key,value in row.items():
                if value is None: assert expected[key] is None
                else:max_score=max(max_score,float(np.max(abs(np.asarray(value)-np.asarray(expected[key])))))
            rows[method]=row
        for kind,ref in [('F',full[2]),('C',partner[2])]:
            baseline=full[1] if kind=='F' else partner[1]
            qm=float(evalstats[name+'_M'+kind][0,0]); qb=float(evalstats[name+'_b'+kind][0])
            quadratic=baseline-2*theta*qb+theta*theta*qm
            max_quadratic=max(max_quadratic,abs(quadratic-ref))
        independent[name]=rows
        delta=pred[2]-pred[1]
        energy=float(np.sum(delta*delta)); partner_energy=float(np.sum((c@delta)**2))
        case_geometry[name]=dict(training_M=m,training_b=rhs,theta=theta,
            unregularized_scalar=rhs/m if m>0 else None,
            full_delta_energy=energy,partner_delta_energy=partner_energy,
            endpoint_energy_fraction=1-partner_energy/energy if energy>0 else None)
    aggregate={}
    for group,selected in [('all5',[t['name'] for t in spec['tasks']]),('exclude_Q4',[t['name'] for t in spec['tasks'] if t['name']!='Q4'])]:
        fs=np.array([sum(independent[n][method]['full_SSE'] for n in selected) for method in names])
        ps=np.array([sum(independent[n][method]['partner_SSE'] for n in selected) for method in names])
        unique=sorted({i for t in spec['tasks'] if t['name'] in selected for i in t['query_ids']})
        eb=[]
        for edge in unique:
            observations=[]
            for t in spec['tasks']:
                if t['name'] in selected and edge in t['query_ids']:
                    observations.append(edge_losses[t['name']][:,t['query_ids'].index(edge)])
            eb.append(np.mean(observations,axis=0))
        eb=np.mean(eb,axis=0)
        aggregate[group]={}
        for j,method in enumerate(names):
            row=dict(pooled_partner_SSE=float(ps[j]),pooled_full_SSE=float(fs[j]),
                partner_gain_vs_M0_pct=float(100*(1-ps[j]/ps[0])),partner_gain_vs_old12_pct=float(100*(1-ps[j]/ps[1])),
                pooled_full_gain_vs_M0_pct=float(100*(1-fs[j]/fs[0])),edge_balanced_full_SSE=float(eb[j]),
                edge_balanced_full_gain_vs_M0_pct=float(100*(1-eb[j]/eb[0])),edge_balanced_full_gain_vs_old12_pct=float(100*(1-eb[j]/eb[1])))
            for key,value in row.items():max_score=max(max_score,abs(value-reported['aggregate'][group][method][key]))
            aggregate[group][method]=row
    for left in spec['tasks']:
        remain=[t['name'] for t in spec['tasks'] if t['name']!=left['name']]
        den=sum(independent[n]['M0']['partner_SSE'] for n in remain)
        for method in names:
            value=100*(1-sum(independent[n][method]['partner_SSE'] for n in remain)/den)
            max_score=max(max_score,abs(value-reported['leave_one_task_out_partner_gain_vs_M0_pct'][left['name']][method]))
    for method in names[1:]:
        gains=[independent[t['name']][method]['partner_gain_vs_M0_pct'] for t in spec['tasks']]
        assert sum(g>0 for g in gains)==reported['direction'][method]['positive_blocks']
        max_score=max(max_score,abs(float(np.median(gains))-reported['direction'][method]['median_partner_gain_pct']))
    checks=dict(coefficient_max_abs=max_coefficient,training_stat_sum_max_abs=max_stat_sum,
        prediction_rebuild_max_abs=max_prediction,score_and_aggregate_max_abs=max_score,
        orthogonal_decomposition_max_abs=max_decomposition,quadratic_SSE_max_abs=max_quadratic)
    assert max(checks.values())<1e-9,checks
    report=dict(status='PASS',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,
        scope='Independent direct fixed-seven-library means, exact one-coefficient ridge, vector scoring, edge-balanced aggregation, quadratic SSE and freeze hashes; no new candidate fit or selection.',
        tasks=independent,aggregate=aggregate,
        posthoc_case_R_geometry=case_geometry,
        interpretation='For one active normalized feature, rho=.1 gives theta=b/(1.1*M), a fixed9.09percent shrinkage of OLS; it does not bound complete-output correction. The geometry diagnostic is post hoc and not used for selection.',
        code_sha256=digest(__file__))
    (ROOT/'audit/INDEPENDENT_VERIFY.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    print(json.dumps(dict(status='PASS',checks=checks,aggregate=aggregate),ensure_ascii=False))

if __name__=='__main__':main()
