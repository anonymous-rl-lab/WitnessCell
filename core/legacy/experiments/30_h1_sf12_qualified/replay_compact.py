"""Self-contained coefficient, partner-vector and aggregate-score replay.
Full RNA loss is checked from archived per-edge losses and quadratic statistics.
Recreating full RNA prediction vectors requires the original source inputs.
"""
import json
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent
S=json.loads((R/'run/SCORES.json').read_text())
P=json.loads((R/'run/PROTOCOL.json').read_text())['spec']
ev=np.load(R/'run/EVALUATION_STATS.npz',allow_pickle=False)
v=np.load(R/'run/COMPACT_VECTORS.npz',allow_pickle=False)
names=list(v['names']);nq=len(P['tasks']);errors=[];coef=[];quad=[]
lookup={int(k):j for j,k in enumerate(v['unique_ids'])}
partner=[]
for i,t in enumerate(P['tasks']):
    name=t['name'];z=np.load(R/'run'/name/'TRAIN_STATS.npz',allow_pickle=False)
    M=z['M'];b=z['b'];eps=len(t['episodes'])
    assert np.allclose(M,sum(z[f'e{k}_M'] for k in range(eps)),atol=1e-15,rtol=1e-12)
    assert np.allclose(b,sum(z[f'e{k}_b'] for k in range(eps)),atol=1e-15,rtol=1e-12)
    # Independent unscaled equivalent to diagonal-normalized ridge.
    scale=np.sqrt(np.maximum(np.diag(M),0));scale=np.maximum(scale,max(scale.max(initial=0)*1e-10,1e-15))
    theta=np.linalg.solve(M+P['case_R']['rho']*np.diag(scale*scale),b)
    coef.append(float(np.max(abs(theta-z['theta']))))
    y=v['unique_truth'][[lookup[int(k)] for k in v['query_ids'][i]]]
    cy=v['contrast']@y
    ps=np.sum((cy[None]-v['partner_predictions'][i])**2,axis=1);partner.append(ps)
    for j,n in enumerate(names):
        row=S['tasks'][name]['methods'][str(n)]
        errors.extend([abs(float(ps[j])-row['partner_SSE']),abs(float(ev['per_edge_SSE'][i,j].sum())-row['full_SSE'])])
    old=S['tasks'][name]['methods']['old12'];cr=S['tasks'][name]['methods']['case_R']
    for suf,metric in [('F','full_SSE'),('C','partner_SSE')]:
        calc=old[metric]-2*float(theta@ev[name+'_b'+suf])+float(theta@ev[name+'_M'+suf]@theta)
        quad.append(abs(calc-cr[metric]))
partner=np.stack(partner)
for key,idx in [('all5',list(range(nq))),('exclude_Q4',[0,1,2,4])]:
    ps=partner[idx].sum(0);fs=ev['per_edge_SSE'][idx].sum(axis=(0,2));rows=[]
    for edge in v['unique_ids']:
        entries=[ev['per_edge_SSE'][i,:,int(np.flatnonzero(v['query_ids'][i]==edge)[0])] for i in idx if edge in v['query_ids'][i]]
        if entries:rows.append(np.mean(entries,axis=0))
    eb=np.mean(rows,axis=0)
    for j,n in enumerate(names):
        want=S['aggregate'][key][str(n)]
        got=dict(pooled_partner_SSE=float(ps[j]),pooled_full_SSE=float(fs[j]),
            partner_gain_vs_M0_pct=float(100*(1-ps[j]/ps[0])),partner_gain_vs_old12_pct=float(100*(1-ps[j]/ps[1])),
            pooled_full_gain_vs_M0_pct=float(100*(1-fs[j]/fs[0])),edge_balanced_full_SSE=float(eb[j]),
            edge_balanced_full_gain_vs_M0_pct=float(100*(1-eb[j]/eb[0])),edge_balanced_full_gain_vs_old12_pct=float(100*(1-eb[j]/eb[1])))
        errors.extend(abs(got[k]-want[k]) for k in got)
assert max(errors)<1e-9 and max(coef)<1e-12 and max(quad)<1e-9
print(json.dumps(dict(status='PASS',coefficient_max_abs=max(coef),score_max_abs=max(errors),quadratic_full_partner_max_abs=max(quad),
    scope='All partner scores directly from saved vectors; full scores from per-edge losses and case_R quadratic statistics. Not raw-cell/full-feature reconstruction.'),indent=2))
