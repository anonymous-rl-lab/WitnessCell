"""Score only after all five task predictions are frozen."""
from run_five import *
def main():
    start=time.monotonic();spec=check_frozen();frozen=json.loads((OUT/'ALL_PREDICTIONS_FROZEN.json').read_text())
    for n,h in frozen['files'].items():assert sha(OUT/n)==h
    assert sha(OUT/'PROTOCOL.json')==frozen['protocol_sha256']
    fac=factory();scores={};sufficient={};contrasts=[];per_edge=[]
    c=np.array(spec['contrast']);p=np.outer(c,c)
    for t in spec['tasks']:
        folder=OUT/t['name']
        with np.load(folder/'PREDICTIONS.npz',allow_pickle=False) as z:pred=z['predictions'];h=z['h_case_R']
        y=fac.fixed_reference(t['query_ids']);err=y[None]-pred
        fs=np.sum(err*err,axis=(1,2));pe=np.einsum('e,meg->mg',c,err);ps=np.sum(pe*pe,axis=1)
        es=np.sum(err*err,axis=2);per_edge.append(es);contrasts.append(np.einsum('e,meg->mg',c,pred));methods={}
        for i,n in enumerate(NAMES):
            dr=np.einsum('e,eg->g',c,pred[i]-pred[0]);r=c@(y-pred[0]);den=np.linalg.norm(r)*np.linalg.norm(dr)
            methods[n]=dict(full_SSE=float(fs[i]),partner_SSE=float(ps[i]),
                full_gain_vs_M0_pct=float(100*(1-fs[i]/fs[0])),partner_gain_vs_M0_pct=float(100*(1-ps[i]/ps[0])),
                full_gain_vs_old12_pct=float(100*(1-fs[i]/fs[1])),partner_gain_vs_old12_pct=float(100*(1-ps[i]/ps[1])),
                partner_direction_cosine_vs_M0=float(r@dr/den) if den>1e-30 else None,per_query_full_SSE=es[i].tolist())
        ro=y-pred[1];ph=np.einsum('ab,bgp->agp',p,h)
        sufficient[t['name']+'_MF']=np.einsum('egp,egq->pq',h,h)
        sufficient[t['name']+'_bF']=np.einsum('egp,eg->p',h,ro)
        sufficient[t['name']+'_MC']=np.einsum('egp,egq->pq',ph,ph)
        sufficient[t['name']+'_bC']=np.einsum('egp,eg->p',ph,p@ro)
        scores[t['name']]=dict(query_ids=t['query_ids'],query_pairs=t['query_pairs'],methods=methods)
    es=np.stack(per_edge);unique=sorted({i for t in spec['tasks'] for i in t['query_ids']});aggregate={}
    for key,indices in [('all5',list(range(5))),('exclude_Q4',[0,1,2,4])]:
        ps=np.array([[scores[spec['tasks'][j]['name']]['methods'][n]['partner_SSE'] for n in NAMES] for j in indices]).sum(0)
        fs=np.array([[scores[spec['tasks'][j]['name']]['methods'][n]['full_SSE'] for n in NAMES] for j in indices]).sum(0)
        edgevalues=[]
        for edge in unique:
            loc=[(j,t['query_ids'].index(edge)) for j,t in enumerate(spec['tasks']) if j in indices and edge in t['query_ids']]
            if loc:edgevalues.append(np.stack([es[j,:,k] for j,k in loc]).mean(0))
        eb=np.mean(edgevalues,axis=0)
        aggregate[key]={n:dict(pooled_partner_SSE=float(ps[i]),pooled_full_SSE=float(fs[i]),
            partner_gain_vs_M0_pct=float(100*(1-ps[i]/ps[0])),partner_gain_vs_old12_pct=float(100*(1-ps[i]/ps[1])),
            pooled_full_gain_vs_M0_pct=float(100*(1-fs[i]/fs[0])),edge_balanced_full_SSE=float(eb[i]),
            edge_balanced_full_gain_vs_M0_pct=float(100*(1-eb[i]/eb[0])),edge_balanced_full_gain_vs_old12_pct=float(100*(1-eb[i]/eb[1]))) for i,n in enumerate(NAMES)}
    loo={}
    for leave in range(5):
        a=np.array([[scores[t['name']]['methods'][n]['partner_SSE'] for n in NAMES] for i,t in enumerate(spec['tasks']) if i!=leave]).sum(0)
        loo[spec['tasks'][leave]['name']]={n:float(100*(1-a[k]/a[0])) for k,n in enumerate(NAMES)}
    direction={n:dict(positive_blocks=sum(scores[t['name']]['methods'][n]['partner_gain_vs_M0_pct']>0 for t in spec['tasks']),
        median_partner_gain_pct=float(np.median([scores[t['name']]['methods'][n]['partner_gain_vs_M0_pct'] for t in spec['tasks']]))) for n in NAMES[1:]}
    decision=dict(pooled_old12_partner_positive=aggregate['all5']['old12']['partner_gain_vs_M0_pct']>0,
        exclude_Q4_old12_partner_positive=aggregate['exclude_Q4']['old12']['partner_gain_vs_M0_pct']>0,
        old12_full_nonworse=aggregate['all5']['old12']['edge_balanced_full_gain_vs_M0_pct']>=0)
    decision['all_prespecified_conditions']=all(decision.values())
    save(OUT/'SCORES.json',dict(scope=spec['scope'],tasks=scores,aggregate=aggregate,direction=direction,
        leave_one_task_out_partner_gain_vs_M0_pct=loo,decision=decision,scoring_seconds=time.monotonic()-start,
        independent_replicates=False,contrast_identity_rank=3,unique_query_edges=8))
    np.savez_compressed(OUT/'EVALUATION_STATS.npz',**sufficient,per_edge_SSE=es,
        partner_SSE=np.array([[scores[t['name']]['methods'][n]['partner_SSE'] for n in NAMES] for t in spec['tasks']]))
    np.savez_compressed(OUT/'COMPACT_VECTORS.npz',query_ids=np.array([t['query_ids'] for t in spec['tasks']]),
        unique_ids=np.array(unique),unique_truth=fac.fixed_reference(unique),partner_predictions=np.stack(contrasts),
        names=np.array(NAMES),genes=fac.genes,pairs=np.array(fac.pairs),contrast=c)
    print(json.dumps({'stage':'score','aggregate':aggregate,'direction':direction,'decision':decision},ensure_ascii=False),flush=True)
if __name__=='__main__':main()
