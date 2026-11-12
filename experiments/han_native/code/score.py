from common import *
import pandas as pd

def main():
    tasks=json.loads((ROOT/'TASKS.json').read_text());locks={}
    for t in tasks:
        o=ROOT/'runs'/t['panel']/t['baseline']/str(t['seed'])/str(t['fold']);lock=json.loads((o/'LOCK.json').read_text())
        for f,h in lock['files'].items():assert digest(o/f)==h,(o,f)
        assert lock['protocol_sha256']==digest(ROOT/'PROTOCOL.md')
        locks[str(o.relative_to(ROOT))]=digest(o/'LOCK.json')
    dump(ROOT/'ALL_PREDICTIONS_LOCKED.json',dict(utc=now(),tasks=len(tasks),locks=locks,query_scoring_started=False))
    rows=[];checks=[];selections=[]
    for t in tasks:
        o=ROOT/'runs'/t['panel']/t['baseline']/str(t['seed'])/str(t['fold']);p=loadnpz(o/'predictions.npz');y=loadnpz(ROOT/'prepared'/t['panel']/str(t['seed'])/str(t['fold'])/'scoring_only.npz');tr=loadnpz(ROOT/'prepared'/t['panel']/str(t['seed'])/str(t['fold'])/'training.npz');assert np.array_equal(y['query_ids'],p['query_ids']);y=y['y'];proj=Project(p['query_edges'],len(tr['genes']));sel=json.loads((o/'SELECTION.json').read_text())
        add=tr['singles'][p['query_edges'][:,0]]+tr['singles'][p['query_edges'][:,1]]
        preds={'additive':add,'baseline':p['B']}
        for f in ('rank1','spectral'):
            selections.append(t|{'family':f}|sel[f])
            for policy in ('original','whole_same','off','controlled','whole_full'):preds[f+'_'+policy]=p[f+'_'+policy]
        for name,pred in preds.items():
            e=y-pred;ep=proj.p(e);ref=y-p['B'];refp=proj.p(ref)
            for r in range(y.shape[1]):rows.append(t|dict(method=name,repeat=r+1,n=len(y),full_sse=float(np.sum(e[:,r]**2)),partner_sse=float(np.sum(ep[:,r]**2)),baseline_full_sse=float(np.sum(ref[:,r]**2)),baseline_partner_sse=float(np.sum(refp[:,r]**2))))
        for f in ('rank1','spectral'):
            e1=y-p[f+'_original'];ea=y-p[f+'_controlled'];df=float(np.sum(e1**2)-np.sum(ea**2));dp=float(np.sum(proj.p(e1)**2)-np.sum(proj.p(ea)**2));gap=abs(df-dp);assert gap<1e-8*max(1.,np.sum(e1**2))
            checks.append(t|dict(check='full_equals_partner_sse_change',family=f,max_abs=gap))
    out=ROOT/'results';out.mkdir(exist_ok=True);df=pd.DataFrame(rows);df.to_csv(out/'FOLD_REPEAT_SCORES.csv',index=False);pd.DataFrame(selections).to_csv(out/'SELECTIONS.csv',index=False);dump(out/'IDENTITY_CHECKS.json',checks)
    cols=['full_sse','partner_sse','baseline_full_sse','baseline_partner_sse'];group=['panel','baseline','method']
    for name,keys in [('POOLED',group),('BY_PARTITION',group+['seed']),('BY_REPEAT',group+['repeat'])]:
        s=df.groupby(keys)[cols].sum().reset_index()
        for met in ('full','partner'):s[met+'_gain_vs_own_baseline_pct']=100*(1-s[met+'_sse']/s['baseline_'+met+'_sse'])
        s.to_csv(out/(name+'.csv'),index=False)
    dump(out/'SCORE_RECORD.json',dict(utc=now(),tasks=len(tasks),rows=len(df),lock_sha256=digest(ROOT/'ALL_PREDICTIONS_LOCKED.json'),scope='fixed-procedure retrospective target refitting; no independent-study replication count from folds'))
    print(pd.read_csv(out/'POOLED.csv').to_string(index=False))
if __name__=='__main__':main()
