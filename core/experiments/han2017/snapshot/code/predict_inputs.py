"""Frozen model inference, including explicitly labelled shared-function adapter."""
from common import *
import sys,importlib.util
sys.path.insert(0,str(ROOT/'frozen_engine/code'))
spec=importlib.util.spec_from_file_location('frozen_predict',ROOT/'frozen_engine/code/predict.py');engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
import torch

def shared_delta(s,family,genes,pairs,singles,knowledge):
    ix={g:i for i,g in enumerate(genes)};edges=np.array([[ix[a],ix[b]] for a,b in pairs]);scale=float(np.std(singles));su=singles*s.h/scale
    beta=s.meta[family]['native_beta']
    if beta==0:return np.zeros(len(pairs)),0.
    if family=='sf12':
        # Recover the fixed source linear embedding map, using source identities only.
        mean=s.knowledge.astype(float).mean(0);source_centered=s.knowledge.astype(float)-mean
        w=np.linalg.lstsq(source_centered,s.x,rcond=None)[0];replay=float(np.max(abs(source_centered@w-s.x)))
        assert replay<1e-8
        coords=(knowledge.astype(float)-mean)@w
        f=engine.sf_features(su,s.rms,coords,s.center,s.spread,edges)
        key=s.phase_sel['sf12_selected']['key'];raw=f@s.sf_state[key+'__theta']
    else:
        kk=knowledge.astype(np.float32).copy();kk-=s.knowledge.astype(np.float32).mean(0)
        kk/=np.maximum(np.sqrt(np.mean(kk*kk,axis=1,keepdims=True)),1e-8);kt=torch.tensor(kk)
        st=torch.tensor(su/s.rms,dtype=torch.float32);raw=np.zeros(len(pairs));replay=0.
        for model,state in s.networks:
            for start in range(0,len(pairs),2048):
                e=edges[start:start+2048];a,b=torch.tensor(e.T.copy(),dtype=torch.long)
                with torch.no_grad():v=model(kt[a],kt[b],torch.zeros((1,64)),st[a],st[b],torch.zeros(())).numpy()*state['scale']
                raw[start:start+len(e)]+=v/len(s.networks)
    return beta*raw*scale/s.h,replay

def main():
    if (ROOT/'data/INFERENCE_LOCK.json').exists():return
    sources={c:engine.Source(c) for c in ['HeLa','A549','293T']};checks={c:s.replay() for c,s in sources.items()}
    k=load(ROOT/'data/target_knowledge.npz');lookup={g:i for i,g in enumerate(k['genes'])};meta={};portable=[]
    for name in ['common_raw','common_published','full_raw','full_published']:
        inp=load(ROOT/f'data/{name}_inputs.npz');genes=inp['genes'];pairs=inp['pairs'];singles=inp['singles'];ix={g:i for i,g in enumerate(genes)}
        B=np.array([singles[ix[a]]+singles[ix[b]] for a,b in pairs]);models=[];bases=[];deltas=[];alphas=[]
        knowledge=np.array([k['embedding'][lookup[g]] for g in genes])
        if name.startswith('common'):
            P=Projection(pairs)
            for c,s in sources.items():
                byrep=[s.predict(genes,pairs,singles[:,r]) for r in range(singles.shape[1])]
                for family in ['rank1','spectral','sf12','v4']:
                    base=np.column_stack([p['M0'] for p in byrep]);native=np.column_stack([p[family+'__native_selected'] for p in byrep]);d=native-base
                    models.append(c+'__'+family);bases.append(base);deltas.append(d);alphas.append(s.meta[family]['alpha'])
                    if family in ['sf12','v4']:
                        direct=np.column_stack([shared_delta(s,family,genes,pairs,singles[:,r],knowledge)[0] for r in range(singles.shape[1])])
                        error=float(np.max(abs(P.P(direct)-P.P(d))));portable.append(dict(panel=name,source=c,family=family,projected_direction_max_abs_error=error));assert error<2e-6
        else:
            for c,s in sources.items():
                for family in ['sf12','v4']:
                    d=np.column_stack([shared_delta(s,family,genes,pairs,singles[:,r],knowledge)[0] for r in range(singles.shape[1])])
                    models.append(c+'__'+family);bases.append(B.copy());deltas.append(d);alphas.append(s.meta[family]['alpha'])
        payload=dict(pairs=pairs,genes=genes,B=B,M0=np.stack(bases,axis=2),D=np.stack(deltas,axis=2),models=np.array(models),source_alpha=np.array(alphas))
        np.savez_compressed(ROOT/f'data/{name}_frozen_predictions.npz',**payload)
        meta[name]=dict(pairs=len(pairs),models=models,readout_repeats=singles.shape[1],kind='original_frozen_models' if name.startswith('common') else 'shared_function_adapter_target_additive_baseline')
        print(json.dumps(dict(stage='frozen_inference_complete',panel=name,**meta[name])),flush=True)
    for rel,h in json.loads((ROOT/'MODEL_FREEZE.json').read_text())['files'].items():assert sha(ROOT/rel)==h
    dump(ROOT/'audit/PORTABLE_FEATURE_CHECKS.json',dict(source_replay=checks,common_projected_direction=portable,source_weights_unchanged=True,target_model_fits=0))
    dump(ROOT/'data/INFERENCE_LOCK.json',dict(created_utc=now(),query_labels_read=False,target_model_fits=0,panels=meta,predictions_sha256={p.name:sha(p) for p in (ROOT/'data').glob('*_frozen_predictions.npz')}))
if __name__=='__main__':main()
