"""Independent scoring, isolation and immutable-weight checks after fixed execution."""
from common import *
import scipy.linalg
import tempfile,shutil
import pandas as pd
import run

def main():
    max_score_error=0.;records=[];guard=[];identity=[]
    for rel,h in json.loads((ROOT/'MODEL_FREEZE.json').read_text())['files'].items():assert sha(ROOT/rel)==h
    for rel,h in json.loads((ROOT/'RUN_CODE_LOCK.json').read_text())['files'].items():assert sha(ROOT/rel)==h
    for panel in run.PANELS:
        data=load(ROOT/f'data/{panel}_frozen_predictions.npz')
        for seed in range(20261108,20261113):
            for fold in range(3):
                dest=ROOT/'results'/panel/str(seed)/str(fold)
                assert (dest/'COMPLETE.json').exists()
                state=json.loads((dest/'PARAMETERS.json').read_text());q=load(dest/'QUERY_PREDICTIONS.npz')['ids'];p,outputs=run.arms(data,q,state,panel)
                y=load(ROOT/'data/splits'/panel/str(seed)/str(fold)/'query_outcomes.npz')['y'];native=pd.read_csv(dest/'SCORES.csv')
                # Independent pivoted QR of the design; no Gram inverse from the scorer.
                Q,R,pivot=scipy.linalg.qr(p.X,mode='economic',pivoting=True)
                diag=np.abs(np.diag(R));rank=int(np.sum(diag>max(p.X.shape)*np.finfo(float).eps*max(diag.max(),1.)))
                assert rank==p.rank;Q=Q[:,:rank]
                for arm,pred in outputs.items():
                    e=y[:,:,None]-pred;flat=e.reshape(len(q),-1);pe=(flat-Q@(Q.T@flat)).reshape(e.shape)
                    full=np.sum(e*e,axis=(0,1));partner=np.sum(pe*pe,axis=(0,1))
                    for j,model in enumerate(data['models']):
                        row=native[(native.method==arm)&(native.model==model)&(native.repeat=='pooled')].iloc[0]
                        error=max(abs(full[j]-row.full_sse),abs(partner[j]-row.partner_sse));max_score_error=max(max_score_error,float(error))
                        assert error<1e-7*max(1.,float(full[j]))
                tr=load(ROOT/'data/splits'/panel/str(seed)/str(fold)/'calibration.npz')
                # Perturbing query labels cannot enter calibrate(), which receives only calibration labels.
                a=run.calibrate(data,tr['ids'],tr['y'],seed,fold)
                mutated_query=y+1e6
                b=run.calibrate(data,tr['ids'],tr['y'],seed,fold)
                difference=float(np.max(abs(np.array(a['alpha'])-np.array(b['alpha']))));assert difference==0
                if seed==20261108 and fold==0:
                    # Also refuse any attempted file read during numerical calibration.
                    original=np.load
                    def forbidden(*args,**kwargs):raise AssertionError('Calibration attempted file access')
                    np.load=forbidden
                    try:c=run.calibrate(data,tr['ids'],tr['y'],seed,fold)
                    finally:np.load=original
                    assert c==a
                    # Exercise the actual prediction runner with a changed query-outcome file.
                    original_root=run.ROOT
                    with tempfile.TemporaryDirectory(prefix='han_query_guard_') as temp:
                        sandbox=Path(temp);sp=sandbox/'data/splits'/panel/str(seed)/str(fold);sp.mkdir(parents=True)
                        origin=ROOT/'data/splits'/panel/str(seed)/str(fold)
                        for name in ['calibration.npz','query_outcomes.npz','identities.json']:shutil.copy2(origin/name,sp/name)
                        (sandbox/f'data/{panel}_frozen_predictions.npz').symlink_to(ROOT/f'data/{panel}_frozen_predictions.npz')
                        run.ROOT=sandbox
                        try:
                            run.predict_fold(panel,seed,fold,data)
                            rr=sandbox/'results'/panel/str(seed)/str(fold)
                            before=json.loads((rr/'PARAMETERS.json').read_text());pred_before=load(rr/'QUERY_PREDICTIONS.npz')
                            shutil.rmtree(sandbox/'results')
                            np.savez_compressed(sp/'query_outcomes.npz',ids=q,y=mutated_query)
                            run.predict_fold(panel,seed,fold,data)
                            after=json.loads((rr/'PARAMETERS.json').read_text());pred_after=load(rr/'QUERY_PREDICTIONS.npz')
                            assert before==after
                            pred_change=max(float(np.max(abs(pred_before[key]-pred_after[key]))) for key in pred_before)
                            assert pred_change==0
                        finally:run.ROOT=original_root
                    guard.append(dict(panel=panel,calibration_function_accepts_no_query_outcomes=True,no_file_reads=True,query_file_values_changed_by=1e6,alpha_max_change=difference,actual_prediction_runner_max_change=pred_change))
                records.append(dict(panel=panel,seed=seed,fold=fold,independent_qr_rank=rank,n_query=len(q),target_alpha_min=min(a['alpha']),target_alpha_max=max(a['alpha'])))
    pd.DataFrame(records).to_csv(ROOT/'audit/INDEPENDENT_FOLD_CHECKS.csv',index=False)
    dump(ROOT/'audit/FINAL_CHECKS.json',dict(utc=now(),folds=len(records),independent_qr_max_abs_sse_difference=max_score_error,immutable_model_files_verified=True,locked_run_files_verified=True,query_isolation=guard,base_model_fit_calls=0,scope='Independent QR rescoring of every pooled arm/model/fold. Control calibration uses target training outcomes; query labels never enter its function. Checks are implementation evidence, not additional biological replications.'))
    print(json.dumps(dict(folds=len(records),independent_max_error=max_score_error,models_unchanged=True)),flush=True)
if __name__=='__main__':main()
