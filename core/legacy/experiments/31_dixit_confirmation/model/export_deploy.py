#!/usr/bin/env python3
"""Export selected existing evaluation fits, with no fitting or tuning."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse,datetime,hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    here=Path(__file__).resolve().parent;original=module('deploy_original',here/'predict_fixed.py');compact=module('deploy_compact',here/'predict_program.py')
    a.out.mkdir(parents=True,exist_ok=True)
    if (a.out/'EXPORT_COMPLETE.json').exists():raise FileExistsError('Use a fresh output directory')
    with np.load(a.source_run/'PREDICTIONS.npz',allow_pickle=False) as z:predictions={key:z[key] for key in z.files}
    report=[];checks=[];rng=np.random.default_rng(89231)
    for fold in range(3):
        source=a.source_run/f'FITTED_STATE_{fold}.npz';meta=json.loads((a.source_run/f'FOLD_{fold}.json').read_text());lam=float(meta['choices']['global'])
        with np.load(source,allow_pickle=False) as z:state={key:z[key] for key in z.files}
        # Training-only coefficients plus a fixed full-universe endpoint offset.
        weights=state['pt']@state[f'c_{lam}']/float(state['scale']) if float(state['scale']) else np.zeros_like(state[f'c_{lam}'])
        offset=(state['plus']@state['ktt'])@weights+state[f'beta_global_{lam}']
        target=a.out/f'Dixit_G2M_fold{fold}.npz'
        np.savez_compressed(target,format_version=np.array(1),nodes=state['nodes'],training_edges=state['et'],endpoint_gram=state['gram'],
            program_endpoint_effects=state['e'][:,state['pi']],program_sigma=state['sigma'][state['pi']],alpha=state['alpha'],
            selected_lambda=np.array(lam),kernel_coefficients=weights,endpoint_offset=offset)
        qp=predictions['pairs'][meta['query']].tolist();universe=state['universe'].tolist()
        pp=compact.predict(target,qp)['prediction'];truthfree_original=original.predict_from_state(source,universe,lam)['prediction'];whole=compact.predict(target,universe)['prediction']
        perm=rng.permutation(len(universe));order=compact.predict(target,[universe[i] for i in perm])['prediction'][np.argsort(perm)]
        single=np.concatenate([compact.predict(target,[pair])['prediction'] for pair in universe]);chunks=np.concatenate([compact.predict(target,universe[j:j+4])['prediction'] for j in range(0,len(universe),4)])
        swapped=compact.predict(target,['+'.join(pair.split('+')[::-1]) for pair in universe])['prediction']
        repeated=compact.predict(target,universe+[universe[0]])['prediction'][:-1]
        differences={'original_frozen_outer_predictions':float(np.max(abs(pp-predictions['global_selected'][meta['query']]))),
            'original_state_all45':float(np.max(abs(whole-truthfree_original))),'singleton':float(np.max(abs(single-whole))),
            'chunks4':float(np.max(abs(chunks-whole))),'order':float(np.max(abs(order-whole))),
            'endpoint_swap':float(np.max(abs(swapped-whole))),'duplicate_request':float(np.max(abs(repeated-whole)))}
        assert max(differences.values())<=1e-12,differences
        checks.append({'fold':fold,'outer_queries':len(qp),'identity_universe_queries':len(universe),'max_abs_differences':differences,'PASS':True})
        report.append({'file':target.name,'bytes':target.stat().st_size,'sha256':sha(target),'selected_lambda':lam,
            'source_fit_sha256':sha(source),'source_decision_sha256':sha(a.source_run/f'FOLD_{fold}.json'),
            'model_status':'outer-fold evaluation fit; retrospective development; not a confirmed production model',
            'reference_domain':'all unordered distinct pairs of the 10 stored measured-single endpoint identities',
            'source_selection':'training-only inner validation global full-program loss; no fitting or selection by export'})
    dump(a.out/'STATE_MANIFEST.json',report)
    receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS','total_npz_bytes':sum(r['bytes'] for r in report),
        'states':3,'outer_prediction_count':45,'all_domain_prediction_count':135,'tolerance':1e-12,'new_fitting_or_tuning':False,
        'query_labels_or_responses_read_for_inference':False,'compact_expression_fields':'only 10 by 52 measured-single endpoint effects and 52 scales',
        'export_script_sha256':sha(__file__),'predictor_script_sha256':sha(here/'predict_program.py'),'checks':checks}
    dump(a.out/'EXPORT_COMPLETE.json',receipt)
    (a.out/'README.txt').write_text('Dixit fixed-reference G2M evaluation models\n\n'
        'These three files are the existing outer-fold evaluation fits. They are\n'
        'retrospective development models, not independently confirmed production\n'
        'models, and have not been refitted on all 45 observations.\n\n'
        'Single-query inference from the package root:\n'
        'python model/predict_program.py --state model/deploy/Dixit_G2M_fold0.npz --pairs CREB1+YY1\n\n'
        'The selected regularization, measured-single node domain, program scales\n'
        'and full identity reference are fixed inside each file. Requesting one\n'
        'pair or many pairs yields the same prediction within floating precision.\n'
        'Only the 10 stored measured-single genes are supported; the two endpoints\n'
        'must be distinct. Input order and endpoint order do not affect the result.\n'
        'No observed query expression, labels or full RNA matrix is needed.\n')
    print(json.dumps({'status':'PASS','total_npz_bytes':receipt['total_npz_bytes'],'checks':checks}),flush=True)

if __name__=='__main__':main()
