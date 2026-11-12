#!/usr/bin/env python3
"""Behavioral checks for the panel-control interface and retained data splits."""
import os
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import itertools,json,tempfile,subprocess,sys
import numpy as np
from witnesscell_partner_control import incidence,project,predict_panel,fit_amplitude
R=Path(__file__).resolve().parents[1]
checks=[]
def check(name,value):
    if not bool(value):raise AssertionError(name)
    checks.append({'name':name,'passed':True})
def rejects(name,fn):
    try:fn()
    except ValueError:check(name,True);return
    raise AssertionError(name)

rng=np.random.default_rng(32);e=np.array(list(itertools.combinations(range(6),2)));b=rng.normal(size=(len(e),3));d=rng.normal(size=b.shape);y=rng.normal(size=b.shape)
out=predict_panel(b,b+d,e,6,.3)
check('alpha one recovers original',np.allclose(predict_panel(b,b+d,e,6,1)['partner_only'],b+d))
check('alpha zero preserves endpoint correction',np.allclose(predict_panel(b,b+d,e,6,0)['partner_only'],b+out['endpoint_correction']))
check('zero correction preserves baseline',np.array_equal(predict_panel(b,b,e,6,.5)['partner_only'],b))
check('endpoint-only correction is invariant',np.allclose(predict_panel(b,b+incidence(e,6)@rng.normal(size=(6,3)),e,6,0)['partner_correction'],0))
_,whole,_=project(e,6,out['whole_shrink']);_,component,_=project(e,6,out['partner_only'])
check('whole and component partner predictions agree',np.allclose(whole,component))
er0,pr0,_=project(e,6,y-b-d);er1,pr1,_=project(e,6,y-out['partner_only'])
check('endpoint error preserved',np.allclose(er0,er1))
check('full and partner SSE changes agree',np.isclose(np.sum((y-b-d)**2)-np.sum((y-out['partner_only'])**2),np.sum(pr0**2)-np.sum(pr1**2)))
perm=rng.permutation(len(e));rev=e[:,::-1]
check('query permutation equivariance',np.allclose(predict_panel(b[perm],(b+d)[perm],e[perm],6,.3)['partner_only'],out['partner_only'][perm]))
check('unordered pair orientation invariance',np.allclose(predict_panel(b,b+d,rev,6,.3)['partner_only'],out['partner_only']))
single=predict_panel(b[:1],(b+d)[:1],e[:1],6,0)
check('zero contrast dimension retains original',int(single['contrast_dim'])==0 and np.array_equal(single['partner_only'],(b+d)[:1]))
rejects('reversed duplicate pairs rejected',lambda:incidence(np.vstack([e,e[0,::-1]]),6))
rejects('self pairs rejected',lambda:incidence(np.array([[1,1]]),6))
rejects('invalid alpha rejected',lambda:predict_panel(b,b+d,e,6,1.1))
bad=b.copy();bad[0,0]=np.nan
rejects('nonfinite predictions rejected',lambda:predict_panel(b,bad,e,6,.3))

ia=np.load(R/'inputs/INNER_GROWTH_REPLAY.npz');oa=np.load(R/'inputs/CALIBRATION_OUTCOMES.npz')
keys=sorted(set(k.split('__')[0] for k in ia.files))
for key in keys:
    train={tuple(sorted(t)) for t in ia[key+'__edges'].tolist()};query={tuple(sorted(t)) for t in oa[key+'__edges'].tolist()}
    check('training-query pair separation: '+key,not train.intersection(query))

with tempfile.TemporaryDirectory() as td:
    td=Path(td);script=Path(__file__).with_name('witnesscell_partner_control.py')
    for cmd in [[sys.executable,str(script),'fit','--input',str(R/'example/training_oof.npz'),'--output',str(td/'calibration.json')],
                [sys.executable,str(script),'predict','--input',str(R/'example/query_predictions.npz'),'--calibration',str(td/'calibration.json'),'--output',str(td/'prediction.npz')]]:
        subprocess.run(cmd,check=True,capture_output=True)
    cal=json.loads((td/'calibration.json').read_text());pred=np.load(td/'prediction.npz')
    saved=np.load(R/'results/FOUR_ARM_PREDICTIONS.npz')
    check('real-data CLI equals all-config evaluation',np.allclose(pred['partner_only'],saved['HeLa_high__partner_only'],rtol=0,atol=1e-12))
    query=np.load(R/'example/query_predictions.npz');np.savez(td/'bad_query.npz',**{k:query[k] for k in query.files},truth=np.zeros_like(query['baseline']))
    run=subprocess.run([sys.executable,str(script),'predict','--input',str(td/'bad_query.npz'),'--calibration',str(td/'calibration.json'),'--output',str(td/'rejected.npz')],capture_output=True)
    check('query-label field refused by CLI',run.returncode!=0 and not (td/'rejected.npz').exists())
    # Copy the whole runnable package and execute through a new path.
    import shutil
    moved=td/'new_location';shutil.copytree(R/'software',moved/'software');shutil.copytree(R/'example',moved/'example')
    subprocess.run([sys.executable,str(moved/'software/witnesscell_partner_control.py'),'fit','--input',str(moved/'example/training_oof.npz'),'--output',str(moved/'cal.json')],check=True,capture_output=True)
    subprocess.run([sys.executable,str(moved/'software/witnesscell_partner_control.py'),'predict','--input',str(moved/'example/query_predictions.npz'),'--calibration',str(moved/'cal.json'),'--output',str(moved/'pred.npz')],check=True,capture_output=True)
    check('new-directory interface replay',np.array_equal(np.load(moved/'pred.npz')['partner_only'],pred['partner_only']))

report={'status':'PASS','checks':len(checks),'details':checks,'scope':'postprocessing interface and archived pair separation; does not certify original model-training provenance'}
(R/'results/INTERFACE_TESTS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks)}))
