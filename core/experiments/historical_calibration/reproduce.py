"""Portable numerical replay in a temporary copy; no biological model training."""
from pathlib import Path
import tempfile,shutil,subprocess,sys,json,os
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parent
scripts=['regress_select.py','sensitivity.py','partner_component.py','verify.py','figures.py']
tables=['REGRESSION_PERFORMANCE.csv','HELD_PREDICTIONS.csv','MODEL_SELECTION_RECORDS.csv','SELECTION_SUMMARY.csv','SELECTION_BY_PANEL.csv','REGRESSION_SENSITIVITY.csv','PARTNER_COMPONENT_SCORES.csv','PARTNER_COMPONENT_BY_CONDITION.csv']
env=os.environ.copy()
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:env[k]='1'
with tempfile.TemporaryDirectory(prefix='witnesscell_gain_replay_') as td:
 work=Path(td)/'analysis';shutil.copytree(ROOT,work,ignore=shutil.ignore_patterns('__pycache__','*.log'))
 for script in scripts:
  p=subprocess.run([sys.executable,str(work/'code'/script)],env=env,capture_output=True,text=True)
  if p.returncode:raise RuntimeError(script+'\n'+p.stdout+'\n'+p.stderr)
 checks=[]
 for name in tables:
  a=pd.read_csv(ROOT/'results'/name);b=pd.read_csv(work/'results'/name);assert a.shape==b.shape and list(a.columns)==list(b.columns),name
  delta=0.
  for col in a:
   if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col]):
    x=a[col].to_numpy(float);y=b[col].to_numpy(float);assert np.allclose(x,y,rtol=1e-10,atol=1e-8,equal_nan=True),(name,col)
    finite=np.isfinite(x)&np.isfinite(y)
    if finite.any():delta=max(delta,float(np.max(abs(x[finite]-y[finite]))))
   else:assert a[col].fillna('<NA>').equals(b[col].fillna('<NA>')),(name,col)
  checks.append(dict(file=name,rows=len(a),max_numeric_difference=delta))
 audit=json.loads((work/'results/VERIFICATION.json').read_text())
 report={'status':'PASS','tables':checks,'independent_verification_checks':audit['checks'],'scope':'Refits historical regressions, reselects candidates, independently recalculates all 55 functional inner calibration coefficients and outer losses, and recomputes query-batch partner-only predictions. Uses supplied sufficient arrays, no original neural training or source matrix preprocessing.'}
 (ROOT/'results/PORTABLE_REPLAY.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
