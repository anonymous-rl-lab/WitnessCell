#!/usr/bin/env python3
"""Rebuild numerical outputs in an isolated directory and compare saved tables."""
from pathlib import Path
import subprocess,tempfile,shutil,sys,json
import pandas as pd,numpy as np
R=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='witnesscell_control_') as td:
    q=Path(td)/'partner_control';shutil.copytree(R,q,ignore=shutil.ignore_patterns('__pycache__'))
    for script in ['run_analysis.py','replay_rna.py','software/test_partner_control.py','independent_check.py']:
        run=subprocess.run([sys.executable,str(q/script)],capture_output=True,text=True)
        if run.returncode:raise RuntimeError(script+'\n'+run.stdout[-2000:]+'\n'+run.stderr[-2000:])
    checks=[]
    for name in ['DECISIONS.csv','FOUR_ARM_SCORES.csv','REPEAT_DECISIONS.csv','REPEAT_SCORES.csv','QUERY_PANEL_SCORES.csv','RNA_FOUR_ARM_SCORES.csv']:
        a=pd.read_csv(R/'results'/name);b=pd.read_csv(q/'results'/name)
        pd.testing.assert_frame_equal(a,b,check_exact=True)
        checks.append({'file':name,'rows':len(a),'exact_dataframe_match':True})
    out={'status':'PASS','tables':checks,'scope':'fresh-directory scalar fitting, scoring and interface replay; no biological response-model retraining'}
    (R/'results/PORTABLE_REPLAY.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out))
