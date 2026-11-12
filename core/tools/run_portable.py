from pathlib import Path
import argparse,subprocess,sys,shutil,json
R=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser();ap.add_argument('--task',choices=['control','interfaces'],required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out.resolve()
if out.exists() or out.is_relative_to(R):ap.error('Choose a fresh output directory outside the release.')
out.mkdir(parents=True)
if a.task=='control':
 p=out/'partner_control';shutil.copytree(R/'paper/partner_control',p,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 subprocess.run([sys.executable,'-B',str(p/'reproduce.py')],check=True)
 report=json.loads((p/'results/PORTABLE_REPLAY.json').read_text())
else:
 p=out/'software';shutil.copytree(R/'paper/software',p,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for f in ['test_interface.py','test_growth_interface.py']:subprocess.run([sys.executable,'-B',str(p/f)],check=True)
 c=out/'partner_control';shutil.copytree(R/'paper/partner_control',c,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 subprocess.run([sys.executable,'-B',str(c/'software/test_partner_control.py')],check=True)
 report=dict(status='PASS',scope='RNA core synthetic interface checks, functional input validation and archived prediction parity, and partner-control input/operation tests. No biological response-model retraining.')
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
