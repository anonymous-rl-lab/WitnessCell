"""Re-fit both bounded development batches in a fresh directory."""
import argparse,shutil,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args();dest=Path(args.output).resolve()
if dest.exists():raise SystemExit('Use a new output directory; existing results are never overwritten.')
dest.mkdir(parents=True)
for folder in ['code','vendor','inputs','references']:
    shutil.copytree(root/folder,dest/folder,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['PLAN.json','PLAN_PHASE2.json','GENEPT_PROVENANCE.json']:shutil.copy2(root/name,dest/name)
(dest/'results').mkdir();(dest/'diagnostics').mkdir()
for script in ['fit_models.py','fit_recovered.py']:
    for cell in ['HeLa','A549','293T','Jurkat']:
        subprocess.run([sys.executable,str(dest/'code'/script),'--cell',cell],check=True)
subprocess.run([sys.executable,str(dest/'code/evaluate_frozen.py')],check=True)
subprocess.run([sys.executable,str(dest/'code/check_and_sensitivity.py')],check=True)
print('Re-fit and scoring completed:',dest)
