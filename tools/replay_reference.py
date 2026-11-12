"""Current v41 sensitivity and Han attribution from integrated assets."""
from pathlib import Path
import argparse,json,shutil,subprocess,sys
R=Path(__file__).resolve().parents[1];S=R/'experiments/reference_sensitivity'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--refit',action='store_true');a=ap.parse_args();out=a.out.resolve()
if out.exists() or out.is_relative_to(R) or out in R.parents:ap.error('Use a fresh output outside the repository')
out.mkdir(parents=True)
for folder in ['code','vendor','inputs','figure_inputs']:
 shutil.copytree(S/folder,out/folder,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
shutil.copy2(S/'ANALYSIS_PLAN.json',out/'ANALYSIS_PLAN.json')
(out/'results').mkdir();(out/'logs').mkdir()
for name in ['MATCHED_REFERENCE_POOLED.csv','GENE_LOADINGS.csv','MODE_DIAGNOSTICS.csv','MATCHED_REFERENCE_SCORES.csv','MODULE_GAIN_ATTRIBUTION.csv']:
 shutil.copy2(S/'results'/name,out/'results'/name)
def run(name,*args):
 with (out/'logs'/(Path(name).stem+'.log')).open('w') as f:
  subprocess.run([sys.executable,'-B',str(out/'code'/name),*map(str,args)],check=True,stdout=f,stderr=subprocess.STDOUT)
 print('PASS',name,flush=True)
if a.refit:run('refit_saturation.py')
else:shutil.copytree(S/'fitted',out/'fitted')
run('analyze_v40.py','--base',R)
run('validate_v40.py','--base',R,'--expected',S/'results')
subprocess.run([sys.executable,'-B',str(R/'tools/check_sensitivity_tables.py'),'--results',str(out/'results'),'--out',str(out/'DISPLAY_CHECKS.json')],check=True)
report=json.loads((out/'results/INDEPENDENT_CHECKS.json').read_text())
report.update(paper='v41',language_revision='R1',four_configurations_refitted=a.refit,han_predictors_refitted=False)
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
