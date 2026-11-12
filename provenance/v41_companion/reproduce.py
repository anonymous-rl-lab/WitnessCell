"""v41 paper: unchanged v40 analyses; saved replay or four exact refits with --refit."""
from pathlib import Path
import argparse,shutil,subprocess,sys,json
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--base',type=Path,required=True,help='Extracted grn_repo_v30 root');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--refit',action='store_true');ap.add_argument('--render',action='store_true');a=ap.parse_args()
base=a.base.resolve();out=a.out.resolve()
subprocess.run([sys.executable,str(ROOT/'verify_manifest.py')],check=True)
assert (base/'experiments/han_native/inputs/full_raw_inputs.npz').exists(),'Choose the extracted grn_repo_v30 directory.'
if out==ROOT or ROOT in out.parents:raise SystemExit('Choose an output outside this supplement.')
if out.exists() and any(out.iterdir()):raise SystemExit('Output must be empty; existing results are preserved.')
out.mkdir(parents=True,exist_ok=True)
with (out/'base_preparation.log').open('w') as f:
    subprocess.run([sys.executable,'-B',str(base/'tools/prepare.py'),'--out',str(out/'base_preparation')],stdout=f,stderr=subprocess.STDOUT,check=True)
print('PASS base archive preparation (saved cache restoration; no model fitting)',flush=True)
for folder in ['code','vendor','inputs','figure_inputs','render','paper_sources']:
    shutil.copytree(ROOT/folder,out/folder,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(ROOT/'ANALYSIS_PLAN.json',out/'ANALYSIS_PLAN.json')
(out/'results').mkdir();(out/'logs').mkdir()
for name in ['MATCHED_REFERENCE_POOLED.csv','GENE_LOADINGS.csv','MODE_DIAGNOSTICS.csv','MATCHED_REFERENCE_SCORES.csv','MODULE_GAIN_ATTRIBUTION.csv']:
    shutil.copy2(ROOT/'results'/name,out/'results'/name)
def run(script,*params):
    with (out/'logs'/(Path(script).stem+'.log')).open('w') as f:
        subprocess.run([sys.executable,str(out/script),*map(str,params)],stdout=f,stderr=subprocess.STDOUT,check=True)
    print('PASS',script,flush=True)
if a.refit:run('code/refit_saturation.py')
else:shutil.copytree(ROOT/'fitted',out/'fitted')
run('code/analyze_v40.py','--base',base)
run('code/validate_v40.py','--base',base,'--expected',ROOT/'results')
run('code/make_figures.py','--paper',out/'paper_sources')
run('code/assemble_figures.py','--paper',out/'paper_sources')
if a.render:run('render/render_documents.py','--paper',out/'paper_sources')
(out/'REPLAY_STATUS.json').write_text(json.dumps({'status':'PASS','four_configurations_refitted':a.refit,'han_predictors_refitted':False,'base_code':'grn_repo_v30_paper_v39.zip','paper_version':'v41'},indent=2)+'\n')
print('Completed',out)
