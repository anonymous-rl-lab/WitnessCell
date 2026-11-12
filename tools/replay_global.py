"""Replay v39 additions from integrated inputs and canonical manuscript sources (historical numerical module)."""
from pathlib import Path
import argparse,shutil,subprocess,sys,json
R=Path(__file__).resolve().parents[1];src=R/'experiments/global_reference'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--refit',action='store_true');a=ap.parse_args();out=a.out.resolve()
assert not out.exists() and not out.is_relative_to(R),'Use fresh output outside this release'
out.mkdir(parents=True)
for name in ['code','vendor','inputs','figure_inputs']:shutil.copytree(src/name,out/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
for name in ['CONFIGURATIONS.json','ANALYSIS_PLAN.json']:shutil.copy2(src/name,out/name)
shutil.copytree(R/'paper/manuscript_sources',out/'paper_sources')
if not a.refit:shutil.copytree(src/'fitted',out/'fitted')
(out/'results').mkdir();(out/'logs').mkdir()
def run(name,*args):
 with (out/'logs'/(Path(name).stem+'.log')).open('w') as log:subprocess.run([sys.executable,'-B',str(out/'code'/name),*map(str,args)],stdout=log,stderr=subprocess.STDOUT,check=True)
 print('PASS',name,flush=True)
if a.refit:run('refit_global.py')
for name in ['score_and_diagnose.py','han_modes.py','check_theory.py','build_tables.py']:run(name)
run('validate_additions.py','--expected',src/'results')
si=(R/'paper/manuscript_sources/WitnessCell_v41_Supplementary_Information.md').read_text();note=si.split('## Supplementary Note 24 |',1)[1].split('## Supplementary figure legends',1)[0];made=(out/'note24.md').read_text();rows=lambda t:[line.strip() for line in t.splitlines() if line.startswith('|')]
from collections import Counter
assert not (Counter(rows(made))-Counter(rows(note))), 'Regenerated matched-product table rows are absent from current S54-S57'
report={'status':'PASS','refitted_predictors':a.refit,'configurations':7,'all_new_table_rows_match_manuscript':True,'scope':'v39 matched-reference analysis, saved Han modes and bounded-model checks; no Han or transcriptomic model fitting.'}
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
