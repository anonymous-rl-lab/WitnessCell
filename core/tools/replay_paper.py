"""Recompute current paper statistics or redraw all current figures, without fitting."""
from pathlib import Path
import argparse,shutil,subprocess,sys,os,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--figures',action='store_true');a=ap.parse_args();out=a.out.resolve()
if out.exists() or out.is_relative_to(R):ap.error('Use a fresh output outside the release.')
out.mkdir(parents=True);p=out/'paper';shutil.copytree(R/'paper',p,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
env=os.environ.copy();env['WITNESSCELL_LEGACY_ROOT']=str(R/'legacy');env['WITNESSCELL_HAN_SNAPSHOT']=str(R/'experiments/han2017/snapshot')
scripts=['make_figures.py','make_control_figures.py','make_method_v36.py','make_han_figures_v36.py'] if a.figures else ['analyze_revision.py','analyze_applications.py','analyze_four_corners.py','analyze_v28.py']
for s in scripts:
 print(s,flush=True);subprocess.run([sys.executable,str(p/'code'/s)],env=env,check=True)
checks=[]
files=['BENCHMARK_DECOMPOSITION.csv','BENCHMARK_PAIRED_INTERVALS.csv','BENCHMARK_DATASET_ABSOLUTE.csv','RNA_ERROR_BUDGET.csv','RNA_ALL_READOUT_CONTRIBUTIONS.csv','RNA_ALL_RECTANGLE_GAINS.csv','PAIRPICK_ALL_28_TASKS.csv','PAIRPICK_ALL_TASK_SUMMARY.csv','JURKAT_FOUR_CORNERS_V28.csv','JURKAT_FOUR_CORNER_CONTRASTS_V28.csv','BASELINE_EXCLUSION_SENSITIVITY.csv']
if not a.figures:
 for name in files:
  ref=pd.read_csv(R/'paper/revision'/name);now=pd.read_csv(p/'revision'/name)
  pd.testing.assert_frame_equal(now,ref,check_exact=False,rtol=1e-12,atol=1e-12)
  checks.append(dict(file=name,rows=len(ref),passed=True))
else:
 for sub in ['figures','supplementary_figures']:
  expected=sorted(x.name for x in (R/'paper'/sub).glob('*.pdf'));actual=sorted(x.name for x in (p/sub).glob('*.pdf'));assert actual==expected,(sub,actual,expected)
  checks.append(dict(directory=sub,figures=len(actual),passed=True))
 for name in ['Fig4_Han_target_calibration.csv','SuppFig5_Han_estimates_ranges.csv','Han_target_strengths.csv']:
  aa=pd.read_csv(R/'paper/figure_source_data'/name);bb=pd.read_csv(p/'figure_source_data'/name);pd.testing.assert_frame_equal(aa,bb,check_exact=False,rtol=1e-12,atol=1e-12);checks.append(dict(file=name,passed=True))
 audit=json.loads((p/'FIGURE_AUDIT.json').read_text());assert len(audit)==10 and all(not x['text_outside_page'] for x in audit)
report=dict(status='PASS',mode='figures' if a.figures else 'paper_scores',checks=checks,base_model_fits=0,scope='Figures redraw from numerical sources; score replay uses frozen per-operation metrics and saved biological predictions, not original single-cell preprocessing or all competitor inference.')
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
