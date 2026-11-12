"""v41 table/figure and diagnostic replay; isolate all generated outputs."""
from pathlib import Path
import argparse,shutil,subprocess,sys,os,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--task',choices=['paper_tables','diagnostics','attribution','figures'],required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out.resolve()
assert not out.exists() and not out.is_relative_to(R),'Use fresh output outside the release'
out.mkdir(parents=True);env=os.environ.copy();env['WITNESSCELL_LEGACY_ROOT']=str(R/'core/legacy');env['WITNESSCELL_HAN_SNAPSHOT']=str(R/'core/experiments/han2017/snapshot')
def call(script,*args):subprocess.run([sys.executable,'-B',str(script),*map(str,args)],check=True,env=env)
def copy(src,dst):shutil.copytree(src,dst,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
if a.task=='paper_tables':
 p=out/'paper';copy(R/'paper',p);call(p/'code/rebuild_tables.py');call(p/'code/verify_revision.py','--paper',p/'manuscript_sources');report=json.loads((p/'verification/V41_CHECKS.json').read_text())
elif a.task=='diagnostics':
 src=R/'experiments/diagnostics';call(src/'code/functional_diagnostics.py','--control',src/'inputs/control','--out',out/'replay');checks=[]
 for name in ['FUNCTIONAL_ORACLES.csv','FUNCTIONAL_REPEATABILITY.csv']:
  x=pd.read_csv(src/'results'/name);y=pd.read_csv(out/'replay'/name);pd.testing.assert_frame_equal(x,y,check_exact=False,rtol=1e-11,atol=1e-11);checks.append({'file':name,'rows':len(x),'passed':True})
 report={'status':'PASS','checks':checks}
elif a.task=='attribution':
 p=out/'analysis';copy(R/'experiments/diagnostics',p);call(p/'code/han_attribution.py','--han',p/'inputs/han','--out',p/'replay');call(p/'code/summarize.py','--results',p/'replay');call(p/'code/check_analysis.py');report=json.loads((p/'replay/VALIDATION.json').read_text())
else:
 p=out/'established';copy(R/'core/paper',p)
 call(p/'code/make_figures.py');call(p/'code/make_control_figures.py')
 call(R/'paper/code/make_figures.py','--out',out/'current')
 gp=out/'global_reference';copy(R/'experiments/global_reference',gp)
 call(gp/'code/make_figures.py','--paper',out/'current')
 sp=out/'reference_sensitivity';copy(R/'experiments/reference_sensitivity',sp)
 call(sp/'code/make_figures.py','--paper',out/'current')
 figs=out/'figures';figs.mkdir()
 for src in (R/'paper/reference_figures').glob('*.pdf'):
  name=src.name;prefix=name.split('_')[0]
  if prefix in ['Fig1','Fig4','SuppFig5','SuppFig6','SuppFig7']:
   found=list((out/'current').rglob(name))
  else:found=[p/('supplementary_figures' if name.startswith('Supp') else 'figures')/name]
  assert len(found)==1 and found[0].exists(),(name,found)
  shutil.copy2(found[0],figs/name)
  svg=found[0].with_suffix('.svg')
  if svg.exists():shutil.copy2(svg,figs/svg.name)
 # Visual equality against authoritative v41 figures, ignoring PDF timestamps.
 import fitz
 comparisons=[]
 for pdf in sorted(figs.glob('*.pdf')):
  with fitz.open(pdf) as x,fitz.open(R/'paper/reference_figures'/pdf.name) as y:
   assert len(x)==len(y)==1
   xx=x[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False);yy=y[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5),alpha=False)
   assert (xx.width,xx.height)==(yy.width,yy.height)
   xx=np.frombuffer(xx.samples,dtype=np.uint8).astype(float);yy=np.frombuffer(yy.samples,dtype=np.uint8).astype(float);mae=float(np.mean(abs(xx-yy)));ratio=float(np.mean(xx!=yy))
   # Recorded platform uses identical rendering; a changed platform must be inspected.
   assert mae<.05 and ratio<.002,(pdf.name,mae,ratio)
   comparisons.append({'file':pdf.name,'pixel_mean_absolute_difference':mae,'changed_channel_fraction':ratio})
 assert len(comparisons)==12
 report={'status':'PASS','main_figures':5,'supplementary_figures':7,'checks':comparisons,'scope':'Scientific data and redraw parity; no model fitting; PDF metadata need not be identical.'}
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
