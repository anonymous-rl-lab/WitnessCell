"""Render current manuscript with optional Pandoc/XeLaTeX dependencies."""
from pathlib import Path
import argparse,shutil,subprocess,sys,json
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out.resolve()
assert not out.exists() and not out.is_relative_to(R)
shutil.copytree(R/'paper/manuscript_sources',out)
subprocess.run([sys.executable,'-B',str(R/'paper/render/render_documents.py'),'--paper',str(out)],check=True)
subprocess.run([sys.executable,'-B',str(R/'experiments/reference_sensitivity/code/assemble_figures.py'),'--paper',str(out)],check=True)
(out/'RESULT.json').write_text(json.dumps({'status':'PASS','paper':'v41','scope':'Canonical text and vectors rendered; no scientific refitting.'},indent=2)+'\n')
