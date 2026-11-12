from pathlib import Path
import subprocess,sys
p=Path(__file__).resolve().parent
for name in ['make_method_v36.py','make_han_figures_v36.py']:
 subprocess.run([sys.executable,str(p/name)],check=True)
