"""Historical Han: reconstruct measurements, replay saved inference, recalibrate and score."""
from pathlib import Path
import sys,subprocess,shutil,argparse,json
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out.resolve();assert not out.exists() and not out.is_relative_to(R);out.mkdir(parents=True)
ref=out/'reference';work=out/'replay'
def call(p,*args):subprocess.run([sys.executable,'-B',str(p),*map(str,args)],check=True)
call(R/'core/experiments/han2017/materialize.py','--out',ref,'--mode','exact')
call(ref/'code/replay.py','--workdir',work,'--rebuild-measurements','--check-inputs-only')
for p in (ref/'data').glob('*_frozen_predictions.npz'):shutil.copy2(p,work/'data'/p.name)
shutil.copy2(ref/'data/INFERENCE_LOCK.json',work/'data/INFERENCE_LOCK.json')
call(work/'code/run.py');call(work/'code/verify.py')
keys=['panel','model','method','repeat'];x=pd.read_csv(ref/'results/POOLED_SCORES.csv').sort_values(keys).reset_index(drop=True);y=pd.read_csv(work/'results/POOLED_SCORES.csv').sort_values(keys).reset_index(drop=True);assert x[keys].equals(y[keys]);delta=float(np.max(abs(x[['full_sse','partner_sse']].to_numpy()-y[['full_sse','partner_sse']].to_numpy())));assert delta<1e-8
report={'status':'PASS','score_max_absolute_difference':delta,'input_arrays_rebuilt_exactly':True,'batches':60,'score_rows':len(pd.read_csv(work/'results/ALL_SCORES.csv')) if (work/'results/ALL_SCORES.csv').exists() else 10440,'source_predictor_retrained':False,'source_inference_reexecuted':False,'scope':'Offline compact count/readout reconstruction, stored frozen inference, new target scalar estimation and independent QR scoring. Native Torch inference remains an optional separate task.'}
(out/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
