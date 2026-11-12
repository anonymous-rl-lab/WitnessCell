"""Replay in a fresh directory without modifying the delivered evidence."""
from pathlib import Path
import argparse,json,shutil,subprocess,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workdir',type=Path,required=True);ap.add_argument('--rebuild-measurements',action='store_true');ap.add_argument('--check-inputs-only',action='store_true');args=ap.parse_args()
    dst=args.workdir.resolve()
    if dst.exists():raise RuntimeError('Choose a fresh replay directory')
    dst.mkdir(parents=True)
    for name in ['code','frozen_engine']:shutil.copytree(ROOT/name,dst/name,ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['MODEL_FREEZE.json','PROTOCOL.json','RUN_CODE_LOCK.json']:shutil.copy2(ROOT/name,dst/name)
    (dst/'audit').mkdir();(dst/'results').mkdir()
    if args.rebuild_measurements:
        src=dst/'source_measurements';shutil.copytree(ROOT/'source_measurements',src)
        han=src/'han';(han/'counts').mkdir()
        with np.load(src/'counts_compact.npz') as counts:
            for key,run in [('plasmid','SRR5258688'),('t0','SRR5258687'),('repeat1','SRR5258686'),('repeat2','SRR5258685')]:counts[key].astype('<u8').tofile(han/'counts'/f'{run}.counts.u64')
        (dst/'data').mkdir()
        subprocess.run([sys.executable,str(dst/'code/prepare.py'),'--han-dir',str(han),'--knowledge',str(src/'GenePT_fixed64.npz')],check=True)
        for p in (ROOT/'data').glob('*_inputs.npz'):
            with np.load(p) as a,np.load(dst/'data'/p.name) as b:
                assert a.files==b.files
                for key in a.files:np.testing.assert_array_equal(a[key],b[key])
        for p in (ROOT/'data').glob('*_outcomes.npz'):
            with np.load(p) as a,np.load(dst/'data'/p.name) as b:
                for key in a.files:np.testing.assert_array_equal(a[key],b[key])
        # Preserve the original receipt after numerical equality; a regenerated timestamp is not a data difference.
        shutil.copy2(dst/'data/PREPARED.json',dst/'audit/REBUILT_INPUT_RECEIPT.json')
        shutil.copy2(ROOT/'data/PREPARED.json',dst/'data/PREPARED.json')
        print('ALL_INPUT_ARRAYS_REBUILT_EXACTLY',flush=True)
    else:
        shutil.copytree(ROOT/'data',dst/'data',ignore=shutil.ignore_patterns('*_frozen_predictions.npz','INFERENCE_LOCK.json'))
    if args.check_inputs_only:return
    subprocess.run([sys.executable,str(dst/'code/predict_inputs.py')],check=True)
    subprocess.run([sys.executable,str(dst/'code/run.py')],check=True)
    subprocess.run([sys.executable,str(dst/'code/verify.py')],check=True)
    import pandas as pd
    a=pd.read_csv(ROOT/'results/POOLED_SCORES.csv');b=pd.read_csv(dst/'results/POOLED_SCORES.csv');keys=['panel','model','method','repeat']
    a=a.sort_values(keys).reset_index(drop=True);b=b.sort_values(keys).reset_index(drop=True)
    assert a[keys].equals(b[keys]);difference=float(np.max(abs(a[['full_sse','partner_sse']].to_numpy()-b[['full_sse','partner_sse']].to_numpy())))
    assert difference<1e-8
    (dst/'REPLAY_RESULT.json').write_text(json.dumps({'score_max_abs_difference':difference,'base_model_fits':0},indent=2)+'\n')
    print('REPLAY_SCORE_MAX_ABS_DIFFERENCE',difference,flush=True)
if __name__=='__main__':main()
