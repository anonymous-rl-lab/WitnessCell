"""Score all frozen model families on the identical Han panel; no fitting or selection."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import argparse, datetime, hashlib, json
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
CELLS = ['HeLa', 'A549', '293T']
FAMILIES = ['rank1', 'spectral', 'sf12', 'v4']

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def projector(pairs):
    genes = sorted(set(pairs.ravel())); ix = {g:i for i,g in enumerate(genes)}
    X = np.zeros((len(pairs), len(genes)))
    for i,(a,b) in enumerate(pairs):
        X[i,ix[a]] = X[i,ix[b]] = 1
    U,s,_ = np.linalg.svd(X, full_matrices=False)
    rank = int(np.sum(s > max(X.shape)*np.finfo(float).eps*s[0]))
    return np.eye(len(pairs))-U[:,:rank]@U[:,:rank].T, rank

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', default='results'); args = ap.parse_args()
    out = ROOT/args.output
    if (out/'SCORES.csv').exists():
        raise RuntimeError('Do not overwrite scored output')
    lock = json.loads((out/'PREDICTIONS_LOCKED.json').read_text())
    for f,h in lock['outputs'].items(): assert sha(out/f) == h
    for f,h in lock['model_sha256'].items(): assert sha(ROOT/f) == h
    rows = []; comparisons = []; identities = []; deletion = []
    for readout in ('raw_repeats','published'):
        with np.load(out/f'{readout}_predictions.npz') as z: pred = {k:z[k] for k in z.files}
        with np.load(ROOT/'outcomes'/f'{readout}.npz') as z:
            y = z['outcomes']; assert np.array_equal(z['pairs'],pred['pairs'])
        pairs = pred['pairs']; P,rank = projector(pairs)
        slices = [('pooled',slice(None))]
        if y.shape[1] > 1: slices = [('repeat1',slice(0,1)),('repeat2',slice(1,2))]+slices
        for c in CELLS:
            keys = [k.split('__',1)[1] for k in pred if k.startswith(c+'__')]
            scores = {}
            for arm in keys:
                e = y-pred[c+'__'+arm]; pe = P@e
                for repeat,col in slices:
                    a = float(np.sum(e[:,col]**2)); b = float(np.sum(pe[:,col]**2)); n = e[:,col].size
                    scores[arm,repeat] = (a,b)
                    rows.append(dict(readout=readout,source_model=c,arm=arm,repeat=repeat,n_pairs=len(pairs),
                                     partner_dimensions=len(pairs)-rank,full_sse=a,partner_sse=b,
                                     full_mse=a/n,partner_mse=b/n))
            for family in FAMILIES:
                for arm in ('native_selected','partner_control','raw_selected_direction'):
                    key = family+'__'+arm
                    for repeat,col in slices:
                        v = scores[key,repeat]
                        for reference in ('M0','single_additive','rank1__native_selected','rank1__partner_control',
                                          family+'__partner_off',family+'__native_selected',family+'__full_calibrated_whole'):
                            r = scores[reference,repeat]
                            comparisons.append(dict(readout=readout,source_model=c,model=family,arm=arm,
                                repeat=repeat,reference=reference,full_gain_pct=100*(1-v[0]/r[0]),
                                partner_gain_pct=100*(1-v[1]/r[1]),full_sse_reduction=r[0]-v[0],partner_sse_reduction=r[1]-v[1]))
                original = pred[c+'__'+family+'__native_selected']
                ctrl = pred[c+'__'+family+'__partner_control']
                change = ctrl-original
                full = np.sum((y-original)**2)-np.sum((y-ctrl)**2)
                part = np.sum((P@(y-original))**2)-np.sum((P@(y-ctrl))**2)
                endpoint = float(np.max(abs(change-P@change)))
                gap = float(abs(full-part)); assert endpoint < 1e-7 and gap < 1e-7*max(1,abs(full),abs(part))
                identities.append(dict(readout=readout,source_model=c,model=family,endpoint_max_abs=endpoint,SSE_identity_error=gap))
                for g in sorted(set(pairs.ravel())):
                    keep = np.all(pairs != g,axis=1); pp,_ = projector(pairs[keep])
                    for arm in ('native_selected','partner_control'):
                        ep = pp@(y[keep]-pred[c+'__'+family+'__'+arm][keep])
                        eb = pp@(y[keep]-pred[c+'__M0'][keep])
                        eo = pp@(y[keep]-pred[c+'__rank1__native_selected'][keep])
                        for repeat,col in slices:
                            v = np.sum(ep[:,col]**2)
                            deletion.append(dict(readout=readout,source_model=c,model=family,arm=arm,repeat=repeat,
                                deleted_endpoint=g,n_pairs=int(keep.sum()),partner_gain_vs_M0_pct=100*(1-v/np.sum(eb[:,col]**2)),
                                partner_gain_vs_rank1_pct=100*(1-v/np.sum(eo[:,col]**2))))
    df = pd.DataFrame(rows); cmp = pd.DataFrame(comparisons)
    df.to_csv(out/'SCORES.csv',index=False); cmp.to_csv(out/'COMPARISONS.csv',index=False)
    pd.DataFrame(identities).to_csv(out/'ALGEBRA_CHECKS.csv',index=False)
    pd.DataFrame(deletion).to_csv(out/'ENDPOINT_DELETION.csv',index=False)
    # Convenient complete table: original native family outputs, separately from controls.
    keep = df['arm'].isin(['M0','single_additive']+[f+'__native_selected' for f in FAMILIES])
    main = df.loc[keep].copy()
    for ix,row in main.iterrows():
        ref = df[(df['readout']==row['readout'])&(df['source_model']==row['source_model'])&(df['repeat']==row['repeat'])&(df['arm']=='M0')].iloc[0]
        main.loc[ix,'full_gain_vs_M0_pct'] = 100*(1-row['full_sse']/ref['full_sse'])
        main.loc[ix,'partner_gain_vs_M0_pct'] = 100*(1-row['partner_sse']/ref['partner_sse'])
    main.to_csv(out/'MODEL_COMPARISON.csv',index=False)
    receipt = dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),score_rows=len(df),
        model_fit_calls=0,target_calibration_calls=0,notes='All families and repeats retained. Overlapping endpoint deletions are sensitivity analyses, not independent replications.',
        files={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    (out/'SCORING_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print(main[(main.readout=='raw_repeats')&(main.repeat=='pooled')][
        ['source_model','arm','full_mse','partner_mse','full_gain_vs_M0_pct','partner_gain_vs_M0_pct']].to_string(index=False))

if __name__ == '__main__': main()
