Dixit fixed-reference G2M evaluation models

These three files are the existing outer-fold evaluation fits. They are
retrospective development models, not independently confirmed production
models, and have not been refitted on all 45 observations.

Single-query inference from the package root:
python model/predict_program.py --state model/deploy/Dixit_G2M_fold0.npz --pairs CREB1+YY1

The selected regularization, measured-single node domain, program scales
and full identity reference are fixed inside each file. Requesting one
pair or many pairs yields the same prediction within floating precision.
Only the 10 stored measured-single genes are supported; the two endpoints
must be distinct. Input order and endpoint order do not affect the result.
No observed query expression, labels or full RNA matrix is needed.
