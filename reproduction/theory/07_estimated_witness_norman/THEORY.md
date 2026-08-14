# Estimated Witness Risk: the deployable extension

## Why `k_t` cannot be learned from training doubles without an inductive map

Training double perturbations identify only the marginal covariance of their
own discrepancies, `K_train`.  For an unmeasured target, two joint covariance
models can have exactly the same `K_train` and training-data distribution but
different cross-covariances `k_t` and target variances `k_tt`, provided the
corresponding Schur complements remain nonnegative.  Their optimal target
predictors are different.  Therefore no estimator based on training outcomes
alone can identify `k_t` for a new target without an additional structural
assumption or a direct target witness.

This is the **cross-covariance identification boundary**.  In this experiment
the required inductive assumption is explicit and testable: pairs close in a
target-safe descriptor built from their two single-perturbation programs share
non-factorized discrepancy.

## Finite-panel Bayes Witness Risk

The exact unbiased theorem in experiment 06 assumes that every target
direction is estimable and that the factorized coefficient is fixed.  Norman
is sparse and rank deficient, so the deployable predictor uses the regularized
finite-panel form

\[
  \beta\sim N(0,B),\qquad
  y=X\beta+u+\varepsilon,
\]

with `Cov(u)=K` and `Cov(epsilon)=D`.  For target feature `t`, define

\[
C=XBX^\top+K+D,
\quad c_t=XBt+k_t,
\quad c_{tt}=t^\top Bt+k_{tt}+d_t.
\]

Then the posterior mean and its exact squared-loss risk are

\[
\widehat z_t=c_t^\top C^{-1}y,
\qquad
W_t=c_{tt}-c_t^\top C^{-1}c_t.
\]

The result follows directly from conditioning a joint Gaussian.  More
generally it is the minimum-MSE linear predictor for the stated second
moments.  When the prior variance in `B` tends to infinity on an estimable
row-space direction, this expression converges to the unbiased predictor in
experiment 06.  With `K=0`, it becomes the ordinary regularized factorized
virtual-cell baseline.

## Training-only estimator

For pair descriptor `psi(e)`, the discrepancy kernel is

\[
q_\ell(e,e')=\exp\{-\|\psi(e)-\psi(e')\|^2/(\ell s)\}.
\]

All descriptor standardization is fitted on the current training edges.  The
mixture fraction `rho`, length factor `ell`, and noise ratio `lambda` are
selected by nested edge holdout using only outer-training double outcomes.
If `a` is the training-only covariance scale estimate, then

\[
\widehat K=a\rho Q_{train,train},\qquad
\widehat k_t=a\rho Q_{train,t},\qquad
\widehat k_{tt}=a\rho.
\]

The pair descriptor uses PCA coordinates of the two single-perturbation
programs and the symmetric sum, absolute difference, and second-order outer
product.  A held-out target contributes no double-perturbation outcome to any
of these quantities.

## Falsifiable interpretation

- If estimated risk tracks oracle/realized risk and prediction improves, the
  training intervention graph contains transferable discrepancy geometry.
- If risk tracks but prediction does not improve, the estimator is a monitor
  but not a constructive AIVC component.
- If prediction improves but risk does not track, the new kernel is a useful
  predictor but not a valid reliability certificate.
- If both fail, the selected descriptor does not identify target
  cross-covariance; direct witnesses or a stronger structural map are needed.

