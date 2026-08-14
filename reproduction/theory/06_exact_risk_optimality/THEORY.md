# Exact geometry–adequacy risk and its optimal predictor

## Model

For observed interventions, let

\[
y=X\beta+u+\varepsilon,
\]

where `X` is the intervention feature/design matrix, `β` is the shared
factorized mechanism, `u` is representation discrepancy not captured by that
factorization, and `ε` is measurement noise. For an unmeasured target,

\[
z=t^\top\beta+u_\star.
\]

Assume only second moments

\[
\operatorname{Cov}(u+\varepsilon)=\Sigma\succ0,\qquad
\operatorname{Cov}(u+\varepsilon,u_\star)=k,\qquad
\operatorname{Var}(u_\star)=k_{\star\star}.
\]

`Σ` contains both measurement noise and training-side representation error;
`k` describes which observed discrepancies carry information about the target
discrepancy. No target outcome is used in these quantities.

## Theorem 1: exact minimum-risk unbiased predictor

Assume `X` has full column rank. Among all linear predictors `aᵀy` that are
unbiased for every `β`, hence satisfy `Xᵀa=t`, the unique minimum-MSE predictor
has

\[
a^\star=\Sigma^{-1}k+\Sigma^{-1}X
(X^\top\Sigma^{-1}X)^{-1}
(t-X^\top\Sigma^{-1}k).
\]

Its exact prediction risk is

\[
\boxed{
R_\star=
\underbrace{k_{\star\star}-k^\top\Sigma^{-1}k}_{A_\star:
\text{unexplained representation discrepancy}}
+
\underbrace{\tilde t^\top
(X^\top\Sigma^{-1}X)^{-1}\tilde t}_{G_\star:
\text{adequacy-conditioned geometry}}
}
\]

with

\[
\tilde t=t-X^\top\Sigma^{-1}k.
\]

Thus geometry and adequacy are neither two arbitrary scores to average nor two
heads to select with a hard switch. Representation covariance changes the
metric `Σ⁻¹`, supplies a residual uncertainty `A★`, and changes the effective
target direction from `t` to `t̃`.

### Proof

For any unbiased linear predictor,

\[
R(a)=a^\top\Sigma a-2a^\top k+k_{\star\star},
\quad X^\top a=t.
\]

The Lagrangian first-order condition is

\[
\Sigma a-k+X\lambda=0.
\]

Solving this equation together with the unbiasedness constraint gives `a★`
above. Because `Σ` is positive definite, `R(a)` is strictly convex, so this
stationary point is the unique global minimum. Substitution gives the boxed
risk. Equivalently, for every other unbiased `a`,

\[
R(a)-R(a^\star)=(a-a^\star)^\top\Sigma(a-a^\star)\ge0.
\]

This is an exact finite-sample optimality statement, not an asymptotic claim.
Joint Gaussianity is unnecessary for optimality among linear unbiased
predictors. Under a joint Gaussian model and the corresponding diffuse-prior
limit for `β`, the same predictor is the posterior mean and is optimal under
squared loss among all measurable predictors.

## Important limits

1. **Correct factorization:** `K=0`, hence `k=0` and `Σ=σ²I`:

   \[
   R_\star=\sigma^2t^\top(X^\top X)^{-1}t,
   \]

   which is the ordinary target-geometry/V-optimal risk.

2. **Independent pair-specific mismatch:** training discrepancies and target
   discrepancy are independent, so `k=0`, `Σ=(τ²+σ²)I`, and

   \[
   R_\star=\tau^2+(\tau^2+\sigma^2)
   t^\top(X^\top X)^{-1}t.
   \]

   The first term is irreducible without directly informative pair witnesses.
   More globally well-conditioned data cannot remove it.

3. **Correlated mismatch:** `k≠0`. Observed partner/hub discrepancies can
   predict part of `u★`; both the adequacy term and effective geometry change.
   This is the regime where a second structured discrepancy expert can create
   new predictive ability.

## Theorem 2: optimal target-conditioned intervention design

For candidate subset `S`, recompute the observed design and covariance and let
`R_t(S)` be Theorem 1's risk for target `t`. For nonnegative target weights
`w_t`, define

\[
S^\star\in\arg\min_{|S|\le b}\sum_{t\in\mathcal T}w_tR_t(S).
\]

If the model and covariance are correct, `S★` minimizes the exact expected
target squared loss among every candidate subset of budget `b`, when each
subset uses its own optimal predictor from Theorem 1. Exhaustive enumeration
therefore gives a finite global certificate. For budget one, choosing the
candidate with the largest exact reduction in the boxed risk is exactly
one-step Bayes-optimal.

This does **not** imply that a greedy multi-step algorithm is globally optimal,
or that an estimated covariance is oracle-correct. Those are separate
approximation and estimation questions to test empirically.

## Consequence for WitnessCell

The theoretically correct reliability head is the predictive risk `R★`, not a
fixed score average. A constructive WitnessCell needs:

1. a shared interaction feature map `X`;
2. an OOF-estimated discrepancy covariance `K` rather than only a scalar bias;
3. target risk computed with the coupled formula;
4. intervention selection minimizing the resulting target risk.

The immediate CPU experiments verify the formula under known `K`, then compare
the globally enumerated design against pure geometry and discrepancy-only
designs. Learning `K` from Norman remains a later, separate gate.

