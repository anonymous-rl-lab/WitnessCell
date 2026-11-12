# Component-specific control of a learned partner correction

## 1. Prediction object and response components

For a fixed panel of unordered query pairs Q, let X have one row per pair and one column per intervention, with ones at its two endpoints. Let H = XX⁺ and P = I − H. Both are orthogonal projectors under the unweighted Euclidean loss. X may be rank deficient; the Moore–Penrose inverse defines the same unique projection onto its column space.

Let M₀ denote an already fitted baseline, and let D denote a learned correction. Both are matrices with query pairs as rows and readouts or measured repeats as columns. Write D_E = HD and D_P = PD. The baseline itself can have a nonzero partner projection. This decomposition concerns the added correction, not a claim that M₀ is additive.

The original predictor is M₀ + D. Whole-correction attenuation gives M₀ + αD. Component-specific control gives

$$\widehat Y_\alpha=M_0+HD+\alpha PD.$$

Prediction uses query identities, M₀, D and a coefficient calculated from allowed training data. It uses no query outcomes. A query panel is part of the prediction specification: changing the panel can change a pair's endpoint/partner decomposition.

## 2. Exact preservation and error identity

**Proposition 1.** For every scalar α,

$$H(\widehat Y_\alpha-M_0)=HD,\qquad P(\widehat Y_\alpha-M_0)=\alpha PD.$$

**Proof.** H² = H, P² = P and HP = PH = 0. Apply H or P to HD + αPD. These identities hold for every outcome matrix and require no stochastic assumptions.

Let e = Y − M₀, e_E = He and e_P = Pe. Orthogonality gives

$$L_{\rm full}(\alpha)=\|e_E-D_E\|_F^2+\|e_P-\alpha D_P\|_F^2.$$

**Corollary 1.** Relative to the original predictor (α = 1),

$$L_{\rm full}(1)-L_{\rm full}(\alpha)=L_{\rm partner}(1)-L_{\rm partner}(\alpha).$$

The equality concerns unnormalized sums of squared error. Complete and partner percentage gains divide by different baseline losses. It does not imply that a training-selected α improves unseen outcomes. It also does not imply that component control dominates whole-correction shrinkage: shrinking the endpoint part can itself be beneficial on some panels.

Whole and component shrinkage have the same partner prediction because P(M₀ + αD) = PM₀ + αPD. They differ only in whether HD is attenuated. At α = 1 both recover the original predictor. At α = 0 component control retains M₀ + HD, whereas whole shrinkage returns M₀.

## 3. Why complete-response and partner objectives can disagree

For a scalar t applied to the entire D, the unconstrained least-squares optimum is

$$t_{\rm full}^*=\frac{\langle e_E,D_E\rangle+\langle e_P,D_P\rangle}{\|D_E\|_F^2+\|D_P\|_F^2}.$$

If both component energies are positive, define t_E* = ⟨e_E,D_E⟩/‖D_E‖² and t_P* = ⟨e_P,D_P⟩/‖D_P‖². Then

$$t_{\rm full}^*=\frac{\|D_E\|_F^2t_E^*+\|D_P\|_F^2t_P^*}{\|D_E\|_F^2+\|D_P\|_F^2}.$$

Thus optimizing complete response need not choose the amplitude preferred by the partner component. The formula is for unconstrained optima; averaging separately clipped component optima generally does not give the bounded complete optimum. The present implementation fixes the endpoint amplitude at its original value one and fits only the partner coefficient. It adds no second hyperparameter search.

## 4. Coefficient from internal held-combination predictions

For each archived inner fold k, let M₀,k and D_k be predictions from the permitted training combinations. Form r_k = P_k(Y_k − M₀,k) and d_k = P_kD_k using that fold's validation-pair graph. Pool

$$A=\sum_k\langle r_k,d_k\rangle,\quad B=\sum_k\|d_k\|_F^2,\quad E=\sum_k\|r_k\|_F^2.$$

The internal partner loss along this correction is E − 2αA + α²B. Its minimizer on [0,1] is

$$\alpha=\operatorname{clip}_{[0,1]}(A/B).$$

The executable rule takes α = 0 when B ≤ 10⁻¹⁴. Original zero-selected predictors remain zero corrections. This coefficient is fitted for the original selected model, not for an outer-performance-selected alternative.

The archived hyperparameters were selected using the same inner folds. Their held-combination predictions exclude each validation combination for fitting a fixed candidate, but the selected internal score can still be optimistic. This extension was designed after historical outer outcomes were known and is retrospective development. No confidence level or prospective false-activation guarantee is attached to α.

For E,B > 0, set ρ = A/√(BE). The unconstrained improvement along this fixed direction is 100ρ² percent of E. With t* = A/B ≠ 0 and k = α/t*, the achieved internal gain is 100ρ²(2k − k²). These are algebraic identities on the internal arrays. They are neither an independent predictive law nor an upper bound on external or dataset-wide learnability. If A ≤ 0, nonnegative control chooses zero; if B = 0, ratio-based diagnostics are undefined.

## 5. Two targeted checks

**Repeat consistency.** The main coefficient pools both measured repeats. Additional coefficients use repeat 1 or repeat 2 of the internal arrays separately. Every coefficient is evaluated on both outer repeats with the same fixed underlying models. The resulting 990 score rows are descriptive repeated measurements, not independent model fits or independent biological confirmations. Shared fitting and experimental effects remain possible.

**Query-panel dependence.** Query pairs are sorted by their canonical unordered endpoint indices. Alternating ranks define two or three subpanels, with no use of outcomes. The same α is applied separately within each subpanel; the assembled predictions are evaluated on the original full query panel and its original projector. The singleton limit returns the original prediction because each one-row graph has zero partner dimension. The number of controllable directions is the sum of the subpanel contrast dimensions. It does not measure a biological sample size.

Each block's change is orthogonal to that block's endpoint matrix; concatenating these changes remains orthogonal to the full endpoint matrix. Consequently the preservation identity still holds globally after blockwise control, although fewer directions may be adjustable. Splitting a panel is not equivalent to reordering it: reordering leaves predictions equivariant, whereas splitting changes the response space being controlled.

## 6. Scope of the metric and possible extensions

The implemented analysis uses a common unweighted squared-error axis. Positive readout weights can be handled by applying the same column scaling to every term. For row weights W, Euclidean P must instead be replaced by the W-orthogonal projection I − X(XᵀWX)⁺XᵀW, with all losses and inner products evaluated in the same W metric. This weighted extension is a mathematical statement, not an additional empirical result here. Rectangle losses can reweight overlapping contrast directions and do not share the full-loss identity automatically.

For a fixed readout matrix S, the more general operation M₀ + HD + (PD)S also preserves HD. A scalar αI is the implementation tested here. Training-derived RNA readout weights provide a related amplitude-control strategy, but the original RNA procedure M₀ + D diag(w) also changes HD and must retain its own identity and validation. No cross-modal efficacy is inferred from the algebra alone.

## 7. Scientific claim supported by this construction

The response decomposition identifies which part of an existing predictor is changed by a proposed adjustment. Training-side fitting supplies a practical coefficient; held-out outcomes assess whether it helps. The preservation theorem, the empirical utility of calibration and the portability of calibration are separate statements. The supplied comparison evaluates all 55 eligible functional configurations and retains the original predictions, favorable outcomes and adverse outcomes.
