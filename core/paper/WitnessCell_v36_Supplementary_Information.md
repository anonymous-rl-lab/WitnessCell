# Supplementary Information

## WitnessCell predicts combinatorial perturbation responses through endpoint and partner learning

**Companion to the v36 main text and Methods.** This supplement supplies derivations, implementation details, complete comparisons and sensitivity analyses for the six Results sections. It preserves the identities and evaluation status of the experiments in the accompanying `grn_repo_v27` archive. The manuscript name **PairPick (R-G)** refers to the unchanged archived model `R_G`.

## Contents

1. Method identities, evidence populations and notation
2. Target identification by joint interventions
3. Released predictor and complete benchmark comparisons
4. Endpoint elimination and the scope of partner contrasts
5. RNA partner construction and matched evaluation
6. Readout shrinkage, concentration and sensitivity
7. Functional growth prediction and procedure transfer
8. Model-based risk and empirical selection
9. PairPick: response constraints, candidate pools and completed comparisons
10. Independent-source checks and acquisition results
11. Statistical units, original outcomes and amendments
12. Source correspondence and computational reproduction
13. Target-specific decomposition, correction budgets, complete application accounting and software choice
14. Evidence status, score provenance, complete decision scope and executable interfaces
15. Additional quantitative results
16. Partner corrections across evaluation settings
17. Component-specific control: derivation and objectives
18. Complete functional comparison and source relationships
19. Repeat consistency and fixed-query-panel dependence
20. RNA candidate replay, historical diagnostics and software
21. Target-side calibration of frozen partner functions in Han
22. Shared functional predictors: features, source fitting and portable inference

Tables S1–S52 are embedded alongside the corresponding notes. The accompanying source bundle contains full-precision tables, relevant protocol records, a claim-to-source index and the script used to assemble this document. Bibliographic numbers refer to the main-text reference list; selected references are repeated at the end for convenience.

## Supplementary Note 1 | Method identities, evidence populations and notation

### 1.1 Scientific and computational correspondence

WitnessCell is the response-prediction method developed in this study. Its evidence has several levels: a controlled causal system establishes an information property; the released transcriptomic implementation supplies benchmark performance; RNA and growth extensions test additional response components; and selection analyses test uses of the predicted responses. These levels share a scientific motivation but have distinct estimators, training sets and outcomes. Table S1 identifies the implementation responsible for each claim. In particular, the formal benchmark rank belongs to the released response predictor. It is not a measured rank of the later RNA partner head or PairPick ranking extensions.

### Supplementary Table S1 | Main-text claims and their experimental identities

| Main-text result | Evaluated implementation | Evidence population | Principal outcome and interpretation |
|---|---|---|---|
| Method construction | Released endpoint-transfer core; explicitly separate partner extensions | Defined design matrices and archived fitted outputs | Algebraic capability of each response component |
| Joint-intervention information | Candidate least squares and contrast GLS in a specified causal family | Four worlds, two degree-matched panels, 1,000 paired noise seeds | Target risk and class-conditional target confidence sets |
| Formal prediction performance | Released WitnessCell core | Four benchmark datasets, three public split seeds; 654 operation–seed records | First of sixteen eligible entries under the specified aggregation |
| RNA partner learning | Matched M0, M1 and reliability-weighted R | 55 pairs, 11 same-fit rectangles, 4,995 non-target readouts | Additional partner predictions on a model-development dataset |
| Functional prediction | Repeat-specific baseline plus shared rank-one or spectral correction | K562 development; procedure-frozen Jurkat; original Shen source evaluation | Original-unit growth error and endpoint-eliminated error |
| Component control | Original correction with its partner part multiplied by a training-derived coefficient | 55 functional configurations and supplementary RNA candidate replay | Exact endpoint preservation; retrospective held-out utility and sensitivity |
| Target-side partner calibration | Han additive baseline plus frozen SF12/V4 functions | Full-screen Han pair splits; all six portable configurations | Two SF12 configurations improve over unchanged and removed corrections |
| Risk selection | Analytic covariance predictor; separate empirical ranking of released predictions | Pair-disjoint calibration/evaluation for the analytic experiment; released split comparisons for empirical ranking | Accepted error at reported coverage, with implementation-specific assumptions |
| Partner recommendation | PairPick using R-G response predictions | GSE220974: 49 detected double identities, 28 selection tasks from eight physical anchors | Observed utility, regret and empirical-best hits |

The Wessels benchmark and the matched CaRPool RNA panel reuse a biological source. K562 and Jurkat are different cell contexts within the same growth-screen study. GSE220974, Dixit and Combi-seq add other sources, but their readouts, interventions and evaluation objectives differ. Accordingly, the number of analyses is not a count of independent biological confirmations.

### 1.2 Common notation and loss conventions

Let $T$ denote permitted training combinations and $Q$ held-out query combinations. Rows of $Y_T$ or $Y_Q$ are condition-level responses; columns are readouts. For RNA, a response is expressed on the experiment's retained control-centred axis. Horlbeck growth errors use control-centred log2 enrichment normalized by population doublings; Shen uses relative log2 count slopes per day. Each source retains its own units and baseline, and their absolute errors are not pooled. An intervention identity includes activation or repression mode when applicable.

The unsigned incidence matrix $X$ has one row per distinct-endpoint pair and one column per endpoint identity, with two unit entries in each row. $X^+$ denotes its Moore–Penrose pseudoinverse. The endpoint-eliminating projector is $P=I-XX^+$. A complete error is $\|\widehat Y-Y\|_F^2$; an endpoint-eliminated error is $\|P(\widehat Y-Y)\|_F^2$. Division by the stated number of coordinates, or by contrast dimensions and readouts, defines the corresponding reported mean. A projected *ratio* can instead divide by $\|PY\|_F^2$, as in the structural tables. These quantities must not be interchanged.

Except where a table explicitly reports a loss ratio, a gain is $100(1-L_{\mathrm{method}}/L_{\mathrm{reference}})$. Negative gain denotes increased error. Pooling squared errors before forming a ratio differs from averaging fold-wise percentage gains. Every table retains its original convention. The historical weighted benchmark loss, raw non-target RNA loss, rectangle loss, growth loss and selection regret measure different outcomes.

## Supplementary Note 2 | Target identification by joint interventions

### 2.1 Compatible mechanisms and target-specific identification

Let $\mathcal M$ be a declared mechanism class and $D_{\mathcal E}(m)$ the population observation law under measured interventions $\mathcal E$. The equivalence class of a mechanism is

$$[m]_{\mathcal E}=\{m'\in\mathcal M:D_{\mathcal E}(m')=D_{\mathcal E}(m)\}.$$

A target $q(m)$ is identified at this observation law exactly when it is constant on $[m]_{\mathcal E}$. This definition concerns the target, rather than identification of all parameters or edges. Adding an intervention can partition the equivalence class, but it helps this target only if the removed ambiguity changes $q$.

For the special linear case with population observations $D\theta$ and target $c^\top\theta$, identification holds if and only if $c\in\operatorname{row}(D)$. To prove sufficiency, write $c=D^\top a$ and obtain $c^\top\theta=a^\top D\theta$. For necessity, if $c$ is not in the row space, it has a nonzero component $v$ in $\ker D$ with $c^\top v\ne0$. Then $\theta$ and $\theta+v$ have identical observations and different targets. This is a fixed-family estimability statement [15]; biological model misspecification remains a separate issue.

### 2.2 Mechanism family and equal-budget panels

The source $X$ is Bernoulli with probability $p=1/2$. Mediators $A,\ldots,F$ equal the source unless directly knocked out; intervening on the source also sets its descendants to zero. The response has affine mediator terms and a direct-source term, together with at most one gate from

$$h_A=A\lor(B\land C),\quad h_B=B\lor(A\land C),\quad h_C=C\lor(A\land B).$$

The intercept is 0.7 and the six mediator coefficients are $(.25,-.15,.4,-.2,.3,-.1)$ in the generating system. The four worlds are direct-only, with $(\beta_X,\beta_h)=(1,0)$, and one world for each gate, with $(0,1)$. Estimation receives the intervention semantics and candidate family, while the coefficients and generating world are hidden. Fitted coefficients are not constrained to the generating signs.

Both panels include control and singles on $X,A,B,C,D,E,F$. Their four additional conditions are $AB,BC,DE,DF$ in the informative panel and $AD,BD,BE,CF$ in the comparison panel. The added endpoint-degree vector is $(1,2,1,2,1,1)$ in both. Each panel therefore has twelve conditions and the same marginal endpoint coverage. The target queries $AC,AF,CE,EF$ are excluded from both, and

$$\Delta=\mu_{AC}-\mu_{AF}-\mu_{CE}+\mu_{EF}.$$

There are 256 observations per condition, split equally by pre-intervention source state, and independent Gaussian response noise with known standard deviation one. Source-state stratification fixes the stratum counts; it is part of the simulation design. The budget of 3,072 observations is per panel and world, not the total across all Monte Carlo repetitions.

### 2.3 Estimator-independent lower bound

Under control, all permitted singles and the comparison doubles, each candidate gate has the same value as the source. The four worlds consequently induce an identical training law. Their targets are $0,-1/2,0,-1/2$. For any possibly randomized estimator $T$ measurable from the comparison observations, its equal-world risk is

$$\begin{aligned}
\bar R(T)&=\tfrac12\mathbb E[T^2]+\tfrac12\mathbb E[(T+\tfrac12)^2]\\
&=\mathbb E[(T+\tfrac14)^2]+\tfrac1{16}\ \ge\ \tfrac1{16}.
\end{aligned}$$

The expectations can be taken under one common training law because that law is identical in all four worlds. The constant estimator $T=-1/4$ attains the bound. The maximum risk over the four worlds is at least their average, so the same value is also a minimax lower bound on this finite subproblem. The result does not require every world's risk to exceed 0.0625; for example, always predicting zero has zero error in two worlds and error 0.25 in the other two. Repeating the comparison conditions with less noise does not remove the population ambiguity.

### 2.4 Information geometry and the strong-estimator reference

After eliminating the affine nuisance terms, the informative $AB$ and $BC$ measurements distinguish three candidate directions. With a consistent contrast orientation, their population gate contributions are proportional to $-p\beta_h(1,0)$, $-p\beta_h(1,1)$ and $-p\beta_h(0,1)$. Nonzero points on these three lines are distinct. The corresponding target is $-p\beta_h$ for gates A and C and zero for gate B. Detecting a gate therefore differs from identifying a nonzero target.

The primary estimator fits the 24 stratum means after multiplying each by $\sqrt{n_s}/\sigma$. Its candidate design contains an intercept, six mediators, source and gate. It has rank eight in the comparison panel and nine in the informative panel. Minimum-norm least squares uses relative pseudoinverse cutoff $10^{-12}$; equal-SSE ties within $10^{-10}$ follow gate order A, B, C. Within-stratum sums of squares recover observation-level SSE exactly, including the noise component shared by the candidate fits.

The GLS reference uses the two measured conditional interaction contrasts, scaled by $p$. Their covariance is

$$V=\frac{\sigma^2}{N}\begin{pmatrix}2&1\\1&2\end{pmatrix}.$$

The off-diagonal element arises from shared control and B-single terms. Treating these contrasts as independent would change their uncertainty. GLS uses their joint covariance to fit the candidate directions. Table S2 shows that this direct use of the informative measurements has target risk close to the primary procedure. The large panel difference therefore measures information supplied by the interventions; it is not an exclusive estimator advantage.

### Supplementary Table S2 | Equal-world target risk and paired estimator comparisons

| Estimator | Mean target MSE | Monte Carlo SE | 95% Monte Carlo interval |
| --- | --- | --- | --- |
| Candidate least squares, comparison panel | 0.065295 | 0.000127 | [0.065046, 0.065543] |
| Candidate least squares, informative panel | 0.004674 | 0.000219 | [0.004244, 0.005105] |
| Contrast GLS, informative panel | 0.005247 | 0.000450 | [0.004363, 0.006131] |
| Always predict zero | 0.125000 | 0.000000 | [0.125000, 0.125000] |

| Paired comparison | Risk difference | 95% Monte Carlo interval |
| --- | --- | --- |
| OLS_high minus OLS_low | -0.060620 | [-0.061103, -0.060137] |
| OLS_high minus contrast_GLS_high | -0.000573 | [-0.001326, 0.000181] |
| OLS_high minus zero_target | -0.120326 | [-0.120756, -0.119895] |

Each risk first averages the four worlds within a noise seed. Intervals are two-sided Student-$t$ Monte Carlo intervals across 1,000 paired seeds, with 999 degrees of freedom. They describe simulation variability, not variability across biological systems. The point-prediction risk includes every retained fit, including those with empty confidence sets.

### 2.5 Simultaneous candidate-family confidence sets

Whiten the stratum observations and project out the eight-dimensional nuisance column space. An orthonormal basis of the remaining gate span yields $Z\sim N(a_k\beta_h,I_2)$ in the informative panel. The ball

$$\{\mu:\|Z-\mu\|^2\le c\},\qquad c=\chi^2_{2,.95}=5.9914645,$$

contains the true projected mean with probability 0.95 under the stipulated model. For each candidate line, let $\widehat b_k=a_k^\top Z/\|a_k\|^2$ and $d_k^2=\|Z-a_k\widehat b_k\|^2$. If $d_k^2\le c$, its retained coefficients are

$$I_k=\left[\widehat b_k-\sqrt{\frac{c-d_k^2}{\|a_k\|^2}},\ \widehat b_k+\sqrt{\frac{c-d_k^2}{\|a_k\|^2}}\right].$$

Otherwise that candidate contributes an empty set. The target set is the union of $-pI_A$, $\{0\}$ when gate B is retained, and $-pI_C$. On the event that the ball contains the true mean, the true candidate coefficient and hence the true target are retained. Target coverage is thus at least 0.95 within this family. This construction includes all candidates; using only the best-fitting gate would invalidate this argument.

Coverage and zero exclusion use the exact union, not its convex hull. An empty union is a coverage failure and does not certify a nonzero target. The comparison panel leaves relevant gate coefficients unidentified and yields an unbounded target set. Table S3 reports all four informative-panel worlds. Gate B illustrates that a zero target can be well supported despite nonzero interaction signals in measured conditions.

### Supplementary Table S3 | Informative-panel world-specific estimation and confidence sets

| World | Target MSE | Mean query MSE | Target coverage % | Empty set % | Nonzero declaration % |
| --- | --- | --- | --- | --- | --- |
| direct_only | 0.006998 | 0.004706 | 99.0 | 0.0 | 1.0 |
| hA | 0.005715 | 0.003933 | 95.0 | 1.6 | 98.3 |
| hB | 0.000234 | 0.004038 | 99.0 | 1.0 | 0.0 |
| hC | 0.005752 | 0.004080 | 95.0 | 1.5 | 98.4 |

Negative lower limits in normal/Student-$t$ Monte Carlo intervals for a nonnegative risk, where present in the full source table, are untruncated uncertainty approximations rather than negative measured losses. Nonzero-certificate rates in gate worlds describe power; the direct-only and gate-B rates describe false nonzero declarations.

### 2.6 Budget calculation and model-class boundary

For the frozen geometry, the whitened nuisance-residualized gate Gram matrix divided by $N$ is

$$\frac1N H^\top M_XH=
\begin{pmatrix}
11/64&3/32&-5/64\\
3/32&3/16&3/32\\
-5/64&3/32&11/64
\end{pmatrix}.$$

At unit gate strength the minimum squared separation from an incorrect candidate line is $N/8$. Put $r=\sqrt{\chi^2_{2,.95}}$ and $s_N=\min\{r,\max(\sqrt{N/8}-r,0)\}$. The event $\|\varepsilon\|\le s_N$ both retains the true mean and excludes competing lines, apart from probability-zero boundary contacts. Its probability is $F_{\chi^2_2}(s_N^2)$. Table S4 records the prespecified menu; the first budget with this conservative probability at least 0.8 is 256. This calculation uses the declared unit-strength geometry and known Gaussian noise. It is not a universal sample-complexity bound for arbitrary gate magnitudes.

### Supplementary Table S4 | Geometry-based observation-budget choice

| Observations per condition | Panel observations | Minimum squared line separation | Conservative probability | Meets 0.8 target |
| --- | --- | --- | --- | --- |
| 64 | 768 | 8.0 | 0.06990 | No |
| 128 | 1536 | 16.0 | 0.70023 | No |
| 256 | 3072 | 32.0 | 0.95000 | Yes |

Allowing all three gates simultaneously changes the conclusion. On the informative training panel, the function

$$f=h_A+h_C-h_B-X$$

vanishes. On the four queries, its contrast is $-2p=-1$. Adding an arbitrary coefficient times $f$ preserves the training law while changing the target. The enlarged eleven-column design has rank ten, and its null direction is target-relevant. This counterexample establishes the boundary of the causal result. Neither the confidence construction nor the 92.84% risk reduction asserts identification under arbitrary Boolean mechanisms or recovery of a complete biological regulatory network.

The formal source is `formal_bridge/` in the SI bundle, restored from experiment 33. Its `THREE_GATE_DESIGN_ENUMERATION.json` concerns this three-candidate experiment. Earlier bridge variants and broader-class exploratory files retain their original identities in the code archive.

## Supplementary Note 3 | Released predictor and complete benchmark comparisons

### 3.1 Operational details beyond the main-text equations

Measured single effects are kept as observed condition means relative to control. The missing-endpoint head combines the equal-condition training-single mean, ten Jaccard-weighted Gene Ontology neighbours and a self-response anchor. Construction of a held-out single's features removes that endpoint from the fitting information. The numerical ridge is $10^{-8}$ and coefficients are clipped to $[0,2]$. The separate on-target amplitude extension uses generic self-response and sign-matched tail and RMS features, relative ridge 0.2 and coefficient bounds $[-1.5,1.5]$.

The amplitude extension is activated only when nested endpoint holdouts give positive one-sided 95% lower bounds for both relative all-gene MSE improvement and top-100 correlation improvement. Fewer than six usable endpoints retains the preceding head. The top-100 MSE diagnostic is recorded but is not a third activation condition. These empirical gates select this head within its development procedure; they are not a general guarantee for new endpoints.

For the released combination predictor, $S_\alpha(z)=z/(1+\alpha|z|)$ is applied coordinatewise to the sum of endpoint effects. Training residuals use the conservative endpoint head, whereas final query backbones use the final head. This distinction matters when an endpoint is unmeasured. The selected saturation and endpoint-transfer parameters are $\alpha\in\{0,.02,.05,.1,.2,.5\}$ and $\eta\in\{.01,.03,.1,.3,1\}$. Conditional on their validation correction $D$ and residual error $E$, the amplitude $\gamma=\operatorname{clip}_{[0,1]}\langle D,E\rangle/\|D\|_F^2$ minimizes validation squared error over that interval. Zero correction energy gives zero amplitude. Numerical ties favour smaller $\alpha$, then smaller $\eta$.

The released validation convention is preserved rather than retrospectively replaced by the nested fitting convention of later experiments. Similarly, a saturation comparator selected on its own is different from the backbone paired with the final residual. Attribution of a residual branch always uses its actual paired backbone.

### 3.2 Full eligible ranking and metric profile

Table S5 contains all sixteen eligible entries in the formal combination track of scPerturBench [9]. The primary operation ranking and the mean metric rank are distinct aggregations; the former defines the reported first-place position. The 654 operation–seed records include 540 double-perturbation records and 114 unseen-single records. This scope explains why the combination-only loss analysis has a different denominator. The double-only archive also gives different scores to `linearModel` and `trainMean`: 529 of 540 saved MSEs differ, with maximum absolute difference 0.0235. Thus inclusion of unseen singles does not explain their entire ranking difference. A branch-level correspondence is not evidence that historical scoring objects are identical; both official score rows are retained (Supplementary Note 14, Table S31).

### Supplementary Table S5 | Complete formal top-100 leaderboard

| Position | Method | Mean operation rank | Mean metric rank | Records |
| --- | --- | --- | --- | --- |
| 1 | WitnessCell | 3.802752 | 4.851682 | 654 |
| 2 | scouter | 5.411315 | 5.780581 | 654 |
| 3 | linearModel | 5.724771 | 5.964832 | 654 |
| 4 | baseReg | 6.709480 | 6.538226 | 654 |
| 5 | trainMean | 6.770642 | 6.582059 | 654 |
| 6 | GEARS | 7.651376 | 7.049949 | 654 |
| 7 | scGPT | 7.805810 | 7.096840 | 654 |
| 8 | scELMo | 7.915902 | 7.227574 | 654 |
| 9 | AttentionPert | 8.058104 | 7.317278 | 654 |
| 10 | GeneCompass | 8.366972 | 7.494393 | 654 |
| 11 | GenePert | 8.571865 | 7.653160 | 654 |
| 12 | CPA | 9.250765 | 7.976300 | 654 |
| 13 | biolord | 9.845566 | 8.483690 | 654 |
| 14 | scFoundation | 11.689602 | 9.425331 | 654 |
| 15 | baseMLP | 13.018349 | 10.216106 | 654 |
| 16 | baseControl | 13.694190 | 10.703874 | 654 |

### Supplementary Table S6 | Complete metric-specific positions in the formal track

| Method | PCC-delta | MSE | E-distance | symmetric KL | Wasserstein | Common-DEGs |
| --- | --- | --- | --- | --- | --- | --- |
| WitnessCell | 1 | 1 | 1 | 1 | 9 | 8 |
| scouter | 2 | 2 | 13 | 6 | 1 | 3 |
| linearModel | 3 | 3 | 2 | 2 | 10 | 5 |
| baseReg | 4 | 5 | 3 | 3 | 11 | 14 |
| trainMean | 5 | 4 | 4 | 4 | 12 | 9 |
| GEARS | 9 | 6 | 5 | 14 | 5 | 11 |
| scGPT | 6 | 7 | 16 | 7 | 2 | 7 |
| scELMo | 10 | 8 | 6 | 12 | 6 | 15 |
| AttentionPert | 8 | 9 | 7 | 15 | 7 | 12 |
| GeneCompass | 12 | 10 | 8 | 13 | 8 | 10 |
| GenePert | 7 | 11 | 9 | 8 | 13 | 13 |
| CPA | 13 | 13 | 14 | 5 | 3 | 4 |
| biolord | 11 | 12 | 10 | 9 | 14 | 6 |
| scFoundation | 14 | 14 | 15 | 16 | 4 | 16 |
| baseMLP | 15 | 16 | 12 | 10 | 16 | 1 |
| baseControl | 16 | 15 | 11 | 11 | 15 | 2 |

Positions in Table S6 rank the mean metric ranks in the source file; they are not raw metric values. The full-precision mean ranks are retained in the source CSV. WitnessCell's first-place aggregate position coexists with weaker Wasserstein and common-differential-expression positions. Distributional scores use the benchmark's adapter combining predicted means with training-control variation. They do not demonstrate a fitted conditional distribution of perturbed cell states.

Equal dataset weighting gives a WitnessCell top-100 mean operation rank of 3.873152 and first position. Aggregation across both recorded gene panels gives 3.396789, and combining both panels with equal dataset weighting gives 3.456997; both remain first. These are robustness summaries of the same predictions, not independent replications. Table S7 retains paired uncertainty against the closest aggregate comparators. Bootstrap units are dataset–perturbation clusters, carrying their repeated seed records together.

### Supplementary Table S7 | Paired cluster-bootstrap differences in mean operation rank

| Comparator | WitnessCell minus comparator | 95% cluster interval |
| --- | --- | --- |
| scouter | -1.608563 | [-2.015861, -1.195782] |
| linearModel | -1.922018 | [-2.276019, -1.559939] |
| baseReg | -2.906728 | [-3.259264, -2.543472] |

The source analysis uses 20,000 bootstrap draws across 384 dataset–perturbation clusters. A bootstrap fraction of one favouring WitnessCell is an empirical resampling frequency, not a posterior probability or a statement that the population $P$ value is exactly zero.

### 3.3 Absolute losses, readout masks and support amount

The combination-only weighted analysis compares with the released training-mean additive reference. This self-contained weighted analysis counts its named additive reference once. It does not substitute that reference for the official `linearModel` score row or establish equality of historical adapted scoring inputs (Table S31). The 24.87% improvement is not an improvement relative to the strongest neural model. Historical weighting uses absolute target-versus-control Welch statistics from one cell half, min–max scaling and squaring, followed by scoring against the other half. Target coordinates are included in that historical endpoint.

Later non-target controls remove the query's perturbed genes, renormalize weights, use full-cell condition means and geometrically aggregate seed-wise ratios. Their pair bootstrap carries repeated records of the same perturbation together. Their numbers therefore need not equal the historical weighted ratios. Table S8 preserves both sets of outcomes with their proper denominators.

The residual-library comparison holds the fitted endpoint system and backbone fixed while changing the amount of available combination residual information. Its quarter, half and full libraries answer how much that branch benefits from support. They do not measure de novo partner contrasts, because the branch remains inside the endpoint span. All nine recorded Wessels MSE order paths and six of nine weighted-MSE order paths are monotone; the mixed weighted result is retained.

### Supplementary Table S8 | Absolute-error controls and residual-library contribution

These comparisons use earlier endpoint or analytic-risk implementations and are not additional partner-learning replications.

**Released benchmark and combination-only weighted analysis**

| Outcome | Value and scope |
| --- | --- |
| Released aggregation | Mean top-100 operation rank 3.803; first among 15 published entrants plus WitnessCell |
| Benchmark population | 654 condition–seed units, including 540 combinations and 114 unseen singles |
| Combination-only WMSE ratio vs released additive training mean | 0.7513 [0.7031, 0.8019] |
| Combination-only weighted direction-$R^2$ increment | +0.1908 [+0.1648, +0.2171] |
| All 654-unit WMSE ratio | 0.6829 [0.6255, 0.7434]; includes unseen singles |

**All-query non-target loss ratios**

| Dataset | Comparator | MSE ratio [95% pair interval] | WMSE ratio [95% pair interval] |
| --- | --- | --- | --- |
| Norman | Released training mean | 0.8701 [0.8304, 0.9070] | 0.9439 [0.8972, 1.0001] |
| Norman | Shared-endpoint additive | 0.9049 [0.8563, 0.9475] | 0.9996 [0.9464, 1.0622] |
| Norman | Independent saturation | 0.9355 [0.9023, 0.9669] | 0.9612 [0.9256, 0.9981] |
| Norman | Paired factorized backbone | 0.9392 [0.9083, 0.9685] | 0.9677 [0.9377, 0.9997] |
| Wessels | Released training mean | 0.6459 [0.6036, 0.6860] | 0.6976 [0.5650, 0.8610] |
| Wessels | Shared-endpoint additive | 0.6459 [0.6036, 0.6860] | 0.6976 [0.5650, 0.8610] |
| Wessels | Independent saturation | 0.6764 [0.6318, 0.7180] | 0.6185 [0.5184, 0.7330] |
| Wessels | Paired factorized backbone | 0.6675 [0.6238, 0.7084] | 0.6234 [0.5228, 0.7390] |

**Fixed-system Wessels residual-library contribution against the same backbone**

| Loss | Quarter library gain % | Half library gain % | Full library gain % | Full-minus-quarter percentage points |
| --- | --- | --- | --- | --- |
| MSE | 19.10 [16.17, 22.32] | 27.42 [23.81, 31.22] | 33.25 [29.28, 37.49] | 14.15 [12.42, 16.06] |
| WMSE | 26.81 [17.07, 34.52] | 33.91 [22.99, 42.74] | 37.66 [26.18, 47.74] | 10.85 [7.52, 14.36] |

Sources: `publication/released_benchmark/`, `publication/Fig2_prediction_ratios.csv`, `publication/ED1_coverage_comparisons.csv` and the `publication/Fig3_ED4_witness_*` tables in the SI bundle. The broader benchmark inventory is not substituted for the sixteen-entry eligibility rule.

## Supplementary Note 4 | Endpoint elimination and the scope of partner contrasts

### 4.1 Endpoint-transfer cancellation

For the released residual branch,

$$\delta_Q=\gamma\frac{X_QX_T^\top}{2}\left(\frac{X_TX_T^\top}{2}+\eta I\right)^{-1}R_T=X_Q\Theta,$$

where the final equality defines $\Theta$. Since $P_QX_Q=0$, $P_Q\delta_Q=0$ for every residual matrix and every selected amplitude. This is a structural identity independent of whether the training residuals contain interaction effects. Those effects can influence the estimated endpoint coefficients, but this output branch still lies in the query endpoint span.

The nonlinear backbone can have $P_QF_{\mathrm{core}}\ne0$. Consequently, cancellation of the residual branch does not imply that the complete predictor is additive. It identifies which component accounts for any new contrast. Across the twelve released splits, projected-error ratios of the final predictor and its paired backbone differed by at most $6.34\times10^{-8}$, consistent with the algebraic identity within saved-array precision. The recorded quantity is the change in projected-error ratio, not a norm of the difference between projected prediction arrays.

The identity distinguishes a particular predictor branch from the whole fitted method. Systema permits custom reference responses, with a perturbed centroid as a principal example, to assess perturbation-specific content.[11] The present operator removes arbitrary additive effects of each endpoint across the query pairs. These are related reference choices with different scientific targets. The method additionally specifies a learner capable of changing responses outside the endpoint span and tests its held-out increment.

### 4.2 Orthogonal error decomposition and weighted losses

The matrix $P_Q$ is symmetric and idempotent. For $E=\widehat Y-Y$,

$$\|E\|_F^2=\|(I-P_Q)E\|_F^2+\|P_QE\|_F^2,$$

because $P_Q(I-P_Q)=0$. Complete-response improvement may therefore occur wholly in the first term. If a common diagonal readout scaling is applied to every query, the decomposition also holds after that column scaling. Query-dependent element weights generally change the relevant inner product, and the same unweighted projector need not supply a weighted orthogonal decomposition. The structural tables use their stated common readout axis; the historical weighted benchmark uses a different weighting and is reported separately.

For a baseline $M_0$, let $e=Y-M_0$ be its error and $D$ an added correction, evaluated on the same query pairs and readout axis. Expanding the squared norm gives

$$\|P_Qe\|_F^2-\|P_Q(e-D)\|_F^2=2\langle P_Qe,P_QD\rangle_F-\|P_QD\|_F^2.$$

An actual change in the partner prediction requires $P_QD\ne0$. It lowers the measured partner error exactly when the right-hand side is positive: alignment with the baseline error must outweigh the squared correction magnitude. This is an identity on the evaluated arrays, not a guarantee that a training-selected correction will generalize. For a rectangle operator $C$, the same expansion holds with $C$ replacing $P_Q$; the reported rectangle and projector losses retain their different weighting of contrast directions.

An explicit partner learner need not output a pure endpoint-eliminated vector. Its correction decomposes as $D=(I-P_Q)D+P_QD$. The RNA kernel's projection is defined on the training pairs, and its continuation to query pairs need not lie entirely in the image of $P_Q$. Likewise, a signed low-rank growth matrix can change both components. Complete-response and projected losses therefore evaluate different contributions of the same fitted correction. The cancellation $P_Q\delta_Q=0$ remains exact for the core endpoint branch. Supplementary Note 17 uses this decomposition operationally: the endpoint part of an added correction is retained while its partner part receives a separately fitted strength.

### 4.3 Query-graph dimension

Take the graph whose vertices are active endpoint identities and whose edges are query pairs. For a connected component, $Xv=0$ requires $v_g=-v_h$ along every edge. A connected bipartite component admits one alternating-sign null vector; a connected non-bipartite component contains an odd cycle, forcing that vector to zero. If the graph has $n$ active vertices, $m$ edges and $b$ bipartite connected components, then

$$\operatorname{rank}(X)=n-b,\qquad \operatorname{rank}(P)=m-n+b.$$

A tree has no contrast dimension. An isolated odd cycle also has none, whereas an even cycle has one. More complex graphs can have endpoint-eliminated directions without an isolated four-edge rectangle; the rank formula, rather than counting named rectangles, defines the full contrast space. The twelve benchmark graphs have 212 such dimensions. These are linear-algebra dimensions of query sets, not 212 biological replicates.

### 4.4 Rectangles and full projection

For four pairs $AB,AD,CB,CD$, put $c=(1,-1,-1,1)^\top$. Then $c^\top X=0$. If these four edges form the entire query graph, $P=cc^\top/4$ and $\|PE\|_F^2=\|c^\top E\|^2/4$. With multiple or overlapping rectangles, summing their contrast losses weights directions according to the chosen rectangle list, while the orthogonal projector weights an orthonormal basis. The two reported endpoints are therefore complementary, rather than numerically interchangeable.

All four predictions in a rectangle must come from the same fitted model. Mixing predictions from different fits could introduce fitting differences into the contrast. Query projection uses identities to define the operator and observed outcomes only for evaluation. It is not used to tune a held-out correction. Both the true response and its measurement noise survive projection; endpoint removal is not itself a denoising operation.

### Supplementary Table S9 | Complete released-split structural decomposition

The zero-contrast reference is one. All-gene scoring is shown. E is the number of combination–seed evaluation records; d is the endpoint-eliminated dimension within that split. Joint backbone denotes the backbone paired with released core, whereas saturation-only was independently selected. A dash denotes an unavailable matching prediction, not an inferred score.

| Dataset / seed | E / d | Released core | Joint backbone | Saturation-only |
| --- | --- | --- | --- | --- |
| Norman / 1 | 73 / 24 | 0.995140 | 0.995140 | 0.995140 |
| Norman / 2 | 66 / 16 | 1.000736 | 1.000736 | 1.003976 |
| Norman / 3 | 69 / 22 | 1.018454 | 1.018454 | 1.009911 |
| Wessels / 1 | 63 / 39 | 1.000000 | 1.000000 | 0.999904 |
| Wessels / 2 | 54 / 33 | 0.993565 | 0.993565 | 0.993565 |
| Wessels / 3 | 64 / 45 | 1.000000 | 1.000000 | 0.997968 |
| Schmidt / 1 | 22 / 8 | 0.999722 | 0.999722 | — |
| Schmidt / 2 | 19 / 4 | 0.997796 | 0.997796 | — |
| Schmidt / 3 | 19 / 5 | 0.996738 | 0.996738 | — |
| Replogle_exp6 / 1 | 35 / 10 | 1.000000 | 1.000000 | — |
| Replogle_exp6 / 2 | 29 / 4 | 0.998739 | 0.998739 | — |
| Replogle_exp6 / 3 | 27 / 2 | 0.998246 | 0.998246 | — |

The six available Norman/Wessels additive references equal one. The complete component table and additional scoring variants remain in the compact archive; the separate matched RNA comparison is Table S10.

**Dataset-level complete-response and projected effects**

Within each dataset, squared errors and truth energies are summed across the three seeds before forming ratios. Gains compare released core with its own jointly selected backbone; projected ratios use projected observed-response energy as denominator.

| Dataset | Records / d | Readouts | Full-response gain % | Projected ratio |
| --- | --- | --- | --- | --- |
| Norman | 208 / 62 | 5025 | 5.7433 | 1.005527 |
| Wessels | 181 / 117 | 5020 | 31.9396 | 0.998160 |
| Schmidt | 60 / 17 | 5005 | 0.0658 | 0.998415 |
| Replogle_exp6 | 91 / 16 | 5019 | 1.2342 | 0.999490 |

Projected gains against the paired backbone are zero within saved-array quantization. Wessels uses 5,020 raw readouts here, versus 4,995 non-target readouts in Table S10.

## Supplementary Note 5 | RNA partner construction and matched evaluation

### 5.1 A positive semidefinite pair kernel

For permitted single-response vectors $e_g\in\mathbb R^p$, write $L_{gh}=e_g^\top e_h/p$ and

$$K_{(g,h),(u,v)}=\tfrac12(L_{gu}L_{hv}+L_{gv}L_{hu}).$$

Define the symmetric tensor feature $\phi(g,h)=(e_g\otimes e_h+e_h\otimes e_g)/(2p)$. Its inner product with $\phi(u,v)$ is exactly $K_{(g,h),(u,v)}$. The kernel is therefore positive semidefinite and invariant to reversing either pair. The construction represents relationships between combinations through measured single-response similarity; it does not require a learned molecular interaction graph.

Training residualization gives $K_T^\perp=P_TK_{TT}P_T$. For any vector $a$, $a^\top K_T^\perp a=(P_Ta)^\top K_{TT}(P_Ta)\ge0$. The projected kernel thus remains positive semidefinite. Its continuation to queries uses

$$K_{QT}^\perp=(K_{QT}-X_QX_T^+K_{TT})P_T.$$

For a training row used as its own query this reduces to $P_TK_{TT}P_T$. The subtraction extends the training endpoint adjustment through the query endpoint identities; it does not project using held-out responses. Both blocks are divided by $s=\operatorname{tr}(K_T^\perp)/|T|$. A degenerate $s$ returns zero correction.

For baseline residual $R_T=Y_T-M_{0,T}$, the correction is the kernel-ridge prediction

$$G_Q=\frac{K_{QT}^\perp}{s}\left(\frac{K_T^\perp}{s}+\lambda I\right)^{-1}R_T.$$

Training endpoint components cannot contribute through this operator. Unlike the released endpoint kernel, its query output is not constrained to $\operatorname{col}(X_Q)$ and can change partner contrasts. This representational capacity is necessary for the added task but does not establish accuracy on a new source.

### 5.2 Nested selection and computational controls

The matched panel contains fourteen singles and 55 combinations. Five outer folds contain 12, 11, 11, 11 and 10 query pairs. Every inner fit reconstructs its baseline from the permitted observations. Three inner folds select $\lambda\in\{.1,1,10\}$ by complete-rectangle error where available and complete-response error otherwise; ties favour stronger regularization. The unweighted correction defines M1, and training-measurement shrinkage defines R. The final readout mask contains 4,995 non-target genes. Neural models retain their 5,020-gene fitting axis and use this common scoring mask.

GEARS and Scouter are each fitted at five folds and three optimization seeds. GEARS compares learning-rate halving intervals one and five with caps twenty and sixty epochs, patience ten and batch size 32. Scouter uses archived GenePT v2 embeddings, batch size 256, loss coefficient .05 and forty epochs, with the predefined possible extension to eighty. Final fits restart for the selected duration. The neural validation uses the three rectangles in inner fold zero, whereas the kernel uses three inner folds. This is a matched test evaluation with different training and tuning budgets, not a compute-matched algorithm comparison.

Table S10 includes the complete-response instability of two Scouter seeds rather than replacing them by the best seed. Primary neural values are the means of seed-specific losses, with seed-level values retained. The best-seed contrasts are descriptive sensitivity comparisons. This averaging convention differs from the GSE220974 prediction ensemble in Note 9.

### Supplementary Table S10 | Complete matched RNA comparison and seed-specific results

All methods are scored on the same held-out combinations, 11 same-fit rectangles and 4,995 readouts. Neural values are means of the three seed-specific scores; fixed predictions have one realization. MSE columns are multiplied by 1,000 for display. Positive gain means lower rectangle error than RNA M0.

| Model | Rectangle MSE × 1,000 | Gain vs M0 % | Response MSE × 1,000 | Projected ratio |
| --- | --- | --- | --- | --- |
| Measured-single additive | 4.082949 | 3.5945 | 3.717907 | 1.000000 |
| Saturation-only | 4.344639 | -2.5844 | 3.254583 | 1.038787 |
| RNA M0 | 4.235183 | 0.0000 | 1.281651 | 1.018972 |
| RNA M1 | 3.835457 | 9.4382 | 1.214510 | 0.918185 |
| RNA R | 3.559154 | 15.9622 | 1.162690 | 0.865165 |
| GEARS | 3.990340 | 5.7812 | 5.154565 | 0.976397 |
| scouter | 4.085221 | 3.5409 | 18.063703 | 1.000400 |

M1/R full-response gains against M0 are 5.2387%/9.2819%. Notes 4 and 5 distinguish the response, rectangle and projected outcomes and describe the neural selection scope.

**Pairwise comparisons and individual optimization seeds**

Neural references are the mean score or the lowest-error seed; best-seed contrasts are descriptive.

| Reference | R gain vs mean/fixed % | R gain vs best seed % |
| --- | --- | --- |
| RNA M0 | 15.9622 | 15.9622 |
| RNA M1 | 7.2039 | 7.2039 |
| Measured-single additive | 12.8289 | 12.8289 |
| GEARS | 10.8058 | 8.3950 |
| scouter | 12.8773 | 12.5895 |

**Individual optimization seeds.** The same 55 biological combinations are reused; seed variation is not biological replication.

| Model / seed | Rectangle MSE × 1,000 | Response MSE × 1,000 | Projected ratio |
| --- | --- | --- | --- |
| GEARS / 0 | 4.027845 | 6.843357 | 0.966473 |
| GEARS / 1 | 3.885328 | 4.261756 | 0.969947 |
| GEARS / 2 | 4.057849 | 4.358583 | 0.992771 |
| scouter / 0 | 4.071770 | 1.717753 | 0.997512 |
| scouter / 1 | 4.090416 | 27.309634 | 1.001689 |
| scouter / 2 | 4.093476 | 25.163722 | 1.001998 |

**Selected final-fit epochs (folds 0–4).**

| Seed | GEARS | scouter |
| --- | --- | --- |
| 0 | 8, 7, 5, 7, 16 | 12, 4, 17, 5, 37 |
| 1 | 14, 14, 5, 17, 9 | 1, 5, 6, 3, 12 |
| 2 | 14, 9, 3, 5, 23 | 1, 5, 7, 11, 5 |

All GEARS fits selected schedule step 5. Scouter upstream commit: `fe2a73bdda02c9ec639bd55dabe758da2f7ff335`; see Note 5.2 for selection scope.

### Supplementary Table S11 | Matched RNA fold composition and fixed-model gains

| Fold | Pairs | Rectangles | Projected dimensions | M1 rectangle gain % | R rectangle gain % | R full gain % |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | 12 | 3 | 5 | 10.2795 | 16.0860 | 12.4161 |
| 1 | 11 | 2 | 4 | 16.8413 | 23.1330 | 10.8403 |
| 2 | 11 | 2 | 2 | -0.8170 | 0.0448 | 6.9471 |
| 3 | 11 | 2 | 2 | 10.1710 | 14.8266 | 10.0805 |
| 4 | 10 | 2 | 2 | 2.1388 | 15.8084 | 4.8530 |

The rows in Table S11 are derived from the seed-zero records of fixed models; copies of those same deterministic predictions stored alongside the other neural seeds are not counted again. Fold gains use their own baseline SSE. The pooled primary gain uses the sum of rectangle errors, rather than the mean of these fold gains.

### Supplementary Table S12 | All eleven same-fit RNA rectangles

| Rectangle | Fold (zero-based) | Four pairs in signed order |
| --- | --- | --- |
| 1 | 0 | FLI1+GSE1; FLI1+INTS1; GSE1+SMARCD1; INTS1+SMARCD1 |
| 2 | 1 | EP300+FLI1; EP300+INTS1; FLI1+GFI1; GFI1+INTS1 |
| 3 | 2 | DOT1L+EP300; DOT1L+MED24; EP300+HDAC3; HDAC3+MED24 |
| 4 | 3 | FLI1+HDAC3; FLI1+MED24; GFI1+HDAC3; GFI1+MED24 |
| 5 | 4 | GFI1+RING1; GFI1+SMARCD1; RCOR1+RING1; RCOR1+SMARCD1 |
| 6 | 0 | GSE1+HDAC3; GSE1+MED24; HDAC3+SMARCD1; MED24+SMARCD1 |
| 7 | 1 | DOT1L+RCOR1; DOT1L+RING1; FLI1+RCOR1; FLI1+RING1 |
| 8 | 2 | DOT1L+HDAC3; DOT1L+INTS1; HDAC3+RING1; INTS1+RING1 |
| 9 | 3 | EP300+GSE1; EP300+RCOR1; GFI1+GSE1; GFI1+RCOR1 |
| 10 | 4 | DOT1L+GFI1; DOT1L+SMARCD1; EP300+GFI1; EP300+SMARCD1 |
| 11 | 0 | GSE1+INTS1; GSE1+RING1; INTS1+MED24; MED24+RING1 |

The source rectangle-definition file uses one-based fold labels. Table S12 converts them to the zero-based fold labels used in the fitting records and Table S11. Its four listed pairs have signs $+,-,-,+$ in that order. The eleven rectangles span a task distinct from the full fifteen-dimensional endpoint-eliminated query space; overlap prevents treating them as eleven independent experiments.

## Supplementary Note 6 | Readout shrinkage, concentration and sensitivity

### 6.1 Derivation and interpretation of training-half shrinkage

For each outer training set, disjoint cell halves give $Z_1=P_TY_T^{(1)}$ and $Z_2=P_TY_T^{(2)}$. For readout $j$, define

$$w_j=\operatorname{clip}_{[0,1]}\frac{Z_{1,\cdot j}^\top Z_{2,\cdot j}}{(\|Z_{1,\cdot j}\|^2+\|Z_{2,\cdot j}\|^2)/2}.$$

Zero denominators yield zero weights. With $m=(Z_1+Z_2)/2$ and $d=(Z_1-Z_2)/2$, the numerator is $\|m\|^2-\|d\|^2$ and the denominator is $\|m\|^2+\|d\|^2$, column by column. The two formulas are an exact identity. Presenting them as two different confirmations would double-count the same construction.

Under the illustrative model $Z_1=s+\varepsilon_1$, $Z_2=s+\varepsilon_2$ with independent centred errors, the expected cross-product contains signal energy, whereas the expected denominator also contains noise energy. This motivates shrinkage but does not make the ratio unbiased or calibrated. Shared guide, capture or biological effects can remain in both halves. Accordingly, $w_j$ is a training-derived coefficient, not the probability that a gene interacts or a certificate of cross-study reproducibility. Prediction is changed to $M_0+G\operatorname{diag}(w)$ before scoring on the common unweighted axis.

### 6.2 Scalar, permutation and random-increment controls

A uniform-amplitude control fits $c\in[0,1]$ using fully refitted inner correction/error products. Its outer-fold coefficients are 0.838645, 1, 1, 1 and 1. It tests whether reducing the size of the whole correction accounts for the benefit. The readout-weighted predictor improves on this control in all eleven rectangles. A second control permutes the fixed weights within sixteen expression/variance strata, retaining the weight distribution but changing its alignment to readouts. Two hundred conditional permutations are retained, including the best one.

Applying the same weights to $M_0$ minus measured additivity tests whether the weights alone improve a different increment. Same-norm random increments test whether any correction of comparable energy could explain the result. None of the 200 weighted random increments improves on M0. These comparisons condition on the measured panel and saved predictions; their counts are not biological-replication probabilities.

### 6.3 Concentration, scale and fixed-readout deletion

Let $\Delta_j$ be baseline rectangle SSE minus candidate rectangle SSE for readout $j$, summed across the eleven rectangles. The positive-gain share of a fixed set $J$ is $\sum_{j\in J}\max(\Delta_j,0)/\sum_j\max(\Delta_j,0)$. Negative cancellation is $\sum_j\max(-\Delta_j,0)/\sum_j\max(\Delta_j,0)$. These denominators distinguish the share of positive contributions from the share of net improvement.

The retained ten-readout set is S100A9, S100A8, S100A12, LYZ, VCAN, VIM, MS4A6A, MYL9, LGALS1 and PPP1R27. It was identified during panel development and then fixed for the diagnostic deletion. It is not an independently specified biological program. Removing these readouts leaves R's gain at −0.3798%. The raw gain is thus concentrated rather than a broad improvement over all transcripts.

Across-condition standardization divides errors by readout variation estimated from outer-training combinations, with the declared floor based on the first percentile of positive training variances. Within-condition cell-variance standardization measures a different scale. R gives −0.1180% under the former and +2.0407% under the latter. Changing the variance floor to the recorded 0.1st and fifth percentiles retains the respective signs. Raw errors remain the predefined primary endpoint; the standardized analyses show what that endpoint emphasizes.

### Supplementary Table S13 | Complete shrinkage controls and scale sensitivity

Positive percentages indicate lower error than M0. All models, folds and readouts are fixed. Randomized rows summarize conditional perturbations of the retained predictions, not additional biological experiments.

| Predictor/control | Rectangle-error reduction % | Complete-response reduction % |
| --- | ---: | ---: |
| M1 | 9.4382 | 5.2387 |
| Training-inner uniform amplitude | 9.6988 | 5.3343 |
| Reliability-weighted R | 15.9622 | 9.2819 |
| R relative to uniform amplitude | 6.9362 | 4.1700 |
| Weight permutation, median of 200 | 3.2895 | See retained draw table |
| Weight permutation, best of 200 | 8.0423 | See retained draw table |
| Same weights on M0 minus measured additivity | 0.2368 | −134.596 |
| Same-norm random increment, weighted, median of 200 | −0.704 | See retained draw table |

R improved on the uniform-amplitude control in all 11 rectangles. None of the 200 permuted-weight predictions matched R; none of the 200 weighted random increments beat M0. These counts are descriptive controls. The scalar coefficients were 0.838645, 1, 1, 1 and 1 for folds 0–4.

**Across-condition variance standardization and fixed-ten-readout deletion.**

All entries below are rectangle-error reductions relative to M0 under the same transformation. Negative values mean deterioration. Standardization uses only outer-training combination variances; deletion uses the previously fixed ten-gene list in Note 6.3.

| Model | Raw gain % | Standardized gain % | Ten-gene deletion gain % |
| --- | --- | --- | --- |
| Measured-single additive | 3.5945 | 0.0244 | -0.3639 |
| Saturation-only | -2.5844 | -0.0999 | -0.6658 |
| RNA M0 | 0.0000 | 0.0000 | 0.0000 |
| RNA M1 | 9.4382 | -8.6708 | -10.0314 |
| RNA R | 15.9622 | -0.1180 | -0.3798 |
| GEARS | 5.7812 | -21.2418 | -9.9918 |
| scouter | 3.5409 | -7.7030 | -0.4520 |

**Concentration of positive gains.** Denominators are defined in Note 6.3; the top-ten share is not a share of net improvement.

| Model | Improved readouts / 4,995 | Ten-gene share of positive gains % | Negative cancellation % |
| --- | --- | --- | --- |
| RNA M1 | 1596 | 82.2127 | 54.6826 |
| RNA R | 1929 | 92.1253 | 9.5238 |

R also gave +2.0407% under within-condition cell-variance standardization. These variances measure a different scale from across-condition response variance. Cross-readout rectangle Spearman/direction agreement were −0.003834/49.35% (M0), 0.027902/50.91% (M1) and 0.026450/50.19% (R).

### 6.4 Interpretation of named readouts

The S100A9 and LYZ illustrations were selected after aggregate scoring to explain particular large changes. Their displayed beneficial coordinates rank first and tenth among 54,945 rectangle–readout coordinates by the archived absolute-benefit criterion. Their examples are not independent validation. The adverse LYZ rectangle involving GFI1/RCOR1 and RING1/SMARCD1 is retained, with a shrinkage weight of 0.879965. Some favourable examples are predicted more closely by M1 than by R, so a plot of those coordinates alone cannot explain the aggregate benefit of readout weighting. The bundle retains all eleven contrasts for the selected readouts and the four-corner values rather than only favourable points.

The modest cross-readout Spearman and direction-agreement results in Table S13 further delimit interpretation. A lower raw squared error supports improved prediction of the measured response scale. It does not, by itself, establish a transcriptome-wide molecular interaction map or causal direction for individual gene pairs.

## Supplementary Note 7 | Functional growth prediction and procedure transfer

### 7.1 Response construction and protected evaluation sets

The growth measurements derive from the Triple Sequencing screen [3]. Guide orientations are averaged within physical guide pair and source repeat, then equally across guide pairs. Processing requires finite measurements and both repeats; conflicting duplicate orientations stop processing. Published whole-map interaction scores are excluded from the fitting target because their fitted baseline could incorporate held-out combinations. The observed target is the source growth phenotype in its original units.

Identity hashing selects eighty K562 genes and their 3,160 pairs. The original registered 632-pair test set is excluded from revised development. A different 632-pair pool supplies three development folds, with 211, 211 and 210 queries. Each fold's high training pool contains 2,317, 2,317 or 2,318 pairs; its nested low pool contains 421, 421 or 422. The same test pairs compare training budgets within a fold. This design separates improvement in the baseline with more data from improvement contributed by the partner correction.

Jurkat retains 79 genes because GINS4 lacks the required single responses. There are 2,464 training and 617 evaluation pairs. The procedure was frozen before Jurkat fitting and scoring, but its coefficients were fitted from Jurkat training measurements. It therefore tests transfer of a procedure to another cellular context within the source study, rather than direct application of K562 weights or a wholly independent experimental platform.

### 7.2 Baseline and penalized shared correction

Each repeat standardizes its single-response effects by their across-gene standard deviation. In those units, the baseline contains symmetric endpoint, slope and curvature terms as specified in the main Methods. Endpoint deviations have zero-sum constraints and ridge penalty $10/s$, with $s$ the supported-gene count; global coefficients are unpenalized. Every inner training subset refits these terms. Predictions are returned to original growth units before model selection and evaluation.

The two repeats have separate baselines. Their mean training residual $\bar r_{gh}$ supplies a shared symmetric correction with objective

$$J(G)=\frac1{|T|}\sum_{(g,h)\in T}(\bar r_{gh}-G_{gh})^2+\frac{\tau}{n}\|G\|_*+\frac{\nu}{n^2}\|G\|_F^2.$$

For the primary signed rank-one family $G=\lambda uu^\top$ with $\|u\|_2=1$, the penalties reduce to $\tau|\lambda|/n+\nu\lambda^2/n^2$. The sign of $\lambda$ is unrestricted. This restriction is nonconvex; six starts and the recorded stationarity criterion support numerical fitting, not a proof of its global optimum. The spectral family instead optimizes over general symmetric matrices. Its nuclear-norm proximal step soft-thresholds signed eigenvalues, with the quadratic penalty included in the declared optimization [19].

The zero-gradient threshold $\tau_0=n\|\nabla\ell_T(0)\|_{\mathrm{op}}$ follows from the spectral-norm dual of the nuclear norm: zero satisfies the first-order condition when $\|\nabla\ell_T(0)\|_{\mathrm{op}}\le\tau/n$. The grid contains $\tau/\tau_0\in\{.3,.1,.03,.01\}$ and $\nu\in\{.1,1\}$, together with the explicit zero correction. Three inner folds select by pooled original-unit SSE across repeats, favouring stronger penalties for numerical ties. Unsupported endpoints receive zero correction.

A shared correction is identical in the two repeats by construction. Its equality is therefore not evidence of independently recovered interactions. The relevant check is whether the same correction improves prediction against each repeat's actual observed responses. The replicate-specific gains and adverse cases are retained below.

### Supplementary Table S14 | K562 budget comparisons, all folds and model selections

All gains are percentage reductions relative to M0 fitted on the same allowed training pairs. K562 combines SSE from three outer development folds; Jurkat uses one held-out evaluation. Each pair has two source replicate measurements.

| Environment / budget | Training pairs | Evaluation pairs | M0 pooled MSE | Rank-one / spectral gain % |
| --- | --- | --- | --- | --- |
| K562 low | 421 / 421 / 422 | 632 total | 0.00425513 | 0.2652 / -0.2879 |
| K562 high | 2,317 / 2,317 / 2,318 | 632 total | 0.00382441 | 11.0292 / 11.3865 |
| Jurkat high | 2,464 | 617 | 0.00092919 | 6.4308 / 8.1606 |

**Endpoint-eliminated error and replicate direction.** Projected gains use projected M0 error as denominator and sum within-fold errors.

| Environment / budget | Contrast dimension | Rank-one gain % | Spectral gain % | Positive raw-error replicate–fold counts |
| --- | --- | --- | --- | --- |
| K562 low | 132 + 131 + 130 | -0.0326 | -0.8591 | 2/6; 3/6 |
| K562 high | 132 + 131 + 130 | 11.6838 | 10.6710 | 6/6; 6/6 |
| Jurkat high | 538 | 6.1750 | 5.8975 | 2/2; 2/2 |

Positive counts are rank one; spectral. High-budget M0 improved 10.1225% over low-budget M0; the tabulated partner gains use their own same-budget baseline.

**K562 fold-specific gains and selected settings**

Fold labels 0–2 match the stored source identifiers. Parentheses give replicate 1 and replicate 2 gains; all entries are percentages versus the same-budget M0. Evaluation counts are 211, 211 and 210.

| Fold | Budget | Rank-one gain (replicates) % | Spectral gain (replicates) % | Projected gains: rank one / spectral % |
| --- | --- | --- | --- | --- |
| 0 | low | 0.0000 (0.0000, 0.0000) | -1.1441 (-0.9392, -1.5961) | 0.0000 / -2.6890 |
| 0 | high | 5.7110 (8.2111, 0.6430) | 5.4163 (6.2223, 3.7825) | 6.5305 / 2.9572 |
| 1 | low | 0.0000 (0.0000, 0.0000) | 0.0159 (0.1854, -0.2272) | 0.0000 / -0.3035 |
| 1 | high | 16.8664 (12.1849, 23.6444) | 14.1850 (11.1525, 18.5756) | 18.6244 / 15.8522 |
| 2 | low | 0.7412 (0.0487, 1.6475) | 0.1415 (0.1715, 0.1021) | -0.0948 / 0.3397 |
| 2 | high | 10.0215 (7.4512, 13.3659) | 13.5532 (10.3977, 17.6591) | 9.8256 / 13.0592 |

**Training-only choices.** Entries are $(\tau/\tau_0,\nu)$; zero means the selected predictor retained M0. The best penalized fit was saved separately when zero was chosen.

| Fold | Low: rank one | Low: spectral | High: rank one | High: spectral |
| --- | --- | --- | --- | --- |
| 0 | Zero | (0.01, 0.1) | (0.03, 0.1) | (0.3, 0.1) |
| 1 | Zero | (0.03, 1) | (0.03, 0.1) | (0.3, 0.1) |
| 2 | (0.3, 1) | (0.01, 1) | (0.03, 0.1) | (0.3, 0.1) |

All 288 inner fits, 12 outer correction fits and 144 objective-order checks met the recorded criteria.

### 7.3 Jurkat results and conditional sensitivity

Rank one is the prespecified primary functional extension and spectral correction is secondary. Both improve pooled error and each source repeat. Projecting each error through the fixed query graph gives 538 contrast dimensions and shows that the benefit is not entirely an endpoint-space change. These dimensions do not count independent biological observations.

The gene-reweighting calculation fixes the predictions and observed outcomes and varies endpoint weights over 10,000 paired draws. Its intervals describe dependence on the observed gene composition; they omit model refitting and do not establish coverage over newly sampled cell lines. The 79 leave-one-gene-out evaluations are overlapping deletions from the same test population. They show that the gain persists when any one observed endpoint's pairs are removed, rather than supplying 79 independent confirmations.

The M0 replicate-mean comparator is scored against both observed repeats. Its deterioration in complete error shows that averaging baseline predictions alone does not explain the partner model's improvement. The same comparator has a different projected outcome, reinforcing the need to retain complete and endpoint-eliminated losses separately.

### Supplementary Table S15 | Complete Jurkat primary, secondary and sensitivity results

The panel contains 79 genes, 2,464 training pairs and 617 test pairs. Rank one is primary and spectral correction is secondary.

| Quantity | Rank one: primary | Spectral: secondary | M0 replicate-mean reference |
| --- | --- | --- | --- |
| Pooled MSE | 0.00086943 | 0.00085336 | 0.00112689 |
| Pooled gain % | 6.4308 | 8.1606 | -21.2774 |
| Endpoint-eliminated gain % | 6.1750 | 5.8975 | 0.9177 |
| Replicate 1 gain % | 4.4542 | 8.0730 | -19.4882 |
| Replicate 2 gain % | 8.8845 | 8.2692 | -23.4982 |
| Conditional reweighting interval % | [1.7335, 12.1290] | [3.0017, 13.1873] | [-40.5760, -8.9395] |
| Minimum gain after one-gene deletion % | 4.9492 | 7.4284 | -22.9565 |
| Positive one-gene deletions | 79/79 | 79/79 | 0/79 |

Both partner families selected $(\tau/\tau_0,\nu)=(0.3,0.1)$. All 48 inner fits, 24 objective-order probes and two final correction fits passed. The 10,000 paired gene-reweighting draws condition on the observed measurements and fitted predictions; their intervals exclude refitting and are sensitivity summaries rather than calibrated coverage guarantees. The 79 overlapping deletions are likewise descriptive.

Named combinations illustrate the aggregate evaluation. For example, TTI2–PGK1 has a mean observed response near −0.321, baseline prediction −0.274 and corrected prediction −0.287, whereas TTI2–OGFOD1 has values near −0.036, −0.069 and −0.058. These are descriptive selected examples. The source bundle retains repeat-specific predictions and all evaluated pair results; these examples are not validation of RNA-based recommendations in Jurkat.

## Supplementary Note 8 | Model-based risk and empirical selection

### 8.1 Derivation of the linear minimum-risk coefficient

Consider $y=X\beta+u+\varepsilon$ and target $y_t=t^\top\beta+u_t$. Put $\Sigma=\operatorname{Cov}(u+\varepsilon)$, $k_t=\operatorname{Cov}(u+\varepsilon,u_t)$ and $k_{tt}=\operatorname{Var}(u_t)$. For a predictor $a^\top y$, unbiasedness for every $\beta$ requires $X^\top a=t$. Under the stated second moments, its risk is

$$R(a)=a^\top\Sigma a-2a^\top k_t+k_{tt}.$$

Assume $\Sigma$ is positive definite and the retained design columns have full rank. Completing the square gives

$$R(a)=k_{tt}-k_t^\top\Sigma^{-1}k_t+(a-\Sigma^{-1}k_t)^\top\Sigma(a-\Sigma^{-1}k_t).$$

Set $b=a-\Sigma^{-1}k_t$ and $\widetilde t=t-X^\top\Sigma^{-1}k_t$. Minimizing $b^\top\Sigma b$ subject to $X^\top b=\widetilde t$ by Lagrange multipliers gives

$$b=\Sigma^{-1}X(X^\top\Sigma^{-1}X)^{-1}\widetilde t.$$

Substitution yields the coefficient and risk $W_t$ in the main Methods [14,15]. A rank-deficient design can first be expressed in a full-rank basis, provided the requested endpoint target is estimable. Substituting a pseudoinverse without checking estimability does not resolve an unsupported target.

### 8.2 What is assumed by a train–query covariance

The training covariance $\Sigma$ does not determine $k_t$ and $k_{tt}$. For example, the same positive definite $\Sigma$ can be extended by many vectors $k$ and variances $v\ge k^\top\Sigma^{-1}k$ to form a positive semidefinite joint covariance. Those extensions have the same training law but different target risks and optimal coefficients. An observed training residual matrix alone therefore cannot identify arbitrary train–query covariance.

The analytic implementation supplies this continuation through standardized single-response programs, a training-fitted principal-component representation, symmetric pair descriptors and an RBF discrepancy component. These are explicit modeling assumptions. Estimating their parameters gives a plug-in risk score. The algebra is exact for the stated moments, whereas its performance with fitted moments is evaluated empirically. It supplies neither distribution-free coverage nor the synthetic causal model's confidence guarantee to the biological predictor.

### 8.3 Pair-balanced calibration and evaluation

Calibration uses five development splits, 130 records and eighty canonical pairs. All those pair identities are removed from evaluation, leaving 213 records from 33 pairs across thirty overlapping split seeds. The fixed risk threshold is the calibration pair-balanced median, 0.0923227147328771. A record of pair $q$ has weight proportional to $1/n_q$, so each canonical pair contributes equal total weight before selection.

Coverage is accepted record weight divided by total weight. Accepted loss uses accepted weights normalized over the accepted records; all-query loss uses weights normalized over all records. The loss ratio compares these two population summaries at the recorded coverage. Uncertainty resamples canonical pairs, carrying all their records and the fixed acceptance labels. Individual records from repeated fits are not treated as independent samples.

The subsequent weighted-error evaluation has 207 complete records from 32 pairs. Six missing weighted endpoints belong to one pair and are omitted without imputation. The original risk threshold and decisions are unchanged for retained records. Its result is a follow-up endpoint on the same prediction process, rather than a second independent selective-prediction experiment.

### 8.4 Empirical ranking of the released predictor

A separate analysis learns relative error rankings from training-only out-of-fold losses of the released predictor. Features include predicted effect and an eight-member condition-subsampling ensemble variance. The categorical-exposure comparator C and support model S share weighted standardization and ridge fitting; S additionally uses inverse training-combination degree. Nominal coverage comparisons use common query populations. Within-exposure reselection fixes quotas before comparing which queries each score accepts.

Table S16 retains the weak Norman comparison as well as the favourable Wessels comparison. Additional effect-quartile restrictions have intervals crossing zero and remain in the source tables. These outcomes test whether support features add ranking information beyond the declared controls. They do not establish that adding a support score is beneficial on every dataset. This ranking, the analytic $W_t$ score and RNA readout shrinkage are three different constructions.

### Supplementary Table S16 | Analytic-risk evaluation and empirical support ranking

**Separate analytic-risk threshold experiment**

| Outcome | Original MSE endpoint | Subsequent WMSE complete case |
| --- | --- | --- |
| Queries / canonical pairs | 213 / 33 | 207 / 32 |
| Accepted records | 110 | 106 |
| Pair-balanced coverage | 49.48% | 48.95% |
| Accepted/all loss ratio | 0.5343 | 0.4775 |
| Pair-cluster percentile interval | [0.3646, 0.7591] | [0.3175, 0.7071] |

The six weighted-unevaluable queries are from one pair; no replacement weights or outcomes were introduced. Pair-balanced coverage and losses use the same record-weight definition in Note 8.3.

**Released-predictor empirical ranking at 50% nominal coverage**

| Dataset / pool | Support vs categorical-exposure MSE improvement % | 95% conditional pair-bootstrap interval |
| --- | ---: | --- |
| Norman, all queries | −2.09 | [−11.23, 10.84] |
| Norman, seen1 reselection | 1.83 | [−18.42, 16.10] |
| Wessels, all queries | 8.27 | [1.75, 15.30] |
| Wessels, seen1 reselection | 17.57 | [5.52, 26.11] |

Exposure-aware controls narrow the interpretation of earlier gains against effect/variance-only ranking. Coarse effect-adjusted comparisons did not establish a general effect-independent advantage. Full retained coverage/loss grids are in machine-readable sources.

## Supplementary Note 9 | PairPick: response constraints, candidate pools and completed comparisons

### 9.1 Cohort construction and the intended decision

PairPick takes an anchor intervention, a fixed transcriptional objective and a prescribed candidate partner pool. It returns response predictions, predicted utilities and the highest-utility candidate. Held-out measurements assess that recommendation. The GSE220974 evaluation uses the source's mixed activation/repression measurements [23], keeping mode identities distinct throughout preprocessing and fitting.

The operational cohort contains all 49 detected double labels involving two distinct regular genes, representing 28 physical pairs. Enhancer and same-gene doubles are excluded. Control identities are reconstructed from guides, removing 860 enhancer-target cells previously assigned to the control category. The source's intended designed-construct count is smaller than the detected-label cohort. The original-table audit obtained the original intended-construct table: 15 labels match gene-mode designs directly, one has only an enhancer-level correspondence, and 33 are not listed as designed pairs. Cell-level guide assignment is not verified by a gene-mode match; incidental co-delivery, assignment error or doublets remain possible (Supplementary Note 13). The cohort was not filtered on response strength to force agreement with the source design count.

Counts are transformed by $\log(1+10^4c/\ell)$, centred against capture-matched non-targeting controls and shifted to the equal D1/D2/D3 control reference. Negative corrected values are retained. Variance selection from controls and singles chooses 5,000 genes; adding physical targets yields 5,006. Evaluation excludes all nineteen physical targets, leaving 4,987 readouts. The intervention axis has 35 mode-specific identities.

The objective $v_a$ uses only the anchor's single-response measurements. Complementary halves are formed within condition/capture strata with seed 20261102 and centred against their own controls. The full single effect is multiplied by clipped half-agreement, target genes are removed and the vector is normalized to unit L1 norm. Utility is $U(a,b)=v_a^\top y_{ab}$. The same objective and candidate pool apply to all predictors. The task rewards reproducing the anchor-associated transcriptional direction; it does not directly measure survival, drug synergy or therapeutic benefit.

### 9.2 Exact one-dimensional selection constraint in R-G

For each $\lambda\in\{.1,1,10\}$ and power $p\in\{1,2\}$, write the candidate increment as $A=G_\lambda\operatorname{diag}(w^p)$. Inner-fold predictions use $M_0+\beta A$, with $\beta\ge0$. Accumulating rectangle and complete-response residual products gives

$$L_P(\beta)=p_0+2p_1\beta+p_2\beta^2,\qquad L_R(\beta)=r_0+2r_1\beta+r_2\beta^2,$$

where $p_1$ and $r_1$ are error–increment inner products and $p_2,r_2$ are squared increment norms. The admissible set is

$$\mathcal B=\{\beta\ge0:L_R(\beta)\le L_R^{\mathrm{ref}}\},$$

using a separately nested original-R reference. When $r_2>0$, nonempty quadratic sublevel sets lie between the roots

$$\beta_\pm=\frac{-r_1\pm\sqrt{r_1^2-r_2(r_0-L_R^{\mathrm{ref}})}}{r_2},$$

intersected with $[0,\infty)$. If the discriminant is negative beyond numerical tolerance, the candidate is infeasible. Zero-curvature cases use the corresponding linear or constant inequality. Because these coefficients arise from squared norms, exact zero $r_2$ also implies exact zero $r_1$; the implementation still handles numerical degeneracies.

For $p_2>0$, the optimum is the projection of $-p_1/p_2$ onto $\mathcal B$. In a flat objective the archived convention prefers an admissible amplitude nearest one. There is no upper amplitude cap. Numerical comparisons use the recorded relative tolerance based on $64\epsilon$; they do not introduce a substantive relaxation chosen from outer outcomes. Candidate ordering prefers the retained original-R fallback, then power one, stronger ridge and an amplitude nearer one under ties.

Thus the guard constrains complete RNA loss on the nested selection observations while optimizing partner loss there. It does not guarantee outer-test improvement in either quantity. R-G uses the native partner kernel; graph-derived ranking features belong to other archived heads. The name PairPick changes the user-facing application identity without replacing the algorithm or claiming a new frozen experiment.

### 9.3 Splits, contrast degrees and candidate-pool dependence

Five balanced, hash-ordered outer folds group all mode variants of a physical gene pair. Three inner folds use the same identity-only grouping rule. Outer query counts are 10, 10, 11, 9 and 9, with partner dimensions 0, 1, 2, 0 and 1. Zero-dimensional folds have no estimable partner-loss endpoint; they are not assigned zero partner error. Complete RNA scoring still includes their predictions.

Fixed-anchor pools are formed within an outer fold and require at least two candidates. Table S19 lists every pool. Twenty-eight tasks arise from eight physical anchors, and five pools have at least four candidates. Tasks sharing an anchor, a candidate or a measured pair are dependent. A target label that changes activation/repression mode is not automatically a new independent physical anchor.

For task $t$, regret is $\max_{b\in B_t}U_t(b)-U_t(\widehat b_t)$. Since the observed maximum is shared across methods, the paired difference in their regrets is the negative of their paired utility difference. Empirical-best hits use the archived deterministic tie rule. They describe the best observed candidate among the included pool; an unmeasured global optimum is not being evaluated.

### 9.4 Full comparisons and interpretation

Table S17 includes every archived result row, including individual neural seeds, historical ranking heads and the PRESAGE adapted reference. The primary comparison uses R-G, standard ridge selected by RNA error, single-response addition and the two neural prediction ensembles. The table does not replace either ensemble by a favourable seed. Standard ridge has the lowest complete RNA MSE in this primary set; GEARS has the lowest partner MSE; R-G has the lowest point estimate of regret and the most hits. These are distinct endpoints.

R-G and the PRESAGE combination adapter select the same candidate in all 28 tasks, giving equal utility and regret despite different response predictions. Scouter seed zero also has lower regret than R-G in the complete historical table, but it was not selected as the primary neural comparator. TargetRank's mean regret is 27.59% higher than R-G's on this source. These observations support retaining recommendation as a demonstrated application of the response system, without attributing universal ranking superiority to a head.

### Supplementary Table S17 | All fifteen GSE220974 result rows

| Archived model | Primary set | RNA MSE | Partner MSE | Mean regret | Hits / 28 | Larger-pool regret |
| --- | --- | --- | --- | --- | --- | --- |
| `Scouter_seed0` | — | 0.004885 | 0.003625 | 0.005685 | 18 | 0.009846 |
| `PRESAGE_Combo_seed24` | — | 0.005039 | 0.003608 | 0.006197 | 19 | 0.009846 |
| PairPick (R-G) / `R_G` | Yes | 0.004903 | 0.003614 | 0.006197 | 19 | 0.009846 |
| `Scouter_seed1` | — | 0.004666 | 0.003629 | 0.006497 | 18 | 0.009846 |
| `StandardRidge_RNA` | Yes | 0.004880 | 0.003631 | 0.006591 | 18 | 0.012052 |
| `Scouter_Ensemble` | Yes | 0.005199 | 0.003627 | 0.006623 | 17 | 0.009846 |
| `GEARS_seed0` | — | 0.006020 | 0.003398 | 0.006859 | 17 | 0.012949 |
| `PairPick_Stable` | — | — | — | 0.007686 | 19 | 0.009846 |
| `PairPick_TargetRank` | — | — | — | 0.007908 | 18 | 0.012052 |
| `GEARS_Ensemble` | Yes | 0.006046 | 0.003528 | 0.010258 | 18 | 0.034809 |
| `GEARS_seed1` | — | 0.006128 | 0.003625 | 0.010551 | 19 | 0.034809 |
| `SingleB_additive` | Yes | 0.006096 | 0.003631 | 0.012220 | 18 | 0.049277 |
| `StandardRidge_Choice` | — | 0.005037 | 0.003631 | 0.012292 | 16 | 0.043387 |
| `Scouter_seed2` | — | 0.010245 | 0.003628 | 0.014407 | 12 | 0.041181 |
| `GEARS_seed2` | — | 0.006991 | 0.003619 | 0.014968 | 16 | 0.052700 |

An em dash in a response-score column denotes an unavailable response-level score for a ranking-only output, not zero error. An em dash in the primary-set column identifies a retained secondary or historical row. Larger-pool regret is evaluated on the same five tasks with at least four candidates.

### Supplementary Table S18 | Paired utility differences for R-G and retained references

| Reference | R-G minus reference utility | 95% anchor-bootstrap interval | Anchor clusters |
| --- | --- | --- | --- |
| PRESAGE_Combo_seed24 | 0.000000000 | [0.000000000, 0.000000000] | 8 |
| StandardRidge_RNA | 0.000393849 | [0.000000000, 0.001140805] | 8 |
| StandardRidge_Choice | 0.006094557 | [0.000091979, 0.019584558] | 8 |
| SingleB_additive | 0.006022491 | [-0.002529016, 0.020191752] | 8 |
| GEARS_Ensemble | 0.004060874 | [-0.006994990, 0.022476783] | 8 |
| Scouter_Ensemble | 0.000425511 | [0.000000000, 0.000882707] | 8 |

Intervals in Table S18 resample eight physical anchors and are conditional on the observed pools and predictions. With few clusters and repeated ties, bootstrap endpoints at zero do not establish a general strict advantage. No candidate set or anchor was removed to improve the result. Complete task-level choices are retained in `gse220974/CHOICES.csv` (420 rows: fifteen models by 28 tasks).

### Supplementary Table S19 | All fixed-anchor candidate pools

| Task | Fold | Mode-specific anchor | Candidates | Partners |
| --- | --- | --- | --- | --- |
| 1 | 0 | CEBPB-a | 2 | FOSL1-i, SPI1-a |
| 2 | 0 | CEBPB-i | 2 | FOSL1-a, SPI1-a |
| 3 | 0 | FOSL1-a | 2 | CEBPB-i, LYL1-i |
| 4 | 0 | GATA1-i | 2 | JUN-a, SPI1-a |
| 5 | 0 | JUN-a | 2 | GATA1-i, SPI1-a |
| 6 | 0 | SPI1-a | 5 | CEBPB-a, CEBPB-i, GATA1-a, GATA1-i, JUN-a |
| 7 | 1 | FOSL1-a | 2 | KLF1-a, KLF1-i |
| 8 | 1 | GATA1-a | 2 | KLF1-a, KLF1-i |
| 9 | 1 | KLF1-a | 6 | CEBPB-i, FOSL1-a, FOSL1-i, GATA1-a, GATA1-i, JUN-i |
| 10 | 1 | KLF1-i | 4 | CEBPB-a, FOSL1-a, GATA1-a, JUN-a |
| 11 | 2 | FOSL1-a | 2 | GATA1-a, GATA1-i |
| 12 | 2 | GATA1-a | 3 | FOSL1-a, FOSL1-i, YY1-a |
| 13 | 2 | GATA1-i | 2 | FOSL1-a, YY1-a |
| 14 | 2 | KLF1-a | 2 | LYL1-a, LYL1-i |
| 15 | 2 | LYL1-a | 3 | KLF1-a, KLF1-i, SPI1-a |
| 16 | 2 | LYL1-i | 3 | CEBPB-a, KLF1-a, SPI1-a |
| 17 | 2 | SPI1-a | 2 | LYL1-a, LYL1-i |
| 18 | 2 | YY1-a | 2 | GATA1-a, GATA1-i |
| 19 | 3 | FOSL1-a | 2 | BTF3-i, YY1-a |
| 20 | 3 | GATA1-i | 2 | CEBPB-a, LYL1-i |
| 21 | 3 | LYL1-i | 2 | GATA1-a, GATA1-i |
| 22 | 3 | YY1-a | 5 | CEBPB-i, FOSL1-a, FOSL1-i, JUN-i, SPI1-a |
| 23 | 4 | FOSL1-a | 2 | JUN-i, SPI1-a |
| 24 | 4 | FOSL1-i | 2 | JUN-a, SPI1-a |
| 25 | 4 | KLF1-a | 2 | SPI1-a, YY1-a |
| 26 | 4 | KLF1-i | 2 | SPI1-a, YY1-a |
| 27 | 4 | SPI1-a | 4 | FOSL1-a, FOSL1-i, KLF1-a, KLF1-i |
| 28 | 4 | YY1-a | 3 | KLF1-a, KLF1-i, LYL1-i |

Suffixes `a` and `i` denote activation and repression. Pool membership is fixed by the archived split, not optimized for a method. The five larger pools are a subset of the 28 tasks and are not counted as additional observations.

### 9.5 Completed neural training and mode adapters

GEARS uses disjoint mode-labelled copies of its perturbation graph and a physical-gene coexpression graph inferred from permitted controls and singles. Scouter uses direct-sum GenePT embeddings, activation $[e,0]$ and repression $[0,e]$. These adapters expose mode identity to the upstream architectures; their results are not an unmodified mixed-mode capability claimed by the upstream papers [8,18]. Standard ridge uses an unpenalized intercept and separate mode-specific endpoint columns, with equal condition weights for controls, singles and training doubles.

GEARS uses Adam at .001, weight decay .0005 and batch size 32; its two learning-rate-halving schedules have intervals one/five and caps twenty/sixty, with patience ten. Scouter uses Adam at .001, batch size 256 and exponential learning-rate factor .9. Its forty-epoch selection run has a predefined possible extension to eighty when the recorded late-improvement criterion is met. Table S20 reports the actual selected durations; no GSE220974 Scouter run took that extension. Each final fit restarts and trains for the selected epoch count.

### Supplementary Table S20 | Actual GSE220974 final-fit durations for all thirty jobs

| Model | Seed | Final-fit epochs, folds 0–4 | LR halving interval, folds 0–4 |
| --- | --- | --- | --- |
| gears | 0 | 12, 5, 20, 30, 10 | 5, 5, 5, 5, 5 |
| gears | 1 | 26, 5, 5, 15, 5 | 5, 1, 1, 1, 5 |
| gears | 2 | 7, 7, 1, 7, 1 | 1, 1, 1, 1, 5 |
| scouter | 0 | 40, 40, 5, 37, 3 | Exponential factor 0.9 |
| scouter | 1 | 40, 39, 15, 26, 3 | Exponential factor 0.9 |
| scouter | 2 | 39, 40, 2, 34, 1 | Exponential factor 0.9 |

The completion record verifies thirty frozen, independently archived and replayed jobs: five folds by three seeds for each model. GEARS folds three and four for all three seeds were completed on GPU; the other final jobs used the recorded CPU environment. GPU execution used FP32 without AMP or TF32. The archive records the CPU PyTorch 2.14.1 and GPU PyTorch 2.1.2+cu118 environments, with the GPU runtime's Python 3.10.8 and NumPy 1.26.4. These differences are part of execution provenance; they do not establish identical stochastic trajectories across environments.

For each neural model, the primary prediction is the arithmetic mean of three saved prediction vectors, followed by RNA, partner and decision scoring. In general, scoring this mean differs from averaging seed-specific losses or recommendations. Some alternative models' outer scores were known before neural completion. Fixed neural budgets and selection rules were retained, but the whole comparison is reported as retrospective source evaluation rather than an untouched preregistered confirmation.

## Supplementary Note 10 | Independent-source checks and acquisition results

### 10.1 Frozen Dixit RNA comparison and subsequent development

The independent RNA check uses ten permitted TF singles, 45 pairs and 23,519 non-target readouts from Dixit [21]. Three outer folds each withhold fifteen pairs, with all baseline and correction selection restricted to the remaining thirty. Each query graph has five endpoint-eliminated dimensions. A predefined G2M program maps 52 of 54 unique symbols and supplies a secondary program endpoint.

The selected uniform correction improves complete RNA error by 0.26980% but increases the primary endpoint-eliminated error by 0.03260%. M0 itself substantially improves complete response relative to measured additivity; that baseline performance does not establish the added partner capability. Table S21A preserves the full frozen comparison, its fold outcomes and conditional sensitivity. The failed primary partner result remains the source-transfer conclusion.

Subsequent G2M development introduced a fixed-reference continuation. A projection defined from the current query batch can change predictions when an unrelated query is added or removed. The fixed full-roster construction instead locks the reference identities, giving a query-batching-invariant mapping. Its 14.2994% projected gain and 2.2671% complete gain are relative to the independently fitted scalar saturation baseline in that later program task. Two implementations with those identical values are mathematically equivalent and are counted once. This development result does not revise the earlier frozen all-RNA decision.

### 10.2 Frozen Combi-seq program evaluation

The Combi-seq evaluation [22] contains 374 drug pairs, 36 endpoints and 51 program genes. Three outer folds have 125, 125 and 124 query pairs; their full query graphs supply 266 contrast dimensions. The principal system adds the fixed-reference program correction to the shared endpoint baseline, keeping the same outer scaling and program units for all methods.

Its complete gain relative to that endpoint baseline is −0.023282%, and its projected gain is −0.030240%. Comparison with scalar saturation alone gives a small complete gain, but the endpoint baseline already accounts for that improvement. The proper added-component comparison therefore remains adverse. Three retained experimental-run targets are all evaluated against the same predictions fitted from pooled training means. Their mixed results, including two adverse complete-error changes, are retained in Table S21D. These targets are repeated measurements of the same pairs, not three independently selected success tests.

### Supplementary Table S21 | Complete independent-source and subsequent program results

**A. Frozen Dixit RNA and fold sensitivity**

All entries are percentage reductions relative to M0; positive values indicate lower error. A negative value for a comparator indicates that it has larger error than M0, not the reciprocal improvement of M0.

| Method | Full response | Projected primary | Full G2M | Projected G2M |
|---|---:|---:|---:|---:|
| M0 / zero | +0.00000% | +0.00000% | +0.00000% | +0.00000% |
| Unit | +0.22204% | -0.12901% | +3.02328% | +6.60838% |
| Uniform | +0.26980% | -0.03260% | +2.56969% | +4.03391% |
| Amplitude-only | +0.25987% | -0.03322% | +2.45074% | +3.60235% |
| Amplitude-noise | +0.07927% | -0.00355% | +0.56885% | +0.39142% |
| Selected | +0.26980% | -0.03260% | +2.56969% | +4.03391% |
| Measured additive | -30.60428% | +0.00473% | +19.14509% | +0.01230% |
| Saturation-only | -30.26808% | +0.00000% | +19.21281% | +0.00000% |

All 45 pairs were withheld once in the three fixed outer folds. The gene program was specified before this evaluation and uses 52 of 54 exact unique mappings. The primary partner result is negative; the secondary application does not outperform measured additivity in complete error.

| Fold | Train / test pairs | Selected c | Core gamma | Full gain | Projected gain |
|---|---|---:|---:|---:|---:|
| 0 | 30 / 15 | 0.50 | 0.77691 | +0.20136% | +0.02930% |
| 1 | 30 / 15 | 0.75 | 0.78138 | +0.41490% | -0.04779% |
| 2 | 30 / 15 | 0.75 | 0.76559 | +0.19528% | -0.07727% |

All folds selected uniform shrinkage, lambda=10, alpha=0.5 and eta=1. Each test graph has rank ten and five contrast dimensions. The complete core is reselected within every allowed nested fit, not selected on these three outer outcomes.

| Sensitivity | Complete gain | Projected gain | Interpretation |
|---|---:|---:|---|
| Positive-node multiplier, central 95% | +0.10351% to +0.44214% | −0.10916% to +0.07366% | Fixed predictions; projections refitted to weights |
| Delete one endpoint, minimum to maximum | +0.24285% to +0.35761% | −0.08717% to +0.05049% | Ten deletions; test projection rebuilt |

These ranges are conditional descriptions rather than biological-replicate confidence intervals. Original protocol, implementation correction and freeze receipts are retained in the source record.

**B. Later Dixit program development**

All rows evaluate the same 45 outer queries and the same program units. Positive gains mean lower error than the independently fitted scalar saturation baseline. The selected full-roster reference is the primary development result; all other rows retain their specified comparator identities.

| Construction | Total SSE | Projected SSE | Total gain % | Projected gain % | Query batching invariant |
|---|---:|---:|---:|---:|---|
| Scalar saturation | 6.538922596 | 2.307093784 | 0.0000 | 0.0000 | Yes |
| Unprotected, matched ridge | 6.695915797 | 1.977192992 | -2.4009 | 14.2994 | Yes |
| Unprotected, own inner selection | 6.550751732 | 1.961337847 | -0.1809 | 14.9866 | Yes |
| Current-query projection | 6.209021805 | 1.977192992 | 5.0452 | 14.2994 | No |
| Fixed full-roster reference | 6.390677159 | 1.977192992 | 2.2671 | 14.2994 | Yes |
| Fixed unmeasured-pool reference | 6.278276738 | 2.046447926 | 3.9861 | 11.2976 | Yes |
| Reference kernel, full residual loss | 6.710468769 | 2.234667612 | -2.6235 | 3.1393 | Yes |
| Reference kernel, training-projected loss | 6.390677159 | 1.977192992 | 2.2671 | 14.2994 | Yes |

The two fixed-reference implementations with 2.2671% total gain are mathematically equivalent and are counted once. Keeping projected gain unchanged at fixed coefficients follows from endpoint cancellation. Conditional sensitivity and adverse endpoint deletion remain in the corresponding fixed-reference source records.

**C. Frozen Combi-seq comparison**

All methods predict the same 374 held-out pairs. Gains in the table are percentages relative to the shared endpoint baseline B+E. Positive values mean lower error. Program units and outer-fold scaling are identical across methods; projected losses pool the 266 fixed test-graph contrast dimensions.

| Method | Complete SSE | Projected SSE | Complete gain % | Projected gain % |
|---|---:|---:|---:|---:|
| Measured-single additive | 177.237788737 | 125.915191495 | -3.006695 | +0.475131 |
| Scalar saturation B | 172.818152532 | 126.516309700 | -0.438100 | +0.000000 |
| Endpoint baseline B+E | 172.064339459 | 126.516309700 | +0.000000 | +0.000000 |
| Previous $B+g_F$ | 172.876243066 | 126.554568339 | -0.471860 | -0.030240 |
| Principal $B+E+g_F$ | 172.104399943 | 126.554568339 | -0.023282 | -0.030240 |

The principal system reduces complete error by 0.413008% relative to scalar saturation, while the endpoint baseline alone reduces it by 0.436189%. The preceding saturation-plus-program assembly gives -0.033614% versus saturation. The principal partner increment is negative, and no alternative comparator changes that decision.

**D. Combi-seq folds and experimental-run targets**

Fold outcomes compare the principal system $B+E+g_F$ with the same endpoint baseline $B+E$. All three folds selected alpha=0.5. Lambda concerns the unchanged partner operator; noise and gamma concern the shared endpoint baseline.

| Fold | Train / test | Lambda | Noise / gamma | Complete % | Projected % |
|---|---|---:|---|---:|---:|
| 0 | 249 / 125 | 10 | 0.3 / 0.997148 | +0.035808 | -0.000757 |
| 1 | 249 / 125 | 0.1 | 0.01 / 0.445541 | -1.798455 | -0.767833 |
| 2 | 250 / 124 | 0.1 | 0.01 / 0.327347 | +0.433845 | -0.500088 |

All three retained experimental-run targets were scored against the same predictions fitted to pooled training means. They were not used to select different models or independent success thresholds. Each row includes the same 374 pair identities.

| Run target | Endpoint SSE | Principal system SSE | Complete gain % |
|---|---:|---:|---:|
| rep1 | 82.344587094 | 82.541384556 | -0.238993 |
| rep2 | 1372.160073814 | 1372.332694951 | -0.012580 |
| rep3 | 57.254705045 | 57.005467898 | +0.435313 |

Only one of three run targets improves. The principal frozen decision is not confirmed. The conditional complete-gain range is -1.098310% to +0.688105%; this is not a biological-replicate confidence interval.

### 10.3 One-measurement acquisition: development and outer evaluation

The acquisition experiment asks which additional measured condition to acquire for a specified target, a different decision from PairPick's choice among candidate interventions based on predicted utility. The target-aligned development procedure improves partner error by 4.68324% and target-AB error by 22.52110% relative to the expected uniform acquisition reference. Its primary development run gate required at least 5% partner improvement and was not met. An explicit subsequent continuation retained that result and evaluated the outer targets rather than redefining the original gate as passed.

On the six outer targets the corresponding improvements are negative: −2.41963% for partner error and −46.91288% for target-AB error. Table S22 keeps both stages. No outer target is removed from the aggregate decision. A favourable geometric-width change or a useful development example is not substituted for measured outer predictive utility. The experiment therefore supplies a boundary on active acquisition, not an additional successful end-to-end application of PairPick.

### Supplementary Table S22 | Development and outer acquisition outcomes

| Stage | Targets | Partner gain % | Target-AB gain % | Complete gain % | Partner gains by repeat % |
| --- | --- | --- | --- | --- | --- |
| development | 6 | 4.68324 | 22.52110 | 13.68849 | 5.00086, 4.38026 |
| outer | 6 | -2.41963 | -46.91288 | -19.51025 | -10.81939, 2.08715 |

All gains use the expected uniform acquisition reference and retain both measured repeats.

### 10.4 Related development analyses and evidence accounting

Earlier fixed-G integration, measurement-reference calibration and shared-feature development are retained in the code archive and the SI source bundle. Measurement-reference calibration concerns better input estimation for M0; it is not evidence of a new partner head. The qualified shared-feature tasks reuse overlapping rectangles, and their case-conditioned comparator has adverse results. Their full tables are provided as `additional_development_tables.md` with their original task labels. These analyses are not pooled with the frozen Dixit or Combi-seq tests or counted as independent confirmations of the matched-panel partner claim.

The biological results establish different levels of transfer: a substantial matched-panel raw RNA gain, a procedure-frozen growth gain in another cell context, a retrospective recommendation demonstration on detected labels and unresolved or adverse source-transfer results for the specified RNA/program extensions. This differentiation defines the current empirical scope of the common response construction.

## Supplementary Note 11 | Statistical units, original outcomes and amendments

### 11.1 Unit of inference

The synthetic experiment uses paired random-noise seeds as the Monte Carlo unit, averaging the four generating worlds inside each seed. The benchmark's paired resampling groups dataset and perturbation identity and keeps repeated split records together. Historical combination-error resampling groups canonical pairs. RNA cell halves estimate training agreement; they are not independent biological replicates. Rectangle overlaps and shared endpoints are reported as dependence, and readouts are not substituted for independent screens.

Growth errors pool the two measured repeats and retain each repeat's direction. Gene-reweighting and endpoint-deletion analyses condition on the observed screen. GSE220974 selection intervals group by physical anchor, not the 28 nominal tasks. The reported intervals therefore have different meanings: Monte Carlo estimation uncertainty, conditional pair/anchor resampling uncertainty or fixed-prediction sensitivity. They are labelled at the point of use.

No universal significance claim is formed by counting all favourable secondary analyses. The primary endpoint and comparator of each experiment define its decision. Descriptive seed comparisons, selected molecular examples, scaling changes and alternative ranking heads remain supplementary evidence. In particular, an improvement in complete RNA error does not establish partner improvement, and neither error automatically establishes lower recommendation regret.

### 11.2 Original plans and local freezes

The original M1 model-development target was at least 10% rectangle improvement; its 9.4382% outcome did not meet that criterion. The later R result is development on the same panel. Original registered K562 experiments selected zero correction at all five budgets. Subsequent revised K562 results use a different evaluation pool and retain model-development status. The Jurkat transfer procedure was locally frozen before fitting and scoring after an availability-only amendment removing GINS4.

Table S23 preserves this chronology. Times are explicitly UTC. A local prediction lock identifies when predictions were fixed relative to scoring. A later external deposit documents those records but cannot retroactively preregister the revised method. Similarly, the GSE220974 completion records prove completed jobs and retained outputs, not that no related outer score had previously been seen.

### Supplementary Table S23 | Original outcomes, amendments and functional prediction locks

Times below are UTC on 18 October 2026; China Standard Time is UTC + 8 hours. Local model freezes are distinct from external registration and later results deposit.

| Study stage | Fixed comparison or event | Outcome and interpretation |
| --- | --- | --- |
| RNA M1 development | Original partner-kernel contrast criterion: at least 10% improvement | 9.4382%; original criterion not met |
| RNA R development and subsequent diagnostic freeze | Same 55-pair panel; readout reliability added, then retained for diagnostics | 15.9622%; model-development result, not independent RNA confirmation |
| Original registered K562 | Separate-repeat ranks 0–3, five budgets, original 632-pair test; 20% primary | Both primary gains zero; 266/270 positive-rank inner fits failed at primary budget |
| Original K562 scoring, 12:43:30–12:43:31 | Formal outcomes evaluated | Known before revised K562 development |
| Revised K562 development | New matrix penalties/shared correction; distinct 632-pair evaluation from non-test pool | Three-fold budget comparison; original formal-test pairs excluded |
| K562 completion plan, 15:09:48 | Finalized smoke procedure applied to remaining folds | Fold 0 previously examined; three-fold results remain development evidence |
| Initial Jurkat freeze, 15:31:12 | Fixed procedure and original 80 gene identities | Availability check stopped for missing GINS4 singles |
| Availability-only amendment, 15:32:20 | Remove GINS4 without replacement and refreeze | 79 genes; before fitting or combination outcome summaries |
| Jurkat prediction lock, 15:32:52 | Predictions and selection records saved and hashed | Trainer had not read sealed evaluation outcomes |
| Jurkat scoring, 15:33:24 | Evaluate 617 held-out pairs | Primary rank-one criterion met; spectral remains secondary |
| Subsequent OSF results deposit | Two-stage archive and disclosure documents | Preserves follow-up, without retroactively preregistering revised methods |

The original five K562 budgets used 158, 316, 632, 1,264 and 2,528 training pairs. All selected the baseline, with incremental gains zero and pooled MSEs 0.00596059, 0.00450399, 0.00396413, 0.00369680 and 0.00356008. Zero denotes identical selected predictions, not a biological effect estimated to be exactly zero. Original outcomes and computational records remain in the registered-plan status documents; subsequent development does not replace them.

The historical pre-amendment Jurkat protocol ZIP is supplied within the complete scientific archive and is additionally exposed at `reproduction/protocol/v1_before_availability_amendment.zip`. The recovered file has 24,655 bytes and passes ZIP CRC checks. It includes the original `FREEZE.json`, `PLAN.json`, 80-gene list and preparation, fitting and scoring code. `reproduction/RECOVERED_JURKAT_PROTOCOL.json` records the complete two-level nested member path, hashes and member inventory. The direct copy is byte-identical to the nested original. This protocol records the original local freeze and subsequent availability amendment; it was not externally preregistered.

### 11.3 Interpretation of the scientific sequence

The controlled system proves that intervention information can distinguish target answers under its model class. The benchmark establishes competitive prediction by the released core. The component analyses identify what the extensions can represent and where their gains occur. Functional and decision evaluations then measure particular uses. This is a sequence of complementary scientific questions, not an experimentally established causal mediation chain from network recovery through prediction to recommendation. Each stage has its own direct observable and comparison.

## Supplementary Note 12 | Source correspondence and computational reproduction

### 12.1 Manuscript components and scientific sources

This supplement provides the derivations, complete comparisons and reproduction instructions for the main manuscript. The full implementation archive is `reproduction/grn_repo_v27.zip`. The current five main and three supplementary figures are generated by `code/make_figures.py`; the archive's legacy `paper` task retains its own historical figure numbering. New-input core and functional interfaces call the original scientific modules. The functional example refits those modules on the archived Jurkat training input and compares the predictions with the saved output.

The package contains the main text, this supplement, figure exports, numerical sources, executable examples and verification scripts. Long comparison tables remain machine-readable, including benchmark metrics, RNA fold and rectangle losses, all GSE220974 recommendations and source-transfer results. The `evidence/` tree preserves source files under their original names; `CLAIM_SOURCE_MAP_V31.csv` connects the current claims to these files. Large raw matrices and neural assets are retrieved according to the full archive's external-asset instructions.

### 12.2 Reproduction levels

Four operations should be distinguished. Score replay starts from saved predictions, residuals or sufficient statistics. Weight replay starts from saved model parameters and required input assets. Compact refitting trains a declared model on retained condition-level data. Full raw-cell retraining additionally reconstructs preprocessing, graphs, embeddings, data splits and stochastic optimization. Table S24 records which operation each code entry actually performs. A hash check validates saved bytes; it does not by itself verify a scientific assumption.

From the extracted code archive, a task is invoked as:

```bash
python tools/reproduce.py --task gse220974_final --out /absolute/path/to/new_output
```

Other task names are listed in Table S24. Each output should be a fresh location. Large raw matrices and neural weights follow their external-asset manifests. The GSE220974 compact increment stores lossless numerical arrays and external weight identities; it does not include every raw-cell input or all model-weight bytes.

### Supplementary Table S24 | Complete code task map and actual reproduction scope

| Task name | Entry script | Recorded scope |
| --- | --- | --- |
| core | `compact/legacy_core/run.py` | Exact released per-operation metrics, paired ablations, rankings and historical bootstrap |
| risk | `compact/legacy_risk/run.py` | Compact moment/row replay of frozen registration, risk selection, metric calibration, bootstrap and intervals |
| support | `experiments/25_support_controls_and_witness_quantity/reproduce.py` | Unchanged experiment25 full endpoint-state, prediction-array and statistics reconstruction |
| growth | `tools/run_growth.py` | Actual K562 and Jurkat model fitting, selection, prediction and scoring from original compact condition inputs |
| step2 | `tools/run_step2.py` | Unchanged compact main RNA/structural score replay |
| rna_shrinkage | `tools/run_rna_shrinkage.py` | 11x4995 signed vectors; M0/M1/uniform/R comparisons |
| examples | `tools/run_examples.py` | Unchanged RNA/Jurkat example tables and four figures |
| calibration | `experiments/29_h1_calibration/run.py` | 168 scores,20 inner selections,80 independently reconstructed calibration operators |
| sf12 | `experiments/30_h1_sf12_qualified/run.py` | 60 SF12 coefficients re-solved; five qualified tasks and rejected case_R reproduced |
| dixit | `experiments/31_dixit_confirmation/run.py` | Primary RNA archived summary audit; separate post-confirmation G2M actual refit |
| combiseq | `experiments/32_combiseq_confirmation/run.py` | Frozen program actual refit from compact count matrices |
| paper | `publication/run.py` | 7 main +10 ED figures redraw;19 supplementary table references; exact scalar checks |
| recent_bridge | `experiments/33_recent_bridge_and_acquisition_20261031/reproduce.py` | Nine 2026-10-31 bridge/acquisition rounds: lossless source restoration, declared numeric replay and exact cache regeneration; per-round scope and verification are recorded inside experiment33. |
| gse220974_final | `experiments/34_gse220974_five_model_final/reproduce.py` | Closed GSE220974 experiment: exact float64 residual/utility replay, five main models, 420 choices, six paired intervals; external full prediction/weight replay identities retained. No training. |

The `recent_bridge` task can restore nine archived rounds and regenerate the declared deduplicated cache families. Merely extracting its stored files does not perform cache regeneration or repeat all fitting. The `dixit` task audits the frozen primary RNA summary and separately refits the later G2M development model; these are different operations. The legacy `paper` task's seven main figures, ten extended-data figures and nineteen supplementary-table references retain their v22 numbering.

### 12.3 Numerical and interface verification

The verification scripts check source hashes, table values and correspondence between figures and numerical inputs. Benchmark checks reconstruct 654 operation ranks from six metric ranks and retain all sixteen entries and fifteen comparator intervals. RNA budgets are calculated from the 11 × 4,995 signed prediction arrays. Application checks include all 28 recommendation tasks and their design-label matches. The manifest records a hash for each delivered file.

The core example fits the original core modules to synthetic input. The functional example fits native modules to archived training measurements and reproduces all 617 query predictions. The fitting interfaces reject query outcomes; separate tests check native-call parity and leakage rejection. Code and licence scope are described in Note 14.4, and public data accessions are listed in Data availability.

## Supplementary Note 13 | Target-specific decomposition, correction budgets and complete application accounting

### 13.1 Frozen benchmark decomposition

The benchmark decomposition retains the released eligible method set, top-100 panel, six metrics, tie convention and operation-weighted aggregation. Each stored operation rank is checked against the rank reconstructed from its six metric ranks. Rows are then divided into doubles and unseen singles using both their identity and released split metadata. This is a retrospective reanalysis, not a newly preregistered benchmark. No scores, model fits or target labels are adjusted.

### Supplementary Table S25 | Complete benchmark separated by target type

| Method | All 654: rank (position) | Doubles 540: rank (position) | Singles 114: rank (position) |
| --- | --- | --- | --- |
| WitnessCell | 3.803 (1) | 3.515 (1) | 5.167 (2) |
| Scouter | 5.411 (2) | 5.150 (2) | 6.649 (4) |
| linearModel | 5.725 (3) | 6.180 (4) | 3.570 (1) |
| baseReg | 6.709 (4) | 6.174 (3) | 9.246 (12) |
| trainMean | 6.771 (5) | 6.315 (5) | 8.930 (11) |
| GEARS | 7.651 (6) | 7.754 (6) | 7.167 (5) |
| scGPT | 7.806 (7) | 7.881 (7) | 7.447 (7) |
| scELMo | 7.916 (8) | 7.957 (8) | 7.719 (9) |
| AttentionPert | 8.058 (9) | 8.050 (9) | 8.096 (10) |
| GeneCompass | 8.367 (10) | 8.531 (10) | 7.588 (8) |
| GenePert | 8.572 (11) | 9.037 (11) | 6.368 (3) |
| CPA | 9.251 (12) | 9.659 (12) | 7.316 (6) |
| biolord | 9.846 (13) | 9.850 (13) | 9.825 (13) |
| scFoundation | 11.690 (14) | 11.826 (14) | 11.044 (14) |
| baseMLP | 13.018 (15) | 12.515 (15) | 15.404 (16) |
| baseControl | 13.694 (16) | 13.919 (16) | 12.632 (15) |


The double subset contains 287 dataset–pair clusters; singles contain 97 dataset–single clusters. Cluster counts differ from operation–seed counts because one identity may occur in multiple splits. These are benchmark units, not independent biological replicates.

### Supplementary Table S26 | Exposure and dataset decomposition with paired intervals

| Population | Records | Clusters | WitnessCell | Scouter | Paired difference [95% interval] |
| --- | --- | --- | --- | --- | --- |
| double | 540 | 287 | 3.515 (1) | 5.150 (2) | -1.635 [-2.117, -1.141] |
| single | 114 | 97 | 5.167 (2) | 6.649 (4) | -1.482 [-2.492, -0.469] |
| combo seen0 | 43 | 43 | 4.279 (1) | 6.116 (3) | -1.837 [-3.233, -0.395] |
| combo seen1 | 374 | 242 | 3.497 (1) | 5.291 (2) | -1.794 [-2.348, -1.232] |
| combo seen2 | 123 | 112 | 3.301 (1) | 4.382 (2) | -1.081 [-1.807, -0.336] |
| Norman double | 208 | 118 | 3.317 (1) | 6.139 (4) | -2.822 [-3.476, -2.156] |
| Wessels double | 181 | 106 | 3.718 (1) | 4.680 (2) | -0.961 [-1.710, -0.172] |
| Schmidt double | 60 | 28 | 6.033 (4) | 2.500 (2) | 3.533 [2.558, 4.661] |
| Replogle exp6 double | 91 | 35 | 1.901 (1) | 5.571 (5) | -3.670 [-4.512, -2.782] |


Values in the method columns are mean operation rank (position among sixteen methods). The last column is WitnessCell minus Scouter, with 20,000 cluster-bootstrap draws. Lower values favor WitnessCell. Resampling carries all repeated seed records of a dataset–identity together and divides each draw's summed ranks by its summed record count. The seen0/seen1/seen2 state follows the benchmark's endpoint-exposure metadata, not a new post hoc definition. The same pair can appear in different exposure states in different splits, so subgroup intervals are descriptive and not independent confirmations. All fifteen paired comparator intervals are retained in `revision/BENCHMARK_PAIRED_INTERVALS.csv`; all sixteen ranks in every population are in `BENCHMARK_DECOMPOSITION.csv`. Dataset-specific MSE and correlation comparisons in `BENCHMARK_DATASET_ABSOLUTE.csv` retain their original scales and include unfavorable outcomes; raw MSEs are not pooled across datasets.

### 13.2 Positive, adverse and net correction contributions

For each readout, error reduction is summed across all eleven same-fit rectangles before its sign is taken. Positive and adverse sums therefore characterize the readout-aggregated budget. Their difference exactly recovers the net SSE reduction. This decomposition does not imply that every rectangle or every coordinate improves, and neither the readout counts nor the eleven overlapping rectangles create independent replication.

### Supplementary Table S27 | Error budget for M1, uniform shrinkage and R

| Model | Readouts | Positive / all M0 (%) | Adverse / all M0 (%) | Net / same-set M0 (%) |
| --- | --- | --- | --- | --- |
| M1 | all | 20.8270 | 11.3887 | 9.4382 |
| M1 | fixed ten | 17.1224 | 0.0000 | 73.1778 |
| M1 | remaining | 3.7046 | 11.3887 | -10.0314 |
| scalar | all | 20.5327 | 10.8339 | 9.6988 |
| scalar | fixed ten | 16.8661 | 0.0000 | 72.0823 |
| scalar | remaining | 3.6666 | 10.8339 | -9.3566 |
| R | all | 17.6425 | 1.6802 | 15.9622 |
| R | fixed ten | 16.2532 | 0.0000 | 69.4628 |
| R | remaining | 1.3893 | 1.6802 | -0.3798 |


The all-readout M0 SSE is 232.70214771649012. The positive and adverse columns use this same denominator in every row; the final column uses the relevant subset's M0 SSE. Thus subset percentages in the positive/adverse columns add to the all-readout values, whereas subset net percentages need not. `scalar` denotes the training-selected uniform shrinkage control. The ten fixed readouts are S100A9, S100A8, S100A12, LYZ, VCAN, VIM, MS4A6A, MYL9, LGALS1 and PPP1R27, identified during model development and fixed for the deletion sensitivity analysis.

For all readouts, M1 has positive/adverse SSE contributions 48.46481740654159/26.501830662818197, and R has 41.054402489037855/3.9099455447293474. R's total positive error reduction is 84.7097% of M1's, and its total adverse increment is 85.2465% lower. Each model determines its own positive and adverse readout sets, so the first ratio compares aggregate magnitudes rather than matched-readout retention. The signs and aggregate budgets are descriptive on the original development panel; the result explains the net improvement without expanding its cross-source scope. `RNA_ALL_READOUT_CONTRIBUTIONS.csv` retains every model–readout row, including negative contributions, and `RNA_ALL_RECTANGLE_GAINS.csv` retains all rectangles.

### 13.3 Complete task outcomes and original-construct audit

### Supplementary Table S28 | All-task decision outcomes and identity qualification

**a. Aggregate task outcomes.** Wins, ties and losses compare R-G's achieved held-out utility with the indicated reference on the identical task; ties use an absolute tolerance of $10^{-12}$. The entire set of 28 tasks is retained.

| Reference | Mean regret | Hits / 28 | R-G wins / ties / losses | Same choice / 28 |
| --- | --- | --- | --- | --- |
| R_G | 0.006197 | 19 | 0 / 28 / 0 | 28 |
| StandardRidge_RNA | 0.006591 | 18 | 1 / 27 / 0 | 27 |
| Scouter_Ensemble | 0.006623 | 17 | 2 / 26 / 0 | 26 |
| GEARS_Ensemble | 0.010258 | 18 | 6 / 17 / 5 | 17 |
| SingleB_additive | 0.012220 | 18 | 5 / 20 / 3 | 20 |
| PRESAGE_Combo_seed24 | 0.006197 | 19 | 0 / 28 / 0 | 28 |


**b. Every task's regret.** Losses are multiplied by $10^3$. Rows are ordered by physical anchor, mode-specific anchor and fold. They are not ordered or selected by model performance. Mode suffixes a/i denote activation/repression. The five principal methods use the same pools; their median regrets are all zero, illustrating the importance of the nonzero tail and complete task distribution.

| Task | Candidates | PairPick | Linear | Scouter | GEARS | Additive |
| --- | --- | --- | --- | --- | --- | --- |
| f0:CEBPB-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f0:CEBPB-i | 2 | 1.020 | 1.020 | 1.020 | 1.020 | 0.000 |
| f0:FOSL1-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f1:FOSL1-a | 2 | 0.000 | 0.000 | 8.381 | 8.381 | 0.000 |
| f2:FOSL1-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f3:FOSL1-a | 2 | 0.000 | 0.000 | 0.000 | 11.130 | 0.000 |
| f4:FOSL1-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f4:FOSL1-i | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f1:GATA1-a | 2 | 5.969 | 5.969 | 5.969 | 0.000 | 5.969 |
| f2:GATA1-a | 3 | 21.339 | 21.339 | 21.339 | 0.000 | 21.339 |
| f0:GATA1-i | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f2:GATA1-i | 2 | 15.453 | 15.453 | 15.453 | 15.453 | 0.000 |
| f3:GATA1-i | 2 | 0.000 | 0.000 | 3.533 | 3.533 | 3.533 |
| f0:JUN-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 15.499 |
| f1:KLF1-a | 6 | 49.230 | 49.230 | 49.230 | 0.000 | 89.451 |
| f2:KLF1-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f4:KLF1-a | 2 | 12.876 | 12.876 | 12.876 | 0.000 | 12.876 |
| f1:KLF1-i | 4 | 0.000 | 11.028 | 0.000 | 0.000 | 0.000 |
| f4:KLF1-i | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f2:LYL1-a | 3 | 27.173 | 27.173 | 27.173 | 27.173 | 27.173 |
| f2:LYL1-i | 3 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f3:LYL1-i | 2 | 31.083 | 31.083 | 31.083 | 0.000 | 0.000 |
| f0:SPI1-a | 5 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f2:SPI1-a | 2 | 0.000 | 0.000 | 0.000 | 37.111 | 0.000 |
| f4:SPI1-a | 4 | 0.000 | 0.000 | 0.000 | 157.352 | 156.676 |
| f2:YY1-a | 2 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| f3:YY1-a | 5 | 0.000 | 0.000 | 0.000 | 16.695 | 0.256 |
| f4:YY1-a | 3 | 9.385 | 9.385 | 9.385 | 9.385 | 9.385 |


**c. Comparison with the published designed-construct table.** The original article is Pacalin et al., doi:10.1038/s41587-024-02213-3. Its Supplementary Table 2, `SuppTable2 - Double gRNA K562`, was obtained from the Europe PMC supplementary archive for PMC11569274 on 7 November 2026. The retained XLSX is the published supplement, SHA-256 `705173fc435f7a596c8cb851abda30a7725e60505ec93b82d7e767faaace54d2`. It is attributed to the original authors and supplied under the article's CC BY 4.0 terms; the audit tables are derived mappings.

| Detected label status | Count | Interpretation |
| --- | --- | --- |
| Exact gene-mode design match | 15 | Matches a listed pair after stripping guide index; per-cell identity still unverified |
| Enhancer-only correspondence | 1 | GATA1-a + KLF1-i only matches an enhancer-specific GATA1 design |
| Not listed as designed pair | 33 | Detected label lacks a corresponding designed pair in Supplementary Table 2 |


Only terminal guide indices are stripped; pair order is canonicalized and enhancer tags are preserved. The gene-mode design set is therefore distinct from the detected-label cohort. Matching 15 labels does not validate their cell-level guide identities; the remaining statuses do not determine the cause of the mismatch. The possibilities include incidental co-delivery, assignment errors and other assay-specific labeling effects. No outcome-based filtering, cohort substitution or reranking is performed. All 49 mappings, their exact source-table rows and the complete designed constructs are retained in `GSE220974_IDENTITY_AUDIT.csv` and `GSE220974_DESIGNED_CONSTRUCTS.csv`. The selection experiment supports an operational retrospective demonstration on detected labels, not independent confirmation of intended combinations or their molecular interactions.

### 13.4 Input-based configuration choice and new-input software

### Supplementary Table S29 | Choosing a configuration from inputs and objective

| Input and objective | Configuration | Required support and delivered output |
| --- | --- | --- |
| Normalized RNA; complete response prediction | Released core | Control, measured training singles, training/validation combinations and annotations for missing endpoints; predicted condition means and baseline |
| Matched RNA; additional partner contrasts | R | Measured endpoint effects, training cell halves and adequate partner overlap; readout-weighted residual correction |
| Fixed anchor and transcriptional objective | PairPick (R-G) | Mode-aware candidate predictions and training-derived objective; utilities and partner choice |
| Repeated functional measurements | Functional rank-one model | Single/double functional support and repeat-specific baseline; growth response plus shared residual; spectral correction is secondary |
| Existing shared functional function; target singles and combination calibration data | Frozen SF12 or V4 adapter | Source transformations and weights; target additive baseline; a scalar strength for each specified source function, with all evaluated configurations retained |
| Existing predictions and internal held-combination fits | Component control | Fixed query panel, M0 and D; training-derived coefficient; endpoint-preserving adjusted prediction |
| Existing shared partner function and target single/combination measurements | Full-screen target adapter | Fixed additive baseline and frozen direction; coefficient from calibration pairs; predicted held-out responses |

The input conditions select the analysis family before query outcomes are seen. Within-family model selection uses permitted training or validation data. This table is a user guide, not a validated performance-optimal routing algorithm. Inadequate contrast geometry produces an undefined partner diagnostic rather than proof of absent biology. The fitting interfaces cover the released core and functional rows; the control interface accepts their saved predictions and internal held-combination outputs. The RNA R and mode-aware R-G fitting/replay entries remain in the complete included archive. The real functional example reuses existing Jurkat inputs and is a software-access check, not a further biological evaluation.

From the package root, run:

```bash
python -m pip install -r software/requirements.txt
python software/witnesscell_predict.py --input software/example_input.json --output my_predictions
python software/test_interface.py
python code/analyze_revision.py
python code/analyze_applications.py
python code/audit_construct_identity.py
python code/make_figures_v36.py
python software/test_growth_interface.py
python code/verify_v36.py
```

The adapter takes JSON condition summaries with explicit train, validation and query identities. Query outcomes are rejected; reversed-pair leakage, malformed vectors and invalid variances are checked. Seven unchanged native modules supply the released `multihead_v14` implementation. The synthetic example has five new queries and a missing-single endpoint; it is an interface example rather than biological evidence. The output includes predicted normalized means, the nonlinear baseline, selected parameters and an input hash. Native-call parity, absence of query truth and a common-input-shift test are recorded in `software/INTERFACE_TESTS.json`. Dependencies and normalization assumptions are documented in the software README. Figures and supplementary audits require the additional dependencies in `code/requirements.txt`.

## Supplementary Note 14 | Evidence status, baseline correspondence and software access

### 14.1 Evidence status and experimental units

### Supplementary Table S30 | Development, freezing and source relationships

| Result | Development | Freeze / separation | Source relationship | Role |
| --- | --- | --- | --- | --- |
| Causal system | Specified simulation | Declared family and layout | Synthetic | Information principle |
| Core benchmark | Public benchmark development | Saved predictions; retrospective split analysis | Four sources | Prediction competition |
| RNA R | Development panel | Outer pairs held out; panel informed design | Wessels shared | Construction and scope |
| K562 growth | Revised development | Original registered failure retained | Same growth study | Support and increment |
| Jurkat growth | Procedure transfer | Before fitting/scoring; coefficients refitted | New cell context, same study | Functional partner gain |
| Original Shen growth | Source evaluation | Original procedure and predictions frozen before scoring | Independent study; three shared-design backgrounds | Complete gain, adverse partner increment |
| Component control | Retrospective development | Coefficients from internal arrays; historical outcomes already visible | Existing Horlbeck and Shen data | Targeted adjustment and usage conditions |
| Han full-screen calibration | Shared-function extension after earlier Han access | Base weights fixed; target coefficients fitted on calibration pairs; query predictions locked before scoring | One target study; two related readout panels | Test target-side strength selection with all six portable configurations |
| Dixit / Combi-seq | Source evaluations | Frozen protocols and predictions | Other studies | Transfer boundaries |
| Risk selection | Separate calibration | Fixed threshold; disjoint pairs | Existing sources | Selective-use outcome |
| GSE220974 | Retrospective comparison | Not globally untouched | Other study; identity limits | Recommendation operation |
| New-input examples | Software checks | Unchanged estimator code | Synthetic/core; archived/growth | Executable access |

The original machine-readable descriptions are in `revision/EVIDENCE_STATUS.csv`; extension status is recorded in `partner_control/PROTOCOL.json` and Table S42. Freezing a procedure, locking predictions and using a new biological source are distinct attributes. Jurkat is the strongest positive functional procedure-transfer result here; it is not an independent-study replication. Dixit and Combi-seq retain their independent-source evaluation status despite nonpositive partner increments. The Jurkat preparation amendment concerned missing GINS4 single measurements and preceded fitting and outcome summaries. The recorded separation is local procedural separation, not independent third-party blinding. 

The RNA task has eleven same-fit rectangles over 4,995 readouts, with shared pairs, endpoints and fitted models. Counts of nine or eleven improved rectangles are descriptive; an independent-binomial sign test is not justified. The twelve released benchmark query graphs' 212 algebraic contrast dimensions are from a different evaluation. Neither those dimensions nor the RNA task's fifteen projected dimensions are biological replicate counts. Readout concentration and annotation are not used to estimate a number of independent biological programs.

### 14.2 Resolve the baseline claim at the level actually observed

### Supplementary Table S31 | Distinct baseline and scoring objects

| Object | Evidence retained | Permitted interpretation |
| --- | --- | --- |
| Self-contained weighted combination reference | Archived training-mean additive predictions and weighted scores | The 24.87% improvement is relative to this named reference only |
| Released implementation-level correspondence | Earlier code audit identified a combination-branch correspondence | Does not establish identical historical adapted inputs or scores |
| Official `linearModel` and `trainMean` rows | 540 common Top-100 double records; 529 different MSEs; maximum difference 0.0235 | Preserve both published benchmark rows; they are not score-identical |
| Exact historical scoring objects | Complete matched cell-level adapted inputs and score-generation states are unavailable here | The origin of the score differences is not adjudicated |

The scoring discrepancy persists in double perturbations alone and cannot be explained by inclusion of unseen singles in the full track. `BASELINE_PROVENANCE_AUDIT.json` records the scoring-file hash, comparison denominator and unresolved historical scoring object. The available record does not identify which adaptation, sampling, preprocessing or file-generation step caused the difference. Both official scoring rows are preserved, and Table S31 evaluates the sensitivity to their inclusion under the unchanged formal ranking rule.

`code/audit_v28.py` generates `BASELINE_SCORE_ALIGNMENT.csv` directly from the Top-100 saved operation rows, using exact dataset and operation keys. It also extracts all fifteen double-only intervals, each below zero. These are pointwise 95% intervals conditional on the existing score archive, not simultaneous intervals or a global multiplicity-adjusted test. Dataset-pair clustering retains repeated seed records, but does not account for every shared-endpoint dependence, the full development history or variation across future studies.

**Comparator-inclusion sensitivity (extension to Table S31).** All 540 Top-100 double-perturbation operation records are retained. For each reduced method set, the six archived metric scores are reranked within each operation using competition (minimum) ranks; their mean is reranked to form the operation rank, which is then averaged across records. The unchanged sixteen-entry case exactly reconstructs all 51,840 saved metric ranks and all 540 WitnessCell-and-comparator operation rank vectors.

| Removed entries | Retained entries | WitnessCell mean operation rank | Position |
| --- | --- | --- | --- |
| None | 16 | 3.514815 | 1 |
| linearModel | 15 | 3.340741 | 1 |
| trainMean | 15 | 3.366667 | 1 |
| Both | 14 | 3.179630 | 1 |

WitnessCell retains first place under each exclusion. Mean ranks from different-sized method sets are not directly comparable measures of improved prediction. This fixed-score sensitivity does not replace the original ranking, resolve the origin of historical scoring differences or add independent biological validation. `code/analyze_v28.py` generates the complete retained-method tables and operation ranks in `revision/BASELINE_EXCLUSION_SENSITIVITY.csv` and `revision/BASELINE_EXCLUSION_OPERATION_RANKS.csv`; `revision/BASELINE_EXCLUSION_AUDIT.json` records the score-source hash, tie rules and reconstruction checks.

### Supplementary Table S32 | Native-scale double-perturbation MSE versus Scouter

| Dataset | WitnessCell MSE | Scouter MSE | Relative MSE reduction (%) | WitnessCell aggregate position |
| --- | --- | --- | --- | --- |
| Norman | 0.056440 | 0.064038 | 11.8650 | 1 |
| Replogle_exp6 | 0.009996 | 0.018109 | 44.8025 | 1 |
| Wessels | 0.062142 | 0.052071 | -19.3415 | 1 |
| Schmidt | 0.016503 | 0.012792 | -29.0163 | 4 |

These raw means are from the same frozen Top-100 operation archive; they are not pooled across studies. Their ordering need not match six-metric, operation-rank aggregation. Four metrics lead in the full-track profile, whereas Wasserstein and common-DE ranks are weaker. The common mean-based adapter does not alone explain method-specific differences. No unverified biological mechanism is assigned to this metric asymmetry.

### 14.3 Recommendation task size and identity qualification

### Supplementary Table S33 | Candidate pools and identity of selected labels

| Candidates per task | Number of tasks |
| --- | --- |
| 2 | 19 |
| 3 | 4 |
| 4 | 2 |
| 5 | 2 |
| 6 | 1 |

Under uniform candidate choice and one designated empirical best per task, the summed expectation is 11.9 hits across 28 tasks. This is an analytical descriptive reference, not an independence-based significance test. The small pools can support meaningful decisions while offering limited resolution among predictors that usually choose the same candidate.

PairPick R-G selected an unlisted-as-designed label in 18 tasks and a direct gene-mode design match in ten tasks; these are task counts, not counts of independent constructs. The full mapping for all primary methods and PRESAGE is in `CHOSEN_LABEL_STATUS.csv`. Absence from the original design list does not prove a detected label is false. Conversely, detected-label prediction and disclosure do not independently establish intended cell-level intervention identity. The original full cohort remains intact. Supplementary Fig. S1 shows all recommendation outcomes and a separately evaluated risk-selection task; no combined pipeline is inferred.

### 14.4 New-input functional fitting with native code

### Supplementary Table S34 | Executable interfaces and verification scope

| Interface | Inputs | Outputs | Verified scope |
| --- | --- | --- | --- |
| Released RNA core | Normalized means, variances, counts, annotations and train/validation/query identities | Predicted means, nonlinear baseline and selected settings | Synthetic 24-readout example; native parity and input isolation |
| Functional procedure | Gene identities, two-repeat singles, training pair responses and outcome-free queries | Two-repeat baseline, rank-one prediction and increment; secondary spectral prediction | Real archived 2,464-pair training input and 617 query identities; exact archived-prediction parity |
| RNA R and PairPick R-G | Task-specific processed inputs in full code archive | Partner corrections and mode-aware objective choices | Original experiment entry points; not replaced by either interface above |

The functional wrapper calls three unchanged native modules: `frozen_models.py`, `spectral_models.py` and `run_smoke.py`. Their hashes are recorded in `software/GROWTH_INPUT_PROVENANCE.json`. Three inner folds are assigned by the original gene-pair identity hash; every inner subset refits the baseline. The unchanged regularization grid, zero-correction alternative and rank-one primary/spectral secondary roles are retained. The wrapper rejects duplicate or reversed train/query pairs, query response objects and nonfinite measurements. Prediction uses only the supplied JSON input and native modules; expected archived predictions are opened by the separate verification script, never by the fitter.

The included functional example contains real already-processed measurements from the existing Jurkat evaluation: 2,464 training pairs and 617 query pairs exhaust the 3,081 unordered pairs among 79 endpoints. The approximately 80% training fraction is a tested setting, not a validated minimum input requirement. Predictions complete a well-supported endpoint set; this example does not test sparse-combination or unseen-endpoint functional extrapolation. Re-executing it reproduced the baseline, rank-one and spectral query predictions with maximum absolute difference zero in the supplied environment. This supports software access and native parity, not a new biological confirmation. Both adapters inherit the input units and preprocessing; the functional wrapper requires two measured repeats and is not a raw-count processor or mixed-mode RNA interface.

```bash
python software/witnesscell_growth.py --input software/growth_example_input.json --output my_growth_predictions
python software/test_growth_interface.py
```

Use a new output directory to preserve existing runs. The fitted baseline and corrected prediction allow a user to inspect the incremental contribution before subsequent experimental measurement.

The full `grn_repo_v27.zip` archive is included under `reproduction/`. It contains the scientific implementations, compact numerical evidence and retrieval indexes for external assets. The interface code and tests use the included MIT licence; upstream scientific modules retain their original rights. No unified redistribution licence is included. `RIGHTS_AND_RELEASE.md` records the scope for each component.

### 14.5 Same-fit functional four-corner illustration

The descriptive four-corner example fixes A = TTI2, B = PGK1 and D = OGFOD1. Candidate backgrounds C are enumerated using pair identities alone. ACTR10 is alphabetically first among ACTR10, MRPL32, NDUFA8 and POLD3, each of which supplies AB, AD, CB and CD in the same held-out Jurkat fit. The additional background is therefore chosen without ranking prediction success. The previously selected anchor/partner example remains descriptive. `revision/FOUR_CORNER_SELECTION_V28.json` records this rule and the exact source hash before the new contrast calculation.

For each repeat and each observed/predicted response, the signed contrast is AB − AD − CB + CD. The four endpoint contributions cancel algebraically. The mean row averages response contrasts over the two measured repeats; it is not a third repeat. Its squared error is not substituted for the pooled-repeat loss used in the primary 617-pair evaluation. No hypothesis test or population interval is attached to these examples. They overlap in AB and AD and are not four independent confirmations.

### Supplementary Table S35 | Every eligible continuation of the fixed Jurkat example

Growth responses use control-centred log2 enrichment normalized by population doublings. The main figure uses ACTR10. All backgrounds use the same fixed TTI2 background and PGK1/OGFOD1 partner exchange. Numerical precision is retained in `revision/JURKAT_FOUR_CORNER_CONTRASTS_V28.csv`; `revision/JURKAT_FOUR_CORNERS_V28.csv` contains all four raw response corners in both repeats.

| Second background | Repeat | Observed contrast | Baseline contrast | Rank-one contrast |
| --- | --- | --- | --- | --- |
| ACTR10 | 1 | -0.059209 | -0.020570 | -0.040740 |
| ACTR10 | 2 | -0.056903 | -0.026049 | -0.046219 |
| ACTR10 | mean | -0.058056 | -0.023309 | -0.043479 |
| MRPL32 | 1 | -0.046414 | -0.020785 | -0.042083 |
| MRPL32 | 2 | -0.052939 | -0.025117 | -0.046416 |
| MRPL32 | mean | -0.049676 | -0.022951 | -0.044250 |
| NDUFA8 | 1 | -0.100244 | -0.037652 | -0.065737 |
| NDUFA8 | 2 | -0.091610 | -0.038361 | -0.066446 |
| NDUFA8 | mean | -0.095927 | -0.038007 | -0.066091 |
| POLD3 | 1 | -0.152039 | -0.043210 | -0.068571 |
| POLD3 | 2 | -0.184849 | -0.046797 | -0.072158 |
| POLD3 | mean | -0.168444 | -0.045003 | -0.070364 |

All four repeat-mean contrasts move towards the observed response after correction. Their role is to illustrate the learned response component; the choice of these anchors was inherited from an already examined example. The primary performance evidence remains the complete, unchanged held-out test.

### 14.6 Relation to measurement controls and established estimators

Technical duplicates in Miller et al. form evaluation positive controls and assess metric sensitivity.[12] The RNA weighting rule here uses only training halves to attenuate a fitted correction, leaving the scored outcome axis unchanged. Shared batch, guide or preparation effects can persist in both halves, so reproducibility of the correction across independent biological preparations is not established by these weights. Functional two-repeat validation tests a separately regularized predictor.

The constituent estimation tools have established precedents: regularized kernel expansions,[39] low-rank matrix completion,[40] spectral thresholding and incomplete-matrix regularization.[41,19] The manuscript specifies the projected kernel and signed symmetric correction actually used; it does not import exact-recovery guarantees from these precedents. The causal result is conditional on its declared mechanism family,[13] and bootstrap summaries retain their explicitly stated sampling units.[42]

## Supplementary Note 15 | Additional quantitative results

This note reports the effect sizes and sensitivity analyses supporting the principal findings. It summarizes the same comparisons documented in the preceding tables, including their denominators, model choices and interpretation.

### 15.1 Controlled-system measurement budget and uncertainty

Two panels add four joint interventions to control and seven singles: the informative panel uses $AB,BC,DE,DF$, and the comparison panel uses $AD,BD,BE,CF$. They match every endpoint degree, the number of conditions and the observation budget. Each condition contains 256 observations, giving 3,072 per panel. The four query pairs $AC,AF,CE,EF$ are absent from both panels; their target is $\Delta=\mu_{AC}-\mu_{AF}-\mu_{CE}+\mu_{EF}$ (Fig. 2b).

Four generating worlds have identical training distributions in the comparison panel but targets $0,-0.5,0,-0.5$. Any estimator based on that panel has equal-world squared-error risk $\mathbb E[(T+0.25)^2]+1/16\geq0.0625$ (Fig. 2c). Repeating those conditions more precisely preserves the ambiguity. The informative panel exposes the relevant gate directions. Across 1,000 paired noise seeds, the same candidate least-squares procedure reduces target MSE from 0.065295 to 0.004674, a 92.84% reduction; the paired difference is −0.060620, with 95% Monte Carlo interval [−0.061103, −0.060137] (Fig. 2d).

Contrast GLS achieves MSE 0.005247 on the informative panel, with an unresolved difference from candidate least squares. The risk reduction therefore demonstrates the value of information layout, not an estimator-exclusive advantage. Allowing arbitrary simultaneous gate contributions introduces a function that vanishes on the informative training panel while changing the query. Identification is lost in that enlarged family (Supplementary Note 2).

### 15.2 Complete benchmark estimates, reference definitions and support analyses

WitnessCell ranked first in the full track, with mean operation rank 3.803, followed by Scouter at 5.411 and the released linear model at 5.725. In the double-only subset, WitnessCell retained first place with mean rank 3.515, followed by Scouter at 5.150; the linear model and GEARS had mean ranks 6.180 and 7.754 (Fig. 3a). Clustering repeated seed records by dataset and pair identity gave a WitnessCell–Scouter difference of −1.635, with a 95% bootstrap interval of −2.117 to −1.141. Corresponding differences were −2.665 against the linear model and −4.239 against GEARS, with intervals of −3.002 to −2.317 and −4.764 to −3.702 (Fig. 3c). All fifteen double-only comparator differences had pointwise 95% cluster-bootstrap intervals below zero. These intervals describe the fixed scoring archive and are not simultaneous, multiplicity-adjusted confidence statements (Supplementary Fig. S2).

Endpoint exposure provided a further decomposition. WitnessCell ranked first by point estimate in each double-perturbation subgroup: neither endpoint seen, one seen and both seen, containing 43, 374 and 123 records, respectively (Fig. 3b). The subgroups share some underlying pair identities across split seeds and are not independent confirmations. In the 114 unseen-single records, WitnessCell ranked second at 5.167 behind the linear model at 3.570. The overall combination result therefore coexists with a distinct ordering for single-response generalization. Full sixteen-method tables, all reference contrasts and dataset-specific raw-loss estimates are provided in Supplementary Tables S25–S26 and their machine-readable sources.

Dataset and metric decompositions delimit the aggregate advantage. WitnessCell ranked first in Norman, Wessels and Replogle_exp6 and fourth in Schmidt. Relative to Scouter, mean native-scale double-perturbation MSE decreased by 11.87% and 44.80% in Norman and Replogle_exp6, but increased by 19.34% and 29.02% in Wessels and Schmidt (Fig. 3d). Mean raw losses and aggregated operation ranks measure different aspects of performance. Across the full track, WitnessCell led the correlation, MSE, energy-distance and symmetric-KL rankings, while ranking ninth for Wasserstein and eighth for common differential expression. The mean-based distribution adapter is shared evaluation machinery; it does not by itself explain this asymmetry or establish full cell-distribution prediction.

A separate weighted-error analysis of 540 combinations gave a 24.87% reduction against the archived training-mean additive reference, with loss ratio 0.7513 and pair-bootstrap interval 0.7031–0.8019. This named reference is used once. The official `linearModel` and `trainMean` rows remain distinct: 529 of 540 double-perturbation MSEs differ, and the historical scoring difference is unresolved. Removing either row or both and recomputing all metric and operation ranks leaves WitnessCell first in the double-only subset (Supplementary Table S31). This comparator-inclusion sensitivity supplements the unchanged formal ranking.

These gains combine several useful operations. The single-response head estimates unavailable endpoints; nonlinear composition accounts for response saturation; and measured combinations refine shared endpoint contributions. After removing manipulated target coordinates, complete non-target MSE ratios against the shared-endpoint additive reference were 0.9049 in Norman and 0.6459 in Wessels, corresponding to reductions of 9.51% and 35.41%. An independently selected saturation comparator still left a 32.36% MSE reduction in Wessels. The advantage therefore extends to downstream readouts rather than being confined to the directly targeted genes.

In Wessels, increasing the retained residual library from one quarter to one half and all available training combinations increased MSE gain over the paired nonlinear backbone from 19.10% to 27.42% and 33.25%; weighted gains were 26.81%, 33.91% and 37.66% (Supplementary Fig. S3). The full-minus-quarter MSE gain was 14.15 percentage points, with a paired interval of 12.42–16.06. This comparison holds the backbone fixed and measures the return from additional combination support.

### 15.3 Functional budgets, repeat-specific results and secondary models

For eighty K562 genes, low and high training budgets contained 421–422 and 2,317–2,318 pairs per fold. The 632 evaluated pairs stayed fixed across the three folds. Increasing support raised the rank-one correction's complete-response gain from 0.27% to 11.03%, while its partner gain rose from −0.033% to 11.68%. The high-budget spectral correction gave complete and partner gains of 11.39% and 10.67%. The baseline itself improved by 10.12% when supplied with additional observations. The data therefore supported separate benefits for the baseline and for the added partner structure (Fig. 4a).

These budgets change condition count, endpoint degree and layout coverage together. Their value lies in demonstrating that more extensive combination measurements can support a transferable residual prediction on fixed queries. The controlled mechanism experiment isolates layout information at equal budget; the growth experiment tests the return from expanding real experimental support. Considering both avoids attributing every biological data-volume effect to a single structural property of the design.

The procedure was frozen before evaluation in Jurkat, using 79 genes and 3,081 possible unordered pairs: 2,464 training pairs (approximately 80%) and 617 held-out pairs (approximately 20%). Coefficients were estimated from Jurkat training data, testing procedure portability to another cellular environment within the same study. The primary rank-one correction reduced pooled growth-response error by 6.43% and endpoint-eliminated partner error by 6.18%. Complete gains were 4.45% and 8.88% in the two measured repeats. The secondary spectral correction achieved gains of 8.16% and 5.90%, respectively (Fig. 4b). Removing each endpoint from scoring left complete rank-one gains between 4.95% and 7.10%.

### 15.4 Complete RNA comparison and correction budget

The unweighted partner extension reduced held-out rectangle error by 9.44% relative to $M_0$. Adding measurement-derived shrinkage, denoted $R$, increased the reduction to 15.96%, with complete-response error reduced by 9.28% (Fig. 5a,b). The shrinkage predictor improved contrast error in nine of eleven rectangles against $M_0$ (Fig. 5c) and in all eleven against both the unit-weight extension and a training-selected uniform shrinkage control. Relative to the latter control's error, the reduction was 6.94%. All five outer folds had positive aggregate gains, including one close to zero.

Matched GEARS and Scouter runs used the same outer combination identities and evaluation mask, each with three optimization seeds.[8,18] The $R$ extension reduced rectangle error by 10.81% against the GEARS seed-mean error and by 12.88% against Scouter. These comparisons aggregate seed-specific losses, whereas the later PairPick application averages prediction vectors before scoring. Inner-selection budgets also differed between the recorded RNA workflows. Optimization seeds characterize fitting variation on this shared outer task.

A signed error budget explains how measurement-informed amplitude control improves the correction (Fig. 5d). Readout-level reductions in squared error, summed over all eleven rectangles, were separated into positive reductions and adverse increments. With the same total $M_0$ error as denominator, $M_1$ contributed 20.83% positive reduction and 11.39% adverse increment, yielding 9.44% net improvement. For R these quantities were 17.64%, 1.68% and 15.96%. R's total positive error reduction was 84.71% of M1's, while its total adverse increment was 85.25% lower. Positive and adverse readout sets are defined separately for each model. This decomposition explains the net benefit on the development panel.

Readout contributions delimit this benefit. Ten readouts accounted for 92.13% of R's positive error reductions, with negative contributions retained in the primary loss. Its gain was −0.118% under training-variance-standardized scoring and −0.380% after deleting the previously locked ten-readout set (Fig. 5e,f). The raw-scale increment is therefore concentrated in a small response subset; neither independent biological programs nor their effective number were established.

Outside the fixed set, gain improved from −10.03% for $M_1$ to −0.38% for R, showing reduced damage in weakly supported readouts. The set includes S100A9, S100A8, S100A12 and LYZ; archived four-corner responses illustrate the corrected variation. Uniform shrinkage, 200 expression/variance-stratified weight permutations and matched-norm random increments provide complementary controls (Supplementary Notes 6 and 13).

### 15.5 Quantitative outcomes of the additional research uses

Two separate retrospective analyses examine further uses (Supplementary Fig. S1). Risk selection retained approximately half the evaluation weight with accepted/all-query MSE ratio 0.5343 on 213 records from 33 calibration-disjoint pair identities. PairPick, using R-G, selected the empirical best in 19 of 28 fixed-anchor tasks, with mean regret 0.006197. Its achieved utility tied linear prediction in 27 tasks, and its choices matched PRESAGE in all 28. The original design table directly matched 15 of 49 detected gene-mode labels, with one enhancer-only correspondence and 33 unlisted pairs.[23] These results illustrate an executable recommendation operation and its limits: response differences often leave choices unchanged, and detected-label evaluation does not confirm intended cell-level constructs. Full outcomes, candidate sizes and identity status are retained in Supplementary Tables S17–S19, S28 and S33.

## Supplementary Note 16 | Partner corrections across evaluation settings

### 16.1 Design, comparator and evaluation status

Table S36 places the principal partner-increment results beside the experimental support used to obtain them. A panel's total number of pairs is different from the number available to a particular fit. Training counts below describe the final outer fit; inner selection can use smaller subsets. The table groups transcriptomic and functional tasks without treating them as points on a common learning curve. Each gain compares a specified correction with its own baseline under the original loss and aggregation rule.

### Supplementary Table S36 | Experimental support and observed partner increments

**A. Permitted training support and held-out measurements.** Counts refer to unique pair identities, not cells, technical half-samples or repeated optimization runs.

| Setting | Eligible pair pool | Endpoints | Training pairs per fit | Query pairs per fold | Readout |
|---|---:|---:|---|---|---|
| CaRPool RNA | 55 | 14 | 43 / 44 / 44 / 44 / 45 | 12 / 11 / 11 / 11 / 10 | 4,995 non-target RNA readouts |
| Dixit RNA | 45 | 10 | 30 / 30 / 30 | 15 / 15 / 15 | 23,519 non-target RNA readouts |
| Combi-seq | 374 | 36 drugs | 249 / 249 / 250 | 125 / 125 / 124 | 51 program genes |
| K562 low | 2,528 | 80 | 421 / 421 / 422 | 211 / 211 / 210 | Growth; two repeats |
| K562 high | 2,528 | 80 | 2,317 / 2,317 / 2,318 | 211 / 211 / 210 | Growth; two repeats |
| Jurkat | 3,081 | 79 | 2,464 | 617 | Growth; two repeats |
| Shen HeLa | 2,619 | 73 | 2,096 | 523 | Growth slope per day; two repeats |
| Shen A549 | 2,625 | 73 | 2,100 | 525 | Growth slope per day; two repeats |
| Shen 293T | 2,573 | 73 | 2,056 | 517 | Growth slope per day; two repeats |

K562's eligible pool is the 3,160 pairs among eighty endpoints after exclusion of the original registered 632-pair test set. Within each revised development fold, the low training pool is nested in the high pool and the queries are unchanged. Low-budget fits leave some eligible measurements unused. The 632 revised query identities are distinct from the original registered test identities. Jurkat uses approximately 80% of its pair universe for fitting; this fraction is the evaluated setting, not an estimated minimum requirement. Fold-specific counts and the identities of measurements are retained in the original sources.

**B. Added prediction, primary outcome and evidence status.** Positive gain means lower error than the stated reference. Complete-response gains, alternative corrections and secondary program analyses remain in their original tables.

| Setting | Correction and reference | Partner gain | Evaluation status |
|---|---|---|---|
| CaRPool RNA | Readout-shrunk product-kernel R versus M0 | +15.9622% raw; −0.1180% variance-standardized; −0.3798% after ten-readout deletion | Development panel; source shared with Wessels benchmark |
| Dixit RNA | Frozen selected RNA correction versus refitted M0 | −0.0326% | Frozen procedure; independent source |
| Combi-seq | Fixed-reference program correction B+E+gF versus B+E | −0.0302% | Frozen procedure; independent drug source |
| K562 low | Primary rank-one correction versus same-budget M0 | −0.0326% | Revised development |
| K562 high | Primary rank-one correction versus same-budget M0 | +11.6838% | Revised development |
| Jurkat | Primary rank-one correction versus refitted M0 | +6.1750% | Frozen procedure; new cell context within the same study |
| Shen HeLa | Original primary rank-one correction versus refitted M0 | −1.0995% | Original frozen procedure; independent study |
| Shen A549 | Original primary rank-one correction versus refitted M0 | −0.1092% | Same independent study; shared-design background |
| Shen 293T | Original primary rank-one correction versus refitted M0 | −0.5265% | Same independent study; shared-design background |

CaRPool's primary endpoint averages squared error over eleven same-fit four-corner contrasts and 4,995 readouts. Its standardized sensitivity uses training-condition variance with the declared 0.01 positive-quantile floor; the deletion analysis retains the previously locked ten-readout set. The other rows use endpoint-eliminated losses on their full query graphs. For growth, squared errors are pooled across the specified folds and measured repeats before taking the ratio. These losses remove endpoint contributions, but their designs, weighting, units and predictive models differ. The table therefore preserves the raw-scale CaRPool benefit alongside its sensitivity results rather than selecting one scoring convention to imply a sample-size rule.

### 16.2 What the support comparison establishes

The revised K562 comparison shows that expanding the measured combinations can make the rank-one correction useful on unchanged queries. Both the baseline and correction are refitted at each budget, so the incremental gain is measured against a baseline receiving the same additional data. The low and high settings jointly change condition count, endpoint degree and layout coverage. Only two budgets were evaluated for this revised procedure; they neither locate a transition nor establish monotonic improvement. The original registered five-budget experiment selected zero correction throughout and remains separately reported in Supplementary Note 11 and Table S23.

Across sources, support quantity also varies with intervention type, readout structure, measurement noise and fitting procedure. A high-dimensional RNA vector can share a low-dimensional response pattern, while a scalar growth phenotype can depend on many partner effects. Output dimension alone therefore supplies no quantitative sample requirement. The table documents where the tested corrections helped and where an increment was not detected; it does not establish impossibility for other estimators or a universal threshold for partner learning. Jurkat tests transfer of a frozen procedure with local re-estimation; Dixit and Combi-seq provide independent-source checks on their distinct targets. None of these rows changes its original development or evaluation status.

Machine-readable values and source hashes are supplied in `revision/PARTNER_SCOPE_V30.csv` and `revision/PARTNER_SCOPE_AUDIT_V30.json`. The original overview is reconstructed from existing results. The added Shen rows reproduce the original new-study evaluation; subsequent component control is separately reported in Tables S37–S40 and does not replace these rows.

## Supplementary Note 17 | Component-specific control: derivation and objectives

### 17.1 Preserve the original endpoint contribution

In this section α denotes the partner-control coefficient, distinct from the saturation parameter in the core predictor.

For a fixed query panel let $H=XX^+$ and $P=I-H$. Decompose an added correction as $D_E=HD$ and $D_P=PD$. The controlled prediction is

$$\widehat Y_\alpha=M_0+D_E+\alpha D_P.$$

Because $HP=PH=0$, $H(\widehat Y_\alpha-M_0)=HD$ and $P(\widehat Y_\alpha-M_0)=\alpha PD$. The baseline can already have partner variation. The operation preserves the endpoint contribution of the added correction, not every component of the baseline error.

Let $e=Y-M_0$. Orthogonality gives

$$L_{\rm full}(\alpha)=\|He-HD\|_F^2+\|Pe-\alpha PD\|_F^2.$$

Subtracting this expression at $\alpha$ from its value at one proves

$$L_{\rm full}(1)-L_{\rm full}(\alpha)=L_{\rm partner}(1)-L_{\rm partner}(\alpha).$$

These are unweighted SSE changes relative to the original correction. Complete and partner percentage gains use different denominators. The identity is exact for the operation, irrespective of the outcomes; improvement for a training-selected coefficient remains empirical. Whole attenuation $M_0+\alpha D$ has the same partner prediction but a different endpoint part. Component control need not outperform whole attenuation when shrinking the endpoint part is itself helpful.

### 17.2 Distinct optima for complete and partner prediction

Write $e_E=He$ and $e_P=Pe$. The unconstrained common-amplitude optimum is

$$t^*_{\rm full}=\frac{\langle e_E,D_E\rangle+\langle e_P,D_P\rangle}{\|D_E\|_F^2+\|D_P\|_F^2}.$$

When both energies are positive, it is the correction-energy-weighted mean of the endpoint and partner optima. Consequently, complete-response selection can retain a correction whose partner contribution is excessive. The relation applies to unconstrained optima, not an average of separately clipped optima. We retain endpoint amplitude one and fit only the partner coefficient, rather than introduce another model search.

For each inner held-combination graph, let $r_k=P_k(Y_k-M_{0,k})$ and $d_k=P_kD_k$. Set $A=\sum_k\langle r_k,d_k\rangle$, $B=\sum_k\|d_k\|_F^2$ and $E=\sum_k\|r_k\|_F^2$. The internal loss is $E-2\alpha A+\alpha^2B$. Minimizing it over $[0,1]$ gives $\alpha=\operatorname{clip}_{[0,1]}(A/B)$; implementation uses zero if $B\le10^{-14}$. Original zero-selected models remain zero corrections.

For positive $B,E$, $\rho=A/\sqrt{BE}$ describes alignment on these internal arrays. The unconstrained scalar improvement along this fixed direction is $100\rho^2$ percent of E. When $A\ne0$, defining $k=\alpha/(A/B)$ gives internal gain $100\rho^2(2k-k^2)$. These are algebraic descriptions, not independent gain forecasts. Shen's current internal directions admit only approximately 0.035%, 0.067% and 0.043% unconstrained scalar improvement in HeLa, A549 and 293T. Attenuation chiefly avoids magnifying a weak direction; it does not learn a better direction.

### 17.3 Metric and panel scope

The proof uses one unweighted squared-error axis. Common readout scaling preserves the argument after the same scaling is applied everywhere. Row-weighted loss requires the corresponding weighted projection; arbitrary element weights do not inherit the unweighted decomposition. Rectangle scores retain their original weighting of overlapping contrasts.

The projector uses query identities, so the operation is panel-dependent. A change in panel membership may change a pair's prediction even without query-label access. A graph with no partner directions retains the original predictor. The full derivation, including weighted notation and readout-specific extensions, is supplied in `partner_control/METHOD_DERIVATION.md`.

## Supplementary Note 18 | Complete functional comparison and source relationships

### 18.1 Original Shen evaluation and subsequent control

Shen et al. studied combinatorial CRISPR–Cas9 perturbations in HeLa, A549 and 293T.[20] The accompanying source archive recovers raw count columns from the authors' official repository, qualifies physical constructs and constructs repeated log2 relative-count growth slopes per day. This unit differs from Horlbeck's population-doubling-normalized phenotype. Training/query pairs after qualification are 2,096/523, 2,100/525 and 2,056/517, respectively, with 73 endpoints in each background. The same physical identity split precedes qualification. Whole-map published interaction scores are not inputs.

The original signed rank-one procedure retains complete-response inner selection and the zero option. It improved complete prediction in all three backgrounds but increased partner error in each. Its equal-background partner gain is −0.5784%; the original conditional gene-weight sensitivity interval includes zero. Original predictions, all source-repeat outcomes and secondary spectral/additive comparisons are retained in `reproduction/WitnessCell_Shen2017_Formal_20261108.zip`.

The control rule was developed after these outcomes and historical Horlbeck analyses had been inspected. It uses training arrays to calculate coefficients but is retrospective development. The scalar decisions are written before the new replay reads external outcomes; that execution order does not retrospectively establish blinded rule development or preregistration.

### 18.2 All eligible functional configurations and four outputs

Eligibility requires matching archived selected-model internal predictions and external predictions. All 55 eligible configurations are included: six revised K562 budget folds, one Jurkat fit, three Shen backgrounds and 45 K562 layout configurations. Twenty-eight selected zero correction; twenty-seven selected a nonzero correction. Every case receives four outputs: baseline M0, original M0+D, whole attenuation M0+αD, and component control M0+HD+αPD. No outer score selects an amplitude or an alternative predictor.

### Supplementary Table S37 | Complete functional four-policy comparison

All gains are percentages relative to the corresponding M0. Baseline gains are zero by definition and are included explicitly in the 220-row source table. Each row below pools the same-condition squared errors before taking a ratio. K562 layout pooling is a descriptive configuration summary over repeated queries, not 45 independent experiments.

| Condition | Policy | Complete gain (%) | Partner gain (%) |
|---|---|---|---|
| K562 low | Original | 0.265197 | -0.032639 |
| K562 low | Whole attenuation | 0.265197 | -0.032639 |
| K562 low | Partner only | 0.265197 | -0.032639 |
| K562 high | Original | 11.029218 | 11.683843 |
| K562 high | Whole attenuation | 10.955798 | 11.592353 |
| K562 high | Partner only | 10.973593 | 11.592353 |
| Jurkat | Original | 6.430835 | 6.175031 |
| Jurkat | Whole attenuation | 6.430835 | 6.175031 |
| Jurkat | Partner only | 6.430835 | 6.175031 |
| K562 layouts | Original | 0.056248 | -0.221688 |
| K562 layouts | Whole attenuation | 0.293178 | 0.192359 |
| K562 layouts | Partner only | 0.278864 | 0.192359 |
| HeLa | Original | 4.161291 | -1.099485 |
| HeLa | Whole attenuation | 0.779209 | -0.037207 |
| HeLa | Partner only | 4.918274 | -0.037207 |
| A549 | Original | 3.649825 | -0.109180 |
| A549 | Whole attenuation | 1.154472 | 0.116419 |
| A549 | Partner only | 3.816872 | 0.116419 |
| 293T | Original | 6.182256 | -0.526543 |
| 293T | Whole attenuation | 1.255438 | 0.060085 |
| 293T | Partner only | 6.595579 | 0.060085 |

`partner_control/results/FOUR_ARM_SCORES.csv` retains every configuration, all four policies and complete/endpoint/partner SSE. Relative to the original predictor, component control improves fifteen configurations, ties thirty-three and worsens seven. These counts use absolute SSE tolerance $10^{-12}$ and are descriptive. Jurkat remains unchanged; the K562 high-support pooled partner gain changes from 11.683843% to 11.592353%. The low-support correction remains unhelpful. Thus the control rule preserves endpoint contributions exactly but does not guarantee partner improvement.

Independent replay of the internal coefficient agrees with the historical calculation to numerical precision. The observed maximum endpoint-preservation error is below 1.1×$10^{-16}$; complete and partner SSE changes agree within 4.5×$10^{-16}$. These are computational checks of the operation, not evidence of external predictive reliability.

## Supplementary Note 19 | Repeat consistency and fixed-query-panel dependence

### 19.1 Repeated measurements of the calibration signal

The primary coefficient pools both internal repeats. Diagnostic coefficients use each repeat separately, with the original fitted biological models held fixed, and are evaluated on both outer repeats. A repeat-specific coefficient does not imply a model trained exclusively on that repeat. All 110 coefficients and 990 policy-by-evaluation-repeat scores are supplied in `REPEAT_DECISIONS.csv` and `REPEAT_SCORES.csv`.

### Supplementary Table S38 | Coefficients estimated from pooled and individual training repeats

| Background | Pooled α | Repeat 1 α | Repeat 2 α |
|---|---|---|---|
| Jurkat | 1.000000 | 1.000000 | 1.000000 |
| HeLa | 0.118745 | 0.000000 | 0.286283 |
| A549 | 0.239196 | 0.291960 | 0.186432 |
| 293T | 0.142838 | 0.144390 | 0.141286 |

Among the twenty-seven nonzero corrections, the median absolute difference between repeat-specific coefficients is 0.155482 and the maximum is one. Including zero-selected models would reduce the median to zero, so the active-model summary is reported separately. HeLa's repeat-one alignment is negative and selects zero, whereas repeat two selects approximately 0.286. The pooled coefficient is approximately 0.119. This inconsistency limits the interpretation of the pooled training signal.

### Supplementary Table S39 | Outer-repeat partner gains using the pooled coefficient

| Condition | Measured repeat | Original gain (%) | Controlled gain (%) |
|---|---|---|---|
| K562 low | 1 | -0.012158 | -0.012158 |
| K562 low | 2 | -0.068831 | -0.068831 |
| K562 high | 1 | 9.305122 | 9.219758 |
| K562 high | 2 | 15.862689 | 15.760437 |
| Jurkat | 1 | 3.814777 | 3.814777 |
| Jurkat | 2 | 9.041476 | 9.041476 |
| K562 layouts | 1 | -0.040058 | 0.216406 |
| K562 layouts | 2 | -0.536966 | 0.150619 |
| HeLa | 1 | -0.840157 | 0.006489 |
| HeLa | 2 | -1.302654 | -0.071440 |
| A549 | 1 | -0.909343 | -0.073242 |
| A549 | 2 | 0.671975 | 0.301575 |
| 293T | 1 | -1.368366 | -0.048565 |
| 293T | 2 | 0.192114 | 0.152839 |

None of the three Shen backgrounds has positive controlled partner gain in both repeats. Their small pooled changes therefore should not be described as stable external partner recovery. Jurkat retains positive gains in both repeats; K562 high support also remains positive in both. The full repeat-specific coefficient-by-outcome matrix is retained without choosing the favorable training repeat.

### 19.2 Outcome-free changes to the query panel

Canonical unordered endpoint indices sort the pairs. Alternating ranks define two and three subpanels. The coefficient, M0 and D stay fixed; only the projection panel changes. The subpanel predictions are assembled and scored on the original full query graph with its original projector, maintaining a common evaluation target. Singleton control retains each original prediction. The total number of adjustable directions is the sum of the subpanel contrast dimensions.

### Supplementary Table S40 | Partner gain under fixed identity-based query partitions

| Condition | Full panel (%) | Halves (%) | Thirds (%) | Singletons (%) |
|---|---|---|---|---|
| K562 low | -0.032639 | -0.032639 | -0.032639 | -0.032639 |
| K562 high | 11.592353 | 11.589502 | 11.664539 | 11.683843 |
| Jurkat | 6.175031 | 6.175031 | 6.175031 | 6.175031 |
| K562 layouts | 0.192359 | -0.071830 | -0.186217 | -0.221688 |
| HeLa | -0.037207 | -0.118462 | -0.769573 | -1.099485 |
| A549 | 0.116419 | -0.243117 | 0.032759 | -0.109180 |
| 293T | 0.060085 | -0.207309 | -0.286356 | -0.526543 |

Some signs change after panel splitting, notably A549 and 293T. Jurkat is invariant because its coefficient equals one. These results establish a usage requirement: specify the query panel when applying control. Reordering the same panel or reversing pair orientation leaves the operation equivariant. Splitting the panel changes the available directions and is not required to preserve predictions. The 220-row panel table includes complete losses, contrast dimensions and prediction changes for all configurations. Supplementary Fig. S4 displays the full analysis without treating shared-query configurations as independent repeats.

## Supplementary Note 20 | RNA candidate replay, historical diagnostics and software

### 20.1 Matched-candidate RNA comparison

All 265 eligible saved RNA candidate records with three matching internal-fold prediction arrays and corresponding external predictions were processed. They come from 35 archived fit/snapshot pools spanning the same five CaRPool outer folds. Repeated snapshots, fallback branches and candidate variants are not independent experiments. We recomputed each candidate's internal cross-product and correction energy, fitted the same bounded coefficient, and evaluated all four policies on the full query projection. This produces 1,060 score rows, separate from the primary eleven-rectangle R comparison in Tables S10–S13.

### Supplementary Table S41 | RNA replay scope and selection identity

| Item | Result |
|---|---|
| Eligible candidate records / outer folds | 265 / 5 |
| Archived fit/snapshot pools | 35 |
| Four-policy score rows | 1,060 |
| Nontrivially attenuated candidate records | 45 |
| Partner SSE improves / ties / worsens | 42 / 220 / 3 |
| Primary rectangle result | Unchanged; separate metric and prediction definition |
| Stage-one R_12 workflow | All five folds select native branch, η = 0 |
| Stage-one R_neural workflow | Four native; one neural blend (fold 2, η = 0.25) |

The counts describe fixed candidates and use SSE tolerance $10^{-10}$. They are not selection-corrected success rates. Later stage-two snapshots also remain in the complete source table and have different selected branches. `RNA_SUFFICIENT_STATISTICS.csv` contains the exact inner and outer endpoint/partner baseline energies, cross-products and correction energies, with source hashes. `replay_rna.py` reconstructs all policies without accessing the historical source directories. Inner and original external gain reconstruction agree with the archive to less than 3×$10^{-14}$ percentage points. These diagnostics do not replace the original readout shrinkage procedure or establish independent RNA transfer.

### 20.2 Diagnostic variables and historical regression

Original projected-response agreement, q_raw, characterizes the stated measurement and preprocessing. Agreement between projected residuals after a baseline, q_residual,A, also depends on the fitted model A. Internal held-combination alignment describes whether the current correction points toward the remaining response, while its energy measures the cost of its magnitude. These quantities remain distinct: repeated measurements can agree without being predictable from the training partners by a particular model.

The retained historical ledger contains 788 matched records from 64 training configurations, eight panels and five studies, including shared queries, model candidates and archive aliases. Its data-plus-internal-score regression has study-held-out R² of approximately −0.084. A rank-one functional subset shows a positive study-held-out relationship but contains only two studies. These analyses motivate diagnostics and expose their transfer limitations; they are not a universal gain law or a calibrated new-study forecast. Original ledgers and regression summaries are included under `partner_control/inputs/ARCHIVED_*` and the accompanying calibration source archive.

### Supplementary Table S42 | Evidence roles in the control extension

| Evidence | Role | Evaluation status |
|---|---|---|
| Endpoint preservation and SSE identity | Exact behavior on a fixed panel and common loss | Algebraic result |
| Original Jurkat result | Positive frozen-procedure functional increment | Same-study refitting in another cell environment |
| Original Shen result | Complete gain with adverse partner increment | New-study evaluation of the prior procedure |
| Functional four-policy comparison | Effect of retrospective scalar/component control | 55 saved configurations; historical outcomes visible |
| Repeat and panel checks | Robustness of this operation | Fixed predictions; no new biological experiments |
| RNA candidate replay | Generality of the arithmetic across saved candidates | Same development panel; separate global projection |
| Historical regression | Diagnostic association and transfer boundary | Shared configurations and aliases; retrospective |
| Han full-screen calibration | Test target-estimated strength with frozen shared functions | Prior Han access; locked full-screen reassessment; all six portable configurations retained |

### 20.3 Executable input and outcome separation

### Supplementary Table S43 | Partner-control interface

| Entry | Allowed input | Output |
|---|---|---|
| fit_amplitude / fit | Internal observations, held-combination baseline and predictions, pairs, folds | Bounded coefficient and internal diagnostics |
| predict_panel / predict | Query baseline and predictions, query pair identities, coefficient | Original, whole-attenuated and component-controlled predictions |
| run_analysis.py | Saved inner and outer arrays; scoring after decisions | All functional comparisons and two checks |
| replay_rna.py | Archived SSE sufficient statistics | All supplementary RNA comparison scores |
| test_partner_control.py | Synthetic behavior cases and real archived pair splits | Interface, preservation, permutation and label-rejection checks |

The fit interface expects the caller to supply genuine internal held-combination predictions. It cannot verify model-training provenance from arrays. The prediction command rejects extra fields, including query outcomes. A real HeLa example fits the coefficient anew and applies it without reading its separate evaluation-truth file. New-directory replay gives identical predictions. This is a runnable adjustment interface, distinct from refitting the underlying biological predictor. The parent package retains the original RNA and functional fitting interfaces.

## Supplementary Note 21 | Target-side calibration of frozen partner functions in Han

### 21.1 Analysis identity, source models and outcome access

The reconstructed panel comprised 20,359 qualified pairs among 206 active endpoints, with two reconstructed repeats. The published panel comprised 20,791 pairs among 205 active endpoints. These panels differ in response construction, permitted single-response inputs and coverage; they do not isolate the effect of response definition alone. Raw-panel eligibility followed the fixed count, plasmid, construct-count and orientation-support rules, without filtering by model performance. ADRBK1 lacked a stored embedding and used the existing global-mean rule; complete endpoint exclusion was examined separately.

All three source backgrounds and both portable model families were retained in the full-screen analysis. The full-screen comparison is the shared-function extension of the archived protocol, alongside the protocol-designated common-identity rank-one primary analysis. Table S49 states the original primary analysis and its outcome alongside the full-screen extension; complete records remain in the source archive. Source configurations are model variants applied to one target study, not independent target-study replications.

### Supplementary Table S44 | Full-screen design and model coverage

| Item | Reconstructed panel | Published panel |
|---|---|---|
| Eligible unordered pairs | 20,359 | 20,791 |
| Active endpoints | 206 | 205 |
| Response measurements | Two end-sample repeats | One published aggregate phenotype |
| Baseline | Corresponding Han single-response addition | Corresponding Han single-response addition |
| Frozen partner functions | SF12 and V4 from HeLa, A549 and 293T | Same six configurations |
| Target-estimated quantities | One clipped scalar per configuration and fold | Same rule |
| Calibration / query pairs per fold | 13,572–13,573 / 6,786–6,787 | 13,860–13,861 / 6,930–6,931 |
| Pair splits | Five partitions, three outer folds, two calibration blocks | Same design |
| Base-predictor refitting | None | None |
| Query access | Predictions locked before query scoring | Predictions locked before query scoring |
| Evaluation status | Frozen reassessment after prior Han access | Related response definition from the same study |

Raw read reconstruction used SRR5258688 (plasmid), SRR5258687 (initial sample), SRR5258686 and SRR5258685 (end repeats).[43] The eight compressed read files total 21,533,213,518 bytes; their source checks and compact counts are preserved in the supplied archive. Eligible constructs require an initial count of at least 50 and a positive plasmid count. For construct c and end repeat r, enrichment is log2[max(n_c,r,1)/n_c,plasmid]. Each repeat is centred by the median eligible safe-control enrichment and scaled by its sample standard deviation. Single responses are medians over gene/safe-control constructs in both orientations. Double responses require at least three eligible ordered constructs and both orientations. Among 21,321 possible pairs from 207 finite single inputs, these fixed rules retain 20,359 pairs; the reasons for all exclusions are supplied in `audit/ALL_RAW_PAIR_QUALIFICATION.csv` inside the archive. The reconstructed response is a declared control-standardized reanalysis rather than an exact reconstruction of published pZ values.

### 21.2 Same-direction control comparisons

### Supplementary Table S45 | Full-screen target-side calibration

 Positive percentages indicate error reduction. “Unchanged” uses the same frozen direction with α = 1. “Removed” uses α = 0 and is identical to the additive baseline in this adapter.

| Readout | Source and model | Partner reduction vs unchanged (%) | Partner reduction vs removed (%) | Complete reduction vs additive (%) |
|---|---|---:|---:|---:|
| Reconstructed | HeLa SF12 | 12.839411 | 0.261676 | 0.186996 |
| Reconstructed | 293T SF12 | 0.054207 | 0.651347 | 0.465458 |
| Reconstructed | A549 SF12 | 0.000000 | 0.000000 | 0.000000 |
| Reconstructed | HeLa V4 | 7.500390 | −0.003666 | −0.002620 |
| Reconstructed | 293T V4 | 1.107888 | −0.003626 | −0.002591 |
| Reconstructed | A549 V4 | −0.001361 | 0.003431 | 0.002452 |
| Published | HeLa SF12 | 14.495292 | 0.143619 | 0.094699 |
| Published | 293T SF12 | 0.288967 | 0.305757 | 0.201609 |
| Published | A549 SF12 | 0.000000 | 0.000000 | 0.000000 |
| Published | HeLa V4 | 6.276387 | −0.008098 | −0.005339 |
| Published | 293T V4 | 0.926410 | −0.001970 | −0.001299 |
| Published | A549 V4 | −0.003073 | −0.003932 | −0.002592 |

### 21.3 Partition and repeat results

For reconstructed outcomes, the HeLa SF12 partner gain over the additive baseline ranged from 0.176070% to 0.351213% across the five partitions; the corresponding range for 293T SF12 was 0.592577% to 0.691432%. Gains over the unchanged correction were also positive in every partition. Estimated strengths across folds ranged from 0.078561 to 0.166958 for HeLa and from 0.582805 to 0.868120 for 293T. The five overlapping partitions describe sensitivity to pair allocation and are not independent biological replicates.

### Supplementary Table S46 | Reconstructed-repeat gains over the additive baseline

 Each entry pools the locked predictions across the five partitions.

| Source and model | Repeat | Partner reduction (%) | Complete reduction (%) |
|---|---|---:|---:|
| HeLa SF12 | 1 | 0.442913 | 0.318944 |
| HeLa SF12 | 2 | 0.080751 | 0.057269 |
| 293T SF12 | 1 | 0.741102 | 0.533672 |
| 293T SF12 | 2 | 0.561746 | 0.398393 |

The two reconstructed repeats share initial and plasmid controls. The published phenotype derives from the same study. These results assess repeat and processing sensitivity within Han rather than adding independent source confirmations. No population-level confidence claim is inferred from the five overlapping partitions.

### Supplementary Table S47 | Pair-partition sensitivity for every full-screen configuration

Ranges cover the five partition-level partner-error reductions versus removing the correction. The last column counts positive partitions; it is descriptive, not an independent-binomial significance test.

| Readout | Source and model | Minimum (%) | Maximum (%) | Positive partitions / 5 |
|---|---|---:|---:|---:|
| Reconstructed | 293T SF12 | 0.592577 | 0.691432 | 5 |
| Reconstructed | 293T V4 | -0.011665 | 0.001466 | 1 |
| Reconstructed | A549 SF12 | 0.000000 | 0.000000 | 0 |
| Reconstructed | A549 V4 | 0.000261 | 0.005208 | 5 |
| Reconstructed | HeLa SF12 | 0.176070 | 0.351213 | 5 |
| Reconstructed | HeLa V4 | -0.012750 | 0.000000 | 0 |
| Published | 293T SF12 | 0.254765 | 0.338189 | 5 |
| Published | 293T V4 | -0.007876 | 0.003819 | 2 |
| Published | A549 SF12 | 0.000000 | 0.000000 | 0 |
| Published | A549 V4 | -0.006640 | -0.000208 | 0 |
| Published | HeLa SF12 | 0.031681 | 0.226332 | 5 |
| Published | HeLa V4 | -0.018471 | 0.004286 | 1 |

### 21.4 Scope of the control comparison

The full-screen adapter fixes the endpoint contribution to that of B and supplies only $P_QD_Q$. Consequently, changing α preserves that contribution and changes complete and partner squared-error sums by the same absolute amount. Their percentage changes use different denominators and are not independent successes. This experiment tests whether the target-estimated strength is useful; the distinction between component control and shrinking a correction that contains both components is evaluated by the original Shen analyses.

### 21.5 Access history and identity checks

Earlier Han aggregate and reconstructed outcomes had been accessed and scored before this extension. The new protocol, source files and run code were locked before its reassessment scores, and each query prediction was locked before scoring, but this does not reset the prior access history. Excluding previously scored pairs from both calibration and query pools retained positive SF12 gains: reconstructed partner reductions were 0.268364% for HeLa and 0.659370% for 293T. Excluding the missing-embedding endpoint also retained the gains. These are post-score, identity-defined sensitivity checks under the same algorithm, not new confirmatory experiments.

### Supplementary Table S48 | Identity-defined sensitivity of the SF12 increments

The exclusions were applied to both calibration and query pools after the main results were visible. Models, coefficient formula and partition seeds were unchanged. Full results for all models accompany the source table.

| Readout | Exclusion | Source | Partner reduction (%) | Complete reduction (%) |
|---|---|---|---:|---:|
| Published | Missing-embedding endpoint | 293T | 0.306129 | 0.202036 |
| Published | Missing-embedding endpoint | HeLa | 0.143527 | 0.094724 |
| Published | Prior scored pairs | 293T | 0.312258 | 0.205762 |
| Published | Prior scored pairs | HeLa | 0.151611 | 0.099904 |
| Reconstructed | Missing-embedding endpoint | 293T | 0.652638 | 0.465389 |
| Reconstructed | Missing-embedding endpoint | HeLa | 0.263656 | 0.188010 |
| Reconstructed | Prior scored pairs | 293T | 0.659370 | 0.471322 |
| Reconstructed | Prior scored pairs | HeLa | 0.268364 | 0.191829 |


### 21.6 Reproduction and numerical checks

`reproduction/WitnessCell_Han2017_Frozen_Reassessment_v3.zip` contains compact counts, the published source table, permitted inputs, frozen shared functions and source feature assets, split identities, per-fold coefficients and locked predictions. Its recorded checks include unchanged source-model hashes, exact replay of the original shared-identity predictions, independent QR rescoring of all 60 fold jobs (maximum absolute SSE difference $1.36 \times 10^{-12}$), and query-outcome substitution with unchanged coefficients and predictions. These are implementation checks, not biological replications. The current manuscript verification also rebuilds the reported relative reductions from saved fold-level SSEs.

From the manuscript package root:

```bash
python code/replay_han.py --workdir /absolute/path/new_han_replay
```

Add `--rebuild-measurements` to reconstruct the input arrays from the retained compact counts and published table. This replays frozen inference and fits the declared scalar coefficients; no base model is trained. The complete source archive also preserves the common-identity primary analysis (Table S49). `--mode exact` requires exact numerical equality. Default compatible replay checks immutable assets exactly, then compares inference, coefficients and scores with precision-aware bounds. It records the earliest observed numerical difference rather than attributing that difference to an environment version without an isolation experiment. The full acceptance rules and measured replay report are in `reproduction/NUMERICAL_REPLAY.md` and `revision/replay_v36/`.


### 21.7 Original primary analysis and full-screen extension

The reassessment protocol distinguishes its common-identity rank-one primary from the shared-function extension. The primary pools the three source rank-one models across the fixed partitions for each readout. The full-screen analysis changes the input adapter: target single-response addition supplies the baseline, and only a projected shared-function correction is added. It therefore addresses a different, narrower use of the frozen predictors. Current emphasis on the full-screen results does not change the original primary designation.

### Supplementary Table S49 | Analysis stages and their evidential roles

Percentages below are pooled partner-error reductions, positive for improvement. The first row is reconstructed from `revision/han_calibration/PRIMARY_COMMON_RANK1.csv`; the protocol is retained alongside it.

| Analysis | Designated role and question | Outcome and current presentation |
|---|---|---|
| Common-identity rank-one reassessment | Protocol primary: can target-estimated strength improve the frozen original-model correction on shared identities? All three source rank-one models retained. | Relative to unchanged correction: +6.6700% reconstructed, +6.3133% published. Relative to removing the partner part: −0.9961% reconstructed, +0.0350% published. The reconstructed primary supports attenuation but not a net gain over removal. Complete source records remain available. |
| Full-screen SF12/V4 adaptation | Protocol-defined extension: can target calibration use a frozen shared direction on eligible new endpoint identities, with a target additive baseline? | Two SF12 configurations improve over unchanged and removed corrections in both readouts; all six configurations appear in Table S45 and Supplementary Fig. S5. This is the current main-text Han result. |
| Identity exclusions | Post-score sensitivity: exclude previously scored pairs or the missing-embedding endpoint from calibration and query pools. | Table S48 retains both readouts. These checks assess the same extension, with no new independent study or revised primary designation. |

The reassessment followed prior access to Han. Its local protocol and per-fold prediction locks document decision order within that reassessment; they do not turn the study into an untouched prospective confirmation. Earlier fixed-source transfer analyses remain historical source records, separate from the target-calibrated reassessment summarized here.

## Supplementary Note 22 | Shared functional predictors: features, source fitting and portable inference

### 22.1 Correspondence to the functional and control interfaces

The indexed rank-one and spectral models predict entries of a fitted gene-by-gene residual matrix. SF12 and V4 evaluate shared functions of single responses and gene annotations. The latter can be evaluated on new identities with suitable annotations without assigning those identities new fitted matrix parameters. The common interface is a baseline $M_0$, a correction $D$ and permitted calibration predictions. The functions have different parameters and source fitting rules; none inherits the released core's benchmark rank by sharing this interface.

The full-screen adapter uses the target additive baseline $B$ and supplies $P_QD_Q$. Source functions are specified before query evaluation. The target coefficient adjusts each given function; the experiment supplies no validated rule for selecting a best source from the six configurations using query results. The original fitting scripts, complete source selections and checkpoints are included in `reproduction/WitnessCell_Shen_Jurkat_Partner_Repair_20261108.zip`. Their retained parameters match the Han frozen-engine copies; `revision/SOURCE_MODEL_CORRESPONDENCE_V36.json` records that correspondence.

### 22.2 Exact SF12 features and source fitting

Let $s_g$ be the source's repeat-mean single effect and $r_s$ its root-mean-square across endpoints. Set $t_g=\tanh(s_g/r_s)$, $t_g^+=\max(t_g,0)$ and $t_g^-=\max(-t_g,0)$. The four local terms are

$$
\ell_{1,gh}=r_st_g^+t_h^+,\qquad
\ell_{2,gh}=r_st_g^-t_h^-,\qquad
\ell_{3,gh}=r_s(t_g^+t_h^-+t_g^-t_h^+),\qquad
\ell_{4,gh}=\ell_{3,gh}\tanh(t_g+t_h).
$$

All terms are invariant to exchanging $g$ and $h$. The stored 64-dimensional GenePT annotation matrix $E_s$ is centred over source endpoints. Its normalized Gram matrix is eigendecomposed; the two leading eigenvectors are scaled by square-root eigenvalues, oriented by the largest absolute entry in each column and jointly normalized to unit root-mean-square. Denote these coordinates by $z_{g1},z_{g2}$. For $j=1,2$, the source training pairs determine the mean $c_j$ and standard deviation $v_j$ of $z_{gj}z_{hj}$. Set $v_j\geq10^{-8}$ and $m_{j,gh}=(z_{gj}z_{hj}-c_j)/v_j$, with $m_{0,gh}=1$.

### Supplementary Table S50 | The twelve SF12 features

Features are ordered by local term first, then modifier; each has its own fitted coefficient. The quantities $r_s,z,c,v$ are fixed by the permitted source inputs and source training pairs.

| Coefficient indices | Local response term | Constant feature | First annotation-product feature | Second annotation-product feature |
|---|---|---|---|---|
| 1–3 | Positive–positive | $\ell_1$ | $\ell_1m_1$ | $\ell_1m_2$ |
| 4–6 | Negative–negative | $\ell_2$ | $\ell_2m_1$ | $\ell_2m_2$ |
| 7–9 | Opposite signs | $\ell_3$ | $\ell_3m_1$ | $\ell_3m_2$ |
| 10–12 | Signed opposite signs | $\ell_4$ | $\ell_4m_1$ | $\ell_4m_2$ |

Let $F_T$ contain these features on source training pairs and $e_T$ be the repeat-mean residual from the fitted baseline. The source projection gives $G=F_T^\top P_TF_T$ and $b=F_T^\top P_Te_T$. With diagonal scale $C_{jj}=\sqrt{G_{jj}}$ floored at $\max(10^{-10}\max_j\sqrt{G_{jj}},10^{-15})$, the fitted coefficients are

$$\theta=C^{-1}(C^{-1}GC^{-1}+\lambda I)^{-1}C^{-1}b.$$

Candidate penalties are $\lambda\in\{0.01,0.1,1,10\}$. A source fit predicts $F_Q\theta-X_QX_T^+F_T\theta$ before its selected source amplitude is applied. This subtraction extends the endpoint fit learned on source training pairs. The full-screen Han adapter omits the source-indexed endpoint subtraction and projects the shared function on the declared query panel instead; endpoint shifts cancel under that projection.

### 22.3 V4 functional head and source selection

V4 uses the actual `PartnerHeadV4` implementation, adapted to a scalar growth response. Source gene annotations are centred over source identities and each row is normalized by its root-mean-square, floored at $10^{-8}$. Single effects are divided by their source root-mean-square. The RNA query and readout mean are zero; extra RNA features are disabled. The upstream RNA `SourceEncoderV4` and its checkpoints are not used for this functional task.

### Supplementary Table S51 | V4 functional computation and fitting settings

| Element | Exact functional setting |
|---|---|
| Input | Two fixed 64-dimensional annotation vectors and two normalized single effects; zero 64-dimensional query and zero readout mean |
| Environment adapter | Shared $64\to4\to64$ layers with GELU; adapted vector equals input plus adapter output |
| Pair attention | Learned query, key and value maps to 16 dimensions; softmax over the two endpoints |
| Ordered head | Concatenate two adapted endpoint vectors, query, attended 16-vector and three scalars (211 inputs); $211\to16\to8\to1$, GELU after the first two layers, dropout 0.1 after the first |
| Symmetry | Average the same ordered head in both endpoint orders |
| Stored model size | 8,517 parameters, including extra-feature modules retained in the class but disabled for functional inference |
| Source fitting | AdamW, learning rate 0.003, weight decay 0.001 or 0.1; gradient norm capped at 5; checkpoints 50, 150 or 300 updates |
| Source target | Endpoint-projected repeat-mean baseline residual, divided by its training RMS (floor $10^{-10}$); prediction is projected on source training pairs during fitting |
| Seeds and output | Seeds 17 and 29; mean their corrections after restoring each fit's residual scale; dropout disabled for evaluation |

For both families, three whole-combination inner folds refit the baseline and the candidate using only the complementary source training pairs. Candidate predictions retain both repeats for scoring. For each candidate, source amplitude $\beta\geq0$ minimizes pooled inner partner SSE, subject to complete-response SSE no greater than the corresponding baseline. If $e$ denotes baseline minus observation and $d$ the candidate direction, the objective is $\sum_k\|P_k(e_k+\beta d_k)\|^2$, with constraint $\sum_k\|e_k+\beta d_k\|^2\leq\sum_k\|e_k\|^2$. The quadratic feasible interval and constrained minimizer are computed by the supplied `solve_beta` function. Zero correction is always a candidate. Each family selects by partner SSE, then complete SSE and candidate key. This source amplitude is not restricted to one. The Han coefficient $\widehat\alpha\in[0,1]$ subsequently attenuates the already source-scaled correction.

### Supplementary Table S52 | Frozen functional functions used in Han

The selected candidates are refitted on the permitted source training pool; stored weights and source amplitudes are retained in Han. V4 entries refer to the average of both seeds.

| Source background | Function | Selected source candidate | $\beta_{\mathrm{source}}$ |
|---|---|---|---:|
| HeLa | SF12 | Ridge 10 | 3.1308238646 |
| HeLa | V4 | Weight decay 0.1; 300 updates | 0.0821636362 |
| A549 | SF12 | Zero correction | 0.0000000000 |
| A549 | V4 | Weight decay 0.1; 50 updates | 0.0029737686 |
| 293T | SF12 | Ridge 10 | 4.0467336920 |
| 293T | V4 | Weight decay 0.1; 300 updates | 0.0399278590 |

### 22.4 Portable feature coordinates, units and zero rules

The SF12 coordinates for new genes use the source-derived map $W=(E_s-\overline E_s)^+Z_s$. The implementation checks reconstruction of the source coordinates to absolute error below $10^{-8}$, then sets $Z_t=(E_t-\overline E_s)W$. Pair-product centring and scaling continue to use the source training values. V4 likewise uses the source annotation mean followed by per-gene RMS normalization. The missing Han annotation for ADRBK1 uses the stored global GenePT mean policy; Table S48 reports the corresponding identity exclusion.

For each target response repeat, $\sigma_t$ is the standard deviation of the permitted single-effect vector. The source reference $h_s$ is the harmonic mean of the two source repeat single-effect standard deviations. Singles are converted to source units as $s_t h_s/\sigma_t$ before applying the fixed source RMS and feature map. Function outputs, including $\beta_{\mathrm{source}}$, are converted back by $\sigma_t/h_s$. These transformations use single responses and annotations, not query combination outcomes. The two V4 seed outputs restore their own stored residual scales before averaging. A source-selected zero returns zero without evaluating the shared function; the target scalar also returns zero when projected calibration energy is at most $10^{-14}$. A positive finite target single-effect scale is required by the evaluated adapter.

The source training entry point is `code/fit_recovered.py` in the source-fit archive. The frozen inference definitions are `frozen_engine/code/predict.py::sf_features` and `Source`; full-screen adaptation is `code/predict_inputs.py::shared_delta` inside the Han archive. Source endpoint shifts are absent from the shared-function output but disappear after projection on known identities, as checked in `PORTABLE_FEATURE_CHECKS.json`. This projection equivalence concerns correction direction; it does not equate the new target additive baseline with the complete original source predictor.


## Supplementary figure legends

## Supplementary Figure S1 | Risk selection and complete retrospective partner-selection outcomes

**a,** PairPick use case. An anchor's single response defines a fixed transcriptional objective; candidate response predictions yield projected utilities and a partner recommendation. Held-out measurements assess achieved utility in the prescribed pool. **b,** A separate analytic-risk evaluation, not a filtering stage applied to PairPick. The fixed calibration threshold gives accepted/all-query MSE ratio 0.5343 at 49.48% pair-balanced coverage, using 213 records from 33 calibration-disjoint pair identities. The weighted complete-case follow-up gives 0.4775 at 48.95% coverage, using 207 records from 32 pairs. Whiskers are the recorded 95% pair-cluster bootstrap intervals. **c,** Regret for all 28 GSE220974 tasks and all five primary methods, with identical candidate pools and physical-pair holdouts. Rows are ordered by physical anchor, mode-specific anchor and fold, independently of outcome; white separators mark physical-anchor groups. Suffixes a/i denote activation/repression, and the number after the slash is the fold. Dots denote zero regret; nonzero values and the color scale use regret multiplied by 1,000. An inverse-hyperbolic-sine color normalization with linear width 5 makes smaller losses visible while retaining the full range. Neural entries score the mean of three seed predictions. The original design-table audit directly matches 15 of the 49 detected gene-mode labels, identifies one enhancer-only correspondence and leaves 33 unlisted as designed pairs. All tasks are retained as a detected-label demonstration; label matching does not verify cell-level constructs. **d,** R-G-minus-reference achieved-utility differences with 95% bootstrap intervals clustering the eight physical anchors. R-G and the archived PRESAGE adapter choose identically in all 28 tasks. Complete regret, hit, win/tie/loss and prediction-error tables are in Supplementary Tables S17–S19 and S28. These descriptive outcomes do not establish superiority across all methods or independent biological confirmation of intended constructs.

## Supplementary Figure S2 | All double-only comparator intervals

Mean operation-rank differences for WitnessCell minus each of the fifteen eligible comparator entries, including simple references. All methods use the same 540 double-perturbation operation records and the original Top-100 six-metric rule. Whiskers are pointwise central 95% intervals from 20,000 resamples over 287 dataset-pair clusters, carrying repeated seed records together. All upper bounds are below zero. The intervals condition on the saved score archive; they do not provide simultaneous, multiplicity-adjusted coverage or incorporate the full development history and future-study variation. These are fifteen comparisons of the same prediction population, not fifteen independent biological confirmations.

## Supplementary Figure S3 | Complete-response benefit from additional combination support

Wessels gains relative to the same fixed nonlinear backbone as the retained training-combination residual library increases from 25% to 50% and 100%. Points and lines show the recorded MSE and weighted-MSE gains. Shading shows archived 95% canonical-pair bootstrap intervals from 5,000 draws. This analysis evaluates endpoint refinement of complete responses; it does not establish a new partner contrast. Lines connect the three evaluated fractions and do not validate optimal acquisition or extrapolate an experiment-saving rate.

## Supplementary Figure S4 | Complete functional control and targeted sensitivity checks

**a,** Original and component-controlled partner gains for all twenty-seven configurations with nonzero selected correction. The twenty-eight zero-selected configurations remain unchanged. The diagonal denotes equal performance. These are configurations, including overlapping K562 layouts, rather than independent biological replications. **b,** Coefficients derived from each internal measured repeat for the same twenty-seven configurations; underlying response fits are held fixed. **c,** Partner gains when the fixed query identities are processed together, as identity-sorted alternating halves or thirds, or separately. All assembled predictions are scored on the original full query projector. Points summarize SSE-pooled conditions and are not uncertainty intervals. Full source tables retain all 55 configurations, both outcome repeats, and all three training-coefficient choices.

## Supplementary Figure S5 | Target-side calibration across all frozen shared partner functions

**a,** Partner-error reduction from target-side calibration relative to the unchanged correction. **b,** Partner-error reduction relative to removing the correction. **c,** Complete-response error reduction relative to Han single-response addition. Every source (HeLa, A549 and 293T) and portable family (SF12 and V4) is shown. Teal filled diamonds denote reconstructed outcomes; blue open circles denote published phenotypes. Horizontal segments span the minimum and maximum over five overlapping pair partitions; large points are pooled-SSE estimates. Segments are sensitivity ranges, not confidence intervals. **d,** Target-estimated correction strengths across the fifteen folds for each readout; small points are fold values and the large point is their median. Source-selected zero correction remains zero. The scientific comparison uses the same baseline and direction within each configuration, changing only the target-estimated scalar. Neither source configurations nor overlapping partitions are independent target studies. Tables S44-S48 specify coverage, full estimates, repeat results and access history.

## Supplementary reference correspondence

Reference numbers correspond to the main manuscript.

**[3]** Horlbeck, M. A. et al. Mapping the genetic landscape of human cells. *Cell* **174**, 953–967.e22 (2018). doi:10.1016/j.cell.2018.06.010.

**[8]** Roohani, Y., Huang, K. & Leskovec, J. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nature Biotechnology* **42**, 927–935 (2024). doi:10.1038/s41587-023-01905-6.

**[9]** Wei, Z. et al. Benchmarking algorithms for generalizable single-cell perturbation response prediction. *Nature Methods* **23**, 451–464 (2026). doi:10.1038/s41592-025-02980-0.

**[11]** Viñas Torné, R. et al. Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation. *Nature Biotechnology* **44**, 1050–1059 (2026). doi:10.1038/s41587-025-02777-8.

**[12]** Miller, H. E. et al. Deep learning perturbation models can outperform baselines on calibrated metrics. *Nature Biotechnology* (2026). doi:10.1038/s41587-026-03307-w.

**[13]** Peters, J., Janzing, D. & Schölkopf, B. *Elements of Causal Inference: Foundations and Learning Algorithms* (MIT Press, 2017). https://mitpress.mit.edu/9780262037310/elements-of-causal-inference/.

**[14]** Robinson, G. K. That BLUP is a good thing: the estimation of random effects. *Statistical Science* **6**, 15–32 (1991). doi:10.1214/ss/1177011926.

**[15]** Rao, C. R. Unified theory of linear estimation. *Sankhyā Series A* **33**, 371–394 (1971).

**[18]** Zhu, O. & Li, J. Scouter predicts transcriptional responses to genetic perturbations with large language model embeddings. *Nature Computational Science* **6**, 21–28 (2026). doi:10.1038/s43588-025-00912-8.

**[19]** Mazumder, R., Hastie, T. & Tibshirani, R. Spectral regularization algorithms for learning large incomplete matrices. *Journal of Machine Learning Research* **11**, 2287–2322 (2010). https://jmlr.org/papers/v11/mazumder10a.html.

**[20]** Shen, J. P. et al. Combinatorial CRISPR–Cas9 screens for de novo mapping of genetic interactions. *Nature Methods* **14**, 573–576 (2017). doi:10.1038/nmeth.4225.

**[21]** Dixit, A. et al. Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. *Cell* **167**, 1853–1866.e17 (2016). doi:10.1016/j.cell.2016.11.038.

**[22]** Mathur, L. et al. Combi-seq for multiplexed transcriptome-based profiling of drug combinations using deterministic barcoding in single-cell droplets. *Nature Communications* **13**, 4450 (2022). doi:10.1038/s41467-022-32197-0.

**[23]** Pacalin, N. M. et al. Bidirectional epigenetic editing reveals hierarchies in gene regulation. *Nature Biotechnology* **43**, 355–368 (2025). doi:10.1038/s41587-024-02213-3.

**[39]** Schölkopf, B., Herbrich, R. & Smola, A. J. A generalized representer theorem. In *Computational Learning Theory*, Lecture Notes in Computer Science **2111**, 416–426 (Springer, 2001). doi:10.1007/3-540-44581-1_27.

**[40]** Candès, E. J. & Recht, B. Exact matrix completion via convex optimization. *Foundations of Computational Mathematics* **9**, 717–772 (2009). doi:10.1007/s10208-009-9045-5.

**[41]** Cai, J.-F., Candès, E. J. & Shen, Z. A singular value thresholding algorithm for matrix completion. *SIAM Journal on Optimization* **20**, 1956–1982 (2010). doi:10.1137/080738970.

**[42]** Efron, B. Bootstrap methods: another look at the jackknife. *The Annals of Statistics* **7**, 1–26 (1979). doi:10.1214/aos/1176344552.

**[43]** Han, K. et al. Synergistic drug combinations for cancer identified in a CRISPR screen for pairwise genetic interactions. *Nature Biotechnology* **35**, 463–474 (2017). doi:10.1038/nbt.3834.
