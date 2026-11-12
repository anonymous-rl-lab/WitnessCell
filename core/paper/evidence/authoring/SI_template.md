# Supplementary Information

## WitnessCell predicts combinatorial perturbation responses using experimental support

**Companion to the v24 main text and Methods.** This supplement supplies derivations, implementation details, complete comparisons and sensitivity analyses for the six Results sections. It preserves the identities and evaluation status of the experiments in the accompanying `grn_repo_v27` archive. The manuscript name **PairPick (R-G)** refers to the unchanged archived model `R_G`.

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

Tables S1–S24 are embedded alongside the corresponding notes. The accompanying source bundle contains full-precision tables, relevant protocol records, a claim-to-source index and the script used to assemble this document. Bibliographic numbers refer to the v24 main-text reference list; selected references are repeated at the end for convenience.

## Supplementary Note 1 | Method identities, evidence populations and notation

### 1.1 Scientific and computational correspondence

WitnessCell is the response-prediction method developed in this study. Its evidence has several levels: a controlled causal system establishes an information property; the released transcriptomic implementation supplies benchmark performance; RNA and growth extensions test additional response components; and selection analyses test uses of the predicted responses. These levels share a scientific motivation but have distinct estimators, training sets and outcomes. Table S1 identifies the implementation responsible for each claim. In particular, the formal benchmark rank belongs to the released response predictor. It is not a measured rank of the later RNA partner head or PairPick ranking extensions.

### Supplementary Table S1 | Main-text claims and their experimental identities

| Main-text result | Evaluated implementation | Evidence population | Principal outcome and interpretation |
|---|---|---|---|
| Method construction | Released endpoint-transfer core; explicitly separate partner extensions | Defined design matrices and archived fitted outputs | Algebraic capability of each response component |
| Joint-intervention information | Candidate least squares and contrast GLS in a specified causal family | Four worlds, two degree-matched panels, 1,000 paired noise seeds | Target risk and class-conditional target confidence sets |
| Formal prediction performance | Released WitnessCell core | Four benchmark datasets, three public split seeds; 654 operation–seed records | First of sixteen eligible entries under the specified aggregation |
| RNA partner learning | Matched M0, M1 and reliability-weighted R | 55 pairs, 11 same-fit rectangles, 4,995 non-target readouts | Added response capability on a developmental panel |
| Functional prediction | Repeat-specific baseline plus shared rank-one or spectral correction | K562 development; procedure-frozen Jurkat evaluation of 617 pairs | Original-unit growth error and endpoint-eliminated error |
| Risk selection | Analytic covariance predictor; separate empirical ranking of released predictions | Pair-disjoint calibration/evaluation for the analytic experiment; released split comparisons for empirical ranking | Accepted error at reported coverage, with implementation-specific assumptions |
| Partner recommendation | PairPick using R-G response predictions | GSE220974: 49 detected double identities, 28 selection tasks from eight physical anchors | Observed utility, regret and empirical-best hits |

The Wessels benchmark and the matched CaRPool RNA panel reuse a biological source. K562 and Jurkat are different cell contexts within the same growth-screen study. GSE220974, Dixit and Combi-seq add other sources, but their readouts, interventions and evaluation objectives differ. Accordingly, the number of analyses is not a count of independent biological confirmations.

### 1.2 Common notation and loss conventions

Let $T$ denote permitted training combinations and $Q$ held-out query combinations. Rows of $Y_T$ or $Y_Q$ are condition-level responses; columns are readouts. For RNA, a response is expressed on the experiment's retained control-centred axis. For growth, original-unit errors use control-centred log2 enrichment normalized by population doublings. An intervention identity includes activation or repression mode when applicable.

The unsigned incidence matrix $X$ has one row per distinct-endpoint pair and one column per endpoint identity, with two unit entries in each row. $X^+$ denotes its Moore–Penrose pseudoinverse. The endpoint-eliminating projector is $P=I-XX^+$. A complete error is $\|\widehat Y-Y\|_F^2$; an endpoint-eliminated error is $\|P(\widehat Y-Y)\|_F^2$. Division by the stated number of coordinates, or by contrast dimensions and readouts, defines the corresponding reported mean. A projected *ratio* can instead divide by $\|PY\|_F^2$, as in the structural tables. These quantities must not be interchanged.

Except where a table explicitly reports a loss ratio, a gain is $100(1-L_{\mathrm{method}}/L_{\mathrm{reference}})$. Negative gain denotes increased error. Pooling squared errors before forming a ratio differs from averaging fold-wise percentage gains. Every table retains its original convention. The historical weighted benchmark loss, raw non-target RNA loss, rectangle loss, growth loss and selection regret measure different outcomes.

## Supplementary Note 2 | Target identification by joint interventions

### 2.1 Compatible mechanisms and target-specific identification

Let $\mathcal M$ be a declared mechanism class and $D_{\mathcal E}(m)$ the population observation law under measured interventions $\mathcal E$. The equivalence class of a mechanism is

$$[m]_{\mathcal E}=\{m'\in\mathcal M:D_{\mathcal E}(m')=D_{\mathcal E}(m)\}.$$

A target $q(m)$ is identified at this observation law exactly when it is constant on $[m]_{\mathcal E}$. This definition concerns the target, rather than identification of all parameters or edges. Adding an intervention can partition the equivalence class, but it helps this target only if the removed ambiguity changes $q$.

For the special linear case with population observations $D\theta$ and target $c^\top\theta$, identification holds if and only if $c\in\operatorname{row}(D)$. To prove sufficiency, write $c=D^\top a$ and obtain $c^\top\theta=a^\top D\theta$. For necessity, if $c$ is not in the row space, it has a nonzero component $v$ in $\ker D$ with $c^\top v\ne0$. Then $\theta$ and $\theta+v$ have identical observations and different targets. This is a fixed-family estimability statement [13]; biological model misspecification remains a separate issue.

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

{{FORMAL_RISK}}

Each risk first averages the four worlds within a noise seed. Intervals are two-sided Student-$t$ Monte Carlo intervals across 1,000 paired seeds, with 999 degrees of freedom. They describe simulation variability, not variability across biological systems. The point-prediction risk includes every retained fit, including those with empty confidence sets.

### 2.5 Simultaneous candidate-family confidence sets

Whiten the stratum observations and project out the eight-dimensional nuisance column space. An orthonormal basis of the remaining gate span yields $Z\sim N(a_k\beta_h,I_2)$ in the informative panel. The ball

$$\{\mu:\|Z-\mu\|^2\le c\},\qquad c=\chi^2_{2,.95}=5.9914645,$$

contains the true projected mean with probability 0.95 under the stipulated model. For each candidate line, let $\widehat b_k=a_k^\top Z/\|a_k\|^2$ and $d_k^2=\|Z-a_k\widehat b_k\|^2$. If $d_k^2\le c$, its retained coefficients are

$$I_k=\left[\widehat b_k-\sqrt{\frac{c-d_k^2}{\|a_k\|^2}},\ \widehat b_k+\sqrt{\frac{c-d_k^2}{\|a_k\|^2}}\right].$$

Otherwise that candidate contributes an empty set. The target set is the union of $-pI_A$, $\{0\}$ when gate B is retained, and $-pI_C$. On the event that the ball contains the true mean, the true candidate coefficient and hence the true target are retained. Target coverage is thus at least 0.95 within this family. This construction includes all candidates; using only the best-fitting gate would invalidate this argument.

Coverage and zero exclusion use the exact union, not its convex hull. An empty union is a coverage failure and does not certify a nonzero target. The comparison panel leaves relevant gate coefficients unidentified and yields an unbounded target set. Table S3 reports all four informative-panel worlds. Gate B illustrates that a zero target can be well supported despite nonzero interaction signals in measured conditions.

### Supplementary Table S3 | Informative-panel world-specific estimation and confidence sets

{{FORMAL_WORLDS}}

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

{{BUDGET}}

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

Table S5 contains all sixteen eligible entries in the formal combination track of scPerturBench [10]. The primary operation ranking and the mean metric rank are distinct aggregations; the former defines the reported first-place position. The 654 operation–seed records include 540 double-perturbation records and 114 unseen-single records. This scope explains why the combination-only loss analysis has a different denominator. It also explains why `linearModel` and `trainMean` can have different overall ranks despite identical predictions in their released combination branches.

### Supplementary Table S5 | Complete formal top-100 leaderboard

{{LEADERBOARD}}

### Supplementary Table S6 | Complete metric-specific positions in the formal track

{{METRIC_POSITIONS}}

Positions in Table S6 rank the mean metric ranks in the source file; they are not raw metric values. The full-precision mean ranks are retained in the source CSV. WitnessCell's first-place aggregate position coexists with weaker Wasserstein and common-differential-expression positions. Distributional scores use the benchmark's adapter combining predicted means with training-control variation. They do not demonstrate a fitted conditional distribution of perturbed cell states.

Equal dataset weighting gives a WitnessCell top-100 mean operation rank of 3.873152 and first position. Aggregation across both recorded gene panels gives 3.396789, and combining both panels with equal dataset weighting gives 3.456997; both remain first. These are robustness summaries of the same predictions, not independent replications. Table S7 retains paired uncertainty against the closest aggregate comparators. Bootstrap units are dataset–perturbation clusters, carrying their repeated seed records together.

### Supplementary Table S7 | Paired cluster-bootstrap differences in mean operation rank

{{BENCHMARK_PAIRED}}

The source analysis uses 20,000 bootstrap draws across 384 dataset–perturbation clusters. A bootstrap fraction of one favouring WitnessCell is an empirical resampling frequency, not a posterior probability or a statement that the population $P$ value is exactly zero.

### 3.3 Absolute losses, readout masks and support amount

The combination-only weighted analysis compares with the released training-mean additive reference. Its `linearModel` combination branch yields the same prediction, so this analysis counts that reference once. The 24.87% improvement is not an improvement relative to the strongest neural model. Historical weighting uses absolute target-versus-control Welch statistics from one cell half, min–max scaling and squaring, followed by scoring against the other half. Target coordinates are included in that historical endpoint.

Later non-target controls remove the query's perturbed genes, renormalize weights, use full-cell condition means and geometrically aggregate seed-wise ratios. Their pair bootstrap carries repeated records of the same perturbation together. Their numbers therefore need not equal the historical weighted ratios. Table S8 preserves both sets of outcomes with their proper denominators.

The residual-library comparison holds the fitted endpoint system and backbone fixed while changing the amount of available combination residual information. Its quarter, half and full libraries answer how much that branch benefits from support. They do not measure de novo partner contrasts, because the branch remains inside the endpoint span. All nine recorded Wessels MSE order paths and six of nine weighted-MSE order paths are monotone; the mixed weighted result is retained.

### Supplementary Table S8 | Absolute-error controls and residual-library contribution

{{BREADTH}}

Sources: `publication/released_benchmark/`, `publication/Fig2_prediction_ratios.csv`, `publication/ED1_coverage_comparisons.csv` and the `publication/Fig3_ED4_witness_*` tables in the SI bundle. The broader benchmark inventory is not substituted for the sixteen-entry eligibility rule.

## Supplementary Note 4 | Endpoint elimination and the scope of partner contrasts

### 4.1 Endpoint-transfer cancellation

For the released residual branch,

$$\delta_Q=\gamma\frac{X_QX_T^\top}{2}\left(\frac{X_TX_T^\top}{2}+\eta I\right)^{-1}R_T=X_Q\Theta,$$

where the final equality defines $\Theta$. Since $P_QX_Q=0$, $P_Q\delta_Q=0$ for every residual matrix and every selected amplitude. This is a structural identity independent of whether the training residuals contain interaction effects. Those effects can influence the estimated endpoint coefficients, but this output branch still lies in the query endpoint span.

The nonlinear backbone can have $P_QF_{\mathrm{core}}\ne0$. Consequently, cancellation of the residual branch does not imply that the complete predictor is additive. It identifies which component accounts for any new contrast. The twelve released splits verify equality of projected final and paired-backbone outputs within saved-array precision, with maximum recorded discrepancy $6.33775\times10^{-8}$.

### 4.2 Orthogonal error decomposition and weighted losses

The matrix $P_Q$ is symmetric and idempotent. For $E=\widehat Y-Y$,

$$\|E\|_F^2=\|(I-P_Q)E\|_F^2+\|P_QE\|_F^2,$$

because $P_Q(I-P_Q)=0$. Complete-response improvement may therefore occur wholly in the first term. If a common diagonal readout scaling is applied to every query, the decomposition also holds after that column scaling. Query-dependent element weights generally change the relevant inner product, and the same unweighted projector need not supply a weighted orthogonal decomposition. The structural tables use their stated common readout axis; the historical weighted benchmark is not silently decomposed by this identity.

### 4.3 Query-graph dimension

Take the graph whose vertices are active endpoint identities and whose edges are query pairs. For a connected component, $Xv=0$ requires $v_g=-v_h$ along every edge. A connected bipartite component admits one alternating-sign null vector; a connected non-bipartite component contains an odd cycle, forcing that vector to zero. If the graph has $n$ active vertices, $m$ edges and $b$ bipartite connected components, then

$$\operatorname{rank}(X)=n-b,\qquad \operatorname{rank}(P)=m-n+b.$$

A tree has no contrast dimension. An isolated odd cycle also has none, whereas an even cycle has one. More complex graphs can have endpoint-eliminated directions without an isolated four-edge rectangle; the rank formula, rather than counting named rectangles, defines the full contrast space. The twelve benchmark graphs have 212 such dimensions. These are linear-algebra dimensions of query sets, not 212 biological replicates.

### 4.4 Rectangles and full projection

For four pairs $AB,AD,CB,CD$, put $c=(1,-1,-1,1)^\top$. Then $c^\top X=0$. If these four edges form the entire query graph, $P=cc^\top/4$ and $\|PE\|_F^2=\|c^\top E\|^2/4$. With multiple or overlapping rectangles, summing their contrast losses weights directions according to the chosen rectangle list, while the orthogonal projector weights an orthonormal basis. The two reported endpoints are therefore complementary, rather than numerically interchangeable.

All four predictions in a rectangle must come from the same fitted model. Mixing predictions from different fits could introduce fitting differences into the contrast. Query projection uses identities to define the operator and observed outcomes only for evaluation. It is not used to tune a held-out correction. Both the true response and its measurement noise survive projection; endpoint removal is not itself a denoising operation.

### Supplementary Table S9 | Complete released-split structural decomposition

{{STRUCTURAL}}

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

{{RNA_MATCHED}}

### Supplementary Table S11 | Matched RNA fold composition and fixed-model gains

{{RNA_FOLDS}}

The rows in Table S11 are derived from the seed-zero records of fixed models; copies of those same deterministic predictions stored alongside the other neural seeds are not counted again. Fold gains use their own baseline SSE. The pooled primary gain uses the sum of rectangle errors, rather than the mean of these fold gains.

### Supplementary Table S12 | All eleven same-fit RNA rectangles

{{RNA_RECTANGLES}}

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

{{RNA_CONTROLS}}

### 6.4 Interpretation of named readouts

The S100A9 and LYZ illustrations were selected after aggregate scoring to explain particular large changes. Their displayed beneficial coordinates rank first and tenth among 54,945 rectangle–readout coordinates by the archived absolute-benefit criterion. Their examples are not independent validation. The adverse LYZ rectangle involving GFI1/RCOR1 and RING1/SMARCD1 is retained, with a shrinkage weight of 0.879965. Some favourable examples are predicted more closely by M1 than by R, so a plot of those coordinates alone cannot explain the aggregate benefit of readout weighting. The bundle retains all eleven contrasts for the selected readouts and the four-corner values rather than only favourable points.

The modest cross-readout Spearman and direction-agreement results in Table S13 further delimit interpretation. A lower raw squared error supports improved prediction of the measured response scale. It does not, by itself, establish a transcriptome-wide molecular interaction map or causal direction for individual gene pairs.

## Supplementary Note 7 | Functional growth prediction and procedure transfer

### 7.1 Response construction and protected evaluation sets

The growth measurements derive from the Triple Sequencing screen [4]. Guide orientations are averaged within physical guide pair and source repeat, then equally across guide pairs. Processing requires finite measurements and both repeats; conflicting duplicate orientations stop processing. Published whole-map interaction scores are excluded from the fitting target because their fitted baseline could incorporate held-out combinations. The observed target is the source growth phenotype in its original units.

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

{{GROWTH_K562}}

### 7.3 Jurkat results and conditional sensitivity

Rank one is the prespecified primary functional extension and spectral correction is secondary. Both improve pooled error and each source repeat. Projecting each error through the fixed query graph gives 538 contrast dimensions and shows that the benefit is not entirely an endpoint-space change. These dimensions do not count independent biological observations.

The gene-reweighting calculation fixes the predictions and observed outcomes and varies endpoint weights over 10,000 paired draws. Its intervals describe dependence on the observed gene composition; they omit model refitting and do not establish coverage over newly sampled cell lines. The 79 leave-one-gene-out evaluations are overlapping deletions from the same test population. They show that the gain persists when any one observed endpoint's pairs are removed, rather than supplying 79 independent confirmations.

The M0 replicate-mean comparator is scored against both observed repeats. Its deterioration in complete error shows that averaging baseline predictions alone does not explain the partner model's improvement. The same comparator has a different projected outcome, reinforcing the need to retain complete and endpoint-eliminated losses separately.

### Supplementary Table S15 | Complete Jurkat primary, secondary and sensitivity results

{{GROWTH_JURKAT}}

Named combinations illustrate the aggregate evaluation. For example, TTI2–PGK1 has a mean observed response near −0.321, baseline prediction −0.274 and corrected prediction −0.287, whereas TTI2–OGFOD1 has values near −0.036, −0.069 and −0.058. These are descriptive selected examples. The source bundle retains repeat-specific predictions and all evaluated pair results; these examples are not validation of RNA-based recommendations in Jurkat.

## Supplementary Note 8 | Model-based risk and empirical selection

### 8.1 Derivation of the linear minimum-risk coefficient

Consider $y=X\beta+u+\varepsilon$ and target $y_t=t^\top\beta+u_t$. Put $\Sigma=\operatorname{Cov}(u+\varepsilon)$, $k_t=\operatorname{Cov}(u+\varepsilon,u_t)$ and $k_{tt}=\operatorname{Var}(u_t)$. For a predictor $a^\top y$, unbiasedness for every $\beta$ requires $X^\top a=t$. Under the stated second moments, its risk is

$$R(a)=a^\top\Sigma a-2a^\top k_t+k_{tt}.$$

Assume $\Sigma$ is positive definite and the retained design columns have full rank. Completing the square gives

$$R(a)=k_{tt}-k_t^\top\Sigma^{-1}k_t+(a-\Sigma^{-1}k_t)^\top\Sigma(a-\Sigma^{-1}k_t).$$

Set $b=a-\Sigma^{-1}k_t$ and $\widetilde t=t-X^\top\Sigma^{-1}k_t$. Minimizing $b^\top\Sigma b$ subject to $X^\top b=\widetilde t$ by Lagrange multipliers gives

$$b=\Sigma^{-1}X(X^\top\Sigma^{-1}X)^{-1}\widetilde t.$$

Substitution yields the coefficient and risk $W_t$ in the main Methods [12,13]. A rank-deficient design can first be expressed in a full-rank basis, provided the requested endpoint target is estimable. Substituting a pseudoinverse without checking estimability does not resolve an unsupported target.

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

{{RISK_TABLES}}

## Supplementary Note 9 | PairPick: response constraints, candidate pools and completed comparisons

### 9.1 Cohort construction and the intended decision

PairPick takes an anchor intervention, a fixed transcriptional objective and a prescribed candidate partner pool. It returns response predictions, predicted utilities and the highest-utility candidate. Held-out measurements assess that recommendation. The GSE220974 evaluation uses the source's mixed activation/repression measurements [20], keeping mode identities distinct throughout preprocessing and fitting.

The operational cohort contains all 49 detected double labels involving two distinct regular genes, representing 28 physical pairs. Enhancer and same-gene doubles are excluded. Control identities are reconstructed from guides, removing 860 enhancer-target cells previously assigned to the control category. The source's intended designed-construct count is smaller than the detected-label cohort. A complete intended-construct whitelist was not verified, so incidental co-delivery, assignment error or doublets remain possible within the operational labels. The cohort was not filtered on response strength to force agreement with the source design count.

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

{{GSE_ALL}}

An em dash in a response-score column denotes an unavailable response-level score for a ranking-only output, not zero error. An em dash in the primary-set column identifies a retained secondary or historical row. Larger-pool regret is evaluated on the same five tasks with at least four candidates.

### Supplementary Table S18 | Paired utility differences for R-G and retained references

{{GSE_PAIRED}}

Intervals in Table S18 resample eight physical anchors and are conditional on the observed pools and predictions. With few clusters and repeated ties, bootstrap endpoints at zero do not establish a general strict advantage. No candidate set or anchor was removed to improve the result. Complete task-level choices are retained in `gse220974/CHOICES.csv` (420 rows: fifteen models by 28 tasks).

### Supplementary Table S19 | All fixed-anchor candidate pools

{{GSE_GROUPS}}

Suffixes `a` and `i` denote activation and repression. Pool membership is fixed by the archived split, not optimized for a method. The five larger pools are a subset of the 28 tasks and are not counted as additional observations.

### 9.5 Completed neural training and mode adapters

GEARS uses disjoint mode-labelled copies of its perturbation graph and a physical-gene coexpression graph inferred from permitted controls and singles. Scouter uses direct-sum GenePT embeddings, activation $[e,0]$ and repression $[0,e]$. These adapters expose mode identity to the upstream architectures; their results are not an unmodified mixed-mode capability claimed by the upstream papers [6,7]. Standard ridge uses an unpenalized intercept and separate mode-specific endpoint columns, with equal condition weights for controls, singles and training doubles.

GEARS uses Adam at .001, weight decay .0005 and batch size 32; its two learning-rate-halving schedules have intervals one/five and caps twenty/sixty, with patience ten. Scouter uses Adam at .001, batch size 256 and exponential learning-rate factor .9. Its forty-epoch selection run has a predefined possible extension to eighty when the recorded late-improvement criterion is met. Table S20 reports the actual selected durations; no GSE220974 Scouter run took that extension. Each final fit restarts and trains for the selected epoch count.

### Supplementary Table S20 | Actual GSE220974 final-fit durations for all thirty jobs

{{GSE_EPOCHS}}

The completion record verifies thirty frozen, independently archived and replayed jobs: five folds by three seeds for each model. GEARS folds three and four for all three seeds were completed on GPU; the other final jobs used the recorded CPU environment. GPU execution used FP32 without AMP or TF32. The archive records the CPU PyTorch 2.14.1 and GPU PyTorch 2.1.2+cu118 environments, with the GPU runtime's Python 3.10.8 and NumPy 1.26.4. These differences are part of execution provenance; they do not establish identical stochastic trajectories across environments.

For each neural model, the primary prediction is the arithmetic mean of three saved prediction vectors, followed by RNA, partner and decision scoring. In general, scoring this mean differs from averaging seed-specific losses or recommendations. Some alternative models' outer scores were known before neural completion. Fixed neural budgets and selection rules were retained, but the whole comparison is reported as retrospective source evaluation rather than an untouched preregistered confirmation.

## Supplementary Note 10 | Independent-source checks and acquisition results

### 10.1 Frozen Dixit RNA comparison and subsequent development

The independent RNA check uses ten permitted TF singles, 45 pairs and 23,519 non-target readouts from Dixit [17]. Three outer folds each withhold fifteen pairs, with all baseline and correction selection restricted to the remaining thirty. Each query graph has five endpoint-eliminated dimensions. A predefined G2M program maps 52 of 54 unique symbols and supplies a secondary program endpoint.

The selected uniform correction improves complete RNA error by 0.26980% but increases the primary endpoint-eliminated error by 0.03260%. M0 itself substantially improves complete response relative to measured additivity; that baseline performance does not establish the added partner capability. Table S21A preserves the full frozen comparison, its fold outcomes and conditional sensitivity. The failed primary partner result remains the source-transfer conclusion.

Subsequent G2M development introduced a fixed-reference continuation. A projection defined from the current query batch can change predictions when an unrelated query is added or removed. The fixed full-roster construction instead locks the reference identities, giving a query-batching-invariant mapping. Its 14.2994% projected gain and 2.2671% complete gain are relative to the independently fitted scalar saturation baseline in that later program task. Two implementations with those identical values are mathematically equivalent and are counted once. This development result does not revise the earlier frozen all-RNA decision.

### 10.2 Frozen Combi-seq program evaluation

The Combi-seq evaluation [18] contains 374 drug pairs, 36 endpoints and 51 program genes. Three outer folds have 125, 125 and 124 query pairs; their full query graphs supply 266 contrast dimensions. The principal system adds the fixed-reference program correction to the shared endpoint baseline, keeping the same outer scaling and program units for all methods.

Its complete gain relative to that endpoint baseline is −0.023282%, and its projected gain is −0.030240%. Comparison with scalar saturation alone gives a small complete gain, but the endpoint baseline already accounts for that improvement. The proper added-component comparison therefore remains adverse. Three retained experimental-run targets are all evaluated against the same predictions fitted from pooled training means. Their mixed results, including two adverse complete-error changes, are retained in Table S21D. These targets are repeated measurements of the same pairs, not three independently selected success tests.

### Supplementary Table S21 | Complete independent-source and subsequent program results

{{EXTERNAL_TABLES}}

### 10.3 One-measurement acquisition: development and outer evaluation

The acquisition experiment asks which additional measured condition to acquire for a specified target, a different decision from PairPick's choice among candidate interventions based on predicted utility. The target-aligned development procedure improves partner error by 4.68324% and target-AB error by 22.52110% relative to the expected uniform acquisition reference. Its primary development run gate required at least 5% partner improvement and was not met. An explicit subsequent continuation retained that result and evaluated the outer targets rather than redefining the original gate as passed.

On the six outer targets the corresponding improvements are negative: −2.41963% for partner error and −46.91288% for target-AB error. Table S22 keeps both stages. No outer target is removed from the aggregate decision. A favourable geometric-width change or a useful development example is not substituted for measured outer predictive utility. The experiment therefore supplies a boundary on active acquisition, not an additional successful end-to-end application of PairPick.

### Supplementary Table S22 | Development and outer acquisition outcomes

{{ACQUISITION}}

### 10.4 Related development analyses and evidence accounting

Earlier fixed-G integration, measurement-reference calibration and shared-feature development are retained in the code archive and the SI source bundle. Measurement-reference calibration concerns better input estimation for M0; it is not evidence of a new partner head. The qualified shared-feature tasks reuse overlapping rectangles, and their case-conditioned comparator has adverse results. Their full tables are provided as `additional_development_tables.md` with their original task labels. These analyses are not pooled with the frozen Dixit or Combi-seq tests or counted as independent confirmations of the v24 partner claim.

The biological results establish different levels of transfer: a substantial matched-panel raw RNA gain, a procedure-frozen growth gain in another cell context, a useful source-specific recommendation result and unresolved or adverse source-transfer results for the specified RNA/program extensions. This differentiation defines the current empirical scope of the common response construction.

## Supplementary Note 11 | Statistical units, original outcomes and amendments

### 11.1 Unit of inference

The synthetic experiment uses paired random-noise seeds as the Monte Carlo unit, averaging the four generating worlds inside each seed. The benchmark's paired resampling groups dataset and perturbation identity and keeps repeated split records together. Historical combination-error resampling groups canonical pairs. RNA cell halves estimate training agreement; they are not independent biological replicates. Rectangle overlaps and shared endpoints are reported as dependence, and readouts are not substituted for independent screens.

Growth errors pool the two measured repeats and retain each repeat's direction. Gene-reweighting and endpoint-deletion analyses condition on the observed screen. GSE220974 selection intervals group by physical anchor, not the 28 nominal tasks. The reported intervals therefore have different meanings: Monte Carlo estimation uncertainty, conditional pair/anchor resampling uncertainty or fixed-prediction sensitivity. They are labelled at the point of use.

No universal significance claim is formed by counting all favourable secondary analyses. The primary endpoint and comparator of each experiment define its decision. Descriptive seed comparisons, selected molecular examples, scaling changes and alternative ranking heads remain supplementary evidence. In particular, an improvement in complete RNA error does not establish partner improvement, and neither error automatically establishes lower recommendation regret.

### 11.2 Original plans and local freezes

The original M1 developmental target was at least 10% rectangle improvement; its 9.4382% outcome did not meet that criterion. The later R result is development on the same panel. Original registered K562 experiments selected zero correction at all five budgets. Subsequent revised K562 results use a different evaluation pool and retain developmental status. The Jurkat transfer procedure was locally frozen before fitting and scoring after an availability-only amendment removing GINS4.

Table S23 preserves this chronology. Times are explicitly UTC. A local prediction lock identifies when predictions were fixed relative to scoring. A later external deposit documents those records but cannot retroactively preregister the revised method. Similarly, the GSE220974 completion records prove completed jobs and retained outputs, not that no related outer score had previously been seen.

### Supplementary Table S23 | Original outcomes, amendments and functional prediction locks

{{CHRONOLOGY}}

One historical pre-amendment protocol ZIP is referenced by the original records but is not supplied as a recovered file in the present compact materials. The available freeze receipts, amendment text and scoring records are retained. The current source map does not convert that historical pointer into a claim that every earlier archive has been recovered.

### 11.3 Interpretation of the scientific sequence

The controlled system proves that intervention information can distinguish target answers under its model class. The benchmark establishes competitive prediction by the released core. The component analyses identify what the extensions can represent and where their gains occur. Functional and decision evaluations then measure particular uses. This is a sequence of complementary scientific questions, not an experimentally established causal mediation chain from network recovery through prediction to recommendation. Each stage has its own direct observable and comparison.

## Supplementary Note 12 | Source correspondence and computational reproduction

### 12.1 Artifact and version correspondence

This SI accompanies `WitnessCell_v24_Main_Methods.md`. It reorganizes existing evidence and adds explicit derivations of the stated constructions; it does not retrain a model, alter a frozen prediction or replace an adverse result. The code distribution is `grn_repo_v27.zip`. Its publication renderer still targets the legacy v22 figures and tables. Running its `paper` task does not produce a new v24 figure set, and code release 27 is not manuscript version 27.

The compact SI support bundle contains the unchanged v24 main text, this supplement, selected source tables and protocols, `CLAIM_SOURCE_MAP.csv`, a source SHA-256 manifest, and the document assembly/check script. Long tables remain machine-readable: complete benchmark panels and metrics, RNA fold and rectangle results, all GSE220974 choices and selections, and the source-specific confirmation comparisons. The bundle is a traceable document-and-source delivery, not a second complete training repository.

### 12.2 Reproduction levels

Four operations should be distinguished. Score replay starts from saved predictions, residuals or sufficient statistics. Weight replay starts from saved model parameters and required input assets. Compact refitting trains a declared model on retained condition-level data. Full raw-cell retraining additionally reconstructs preprocessing, graphs, embeddings, data splits and stochastic optimization. Table S24 records which operation each code entry actually performs. A hash check validates saved bytes; it does not by itself verify a scientific assumption.

From the extracted code archive, a task is invoked as:

```bash
python tools/reproduce.py --task gse220974_final --out /absolute/path/to/new_output
```

Other task names are listed in Table S24. Each output should be a fresh location. Large raw matrices and neural weights follow their external-asset manifests. The GSE220974 compact increment stores lossless numerical arrays and external weight identities; it does not include every raw-cell input or all model-weight bytes.

### Supplementary Table S24 | Complete code task map and actual reproduction scope

{{REPRODUCTION}}

The `recent_bridge` task can restore nine archived rounds and regenerate the declared deduplicated cache families. Merely extracting its stored files does not perform cache regeneration or repeat all fitting. The `dixit` task audits the frozen primary RNA summary and separately refits the later G2M development model; these are different operations. The legacy `paper` task's seven main figures, ten extended-data figures and nineteen supplementary-table references retain their v22 numbering.

### 12.3 New-document verification

The SI assembly checks source presence and hashes, embedded table completeness, table/section numbering, all thirty neural selection records and their final-fit durations, the fifteen-row/420-choice GSE220974 correspondence, and the unchanged main-text hash. It independently recalculates the formal risk reduction from the saved summary, the paired-seed formal interval, matched-RNA pooled gains from fold SSE and GSE220974 choice summaries from task-level outputs. These checks concern retained evidence and document consistency. They do not claim new model training, new biological experiments or a newly public software release.

Public-data accessions and bibliographic provenance are given in the main Data availability statement. Public repository registration, a verified software licence and a documented new-input interface remain distinct from delivery of this reviewer archive. No placeholder public DOI or unverified licence is inserted here.

## Supplementary reference correspondence

Numbers below match the main-text references; the omitted numbers remain in that list.

{{REFERENCES}}
