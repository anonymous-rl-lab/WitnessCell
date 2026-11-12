from pathlib import Path
import pandas as pd,numpy as np
W=Path(__file__).resolve().parents[1];R=W/'results'
pool=pd.read_csv(R/'MATCHED_REFERENCE_POOLED.csv');scores=pd.read_csv(R/'MATCHED_REFERENCE_SCORES.csv');modes=pd.read_csv(R/'MODE_DIAGNOSTICS.csv');att=pd.read_csv(R/'GLOBAL_ATTRIBUTION.csv');mod=pd.read_csv(R/'MODULE_GAIN_ATTRIBUTION.csv');han=pd.read_csv(R/'HAN_MODE_DIAGNOSTICS.csv')
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,row))+' |' for row in rows])
def f(x,n=4):return 'Inactive' if pd.isna(x) else f'{x:.{n}f}'
parts=[]
parts.append(r'''## Supplementary Note 24 | Global-reference sensitivity and bounded-intermediate interpretation

### 24.1 Three distinct response questions

Endpoint elimination tests whether a response can be expressed as a sum of freely varying endpoint contributions. It does not remove every function of single-effect strength. If $Y_{gh}=F(s_g+s_h)$, curvature of $F$ can generate a nonzero four-corner contrast even though the same rule applies to every gene pair. Such a global effect remains a legitimate part of the paper's operational partner response. A further comparison is required to ask what an added learner predicts beyond a specified global rule. Global epistasis models provide precedents for this distinction.[44]

The new analyses separate (i) complete versus endpoint-eliminated prediction; (ii) a single-effect-related global reference versus an additional learned correction; and (iii) a descriptive gene-associated pattern versus a molecular mechanism. The first two are prediction questions. The third cannot be settled by a rank, loading correlation or prediction gain alone. All original Horlbeck and Han results retain their reported baselines and evaluation status.

### 24.2 Matched global-reference refitting

For repeat $r$, extend the original standardized baseline to

$$B^\times_{gh,r}=B_{gh,r}+\gamma_r a_{g,r}a_{h,r}.$$

The single inputs $a$ retain the original across-gene standardization. The new coefficient is unpenalized and jointly fitted with all original baseline coefficients; sum-zero endpoint, slope and curvature deviations retain ridge $10/s$. This differs from merely adding a product to frozen predictions. The original global coefficients in $B$ are unpenalized, but a cross-product is otherwise represented indirectly through penalized endpoint-specific slopes. The extension gives that direction an explicit unpenalized coefficient without changing the remaining regularization.

Seven original configurations are retained: three low-support K562 folds, three high-support folds and Jurkat. Every inner subset refits $B^\times$ and the matrix learner. The original rank-one/spectral families, fractions $\{.3,.1,.03,.01\}$, penalties $\{.1,1\}$, convergence checks and zero-correction option are unchanged. Pooled internal complete SSE selects the correction, with the original stronger-penalty tie rule. All seven selections and outer predictions were persisted before the scoring pass. The analysis was motivated by already known results and is therefore retrospective baseline sensitivity. It adds no new source or untouched confirmation.

### Supplementary Table S54 | Matched global-reference sensitivity

**A. Pooled gains relative to each named baseline.** Original and expanded-reference percentages have different denominators. Control uses the same internal projected-error rule as Note 17.
''')
rows=[]
for group,label in [('K562_low','K562 low'),('K562_high','K562 high'),('Jurkat_high','Jurkat')]:
 a=pool[pool.group==group].set_index('method')
 for fam,ff in [('rank1','Rank one'),('spectral','Spectral')]:
  o=a.loc['native_'+fam];n=a.loc['global_'+fam];c=a.loc['global_'+fam+'_control']
  rows.append([label,ff,f(o.full_gain_pct),f(o.partner_gain_pct),f(n.full_gain_pct),f(n.partner_gain_pct),f(c.partner_gain_pct)])
parts.append(table(['Setting','Learner','Original full (%)','Original partner (%)','Expanded full (%)','Expanded partner (%)','Expanded controlled partner (%)'],rows))
parts.append('\n**B. Expanded baseline versus original baseline.** These gains evaluate the baseline alone, before matrix correction.\n')
rows=[]
for group in ['K562_low','K562_high','Jurkat_high']:
 a=pool[pool.group==group].set_index('method');o=a.loc['native_B'];b=a.loc['global_B']
 rows.append([group.replace('_',' '),f(o.full_sse,6),f(b.full_sse,6),f(o.partner_sse,6),f(b.partner_sse,6),f(100*(1-b.partner_sse/o.partner_sse))])
parts.append(table(['Setting','Original full SSE','Expanded full SSE','Original partner SSE','Expanded partner SSE','Baseline partner gain (%)'],rows))
parts.append('\n**C. All seven configurations after global-reference refitting.** Inner selection remains based on complete SSE; a selected correction may increase outer partner error.\n')
rows=[]
for c in modes.itertuples():
 a=scores[(scores.configuration==c.configuration)&(scores.repeat=='pooled')].set_index('method')
 rows.append([c.configuration,f(a.loc['global_rank1','partner_gain_pct']),f(a.loc['global_spectral','partner_gain_pct']),c.global_refit_rank1_selected,c.global_refit_spectral_selected])
parts.append(table(['Configuration','Rank-one partner (%)','Spectral partner (%)','Rank-one selection','Spectral selection'],rows))
parts.append(r'''
The full per-repeat losses, absolute denominators, original outputs, control outputs and internal selection grids are supplied in the revision supplement. Low-support nonpositive outcomes are retained. The two K562 support levels share query identities within folds; the seven configurations are not seven independent biological studies.

### 24.3 Product-direction diagnostics and exact gain attribution

For internal held-combination fold $k$ and repeat $r$, let $t_{k,r}=P_k(s_gs_h)$ and $e_{k,r}=P_k(Y-B_k)$. A scalar diagnostic uses

$$\widehat\beta_r=\frac{\sum_k t_{k,r}^{\mathsf T}e_{k,r}}{\sum_k\|t_{k,r}\|^2},$$

returning zero for energy at most $10^{-14}$. It retains the original baseline and does not refit a matrix learner. Both repeat-mean singles (the same input vector for both repeats) and each repeat's own singles are reported. The coefficient is not constrained to a sign. Outer scoring uses the fixed projected query template, with no query outcome in fitting.

For retrospective attribution, write the original projected correction as $d=d_t+d_\perp$, where $d_t$ is its orthogonal projection onto the single-product template separately for each repeat. For original baseline error $e=P_Q(Y-B)$,

$$\|e\|^2-\|e-d\|^2=[2\langle e,d_t\rangle-\|d_t\|^2]+[2\langle e,d_\perp\rangle-\|d_\perp\|^2].$$

The geometric split uses query identities, singles and predictions; the attributed error reductions use query outcomes. Negative contributions are retained. This is a descriptive decomposition of saved predictions, not a training-side forecast. Projection onto one product direction does not remove every possible global nonlinear response.

### Supplementary Table S55 | Original-model modes and single-product attribution

**A. Mode identities.** The fitted mode is the eigenvector of the selected archived model matrix. The residual mode is the largest-absolute-eigenvalue eigenvector after projecting repeat-mean training residuals and placing them symmetrically on the observed graph, with unobserved entries and the diagonal set to zero. Correlations are absolute centered Pearson correlations; no mode is attributed to an inactive zero-selected model.
''')
rows=[]
for q in modes.itertuples():rows.append([q.configuration,f(q.fitted_mode_single_abs_r),f(q.residual_mode_single_abs_r),f(q.fitted_residual_mode_abs_r),q.selected])
parts.append(table(['Configuration','Fitted mode vs singles','Residual mode vs singles','Fitted vs residual mode','Original selection'],rows))
parts.append('\n**B. Query gains on the original baseline.** The last two columns use the original complete partner-error denominator and sum to the original rank-one gain.\n')
rows=[]
for q in att[att.configuration.str.contains('high')].itertuples():rows.append([q.configuration,'Mean' if q.template=='mean_singles' else 'Repeat-specific',f(q.product_gain_pct),f(q.original_rank1_gain_pct),f(q.global_direction_gain_pp),f(q.orthogonal_gain_pp)])
parts.append(table(['Configuration','Single inputs','Scalar product gain (%)','Original rank-one gain (%)','Along product (pp)','Orthogonal remainder (pp)'],rows))
parts.append(r'''
The Jurkat scalar product recovers most of the original increment. The matched reference in Table S54 is the stronger check because it jointly refits baseline and matrix correction. A low loading correlation alone would not exclude a nonlinear relationship with single effects. No empirical boundary density or causal bottleneck is estimated from the scalar coefficient.

### 24.4 Descriptive module-associated structure in K562

The original high-support fitted mode is stable in its broad loading pattern across the three folds. Its positive orientation is chosen so that the designated mitochondrial group has the larger mean loading; a rank-one eigenvector's sign itself is arbitrary. This group, specified in the supplied retrospective proposal and fixed before this replay, contains ATP5J2, MRPL3, MRPL13, MRPL18, MRPL27, MRPL32, MRPL33, MRPL43, MRPS14, MRPS27, MRPS30, MRPS35, NDUFA8, PET112, PTCD1 and SSBP1. It is a descriptive annotation, not a prospectively tested enrichment hypothesis.

In fold 0, high positive fitted loadings include MRPL33, MRPL3, MRPL43, PET112, MRPS35 and NDUFA8; high negative loadings include ACTL6A, PSMD12, CACTIN, MAX and SNRPG. All gene loadings for every configuration are retained in `GENE_LOADINGS.csv`. No gene is removed from prediction or scoring. The result motivates a module-associated interpretation, consistent with previous work on modular interaction patterns,[45] without assigning a physical mechanism to each inferred relationship.

### Supplementary Table S56 | Query coverage and gain associated with the mitochondrial set

For each model, the numerator sums per-pair projected error reductions, including negative contributions, over pairs containing at least one listed gene. The denominator is that model's net projected SSE reduction on all queries. These shares are descriptive; endpoint sharing prevents an independent-pair interpretation.
''')
rows=[]
for c in ['K562_high_fold0','K562_high_fold1','K562_high_fold2']:
 a=mod[(mod.configuration==c)&(mod.pair_class=='any')].set_index('model');o=a.loc['original'];g=a.loc['global_refit']
 rows.append([c,int(o.pairs),int(o.total_pairs),f(100*o.query_share,2),f(100*o.net_gain_share,2),f(100*g.net_gain_share,2)])
parts.append(table(['Configuration','Pairs with set member','All query pairs','Query share (%)','Original net gain share (%)','Expanded-reference net gain share (%)'],rows))
parts.append(r'''
### 24.5 Han model modes beyond the tested nonlinear reference

The diagnostic reads the saved final matrices for all sixty native-analysis batches and both correction families. The largest-absolute-eigenvalue component defines the leading mode. A numerical eigenvalue exceeds $10^{-8}\max(1,\max_j|\lambda_j|)$ to count toward fitted rank, matching the original numerical definition. A zero-selected candidate remains unused even if its diagnostic fitted matrix is nonzero. Correlations use mean singles and the same centered Pearson convention as Table S55.

### Supplementary Table S57 | Saved Han mode diagnostics

Ranges span fifteen fits per readout/baseline/family and describe repeated partitions of one study. Correlation ranges exclude zero-selected fits. The complete top-ten loadings of up to three leading modes per selected fit are supplied in `HAN_MODE_LOADINGS.csv`.
''')
rows=[]
for (panel,base,fam),a in han.groupby(['panel','baseline','family']):
 active=a[a.used_rank>0];rows.append(['Reconstructed' if panel=='full_raw' else 'Published',base.title(),'Rank one' if fam=='rank1' else 'Spectral',f'{len(active)}/15',f'{a.used_rank.min()}–{a.used_rank.max()}',('Inactive' if active.empty else f'{active.leading_single_abs_r.min():.3f}–{active.leading_single_abs_r.max():.3f}')])
parts.append(table(['Readout','Baseline','Learner','Nonzero selections','Used numerical rank','Absolute mode correlation with singles'],rows))
parts.append(r'''
The secondary spectral gains over the saturation reference remain 4.8245% and 7.0357% in the reconstructed and published readouts (Table S45). These results establish an increment beyond that fitted reference. Neither the rank range nor the loading correlations identify how many biological pathways generated the response. The saved native rank-one modes also have low single-effect correlations; loss of a gain after baseline adjustment should therefore not be equated with proof that its fitted mode was proportional to singles. No additional Han model was fitted for this diagnostic.

### 24.6 A bounded-intermediate explanation under explicit assumptions

Let a latent process $k$ have wild-type level $b_k$ within $[L_k,U_k]$, with headrooms $h_k^-=b_k-L_k$ and $h_k^+=U_k-b_k$. Define its change under input $s$ as

$$f_k(s)=\operatorname{clip}(s,-h_k^-,h_k^+).$$

Assume independent latent processes, additive intervention inputs $\delta_{gk}$ and a linear measured readout. Then

$$Y_g-Y_0=\sum_k w_k f_k(\delta_{gk}),\qquad
\varepsilon_{gh}=\sum_k w_k[f_k(\delta_{gk}+\delta_{hk})-f_k(\delta_{gk})-f_k(\delta_{hk})].$$

If all three arguments remain in the linear range, that process contributes zero interaction. For $d,c\geq0$, its contribution before readout weighting is exactly

$$\varepsilon_k=-[\min(d,c,h_k^+,d+c-h_k^+)]_+.$$

This follows by considering whether neither, one or both singles reach the upper bound. The negative-input formula follows by reversing signs and exchanging headrooms. Since clipping is monotone and 1-Lipschitz with $f_k(0)=0$,

$$|\varepsilon_k|\leq\min(|d|,|c|).$$

For example, when $d\geq0$, the increment $f_k(c+d)-f_k(c)$ and $f_k(d)$ both lie in $[0,d]$, bounding their difference by $d$; exchange the inputs for the other bound. This inequality concerns latent input shifts, not arbitrary measured single-effect magnitudes. With $f(s)=\operatorname{clip}(s,-1,2)$, $d=3$ and $c=-10$, the interaction has magnitude 2 while the smaller observed single magnitude is 1.

For nonnegative weights and nonpositive inputs to every process, the aggregate single changes $s_g,s_h$ obey the conditional envelope

$$s_g+s_h\leq Y_{gh}-Y_0\leq\min(s_g,s_h).$$

Each clipped negative input is monotone and subadditive in its nonnegative loss magnitude, giving both inequalities after weighting and summation. Observing negative total single responses does not establish the required componentwise signs.

A cooperative recurrent extension has loss state $q=\min(h,p+Wq)$, where $p\geq0$, $W\geq0$ and $\|W\|_\infty<1$. Contraction gives a unique fixed point. Iteration from zero establishes monotonicity and $q(p_1+p_2)\leq q(p_1)+q(p_2)$: use $\min(h,x+y)\leq\min(h,x)+\min(h,y)$ at each step. A nonnegative linear readout inherits the same loss envelope. This is a theorem for the declared model class. Agreement with an envelope does not identify that class, and disagreement can reflect a failure of any of its assumptions, including the readout or input model.

Heterogeneous boundaries can smooth a hard response. For a random upper headroom $K$, $F(s)=\mathbb E[\min(s,K)]$ has $F'(s)=\Pr(K>s)$ and $F''(s)=-p_K(s)$ at differentiable interior points. A local expansion yields a product term $F''(0)dc$ when that expansion is valid. Latent scales and readout weights remain unspecified, so an empirical product coefficient is a composite curvature parameter, not an identified biological boundary density.

### 24.7 Why response rank does not count biological intermediates

A local second-order expansion can yield $G\approx D\operatorname{diag}(w_k\kappa_k)D^{\mathsf T}$ and hence rank at most the number of retained curved latent directions. This statement requires the approximation and representation assumptions.[46] It is not a general rank identity for bounded response functions.

For a counterexample, take one intermediate $f(x)=\min(x,1)$ and two groups of three distinct interventions with strengths $(.2,.5,.8)$ and $(.3,.6,.9)$. Their cross-group interaction matrix is

$$G=\begin{pmatrix}0&0&-.1\\0&-.1&-.4\\-.1&-.4&-.7\end{pmatrix}.$$

It has determinant $0.001$ and rank 3. Removing additive row and column contributions gives rank 2. Thus, one bounded intermediate can already generate more than one partner direction. Conversely, multiple latent contributions can be redundant or cancel. Neither a fitted numerical rank nor one observed growth readout identifies latent-process count, loading uniqueness or molecular identity.

### 24.8 Endpoint-equivalent module descriptions

Let $m_g$ indicate a gene set. Suppose a block description assigns interaction $a$ inside the set, $b$ across it and $c$ outside it. Then

$$\varepsilon_{gh}=c+(b-c)(m_g+m_h)+(a-2b+c)m_gm_h.$$

The constant and endpoint terms disappear under $P_Q$, leaving

$$P_Q\varepsilon=(a-2b+c)P_Q(m_gm_h).$$

When the last direction is nonzero, the projected contrast concerns the block difference $a-2b+c$. Stronger within-set buffering and weaker cross-set buffering can generate the same projected sign pattern as other absolute interaction descriptions. A negative projected entry therefore does not establish an aggravating interaction relative to the biological additive reference. This identity connects endpoint elimination to module interpretation without changing the prediction loss or asserting a causal mechanism.

### 24.9 Verification and empirical scope

The supplied `check_theory.py` checks 100,000 independently drawn clipping configurations, 1,000 eight-node cooperative networks, the block-contrast identity and the explicit rank counterexample. These finite checks support implementation consistency; the mathematical conclusions rest on the stated arguments. No empirical three-standard-error classification, pathway-count estimate, inferred network-sign map, clipping-based predictor or new repeatability gate is claimed. The two-repeat measurements and shared controls require their own error model before an envelope-frequency test could be interpreted statistically. The new diagnostics retain the empirical gains and delimit the biological explanations they can support.
''')
(W/'note24.md').write_text('\n\n'.join(parts)+'\n')
print('note24 written')
