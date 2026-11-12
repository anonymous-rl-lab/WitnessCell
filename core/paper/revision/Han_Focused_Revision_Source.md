# WitnessCell Han 定向论证与论文修订稿

本稿集中补强一个结论：**在Han单扰动相加的基线上，保持已有伙伴预测器不变，只用Han校准组合决定修正强度，两个SF12来源配置在留出组合上取得了正向增量。**

它延伸论文已有的“端点与伙伴响应分开控制”主线。原有预测榜单、因果信息、Jurkat功能增量和RNA伙伴结果继续承担各自的结论。本稿提供可直接并入论文的Results、Methods、Discussion与SI内容；未改动v34原文件，也未启动新实验。

## 论证只围绕三个问题

| 需要回答的问题 | 直接对照 | Han全量结果说明什么 |
|---|---|---|
| 训练数据选出的强度有用吗？ | 目标侧校准 vs 同一伙伴修正原样使用 | HeLa和293T来源SF12均降低留出伙伴误差 |
| 是不是把伙伴全部关掉就够了？ | 目标侧校准 vs 完全关闭伙伴 | 两个配置仍优于关闭，保留了小幅有用信息 |
| 完整预测是否实际变好？ | 校准后的完整预测 vs Han单扰动相加 | 完整误差分别降低0.187%和0.465%，端点贡献保持不变 |

HeLa的12.84%表示相对过强修正的误差降低，不是相对简单相加的最终净收益。正文同时交代“修正损伤减少”和“基础预测之上的净增量”，避免用较大的前一个数字替代后一个数字。

## Results英文正文

### Target-side calibration improves held-out prediction with frozen partner functions

We next asked whether measurements from a target screen could determine how strongly to use an existing partner predictor. In a full-screen extension of the Han analysis, we fixed the baseline to the sum of the measured single-perturbation responses and supplied partner corrections from frozen SF12 and V4 functions trained on the three Shen cell backgrounds. Only the correction strength was estimated from Han calibration combinations. The reconstructed panel contained 20,359 eligible pairs; whole-pair splits separated calibration and query outcomes, and predictions were saved before query scoring. All six source–model configurations were retained.

Target-side calibration improved the HeLa- and 293T-derived SF12 corrections relative to both applying the correction unchanged and removing it entirely. On the reconstructed readout, partner error fell by 12.84% and 0.054%, respectively, relative to the unchanged correction. Relative to removing the correction, the reductions were 0.262% and 0.651%, yielding complete-response error reductions of 0.187% and 0.465% over the additive baseline. The endpoint contribution was unchanged. Thus, attenuation removed substantial damage from the HeLa correction while retaining a small useful increment, whereas the already beneficial 293T correction received a smaller further improvement.

Both SF12 configurations improved over the unchanged and removed corrections in each of five pair partitions, and their pooled gains over the additive baseline were positive in both reconstructed repeats. The published-phenotype panel also showed positive partner and complete-response increments. The other configurations did not show the same consistent benefit. These results provide evidence that target-side calibration can retain useful information in a frozen partner predictor while reducing excessive correction, with modest net gains beyond the additive baseline.

## Methods英文补充

### Target-side calibration of frozen shared partner functions

The full-screen extension used a target additive baseline, B, computed from the permitted Han single-perturbation responses. Existing SF12 and V4 shared functions from HeLa, A549 and 293T supplied the correction D. Their fitted coefficients, network weights, source-selected amplitudes and source feature transformations were retained. A source-selected zero correction remained zero. This interface extends the shared functions to eligible Han identities; it does not transfer the complete gene-indexed source baseline. Rank-one and spectral matrices were evaluated separately on the common-identity panel.

For each of five fixed pair partitions, three outer folds separated calibration combinations from query combinations. Both reconstructed repeats followed the same unordered-pair identity. Within each calibration pool, two blocks defined the endpoint-removal projections used to estimate a single coefficient for each source–model configuration. For block k, let r_k = P_k(Y_k − B_k) and d_k = P_kD_k. We used

\[
\widehat\alpha=\operatorname{clip}_{[0,1]}
\frac{\sum_k\langle r_k,d_k\rangle}
{\sum_k\|d_k\|_F^2},
\]

with zero returned when correction energy was at most 10⁻¹⁴. The coefficients were shared across the reconstructed repeats. Query predictions were B_Q + α̂P_QD_Q. The unchanged and removed-correction comparators used α = 1 and α = 0 on the same query panel. P_Q was determined by query pair identities, without query outcomes; predictions therefore depend on the specified query panel. Since only the partner component was added, its scalar shrinkage and whole-vector shrinkage are identical in this adapter.

Coefficients and prediction hashes were saved before each fold’s query outcomes were read for scoring. Squared errors were pooled over folds before calculating relative reductions. The published-phenotype panel followed the same procedure with its corresponding permitted single-response inputs. The base predictors were not retrained; estimating α̂ did use Han calibration outcomes. This extension followed earlier access to Han results and is a frozen reassessment, not an untouched prospective confirmation.

## Discussion英文替换段

The Han extension connects component control to a concrete use: retaining an existing partner predictor while estimating its strength from measurements available in a target screen. Two SF12 configurations benefited both over using the full correction and over removing it, providing evidence that useful partner information can survive attenuation. Their modest complete-response gains distinguish this result from merely repairing a poor uncalibrated predictor. Together with the Shen analysis, which separates endpoint retention from whole-correction shrinkage, these results support using the response decomposition to control a learned correction as well as to evaluate it.

## SI补充

### Analysis identity and complete model reporting

The reconstructed panel comprised 20,359 qualified pairs among 206 active endpoints, with two reconstructed repeats. The published panel comprised 20,791 pairs among 205 active endpoints. These panels differ in response construction, permitted single-response inputs and coverage; they do not isolate the effect of response definition alone. Raw-panel eligibility followed the fixed count, plasmid, construct-count and orientation-support rules, without filtering by model performance. ADRBK1 lacked a stored embedding and used the existing global-mean rule; complete endpoint exclusion was examined separately.

All three source backgrounds and both portable model families were retained in the full-screen analysis. The frozen protocol designated the common-panel rank-one comparison as primary and the other model families as secondary. The positive full-screen SF12 findings therefore do not replace the common-panel primary result. Source configurations are model variants applied to one target study, not independent target-study replications.

**Table H1. Full-screen target-side calibration.** Positive percentages indicate error reduction. “Unchanged” uses the same frozen direction with α = 1. “Removed” uses α = 0 and is identical to the additive baseline in this adapter.

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

### Partition and repeat results

For reconstructed outcomes, the HeLa SF12 partner gain over the additive baseline ranged from 0.176070% to 0.351213% across the five partitions; the corresponding range for 293T SF12 was 0.592577% to 0.691432%. Gains over the unchanged correction were also positive in every partition. Estimated strengths across folds ranged from 0.078561 to 0.166958 for HeLa and from 0.582805 to 0.868120 for 293T. The five overlapping partitions describe sensitivity to pair allocation and are not independent biological replicates.

**Table H2. Reconstructed-repeat gains over the additive baseline.** Each entry pools the locked predictions across the five partitions.

| Source and model | Repeat | Partner reduction (%) | Complete reduction (%) |
|---|---|---:|---:|
| HeLa SF12 | 1 | 0.442913 | 0.318944 |
| HeLa SF12 | 2 | 0.080751 | 0.057269 |
| 293T SF12 | 1 | 0.741102 | 0.533672 |
| 293T SF12 | 2 | 0.561746 | 0.398393 |

The two reconstructed repeats share initial and plasmid controls. The published phenotype derives from the same study. These results assess repeat and processing sensitivity within Han rather than adding independent source confirmations. No population-level confidence claim is inferred from the five overlapping partitions.

### Scope of the control comparison

The full-screen adapter fixes the endpoint contribution to that of B and supplies only P_QD_Q. Consequently, changing α preserves that contribution and changes complete and partner squared-error sums by the same absolute amount. Their percentage changes use different denominators and are not independent successes. This experiment tests whether the target-estimated strength is useful; the distinction between component control and shrinking a correction that contains both components is evaluated by the original Shen analyses.

The common-panel rank-one primary analysis retained a different result: on reconstructed outcomes, target-side control reduced pooled partner error by 6.670017% relative to the original correction, but increased it by 0.996058% relative to removing the correction. That outcome remains part of the evidence. No success of the secondary SF12 adapter is transferred to the rank-one model.

### Access history and identity checks

Earlier Han aggregate outcomes and 88 reconstructed pairs had been scored before this extension. The new protocol, source files and run code were locked before its reassessment scores, and each query prediction was locked before scoring, but this does not reset the prior access history. Excluding previously scored pairs from both calibration and query pools retained positive SF12 gains: reconstructed partner reductions were 0.268364% for HeLa and 0.659370% for 293T. Excluding the missing-embedding endpoint also retained the gains. These are post-score, identity-defined sensitivity checks under the same algorithm, not new confirmatory experiments.

## 定向并入方式

1. **Results：** 以上三段承接现有Shen分量控制结果。正文交代训练侧调节、两个对照和最终净收益；不再展开多个调试方向。
2. **Methods：** 增补全量共享函数适配与Han强度估计，明确它与v34固定源强度迁移的区别。旧实验方法保留在SI，不能把新分析追溯写成旧协议。
3. **Discussion：** 用上述一段解释新增用途；原有完整预测、因果信息与伙伴学习贡献不动。
4. **Fig.4：** 保留既有功能结果。Han位置只呈现“原样使用、全部关闭、训练侧调节”三个对照；所有来源及读出列入SI。图尚未在本稿中重绘。
5. **SI：** H1和H2是本定向稿的局部表号，正式合稿时连续编号。完整既有配置、原88对结果、访问历史和敏感性分析保留。

## 已核对的来源

- `han2017_frozen_reassessment_v3/PROTOCOL.json`：主次分析身份、六种全量配置、固定校准公式与访问历史。
- `han2017_frozen_reassessment_v3/code/run.py`：校准只接收校准结局，查询预测锁定后才进入评分。
- `han2017_frozen_reassessment_v3/results/POOLED_SCORES.csv`：本文主百分比可由保存SSE直接复算。
- `han2017_frozen_reassessment_v3/results/FULL_HAN_STABILITY.csv`：五套划分和控制系数范围。
- `han2017_frozen_reassessment_v3/results/PRIMARY_COMMON_RANK1.csv`：原共同身份面板主分析。
- `han2017_frozen_reassessment_v3/audit/FINAL_CHECKS.json`：归档的查询结局隔离、独立QR评分与冻结文件检查。本稿未将这些实现检查计作生物学证据。

本轮仅核对已有代码与结果、复算保存SSE并撰写定向文本；未重训模型，未改变评分口径，未生成新的实验结果。
