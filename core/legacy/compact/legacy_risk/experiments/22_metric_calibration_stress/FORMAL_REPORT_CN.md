# Experiment 22 v2 正式结果与审计报告

**状态：** 正式执行完成  
**执行时间（UTC）：** 2026-08-11  
**冻结清单 SHA-256：** `0bf032f88118523f05e11751c19cfee0201744b340b14357136090161041a756`  
**正式结果清单 SHA-256：** `db697b9a91fa6bbb7fd00f9aadb9e2be6d9554a1135f4d3e186ea517255ebc1a`

## 1. 结论

Experiment 22 在不改变 WitnessCell v14、Gate 19 或 Gate 21 科学参数的条件下完成。结果支持三个有边界的结论：

1. 在四数据集、12 个 split、654 个操作上，WitnessCell v14 对 canonical signal-weighted WMSE 与 weighted ΔR² 的优势不依赖未加权 MSE/Pearson 口径；相对 trainMean 和 baseControl 的双主指标均通过冻结检验。
2. Norman Gate 21 的冻结选择规则在 207/213 个具有 canonical 基因权重的完整案例中显著富集低 WMSE query；在匹配覆盖率下，单纯几何风险对照未通过。
3. 训练侧加权双指标 shadow gate 与原 Gate 19 在 10/12 个 split 上一致；两次不一致均为原 active 转为保守 inactive，没有 inactive 转 active。该分析是描述性反事实，不改写 v14。

本实验不支持以下越界声明：对 raw linearModel 的重新评分、15 个已发表方法的完整 weighted leaderboard、跨数据集选择性弃权、细胞分布预测，或全部 213 个 Gate 21 query 的 canonical WMSE 判决。

## 2. 冻结与执行完整性

- 预冻结 checklist 的 26 个 S/A/C/G/Q/F/D 项全部 `PASS`。
- 16 个源代码一致性、边界和负向单元测试全部通过。
- Phase M 在候选预测不可读的防火墙内先行执行；`candidate_predictions_read=false`。
- 冻结时正式输出目录不存在；76 个协议、代码、审计和合同对象随后一次性写入冻结清单。
- 正式结束后，28 个结果对象由 `FORMAL_RESULTS_MANIFEST.sha256` 逐项复核，全部 `OK`。
- 没有使用二元 DRF 有效性阈值；DRF 仅连续报告。

## 3. Phase M：候选盲的指标有效性

654/654 个正式操作可评估。以 technical duplicate 为正控制、mean prediction 为负控制时，canonical WMSE 的中位 DRF 为：Norman 0.967、Wessels 0.845、Schmidt 0.533、Replogle 0.878。对应 plain MSE 的中位 DRF 分别为 0.708、0.411、-0.372 和 0.251。

该结果表明，在此四数据集协议中，signal-weighted loss 通常比未加权 loss 更能分离技术重复上界与无信息均值基线。weighted ΔR² 与 WMSE 给出近似相同的 DRF 排序，这是其共享加权平方误差结构的预期结果，不应被解释为独立的第二份动态范围证据。

NIR 的 technical-duplicate、interpolated-duplicate 和 mean/control-negative 总体呈预期排序，但锁定来源未提供 NIR 的 DRF 认证，因此 NIR 只作描述性结果。control-referenced Pearson 的 baseControl 方向为常数，相关系数未定义，`valid_n=0`，不进入主判决。

## 4. Phase P：预测稳健性

四个候选/基线归档（v14、v13、trainMean、baseControl）均覆盖 654/654 个操作。20,000 次按数据集分层、以 condition 为 cluster 的 bootstrap 结果如下：

| 比较 | dataset-equal WMSE ratio | 95% CI | weighted ΔR² difference | 95% CI | 冻结判决 |
|---|---:|---:|---:|---:|---|
| v14 / trainMean | 0.683 | 0.626–0.743 | +0.227 | +0.195–+0.261 | 双主指标 PASS |
| v14 / baseControl | 0.398 | 0.354–0.443 | +0.649 | +0.590–+0.709 | 双主指标 PASS |
| v14 / v13 | 0.988 | 0.980–0.997 | +0.0076 | +0.0037–+0.0125 | 双主指标 PASS |

相对 trainMean，v14 的 dataset-equal WMSE 降低 31.7%；相对 baseControl 降低 60.2%。v14 相对 v13 的增量较小（约 1.17%），但冻结 cluster bootstrap 的双主指标均通过。四个数据集内 v14/trainMean WMSE ratio 均小于 1。

由于冻结前未找到 raw linearModel condition-mean predictions，且原始 R/runtime 状态不足以精确重建，该比较按协议判为 `PRED_LINEAR=NOT_ADJUDICATED`，不得由已发表汇总分替代。比较范围固定为 `MANDATORY_BASELINES_ONLY`，因此本实验不是 15 方法 weighted leaderboard。

## 5. Phase D：冻结选择性预测

Gate 21 原面板含 213 行、33 个 pair，冻结阈值为 `0.0923227147328771`，原始 accepted/rejected 身份保持 110/103。canonical truth-side weights 可为 207 行、32 个 pair 构造；缺失的 6 行均属于 `ELMSAN1+MAP2K3`。协议禁止填补，因此：

- 全 213 行判决：`GATE21_WMSE=NOT_ADJUDICATED_FULL_213`；
- 207 行预先规定 complete-case sensitivity：`GATE21_WMSE_COMPLETE_CASE=PASS`。

complete-case 主结果：

- pair-balanced coverage：0.489；
- accepted WMSE：0.0863；all WMSE：0.1807；rejected WMSE：0.2712；
- accepted/all ratio：0.477（95% cluster-bootstrap CI 0.317–0.707）；
- accepted/rejected ratio：0.318；
- frozen score 与 WMSE 的 Spearman：0.732；
- within-seed random-selection 单侧检验：`p=5.0×10^-5`；
- 五项冻结判据全部通过。

匹配覆盖率的 geometry control 得到 accepted/all ratio 1.168（95% CI 0.844–1.448；`p=0.701`），五项中仅覆盖率条件通过。由此，低损失富集来自冻结 Witness score，而非仅由几何覆盖率排序解释。

`decision_verdict.json` 还保留一个次级字段 `risk_wmse_spearman=NaN`：它对原 213 行数组直接计算，因 6 行缺失权重而未定义。完整案例主字段 `primary.score_loss_spearman=0.732` 有效；该次级字段不参与任何冻结判据，也未被事后替换。

## 6. Training-only shadow gate

shadow gate 的权重、拟合、置信下界和决定均只使用训练单扰动 LOO；未读取 validation 或 test outcome。结果为：

- 12 个 split 中 10 个与原 Gate 19 一致（83.3%）；
- 6 个原 active 仍 active；4 个原 inactive 仍 inactive；
- Schmidt seed 1/2 从原 active 转为 shadow inactive；
- inactive→active 翻转为 0。

因此，指标校准没有诱发更激进的启用；但因精确一致率为 83.3%，不能声称“门控决定完全不受指标选择影响”。正确表述是：加权双指标保留了全部原回退，并在两个 Schmidt split 上更加保守。

## 7. 执行期修正

Phase D 在生成任何 D 结果前两次触发仅与 complete-case 索引对齐有关的工程异常。冻结文件未修改：

- Amendment A1：将 207 行 complete-case DataFrame 重置为稠密位置索引；
- Amendment A2：在 geometry quantile 中应用冻结合同内同一个 `scoring_evaluable` mask，并取代 A1 作为正式恢复入口。

两次修正均未改变 query、预测、真值、权重、阈值、随机种子、重复次数或判据；Phase M/P 未重跑。A1/A2 文档、wrapper 和 receipt 均独立哈希保留。

## 8. 可用于论文的有边界表述

> In a candidate-blind metric-validity phase, signal-weighted losses more clearly separated technical duplicates from uninformative mean predictions than their unweighted counterparts. Across 654 operations in four datasets, WitnessCell reduced dataset-balanced weighted MSE by 31.7% relative to the training-mean baseline (95% cluster-bootstrap CI, 25.7–37.4%) while increasing weighted ΔR² by 0.227 (0.195–0.261). A previously frozen selective rule enriched low-loss queries in the 207-of-213-query complete-case Norman panel (accepted-to-all weighted-MSE ratio, 0.477; 95% CI, 0.317–0.707; within-seed random-selection P=5.0×10⁻⁵), whereas a matched-coverage geometry control did not. The full 213-query weighted analysis and the raw linear-model comparison were not adjudicated because the required canonical weights and raw predictions, respectively, were unavailable.

## 9. 最终状态矩阵

| 项目 | 状态 |
|---|---|
| SOURCE_PARITY | PASS |
| METRIC_VALIDITY | PASS_PHASE_M_EXECUTION_AND_SEAL |
| PRED_UNINFORMATIVE | PASS |
| ENDPOINT_COMPATIBILITY | STRONG |
| PRED_LINEAR | NOT_ADJUDICATED |
| GATE21_WMSE（全 213） | NOT_ADJUDICATED_FULL_213 |
| GATE21_WMSE（207 complete-case） | PASS |
| SHADOW_GATE | DESCRIPTIVE_COMPLETE；10/12 一致，0 个冒进翻转 |

