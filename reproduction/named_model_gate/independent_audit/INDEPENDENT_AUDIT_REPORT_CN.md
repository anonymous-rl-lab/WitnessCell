# GSE146194 具名 AIVC 正式门 - 独立审计报告

## 审计结论

**PASS_INDEPENDENT_AUDIT。** 本次审计从 18 个原始具名模型预测档案与 3 个
CPU truth anchors 独立重建六种策略，复算 756 行逐靶指标。与 AutoDL 包内正式
结果相比，四项指标的最大绝对差为 `1.11e-16`，等于浮点机器误差。

这道门支持论文 B 的建设性中心结论：**target-conditioned double-perturbation
witnesses 可以形成一个独立的组合扰动虚拟细胞，并在目标干预几何被见证时，
提高 strict-seen1 组合响应预测。**

## 1. 冻结与完整性

- 科学协议 SHA-256：`459bb2b1663a78fefb4231d101941d2fe54a41c6ad0ceca1198a6d1dbe9e66b0`
- 训练预算 SHA-256：`6d64b639a60405df0cd2fbd0cc2f6869df256d4b22d95552abaeb60e8ebda7c7`
- 数据拆分：3 个 independent cell-split seeds（41, 47, 59）
- 模型随机种子：每个数据拆分 3 个（0, 1, 2）
- 具名模型：GEARS 9 fits + CPA 9 fits，共 18 fits
- 最终靶对：14 个 strict-seen1 pair；同一靶对跨 fit 重复，不能把 9 fits 当成 9 个生物重复
- 六策略：Saturation、GEARS、CPA、GEARS+Witness、CPA+Witness、WitnessCell

原始结果包：`GSE146194_named_AIVC_formal_results.zip`，9,855,248 bytes，
SHA-256 `316b3e360ee368d37ac021336d235504c1c94e34963807526028708c9c9a6eb8`。

## 2. 泄漏与工程审计

- 9/9 GEARS fits 完成冻结的 40 epochs；训练审计均记录
  `target_expression_used=false`，最终查询行来自 source controls。
- 9/9 CPA fits 完成冻结的 80 epochs 与 16 次验证检查；checkpoint 仅按 validation
  `cpa_metric` 选择，且只能从 epoch 65/70/75/80 的 post-warm 候选中选择。
- 9/9 CPA 反事实查询审计均记录 `test_expression_used=false`、
  `test_condition_cell_count_used=false`、registry 已重编码且预测未塌缩。
- 1000 个评价基因由 source-half-A 的 12 个训练双扰动选择；14 个 final target
  outcomes 未参与基因选择、模型训练、checkpoint 选择、gamma 或 fusion weight 选择。
- source-half-A 的 5 个 validation pairs 只用于 checkpoint/Witness hyperparameter；
  destination-half-B 的 5 个 calibration pairs 只用于闭式 gamma 与 model-specific lambda；
  14 个 destination-half-B final pairs 只用于评分。

未发现目标表达泄漏、测试集选 epoch、object/pickle NPZ、档案 schema 漂移或缺失 fit。

## 3. 六策略绝对表现

数值先在 9 个重复中按 pair 聚合，再以 14 个生物学 pair 为统计单位：

| 策略 | Residual MSE | Full-effect cosine |
|---|---:|---:|
| Saturation | 0.020748 | 0.834122 |
| GEARS | 0.042356 | 0.720165 |
| CPA | 0.197577 | 0.084906 |
| GEARS+Witness | **0.011710** | **0.906217** |
| CPA+Witness | 0.011876 | 0.905139 |
| WitnessCell | 0.011933 | 0.904838 |

在这个稀疏外部 pathway panel 和冻结训练预算下，原生 GEARS/CPA 都不应被当作
强于 saturation 的基线。最能证明 Witness 方向本身价值的比较是 WitnessCell vs
saturation，而不是对较弱 CPA 报 94% 的夸张头条。

## 4. 预冻结四个主比较

| Candidate vs comparator | Residual MSE 改善 | Pair bootstrap 95% CI | Pair wins | Exact p | Full cosine delta |
|---|---:|---:|---:|---:|---:|
| GEARS+Witness vs GEARS | 72.35% | 70.19%-74.35% | 14/14 | 6.10e-5 | +0.1861 |
| CPA+Witness vs CPA | 93.99% | 93.50%-94.46% | 14/14 | 6.10e-5 | +0.8202 |
| WitnessCell vs GEARS | 71.83% | 69.40%-74.14% | 14/14 | 6.10e-5 | +0.1847 |
| WitnessCell vs CPA | 93.96% | 93.45%-94.45% | 14/14 | 6.10e-5 | +0.8199 |

四个比较均在 3/3 data seeds 与 9/9 runs 上为正，并通过全部冻结条件。

## 5. 保守的论文级比较

| Candidate vs comparator | Residual MSE 改善 | Pair bootstrap 95% CI | Pair wins | Exact p |
|---|---:|---:|---:|---:|
| **WitnessCell vs Saturation** | **42.49%** | **40.31%-44.91%** | **14/14** | **6.10e-5** |
| GEARS+Witness vs Saturation | 43.56% | 41.76%-45.56% | 14/14 | 6.10e-5 |
| CPA+Witness vs Saturation | 42.76% | 40.60%-45.11% | 14/14 | 6.10e-5 |
| GEARS+Witness vs WitnessCell | 1.87% | 0.45%-3.23% | 10/14 | 0.0212 |
| CPA+Witness vs WitnessCell | 0.48% | -0.11%-1.08% | 8/14 | 0.0870 |

最稳妥的头条是：**独立 WitnessCell 相对无交互的 saturation 构造对照将残差
MSE 降低 42.49%，且 14/14 靶对一致获益。** GEARS 给 WitnessCell 贡献了统计上
可见但很小的 1.87% 额外信息；CPA 的额外贡献不稳健。

## 6. 融合权重揭示的机制

- `lambda_GEARS` 平均 0.9183，范围 0.9096-0.9354。
- `lambda_CPA` 平均 0.9790，范围 0.9725-0.9892。
- `gamma_Witness` 跨三数据拆分平均 0.6736。

因此 `GEARS+Witness` 与 `CPA+Witness` 的优异表现主要来自 Witness。数据允许
如下机制性表述：

1. Witness geometry 本身已经构成高精度的组合 AIVC；
2. GEARS 残差中存在少量、可校准地保留的互补信息；
3. CPA 在本门中基本被 Witness 替代，不能写成“CPA 被普适修复”；
4. 目标条件 calibration 是对信息源做凸组合，而不是重新训练具名模型。

## 7. 可用主张与禁止升级

### 可用

> On a frozen strict-seen1 pathway panel, target-conditioned witness geometry
> forms an independent combinatorial virtual cell that reduces interaction-
> residual error relative to saturation and outperforms GEARS and CPA trained
> under the same sparse information budget.

### 暂不可用

- 不可声称所有 AIVC、所有数据集或所有基因对均被普适修复。
- 不可声称已解决 seen0 冷启动。
- 不可把 18 fits 写成 18 个独立生物学实验。
- 不可把本结果直接升级为主动 acquisition 或实验预算节省结论。
- 不可隐去全部 14 targets 共享 FDPS 这一结构限制。

## 8. 审计复现

```bash
python independent_audit/independent_audit.py \
  --experiment-root . \
  --out independent_audit
```

输出包括 `independent_target_rows.csv`、`strategy_metrics.csv`、
`frozen_comparisons.csv`、`conservative_comparisons.csv`、`run_metadata.csv` 与
`audit_summary.json`。

