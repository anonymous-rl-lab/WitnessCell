# WitnessCell Gate 21：冻结选择性预测正式报告

## 结论

**PASS_FROZEN_SELECTIVE_PREDICTION_GATE**

WitnessCell 的训练侧 estimated Witness Risk 在 calibration 从未出现的
目标组合上支持一个冻结的 accept/abstain 决策。它能够在接近 50% 覆盖率时
选择出误差显著更低的一组目标，但不支持把高风险目标自动路由到 geometry-only
预测器。

## 冻结设计

- calibration：Gate 07 工程 seeds 0–4，共 130 rows、80 个目标组合；
- formal test：Gate 07 seeds 100–129 中从未在 calibration 出现的目标身份，
  共 213 rows、33 个组合；pair identity overlap 为 0；
- score：`estimated_witness_risk`，只由训练双扰动和 target-safe 单扰动描述估计；
- primary threshold：`0.0923227147328771`；
- threshold 来源：calibration pair-balanced 50% risk quantile；
- test outcome、test reliability 和 oracle risk 在协议、代码、阈值与 manifest
  哈希冻结前均不可读取；
- 冻结协议 SHA-256：
  `099e84619bd7f235845de5297138520a6b0f06558fa8db1f3fc58b2b2c622137`。

## 主结果

| Quantity | Frozen result |
|---|---:|
| Pair-balanced coverage | 49.48% |
| Accepted / rejected rows | 110 / 103 |
| All-target MSE | 0.12385 |
| Accepted MSE | 0.06617 |
| Rejected MSE | 0.18035 |
| Accepted / all MSE ratio | 0.5343 |
| Accepted / rejected MSE ratio | 0.3669 |
| Accepted/all pair-cluster bootstrap 95% CI | 0.3646–0.7591 |
| Within-seed random-selection one-sided p | 0.000050 |
| Test risk–error Spearman | 0.6357 |
| Risk–coverage curve Spearman | 1.000 |

相对于在全部目标上无条件输出，冻结 gate 在保留约一半 query 时使选择性 MSE
下降 **46.57%**。五项预声明判据——非退化覆盖率、至少 20% 实用效应、
pair-cluster CI 排除 1、胜过随机选择以及 accepted/rejected 分离——全部通过。

## 风险–覆盖曲线

| Calibration target coverage | Actual pair-balanced coverage | Accepted/all MSE ratio |
|---:|---:|---:|
| 30% | 24.04% | 0.5209 |
| 50% | 49.48% | 0.5343 |
| 70% | 84.14% | 0.8569 |
| 90% | 95.56% | 0.9670 |
| 100% | 100% | 1.0000 |

阈值均由 calibration 冻结。较高覆盖率的实际 coverage 偏离 calibration target，
因此横轴必须报告实际 coverage，不能把 target coverage 冒充 test coverage。

## 稳健性与负对照

- 28 个同时含 accept/reject 的 formal seed 中，28/28 的 accepted MSE 低于
  该 seed 全目标 MSE；median accepted/all ratio 为 0.667；
- pair-median decision 的 accepted/all ratio 为 0.561；
- leave-one-pair-out ratio 范围为 0.502–0.630；
- outcome-blind 等覆盖率 geometry-risk 对照的 accepted/all ratio 为 1.152，
  即选择后反而恶化 15.2%；
- calibration isotonic mapping 的 test mean-MSE 预测为 0.1050，实际为 0.1239，
  绝对风险低估约 15.2%。因此排序与选择成立，但不声称精确概率校准；
- accepted 用 Witness、rejected 自动改用 geometry-only 的整体 MSE 为
  always-Witness 的 1.083 倍。fallback routing 未通过。

## 可支持的论文主张

> A risk threshold frozen on five engineering splits accepted 49.5% of
> calibration-unseen Norman target queries. Their pair-balanced MSE was 46.6%
> lower than that over the full target panel (accepted/all ratio, 0.534;
> pair-cluster bootstrap 95% CI, 0.365–0.759; within-seed random-selection
> p=5.0×10−5). The ordering was stable to pair deletion and failed under a
> geometry-only risk control.

## 不能支持的主张

1. 这不是 prospective biological validation；正式 rows 在历史实验中已经存在，
   本门是新冻结的 retrospective heldout reanalysis。
2. 这不证明 Witness Risk 是模型无关 uncertainty score。
3. 这不提供任意 coverage 下的严格风险保证；绝对 isotonic risk 仍有偏差。
4. 这不支持自动 geometry fallback。当前已验证的 decision 是 accept/abstain，
   不是多模型 router。
5. 结论限定于 Norman development-panel pseudobulk residual/full-effect MSE；
   不自动外推到 scPerturBench 四数据集或分布级指标。
