# Wessels 固定端点的见证库数量实验

协议：`WitnessCell-step1-followup-v1-20261017`。

本目录只改变原始训练组合中用于残差回归的见证行。每个原 seed 只从原始训练 singles 恢复一次最终 v14 端点头，然后固定端点、factorized 骨架、alpha、noise、gamma、控制均值、基因轴和测试查询。删除组合不会进入端点重拟合或参数重选。

## 冻结的实施定义

原 core 的训练残差使用未做增量幅度修正的 `baseline_head`：

`training effect − saturate(sum of fixed baseline endpoint effects, fixed alpha)`。

正式目标骨架使用最终 corrected v14 端点；每个子集只重算

`fixed backbone + fixed gamma × solve(K_subset + fixed noise × I, K_subset,target).T × fixed residual_subset`。

核保持原始 `incidence @ incidence.T / 2`。因此这里保留了原实现的训练残差定义，不把 corrected 端点的训练残差替换进去。

三个原 seed 的完整见证库为 28、42、36 个训练组合。三个预定 SHA256 排序生成嵌套前缀；25/50/100% 分别取 7/14/28、10/21/42、9/18/36 个组合。100% 每 seed 保存一次，由三个 order 共同引用。

所有测试查询保持固定。零直接见证定义为：所选训练组合中没有任何一个与该测试组合共享端点。它是核的直接连接数，不是额外训练的风险分数。

## 执行与身份核验

`protocol.json`、`subsets.json` 和 `protocol_freeze.json` 在 pilot 前冻结。`predictions_frozen.json` 在所有预测完成后冻结，之后才由独立 `evaluate` 子命令解码原始 `evaluation_only.npz` 并计算损失。

三个 seed 的 100% 完整输出及 factorized 输出均与原始 v14 的 float32 数组完全一致。预先允许的阈值为所有坐标 `atol=1e-6, rtol=1e-6`；实测最终预测最大误差为 0。共享加性参照也在 float32 下完全一致。

共 27 个逻辑 seed×order×fraction 设置、1,629 条逐 query 预测／损失记录、21 份实际预测 NPZ。三个 order 共用同一份 100% 输出，不应重复计算拟合次数或当作独立复制。

pilot 的拟合阶段墙钟时间为 2.552 秒；正式阶段为 3.949 秒，总计 6.501 秒。端点恢复计算合计 0.900 秒；21 次子集闭式残差求解合计 0.032 秒。较短的实际用时来自固定端点后的 7–42 行小型闭式核回归；没有重跑 ensemble。

`audit.json` 记录输出核验：完整结果与原 Step1 的逐 query 基线 MSE 最大差 2.61e−18，WMSE 最大差 8.89e−16；嵌套见证计数、冻结哈希、阶段顺序通过检查。

## 结果口径

每个 seed 和 fraction 先平均三个 order 的逐 query 损失，再平均测试 queries，除以同 seed 的固定 factorized 平均损失。三个 seed 的比值取几何平均为 rho；增量 G=1−rho。主对比为 G100−G25。

以下括号为 5,000 次共同 pair cluster bootstrap 的点态 95% 区间。前三个 G 列以百分比表示；最后一列为百分点。

| 损失／固定参照 | G25 [95% CI] | G50 [95% CI] | G100 [95% CI] | G100−G25 [95% CI] |
|---|---:|---:|---:|---:|
| 非靶 MSE／factorized，主要 | 19.10 [16.17, 22.32] | 27.42 [23.81, 31.22] | 33.25 [29.28, 37.49] | **14.15 [12.42, 16.06]** |
| 非靶 MSE／共享加性 | 21.71 [18.87, 24.74] | 29.76 [26.31, 33.37] | 35.41 [31.56, 39.49] | 13.69 [11.99, 15.53] |
| 非靶 WMSE／factorized | 26.81 [17.07, 34.52] | 33.91 [22.99, 42.74] | 37.66 [26.18, 47.74] | 10.85 [7.52, 14.36] |
| 非靶 WMSE／共享加性 | 18.10 [2.86, 29.42] | 26.05 [10.59, 38.14] | 30.24 [13.54, 43.64] | 12.14 [8.64, 15.97] |

九个 seed×order 路径的普通 MSE 均在 25→50→100% 见证库增加时改善。这个方向一致性是描述性结果，九条路径共享测试对和部分训练组合，不能作为九个独立样本。

WMSE 在九条路径中，100% 优于 25% 为 8/9，50% 优于 25% 为 7/9，100% 优于 50% 为 8/9。平均三条 order 后，三个 seed 的 MSE 和 WMSE 都在 25→50→100% 时改善。

零 incident witness 的 query×order 比例从 25% 见证库的 32.04% 降至 50% 的 20.07%，再降至 100% 的 16.57%。这些分母包含同一 query 的不同 order 预测，不能解释为独立生物样本；各 seed/order 的原始数量保存在诊断 CSV。

统计使用 `default_rng(2026101702)`，对按字典序排列的 106 个 canonical pair 并集共同生成 5,000 次 multinomial multiplicity，每次携带同一 pair 的全部 seed/order/fraction 记录。先在每个 query 内平均三条 order 的损失，再按对应 seed 的 pair multiplicity 求均值，然后对三个 seed 的损失比取几何平均。主差值区间直接来自每次同一抽样内的 `rho25−rho100`，没有相减两个端点的置信区间。全部报告均为 5,000 次有效、0 次无效；没有重抽训练集或把九条顺序当作独立 n。

统计脚本在运行前另行记录 source SHA256 与输入哈希；此前已查看点估计，这一事实在 `statistics/analysis_freeze.json` 显式记录。它是原已冻结统计口径的回顾性计算，不是新的独立验证。中央协议使用简写命名空间 `dataset`，实际见标签前冻结的子集使用 `wessels_witness_quantity`；原子集原样保留，没有根据结果重抽。

本实验是对已有系统的回顾性受控消融。固定参数来自完整原始模型；结果不能直接换算为前瞻性实验预算节省，也不证明重新训练、重新调参后的少见证系统会达到相同性能。

## 文件

- `witness_quantity.py`：冻结、pilot、正式预测、预测冻结和评价的可执行入口。
- `protocol.json`、`subsets.json`、`protocol_freeze.json`：协议、嵌套见证身份和冻结记录。
- `seed*/frozen_state.npz`：固定端点对应的骨架、训练残差、固定设计矩阵和参数。
- `seed*/restoration.json`：端点恢复参数、用时及 100% 身份验证。
- `seed*/prediction*.npz`、同名 JSON：每份实际预测、query 直接见证数、所选见证身份、用时与哈希。
- `predictions_frozen.json`：全部实际文件哈希与 27 个逻辑设置到文件的映射。
- `evaluation/query_losses.csv`：1,629 条 query×seed×order×fraction 损失；包含预测文件和行号。
- `evaluation/fixed_comparator_losses.csv`：原 final v14、factorized、共享加性三种固定参照的 543 条逐 query 损失。
- `evaluation/witness_counts.csv`：每个设置的见证库大小、零直接见证比例和求解时间。
- `evaluation/point_summary.json`：按冻结协议聚合的点估计。
- `evaluation/evaluation_manifest.json`：评价时间、上游预测冻结哈希、输出哈希。
- `pilot_report.json`、`formal_fit_report.json`、`audit.json`：执行报告与核验。
- `bootstrap_witness_quantity.py`：统计脚本，先锁定自身和输入哈希，再运行共同 pair bootstrap。
- `statistics/analysis_freeze.json`、`statistics/run_manifest.json`：统计锁定与执行元数据。
- `statistics/summary.json`、`statistics/fraction_results.csv`、`statistics/contrasts.csv`：各比例 G、点态区间及同次抽样差值区间。
- `statistics/seed_results.csv`、`statistics/seed_order_stability.csv`：各 seed／order 的方向稳定性。
- `statistics/zero_witness_by_setting.csv`、`statistics/zero_witness_summary.csv`：零 incident witness 比例和直接见证数。
- `statistics/bootstrap_draws.npz`：106 pair 顺序、全部共同 multiplicity、rho 与直接差值抽样数组，便于不重训核对区间。
- `evaluation_inputs/Wessels_seed*.evaluation_only.npz`：三个原始评价输入的逐字节副本，合计 9,819,542 字节。
- `evaluation_inputs/source_manifest.json`：评价输入来源路径、大小和 SHA256；与最初实验协议所锁的原输入 SHA 一致。
- `replay_packaged.py`、`replay_packaged_report.json`：仅用包内冻结状态重新计算预测及损失的验证入口与实测报告。

## 从冻结状态开始的可移交复现

下面的入口只依赖 Python 和 NumPy。它只读取本目录内的冻结状态、子集身份、预测、评价输入与结果表，不需要旧仓库、GO 资产、原端点训练代码或网络：

```text
python replay_packaged.py --root PACKAGED_WESSELS_FOLDER
```

默认只输出验证结果，不写文件；如需保存新的 JSON 报告，可加 `--out PACKAGED_WESSELS_FOLDER/new_replay_report.json`，已有文件不会被覆写。验证通过返回 0；验证失败返回 1。

它先核对包内冻结哈希，再从 `frozen_state.npz` 的固定骨架、训练残差、设计矩阵和参数重算 21 份预测；逐元素核对原预测、见证身份和直接见证数，然后用随包保留的原 `evaluation_only.npz` 复算所有非靶 MSE／WMSE。允许的逐元素预测误差沿用原协议 `atol=1e-6, rtol=1e-6`，报告同时记录是否 float32 完全一致。

本次实测 **PASS**：21 份预测、6,360,340 个数值全部在 float32 下完全一致；1,629 条 query 损失和 543 条固定参照损失逐行复算后，MSE 与 WMSE 最大差均为 **0**。耗时约 0.47 秒，未改动任何冻结预测或原科学结果。

这是**从冻结数值状态开始的预测和评价复现**。它没有重新训练原始端点头，也没有替代原始端点训练的复现；从原训练数据恢复端点头仍需要下节所列旧资产。

## 从原端点训练资产开始重现

需要原 Step1 的 Wessels 三个 `results/Wessels_seed*` 目录，以及含 E19 source、GO 资产及正式 manifest 的原仓库。脚本对原始输入只读；输出请选择一个新的目录。以下每条命令使用同一组参数：

```text
python witness_quantity.py freeze --repo REPO --step1-root STEP1_RUN --out NEW_OUTPUT
python witness_quantity.py pilot --repo REPO --step1-root STEP1_RUN --out NEW_OUTPUT
python witness_quantity.py predict --repo REPO --step1-root STEP1_RUN --out NEW_OUTPUT
python witness_quantity.py freeze_predictions --repo REPO --step1-root STEP1_RUN --out NEW_OUTPUT
python witness_quantity.py evaluate --repo REPO --step1-root STEP1_RUN --out NEW_OUTPUT
python bootstrap_witness_quantity.py lock --root NEW_OUTPUT --analysis-plan ANALYSIS_PLAN_JSON
python bootstrap_witness_quantity.py run --root NEW_OUTPUT --analysis-plan ANALYSIS_PLAN_JSON
```

已冻结的协议、已冻结预测及已完成评价不能被这些入口覆写。代码使用原始 core 及 Step1 已验证的 amplitude bank 优化，关闭向原 source 写入 Python 字节码。
