# Gate 14 - GSE146194 named-AIVC confirmation

本目录冻结并独立复核 WitnessCell 在 GSE146194 上的六策略正式比较：

`Saturation, GEARS, CPA, GEARS+Witness, CPA+Witness, WitnessCell`。

## 一句话结论

在同一批 14 个 strict-seen1、共享 FDPS 的最终靶对上，独立 WitnessCell 相对
saturation 将 interaction-residual MSE 降低 **42.49%**（pair-bootstrap 95% CI
40.31%-44.91%，14/14 靶对获胜，exact one-sided sign-flip
`p=6.10e-5`）；它也显著优于本门训练预算下的 GEARS 与 CPA。

GEARS+Witness 相对独立 WitnessCell 仍有 **1.87%** 的小幅增益（95% CI
0.45%-3.23%），说明 GEARS 保留少量互补信息。CPA+Witness 相对 WitnessCell 的
0.48% 增益区间跨零，不能主张 CPA 的稳健互补贡献。

## 目录

- `FROZEN_SCIENTIFIC_PROTOCOL.json`：具名模型预测生成前冻结的科学协议。
- `code/`：正式六策略评分代码、汇总器与冻结训练预算。
- `raw_results/`：18 个具名模型预测、18 份训练审计、18 份日志、9 个六策略评分和 3 个 CPU truth anchors。
- `independent_audit/`：不导入原评分器的独立重建、逐靶复算、统计表和审计报告。
- `figures/`：论文主图、诊断图、作图脚本和图注。

## 独立复核

```bash
python independent_audit/independent_audit.py \
  --experiment-root . \
  --out independent_audit
```

审计器会从原始 GEARS/CPA prediction archives 和 CPU anchors 重建全部六策略，
逐项比对 756 行目标指标，并以 14 个生物学靶对为统计单位重新计算效应量、
bootstrap 区间和精确 sign-flip 检验。

## 重画论文图

```bash
MPLCONFIGDIR=/tmp/witnesscell-mpl \
python figures/make_paper_figures.py
```

## 结论边界

- 14 个目标都共享 FDPS；这是一个 pathway/hub panel，不代表任意基因对总体。
- 本门测试 strict-seen1，不测试 seen0。
- 3 个 data seeds 和 3 个 model seeds 是稳定性重复，不是独立生物学重复。
- 本协议在生成具名模型预测前冻结，但 14 个靶对的生物学结果已用于先前的
  CPU WitnessCell 实验；因此它是 model-blind named-AIVC confirmation，
  不是从未查看过的全新生物学队列。
- 正式自动状态名中的 `UNIVERSAL_AUGMENTATION` 只表示四个预冻结比较均通过；
  论文不使用“万能 wrapper”措辞。融合权重显示 GEARS+Witness 和 CPA+Witness
  分别平均使用 91.83% 和 97.90% Witness，因此应解释为独立 WitnessCell 已成立、
  GEARS 有小幅互补、CPA 在本门主要被替代。

