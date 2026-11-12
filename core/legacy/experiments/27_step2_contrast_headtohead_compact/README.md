# Experiment 27 — 第二步精简结果复现

本目录仅新增独立的小于 10 MB 的第二步复现包及其身份记录；不重复嵌入原大档案。原始逐预测数组、AutoDL 完整结果和训练输入仍在原实验档案中单独保留。

小包覆盖：12 个发布划分的端点消去分析，以及 55 个 RNA 组合的 GEARS/scouter/M0/M1/R 同台比较。它保存无损 float64 充分统计量、固定划分与基因身份、协议、复算代码和参考结果，能重新计算既定主表、敏感性及贡献分析。数值没有降精度或量化。

在总包根目录执行：

```bash
python tools/run_step2_v22.py --mode verify
python tools/run_step2_v22.py --mode reproduce --out /tmp/witnesscell_v22_step2_run1
```

也可单独解压本目录中的 ZIP，在其顶层运行：

```bash
python verify.py
python reproduce.py --output /tmp/witnesscell_step2_run1
```

输出目录必须新建。默认评分仅需现有 Python/NumPy；若显式加 `--figure`，则还需要 Matplotlib。没有安装、下载或 GPU 训练操作。

这里的精确复现是“用保留的数值统计量重新计算指定结果”，并有原始数组到统计量的独立核对记录；它不等于恢复全部逐基因预测或从头重训。需要完整预测数组进行新指标评分时，应回到原 `WitnessCell_Prediction_Contrast_Audit_20261018.zip` 与 `WitnessCell_HeadToHead_Results_20261018T193222.zip`。小包内可选逐数组评分脚本也需要显式提供这些外部原档案。

12 个发布划分和 55 组合任务保留各自评分范围。新 GEARS/scouter 结果不是旧 12 划分的匹配输出，不能填旧表缺项。神经模型种子代表优化重复；55 组合面板仍属于开发分析。R 收益集中于少数读出、方差标准化和固定十基因剔除等结果全部保留。

当前小包的字节数、SHA-256 和顶层目录见 `COMPACT_ARCHIVE.json`。本次仅重组并复算原结果，不修改模型、协议或科学结论。
