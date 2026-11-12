# WitnessCell v24 补充材料

本包配套 `WitnessCell_v24_Main_Methods.md`，按照正文六节 Results 重构，正文原文件未修改。主交付为英文 Markdown SI，包含12个补充说明和24张补充表。

## 文件

- `WitnessCell_v24_Supplementary_Information.md`：完整SI正文，含理论证明、算法细节、完整对照、敏感性分析和复现范围。
- `WitnessCell_v24_Main_Methods.md`：原v24主文与Methods，原样收录。
- `CLAIM_SOURCE_MAP.csv`：正文结果—SI章节／表格—证据文件逐项对应。
- `NEURAL_TRAINING_COMPLETION.csv`：GSE220974全部30项神经模型训练的实际选中轮数与完成记录。
- `SOURCE_MANIFEST.json`：164个来源文件的路径、字节数和SHA-256。
- `SI_VERIFICATION.json`：本次文档检查及已保存结果复算记录。
- `additional_development_tables.md`：旧版的系统集成、测量校准和共享特征开发表，保留其原编号，避免将其误认作新增独立确认。
- `source_data/`：上述来源文件的副本；未更改原始数值或不利结果。
- `authoring/`：SI文本模板和组装／检查脚本。

## 与正文的互补关系

| 正文 | SI支撑 |
|---|---|
| 方法与流程 | Note 1方法身份；Note 3实施细节；Note 4端点空间证明；Note 5伙伴核构造 |
| 联合干预辨识信息 | Note 2：不可辨识风险下界、置信集、强GLS对照、预算依据与机制家族边界 |
| 正式组合预测性能 | Note 3：完整16方法榜单、六指标排名、配对区间、掩码与支持量对照 |
| 伙伴响应 | Notes 4–6：全部结构分解、55对面板、11矩形、种子／折结果、权重控制与读出集中度 |
| 功能表型 | Note 7：K562训练预算与全部折、Jurkat主／次终点及条件敏感性 |
| 风险与PairPick应用 | Notes 8–9：风险推导与校准、R-G精确约束、全部15行结果／28候选池／30项训练 |
| 讨论中的迁移边界 | Note 10：Dixit、Combi-seq完整结果与单步采集外层结果；Note 11原计划与修订 |
| 数据与代码可用性 | Note 12及来源清单：区分评分复算、权重重放、紧凑重拟合与完整重训 |

## 本次实际核验

已核对164个来源文件的哈希、24表与12节编号；从保存的1000个种子结果重新计算因果配对风险差和区间；从RNA逐折SSE重算M1／R增益；从420条决策记录重新计算全部15模型的效用、错选损失与命中数；核对全部30项神经训练的选择轮数和最终训练轮数一致。主文哈希与原文件一致。

这些检查没有重新训练模型。GSE220974的数值差为浮点汇总量级，完整记录见 `SI_VERIFICATION.json`。

## 复建本SI

在已安装NumPy和SciPy的Python环境，从本包根目录执行：

```bash
python authoring/build_si.py
```

该命令仅从包内来源表重建SI并检查记录，不联网，不训练模型。完整实验仍使用 `grn_repo_v27` 的各任务入口。代码包内 `paper` 任务对应旧v22图表，不代表已生成v24新图。本包没有新增插图、公开代码地址或未经核实的许可证。

PairPick (R-G) 与归档 `R_G` 一一对应；PairPick-Stable、PairPick-TargetRank保留各自身份。v24的新文稿组织不改变原实验的开发、局部冻结或外部比较状态。
