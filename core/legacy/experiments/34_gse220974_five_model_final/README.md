# GSE220974 五模型最终封存包

状态：实验截至 2026-11-07 结束。GEARS、Scouter 各五折×三种子，共30项全部冻结、原模型权重重放通过并独立存档。此包只整理和复现现有结果，不启动训练或调参。主表只含 R-G、GEARS、Scouter、标准线性、单扰动相加；两神经模型使用三个种子预测的均值后评分。

## 默认离线复算（NumPy即可，无GPU、无联网）

```bash
python reproduce.py --out /tmp/witnesscell_gse220974_replay
```

输出目录必须是包外的新目录。重新计算五模型主表全部字段、15个已归档模型/种子的420条逐任务选择、6组物理锚点聚类配对区间。RNA误差与伙伴误差从完整49×4987逐基因误差矩阵重算，伙伴指标同时用SVD与投影矩阵独立求得；不会读取已保存的MSE数值作为计算输入。参考表只在计算后用于核验。26个压缩数组逐一按原dtype、shape和SHA256恢复，无降精度、舍入或有损压缩。真实程序得分和预测程序得分是已冻结的输入读出；默认入口不重新计算原始表达矩阵的目标投影。

## 完整原始预测复算

先下载 EXTERNAL_ASSETS.json 中的原始完整比较ZIP（约26.3MB），然后：

```bash
python reproduce.py --out /tmp/witnesscell_gse220974_full --reference-zip /path/WitnessCell_Completed_Comparison_20261106.zip
```

此入口验证整个外部ZIP及关键文件SHA，重新读取所有原始预测/真值/目标向量，核对小包全部26个数组（直接存储与减法形成的数组逐位核对；目标投影矩阵乘法允许已声明的BLAS末位舍入容差），再执行未经改动的原始score.py，复算全部15行（含各单种子与历史参照）的RNA及伙伴指标。这些额外模型保留为来源记录，不增加到五模型主表。

## 已训练权重重放与重训边界

reference/FINAL_WEIGHTS_ARCHIVE_INDEX.json列出30个独立权重ZIP、准确文件名、SHA和下载身份；完整数据/基础代码/GPU适配包见EXTERNAL_ASSETS.json。先恢复基础代码与细胞数据，再叠加30个完成ZIP。原恢复入口为 source/run/restore_completed_batches.py，原核验入口为 source/run/verify_neural_predictions.py；都应在恢复的原实验目录执行。FULL_INPUT_HASHES.json锁定实际输入文件。原训练代码、损失、批次、划分、种子、预算、内层选轮、最终轮数、曲线及CPU/GPU来源均保留。

30个权重归档合计超过1GB，无法与原始细胞矩阵一起放入10MB；小包不宣称包含这些外部资产，也不把冻结结果复算称为从细胞数据全量重训。外部资产均已有独立存档，完整重放必须按SHA核验后使用。CPU/GPU和Torch版本的训练轨迹不保证逐位相同；最终已冻结权重才是本次成绩的精确模型身份。本次封存不重新训练模型。

## 结果与适用范围

reference/MAIN_FIVE_MODELS.csv是最终主表，reference/PAIRED_COMPARISONS.csv给出描述性区间。R-G在主表的错选损失和命中数点估计领先；RNA误差最低为标准线性，伙伴误差最低为GEARS。不得宣称所有指标均领先或全面显著获胜。

保持原协议限制：49种检出双干预标签尚非已核实设计构建体白名单；28个选择任务来自8个物理锚点，只有5个大候选池，伙伴自由度合计4。最优指固定候选池中RNA程序得分最高。既往其他模型外层结果已经查看，本轮GEARS/Scouter未依据外层结果改预算或选轮规则。没有增加实验、删除不利结果或改变研究主张。

## 总包集成

本实验并入总代码 experiments/34_gse220974_five_model_final/，统一入口：

```bash
python tools/verify.py
python tools/reproduce.py --task gse220974_final --out /tmp/witnesscell_final_replay
```

10MB限制适用于本实验独立增量ZIP，按严格小于10,000,000字节验收。总代码保留原v26全部科学内容，体积不受本增量限制。
