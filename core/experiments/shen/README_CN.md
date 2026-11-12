# WitnessCell Shen2017：正式实验已完成

下载并恢复了作者历史仓库中的原始计数；三细胞系均完成真实拟合、预测冻结和正式评分。

主秩一模型相对M0：HeLa完整误差降低4.16%、伙伴误差增加1.10%；A549降低3.65%、伙伴增加0.11%；293T降低6.18%、伙伴增加0.53%。预设伙伴主终点未改善，三细胞等权平均增益−0.5784%。

详细结论、下载过程及来源限定见 `WitnessCell_Shen2017_Formal_20261108_Report.md`。

## 核查与复现

先运行 `python code/verify_bundle.py` 核查归档文件。

安装 `requirements.txt` 后：

```bash
python code/independent_audit.py --out /tmp/shen_audit_fresh
python code/replay_formal.py --out /tmp/shen_replay_fresh
```

两个输出目录均须不存在。第二个命令从原始计数的规范表重新准备测量、拟合三细胞模型、统一冻结预测后评分；不会覆写本包的正式预测。已执行一次这样的数值重放，9个预测数组最大差均为0。

`protocol/PROTOCOL.json`及`protocol/ORIGINAL_PROTOCOL_CN.md`保留2026-11-07的原始本地锁定记录，文中“源文件尚未下载”描述当时状态。当前完成状态以`RUN_STATUS.json`、`audit/SOURCE_AND_FIT_LOCK.json`及正式报告为准。来源核对落实了计数映射，模型、阈值、划分和主终点不变。

`runs/HeLa_source_smoke`是真实HeLa源内烟测；`runs/formal`是真实三细胞正式拟合。数值重放不计为新的外部研究。没有外部预注册，没有重新训练GEARS或Scouter。

主要结果和图形由`results/formal/SCORES.json`和`audit/INDEPENDENT_SCORES.csv`提供；完整原始计数、模型权重、内层成绩、真值与预测均已包括。所有原文件与第三方代码保持原许可。
