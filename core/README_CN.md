# WitnessCell 复现代码总包 v28（论文 v36）

本发行版以论文 v36 为固定基准。论文版本与代码发行版本分别编号：**paper v36 / code 28.0.0**。原始论文正文、SI、数据、模型权重和科学结果未改写。轻便论文另见 `WitnessCell_v36_Light.zip`。

## 快速开始

使用 Python 3.12 和 `requirements-tested.txt` 所列依赖。默认流程不重新训练生物学预测模型；它复算保存结果、进行冻结权重推断、重新估计 Han 的目标控制系数，并重建图件。RNA 接口测试包含小型人工输入拟合。

```bash
python -B tools/verify.py
python -B tools/reproduce.py --out ../v36_reproduction_new
```

输出必须是包外尚不存在的新目录，避免覆盖参考数据。总入口有逐任务日志、失败状态、超时和最终 JSON 报告。可按任务复现：

```bash
python -B tools/reproduce.py --task han --out ../han_replay_new
python -B tools/reproduce.py --task paper_scores,theory_growth,component_control,figures --out ../paper_replay_new
```

完整任务与实际入口见 `REPRODUCTION_INDEX.json`，各图、表对应见 `docs/PAPER_COVERAGE.csv`。默认九组任务均有本轮实际执行记录，见 `qa/VALIDATION_SUMMARY.json`。环境准备时间不包含在运行记录中。

## 目录

| 路径 | 内容 |
| --- | --- |
| `paper/` | v36 Markdown 主文、SI、图注、当前图件源数据、图件脚本和参考矢量图；无整套论文 PDF 或 PNG 副本 |
| `legacy/` | 原 grn_repo_v27 全部文件原样展开，涵盖核心预测、理论、RNA、功能、风险和选择等早期成果；其中版本说明属于历史发行记录 |
| `experiments/shen/` | Shen 正式实验、保存预测和原协议 |
| `experiments/source_fits/` | SF12/V4 原训练实现、候选配置与冻结权重；默认不重训 |
| `experiments/historical_calibration/` | 训练侧历史诊断与校准证据 |
| `experiments/han2017/` | 与论文直接相关的 Han 全量实验及原主分析摘要、数据、模型、协议、重建器与压缩清单 |
| `paper/partner_control/` | 分量控制接口、四组对照、针对性检查和参考输出 |
| `paper/software/` | RNA 核心与功能接口及示例 |
| `tools/` | 当前统一核验、复算、图件重建入口 |
| `qa/` | 本轮实际运行报告和数值差异记录 |

## “精准复现”的具体范围

保存文件按 SHA256 校验；Han 被剥离的 180 个重复 NPZ 在本轮全部逐字节恢复。原 SF12 主要结果完全重现。V4 在本轮 PyTorch 2.10.0+cpu 环境存在末位浮点差异，按 v36 已有的逐阶段误差传播规则检查通过，**不称为 V4 全数值逐位相同**。

默认复算覆盖论文报告的保存指标、保存预测及当前图表；Han 还从压缩计数/发表工作表重建输入并运行冻结模型。原始 FASTQ 到计数、全部竞争模型的原始单细胞预处理和重新训练不属于默认离线入口。原始数据索引、获取说明、权重索引及历史训练实现保留。详见 `docs/REPRODUCTION_LEVELS.md`。

旧稿中的 `reproduction/*.zip` 路径现按 `docs/PATH_MIGRATION.csv` 映射，不必重新下载这些被整合的归档。科学源码的许可边界保持原状，见 `RIGHTS_AND_RELEASE.md`；此打包不新增或扩大授权。
