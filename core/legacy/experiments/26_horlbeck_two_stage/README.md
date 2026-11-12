# Experiment 26 — Horlbeck two-stage reproducibility

本目录把已完成的 K562 开发比较与 Jurkat 保留环境确认并入 `grn_repo_v21`。原精简包保持原字节，统一入口位于仓库根目录的 `tools/run_horlbeck_v21.py`。

## 内容与身份

| 文件 | 内容 |
|---|---|
| `WitnessCell_Horlbeck_TwoStage_Repro_20261018.zip` | 181 个文件，包含两阶段代码、派生输入、划分、保存模型/预测、结果与复现入口 |
| `registration/ANALYSIS_AMENDMENT.md` | 方法调整、数据访问和冻结时间线 |
| `registration/REGISTERED_PLAN_STATUS.md` | 原注册假设结果与后续两阶段分析的关系 |

精简 ZIP 为 **2,491,533 B**，展开文件合计 **5,317,664 B**。SHA-256：

```text
b60e01230f3f43f156271697fe1fdc52b6c4ce5f68b522c3664c6e4aa6ff4082
```

K562 部分包括三个外层开发折和两个训练预算；Jurkat 部分是一个外层测试集、三个训练内部选择折。原先两个较大实验包没有作为嵌套 ZIP 重复装入，但复现所需阶段内容已纳入本精简包。其来源对照见内包 `SOURCE_ARCHIVES.json`。

## 从总代码仓库根目录运行

使用已有兼容 Python、NumPy、SciPy；无需 CUDA、PyTorch 或重建环境。

```bash
python tools/run_horlbeck_v21.py --mode verify
python tools/run_horlbeck_v21.py --mode reproduce --stage all \
  --out /tmp/witnesscell_v21_horlbeck_run1
```

核验模式复算保存结果，不重新训练。重新拟合模式在新仓库外目录展开工作副本，按原流程重新拟合、训练内选择和评价；主要检查报告为：

```text
/tmp/witnesscell_v21_horlbeck_run1/reproduced/REPRODUCTION_CHECK.json
```

`--stage K562` 或 `--stage Jurkat` 可单独运行一个阶段。任何输出目录必须是新目录；不覆盖当前保存资产。该实验不依赖历史大数据恢复；只有调用早期实验时才需执行总仓库的 `tools/restore_data_v21.py --verify`。

## 结果与研究身份

| 比较 | 共享秩一增益 | 对称谱收缩增益 |
|---|---:|---:|
| K562 低预算，相对同预算 M0 | 0.2652% | −0.2879% |
| K562 高预算，相对同预算 M0 | 11.0292% | 11.3865% |
| K562 高预算，消去自由端点后 | 11.6838% | 10.6710% |
| Jurkat 测试误差 | **6.4308%，主要方法** | 8.1606%，次要方法 |
| Jurkat 消去自由端点后的误差 | 6.1750% | 5.8975% |

M0 按各训练子集、各重复重新拟合；共享修正矩阵 G 按相同信息和原单位残差拟合。评价时所有模型使用同一测试集合。两次实验重复分别报告，不能当作相互独立的组合数。端点消去增益的分母为同空间的 M0 误差。

原 OSF 注册 K562 主检验未达到正面证据标准，其原始结论保留。后续 K562 三折属于开发；Jurkat 在修订流程本地冻结后评分，属于同一原研究中另一细胞环境的确认。新方法未被追溯写成原注册方法；该结果也不被表述为新增湿实验、RNA 预测确认或独立研究复现。

- Project 匿名访问：<https://osf.io/psg7v/overview?view_only=db54056fdde34bc7b0e5ec5f12135bb4>
- 原注册匿名访问：<https://osf.io/2uwqk/overview?view_only=58325479562749a1a25bf1002fd13cfe>

本目录展示的是已经完成实验的保存结果；本次整合后的实际检查与重新拟合状态以根目录 `RELEASE_STATUS.json` 及 v21 验收报告为准。
