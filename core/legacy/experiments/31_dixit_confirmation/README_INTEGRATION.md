# Dixit：冻结 RNA 结局与后续 G2M 开发

本模块同时保留两个不同实验，不能将后者的正结果替代前者的阴性主终点。

| 证据 | 结果 | 本包复现层级 |
| --- | --- | --- |
| 冻结全 RNA，45 对、10 TF、23,519 非靶读出 | 原 M0 相对单扰动加性完整误差改善 23.4328%；新增项完整改善 0.2698%，伙伴改善 −0.0326% | 从逐组合完整 MSE 与逐折 SSE 重新汇总，核对原始协议和结果；不包含原始全 RNA 预测数组及其重新拟合 |
| 后续固定参考 G2M 开发 | 完整改善 2.2671%，伙伴改善 14.2994% | 从双精度条件均值重新拟合、选模、预测、独立重构、导出查询模型；237 个命名数组逐位核验 |

统一复算入口（从本目录运行，输出目录须为空且位于模块外）：

```bash
python run.py --out /tmp/dixit_confirmation_replay_new
```

以下是保留的分项入口：

```bash
python replay_primary_summary.py --out /tmp/dixit_primary_audit.json
python REPRODUCE.py --output /tmp/dixit_g2m_replay_new
python model/predict_program.py --state model/deploy/Dixit_G2M_fold0.npz --pairs CREB1+YY1
```

`REPRODUCE.py` 的输出目录必须为空；正常运行约数秒，不需要网络。Python 3、NumPy；绘图另需 Matplotlib。原始程序、冻结记录、来源哈希、失败候选全部保留。新增文件不改变 `SOURCE_MANIFEST.json` 中冻结源码的身份。

`primary_rna_archive/` 是 v22 论文已使用的冻结主终点记录。`replay_primary_summary.py` 验证其来源哈希，重新汇总 8 方法×3 外折、45 对的完整误差，检查端点分解及增益分母。这里的复算层级明确低于完整 RNA 重新训练。

根目录原始 `README_CN.txt` 与 `REPRODUCE.py` 对应后续 G2M 开发。其精简输入包括全 RNA 条件均值，保留训练侧标准化所需信息；不含逐细胞矩阵、半样本或方差。重跑得到的是已揭示数据上的复现，既不是新确认，也不是前瞻验证。

`validation/INTEGRATION_REPLAY.json` 记录本次集成验证。`INTEGRATION_PROVENANCE.json` 记录原 ZIP 的 SHA-256。原 `MANIFEST_SHA256.json` 对应原始归档成员；新总包的清单覆盖全部新增文件。
