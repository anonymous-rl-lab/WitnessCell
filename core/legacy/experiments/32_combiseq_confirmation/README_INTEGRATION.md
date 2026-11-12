# Combi-seq 冻结程序确认：保留阴性结果

结果为 **NOT_CONFIRMED**。374 个药物组合、36 个端点、3 次实验运行，在固定 51 基因程序上，相对同一个完整端点参照：完整误差改善 **−0.0232823%**，去端点伙伴误差改善 **−0.0302401%**。正数才表示降低误差。

统一复算入口（从本目录运行，输出目录须为空且位于模块外）：

```bash
python run.py --out /tmp/combiseq_confirmation_replay_new
```

以下是保留的分项入口：

```bash
python REPRODUCE.py --out /tmp/combiseq_replay_new
```

运行从包内 2.25 MB 压缩计数矩阵与固定元数据开始，重新生成条件均值、训练选择、拟合与预测；检查 108 个命名数组的 dtype/shape/字节哈希、全部评分表和科学结论。原环境约数秒完成，不需要网络或旧工作区。Python 3.12、NumPy 2.3.5、Pandas 2.2.3；绘图另需 Matplotlib。

这是一项新公开来源、跨干预类型的回顾性确认：它研究药物组合，不是新增双基因 RNA 实验。全部运行及组合方向共同留出；来源、输入、预定义程序与冻结记录见原 `README.txt`、`CONFIRMATION_PROTOCOL.txt`、`FREEZE.json`。

阴性结论保持不变。相对更弱的饱和模型的完整增益，不能冒充新增伙伴分量的收益。复现本身也不构成第二次独立确认。

`validation/INTEGRATION_REPLAY.json` 保存本次 108 数组精确核验记录。`INTEGRATION_PROVENANCE.json` 记录原 ZIP SHA-256。原始文件逐字节保留；原 `PACKAGE_MANIFEST.json` 是原归档成员清单，新代码总包清单覆盖本次新增文档。
