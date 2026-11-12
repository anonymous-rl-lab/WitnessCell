# WitnessCell 精准复现代码总包 v27

本版保留v26全部科学内容，并加入已结束的GSE220974五模型实验：`experiments/34_gse220974_five_model_final/`。GEARS、Scouter各五折×三种子共30项训练全部完成、核验并逐项存档。此次仅封存与复现入口整合，不启动新训练、不改模型或结果。论文版本仍为v22，代码发行版本为27.0.0。

## 入口

```bash
python tools/verify.py
python tools/reproduce.py --task gse220974_final --out /tmp/witnesscell_v27_final
```

默认离线、无GPU；新增实验只需NumPy。输出须为包外的新目录。参考环境见requirements-reproduction.txt。完整运行原有任务及新增任务可用 `python tools/reproduce.py --out /tmp/witnesscell_v27_all`，仍遵循各任务原复现层次与限时。

本次新增任务从无损保存的逐基因预测误差和候选评分重新计算五模型主表、420条选择及6组区间。所有26个存储数组按原dtype和SHA256逐位恢复；主表及420条选择在实际验证中与原CSV字节一致。完整15模型原始预测复算也已执行，平台/BLAS末位误差以QA实测记录为准，不用四舍五入隐藏。

全部原始预测、真值、目标向量、30个模型权重和细胞数据仍通过已存档的SHA索引恢复。新增实验README提供原始预测复算与原模型权重重放入口；不能把小包离线评分说成全量重训，或承诺跨硬件随机训练轨迹逐位相同。

## 文件与尺寸

本次实验独立增量ZIP严格小于10,000,000字节；总代码原v26约102MB内容保留，因此总包不受10MB限制。发行尺寸、原文件保留情况与测试范围见RELEASE_STATUS.json及qa/V27_FINAL_INTEGRATION.json。

v26原有1053个文件全部可追溯：除6个发行导航文件外，其余保留原相对路径和原字节；6个原文在provenance/base_release_v26。保护实验25–32的文件仍全部保持。旧任务的历史验证不冒充本次重跑；本次只运行新增实验及一个原有入口烟测。v26完整说明保留在provenance/base_release_v26/README.md。

## 结果边界

主表仅R-G、GEARS、Scouter、标准线性、单扰动相加。R-G错选损失及命中数点估计最好，标准线性RNA误差最低，GEARS伙伴误差最低。49种检出标签、8个物理锚点、28任务、5个大池及4个伙伴自由度的原限制完整保留；不宣称全面显著优势，不删除负结果。
