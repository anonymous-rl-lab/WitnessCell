# WitnessCell v14 Gate 19 正式报告

## 正式结论

v14 已作为独立完整版本并入 v13 全部 Gate 00–18 资产。算法不是替换 v13，而是在
冻结 v13 背景、因子化骨架、交互见证核、split 和 validation 选择后，只对 unseen
endpoint 追加训练侧幅度残差头。

12 个正式 split 中 8 个通过双升级门；Wessels 三个 seed 与 Schmidt seed 3 自动
逐元素回退 v13。654 个 condition×seed 单元相对 v13：

| 指标 | v13 | v14 | v14 相对变化 |
|---|---:|---:|---:|
| 全基因 MSE | 0.00309401 | 0.00309067 | -0.108% |
| top-100 MSE | 0.0471909 | 0.0470472 | -0.305% |
| top-100 PCC | 0.787607 | 0.788908 | +0.00130 |

相对冻结 v1，v14 达到全基因 MSE `-5.157%`、top-100 MSE `-9.435%`、PCC
`+0.08357`。top-100/top-5000 方向二指标榜均保持第 1，published-best PCC/MSE
单元保持 `18/24` 与 `14/24`。

## 底层机制

v13 已经给出可靠的方向背景；剩余误差主要是 unseen endpoint 自响应幅度不恒定。
v14 将目标基因对已知扰动的响应指纹压成两个稳定标量，与通用 self response 一起
拟合 v13 的 top-100 残差。高维 reciprocity、GO 稀疏程序与 response-kNN 在严格
outer-LOO 中均表现为方差过大或指标交换，因此没有进入正式头。

最终相对 ridge 固定为 `0.2`，系数限制在 `[-1.5, 1.5]`。每个 outer endpoint
都会先从 v13 refit 和 correction fit 中移除；只有 all-gene MSE 改善 LCB 与
top-100 PCC 改善 LCB 同时大于 0 才激活。top-100 MSE 只记录，不充当第三道门。

## 正式复现与发布边界

- 12 个 split 已从 v13 条件矩缓存重新运行；除可移植化缓存路径外，manifest 与
  候选运行逐字段一致。
- 654 单元配对结果及方向榜与候选结果逐值一致。
- 正式部署目录只含 target-free archives，不含 `truth` 或 `truth_variance`。
- 4 个拒绝升级的 split 与 v13 prediction 逐元素完全相同。
- Gate 18 的官方六指标 SOTA 证据完整保留，仍证明 v13 冻结协议总榜第 1。
- 本次没有原始 h5ad 分布矩阵，故不重算 v14 Wasserstein、E-distance、KL 与
  Common-DEGs；这些 v13 表不会被改名成 v14 新结果。

因此 v14 是“正式均值响应算法与部署版本”，不是“已经新跑完分布六指标的版本”。
这一边界由 `audit.py` 和仓库回归测试强制检查。
