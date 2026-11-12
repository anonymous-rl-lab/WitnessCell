# WitnessCell 纯计算预注册 v2 最终审计报告

## 总体裁决

`PASS_READY_FOR_NEW_OSF_REGISTRATION`

本包已达到重新注册条件。该裁决只说明注册文本、冻结合同、一键执行链和审计器已经通过注册前质量控制，不是注册后的正式确认性结果。

## 审计范围

1. 新旧注册边界与既往运行披露；
2. 纯计算研究口径及禁止性表述；
3. v14 源码包 SHA-256 与内部 `MANIFEST.sha256`；
4. Python 3.12 精确环境从零重建；
5. 一键入口从解压到结果归档的端到端运行；
6. 四数据集、三种子、12 个 split 的重新生成；
7. 156 个部署数组的 schema、dtype、shape 与逐元素一致性；
8. 八个激活 split 与四个精确回退 split；
9. 654 个 condition--seed 单元的 v14--v13 配对结果；
10. Gate18 继承榜单聚合的完整性；
11. 错误源码哈希与预测数组篡改的 fail-closed 负向测试；
12. 压缩包清洁度、体积与审计材料完整性。

## 端到端结果

- 精确环境：PASS；
- 源码内部 manifest：PASS；
- Gate19、Gate18 与 repository 回归入口：全部 PASS；
- 重新生成 split：12/12；
- 比较数组：156/156；
- 全局最大绝对差：0；
- 门控：8 个激活、4 个精确回退；
- 配对单元：654；
- C1--C5：全部 PASS；
- 官方继承聚合：平均 operation rank 3.871559633，排名第 1。

## 审计中发现并修复的问题

### 1. 删除无关强制依赖

流程不读取 h5ad，因此删除了 `h5py` 强制依赖，避免运行前无意义下载或环境失败。

### 2. 写死 Gate18/Gate19 边界

Gate19/v14 重新生成的是均值预测与增量端点头。六指标表来自 Gate18 的冻结继承证据。本包已将 C5 改为“继承 Gate18 聚合完整性审计”，禁止写成新 v14 raw-h5ad 六指标重算。

### 3. 缩减结果包

原脚本会把解压后的源码副本和内部 `predictions.npz` 一并归档，沙盒结果包约 402 MB。修复后排除这两类文件，测试归档约 38 MB，同时保留 12 个 target-free deploy 档案、12 个 manifest、逐 split 日志、环境、配对结果和最终 verdict。

### 4. 加强环境与日志冻结

直接依赖和会影响执行的传递依赖均固定版本；每个 dataset--seed 执行单独保存日志；记录 UTC 开始与结束时间。

## Fail-closed 测试

- 输入错误源码包：在 SHA-256 门被拒绝，退出码 4；
- 人为修改一个预测值 `1e-6`：审计器定位 Norman seed 1 的 `prediction` 数组并拒绝 PASS；
- 恢复原数组：C1--C5 再次全部 PASS。

## 注册后允许写入论文的表述

仅当新 OSF 注册时间戳之后的一次完整运行通过 C1--C5，才可写：

> A preregistered confirmatory computational replication reproduced the frozen WitnessCell v14 deployment archives, evidence-gate branches and paired mean-response results, together with the integrity of the explicitly inherited Gate18 benchmark aggregation.

不得写为前瞻性生物验证、新扰动实验、135 对真值验证或 v14 raw-h5ad 六指标重新评分。

## 最终操作要求

1. 只上传本次最终 zip，不使用此前任何 v1 或较早 v2 文件；
2. OSF 接受新注册并生成时间戳之后再运行；
3. 不修改源码、注册包、阈值或一键脚本；
4. 数值失败不得选择性重跑；基础设施在产生数值判据前中断，才允许从零完整重启一次；
5. 上传完整正式结果包及其 `.sha256` 文件。
