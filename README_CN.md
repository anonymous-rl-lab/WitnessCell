# WitnessCell 匿名 GitHub 复现仓库

这是由 `grn_repo_v18` 整理得到的双盲公开版本，包含：

- 工业化 `witnesscell 0.1.0` Python API、CLI、测试与发布工作流；
- WitnessCell v14 最终条件均值算法及四数据集三种子复现入口；
- scPerturBench 正式评分表、Gate 21 选择性预测与 Gate 22 指标压力测试；
- 精确 Witness Risk、estimated Witness Risk 与几何规律验证代码。

原始 v18 证据包约 430 MiB。本仓库删除了历史嵌套压缩包、重复版本预测、原始
公开单细胞矩阵、缓存、编译产物和包含机器绝对路径的日志；保留了最终算法、冻结
协议、紧凑结果、必要派生输入和原资产校验值。

## 1. 快速验证

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-repro-v14.txt
python -m pip install -e ".[dev]" --no-deps
bash scripts/run_smoke.sh
```

## 2. 单个正式 split 复现

```bash
DATASETS=Norman SEEDS=1 bash scripts/reproduce_v14.sh
```

## 3. 四数据集 × 三种子完整均值预测

```bash
bash scripts/reproduce_v14.sh
```

默认只保留不含测试真值的 `deploy_predictions.npz`，并按数组内容而非 ZIP 时间戳
核验冻结语义摘要。

## 4. 完整细胞级评分

```bash
bash scripts/download_full_data.sh
```

该步骤会下载约 0.9 GiB 的四个公开 `.h5ad.gz`。详细环境隔离、预计耗时、数据
校验值和完整命令见 `REPRODUCIBILITY.md` 与 `DATA.md`。

公开前请按 `docs/GITHUB_UPLOAD.md` 操作。仓库保持双盲匿名；不要补入作者、单位、
个人 GitHub 账号、项目 DOI、邮箱或 OSF 身份信息。
