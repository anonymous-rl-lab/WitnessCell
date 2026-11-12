"""Finish independent verification and package the complete five-model comparison.

Run only after finish_and_backup.py has exited successfully. This script does
not fit models, select checkpoints or change the experimental protocol.
"""
from pathlib import Path
import csv, hashlib, json, shutil, subprocess, sys, zipfile
import numpy as np

R = Path(__file__).resolve().parent
E = R / 'results/detected49/evaluation_main_20261105'
OUT = Path('/workspace/scratch/220c4aa7ee53/WitnessCell_Completed_Comparison_20261106')

def readcsv(p):
    with p.open() as f:
        return {r['model']: r for r in csv.DictReader(f)}

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        while b := f.read(8 * 1024 * 1024):
            h.update(b)
    return h.hexdigest()

def main():
    required = [R / 'results/detected49' / f'{k}_seed{s}' / f'fold_{f}' / 'FREEZE.json'
                for k in ('gears', 'scouter') for s in range(3) for f in range(5)]
    assert all(p.exists() for p in required), 'All thirty jobs must be frozen first.'
    assert json.loads((E / 'EVALUATION.json').read_text())['status'] == 'COMPLETE_MAIN_ROSTER_15_ROWS'
    subprocess.run([sys.executable, str(R / 'verify_neural_predictions.py')], cwd=R, check=True)
    OUT.mkdir(exist_ok=True)
    shutil.copyfile(R / 'portable_score_completed.py', OUT / 'score.py')
    subprocess.run([sys.executable, str(OUT / 'score.py'), '--source', str(R)], check=True)
    official, portable = readcsv(E / 'SUMMARY.csv'), readcsv(OUT / 'SUMMARY.csv')
    assert official.keys() == portable.keys()
    checks = []
    mapping = [('mean_U', 'mean_U'), ('mean_regret', 'mean_regret'),
               ('RNA_MSE', 'RNA_MSE'), ('partner_projection_MSE', 'partner_MSE'),
               ('half0_mean_U', 'half0_U'), ('half1_mean_U', 'half1_U')]
    for name, row in official.items():
        for a, b in mapping:
            if not row[a]:
                assert not portable[name][b]
                continue
            x, y = float(row[a]), float(portable[name][b])
            np.testing.assert_allclose(x, y, rtol=1e-10, atol=1e-12)
            checks.append(dict(model=name,metric=a,absolute_difference=abs(x-y)))
        assert round(float(row['top_choice_accuracy']) * int(row['choice_tasks'])) == int(portable[name]['hits'])
    official_choices = {}
    with (E / 'CHOICES.csv').open() as f:
        for row in csv.DictReader(f):
            official_choices[(row['model'], row['group'])] = row['pair']
    with (OUT / 'CHOICES.csv').open() as f:
        for row in csv.DictReader(f):
            assert official_choices[(row['model'], row['group'])] == row['pair']
    (OUT / 'INDEPENDENT_SCORE_CHECK.json').write_text(json.dumps(dict(
        status='PASS', all_model_metrics=True, all_choices=True, checks=checks), indent=2)+'\n')
    for name in ('NATIVE_WEIGHT_REPLAY.json', 'PAIRED_COMPARISONS.csv', 'EVALUATION.json', 'REPORT_CN.txt', 'GROUPS.json'):
        shutil.copyfile(E / name, OUT / name)
    for name in ('PROTOCOL.md', 'RUNTIME_NOTES_20261105.txt', 'requirements_runtime_20261105.txt',
                 'RUNTIME_EQUIVALENCE.json', 'RUNTIME_EQUIVALENCE_FP64.json', 'DATA_INTERFACE_EQUIVALENCE.json',
                 'AUTODL_DEPLOYMENT_HANDOFF_20261107.json', 'EXTERNAL_JOBS_20261107.json'):
        if (R / name).exists():
            shutil.copyfile(R / name, OUT / name)
    code = OUT / 'code'; code.mkdir(exist_ok=True)
    for name in ('common.py', 'train_cell_models.py', 'train_cell_models_runtime.py',
                 'lazy_cell_data.py', 'fast_gears_ops.py', 'score_main_20261105.py',
                 'verify_neural_predictions.py', 'finalize_comparison_20261105.py'):
        shutil.copyfile(R / name, code / name)
    readme = '''WitnessCell：完整五模型同台比较

状态：GEARS、Scouter 各五折、每折三个种子，共三十项主训练完成。
MAIN_FIVE_MODELS.csv 为主表；SUMMARY.csv 包括逐种子及历史参照。
主表使用三个种子预测的均值再评分，不是选择最好的种子。

轻量复算：安装 NumPy 后在本目录执行 python score.py。
scoring_inputs.npz 包含冻结预测、评价真值、划分和目标读出。
其结果已与原评分器逐项核对；每个神经模型的预测另由原始上游
模型类加载最终权重重新生成并核对。详见两份验证 JSON。

数据范围：GSE220974 的 49 种检出双干预标签，不是已经逐一核实
的设计构建体白名单。五折按物理基因对留出，模式变体同折。
28 个固定 A 的候选池选择任务，来自 8 个物理锚点；只有 5 个池
含至少 4 个候选。伙伴指标合计只有 4 个独立条件对比。
最优指当前候选池中 A 定义的 RNA 程序得分最高，不代表生长、
治疗效果或全图全局最优。此前其他模型的外层结果已经查看，
本轮 GEARS/Scouter 未据外层成绩改变已固定的选轮和训练规则。

RNA_MSE、partner_MSE 和 mean_regret 越低越好，hits 越高越好。
partner_MSE 消除自由端点线性主效应后评分；均值预测改善与
伙伴预测改善须分别解释。配对区间按物理锚点聚类，属于该
候选池的描述性不确定性，不是跨实验研究的独立验证。

完整重训需原续跑代码包、细胞数据包及最终权重索引。
本轻量包的 score.py 可独立复算已完成的全部结果。
GEARS 的 CPU 加速保持数学运算等价，FP32 存在已记录的舍入
差异；FP64 梯度和前向均验证。没有因此更换科学模型。
GEARS fold3/4 的六项由已授权 GPU 环境续跑：Python3.10.8、
Torch2.1.2+cu118、NumPy1.26.4；此前 CPU 环境为 Torch2.14.1。
保留 FP32、原划分、预算和内层选轮规则，未使用 AMP/TF32。
跨环境重放按预定容差核验，不宣称 CPU/GPU 训练逐位等价。
各项 WORKER_RUNTIME.json、执行指纹及适配源码保存在独立完成存档中。
'''
    (OUT / 'README_CN.txt').write_text(readme)
    manifest = {str(p.relative_to(OUT)): dict(bytes=p.stat().st_size, sha256=digest(p))
                for p in sorted(OUT.rglob('*')) if p.is_file() and p.name != 'MANIFEST.json'}
    (OUT / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    archive = OUT.with_suffix('.zip')
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file():
                z.write(p, str(p.relative_to(OUT)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
    print(json.dumps(dict(status='COMPLETE_VERIFIED', artifact=str(archive),
                          bytes=archive.stat().st_size, sha256=digest(archive))), flush=True)

if __name__ == '__main__':
    main()
