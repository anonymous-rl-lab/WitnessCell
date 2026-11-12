from pathlib import Path
import sys,json,shutil
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT if (ROOT/'fitted').exists() else ROOT/'delivery/Han_Native_Refit_v1'

def md(columns,rows):
    return '\n'.join(['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
def main():
    p=pd.read_csv(ROOT/'results/POOLED.csv');part=pd.read_csv(ROOT/'results/BY_PARTITION.csv');rep=pd.read_csv(ROOT/'results/BY_REPEAT.csv');ctrl=pd.read_csv(ROOT/'results/CONTROL_COMPARISONS.csv');sel=pd.read_csv(ROOT/'results/SELECTIONS.csv')
    keys=[('native','rank1_original','原流程／主秩一'),('native','rank1_controlled','原流程／秩一＋分量控制'),('native','spectral_original','原流程／预设谱修正'),('native','spectral_controlled','原流程／谱修正＋分量控制'),('saturation','rank1_original','饱和参照／秩一'),('saturation','rank1_controlled','饱和参照／秩一＋分量控制'),('saturation','spectral_original','饱和参照／预设谱修正'),('saturation','spectral_controlled','饱和参照／谱修正＋分量控制')]
    rows=[];machine=[]
    for base,method,label in keys:
        vals=[]
        for panel in ('full_raw','full_published'):
            x=p[(p.panel==panel)&(p.baseline==base)&(p.method==method)].iloc[0];vals.extend([f'{x.partner_gain_vs_own_baseline_pct:.4f}%',f'{x.full_gain_vs_own_baseline_pct:.4f}%'])
            machine.append(dict(panel=panel,baseline=base,method=method,partner_gain_pct=x.partner_gain_vs_own_baseline_pct,full_gain_pct=x.full_gain_vs_own_baseline_pct))
        rows.append([label]+vals)
    table=md(['配置','重建：伙伴误差降低','重建：完整误差降低','发表：伙伴误差降低','发表：完整误差降低'],rows)
    stable=[]
    for base,method,label in keys:
        for panel,pn in [('full_raw','重建'),('full_published','发表')]:
            a=part[(part.panel==panel)&(part.baseline==base)&(part.method==method)];v=a.partner_gain_vs_own_baseline_pct
            stable.append([label,pn,f'{v.min():.4f}%—{v.max():.4f}%',f'{sum(v>0)}/5'])
    st=md(['配置','读出','五套划分的伙伴降幅范围','正向划分'],stable)
    repeatrows=[]
    for base,method,label in keys:
        a=rep[(rep.panel=='full_raw')&(rep.baseline==base)&(rep.method==method)].sort_values('repeat')
        repeatrows.append([label]+[f'{x:.4f}%' for x in a.partner_gain_vs_own_baseline_pct])
    rt=md(['配置','重复1：伙伴降幅','重复2：伙伴降幅'],repeatrows)
    pd.DataFrame(machine).to_csv(OUT/'results/KEY_RESULTS.csv',index=False)
    ind=json.loads((ROOT/'verification/INDEPENDENT_CHECKS.json').read_text());replay=json.loads((OUT/'replay_results/REPLAY_CHECKS.json').read_text())
    text=f'''# Han2017 原生重拟合：完成报告

日期：2026-11-09（北京时间）。状态：60/60 正式拟合批次、统一锁定评分及独立复算全部完成。

## 核心结论

本轮补到了独立研究数据上的伙伴学习正证据。复用冻结的 Jurkat 功能流程，在 Han 允许训练组合上重新估计参数后，主秩一模型在重建读出上的伙伴误差降低 **1.9105%**、完整误差降低 **9.0885%**；内部留出估计的分量控制进一步达到 **2.4431%**、**9.5212%**。发表聚合读出分别得到同方向的 **2.3017%→2.7695%** 伙伴改善。

更有力的归因补强来自预先指定的谱修正：在原功能基线上，重建／发表读出的伙伴误差降低 **5.6538%／7.8806%**；改用 v37 的饱和参照后，仍降低 **4.8245%／7.0357%**。因此，本轮的谱修正增量能够超出该非线性参照已经解释的变化。

这是在另一项研究中的冻结流程再评价。Han 先前已经接触过，现有划分也已经使用过；本次固定了新拟合规则和比较顺序、先锁预测再评分，但不据此改称完全未触碰的前瞻性确认。两种读出、两个重复及五套划分属于同一 Han 研究。

## 1. 实际执行的实验

- 重建读出：20,359 个组合，两个重建重复；每折训练 13,572—13,573 个组合，评价 6,786—6,787 个组合。保留输入中的 207 个单扰动身份，其中 206 个出现在组合面板。
- 发表读出：20,791 个组合，205 个端点，单列聚合响应；作为读出敏感性分析，不人为复制为两个生物重复。
- 两套基线：原 Jurkat 功能基线，以及 v37 饱和单效应参照。后者的选择与重拟合完全嵌入相应训练层，伙伴内层验证组合不参与基线拟合。
- 两种既有学习器：主秩一模型、预设次要谱模型。原惩罚网格、三内层折、零修正、完整误差选模、数值收敛要求均保留。
- 两基线 × 两读出 × 五划分 × 三外层折，共 60 个批次、2,880 次内层候选拟合和 120 次最终拟合。两次内层候选未达到原收敛阈值，按原规则排除；120 个最终拟合均收敛，没有超时任务。
- 每个批次保存基线、修正、内部预测、参数、收敛记录及时间戳；全部 60 批锁定后，统一执行外层评分。没有按外层胜负更改模型或网格。

## 2. 全部主要结果

表中每个降幅均相对**同一行自己的拟合基线**，正数表示误差降低。两种基线的百分比分母不同，不能直接把不同基线的百分比当成绝对误差比较。

{table}

解释需要区分三个层次：

1. **主分析成立。** 原流程的秩一模型在两种读出上均有正伙伴增量；完整响应改善也保留。
2. **超过饱和参照的证据以预设谱模型最清楚。** 谱修正在两种读出上都保留增量。饱和参照之上的秩一模型只在重建读出保留较小收益，发表读出的 15 折全部选择零修正。这一结果保留，不以谱模型替换原先指定的主模型身份。
3. **控制收益依配置而异。** 原流程的秩一模型经内部强度控制后进一步改善。谱模型多数内部强度已接近 1；饱和参照上的调幅略降了原谱模型收益。因此本轮支持学习与有条件的控制效用，不支持每一种配置都因控制而提高。

## 3. 五套划分与两个重复

{st}

重建数据分别按两次重复计算：

{rt}

五套划分重用同一批组合；这里报告的是划分稳定性，不把五套划分视为五项独立生物验证。重复共享部分初始及质粒对照，重复相关也不直接等同于噪声天花板。

## 4. 控制层具体增加了什么

在原功能基线的主秩一配置中，分量控制相对未调幅模型的伙伴误差再降低 **0.5429%（重建）／0.4788%（发表）**。这两个百分比以未调幅模型误差为分母，与表中相对基线的降幅区分。

相对同系数整体衰减，伙伴预测相同，而分量控制保留端点贡献，使完整误差额外降低 **2.3229%／3.1683%**。相对另行按内部完整误差优化系数的整体衰减，分量控制的完整误差仍降低 **0.5513%／0.4454%**，伙伴误差降低 **0.3386%／0.3623%**。这些结论属于该预设主配置，不外推到所有学习器。

完全关闭伙伴分量会丢掉这项非零增量：原流程主秩一控制后的伙伴降幅仍为 **2.4431%／2.7695%**。因此这项主结果包含保留有用伙伴预测的价值。

## 5. 对论文的直接意义

建议新增的准确结论是：

> 冻结的功能学习流程在 Han 的允许训练组合上重新估参后，能够改善留出组合的伙伴预测。预设谱修正在两种读出中均超出饱和参照；原秩一流程的训练侧分量控制进一步改善了伙伴误差，并保留了端点贡献。

本轮加强的是**方法流程在独立研究中的适用性**。此前冻结 Shen SF12 权重的迁移、随机方向检验与饱和参照归因仍保留原结论；新的本地参数拟合不能追溯证明旧来源权重已经携带可迁移知识。原有理论、正式榜单、Jurkat 和 RNA 证据的归属不变。

目前可以结束这一轮补实验，整合已有的主结果和完整敏感性表。无需继续为了更大百分比扩大模型搜索。

## 6. 复算与工程记录

- 首折原功能流程与原始 `run_budget` 程序的基线、秩一、谱预测逐值一致，最大差均为 0。
- 饱和基线与 v37 已保存的 30 折成绩一致，完整 SSE 最大差 4.55×10⁻¹³、伙伴 SSE 最大差 2.27×10⁻¹³。
- 另用 QR 最小二乘残差化独立核对全部 60 批：伙伴 SSE 最大差 {ind['sse_max_abs']['partner_sse']:.3g}，完整 SSE 最大差 {ind['sse_max_abs']['full_sse']:.3g}；内部强度最大差 {ind['max_alpha_abs']:.3g}。
- 从精简包中保存的模型参数重新生成全部预测、重估训练侧强度并重评分，1,080 条成绩全部通过；SSE 最大差 {max(replay['score_max_abs'].values()):.3g}，选中内层修正逐值重放差为 0。
- 首轮工程检查曾发现缓存文件名中的小数点导致候选／折次同名。在任何外层评分前作废该轮全部预测，修复后按不变科学协议重跑。故障及原冻结记录保存在 `engineering_notes/`，不混入正式结果。

完整依赖、输入、已拟合权重、全部候选选择分数、复算脚本与结果一并交付。完整重拟合和只重放已保存权重是两种明确区分的入口。
'''
    (OUT/'REPORT_CN.md').write_text(text)
    readme='''# Han native refit v1 — completed experiment

Start with REPORT_CN.md and results/KEY_RESULTS.csv. All 60 prespecified batches completed; the manuscript itself has not been overwritten.

## Replay saved fitted parameters (no fitting of predictor weights)

Python 3.12; install requirements.txt. From this directory run:

```bash
python verification/compact_replay.py
```

This reconstructs selected inner and final corrections from saved parameters, refits no predictor weights, recomputes the scalar control from allowed inner predictions, and verifies all 1,080 reference score rows. The expected precision is given in replay_results/REPLAY_CHECKS.json. Numerical differences in another BLAS environment are assessed in original squared-error units (absolute tolerance 1e-8); bitwise equality of every intermediate is not asserted across different environments.

## Refit the frozen scientific procedure

Work in a separate copy to preserve the delivered records:

```bash
python code/engineering_checks.py
python code/prepare.py
python code/launch.py --workers 6
python code/score.py
python verification/independent.py
```

Preparation regenerates outcome-free query inputs and training-only files, plus a new execution freeze record. Workers checkpoint every candidate and completed outer batch. Interrupted workers resume completed artifacts; scientific parameters do not change. The scorer refuses to run until all 60 prediction locks exist. Candidate convergence failures remain distinct from selecting zero correction.

## Contents and lossless reduction

- inputs/: unchanged Han panel inputs, outcomes, and original physical-pair partitions.
- code/: frozen fitting, selection and scoring scripts; unchanged native scientific engines under vendor/.
- fitted/: final and selected inner model parameters, exact selected inner predictions, baseline fits and control coefficients. Symmetric spectral matrices store the upper triangle losslessly; rank-one matrices store factors and strength, checked to reconstruct the original matrix exactly.
- results/: all comparisons, per-repeat/per-partition results, all candidate inner selection scores and convergence diagnostics.
- verification/: original-engine equality and independent residualization checks; compact model replay.
- engineering_notes/: the rejected pre-scoring cache-collision run's incident and original freeze. Invalid model outputs are not used as evidence.

The redundant unselected inner matrices and expanded copies of six algebraically derived query predictions are omitted. Full refitting recreates all candidates; saved selected weights recreate every reported outer score. Original LOCK.json files identify the full run before compaction, so some hashed expanded files are intentionally absent here. MANIFEST.json describes this deliverable; replay_results verifies its actual sufficient contents.

## Statistical and usage scope

The primary analysis is reconstructed readout, native baseline, rank-one model. Spectral models, saturation baselines, published readout and component control retain their prespecified secondary roles. Gains are relative to the named baseline. Repeated partitions do not count as independent studies. This is a frozen-procedure reanalysis of previously accessed Han data, not an untouched prospective validation. Published aggregate outcomes are one readout, not biological duplicates.

Scientific-source rights follow the author-provided WitnessCell code archive and its rights statement; this experiment does not assign a new blanket licence. NumPy, SciPy, pandas and optional plotting dependencies retain their own licences.
'''
    (OUT/'README.md').write_text(readme)
    # Standalone exact numerical figure. Partition dots are descriptive.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,2,figsize=(12,4.7),sharey=True)
    colors={'full_raw':'#257A9B','full_published':'#CC7B35'};methods=['rank1_original','rank1_controlled','spectral_original','spectral_controlled'];labels=['Rank one','Rank one\n+ control','Spectral','Spectral\n+ control']
    for ax,base,title in zip(axes,['native','saturation'],['a  Frozen native functional baseline','b  Saturating reference baseline']):
        for j,panel in enumerate(['full_raw','full_published']):
            x=np.arange(4)+(j-.5)*.34
            vals=[p[(p.panel==panel)&(p.baseline==base)&(p.method==m)].iloc[0].partner_gain_vs_own_baseline_pct for m in methods]
            ax.bar(x,vals,width=.31,color=colors[panel],alpha=.83,label='Reconstructed repeats' if j==0 else 'Published aggregate')
            for xx,m,v in zip(x,methods,vals):
                yy=part[(part.panel==panel)&(part.baseline==base)&(part.method==m)].sort_values('seed').partner_gain_vs_own_baseline_pct.to_numpy()
                ax.scatter(xx+np.linspace(-.07,.07,5),yy,s=12,color='#222222',zorder=4,linewidths=.3,edgecolors='white')
                ax.text(xx,max(v,float(yy.max()))+.16,f'{v:.2f}',ha='center',va='bottom',fontsize=8)
        ax.axhline(0,color='#333333',lw=.7);ax.set_xticks(range(4),labels);ax.set_title(title,loc='left',fontsize=11);ax.spines[['top','right']].set_visible(False);ax.set_ylim(-.1,9.2);ax.grid(axis='y',alpha=.14);ax.set_axisbelow(True)
    axes[0].set_ylabel('Partner squared-error reduction (%)');axes[0].legend(frameon=False,loc='upper left',fontsize=8)
    fig.text(.07,.02,'Bars: pooled error reductions relative to each panel’s own baseline. Dots: five overlapping partitions, not independent studies.',fontsize=8)
    fig.tight_layout(rect=[0,.08,1,1]);(OUT/'figures').mkdir(exist_ok=True)
    for ext in ('png','pdf','svg'):fig.savefig(OUT/'figures'/('Han_Native_Partner_Gain.'+ext),dpi=220,bbox_inches='tight')
    plt.close(fig)
    if Path(__file__).resolve()!=(OUT/'verification/report.py').resolve():shutil.copy2(__file__,OUT/'verification/report.py')
    print('Report and figures completed.')
if __name__=='__main__':main()
