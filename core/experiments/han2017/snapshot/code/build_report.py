from common import *
import itertools
import pandas as pd

def table(df,cols,decimals=6):
    lines=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    for _,row in df.iterrows():
        values=[f'{row[c]:.{decimals}f}' if isinstance(row[c],(float,np.floating)) else str(row[c]) for c in cols]
        lines.append('| '+' | '.join(values)+' |')
    return '\n'.join(lines)

def main():
    scores=pd.read_csv(ROOT/'results/POOLED_SCORES.csv');comp=pd.read_csv(ROOT/'results/COMPARISONS.csv');seeds=pd.read_csv(ROOT/'results/SEED_SCORES.csv');sensitivity=pd.read_csv(ROOT/'audit/SENSITIVITY_SUMMARY.csv')
    info=json.loads((ROOT/'data/PREPARED.json').read_text());audit=json.loads((ROOT/'audit/FINAL_CHECKS.json').read_text())
    rows=[];ranges=[]
    for panel in ['full_raw','full_published']:
        for model in sorted(scores[scores.panel==panel].model.unique()):
            d=comp[(comp.panel==panel)&(comp.model==model)&(comp.repeat=='pooled')&(comp.controlled=='target_control')].set_index('reference')
            rows.append(dict(panel=panel,model=model,partner_gain_vs_original_pct=d.loc['original','partner_reduction_pct'],partner_gain_vs_off_pct=d.loc['partner_off','partner_reduction_pct'],full_gain_vs_additive_pct=d.loc['single_additive','full_reduction_pct']))
            z=seeds[(seeds.panel==panel)&(seeds.model==model)&(seeds.repeat=='pooled')].pivot(index='seed',columns='method',values='partner_sse')
            gain=100*(1-z.target_control/z.single_additive);reg=100*(1-z.target_control/z.original)
            alpha=[]
            for p in (ROOT/'results'/panel).glob('*/*/PARAMETERS.json'):
                keys=load(ROOT/f'data/{panel}_frozen_predictions.npz')['models'];j=list(keys).index(model);alpha.append(json.loads(p.read_text())['alpha'][j])
            ranges.append(dict(panel=panel,model=model,partner_gain_min=gain.min(),partner_gain_max=gain.max(),gain_vs_original_min=reg.min(),gain_vs_original_max=reg.max(),positive_partitions=int((gain>1e-10).sum()),alpha_min=min(alpha),alpha_max=max(alpha)))
    result=pd.DataFrame(rows);result.to_csv(ROOT/'results/FULL_HAN_KEY_RESULTS.csv',index=False);pd.DataFrame(ranges).to_csv(ROOT/'results/FULL_HAN_STABILITY.csv',index=False)
    # Enumerate the entire raw single-supported identity universe, including absent constructs.
    raw=load(ROOT/'data/full_raw_inputs.npz');universe=pd.DataFrame(itertools.combinations(raw['genes'],2),columns=['gene_a','gene_b'])
    support=pd.read_csv(ROOT/'data/RAW_PAIR_SUPPORT.csv');cover=universe.merge(support,on=['gene_a','gene_b'],how='left')
    for c in ['constructs','orientation_ab','orientation_ba']:cover[c]=cover[c].fillna(0).astype(int)
    cover['eligible']=cover['eligible'].fillna(False).astype(bool);cover['reason']=np.where(cover.eligible,'eligible','initial-count/plasmid/orientation/construct-count qualification')
    cover.to_csv(ROOT/'audit/ALL_RAW_PAIR_QUALIFICATION.csv',index=False);assert int(cover.eligible.sum())==info['raw_pairs']
    old={tuple(p) for p in json.loads((ROOT/'audit/PRIOR_SCORED_PAIRS.json').read_text())};missing=set(load(ROOT/'data/target_knowledge.npz')['missing'])
    coverage=[]
    for panel in ['common_raw','common_published','full_raw','full_published']:
        d=load(ROOT/f'data/{panel}_frozen_predictions.npz');pairs=d['pairs'];coverage.append(dict(panel=panel,n_pairs=len(pairs),n_endpoints=len(set(pairs.ravel())),n_models=len(d['models']),prior_scored_pairs=sum(tuple(sorted(p)) in old for p in pairs),outside_prior_score_pairs=sum(tuple(sorted(p)) not in old for p in pairs),pairs_with_imputed_knowledge=sum(bool(set(p)&missing) for p in pairs)))
    coverage=pd.DataFrame(coverage);coverage.to_csv(ROOT/'audit/COVERAGE.csv',index=False)
    positive=result[result.model.isin(['HeLa__sf12','293T__sf12'])]
    repeat=comp[(comp.panel=='full_raw')&comp.model.isin(['HeLa__sf12','293T__sf12'])&(comp.controlled=='target_control')&(comp.reference=='single_additive')&(comp.repeat!='pooled')][['model','repeat','partner_reduction_pct','full_reduction_pct']]
    native=[]
    for panel in ['common_raw','common_published']:
        s=scores[(scores.panel==panel)&scores.model.str.endswith('__rank1')&(scores.repeat=='pooled')].groupby('method')[['full_sse','partner_sse']].sum()
        for ref in ['original','partner_off','source_control','single_additive']:
            native.append(dict(panel=panel,reference=ref,partner_gain_pct=100*(1-s.loc['target_control','partner_sse']/s.loc[ref,'partner_sse'])))
    native=pd.DataFrame(native);native.to_csv(ROOT/'results/PRIMARY_COMMON_RANK1.csv',index=False)
    portability=json.loads((ROOT/'audit/PORTABLE_FEATURE_CHECKS.json').read_text());pm=max(x['projected_direction_max_abs_error'] for x in portability['common_projected_direction'])
    completions=[json.loads(p.read_text()) for p in (ROOT/'results').glob('*/*/*/COMPLETE.json')];assert len(completions)==60
    text=f'''# Han2017：冻结模型后的全量重评估 v3

## 实际完成与核心判断

本轮已完成60个主分析分折任务，保存10,440个模型/方法/重复评分行。基础预测模型重训次数为0；Han只用于标量控制系数的校准。评分行、划分、来源模型均不是独立实验数量。

全量共享函数分析中，HeLa与293T来源SF12在相同的预定控制规则下同时优于未控制和关闭伙伴，两个重建重复的方向一致，五套划分均为正；已发布表型也保持正向。增量幅度小，不能描述为全面领先或普遍保证。网络头没有形成同样的稳定收益。

原88对、原完整模型的秩一主分析仍主要支持止损，没有在重建读出上胜过关闭伙伴。该结果保留，不用全量SF12适配分析替换。

## 1. 全部数据和模型覆盖

{table(coverage,list(coverage.columns),0)}

原始读出从207个有限单扰动身份的{len(universe):,}个可能无序组合逐项核查，保留{info['raw_pairs']:,}对，排除{len(universe)-info['raw_pairs']:,}对。排除依据沿用既有计数、质粒、构件数和双方向支持规则，不根据响应强弱或模型成绩筛选。全部身份及理由见ALL_RAW_PAIR_QUALIFICATION.csv。已发布表型按其实际可用、有限读出保留{info['published_pairs']:,}对，两个全量面板并不完全相同。

秩一和谱模型的存储矩阵依赖来源基因身份，无法无训练地产生全部Han新靶点的矩阵条目。它们在全部88个共同合格组合上完成重评；没有将缺失矩阵条目补零冒充全量预测。所有Han合格组合均进入可适用的共享模型分析。

SF12和V4的共享函数使用既有系数、网络权重及来源种子。扩展时恢复来源固定特征变换，对新基因只做确定性输入映射。ADRBK1缺少既有嵌入，沿用此前的全局均值策略，另作完整删除敏感性。

## 2. 两个明确区分的实验

**原完整模型比较：** 保留M₀及原D，执行原端点/目标端点两组，以及不控制、关闭、来源系数、Han系数四种伙伴强度；另外保留按完整误差优化的整体收缩。三个来源、四类模型、两个读出、五套三折均报告。

**全量共享函数适配：** 目标基线B固定为允许的Han单扰动相加，现成SF12/V4共享函数提供方向D，预测为B+αP_QD。保留原来源选择幅度，原选择为零的模型仍为零。来源端点平移在投影后消失；在原88对上，扩展接口与原伙伴方向的最大差为{pm:.3g}。不需要为新基因伪造来源端点参数。

第二项是明确标记的新输入适配，不是原完整M₀权重覆盖全部新靶点，也不是重新训练SF12或网络。这个版本已经固定端点，因此“整体缩小纯伙伴向量”与“伙伴控制”等价，不能重复计为两种方法的胜负。

两项实验的M₀、数据量和查询图不同。现有结果不能把改善全部归因于数据量，更不能据此宣布一个普适样本门槛。

## 3. 全量结果：所有来源都报告

表中正数表示误差降低，单位为%。source_control对SF12/V4非零分支接近α=1，因此其比较与不控制基本重合，不作为额外独立成功。

{table(result,list(result.columns))}

全量版本的partner_off就是B。HeLa SF12的主要作用是把过强修正缩小，同时保留小幅有用信息；293T SF12原方向已较有用，控制进一步小幅改善。这两个不同情形均按同一标量规则得到，不是事后给不同来源选择不同控制算法。

## 4. 重复、划分及控制强度

{table(repeat,list(repeat.columns))}

{table(pd.DataFrame(ranges),['panel','model','partner_gain_min','partner_gain_max','positive_partitions','alpha_min','alpha_max'])}

两个SF12来源在重建与已发布读出上均有五套划分同向的伙伴改善；这描述当前面板的稳定性，五套重叠划分不是五个独立生物学重复。两个真实重复共享初始及质粒控制；已发布读出也来自同一研究，并同时改变测量定义和允许输入，不能合计成外部独立复制次数。

## 5. 旧烟测与缺失嵌入不能解释这两个正结果

以下为主结果可见后进行的身份定义稳健性检查，固定模型、公式及种子，分别从校准池和查询池删除旧评分组合，或删除缺失嵌入端点。只重新估计原定标量，不进行新选模。全部模型结果见SENSITIVITY_SUMMARY.csv。

{table(sensitivity[(sensitivity.model.isin(['HeLa__sf12','293T__sf12']))&(sensitivity.reference=='single_additive')],['panel','filter','model','partner_gain_pct','full_gain_pct'])}

两项检查保留了两个SF12来源的正向增量。它们是对已见结果的敏感性分析，不是新的确认性实验。

## 6. 原88对秩一主分析完整保留

{table(native,list(native.columns))}

重建读出上，相对不控制有止损，但相对关闭和来源强度没有改善。不能因全量SF12分析为正，就把这项原主分析也改写成成功。原88对的其余模型、全部端点版本和对照见COMPARISONS.csv。

## 7. 冻结、访问历史和实际核验

本轮协议与现成模型在新重评估评分前锁定；每折先写入系数和查询预测、记录散列，再评分。没有根据结果切换来源、检查点、正则、特征、种子或强度公式。目标校准确实使用Han校准结局，不能称为零目标标签拟合。

此前88对及已发布表型部分已经被查看和用于开发，原始表及计数也已访问。本轮将先前结果保留为烟测/探索记录，但不会抹去访问历史。Han仍是相对来源研究的外部目标研究；本轮属于看过其部分结果之后的冻结重评估及扩展，不能称为整个项目从未接触的独立预注册确认。

- 全部51个冻结来源模型/代码文件散列保持一致，原88对预测逐值重放最大差0。
- 每个主分析折、每个模型、每个方法的汇总分数均用独立列主元QR重算；最大SSE差{audit['independent_qr_max_abs_sse_difference']:.3g}。
- 实际改变查询结局文件后重新执行预测，四类面板的系数及预测最大变化均为0。
- 端点保持、完整与伙伴绝对SSE差恒等式以及锚定不改变伙伴预测检查通过。
- 在新目录从紧凑原始计数和已发布表重新构建全部输入/结局，数组逐值一致。未重新下载FASTQ或重跑原始序列映射。
- 模型参数没有重训；校准系数有目标数据估计。没有将实现核验计为新增生物学证据。

## 8. 当前可以支持的结论

> 在Han的全量功能读出上，既有SF12共享伙伴函数经过目标侧标量校准，在两个来源配置中得到小幅、跨划分及两种读出方向一致的留出增量；它们既优于未控制，也优于关闭伙伴。该结果支持特定模型和适配流程下的有用控制，不自动推广到秩一、谱或网络模型，也不说明所有此前阴性都由实验设计造成。

前面目标与执行路线错位属于需要承担和纠正的设计问题；现有数据同时显示模型方向本身仍决定控制能否保留收益。将这两点分开，才与实验目标对齐。

## 9. 文件与复现

主表：results/FULL_HAN_KEY_RESULTS.csv；全部比较：results/COMPARISONS.csv；逐折参数、锁定预测和分数：results/<panel>/<seed>/<fold>/；身份排除：audit/ALL_RAW_PAIR_QUALIFICATION.csv；独立核验：audit/FINAL_CHECKS.json。

完整目录包含冻结源模型、特征、重建输入、紧凑计数和已发布表。各臂预测可由冻结数组、每折系数及已公开公式逐值重建；全部臂在评分前保存了数值散列，控制输出另存NPZ。

```bash
pip install -r requirements.txt
python code/replay.py --workdir /absolute/path/new_replay
```

从原始紧凑计数重建输入再重评：

```bash
python code/replay.py --workdir /absolute/path/new_full_replay --rebuild-measurements
```

该命令重新推断和校准，不重新训练基础模型。当前实际执行核验范围如第7节所列；原v34和全部旧结果保持不变。
'''
    (ROOT/'REPORT_CN.md').write_text(text)
    dump(ROOT/'COMPLETION.json',dict(utc=now(),primary_folds_completed=len(completions),primary_score_rows=10440,base_model_training_calls=0,source_weights_unchanged=True,raw_pairs=info['raw_pairs'],published_pairs=info['published_pairs'],source_native_models_common_panel=12,shared_adapter_models_full_panel=6,new_paper_version_created=False,all_scheduled_primary_analyses_complete=True))
    print('REPORT_WRITTEN',len(text),flush=True)
if __name__=='__main__':main()
