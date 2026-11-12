from pathlib import Path
import json,pandas as pd,platform,importlib.metadata as md,shutil
R=Path(__file__).resolve().parents[1];D=R/'results'
def table(df, digits=4):
 def val(v):
  return f'{v:.{digits}f}' if isinstance(v,float) else str(v)
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(val(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
tab=pd.read_csv(D/'PARTNER_GAIN_TABLE.csv');sm=pd.read_csv(D/'MODEL_COMPARISON.csv')
full=sm[(sm.readout=='raw_repeats')&(sm['repeat']=='pooled')][['source_model','arm','full_mse','partner_mse','full_gain_vs_M0_pct','partner_gain_vs_M0_pct']]
checks=json.loads((D/'INDEPENDENT_CHECK.json').read_text())
report=f'''# Han2017：现成模型的同面板比较

## 结论

上一轮确实使用已保存的Shen来源秩一功能修正。本轮完成秩一、谱模型、12参数SF12、V4神经伙伴头四类现成模型的比较，另保留M0和单扰动相加。没有重训任何模型，没有在Han上估计参数或选择模型。全部比较使用同一88个组合、14个查询端点和两个重建重复；另保留同88组合的已发布表型敏感性。

12参数SF12提供了有条件的伙伴收益：293T来源模型在重建重复汇总下相对同一M0降低伙伴误差2.2782%，两个重复分别4.5107%和0.6336%。但完整误差增加6.1869%，已发布表型的伙伴收益为−0.3029%。这是可追踪的模型差异，不是稳定外部突破。

## 伙伴增益主表

单位为百分比，正数代表伙伴误差降低，负数代表增加。每列均比较相同来源M0，列间不是同一分母。SF12和V4保留此前来源内部选择的幅度；秩一和谱模型保留此前选中版本。额外分量控制使用对应模型的来源内部预测计算系数。

{table(tab)}

A549的SF12为0，因为来源内部原选择就是零修正；本轮没有绕过选择去挑Han上更好的候选。V4每个来源固定平均两个已保存种子。

## 完整预测与绝对误差

{table(full,6)}

单扰动相加完整MSE为0.261449，低于本轮所有迁移后的完整模型。SF12在293T的伙伴改善没有补偿其端点部分损伤。因而本轮没有建立外部完整预测的竞争优势，不能只展示伙伴正数而隐去完整误差。

## 重复与读出敏感性

SF12（293T来源）伙伴收益：重复1为4.5107%，重复2为0.6336%，汇总2.2782%；同88对的已发布表型为−0.3029%。固定预测后逐一删除14个查询端点并重新投影，汇总伙伴收益在11/14次删除中为正，范围−4.9026%至7.4633%。删除分析高度重叠，不是14个独立实验，不报告二项显著性。

V4在A549为双重复同向、但幅度仅0.0465%的伙伴改善；HeLa与293T未形成双重复一致正向增量。谱模型原版本三组伙伴误差均增加；控制显著减轻损伤，其中293T重建读出有0.7860%伙伴收益，但完整误差仍增加8.4053%。全部重复、来源和控制版本均保留。

## 方法身份与输入处理

SF12是已保存的功能读出适配：原12特征构造、GenePT64输入及共享系数，不能称为旧RNA权重直接预测fitness。V4是已训练的功能伙伴头，不是新训练网络，也不是原RNA编码器零适配重放。秩一和谱模型均使用原保存矩阵及来源选择，没有重新运行优化。

四类模型共享同一来源M0、同一Han查询身份。来源基因知识坐标、源单扰动归一化和训练端点偏移固定。沿用上一轮单位转换，令h=1/mean(1/std(source single repeats))，Han单扰动映射为s_target*h/std(target all singles)，新增修正按std(target all singles)/h转回Han单位。这是固定推断适配，没有Han参数拟合。

SF12来源幅度beta可大于1，这是旧选择的实际结果；额外分量控制alpha限制在[0,1]，二者不能混称。SF12/V4的非零原生选择已在来源内部按伙伴目标调幅，所以本轮追加alpha接近1，不是额外取得一次控制成功。

没有用Han结果选择来源模型、正则、检查点、种子或幅度。published与raw_repeats用各自定义下的单扰动输入及结局；两种定义的比较并未隔离归一化、输入和结局变化的单独作用。

## 核验

来源预测重放误差为浮点舍入量级；原SF12特征映射逐值相同。上一轮Han的M0、秩一原预测、秩一分量控制和单扰动相加均重放一致，最大差1.78×10⁻¹⁵。来源文件散列未改变。

独立最小二乘残差化重新计算{checks['rows_recomputed']}个评分行，与SVD投影评分的MSE最大差{checks['max_abs_mse_error']:.3g}。分量控制的端点保持及绝对SSE恒等式检查通过。312行是模型、控制、读出和重复的评分行，不是312个独立实验。

本轮在先前Han结果已知之后开展，属于回顾性探索性模型比较。锁定用于保全过程和避免覆盖，不能追认为外部预注册或独立确认。没有修改v34和上一轮正式结果。

## 未纳入的模型

GEARS、Scouter、R-G及FullFeedback现有RNA检查点需要不同输入和输出，没有当前可直接重放的标量功能接口；强行列入会变成新建适配或重新训练。Jurkat原功能M0为基因索引模型，对当前Han面板的端点重叠不足，也不能构成同面板比较。本次未进行这些新训练或新接口开发。

## 复现

CPU即可，Python依赖见requirements.txt。

```bash
python code/predict.py --output results_replay
python code/score.py --output results_replay
```

第一条只读现成权重、来源内部归档和目标单扰动/身份，保存全部预测；第二条才读取目标结局。默认results含本次结果，不允许覆盖。

- results/MODEL_COMPARISON.csv：原生模型的完整/伙伴误差及相对M0收益。
- results/COMPARISONS.csv：所有控制与对照比较。
- results/SCORES.csv：全部312行原尺度成绩。
- results/SOURCE_SELECTION_TABLE.csv：实际来源配置、beta、alpha、gamma。
- results/ENDPOINT_DELETION.csv：逐端点删除敏感性。
- results/INDEPENDENT_CHECK.json：独立评分复核。
- pretrained/：实际现成参数、网络权重和来源内部预测。

既有科学代码权利边界沿用RIGHTS_AND_RELEASE.md，本包不自动扩大许可。
'''
(R/'REPORT_CN.md').write_text(report)
(R/'README.md').write_text('# Han2017 saved-model comparison\n\nSee REPORT_CN.md and PROTOCOL.json. No new fitting.\n\nRun `python code/predict.py --output results_replay` then `python code/score.py --output results_replay`.\n\nAll four families, all three sources, both phenotype definitions and both real repeats retained. Exploratory analysis after earlier Han outcomes were known.\n')
(R/'requirements.txt').write_text('\n'.join(f'{x}=={md.version(x)}' for x in ['numpy','scipy','pandas','torch'])+'\n')
(R/'ENVIRONMENT.json').write_text(json.dumps({'python':platform.python_version(),'platform':platform.platform()},indent=2))
print('REPORT_COMPLETE',len(report))
