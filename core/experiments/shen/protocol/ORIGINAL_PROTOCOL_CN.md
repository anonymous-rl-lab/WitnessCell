# WitnessCell × Shen 2017：独立研究来源功能复现方案 v1

日期：2026-11-07。任务：设计实验与工程烟测；不修改 v31 的科学估计器。

## 0. 状态、问题与边界

本地方案在首次工程烟测前记录时间与SHA256，见 `protocol/LOCK_BEFORE_SMOKE.json`。这是本地设计记录，不是已经在OSF等公开平台登记的预注册。Shen真实工作簿尚未下载，真实来源字段映射、缺失值规则与测量资格仍待核实。本包不包含Shen真实预测成绩。

唯一主要问题：原v31功能流程的有符号秩一修正，能否在另一研究来源中，相对使用同一训练集重新拟合的M0，改善留出组合的端点消去响应预测？本次不测试RNA、TargetRank或普遍样本量门槛，不在Shen结果出现后更换主模型。

来源：Shen等，*Combinatorial CRISPR–Cas9 screens for de novo mapping of genetic interactions*，Nature Methods 14, 573–576 (2017)，doi:10.1038/nmeth.4225。

## 1. 为什么不能直接拿作者的GI表或全图单效应来训练

原文研究73个基因、三个细胞背景（HeLa、A549、293T）、两个生物学重复，分别在感染后3、14、21、28天测序。两个基因的3×3 guide组合形成九个物理构件，不能算成九个独立基因对。

原文Methods显示：构件fitness的单位为相对于双NTC的倍增/天；guide-level单效应由整个构件图的稳健回归估计；作者按这些效应对guide排序加权；最终互作分析还联合两次重复并进行Bayesian ensemble处理。这些全图估计不宜直接充当外层预测任务的单效应输入或独立双重复标签。已发表的显著NDEx边表也不能替代包含弱效应与阴性的连续测量。

本实验因此重新定义并预先固定测量对象：**在原测量时间序列中，逐物理构件、逐生物重复计算控制相对增长斜率，然后在合格且两个重复共有的构件上等权平均。** 其含义是实际guide-library平均响应，不是完美敲除的潜在gene-level因果效应，也不是精确复制作者的Bayesian GI估计量。

## 2. 输入恢复与质控（尚需源表字段审查）

规范化CSV字段：

`cell_line,replicate,day,construct_id,gene_a,guide_a,gene_b,guide_b,count`

这些是本项目中间格式，**不是对尚未读到的XLSX列名的猜测**。必须把真实工作表、列、repeat、NTC身份与源SHA写入映射文件后才能调用 `prepare_shen.py`。

每个 cell × repeat × time 的参考只由双NTC构件提供；以 log2(count) 减去该时点双NTC的平均log2(count)，对时间拟合斜率。共享测序深度因子在此差值中抵消；不以全部组合数据估计guide效应。单扰动输入只从 target+NTC 构件恢复。

预飞默认门：计数至少20、每构件至少3个有效时间点、每时点至少3个双NTC参考、每基因对至少3个两个重复共有的物理构件。缺失/低计数不加伪计数参与拟合。以上是**预先写下的分析者默认规则，不是论文原作者实际阈值**；原文说明其阈值随时点测序深度而变。读到源表的掩码与计数布局后，任何必要的测量适配应形成版本化修订，在第一次Shen模型拟合前重新锁定，不依据模型输赢选择阈值。

同一基因双guide、无单扰动支持端点、不可恢复重复身份等情况单独列出，不默默补值。必须报告每个细胞背景的构件数分布、缺失原因、实际基因对数、guide覆盖差异。不同组合保留的guide组成不同时，不能把构件平均差异无条件解释为纯gene-level互作。

## 3. 划分与烟测隔离

以规范化无序物理基因对为单位，**身份先划分，随后应用逐构件测量QC**。各guide构件、顺反方向、两重复及三个细胞背景中的同一基因对使用同一划分。哈希盐固定在 `PROTOCOL.json`。

完整73基因面板名义上共有2628对：80%训练2102对，20%查询526对；这不是2628个训练对。源数据资格过滤后数量可以下降，不补抽查询，不选择更容易预测的分割。

真实来源烟测只在HeLa的外层训练池内进行：按另一个固定身份哈希取约10%作为开发验证。完整面板名义上为1891个烟测拟合对、211个烟测验证对，正式526个查询不评分。A549和293T不用于方法选择。

正式三背景拟合完成、预测文件全部写入SHA冻结清单之后，才允许统一评分。测量重建可能读取原始查询计数以生成评价标签，但预测器不读取它们；这是软件访问隔离，不是第三方盲法，也不是加密保证。

## 4. 方法保持不变

- M0：两个重复分别拟合，单效应标准化，端点/斜率/曲率的零和偏差与原ridge=10/s保留。
- 主要方法：原有符号秩一修正，G=lambda uu^T，lambda正负均可。
- 次要方法：原谱收缩。
- 配套参照：单扰动响应相加。
- 与v31完全相同的三个内折、tau/tau0=[.3,.1,.03,.01]、nu=[.1,1]、零修正候选；六起点、3000迭代上限、1e-6容差。
- 内层仍按两个重复合并的**完整原单位SSE**选参，不改成伙伴损失，不因外层效果选秩一或谱。
- 各内部子集都重拟M0；新细胞背景重新估参，不直接迁移Horlbeck权重。

原生模块逐文件SHA见 `audit/V31_SOURCE_HASHES.json`，本包不授予其新的软件许可；保留来源许可文件和作者权利边界。

## 5. 主要评价与结果解释

对每个背景的正式查询构造端点发生矩阵X，令P=I-XX+。主要值为

`gain_partner = 1 - sum_r ||P(Y_r - M0_r - D)||² / sum_r ||P(Y_r - M0_r)||²`。

先合并两个重复的SSE，再取比值。配套报告完整误差、两个重复的各自增量、全部方法和零修正选择。投影维度实算；没有可评价对比或基线能量近零时记为不可评价，不伪造增益。

跨背景摘要取三个背景相对增量的等权平均，并保留逐背景结果；不混合不同背景的原尺度SSE。一个新研究里的三个细胞背景不是三个独立研究，九个guide构件也不是九个独立基因重复。

**没有“2/3为正就通过”的规则。** 看效应量、重复一致性、完整误差是否恶化及条件敏感性。不能预先保证成功；零修正属于有效的阴性结果。

采用固定预测、固定原始查询投影P的10,000次端点重加权条件敏感性：每个物理基因独立抽取Exp(1)权重，每条配对取两端权重乘积，在投影后的逐配对误差上重新计算比值；同一基因在三个背景共享同次权重。报告每背景及等权背景均值的中心95%区间，明确不是总体置信保证。固定种子20261107，不重拟模型或投影。实现已加入 `score_prepared.py` 并在合成夹具上验证，尚未对真实Shen执行。该统计细节在原工程烟测后、读取任何Shen工作簿前写入补充锁定记录。

cos(theta)和oracle beta只在冻结后诊断，绝不据其修改正式预测。两档或多个预算也不作为本实验新的变量；本次只做高支持复现，不寻找普遍阈值。

## 6. 停止与解释规则

若原始表不能恢复单独重复、NTC或合法单扰动输入，停止来源确认，不将 pooled GI 拆成两份。若不收敛，保留失败并记录，不挑种子或换模型。若真实烟测发现实现错误，可修工程问题并重新锁定，但不可根据其胜负调整方法再称整项完全未开发。

一个合格阳性支持功能伙伴修正的跨研究来源用途，不自动证明RNA互作恢复。阴性说明当前冻结构造没有在此来源获得收益，不推翻Horlbeck阳性，也不证明所有方法均不可能。

## 7. 执行入口

1. `python code/fetch_sources.py`：获取官方S1/S2并检查真实ZIP结构；本环境当前DNS失败，日志保留。
2. `python code/inspect_source_workbook.py --help`：只查看工作簿结构，尚未对真实Shen表运行。
3. 将真实字段审计映射成规范化CSV和 provenance JSON；不能直接拿源XLSX当CSV输入。
4. `python code/prepare_shen.py --counts raw_counts_canonical.csv --provenance SOURCE_MAPPING.json --out prepared`
5. `python code/fit_prepared.py --prepared prepared --out source_smoke_fit --stage smoke`
6. `python code/score_prepared.py --prepared prepared --fits source_smoke_fit --out source_smoke_score`
7. 源资格/协议锁定后，`python code/fit_prepared.py --prepared prepared --out formal_fit --stage formal`
8. 完整冻结清单生成后，`python code/score_prepared.py --prepared prepared --fits formal_fit --out formal_score`

上述完整准备/拟合/评分协调入口已经在明确标注的合成三背景软件夹具上运行验证，尚未在Shen上运行。

## 8. 来源

- 原文：https://www.nature.com/articles/nmeth.4225
- 作者全文：https://idekerlab.ucsd.edu/wp-content/uploads/publications/Shen_NatureMethods2017.pdf
- Supplementary Table 2：https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4225/MediaObjects/41592_2017_BFnmeth4225_MOESM42_ESM.xlsx
- Supplementary Table 1：https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.4225/MediaObjects/41592_2017_BFnmeth4225_MOESM41_ESM.xlsx

本方案的QC数值、拆分哈希与结果组织是新提议，不能归因于原论文；原生科学估计器与超参数规则来自用户提供的v31归档。
