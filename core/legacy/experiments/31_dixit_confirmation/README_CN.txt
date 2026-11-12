WitnessCell Dixit Fixed Reference — 2026-10-20

主要结果
在固定10个已测单扰动端点的空间内，G2M程序的伙伴校正可逐对交付；同一拟合
模型对同一组合的预测不随同行查询变化。当前开发数据总体SSE改善2.2671%，
去端点SSE改善14.2994%；总体收益的条件敏感性范围跨零。原全转录组主终点及
独立确认状态不变。完整论证见WitnessCell_Dixit_FixedReference_20261020_Report.txt。

立即运行单查询（不用原始RNA矩阵，不输入测试真值）
python model/predict_program.py --state model/deploy/Dixit_G2M_fold0.npz --pairs CREB1+YY1

单批多个查询
python model/predict_program.py --state model/deploy/Dixit_G2M_fold0.npz --pairs CREB1+YY1 E2F4+NR2C2

三个小模型对应原来的三个外层训练划分，合计约22KB。它们是回顾性开发评价
模型，未在全部45对上重训，也不是已独立确认的生产模型。程序尺度随训练折
固定；不能将不同拟合模型的分值差异解释为请求批次依赖。

完整离线复现（输出目录必须为空）
python REPRODUCE.py --output /tmp/witnesscell_fixed_reference_reproduced

需要Python 3和NumPy。绘图另外需要Matplotlib。实际执行包括两个模型分支、
纯代数检查、真实输出独立重构、条件敏感性和小模型导出。全程不访问网络或旧
工作区。保存逐数组哈希、结果表和精简输入的原始来源哈希。

复现用compact_assets内6.50MB原始条件均值与身份。保留双精度，训练局部尺度
从允许均值重算。未纳入大容量逐细胞数据及本轮完全未使用的半均值/方差。
冗余大FITTED_STATE文件由复现命令生成；EXPECTED_REPRODUCTION.json保留每个
数组哈希。model/predict_fixed.py可读取复现后完整状态；直接使用小模型则无需
先运行复现。原始运行源码和冻结记录保留，精简适配后的结果逐位相同。

目录
model/             固定参考主实现、各对照、单查询CLI和小模型
kernel_probe/      等价PSD核实现及完整残差阴性候选
theory_audit/      代数、闭式参考、独立预测重构、条件敏感性
compact_assets/    无损充分中间数据、冻结划分和预定义程序
source_helpers/    复用的原模型/划分/读取辅助源码

方法选择
正式主结果使用各自训练内选择lambda0.1/10/10。所有固定参数与未测池参考
结果都保留。未测池参考在本面板与外层查询池重合，不能用它重新宣称超越批次。
本轮没有改旧论文、旧原始预测或原生RNA主终点。
