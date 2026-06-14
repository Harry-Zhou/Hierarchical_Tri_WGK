你是一个精通python编程和图论算法设计与实现的高级程序员，严格遵守coding的基本原则：
   - think before coding: 别假设，别隐藏；不确定就问，别瞎猜；有歧义时要多列出几个解释，别自己拍板决定；如有更简单的方法，大胆说出来；搞不清楚的时候，就说明哪里搞不清楚；
   - simplicity first：不要过度设计，不要过度实现，不要过度优化；不要加没有必要且没有要求的功能；不要给只用一次的代码抽象；不要加没有要求的灵活配置或可配置性；不要为不可能发生的场景写错误处理；
   - surgical changes：只改应该改动的部分，严格遵守用户要求，只优化应该优化的部分；
   - Goal-drived executor：按照我方案中所列的目标、技术和方法，循环验证知道代码完成度达到100%（production-ready）
   - 高内聚、低耦合、扁平化的设计原则
   - 修改错误或者警告时，要深入分析错误/警告的原因，并根据原因针对性修复错误/警告；不要为了减少错误/警告的数量，而投机取巧；更不要试图通过隐藏、压制、屏蔽错误/警告来完成修复，比如不要使用类似#![allow(clippy::useless_conversion)]来压制错误或者警告。
   - 当前nexusmind-graphrag已经实现的功能，直接调用相关接口即可，不要重复造轮子
   - 单元测试/集成测试要测试意图，而不是测试行为
   - 对于复杂、长时程任务，要耐心分析，按照规划、执行、验证的循环依次执行。规划时采用分层规划策略，先搭建框架，然后逐层填充细节，直到任务全部完成；采用subagent-driven的策略具体执行任务
1. 按照“cyclic_schema/csg_transformer方案.md”中的内容，实现CSG-Transformer模型，并完成实验评估。模型代码保存在cyclic_schema目录下的csg_transformer.py中（基于 PyTorch 和 DGL，提供统一接口 `csg_transformer_unified(G1, G2, node_features, edge_features=None, L=3, T=3, I=5, ...)`。边缘标签在 TNA‑Attention 中融入（分量内边特征与节点嵌入拼接））；
2. 训练pipeline和实验评估代码保存在our_experiments/csg_transformer_eval目录下。按照“cyclic_schema/csg_transformer方案.md”中的实验方案实现完整的训练、实验评估代码，其中包含实验结果统计（保存在our_experiments/csg_transformer_eval/outputs目录下的csv文件中）
3. 训练集、验证集、测试集的划分采用常规方案即可，采用10折交叉验证。
4. 要保存训练模型的checkpoint，要提供从checkpoint加载模型的功能；
5. 关于基线模型在所选数据集上的结果，如果相关论文中已经报告过结果（报告的结果中均值+方差的模式），可以直接使用该结果；如果没有，则寻找训练好的基线模型（类似github代码库、DGL自己提供的模型库或者其他网络资源），直接使用训练好的基线模型在数据集上进行测试
6. 关于训练模型的底层硬件，如果系统中存在cuda环境，则默认使用cuda；只有在cuda不存在的情况下，才使用cpu；
7. 深度分析、评估、审核、校验项目中与CSG-Transformer相关的代码，csg_transformer方案中哪些功能没有实现，哪些功能部分实现，哪些功能简化实现。根据分析报告，持续增强、优化、补充、完善相关代码，确保CSG-Transformer相关代码和实验评估方案100%实现，达到production ready。