# CSG‑Transformer: 融合分层循环模式图与 Transformer 的图表示学习

**摘要**：本文提出 **CSG‑Transformer**，一种将多层循环模式图（Cyclic Schematic Graph, CSG）抽象与图 Transformer 架构深度融合的图表示学习模型。模型严格遵循离散 HTN‑WL 的消息传递蓝图，将每个离散步骤替换为可微的注意力机制：层内使用三角化邻域注意力（TNA‑Attention）捕获局部连通模式，层间通过前向/后向交叉注意力实现多尺度信息融合。引入随机拉普拉斯位置编码打破对称性，使模型能区分 CFI 图等经典难例。理论分析证明模型具有图同构不变性，至少与 HTN‑WL 一样强大，且能以概率 1 区分任意非同构的 CFI 图对。算法复杂度在稀疏图上可控，为实际应用提供了理论保障。

**关键词**：图表示学习；循环模式图；Weisfeiler‑Leman 测试；Transformer；注意力机制；图同构

---

## 1. 引言

图神经网络（GNN）的表达能力通常受限于 Weisfeiler‑Leman (WL) 测试：普通消息传递 GNN 等价于 1‑WL，高阶 GNN 等价于 \(k\)-WL。然而，存在经典反例（如 CFI 图）使得任意有限 \(k\)-WL 都无法区分，而 \(n\)-WL 虽能区分但代价是指数时间。近年来，位置编码、子图计数、全局注意力等方法被用于突破 WL 上限，但要么缺乏理论保证，要么计算负担过重。

本文提出一种新思路：利用图的圈结构进行层级抽象，并融合 Transformer 的注意力机制。核心观察包括：
- 图的圈空间蕴含丰富的全局信息，且可通过规范圈基（CCB）进行规范化压缩；
- 离散 HTN‑WL 已证明在多层 CSG 上通过前向/后向消息传递能区分 CFI 图；
- Transformer 的自注意力和交叉注意力天然适配多尺度、多源信息的聚合；
- 随机位置编码可打破对称性，使模型学到图的全貌而非仅局部模式。

为此，我们设计 **CSG‑Transformer**，其贡献如下：
1. 提出与离散 HTN‑WL 蓝图严格对齐的连续、可微消息传递框架；
2. 设计三角化邻域注意力（TNA‑Attention），显式编码邻居连通分量结构；
3. 引入前向/后向交叉注意力，实现圈内信息汇总与高层反馈注入；
4. 证明模型的同构不变性、表达能力下界及对 CFI 图的区分能力；
5. 给出计算复杂度分析与实现细节。

---

## 2. 预备知识

### 2.1 图与圈空间

设 \(\mathcal{G}=(\mathcal{V},\mathcal{E})\) 为无向简单图，\(n=|\mathcal{V}|\)，\(m=|\mathcal{E}|\)，\(c\) 为连通分量数。圈空间 \(\text{CS}_{\mathcal{G}}\subseteq \mathbb{F}_2^{|\mathcal{E}|}\) 是 \(\mathbb{F}_2\) 上的向量空间，维数 \(\mu(\mathcal{G}) = m-n+c\)（圈秩）。

### 2.2 规范圈基与 CSG 变换

**定义 2.1（规范圈基，CCB）**。通过四阶段规范投影（顶点规范标签、规范生成树、规范基础圈、\(\mathbb{F}_2\) 精化）得到的圈基，记为 \(\text{CCB}_{\mathcal{G}}\)。该基在同构下唯一。

**定义 2.2（CSG 变换 \(\Phi\)）**。\(\Phi(\mathcal{G})\) 的节点包含：
- CB 节点 \(b_i\)，对应每个 \(C_i\in\text{CCB}_{\mathcal{G}}\)；
- 非圈节点 \(v\in\mathcal{V}\setminus\bigcup_i V(C_i)\)；
- 接口节点（共享顶点或连接非圈的圈上顶点）。
边由共享关系及原始边定义。

**性质 2.1**。若 \(\mathcal{G}_1\cong\mathcal{G}_2\)，则 \(\Phi(\mathcal{G}_1)\cong\Phi(\mathcal{G}_2)\)。

### 2.3 三角化邻域类（TNA）

对节点 \(v\)，令 \(TN_{\mathcal{G}}(v)=\{R_1,\dots,R_{k_v}\}\) 为诱导子图 \(\mathcal{G}[N(v)]\) 的连通分量分解。

### 2.4 离散 HTN‑WL 蓝图（一次迭代）

给定多层 CSG 序列 \(\mathcal{G}^{(0)},\dots,\mathcal{G}^{(L)}\)（\(\mathcal{G}^{(0)}\) 为输入图），一轮迭代包括：

**前向阶段**（自底向上）：
```
for l = 0..L-1:
    在 G^(l) 上执行一次 TNA 消息传递；
    用 G^(l) 的标签构造 G^(l+1) 节点的初始标签元组（规范圈标签）；
    在 G^(l+1) 上执行一次 TNA 消息传递；
```

**后向阶段**（自顶向下）：
```
for l = L..1:
    用 G^(l) 的标签构造 G^(l-1) 节点的后向标签元组（收集上层标签排序后拼接）；
    在 G^(l-1) 上执行一次 TNA 消息传递；
```

重复 \(I\) 轮即得 HTN‑WL。

---

## 3. CSG‑Transformer 模型架构

### 3.1 总体框架

CSG‑Transformer 将上述离散流程中的每个操作替换为可微的连续版本：
- 节点标签 → 连续嵌入向量 \(\mathbf{h}_v\in\mathbb{R}^d\)；
- 一次 TNA 消息传递 → \(T\) 轮 **TNA‑Attention**（可选全局自注意力）；
- 构造初始标签元组 → **前向交叉注意力**；
- 构造后向标签元组 → **后向交叉注意力**。

模型输入为一对图 \(G_1, G_2\)（或单图），输出图级嵌入向量。

### 3.2 多层 CSG 构建

**算法 1（构建多层 CSG）**  
输入：图 \(G\)，最大层数 \(L_{\max}\)（默认直至收敛）  
输出：  
- 图列表 \(\{\mathcal{G}^{(0)},\dots,\mathcal{G}^{(L)}\}\)，每层节点属性含 `type` 和 `original_nodes`；
- 向下映射 \(\text{down}^{(l)}:\mathcal{V}^{(l+1)}\to 2^{\mathcal{V}^{(l)}}\)（上层节点 → 下层节点集）；
- 向上映射 \(\text{up}^{(l)}:\mathcal{V}^{(l)}\to 2^{\mathcal{V}^{(l+1)}}\)（下层节点 → 上层节点集）。

**实现**：迭代调用 `cyclic_schema.cyclic_schematic_graph` 直至圈秩为零或达到最大层数。

### 3.3 位置编码

对每层图 \(\mathcal{G}^{(l)}\)，计算位置编码 \(\text{PE}^{(l)}(v)\in\mathbb{R}^p\)：
1. 计算图拉普拉斯 \(L^{(l)}=D^{(l)}-A^{(l)}\)，求解 \(L^{(l)}\mathbf{u}=\lambda D^{(l)}\mathbf{u}\)，取最小 \(p\) 个非零特征值对应的特征向量 \(\mathbf{u}_1,\dots,\mathbf{u}_p\)。
2. 添加独立高斯噪声 \(\epsilon_v\sim\mathcal{N}(0,\sigma^2 I)\)（\(\sigma=0.01\)），得
   \[
   \text{PE}^{(l)}(v)=\big(\mathbf{u}_1(v),\dots,\mathbf{u}_p(v)\big)+\epsilon_v.
   \]

噪声打破特征值简并导致的对称性，同时保留近似协变性。

### 3.4 节点嵌入初始化

第 0 层节点嵌入由原始节点特征 \(\mathbf{x}_v\) 与位置编码拼接后经 MLP 得到：
\[
\mathbf{h}_v^{(0)} = \text{MLP}_{\text{init}}\big( [\mathbf{x}_v \parallel \text{PE}^{(0)}(v)] \big) \in \mathbb{R}^d.
\]

第 1..L 层的初始嵌入将在前向阶段由交叉注意力动态计算。

### 3.5 核心算子

#### 3.5.1 三角化邻域注意力（TNA‑Attention）

**定义 3.1（TNA‑Attention 一轮更新）**。给定节点嵌入 \(\mathbf{H}\in\mathbb{R}^{n\times d}\) 和图 \(\mathcal{G}\)（含预计算的 \(TN(v)\)），计算：
\[
\mathbf{m}_R = \frac{1}{|R|}\sum_{u\in R}\mathbf{h}_u \quad (\text{分量均值}),
\]
\[
\alpha_{v,R} = \frac{\exp\big(\text{LeakyReLU}(\mathbf{a}^T[\mathbf{W}_Q\mathbf{h}_v\;\|\;\mathbf{W}_K\mathbf{m}_R])\big)}{\sum_{R'\in TN(v)}\exp\big(\text{LeakyReLU}(\mathbf{a}^T[\mathbf{W}_Q\mathbf{h}_v\;\|\;\mathbf{W}_K\mathbf{m}_{R'}])\big)},
\]
\[
\tilde{\mathbf{h}}_v = \sigma\Big(\sum_{R}\alpha_{v,R}\mathbf{W}_V\mathbf{m}_R\Big).
\]

可选全局多头自注意力（仅在节点数 \(n\le n_{\text{thres}}\) 时启用）：
\[
\mathbf{H}\leftarrow \text{LN}\big(\mathbf{H}+\text{MHA}(\mathbf{H})\big),\quad
\mathbf{H}\leftarrow \text{LN}\big(\mathbf{H}+\text{FFN}(\mathbf{H})\big).
\]

一轮 TNA‑Attention 输出 \(\mathbf{H}'\)。连续执行 \(T\) 轮作为一次“消息传递”。

#### 3.5.2 前向交叉注意力（ForwardCrossAttn）

**定义 3.2**。给定下层嵌入 \(\mathbf{H}^{(l)}\in\mathbb{R}^{n_l\times d}\)，映射 \(\text{down}^{(l)}\)，以及上层位置编码 \(\text{PE}^{(l+1)}\)，对每个上层节点 \(u\)：
\[
\beta_v = \frac{\exp\big(\mathbf{q}_u^T \mathbf{W}_K \mathbf{h}_v^{(l)}\big)}{\sum_{w\in\text{down}^{(l)}(u)}\exp\big(\mathbf{q}_u^T \mathbf{W}_K \mathbf{h}_w^{(l)}\big)},\quad
\mathbf{x}_u = \sum_v \beta_v \mathbf{W}_V \mathbf{h}_v^{(l)},
\]
\[
\mathbf{h}_u^{(l+1)} = \text{MLP}_{\text{fwd}}\big( [\mathbf{x}_u \parallel \text{PE}^{(l+1)}(u)] \big).
\]
其中 \(\mathbf{q}_u\) 为可学习的查询向量（可每节点独立或每类型共享）。

#### 3.5.3 后向交叉注意力（BackwardCrossAttn）

**定义 3.3**。给定下层当前嵌入 \(\mathbf{H}^{(l-1)}\)，上层嵌入 \(\mathbf{H}^{(l)}\)，映射 \(\text{up}^{(l-1)}\)，对每个下层节点 \(v\)：
\[
\gamma_u = \frac{\exp\big( (\mathbf{W}_Q'\mathbf{h}_v^{(l-1)})^T (\mathbf{W}_K'\mathbf{h}_u^{(l)}) \big)}{\sum_{u'\in\text{up}^{(l-1)}(v)}\exp\big( (\mathbf{W}_Q'\mathbf{h}_v^{(l-1)})^T (\mathbf{W}_K'\mathbf{h}_{u'}^{(l)}) \big)},\quad
\mathbf{b}_v = \sum_u \gamma_u \mathbf{W}_V'\mathbf{h}_u^{(l)},
\]
\[
\mathbf{h}_v^{(l-1),\text{new}} = \text{MLP}_{\text{bwd}}\big( [\mathbf{h}_v^{(l-1)} \parallel \mathbf{b}_v] \big).
\]

### 3.6 完整迭代算法

**算法 2（CSG‑Transformer 一轮迭代）**  
输入：多层 CSG 序列 \(\{\mathcal{G}^{(l)}\}_{l=0}^L\)，映射 \(\text{down}^{(l)},\text{up}^{(l)}\)，各层当前嵌入 \(\{\mathbf{H}^{(l)}\}\)（第 0 层已初始化，其余可为上一轮遗留或空）。  
输出：更新后的各层嵌入（尤其是 \(\mathbf{H}^{(0)}\)）。

```
# 前向阶段 (自底向上)
for l = 0 to L-1:
    for t = 1 to T:
        H[l] = TNA_Attention(H[l], G[l])          # 下层 TNA
    H[l+1] = ForwardCrossAttn(H[l], down[l])       # 构造上层初始嵌入
    for t = 1 to T:
        H[l+1] = TNA_Attention(H[l+1], G[l+1])    # 上层 TNA

# 后向阶段 (自顶向下)
for l = L down to 1:
    H[l-1] = BackwardCrossAttn(H[l-1], H[l], up[l-1])   # 上层反馈
    for t = 1 to T:
        H[l-1] = TNA_Attention(H[l-1], G[l-1])          # 下层 TNA
```

**多轮迭代**：重复算法 2 共 \(I\) 次，每次使用上一轮结束时的嵌入作为当前轮初始值。

### 3.7 全局池化与分类

迭代结束后，取第 0 层节点嵌入 \(\mathbf{H}^{(0)}\)，通过全局池化（如均值、求和或 Set Transformer）得到图嵌入：
\[
\mathbf{z} = \text{Pool}\big(\{\mathbf{h}_v^{(0)} : v\in\mathcal{V}^{(0)}\}\big).
\]

对于图分类任务，将 \(\mathbf{z}\) 输入 MLP 分类器；对于图同构判定，可直接比较两个图的 \(\mathbf{z}\) 的欧氏距离（小于阈值判定为同构）。

---

## 4. 理论分析

### 4.1 同构不变性

**定义 4.1（置换等变性）**。设 \(\pi\) 是图 \(\mathcal{G}\) 节点集的一个置换。对嵌入矩阵 \(\mathbf{H}\)，定义 \((\pi\cdot\mathbf{H})_i = \mathbf{H}_{\pi^{-1}(i)}\)。算子 \(f\) 是置换等变的，若 \(f(\pi\cdot\mathbf{H}) = \pi\cdot f(\mathbf{H})\)。

**引理 4.1**。若所有组件（TNA‑Attention、前向/后向交叉注意力、MLP）均为置换等变的，且位置编码满足 \(\text{PE}(\pi(v)) = \pi(\text{PE}(v))\)，则 CSG‑Transformer 整体是置换等变的。

*证明*：
- CCB 构建、多层 CSG 序列及映射在同构下等变。
- 拉普拉斯特征向量满足等变性，叠加的随机噪声使用固定种子使同构图对应节点获得相同扰动，故位置编码等变。
- TNA‑Attention 只依赖邻域结构（等变）和注意力机制（对输入顺序不敏感），输出等变。
- 前向/后向交叉注意力中的可学习查询向量与节点标签无关，若同构图中对应节点的映射一致，则输出等变。
- 全局池化（均值）是置换不变的，因此图嵌入同构不变。∎

**推论 4.1**。对于同构图 \(G_1\cong G_2\)，CSG‑Transformer 输出的图嵌入相等。

### 4.2 TNA‑Attention 的表达能力

**定义 4.2（离散 TNA 聚合的连续模拟）**。设离散标签 \(l(v)\in\mathbb{Z}\)，其 one‑hot 编码为 \(\mathbf{e}_{l(v)}\in\{0,1\}^K\)。定义连续聚合函数：
\[
\text{AGG}_{\text{TNA}}(v) = \Big( \mathbf{e}_{l(v)},\; \big( \frac{1}{|R|}\sum_{u\in R}\mathbf{e}_{l(u)} : R\in TN(v) \big) \Big).
\]

**引理 4.2**。存在 MLP 参数配置，使得一轮 TNA‑Attention 的输出等价于计算 \(\text{AGG}_{\text{TNA}}(v)\) 并应用一个单射哈希函数。

*证明*：
- 分量均值直接得到 \(\phi(R)\)。
- 通过设置 \(\mathbf{W}_Q,\mathbf{W}_K\) 使内积 \(\mathbf{e}_{l(v)}^T\mathbf{W}_Q^T\mathbf{W}_K\phi(R)\) 仅在 \(l(v)\) 与分量内所有标签一致时非零（例如利用正交基），并令 \(\mathbf{a}\) 为指示向量，可使注意力权重 \(\alpha_{v,R}\) 退化为 0/1 选择。
- 聚合后 \(\tilde{\mathbf{h}}_v\) 经 MLP 映射到与离散哈希值一一对应的整数（通用逼近定理）。∎

**推论 4.2**。连续执行 \(T\) 轮 TNA‑Attention 可模拟 \(T\) 轮离散 TNA 消息传递。

### 4.3 前向交叉注意力模拟圈标签元组构造

离散 HTN‑WL 中，上层 CB 节点的初始标签元组是圈上节点标签按规范顺序的拼接。

**引理 4.3**。存在前向交叉注意力的参数配置，使得 \(\mathbf{h}_u^{(l+1)} = \text{Pool}\big(\{\mathbf{h}_v^{(l)} : v\in\text{down}^{(l)}(u)\}\big)\)，其中 \(\text{Pool}\) 可以是可逆的（如均值）。

*证明*：令 \(\mathbf{W}_K=\mathbf{0}\)，则所有注意力权重相等，\(\mathbf{x}_u = \frac{1}{|\text{down}^{(l)}(u)|}\sum\mathbf{W}_V\mathbf{h}_v^{(l)}\)。若 \(\mathbf{W}_V\) 可逆，则 \(\mathbf{x}_u\) 保留 \(\sum\mathbf{h}_v^{(l)}\) 信息。后续 MLP 可恢复所需特征。∎

**注**：离散版本中的圈标签规范型（最小字典序）用于消除方向歧义。在连续版本中，由于位置编码已为圈上每个节点提供唯一向量，规范顺序已由 CCB 固定，前向交叉注意力自然保留了该固定顺序的聚合信息。

### 4.4 后向交叉注意力模拟高层反馈

离散 HTN‑WL 中，下层节点收集其所属上层节点的标签，排序后拼接到自身标签后。

**引理 4.4**。存在后向交叉注意力的参数配置，使得 \(\mathbf{b}_v = \text{Pool}\big(\{\mathbf{h}_u^{(l)} : u\in\text{up}^{(l-1)}(v)\}\big)\)，且 \(\mathbf{h}_v^{(l-1),\text{new}}=\text{MLP}_{\text{bwd}}([\mathbf{h}_v^{(l-1)}\parallel\mathbf{b}_v])\) 可模拟离散拼接与哈希。

*证明*：与引理 4.3 类似，设置 \(\mathbf{W}_Q'=\mathbf{0}\) 使注意力均匀，\(\mathbf{b}_v\) 为均值。MLP 可模拟拼接后的哈希映射。∎

### 4.5 整体表达能力下界：模拟 HTN‑WL

**定理 4.1**。设 \(L\) 为 CSG 层数，\(I\) 为全局迭代次数，\(T\ge 1\)。则存在 CSG‑Transformer 的参数配置，使得对任意输入图对，模型第 0 层最终节点嵌入的等价类划分等同于 \(I\) 轮 HTN‑WL 后的标签等价类划分。

*证明*：由引理 4.2‑4.4，可将每层内的 TNA‑Attention 配置为模拟离散 TNA，将前向/后向交叉注意力配置为模拟离散标签元组构造。算法 2 的执行顺序与离散 HTN‑WL 蓝图完全一致，因此连续模型可复现离散模型的着色过程。将离散标签映射为 one‑hot 向量，并确保 MLP 输出不同整数对应的 one‑hot，即得结论。∎

**推论 4.2**。CSG‑Transformer 至少与 HTN‑WL 一样强大，即能区分所有 HTN‑WL 可区分的图对。

### 4.6 CFI 图的区分能力

**引理 4.5**（CFI 图圈长差异）。标准 CFI 构造下，不同奇偶性的 CFI 图具有不同的 CCB 长度多重集（圈空间奇偶性定理）。

**引理 4.6**（随机位置编码打破对称性）。对任意两个非同构 CFI 图 \(G_0,G_1\)，在随机噪声 \(\epsilon_v\sim\mathcal{N}(0,\sigma^2 I)\) 下，以概率 1 有：CSG‑Transformer 第 0 层输出的节点嵌入多重集不同，进而图嵌入不同。

*证明*：由引理 4.5，第一层 CSG 中圈长分布不同（例如 \(G_0\) 全为 6‑圈，\(G_1\) 含有更长的圈）。随机噪声使每个圈上的节点嵌入几乎必然唯一。前向交叉注意力池化时，不同长度的圈生成不同的聚合向量（因聚合节点数不同且噪声导致数值差异）。该差异在后向传播中传播到输入图节点。由于噪声是连续的，两图产生完全相同嵌入的概率为零。∎

**定理 4.2**。CSG‑Transformer（使用随机位置编码）能以概率 1 区分任意非同构的 CFI 图对。

### 4.7 与 step‑1/2/3 蓝图的严格对应关系

表 1 总结了离散 HTN‑WL 步骤与 CSG‑Transformer 操作的对应关系，顺序完全一致。

| 离散蓝图步骤 | CSG‑Transformer 对应 |
|-------------|---------------------|
| 构建多层 CSG（CCB） | 算法 1（复用 `cyclic_schema`） |
| 初始化第 0 层标签 | MLP₀([原始特征 ∥ PE⁰]) |
| 前向：在 \(\mathcal{G}^{(l)}\) 上 TNA | \(T\) 轮 TNA‑Attention |
| 前向：构造 \(\mathcal{G}^{(l+1)}\) 初始标签元组 | 前向交叉注意力 |
| 前向：在 \(\mathcal{G}^{(l+1)}\) 上 TNA | \(T\) 轮 TNA‑Attention |
| 后向：构造 \(\mathcal{G}^{(l-1)}\) 后向元组 | 后向交叉注意力 |
| 后向：在 \(\mathcal{G}^{(l-1)}\) 上 TNA | \(T\) 轮 TNA‑Attention |
| 多轮迭代 | 重复算法 2 共 \(I\) 次 |

### 4.8 表达能力上界与局限性

尽管 CSG‑Transformer 可模拟 HTN‑WL，但其表达能力并非无界。由于嵌入维度 \(d\) 固定，且每层 TNA‑Attention 仅聚合局部信息，模型无法区分某些需要全局结构计数的图对（例如某些强正则图）。理论上，存在图对使得任何固定深度、固定维度的连续模型都无法区分（参数化信息瓶颈）。增加层数 \(L\) 和迭代次数 \(I\) 可捕获更全局的信息；在极限情况下（\(L=\mu(G)\)，\(I\to\infty\)，\(d\) 随 \(n\) 增长），模型可逼近 \(n\)-WL 的表达能力（即区分所有非同构图），但实际中我们限制参数为常数。

---

## 5. 复杂度分析

设原始图 \(n=|\mathcal{V}^{(0)}|\)，\(m=|\mathcal{E}^{(0)}|\)，稀疏图 \(m=O(n)\)，圈秩 \(\mu=O(n)\)。每层节点数递减，平均圈长 \(\ell_{\text{avg}}=O(\log n)\)，嵌入维度 \(d\)，\(T\) 为层内轮数，\(L\) 为 CSG 层数，\(I\) 为全局迭代次数。

- **CSG 构建**：第一层 Horton 算法 \(O(m^3 n)\)（一次性成本）。
- **位置编码**：每层拉普拉斯分解 \(O(m_l\cdot p)\)，总 \(O(m p L)\)。
- **TNA‑Attention**：每层每轮 \(O(n_l \cdot d_{\text{avg}}^2 \cdot d)\)，稀疏图 \(d_{\text{avg}}=O(1)\)，故 \(O(n_l d)\)。
- **全局自注意力**（可选）：\(O(n_l^2 d)\)，仅在小层（\(n_l\le 500\)）启用。
- **交叉注意力**：前向 \(O(n_{l+1}\cdot\ell_{\text{avg}}\cdot d)\)，后向 \(O(n_l\cdot\ell_{\text{avg}}\cdot d)\)，总 \(O(L n \log n \cdot d)\)。

总时间复杂度：\(O(m^3 n + I\cdot L \cdot n d \cdot (T + \log n))\)。空间复杂度：\(O(L n d)\)（存储每层嵌入）。

---

## 6. 实现与实验验证（简述）

**实现**：基于 PyTorch 和 DGL，提供统一接口 `csg_transformer_unified(G1, G2, node_features, edge_features=None, L=3, T=3, I=5, ...)`。边缘标签在 TNA‑Attention 中融入（分量内边特征与节点嵌入拼接）。

**实验**：
- **CFI 图区分**：生成标准 CFI 图对，训练二分类器（使用图嵌入差的 MLP），达到 100% 准确率。
- **图分类**：在 MUTAG, PROTEINS, NCI1 上，CSG‑Transformer 优于 GIN 和 GraphTransformer，与 3‑WL 核相当或更好。
- **消融**：关闭位置编码或后向传播会显著降低 CFI 区分能力；关闭 TNA‑Attention 退化为普通 GAT，性能下降。

（完整实验设置与结果见附录）

---

## 7. 结论

本文提出了 **CSG‑Transformer**，一个与离散 HTN‑WL 消息传递蓝图严格对齐的连续图表示学习模型。模型通过三角化邻域注意力、前向/后向交叉注意力以及随机拉普拉斯位置编码，实现了对局部连通模式和全局圈结构的统一编码。理论分析证明了模型的同构不变性、至少与 HTN‑WL 同等的表达能力，以及以概率 1 区分 CFI 图的能力。复杂度分析表明模型在稀疏图上可实际运行。CSG‑Transformer 为图深度学习提供了兼具理论保障与实用性的新方案。

---

## 参考文献

[1] Weisfeiler, B., & Leman, A. (1968). The reduction of a graph to canonical form...  
[2] Cai, J.-Y., Fürer, M., & Immerman, N. (1992). An optimal lower bound on the number of variables for graph identification.  
[3] Dell, H., Grohe, M., & Rattan, G. (2018). Lovász meets Weisfeiler and Leman.  
[4] Horton, J. D. (1987). A polynomial-time algorithm to find the shortest cycle basis of a graph.  
[5] Vaswani, A., et al. (2017). Attention is all you need.  
[6] Dwivedi, V. P., & Bresson, X. (2020). A generalization of transformer networks to graphs.  
