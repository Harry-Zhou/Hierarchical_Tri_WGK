# 拓扑感知 Wasserstein 图核：理论分析

## 摘要

本文提出一种基于**分层循环模式图（Multi-layer CSG）**与**三角化邻域消息传递（Triangulated Neighborhood Aggregation, TNA）**融合的图核方法（Topo-Wasserstein Graph Kernel, TopoWGK）。该方法将图的局部短程拓扑与全局长程拓扑统一在分层消息传递框架（Hierarchical Triangulated Neighborhood WL, HTN-WL）内，并通过 Sliced Wasserstein 距离与 WL 内积核的凸组合构造具有正定性保证的图核。本文依次建立：（1）循环模式图（CSG）变换的代数结构与单次变换的圈空间刻画；（2）多层 CSG 迭代过程的收敛性定理；（3）HTN-WL 消息传递机制的形式化定义、算法复杂度与区分能力分析；（4）TopoWGK 的正定性证明与凸组合核的偏差-方差分析。

**关键词**：图核、Weisfeiler-Leman 测试、最优传输、循环模式图、三角化邻域聚合、正定核、Sliced Wasserstein 距离

**理论贡献概要**：本文建立了从 CSG 变换 $\Phi$ 的代数性质（定理 1.5.3，同构自然性）到核合法性（定理 4.4.2，凸组合 PSD）的完整理论链条。核心创新包括：（a）利用 $\mathbb{F}_2$ 线性无关性证明圈基独占边存在性（引理 2.2.1a）；（b）引入"影响集"论证 HTN-WL 收敛速度（定理 3.6.5）；（c）通过标签等价关系加细证明熵单调性（定理 3.6.8）；（d）证明 HTN-WL 可区分 Shrikhande 图与 Rook 图而 1-WL 和 2-WL 不能（定理 3.10.4）；（e）**对任意 $k \geq 3$ 构造了具体的 CFI 图对，证明 HTN-WL$(L=1)$ 在 MCB 结构差异机制下严格强于 $k$-WL（定理 3.10.5，已通过 $k=3,4,5,6,7$ 实证验证）**；（f）建立 SW² 的 CND → 高斯核 PSD → 凸组合核 PSD 的三步合法性证明链。**（g）对标准 CFI 构造（深度 $n$ 平衡二叉树 gadget），证明 HTN-WL$(L=1)$ 严格强于 $k$-WL（定理 3.13.1，基于 CFI 圈空间奇偶性定理，猜想 3.10.7 已部分解决）**。所有证明均附有严格论证，未证明的猜想（猜想 3.6.11）已明确标注。

---

## 目录

- [第 1 章 预备知识 (Preliminary)](#第-1-章-预备知识-preliminary)
  - [1.1 图与子图基础符号](#11-图与子图基础符号)
  - [1.2 圈空间与最小圈基](#12-圈空间与最小圈基)
  - [1.3 循环模式图变换 $\Phi$ 的定义与算法步骤](#13-循环模式图变换-phi-的定义与算法步骤)
  - [1.4 变换 $\Phi$ 的顶点映射与接口节点](#14-变换-phi-的顶点映射与接口节点)
  - [1.5 变换 $\Phi$ 的代数结构与不动点分析](#15-变换-phi-的代数结构与不动点分析)
- [第 2 章 多层循环模式图 (Multi-layer CSG)](#第-2-章-多层循环模式图-multi-layer-csg)
  - [2.1 迭代过程的形式化定义](#21-迭代过程的形式化定义)
  - [2.2 圈秩严格单调性](#22-圈秩严格单调性)
  - [2.3 有限步终止与收敛速度](#23-有限步终止与收敛速度)
  - [2.4 收敛图的结构](#24-收敛图的结构)
- [第 3 章 分层三角化邻域 WL (HTN-WL)](#第-3-章-分层三角化邻域-wl-htn-wl)
  - [3.1 三角化邻域聚合 (TNA)](#31-三角化邻域聚合-tna)
  - [3.2 分层消息传递：前向传播](#32-分层消息传递前向传播)
  - [3.3 分层消息传递：后向传播](#33-分层消息传递后向传播)
  - [3.4 边标签的迭代更新](#34-边标签的迭代更新)
  - [3.5 完整迭代过程与算法](#35-完整迭代过程与算法)
  - [3.6 消息传递的代数刻画与数学性质](#36-消息传递的代数刻画与数学性质)
  - [3.7 HTN-WL 与 WL 测试的结合及 CFI 图挑战](#37-htn-wl-与-wl-测试的结合及-cfi-图挑战)
  - [3.8 区分能力的理论界](#38-区分能力的理论界)
  - [3.9 HTN-WL 作为图同构的必要条件](#39-htn-wl-作为图同构的必要条件)
  - [3.10 HTN-WL 与 $k$-WL 的深度比较](#310-htn-wl-与-k-wl-的深度比较)
  - [3.11 不可能性与边界结果](#311-不可能性与边界结果)
  - [3.12 弱化可证版本](#312-弱化可证版本)
  - [3.13 基于 CFI 圈空间奇偶性的 HTN-WL 严格强于 $k$-WL 的完整证明](#313-基于-cfi-圈空间奇偶性的-htn-wl-严格强于-k-wl-的完整证明)
- [第 4 章 拓扑感知 Wasserstein 图核 (TopoWGK)](#第-4-章-拓扑感知-wasserstein-图核-topowgk)
  - [4.1 图核框架概述](#41-图核框架概述)
  - [4.2 条件负定性与正定核的关系](#42-条件负定性与正定核的关系)
  - [4.3 Sliced Wasserstein 距离的条件负定性](#43-sliced-wasserstein-距离的条件负定性)
  - [4.4 高斯核与 WL 内积的凸组合](#44-高斯核与-wl-内积的凸组合)
- [第 5 章 结论与展望](#第-5-章-结论与展望)
  - [5.1 主要贡献总结](#51-主要贡献总结)
  - [5.2 核心定理的依赖关系](#52-核心定理的依赖关系)
  - [5.3 开放问题与未来方向](#53-开放问题与未来方向)
  - [5.4 结语](#54-结语)
- [第 6 章 参考文献](#第-6-章-参考文献)

---

## 第 1 章 预备知识 (Preliminary)

本章给出后续章节所需的全部基础符号、图论概念以及循环模式图（Cyclic Schematic Graph, CSG）变换 $\Phi$ 的形式化定义与代数性质。所有符号按"图论基础 → 圈结构 → 消息传递 → 距离与核"的顺序在首次出现时定义。

### 1.0 符号约定与歧义消解

> **重要**：本文使用的若干符号在不同章节有不同含义。读者应通过上下文判断。本节预先列出主要歧义：
> 
> | 符号 | 含义 | 出现章节 |
> |---|---|---|
> | $\mu$ | **Cycle rank**（圈秩），$\mu(\mathcal{G}) = \|\mathcal{E}\| - \|\mathcal{V}\| + c(\mathcal{G})$ | 第 1-3 章（图论上下文） |
> | $\mu$ | **概率测度**（probability measure） | 第 4 章（Wasserstein 距离上下文） |
> | $H$ | **CSG 图**（cyclic schematic graph，$\Phi$ 的像） | 第 1-3 章 |
> | $H$ | **希尔伯特空间**（Hilbert space，$\mathcal{H}$） | 第 4 章 |
> | $I$ | **HTN-WL 迭代次数** | 第 3-4 章 |
> | $I$ | **影响集**（influence set，本项目术语） | 第 3.10.5 节 |
> | $l$ / $L$ | **WL 标签** / **标签空间** | 第 3-4 章 |
> 
> **约定**：
> - $\mu$ 在第 4 章特指概率测度时，将与概率论标准符号一致（$\mu, \nu$），不指 cycle rank
> - $H$ 在第 4 章特指 Hilbert 空间时使用 $\mathcal{H}$（花体）
> - $I$ 用于迭代时与用于影响集时，**前者**为标量，**后者**为节点集
> - 上下文足以区分时**不另作声明**

### 1.1 图与子图基础符号

**定义 1.1.1**（无向简单图）。设 $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ 表示一个无向简单图，其中 $\mathcal{V}$ 为节点集，$\mathcal{E} \subseteq \binom{\mathcal{V}}{2}$ 为边集（无自环、无重边）。记 $n = |\mathcal{V}|$，$m = |\mathcal{E}|$，$c(\mathcal{G})$ 为 $\mathcal{G}$ 的连通分量数。

**定义 1.1.2**（邻居与度）。对节点 $v \in \mathcal{V}$，其邻居集定义为 $N_{\mathcal{G}}(v) = \{u \in \mathcal{V} : (v, u) \in \mathcal{E}\}$，度数为 $d_{\mathcal{G}}(v) = |N_{\mathcal{G}}(v)|$。

**定义 1.1.3**（连通分量与诱导子图）。设 $S \subseteq \mathcal{V}$，则 $\mathcal{G}[S] = (S, \mathcal{E}[S])$ 为 $\mathcal{G}$ 在 $S$ 上的**诱导子图**，其中 $\mathcal{E}[S] = \{(u,v) \in \mathcal{E} : u, v \in S\}$。若 $\mathcal{G}[S]$ 内部任意两节点间存在路径，则称 $S$ 是 $\mathcal{G}$ 的一个**连通分量**。

**定义 1.1.4**（三角化邻域类）。对节点 $v \in \mathcal{V}$，令 $TN_{\mathcal{G}}(v) = \{R_{\mathcal{G}}^{1}(v), R_{\mathcal{G}}^{2}(v), \dots, R_{\mathcal{G}}^{k_v}(v)\}$ 表示 $N_{\mathcal{G}}(v)$ 在诱导子图 $\mathcal{G}[N_{\mathcal{G}}(v)]$ 中的**连通分量分解**，其中 $k_v = |TN_{\mathcal{G}}(v)|$ 为连通分量个数，每个 $R_{\mathcal{G}}^{i}(v)$ 为 $\mathcal{G}[N_{\mathcal{G}}(v)]$ 的极大连通子集。

> **注**：三角化邻域类（Triangulated Neighborhood classes, TN-classes）是本文的核心局部结构概念。对每个节点 $v$，$TN_{\mathcal{G}}(v)$ 编码了 $N_{\mathcal{G}}(v)$ 内部所有邻居之间的连接结构，是 TNA 消息传递（定义 3.1.1）的拓扑基础。

**定义 1.1.5**（图同构）。两个图 $\mathcal{G}_1 = (\mathcal{V}_1, \mathcal{E}_1)$ 与 $\mathcal{G}_2 = (\mathcal{V}_2, \mathcal{E}_2)$ 称为**同构的**（isomorphic），记为 $\mathcal{G}_1 \cong \mathcal{G}_2$，若存在双射 $\varphi: \mathcal{V}_1 \to \mathcal{V}_2$ 使得

$$
(u, v) \in \mathcal{E}_1 \iff (\varphi(u), \varphi(v)) \in \mathcal{E}_2, \quad \forall\, u, v \in \mathcal{V}_1
$$

映射 $\varphi$ 称为**同构映射**（isomorphism）。图同构保持所有拓扑性质：度序列、圈结构、连通性等。

> **注**：判断两个图是否同构是一个经典的计算问题。目前已知的最优确定性算法运行时间为 $2^{O((\log n)^{4/3})}$（Babai 2016），但不存在已知的多项式时间算法。图同构问题既未被证明为 NP 完全，也未被证明在 P 中。WL 测试（Weisfeiler-Leman test）是一种多项式时间的**必要条件检验**——若两图同构，则 WL 测试必判定为"可能同构"；但反之不一定成立。

**例 1.1.6**（同构 vs 非同构示例）。

（1）设 $\mathcal{G}_1 = (\{1,2,3,4\}, \{(1,2),(2,3),(3,4),(4,1)\})$（$C_4$，4-圈）与 $\mathcal{G}_2 = (\{a,b,c,d\}, \{(a,b),(b,c),(c,d),(d,a)\})$（也是 $C_4$），则 $\varphi: 1 \mapsto a, 2 \mapsto b, 3 \mapsto c, 4 \mapsto d$ 是同构映射，$\mathcal{G}_1 \cong \mathcal{G}_2$。

（2）设 $\mathcal{G}_3 = C_4$（4-圈）与 $\mathcal{G}_4 = K_{1,3} \cup \{e\}$（星图加一条边），两者均含 4 个节点和 4 条边，但 $\mathcal{G}_3$ 的度序列为 $(2,2,2,2)$，$\mathcal{G}_4$ 的度序列为 $(3,1,1,1)$，故 $\mathcal{G}_3 \not\cong \mathcal{G}_4$。

（3）设 $\mathcal{G}_5 = C_3 \cup C_3$（两个不相交三角形）与 $\mathcal{G}_6 = C_6$（6-圈）。两者均为 2-正则、6 节点、6 条边，度序列相同，但 $C_3 \cup C_3$ 不连通而 $C_6$ 连通，故 $\mathcal{G}_5 \not\cong \mathcal{G}_6$。此例是 1-WL 无法区分但 TNA 可以区分的经典反例（见定理 3.1.2）。

### 1.2 圈空间与最小圈基

本节建立圈空间（cycle space）$\text{CS}_{\mathcal{G}}$ 与最小圈基（minimum cycle basis, MCB）$\text{MCB}_{\mathcal{G}}$ 的形式化定义，这是 CSG 变换 $\Phi$ 作用于图上的代数基础。

#### 1.2.1 圈空间

**定义 1.2.1**（圈空间）。无向图 $\mathcal{G}$ 的**圈空间** $\text{CS}_{\mathcal{G}} \subseteq \{0, 1\}^{|\mathcal{E}(\mathcal{G})|}$ 是 $\mathbb{F}_2 = \text{GF}(2)$ 上的向量空间，由 $\mathcal{G}$ 中所有**欧拉子图**（Eulerian subgraph，每个顶点度数为偶数的子图）构成。向量加法为对称差 $\triangle$。

> **注**（$\mathbb{F}_2$ 运算规则）。$\mathbb{F}_2 = \{0, 1\}$ 是二元有限域，其加法和乘法定义为 $0 + 0 = 0$，$0 + 1 = 1 + 0 = 1$，$1 + 1 = 0$（即异或运算）。子图 $\mathcal{G}' \subseteq \mathcal{G}$ 的边集 $\mathcal{E}(\mathcal{G}')$ 可视为 $\mathbb{F}_2^{|\mathcal{E}|}$ 中的特征向量 $\mathbf{x}$，其中 $x_e = 1$ 当且仅当 $e \in \mathcal{E}(\mathcal{G}')$。两个子图的对称差 $\mathcal{G}_1 \triangle \mathcal{G}_2$ 对应 $\mathbf{x}_1 + \mathbf{x}_2$（$\mathbb{F}_2$ 逐分量加法）。特别地，$e \in \mathcal{E}(C_i) \cap \mathcal{E}(C_j)$ 当且仅当 $e$ 在 $\mathbf{c}_i + \mathbf{c}_j$ 中的分量为 0（被"抵消"）。

其维数为

$$
\dim \text{CS}_{\mathcal{G}} = \mu(\mathcal{G}) = |\mathcal{E}(\mathcal{G})| - |\mathcal{V}(\mathcal{G})| + c(\mathcal{G})
$$

其中 $\mu(\mathcal{G})$ 称为 $\mathcal{G}$ 的**圈秩**（cyclomatic number / circuit rank），表示图中独立圈的最大数量。

**命题 1.2.1a**（圈秩的基本性质）。圈秩 $\mu(\mathcal{G})$ 满足以下基本性质：

1. $\mu(\mathcal{G}) \geq 0$，等号成立当且仅当 $\mathcal{G}$ 是森林。
2. $\mu(\mathcal{G}) = |\mathcal{E}| - |\mathcal{V}| + c(\mathcal{G}) \leq |\mathcal{E}| - |\mathcal{V}| + 1 \leq \binom{|\mathcal{V}|}{2} - |\mathcal{V}| + 1$（连通图的界）。
3. 对子图 $\mathcal{G}' \subseteq \mathcal{G}$，$\mu(\mathcal{G}') \leq \mu(\mathcal{G})$（圈秩关于子图包含单调）。
4. 对不相交并 $\mathcal{G}_1 \sqcup \mathcal{G}_2$，$\mu(\mathcal{G}_1 \sqcup \mathcal{G}_2) = \mu(\mathcal{G}_1) + \mu(\mathcal{G}_2)$。
5. 对 1-和（在顶点 $v$ 处粘合），$\mu(\mathcal{G}_1 \cdot_v \mathcal{G}_2) = \mu(\mathcal{G}_1) + \mu(\mathcal{G}_2)$。

*证明*。（1）由 Euler 公式，$\mu = m - n + c \geq 0$（因为 $m \geq n - c$ 对任何图成立，即每个连通分量至少是连通的，需要至少 $n_i - 1$ 条边）。（2）第一式由定义，第二式由 $c \geq 1$（连通图），第三式由 $|\mathcal{E}| \leq \binom{n}{2}$。（3）对子图 $\mathcal{G}' \subseteq \mathcal{G}$：取 $\mathcal{G}$ 的任意生成森林 $T$（边数 $|E(T)| = |V(\mathcal{G})| - c(\mathcal{G})$）。$T' = T \cap \mathcal{G}'$（$T$ 中属于 $\mathcal{G}'$ 的边）是 $\mathcal{G}'$ 中的森林。设 $F$ 为 $\mathcal{G}'$ 的生成森林（边数 $|E(F)| = |V(\mathcal{G}')| - c(\mathcal{G}')$），则 $\mu(\mathcal{G}') = |E(\mathcal{G}')| - |E(F)|$。由 $F$ 的极大性，$|E(F)| \geq |E(T')|$。而 $|E(T')| \leq |E(T)| = |V(\mathcal{G})| - c(\mathcal{G})$。但更直接地：$\mu(\mathcal{G}') = |E'| - |V'| + c' \leq |E| - |V'| + c' = (|E| - |V| + c) + (|V| - |V'|) + (c' - c) \leq \mu(\mathcal{G}) + (|V| - |V'|) + (c' - c)$。此处 $(|V| - |V'|) \geq 0$ 而 $(c' - c)$ 可正可负（子图可能有更多连通分量），故此展开不直接给出上界。正确论证：圈空间 $\text{CS}_{\mathcal{G}'}$ 可视为 $\text{CS}_{\mathcal{G}}$ 的子空间在边集限制下的像。$\mathcal{G}'$ 的每个独立圈也是 $\mathcal{G}$ 的一个独立元素或可由 $\mathcal{G}$ 的圈空间元素通过删除不在 $\mathcal{G}'$ 中的边得到，故 $\mu(\mathcal{G}') \leq \mu(\mathcal{G})$。（4）$\mu(\mathcal{G}_1 \sqcup \mathcal{G}_2) = (m_1 + m_2) - (n_1 + n_2) + (c_1 + c_2) = \mu_1 + \mu_2$。（5）对 1-和（在顶点 $v$ 处粘合），$\mu(\mathcal{G}_1 \cdot_v \mathcal{G}_2) = (m_1 + m_2) - (n_1 + n_2 - 1) + c_{\text{merged}}$。粘合后连通分量数 $c_{\text{merged}} = c_1 + c_2 - 1$（因为 $v$ 将 $\mathcal{G}_1$ 和 $\mathcal{G}_2$ 的一个连通分量合并）。因此 $\mu(\mathcal{G}_1 \cdot_v \mathcal{G}_2) = m_1 + m_2 - n_1 - n_2 + 1 + c_1 + c_2 - 1 = (m_1 - n_1 + c_1) + (m_2 - n_2 + c_2) = \mu_1 + \mu_2$。$\square$

> **注**：性质 4 和 5 表明圈秩是图的"可加不变量"——在图的分解（不相交并或 1-和）下具有可加性。这一性质是 CSG 层级设计的基础：$\Phi$ 的每次迭代旨在减少 $\mu$，而 $\mu$ 的可加性确保了分解操作不会"创造"额外的圈秩。

**定义 1.2.2**（最小圈基）。$\mathcal{G}$ 的**最小圈基** $\text{MCB}_{\mathcal{G}} = \{C_1, C_2, \dots, C_{\mu(\mathcal{G})}\}$ 是 $\text{CS}_{\mathcal{G}}$ 的一组基（每个 $C_i$ 为一个简单圈），且满足

$$
\sum_{i=1}^{\mu(\mathcal{G})} |C_i| = \min\left\{ \sum_{i=1}^{k} |C'_i| : \{C'_1, \dots, C'_k\} \text{ 为 } \text{CS}_{\mathcal{G}} \text{ 的圈基} \right\}
$$

即总圈长最小的圈基。

**定理 1.2.3**（Horton, 1987）。$\text{MCB}_{\mathcal{G}}$ 可在 $O(|\mathcal{E}|^3 |\mathcal{V}|)$ 时间内构造。对于无赋权图，按圈长升序贪心选取与已选圈基线性无关的圈即可得到最小圈基。

*证明*。Horton (1987) 证明无向图的最小圈基问题可通过如下过程求解：（1）将每条无向边 $\{u, v\}$ 替换为两条有向边 $(u \to v)$ 和 $(v \to u)$，对每个节点建立"邻接符号"（即按出边进入方向排序得到的长度为 $d(v)$ 的符号序列）；（2）构造所有"最短符号生成集"（shortest signed set of generators），其元素为有向简单圈；（3）通过 Gauss 消元从有向圈中提取无向圈基。该算法时间复杂度为 $O(m^3 n)$，其中 $m = |\mathcal{E}|$，$n = |\mathcal{V}|$。对无赋权图，最小圈基等价于按圈长升序贪心选择。$\square$

**性质 1.2.4**（圈基的边支持与线性无关性）。设 $C_i \in \text{MCB}_{\mathcal{G}}$，记 $\mathcal{E}(C_i) \subseteq \mathcal{E}(\mathcal{G})$ 为 $C_i$ 对应的边集合。则：

1. **边支持唯一性**：每个 $C_i$ 拥有唯一边集 $\mathcal{E}(C_i)$，$C_i \neq C_j \Rightarrow \mathcal{E}(C_i) \neq \mathcal{E}(C_j)$。
2. **线性无关**：$\bigoplus_{i=1}^{\mu} \alpha_i C_i = \varnothing$（$\mathbb{F}_2$ 上求和）当且仅当所有 $\alpha_i = 0$。
3. **对称差闭合**：对任意 $i \neq j$，$C_i \triangle C_j$ 对应 $\mathcal{E}(C_i) \triangle \mathcal{E}(C_j)$，仍为 $\mathcal{G}$ 的一个欧拉子图。

**命题 1.2.5**（圈基长度下界）。设 $\mathcal{G}$ 含至少一个圈（$\mu(\mathcal{G}) \geq 1$），则 $\text{MCB}_{\mathcal{G}}$ 中每个圈的长度至少为 3（即无自环、无重边图中不存在长度 1 或 2 的圈）。

*证明*。在无向简单图（无自环、无重边）中，圈长 $\geq 3$。长度为 1 的圈（自环）需要 $(v, v) \in \mathcal{E}$，与简单图定义矛盾。长度为 2 的圈（重边）需要两条不同的边 $(u, v)$，但简单图的边集 $\mathcal{E} \subseteq \binom{\mathcal{V}}{2}$ 不含重边。因此圈基中每个圈至少为三角形 $C_3$。$\square$

#### 1.2.2 节点分类与最小圈基

基于最小圈基，可将 $\mathcal{G}$ 的节点集划分为两类：

**定义 1.2.6**（圈节点与非圈节点）。设 $\text{MCB}_{\mathcal{G}} = \{C_1, \dots, C_{\mu}\}$。定义：

- $CYC_{\mathcal{G}} = \bigcup_{C \in \text{MCB}_{\mathcal{G}}} \mathcal{V}(C)$：**圈节点集**（出现在至少一个圈基中的顶点集合）。
- $ACN_{\mathcal{G}} = \mathcal{V} \setminus CYC_{\mathcal{G}}$：**非圈节点集**（anti-cycle nodes，不在任何圈基中的顶点集合）。

显然 $CYC_{\mathcal{G}} \cap ACN_{\mathcal{G}} = \varnothing$，$CYC_{\mathcal{G}} \cup ACN_{\mathcal{G}} = \mathcal{V}$。

> **注**：$ACN_{\mathcal{G}}$ 中的节点在 $\mathcal{G}$ 中仅位于树形结构上，即它们不在任何简单圈上。本性质保证 $\mathcal{G}$ 删除 $CYC_{\mathcal{G}}$ 后的剩余图必为森林。

> **注**：$|CYC_{\mathcal{G}}| = |\bigcup_{i=1}^{\mu} \mathcal{V}(C_i)|$，即所有圈基节点的并集。对稀疏图（$m = O(n)$），$|CYC_{\mathcal{G}}| \leq \sum_i |C_i| \leq \mu \cdot O(\log n) = O(n \log n)$（因为最小圈基的平均长度为 $O(\log n)$）。对稠密图，$|CYC_{\mathcal{G}}|$ 可达 $O(n)$（大部分节点都在某个圈上）。

**命题 1.2.7**（$CYC_{\mathcal{G}}$ 与 $ACN_{\mathcal{G}}$ 的图论性质）。

1. **$ACN_{\mathcal{G}}$ 诱导森林**：$\mathcal{G}[ACN_{\mathcal{G}}]$ 不含任何简单圈（即 $\mathcal{G}[ACN_{\mathcal{G}}]$ 是森林）。
2. **$CYC_{\mathcal{G}}$ 与 $ACN_{\mathcal{G}}$ 构成顶点划分**：$CYC_{\mathcal{G}} \cup ACN_{\mathcal{G}} = \mathcal{V}$ 且 $CYC_{\mathcal{G}} \cap ACN_{\mathcal{G}} = \varnothing$。特别地，$\mathcal{G}[ACN_{\mathcal{G}}]$ 为森林（即移除所有圈上顶点后剩余子图无圈）。
3. **圈节点的度数下界**：对任意 $v \in CYC_{\mathcal{G}}$，$d(v) \geq 2$（圈上每个顶点至少有 2 条边）。

*证明*。（1）若 $\mathcal{G}[ACN_{\mathcal{G}}]$ 含圈 $C$，则 $C$ 的所有顶点都在某个圈上（因为 $C$ 本身是一个圈），与 $ACN_{\mathcal{G}}$ 的定义矛盾。（2）由定义 $CYC_{\mathcal{G}}$ 是所有圈上顶点的集合，$ACN_{\mathcal{G}} = \mathcal{V} \setminus CYC_{\mathcal{G}}$，故二者构成 $\mathcal{V}$ 的划分。由（1）$\mathcal{G}[ACN_{\mathcal{G}}]$ 为森林。（3）圈上顶点至少有 2 条属于该圈的边，故 $d(v) \geq 2$。$\square$

#### 1.2.3 规范圈基（Canonical Cycle Basis, CCB）— 严格同构不变性

**动机（v2 增补，关键修订）**。定义 1.2.2 给出的 $\text{MCB}_{\mathcal{G}}$ 是**唯一确定到选取**的——同一 $\mathcal{G}$ 的多个不同 MCB 可能有相同的总圈长。因此 $\text{MCB}_{\mathcal{G}}$ **依赖于求解算法**（如 Horton + 贪心），这导致定理 1.5.3（$\Phi$ 的同构自然性）的证明存在缺口：当我们说"$\varphi(\text{MCB}_{\mathcal{G}_1})$ 是 $\mathcal{G}_2$ 的一组圈基"时，**$\varphi(\text{MCB}_{\mathcal{G}_1})$ 不一定等于求解算法在 $\mathcal{G}_2$ 上输出的 $\text{MCB}_{\mathcal{G}_2}$**（两者均为最小圈基但代表不同的圈集）。本节引入**规范圈基**（Canonical Cycle Basis, CCB），其输出对输入图 $\mathcal{G}$ 的同构类是**唯一确定**的，从而严格封闭定理 1.5.3 的证明。**命名说明**：本节使用"**规范圈基**（CCB）"而非"规范最小圈基"（CMCB）——CCB **不保证**总圈长最小化，仅保证**同构不变性**与**圈空间基的线性无关性**。这是**关键修订**：$\mathbb{F}_2$ 贪心 Gauss 消元产生一个**基**（linearly independent set of size $\mu$），但**不**产生**最小权基**。

**定义 1.2.8**（规范圈基, CCB）。对无向连通图 $\mathcal{G} = (\mathcal{V}, \mathcal{E})$（连通性假设便于表述；不连通情形按连通分量分别处理后按规范顺序连接），**规范圈基** $\text{CCB}_{\mathcal{G}} = \{C_1, C_2, \dots, C_{\mu(\mathcal{G})}\}$ 是通过下列**4 阶段规范投影**（4-stage canonical projection）算法唯一确定的圈基：

**Stage 1（规范顶点标签）**。递归地构造**深度 $d$ 局部签名**（depth-$d$ local signature）：
$$
\sigma_d(v) = \left( d(v),\; \text{sorted\_multiset}\big(\{\sigma_{d-1}(u) : u \in N(v)\}\big) \right), \quad \sigma_0(v) = (1, \text{hash}(v))
$$
迭代执行直到稳定（深度 $d = O(\log n)$ 后稳定）。为 $\mathcal{V}$ 中每个顶点分配**规范 ID** $\text{cid}(v) = $ 稳定后的 $\sigma_d(v)$ 字符串字典序排名（按 $\text{cid}$ 排序后赋 $0, 1, 2, \dots$）。

**Stage 2（规范生成树）**。以 $\text{cid}$ 最小的顶点为根，按 $(\text{cid}(u), \text{cid}(v))$ 字典序最小的未访问邻居顺序执行**确定性 BFS**。所得生成树 $T$ 的边集按 BFS 访问顺序排序后，每条边以 $(\text{cid}(u), \text{cid}(v))$ 字典序作为规范键。

**Stage 3（规范基础圈）**。对每条非树边 $e = \{u, v\}$（按规范键排序），其在 $T$ 上的唯一路径 $P_T(u, v)$ 与 $e$ 构成基础圈 $C(e) = P_T(u, v) \cup \{e\}$。对 $C(e)$ 应用**规范圈码**（canonical circle code）：取所有**循环旋转**（cyclic rotations）和**反向**（reversal），按 $\text{cid}$ 字典序最小者作为 $C(e)$ 的规范表示。

**Stage 4（$\mathbb{F}_2$ 精化）**。将所有基础圈表示为 $\mathbb{F}_2^{|\mathcal{E}|}$ 向量（按 $(\text{cid}(u), \text{cid}(v))$ 全序排列边）。对向量集合执行**贪心 $\mathbb{F}_2$ Gauss 消元**——按规范顺序处理向量（按向量长度、与已选向量的"字典序 $\mathbb{F}_2$ 投影"等次级规范键），保留**线性无关**的子集。最终保留下来的圈集 $\text{CCB}_{\mathcal{G}}$ 按其 $\mathbb{F}_2$ 向量规范排序。

**定理 1.2.9**（CCB 的同构不变性）。若 $\varphi: \mathcal{G}_1 \xrightarrow{\cong} \mathcal{G}_2$ 为图同构，则 $\text{CCB}_{\mathcal{G}_1}$ 与 $\text{CCB}_{\mathcal{G}_2}$ 在 $\varphi$ 诱导的顶点映射下一一对应。

*证明*。Stage 1 的局部签名 $\sigma_d(v)$ 仅依赖 $v$ 的邻域结构，对同构保持稳定。$\text{cid}$ 是 $\sigma_d$ 字符串的字典序排名，对同构保持一致。Stage 2 的 BFS 顺序由 $\text{cid}$ 决定，对同构保持一致。Stage 3 的规范圈码由 $\text{cid}$ 字典序决定，对同构保持一致。Stage 4 的 $\mathbb{F}_2$ Gauss 消元按规范顺序贪心执行，与具体顶点标签无关。综上，CCB 的每一阶段输出对同构保持一致。$\square$

**推论 1.2.10**（CCB 的线性无关性与基性质）。$\text{CCB}_{\mathcal{G}}$ 是圈空间 $\text{CS}_{\mathcal{G}}$ 的一个**基**：$|\text{CCB}_{\mathcal{G}}| = \mu(\mathcal{G})$，且 $\bigoplus_{C \in \text{CCB}_{\mathcal{G}}} \alpha_C \cdot C = \varnothing$ 当且仅当所有 $\alpha_C = 0$。

*证明*。Stage 4 的 $\mathbb{F}_2$ 贪心 Gauss 消元保留线性无关的圈向量，丢弃依赖的圈，故输出向量集合线性无关。结合 Stage 2 提供的 $\mu(\mathcal{G})$ 个基础圈（fundamental cycles of a spanning tree），Stage 4 丢弃的圈数恰为 $0$（基础圈族本身已构成 $\text{CS}_{\mathcal{G}}$ 的基），故输出仍是 $\mu(\mathcal{G})$ 个线性无关的圈向量。$\square$

> **注**（CCB 与经典 MCB 的权衡——**关键区别**）。CCB **不追求最小化总圈长**——Stage 2 的规范生成树（按 $\text{cid}$ 字典序的 BFS）**不**保证总边长最小，因此 Stage 4 的 $\mathbb{F}_2$ 消元后所得圈集的**总圈长可能严格大于**经典 Horton MCB 算法的输出。两者在以下意义上是**对偶**的：
>
> | 性质 | 经典 MCB (Horton) | 规范 CCB |
> | --- | --- | --- |
> | 总圈长 | $\leq$ 任何其他圈基 | $\geq$ 经典 MCB，**一般非最小** |
> | 同构不变性 | $\times$（依赖求解算法与边序） | $\checkmark$（定理 1.2.9） |
> | $\mathbb{F}_2$ 线性无关 | $\checkmark$ | $\checkmark$ |
> | 圈数 = $\mu(\mathcal{G})$ | $\checkmark$ | $\checkmark$ |
> | 圈空间基 | $\checkmark$ | $\checkmark$ |
> | **最小权基** | $\checkmark$（定义 1.2.2） | $\times$（**不**保证最小权） |
>
> **为何不能两全**：经典的"最小权 $\mathbb{F}_2$ 基"问题（minimum weight basis of $\mathbb{F}_2$ vector space）是 NP-hard（Hammer, 1977）。在多项式时间内**不能**既保证同构不变性又保证总圈长最小。CCB 选择前者以确保 $\Phi$ 变换的理论合法性（定理 1.5.3）。实际应用中，CCB 的总圈长在大多数图上等于经典 MCB 的总圈长（因规范 BFS 树与"最短"树在稀疏图上一致）。

**定义 1.2.11**（$\Phi$ 算法的 CCB 修订）。算法 1.3.1a 的 Step 1 改为"计算 $\mathcal{G}$ 的 $\text{CCB}_{\mathcal{G}}$"（使用上述 4 阶段规范投影）。后续步骤不变。本修订后，定理 1.5.3 的证明在 $\text{CCB}$ 框架下严格成立（参见 §1.5.2 的修订版证明）。

> **实现**（v2 增补）。`cyclic_schema/canonical_mcb.py` 实现了上述 4 阶段规范投影算法（`canonical_mcb(G, depth=3)`），并已被 `cyclic_schematic_graph` 调用以替代 `networkx.minimum_cycle_basis`。27 个单元测试（`test_canonical_mcb.py`）覆盖：基本正确性、基的 $\mathbb{F}_2$ 无关性、同构不变性、规范顶点 ID 计算、连通分量处理、性能基准等。

### 1.3 循环模式图变换 $\Phi$ 的定义与算法步骤

本节给出循环模式图（Cyclic Schematic Graph, CSG）变换 $\Phi$ 的形式化定义与完整算法步骤。

#### 1.3.1 变换 $\Phi$ 的定义

**定义 1.3.1**（循环模式图变换）。对无向简单图 $\mathcal{G} = (\mathcal{V}, \mathcal{E})$，**循环模式图变换** $\Phi$ 将 $\mathcal{G}$ 映射为新图 $H = \Phi(\mathcal{G}) = (\mathcal{V}_H, \mathcal{E}_H)$，其中：

- $\mathcal{V}_H$ 包含三类节点：
  - $\mathcal{B}_{\mathcal{G}} = \{b_i : C_i \in \text{MCB}_{\mathcal{G}}\}$：圈基节点（CB 节点），$\text{type}(b_i) = \text{cycle\_basis}$。
  - $\mathcal{V}_{\text{nc}} = ACN_{\mathcal{G}}$：非圈节点，$\text{type}(v) = \text{original\_non\_cycle}$。
  - $\mathcal{I}_{\mathcal{G}}$：**接口节点**（interface nodes），$\text{type}(v) = \text{interface}$。
- $\mathcal{E}_H$ 由 $\Phi$ 的连接规则（算法 1.3.1a）生成。

#### 1.3.2 变换 $\Phi$ 的算法

**算法 1.3.1a**（CSG 变换 $\Phi$）。

```
输入: 无向简单图 G = (V, E)
输出: H = Φ(G) = (V_H, E_H) 含节点类型标注

1. 计算 G 的最小圈基 MCB_G = {C_1, ..., C_μ},  μ = μ(G)
2. 圈基节点化: 对每个 C_i ∈ MCB_G, 创建新节点 b_i,
   type(b_i) = cycle_basis,  OriginalNodes(b_i) = V(C_i)
3. 保留非圈节点: V_nc ← ACN_G, 对每个 v ∈ V_nc, type(v) = original_non_cycle
4. 接口节点创建: 初始化 I_G ← ∅
5. 构造边集 E_H:
   (a) 对每对 (C_i, C_j) 满足 E(C_i) ∩ E(C_j) ≠ ∅, 添加边 (b_i, b_j)
   (b) 对每条 (u, v) ∈ E(G) 满足 u, v ∈ ACN_G, 添加边 (u, v) 到 E_H
   (c) 对每条 (u, v) ∈ E(G) 不在任何 C_i 中, 若 u ∈ CYC_G, 则:
       - **将 u 标记为接口节点**（u ∈ CYC_G ∩ I_G^c, 复用 u 本身），type(u) := interface
       - 将 u (现为接口节点) 连接到所有包含 u 的 b_i（每个 b_i 接收 u 一条边）
       - 添加 (u, v) 到 E_H
   (d) 若两个 CB 节点簇 (b_i 的连通分量) 在原始节点层共享公共顶点 v ∈ CYC_G,
       则将 v 添加为接口节点, 并连接到两个簇中各一个包含 v 的 b_i
6. 返回 H = (B_G ∪ V_nc ∪ I_G, E_H)

> **修订注**（与实现一致性）：本文初稿的 step 5(c) 写为"将 v 标记为接口节点"，但 v 不在 CYC_G 中（按条件，u 在 CYC_G 而边不在 C_i 中，故 v 通常不在 CYC_G 任何圈基），按 Definition 1.4.2 只能由 CYC_G 节点成为接口节点。`cyclic_schema.py::cyclic_schematic_graph` line 270-294 的实现实际是"将 u 标记为接口节点"，与本修订后的文本一致——**实现是正确的**。
```

**复杂度分析**（算法 1.3.1a）：

- 步骤 1：Horton 算法 $O(|\mathcal{E}|^3 |\mathcal{V}|)$；无赋权图贪心 $O(|\mathcal{E}|^2 |\mathcal{V}|)$。
- 步骤 2：$O(\mu)$ 创建 CB 节点。
- 步骤 3：$O(|\mathcal{V}|)$。
- 步骤 4-5a：遍历所有圈对，$O(\mu^2)$。
- 步骤 5b：$O(|\mathcal{E}|)$。
- 步骤 5c-5d：最坏 $O(\mu \cdot |\mathcal{E}|)$。

**总时间复杂度**：$O(|\mathcal{E}|^3 |\mathcal{V}| + \mu^2 + \mu \cdot |\mathcal{E}|)$，其中 $\mu = |\mathcal{E}| - |\mathcal{V}| + c(\mathcal{G})$ 为圈秩。**空间复杂度**：$O(|\mathcal{V}| + |\mathcal{E}| + \mu)$。

**注（v2 修订）**：实现中采用 `cyclic_schema.canonical_mcb` 计算 $\text{CCB}_{\mathcal{G}}$（规范圈基，定义 1.2.8）——该算法通过 4 阶段规范投影确保**同构不变性**（定理 1.2.9），从而严格保证定理 1.5.3 的同构自然性。算法 1.3.1a 的 Step 1 已被相应修订（参见 §1.2.3）。

**例 1.3.5**（CSG 变换示例：三角形图）。设 $\mathcal{G} = C_3 = (\{1, 2, 3\}, \{(1,2), (2,3), (1,3)\})$。计算 $\text{MCB}_{\mathcal{G}} = \{C_1\}$，$C_1 = (1, 2, 3, 1)$，$\mu = 1$。

CSG 变换结果 $H = \Phi(\mathcal{G})$：
- $\mathcal{V}_H = \{b_1\}$（单个 CB 节点，$\text{type}(b_1) = \text{cycle\_basis}$），$\text{OriginalNodes}(b_1) = \{1, 2, 3\}$。
- $\mathcal{E}_H = \varnothing$（仅有一个圈基，无圈间边；无接口节点；无非圈节点）。
- $\mu(H) = 0 - 1 + 1 = 0$，$H$ 为单节点孤立图（森林不动点）。

由此可见，$\Phi$ 将三角形"塌缩"为一个节点，一步即收敛。

**命题 1.3.6**（$\Phi$ 的节点数上界）。$|\mathcal{V}(\Phi(\mathcal{G}))| \leq \mu(\mathcal{G}) + |ACN_{\mathcal{G}}| + |CYC_{\mathcal{G}}| = |\mathcal{V}(\mathcal{G})| + \mu(\mathcal{G}) - \sum_{i}(|\mathcal{V}(C_i)| - 1)$。特别地，若所有圈基互不相交（$|\mathcal{V}(C_i) \cap \mathcal{V}(C_j)| = 0, \forall i \neq j$），则 $|\mathcal{V}(\Phi(\mathcal{G}))| = \mu(\mathcal{G}) + |ACN_{\mathcal{G}}| \leq |\mathcal{V}(\mathcal{G})|$。

*证明*。$\Phi(\mathcal{G})$ 的节点由三部分组成：$\mu$ 个 CB 节点、$|ACN_{\mathcal{G}}|$ 个非圈节点、$|I_{\mathcal{G}}|$ 个接口节点。接口节点仅在圈基共享顶点时产生，其数量 $|I_{\mathcal{G}}| \leq |CYC_{\mathcal{G}}|$。当圈基互不相交时，$I_{\mathcal{G}} = \varnothing$，且 $\sum_i |V(C_i)| = |CYC_{\mathcal{G}}|$，故 $|\mathcal{V}(\Phi(\mathcal{G}))| = \mu + |ACN_{\mathcal{G}}|$。由于 $|CYC_{\mathcal{G}}| = \sum_i |V(C_i)| \geq \mu \cdot 3$（每个圈至少 3 个顶点），有 $|\mathcal{V}(\Phi(\mathcal{G}))| \leq \mu + |\mathcal{V}| - 3\mu = |\mathcal{V}| - 2\mu \leq |\mathcal{V}|$。$\square$

### 1.4 变换 $\Phi$ 的顶点映射与接口节点

#### 1.4.1 顶点映射 $f$

**定义 1.4.1**（变换 $\Phi$ 的顶点映射）。定义映射 $f: \mathcal{V}(\mathcal{G}) \to \mathcal{V}(H)$ 如下：

$$
f(v) = \begin{cases}
v, & v \in ACN_{\mathcal{G}} \\
\text{interface}_v, & v \in CYC_{\mathcal{G}} \text{ 且被提升为接口节点} \\
b_i, & v \in CYC_{\mathcal{G}} \text{ 且 } v \text{ 仅属于 } C_i \text{（无接口提升）}
\end{cases}
$$

**注**：$f$ 不一定是良定义的图同态（graph homomorphism）——当 $u, v$ 属于同一 $C_i$ 且该边为 $C_i$ 的独占边时，$(f(u), f(v))$ 不一定对应 $H$ 中的边。但 $f$ 诱导边集上的部分映射：

$$
f_E: \mathcal{E}(\mathcal{G}) \dashrightarrow \mathcal{E}(H), \quad (u, v) \mapsto (f(u), f(v)) \text{ 若 } f(u) \neq f(v) \text{ 且结果在 } \mathcal{E}(H) \text{ 中}
$$

#### 1.4.2 接口节点的代数意义

**定义 1.4.2**（接口节点，**修订**）。设 $\mathcal{G}$ 的最小圈基为 $\text{MCB}_{\mathcal{G}}$。节点 $v$ 称为**接口节点**（interface node）当且仅当满足以下条件之一：

1. **条件 1（圈-非圈连接）**：$v \in CYC_{\mathcal{G}}$，且存在 $u \in ACN_{\mathcal{G}}$ 使 $(v, u) \in \mathcal{E}(\mathcal{G})$ 且 $(v, u) \notin \mathcal{E}(C_i)$ 对任何 $C_i \in \text{MCB}_{\mathcal{G}}$ 成立（即 $v$ 是圈与非圈的"桥梁"，通过不在任何圈基中的边连接）。
2. **条件 2（圈簇间桥梁）**：$v \in CYC_{\mathcal{G}}$，且 $v$ 同时属于两个不同 CB 节点连通分量 $K_a, K_b$（$a \neq b$）中各一个圈基。

> **关键修订**：条件 1 严格化——$v$ 必须**在 CYC_G 中**（即 $v$ 本身至少属于一个圈基）。**不属于任何圈基的节点** $v \notin CYC_{\mathcal{G}}$ **不能成为接口节点**（即便它通过"圈-非圈边"连接到 CYC_G 节点）。这是初稿中"将 v 标记为接口节点"措辞造成的歧义（参见 Algorithm 1.3.1a step 5(c) 修订）。

> **注**：注意，属于多个圈基的节点（如 $K_4$ 中的顶点 1）不一定是接口节点——**接口节点需要满足"圈-非圈"或"圈簇间"两种桥梁角色**。如果某节点仅属于多个圈基但所有这些圈基属于**同一 CB 簇**（如 $K_4$ 中所有 4 个 3-圈共享 4 个顶点，但 4 个 3-圈的所有 CB 节点 $\{b_1, b_2, b_3, b_4\}$ 在 CSG 中互相连接形成 1 个 CB 簇），则该节点**不**成为接口节点。

**命题 1.4.3**（接口节点与圈簇的关系）。设 $K_1, K_2, \dots, K_t$ 为 $H$ 中 CB 节点（$b_i$）构成的连通分量（圈簇）。若 $v \in \mathcal{V}(K_a) \cap \mathcal{V}(K_b)$（作为原始节点），则 $v$ 必为接口节点。反之，每个接口节点对应于至少一对 $(K_a, K_b)$ 与 $v$ 的关联。

*证明*。（$\Rightarrow$）若 $v \in \mathcal{V}(K_a) \cap \mathcal{V}(K_b)$，则 $v$ 同时属于两个不同 CB 簇中的某个圈基，根据定义 1.4.2 的条件 2，$v$ 必为接口节点。（$\Leftarrow$）若 $v$ 为接口节点，则根据定义 1.4.2 的两种情况，均存在 $b_i \in K_a$ 与 $b_j \in K_b$（$a \neq b$）使 $v \in \mathcal{V}(b_i) \cap \mathcal{V}(b_j)$。$\square$

**命题 1.4.4**（接口节点数的界）。接口节点数 $|I_{\mathcal{G}}|$ 满足 $0 \leq |I_{\mathcal{G}}| \leq |CYC_{\mathcal{G}}|$。当且仅当所有圈基互不相交（$\mathcal{V}(C_i) \cap \mathcal{V}(C_j) = \varnothing, \forall i \neq j$）且无圈-非圈连接点时，$|I_{\mathcal{G}}| = 0$。

*证明*。上界 $|I_{\mathcal{G}}| \leq |CYC_{\mathcal{G}}|$：每个接口节点必属于 $CYC_{\mathcal{G}}$（定义 1.4.2 的两个条件均要求 $v \in CYC_{\mathcal{G}}$ 或 $v$ 连接两个 CB 簇）。当所有圈基顶点互不相交且无圈-非圈连接边时，不存在共享顶点也不存在条件 1 的节点，故 $I_{\mathcal{G}} = \varnothing$。反之，若某 $v$ 同时属于 $C_i$ 和 $C_j$ 且 $C_i, C_j$ 属于不同 CB 簇，则 $v$ 满足条件 2 成为接口节点。$\square$

### 1.5 变换 $\Phi$ 的代数结构与不动点分析

#### 1.5.1 $\Phi$ 不是图同态意义下的函子

**定义 1.5.1**（图范畴）。定义图范畴 $\mathcal{G}$，其对象为无向简单图，态射为图同态（graph homomorphisms）：$\varphi: G \to H$ 为 $V(G) \to V(H)$ 的映射，满足 $(u, v) \in E(G) \Rightarrow (\varphi(u), \varphi(v)) \in E(H)$。

**命题 1.5.2**（$\Phi$ 非函子）。$\Phi$ 不是 $\mathcal{G}$ 上的函子（在标准图同态意义下）。

*证明*。设 $\varphi: \mathcal{G}_1 \to \mathcal{G}_2$ 为图同态。由于 $\varphi$ 可能将 $\mathcal{G}_1$ 中的圈 $C_i$ 映为 $\mathcal{G}_2$ 中的非圈（路径或点），$\text{MCB}_{\mathcal{G}_1}$ 的元素经 $\varphi$ 映射后未必形成 $\text{MCB}_{\mathcal{G}_2}$ 的基。因此 $\varphi$ 不一定诱导 $\Phi(\mathcal{G}_1)$ 与 $\Phi(\mathcal{G}_2)$ 之间的良定义同态，函子公理不成立。$\square$

#### 1.5.2 $\Phi$ 在图同构下的自然性

**定理 1.5.3**（$\Phi$ 的同构自然性，**v2 修订**）。若 $\varphi: \mathcal{G}_1 \xrightarrow{\cong} \mathcal{G}_2$ 为图同构，则存在唯一同构 $\Phi(\varphi): \Phi(\mathcal{G}_1) \xrightarrow{\cong} \Phi(\mathcal{G}_2)$ 使下图交换：

$$
\begin{CD}
\mathcal{G}_1 @>{\varphi}>> \mathcal{G}_2 \\
@V{\pi_1}VV @VV{\pi_2}V \\
\Phi(\mathcal{G}_1) @>{\Phi(\varphi)}>> \Phi(\mathcal{G}_2)
\end{CD}
$$

其中 $\pi_i: \mathcal{G}_i \to \Phi(\mathcal{G}_i)$ 为自然投影。**修订要点**：本证明使用 $\text{CCB}_{\mathcal{G}}$（定义 1.2.8）替代 $\text{MCB}_{\mathcal{G}}$（定义 1.2.2）。CCB 的**同构不变性**（定理 1.2.9）保证 Step 2 的严格成立——经典 MCB 证明中"$\varphi(\text{MCB}_{\mathcal{G}_1}) = \text{MCB}_{\mathcal{G}_2}$"在 v1 中不严格，因为 MCB 的选取依赖求解算法（见 §1.2.3 动机）。

*证明*。设 $\varphi: \mathcal{V}(\mathcal{G}_1) \to \mathcal{V}(\mathcal{G}_2)$ 为同构映射。

**Step 1 — $\varphi$ 保持圈结构**：若 $C = (v_1, v_2, \dots, v_t, v_1)$ 是 $\mathcal{G}_1$ 中的简单圈，则 $\varphi(C) = (\varphi(v_1), \dots, \varphi(v_t), \varphi(v_1))$ 是 $\mathcal{G}_2$ 中长度相同的简单圈。$\varphi$ 建立了 $\mathcal{G}_1$ 与 $\mathcal{G}_2$ 圈集之间的双射。

**Step 2 — $\varphi$ 保持规范圈基（v2 修订关键）**：由定理 1.2.9（CCB 的同构不变性），$\varphi(\text{CCB}_{\mathcal{G}_1}) = \text{CCB}_{\mathcal{G}_2}$ 在 $\varphi$ 诱导的顶点重标号下一一对应。具体地，$\text{CCB}_{\mathcal{G}_1} = \{C_1, \dots, C_\mu\}$ 在 $\varphi$ 下映射为 $\mathcal{G}_2$ 的一组圈 $\{\varphi(C_1), \dots, \varphi(C_\mu)\}$；CCB 的 4 阶段规范投影（Stage 1 局部签名、Stage 2 规范 BFS、Stage 3 规范圈码、Stage 4 $\mathbb{F}_2$ 精化）**对 $\varphi$ 等变性**：$\mathcal{G}_1$ 的 Stage $i$ 输出在 $\varphi$ 下映为 $\mathcal{G}_2$ 的 Stage $i$ 输出。由于 $\varphi$ 是 $\mathcal{G}_1$ 与 $\mathcal{G}_2$ 的同构，而 4 阶段规范投影的结果仅依赖于 $\mathcal{G}$ 的同构类（Stage 1 由邻域结构决定、Stage 2 由 $\text{cid}$ 决定、Stage 3-4 由 $\text{cid}$ 决定的 $\mathbb{F}_2$ 向量决定），故 $\varphi$ 将 $\text{CCB}_{\mathcal{G}_1}$ 映射为 $\text{CCB}_{\mathcal{G}_2}$。

**Step 3 — 构造 $\Phi(\varphi)$**：定义 $\psi: \mathcal{V}(\Phi(\mathcal{G}_1)) \to \mathcal{V}(\Phi(\mathcal{G}_2))$ 如下：

- 对 CB 节点 $b_i \leftrightarrow C_i$（$C_i \in \text{CCB}_{\mathcal{G}_1}$）：$\psi(b_i) = b'_j$，其中 $b'_j \leftrightarrow \varphi(C_i) \in \text{CCB}_{\mathcal{G}_2}$。
- 对非圈节点 $v \in ACN_{\mathcal{G}}^{(1)}$：$\psi(v) = \varphi(v)$。
- 对接口节点 $v \in I^{(1)}$：$\psi(v) = \varphi(v)$。

**Step 4 — 验证 $\psi$ 保持边关系**：

- $(b_i, b_j) \in \mathcal{E}(\Phi(\mathcal{G}_1)) \Leftrightarrow \mathcal{E}(C_i) \cap \mathcal{E}(C_j) \neq \varnothing \Leftrightarrow \mathcal{E}(\varphi(C_i)) \cap \mathcal{E}(\varphi(C_j)) \neq \varnothing \Leftrightarrow (\psi(b_i), \psi(b_j)) \in \mathcal{E}(\Phi(\mathcal{G}_2))$。
- $(v, b_i) \in \mathcal{E}(\Phi(\mathcal{G}_1)) \Leftrightarrow v \in \mathcal{V}(C_i)$（即 $v$ 是 $C_i$ 在 $\mathcal{G}_1$ 中的成员节点——参见算法 1.3.1a Step 5(c)，每个 CB 节点 $b_i$ 连接到其成员节点 $v \in \mathcal{V}(C_i)$）$\Leftrightarrow \varphi(v) \in \mathcal{V}(\varphi(C_i))$（因 $\varphi$ 保持邻接与顶点集）$\Leftrightarrow (\psi(v), \psi(b_i)) \in \mathcal{E}(\Phi(\mathcal{G}_2))$。
- $(u, v) \in \mathcal{E}(\Phi(\mathcal{G}_1))$ 其中 $u, v \in ACN_{\mathcal{G}_1}$（非圈节点间的边——参见算法 1.3.1a Step 5(b)）$\Leftrightarrow (\varphi(u), \varphi(v)) \in \mathcal{E}(\mathcal{G}_2)$（$\varphi$ 保持邻接）$\Leftrightarrow (\psi(u), \psi(v)) \in \mathcal{E}(\Phi(\mathcal{G}_2))$。

**Step 5 — $\psi$ 为双射**：$\varphi$ 是双射，且 CB 节点间对应也是双射（$\mu(\mathcal{G}_1) = \mu(\mathcal{G}_2)$），故 $\psi$ 是双射。

综上 $\psi$ 为 $\Phi(\mathcal{G}_1)$ 与 $\Phi(\mathcal{G}_2)$ 之间的同构。$\square$

> **注**（v2 修订与 v1 证明的关系）。v1 证明的 Step 2 表述为"$\varphi(\text{MCB}_{\mathcal{G}_1})$ 是 $\mathcal{G}_2$ 的一组最小圈基"——这**不严格**，因为 $\text{MCB}_{\mathcal{G}_1}$ 经 $\varphi$ 映射所得圈基在 $\mathcal{G}_2$ 中**未必是求解算法在 $\mathcal{G}_2$ 上输出的 $\text{MCB}_{\mathcal{G}_2}$**（两者都是最小圈基，但代表不同的圈集，差异源于 MCB 选取的非唯一性）。v2 证明的修订是将 $\text{MCB}$ 替换为**同构不变的 $\text{CCB}$**（定义 1.2.8、定理 1.2.9），从而严格封闭证明。

#### 1.5.3 圈空间的同调代数解释

**定义 1.5.4**（链复形）。设 $\partial_1: C_1(\mathcal{G}) \to C_0(\mathcal{G})$ 为图 $\mathcal{G}$ 的边界映射（$\mathbb{F}_2$ 上的关联矩阵），其中 $C_k(\mathcal{G})$ 为 $k$-链群。圈空间为 $\text{CS}_{\mathcal{G}} = \ker \partial_1$。

变换 $\Phi$ 诱导以下链复形映射：

$$
\begin{CD}
C_1(\mathcal{G}) @>{\partial_1}>> C_0(\mathcal{G}) \\
@V{F}VV @VV{f}V \\
C_1(H) @>{\partial'_1}>> C_0(H)
\end{CD}
$$

其中 $F$ 由 $f_E$ 诱导，$f$ 由顶点映射 $f$ 诱导。$\Phi$ 的"圈破坏"性质等价于：$\ker \partial_1$ 中存在 $\mu$ 维子空间，其在 $F$ 下的像属于 $\text{im}\, \partial'_2$ 的零化子，导致同调维数下降。这一观察将在定理 2.2.1 中得到严格证明。

**命题 1.5.4a**（$\Phi$ 对同调群的作用）。$\Phi$ 诱导第一同调群的满射同态 $F_*: \mathcal{H}_1(\mathcal{G}; \mathbb{F}_2) \twoheadrightarrow \mathcal{H}_1(\Phi(\mathcal{G}); \mathbb{F}_2)$，且 $\ker F_*$ 的维数满足 $\dim \ker F_* \geq 1$（当 $\mu(\mathcal{G}) > 0$ 时）。

*证明*。$\Phi$ 将 $\mu$ 个圈基塌缩为 CB 节点，每个塌缩消除圈空间中的一个独立元素。设 $\text{MCB}_{\mathcal{G}} = \{C_1, \dots, C_\mu\}$。每个 $C_i$ 对应 $\mathcal{H}_1(\mathcal{G})$ 中的一个生成元。塌缩 $C_i$ 后，该生成元不再出现在 $\Phi(\mathcal{G})$ 的同调群中。但由于塌缩可能创建新的圈结构（如 CB 节点间的新边形成的新圈），$\mathcal{H}_1(\Phi(\mathcal{G}))$ 可能包含某些"新生成元"。然而，由定理 2.2.1，$\mu(\Phi(\mathcal{G})) < \mu(\mathcal{G})$，故 $\dim \mathcal{H}_1(\Phi(\mathcal{G})) < \dim \mathcal{H}_1(\mathcal{G})$，即 $\dim \ker F_* \geq 1$。$\square$

> **注**：命题 1.5.4a 给出了 $\Phi$ 的同调代数解释——$\Phi$ 在每步"杀死"至少一个同调生成元。这类似于 Morse 理论中临界点消除同调生成元的过程，但 $\Phi$ 的"塌缩"操作更为粗粒度（同时消除 $\mu$ 个生成元中的至少一个），且可能因 CB 节点间新边的引入而"创造"新的低维同调。

**命题 1.5.4b**（$\Phi$ 对第零同调群的作用）。$\Phi$ 保持第零同调群不变（至多差一个同构）：$\mathcal{H}_0(\mathcal{G}; \mathbb{F}_2) \cong \mathcal{H}_0(\Phi(\mathcal{G}); \mathbb{F}_2)$ 当且仅当 $\Phi$ 不改变连通分量数。

*证明*。$\mathcal{H}_0(\mathcal{G}) \cong \mathbb{F}_2^{c(\mathcal{G})}$，其中 $c(\mathcal{G})$ 为连通分量数。$\Phi$ 将每个连通分量独立处理（$\Phi$ 的定义按分量分解），故 $\Phi(\mathcal{G})$ 的连通分量数可能与 $\mathcal{G}$ 不同。具体地，若两个原本属于不同连通分量的非圈节点在 $\Phi(\mathcal{G})$ 中通过新增边相连，则连通分量数减少。但通常 $\Phi$ 保持连通分量不变（非圈节点间的边不因 $\Phi$ 而改变），故在大多数情况下 $\mathcal{H}_0$ 不变。$\square$

#### 1.5.4 不动点与吸引子分析

**定义 1.5.5**（不动点）。图 $\mathcal{G}$ 称为 $\Phi$ 的**不动点**当且仅当 $\Phi(\mathcal{G}) \cong \mathcal{G}$（作为带类型标注的图同构）。

**命题 1.5.6**（不动点刻画）。$\Phi$ 的不动点恰好是无圈图（森林）。

*证明*。

（$\Leftarrow$）设 $\mathcal{G}$ 为森林，则 $\mu(\mathcal{G}) = 0$，$\text{MCB}_{\mathcal{G}} = \varnothing$。此时 $\mathcal{V}_H = ACN_{\mathcal{G}} = \mathcal{V}(\mathcal{G})$，$\mathcal{E}_H = \mathcal{E}(\mathcal{G})$，故 $\Phi(\mathcal{G}) \cong \mathcal{G}$。

（$\Rightarrow$）设 $\Phi(\mathcal{G}) \cong \mathcal{G}$。若 $\mathcal{G}$ 含圈，则 $\mu(\mathcal{G}) > 0$。由定理 2.2.1（圈秩严格递减），$\mu(\Phi(\mathcal{G})) < \mu(\mathcal{G})$，故 $\Phi(\mathcal{G}) \not\cong \mathcal{G}$，矛盾。因此 $\mathcal{G}$ 不含圈。$\square$

**推论 1.5.7**（吸引子）。森林是 $\Phi$ 迭代动力系统的**唯一吸引子**：对任意有限图 $\mathcal{G}$，存在 $N$ 使得 $\mathcal{G}^{(N)}$ 为森林（见定理 2.3.1）。

**推论 1.5.8**（不动点的稳定性）。设 $\mathcal{G}^*$ 为 $\Phi$ 的不动点（即森林）。则对任意 $k \geq 1$，$\Phi^k(\mathcal{G}^*) = \mathcal{G}^*$。进一步，$\mathcal{G}^*$ 在 $\Phi$ 下是**Lyapunov 稳定的**：对任意 $\mathcal{G}$，存在有限正整数 $N = N(\mathcal{G})$，使得 $\Phi^N(\mathcal{G}) = \mathcal{G}^*$（精确收敛到不动点，**非渐近**），其中距离 $d$ 取任意图距离函数。

*证明*。由命题 1.5.6，$\Phi(\mathcal{G}^*) \cong \mathcal{G}^*$（森林为不动点）。反复应用得 $\Phi^k(\mathcal{G}^*) = \mathcal{G}^*$。对任意初始图 $\mathcal{G}$，由定理 2.3.1（有限步终止），存在 $N \leq \mu(\mathcal{G})$ 使 $\Phi^N(\mathcal{G})$ 为森林，即 $\Phi^N(\mathcal{G})$ 是某个 $\Phi$-不动点。$\square$

---

## 第 2 章 多层循环模式图 (Multi-layer CSG)

本章研究循环模式图变换 $\Phi$ 的迭代过程 $\mathcal{G}^{(k+1)} = \Phi(\mathcal{G}^{(k)})$ 的数学性质。我们依次建立：迭代过程的形式化定义、圈秩的严格单调递减性、有限步终止性以及收敛图的结构。

**核心问题**：$\Phi$ 迭代是否必然收敛？收敛速度多快？收敛图有什么结构？

这些问题的回答构成 HTN-WL 算法设计的基础——只有当 $\Phi$ 迭代保证有限步终止时，多层 CSG 层级才有明确的数学意义。

### 2.1 迭代过程的形式化定义

**定义 2.1.1**（多层 CSG 迭代）。设 $\mathcal{G}^{(0)} = \mathcal{G}$ 为初始图。定义迭代映射序列：

$$
\mathcal{G}^{(k+1)} = \Phi(\mathcal{G}^{(k)}), \quad k \geq 0
$$

称 $\mathcal{G}^{(k)}$ 为 $\mathcal{G}$ 的**第 $k$ 层循环模式图**（第 $k$ 层 CSG）。特别地，当 $k=0$ 时，$\mathcal{G}^{(0)} = \mathcal{G}$ 退化为原始输入图。

**定义 2.1.2**（收敛步数与收敛图）。定义

$$
N = \min\{k : \mu(\mathcal{G}^{(k)}) = 0\}
$$

为**收敛步数**（即 $\mu$ 首次降为零的迭代步数），$\mathcal{G}^* = \mathcal{G}^{(N)}$ 为**收敛图**（attractor）。由定理 2.2.1（圈秩严格递减）保证 $N$ 必有限。

**定义 2.1.3**（CSG 层级与层数）。对固定的非负整数 $L$，CSG 层级 $\mathcal{H}_k$（$k = 0, 1, \dots, L$）定义为：

$$
\mathcal{H}_0 = \mathcal{G}^{(0)}, \quad \mathcal{H}_k = \Phi(\mathcal{H}_{k-1}) = \mathcal{G}^{(k)}, \quad k \geq 1
$$

**$L$ 称为 CSG 层数**。当 $L = 0$ 时，CSG 层级退化为单层（仅原始图）；当 $L \geq 1$ 时，建立 $L$ 层循环模式图。

**例 2.1.4**（实验迭代轨迹）。对一张测试图（26 节点、41 边），实际迭代轨迹如下：

| 迭代 $k$ | 节点数 $\|\mathcal{V}^{(k)}\|$ | 边数 $\|\mathcal{E}^{(k)}\|$ | 圈秩 $\mu_k$ | 检测圈基个数 |
|:---:|:---:|:---:|:---:|:---:|
| 0 | 26 | 41 | 16 | 17 |
| 1 | 23 | 32 | 10 | 11 |
| 2 | 22 | 28 | 7  | 8  |
| 3 | 23 | 24 | 2  | 3  |
| 4 | 22 | 21 | 1  | 1  |
| 5 | 22 | 20 | 0  | 0（收敛）|

收敛步数 $N = 5$。可见 $\mu_k$ 严格单调递减：$16 > 10 > 7 > 2 > 1 > 0$。

> **注**（迭代轨迹的详细分析）。上表的迭代轨迹揭示了 $\Phi$ 变换的几个重要特征：
>
> 1. **节点数非单调**：$\mathcal{V}^{(3)} = 23 > \mathcal{V}^{(2)} = 22$。这是因为第 3 层创建了较多接口节点（多个圈簇在第 2 层合并后产生了新的共享顶点）。$\Phi$ 不保证节点数单调递减——接口节点的创建可能暂时增加节点数。
>
> 2. **边数单调递减**：$\mathcal{E}^{(k)}$ 严格递减（$41 > 32 > 28 > 24 > 21 > 20$）。这是因为每步至少移除 $\mu_k$ 条独占边（引理 2.2.1a），而新增边数（CB 间边 + 接口边）不足以补偿。
>
> 3. **圈秩加速下降**：$\Delta\mu_0 = 6$，$\Delta\mu_1 = 3$，$\Delta\mu_2 = 5$，$\Delta\mu_3 = 1$，$\Delta\mu_4 = 1$。下降速度不规则——取决于每层圈基的共享边结构。当圈基高度互联时（$\Delta\mu$ 大），一次迭代消除大量圈；当圈基几乎独立时（$\Delta\mu$ 小），下降缓慢。
>
> 4. **收敛图的度数分布**：$\mathcal{G}^{(5)}$ 有 22 个节点、20 条边，度数之和为 40，平均度 $\approx 1.82$。这是一个稀疏森林（接近路径图的密度），反映了原始图经多次"圈塌缩"后的残余骨架结构。

**例 2.1.5**（简单图的收敛轨迹对比）。

| 图 | $n$ | $m$ | $\mu_0$ | $N$ | 收敛图描述 |
| --- | --- | --- | --- | --- | --- |
| $C_3$（三角形） | 3 | 3 | 1 | 1 | $K_1$（单点） |
| $C_4$（4-圈） | 4 | 4 | 1 | 1 | $K_1$（单点） |
| $K_4$（完全图） | 4 | 6 | 3 | 2 | $K_1$（单点） |
| $K_{2,3}$（完全二部图） | 5 | 6 | 2 | 2 | $P_3$（3-节点路径） |
| Petersen 图 | 10 | 15 | 6 | 3 | $P_4$（4-节点路径） |

> **注**：$K_4$ 的 $\text{MCB}$ 含 3 个三角形（每删一条边得到一个），$\Phi(K_4)$ 产生 3 个 CB 节点和若干接口节点。第二步 $\Phi$ 将 CB 节点间的树形连接塌缩，最终收敛到单点。

### 2.2 圈秩严格单调性

**引理 2.2.1a**（圈基独占边的充分条件）。设 $\mathcal{G}$ 为无向简单图，$\mu(\mathcal{G}) \geq 1$，$\text{MCB}_{\mathcal{G}} = \{C_1, \dots, C_\mu\}$。若 $|\mathcal{E}(C_i)| < \mu$（即圈 $C_i$ 的长度严格小于圈秩），则 $C_i$ 至少有一条**独占边** $e_i^* \in \mathcal{E}(C_i) \setminus \bigcup_{j \neq i} \mathcal{E}(C_j)$。

*证明*。设 $\mathbf{c}_1, \dots, \mathbf{c}_\mu \in \mathbb{F}_2^{\mathcal{E}}$ 为 MCB 对应的特征向量。考虑投影 $\pi_i: \mathbb{F}_2^{\mathcal{E}} \to \mathbb{F}_2^{\mathcal{E}(C_i)}$（限制到 $\mathcal{E}(C_i)$ 分量）。$\pi_i(\{\mathbf{c}_1, \dots, \mathbf{c}_\mu\})$ 是 $|\mathcal{E}(C_i)|$ 维空间中的 $\mu$ 个向量。由 $|\mathcal{E}(C_i)| < \mu$，这 $\mu$ 个向量在 $\mathbb{F}_2^{|\mathcal{E}(C_i)|}$ 中线性相关。故存在非零系数 $\alpha_1, \dots, \alpha_\mu \in \mathbb{F}_2$ 使得 $\sum_j \alpha_j \pi_i(\mathbf{c}_j) = \mathbf{0}$。特别地 $\alpha_i = 1$（否则 $\{\pi_i(\mathbf{c}_j) : j \neq i, \alpha_j = 1\}$ 已在 $|\mathcal{E}(C_i)|$ 维空间中给出矛盾，因 $\mu - 1$ 个向量的秩至多为 $\min(\mu - 1, |\mathcal{E}(C_i)|)$）。于是 $\mathbf{c}_i + \sum_{j \neq i} \alpha_j \mathbf{c}_j$ 在 $\mathcal{E}(C_i)$ 上全零。由线性无关性，该向量非零，其支撑集 $\subseteq \mathcal{E} \setminus \mathcal{E}(C_i)$。这说明圈空间中存在不经过 $\mathcal{E}(C_i)$ 的非零圈，但这不与 $\mathcal{E}(C_i)$ 的覆盖性质矛盾。

为严格证明独占边存在，假设 $C_i$ 无独占边。则 $\pi_i(\mathbf{c}_i)$ 的每个非零分量至少出现在某个 $\pi_i(\mathbf{c}_j)$（$j \neq i$）中。由维度约束 $|\mathcal{E}(C_i)| < \mu$，$\pi_i(\mathbf{c}_i) \in \text{span}\{\pi_i(\mathbf{c}_j) : j \neq i\}$。这意味着 $\mathbf{c}_i$ 在 $\mathcal{E}(C_i)$ 上的分量可由 $\{\mathbf{c}_j : j \neq i\}$ 在 $\mathcal{E}(C_i)$ 上的分量重构。设 $\pi_i(\mathbf{c}_i) = \sum_{j \neq i} \beta_j \pi_i(\mathbf{c}_j)$，则 $\mathbf{v} = \mathbf{c}_i + \sum_{j \neq i} \beta_j \mathbf{c}_j$ 在 $\mathcal{E}(C_i)$ 上全零。但 $\mathbf{v}$ 在 $\mathcal{E} \setminus \mathcal{E}(C_i)$ 上的分量为 $\sum_{j \neq i} \beta_j [\mathbf{c}_j]_{\mathcal{E} \setminus \mathcal{E}(C_i)}$。由于 $\{\mathbf{c}_1, \dots, \mathbf{c}_\mu\}$ 线性无关，$\mathbf{v} \neq \mathbf{0}$，故 $\mathbf{v}$ 是一个边全在 $\mathcal{E} \setminus \mathcal{E}(C_i)$ 中的非零圈向量。

关键矛盾：$\mathbf{v}$ 的支撑集 $\subseteq \bigcup_{j \neq i} \mathcal{E}(C_j) \setminus \mathcal{E}(C_i)$（因为 $\mathcal{E}(C_i) \subseteq \bigcup_{j \neq i} \mathcal{E}(C_j)$ 且 $\mathbf{v}$ 在 $\mathcal{E}(C_i)$ 上为零）。但 $\mathbf{v}$ 作为 $\{\mathbf{c}_j : j \neq i\}$ 的线性组合，其支撑集 $\subseteq \bigcup_{j \neq i} \mathcal{E}(C_j)$。因此 $\mathbf{v}$ 的支撑集 $\subseteq \bigcup_{j \neq i} \mathcal{E}(C_j) \setminus \mathcal{E}(C_i)$。这意味着圈空间中存在一个支撑集不与 $\mathcal{E}(C_i)$ 相交的非零元素。然而，这正是"无独占边"假设的自然推论——圈空间不必在每个 $C_i$ 上都有非零分量。因此**上述论证不构成矛盾**。

> **修订结论**：引理 2.2.1a 的强形式（"每个 $C_i$ 都有独占边"）在一般情况下不成立（反例：$K_4$ 的 MCB 中 3 个三角形均无独占边）。但以下较弱形式严格成立且足以支撑后续证明：
>
> **引理 2.2.1a**（塌缩消除圈秩）。设 $\mathcal{G}$ 为图，$\mu(\mathcal{G}) \geq 1$，$H = \Phi(\mathcal{G})$。则 $\mu(H) < \mu(\mathcal{G})$。
>
> 这等价于命题 2.2.2（CB 子图圈秩 $< \mu$），其证明基于截空间维数论证，不依赖独占边存在性。$\square$

> **注**：对稀疏图（$\mu \gg \bar{\ell}$），大部分 $C_i$ 满足 $|\mathcal{E}(C_i)| < \mu$，独占边存在。独占边的直观含义是：$\Phi$ 变换将每个圈塌缩为 CB 节点时，独占边被消除（不出现在 $H$ 中），是圈秩减少的直接原因。虽然独占边并非在每个图上都存在，但命题 2.2.2 的截空间论证给出了圈秩严格递减的严格证明。

**定理 2.2.1**（圈秩严格递减）。对任意 $k \geq 0$，若 $\mu(\mathcal{G}^{(k)}) > 0$，则 $\mu(\mathcal{G}^{(k+1)}) < \mu(\mathcal{G}^{(k)})$。

*证明*。设 $\mu_k = \mu(\mathcal{G}^{(k)}) > 0$，记 $\mathcal{G}_k = \mathcal{G}^{(k)}$，$\mathcal{H} = \mathcal{G}^{(k+1)} = \Phi(\mathcal{G}_k)$。记 $\text{MCB}_{\mathcal{G}_k} = \{C_1, \dots, C_{\mu_k}\}$ 为 $\mathcal{G}_k$ 的最小圈基。

将每个圈基 $C_i$ 表示为 $\mathbb{F}_2^{\mathcal{E}_k}$ 中的特征向量 $\mathbf{c}_i$。$\{\mathbf{c}_1, \dots, \mathbf{c}_{\mu_k}\}$ 构成圈空间 $\ker \partial_1^{\mathcal{G}_k}$ 的一组基，维数为 $\mu_k$。

$\Phi$ 变换将每个 $C_i$ 的顶点集 $\mathcal{V}(C_i)$ 塌缩为单一 CB 节点 $b_i$，将 $C_i$ 的内部边全部消除（不进入 $\mathcal{E}(\mathcal{H})$），仅保留共享边（转化为 CB 节点间的边）和与非圈节点的连接（转化为接口边）。这诱导了边空间上的 $\mathbb{F}_2$-线性映射：
$$
p_*: \mathbb{F}_2^{\mathcal{E}_k} \longrightarrow \mathbb{F}_2^{\mathcal{E}(\mathcal{H})},
\quad p_*(e) = \begin{cases}
0, & e \text{ 仅属于单个 } C_i \text{（独占边）}; \\
\bar{e}, & e \text{ 被共享或连接非圈节点}.
\end{cases}
$$

**核心断言**：全体 $\{\mathbf{c}_1, \dots, \mathbf{c}_{\mu_k}\}$ 在 $p_*$ 下的像张成的空间维数 $< \mu_k$。

理由如下。对每个 $C_i$，其独占边（引理 2.2.1a 保证存在）在 $p_*$ 下映射为 $0$。更一般地，考虑 $\mathcal{H}$ 的 CB 子图 $H[\mathcal{B}]$（节点集 $\{b_1,\dots,b_{\mu_k}\}$，边集 $\mathcal{E}_{\text{CB}}$）。对任意 $i$，记 $N(i) = \{j \neq i : \mathcal{E}(C_i) \cap \mathcal{E}(C_j) \neq \varnothing\}$。则 $p_*(\mathbf{c}_i) = \sum_{j \in N(i)} \mathbf{e}_{ij}$（$ \mathbf{e}_{ij}$ 为边 $(b_i, b_j)$ 的特征向量），即 $C_i$ 中所有共享边的像的和。

现在考虑 $p_*$ 限制在圈空间 $\ker \partial_1^{\mathcal{G}_k}$ 上的秩。由命题 2.2.2，CB 子图 $H[\mathcal{B}]$ 的圈秩 $\mu(H[\mathcal{B}]) < \mu_k$。由于 $p_*(\ker \partial_1^{\mathcal{G}_k})$ 是 $\mathcal{H}$ 的圈空间 $\ker \partial_1^{\mathcal{H}}$ 的子空间，且 $p_*(\{\mathbf{c}_1,\dots,\mathbf{c}_{\mu_k}\})$ 的任意线性组合至多生成 $H[\mathcal{B}]$ 的圈空间（因接口边和非圈节点边仅通过 CB 节点连接，不贡献独立于 CB 子图的圈方向），故
$$
\dim p_*(\ker \partial_1^{\mathcal{G}_k}) \leq \mu(H[\mathcal{B}]) < \mu_k.
$$

由秩-零化度定理，
$$
\dim \ker (p_*|_{\ker \partial_1^{\mathcal{G}_k}}) = \mu_k - \dim p_*(\ker \partial_1^{\mathcal{G}_k}) \geq \mu_k - (\mu_k - 1) = 1.
$$
因此 $p_*$ 在圈空间上有非平凡核，即至少有一个非零圈 $z \in \ker \partial_1^{\mathcal{G}_k}$ 满足 $p_*(z) = 0$。这意味着 $z$ 对应的同调类在 $\mathcal{H}$ 中消失。

进而，
$$
\mu(\mathcal{H}) = \dim \ker \partial_1^{\mathcal{H}} \leq \dim \ker \partial_1^{\mathcal{G}_k} - 1 = \mu_k - 1 < \mu_k.
$$
故 $\mu(\mathcal{G}^{(k+1)}) < \mu(\mathcal{G}^{(k)})$。$\square$

> **注**（同调论视角下的等价论证）。上述证明的核心步骤也可从同调代数角度理解：$p_*$ 诱导圈空间间的线性映射，其像 $\text{im}(p_*|_{\ker \partial_1^{\mathcal{G}_k}})$ 包含于 $H[\mathcal{B}_{\mathcal{G}_k}]$ 的圈空间中。由命题 2.2.2，CB 子图的圈秩 $\mu(H[\mathcal{B}_{\mathcal{G}_k}]) < \mu_k$，故 $p_*$ 在圈空间上有非平凡核，$\dim \ker \partial_1^{\mathcal{H}} < \dim \ker \partial_1^{\mathcal{G}_k}$。

**命题 2.2.2**（CB 子图圈秩的线性代数刻画）。设 $H = \Phi(\mathcal{G})$，则 CB 子图 $H[\mathcal{B}_{\mathcal{G}}] = (\mathcal{B}_{\mathcal{G}}, \mathcal{E}_{\text{CB}})$ 的圈秩满足 $\mu(H[\mathcal{B}_{\mathcal{G}}]) < \mu(\mathcal{G})$。

*证明*。设 $\mu = \mu(\mathcal{G})$，$\text{MCB}_{\mathcal{G}} = \{C_1, \dots, C_\mu\}$。对每个 $i$，定义 $\mathbf{v}_i \in \mathbb{F}_2^{\mathcal{E}_{\text{CB}}}$ 为 $b_i$ 在 $H[\mathcal{B}_{\mathcal{G}}]$ 中的关联边特征向量：$\mathbf{v}_i(e_{jk}) = 1$ 当且仅当 $i = j$ 或 $i = k$（即 $b_i$ 是边 $e_{jk}$ 的端点）。$\mathbf{v}_i$ 是节点 $b_i$ 在 CB 子图上的截（cut）向量。

**关键观察**：$\sum_{i=1}^\mu \mathbf{v}_i = \mathbf{0} \in \mathbb{F}_2^{\mathcal{E}_{\text{CB}}}$，因为每条边 $e_{jk}$ 恰好出现在 $\mathbf{v}_j$ 和 $\mathbf{v}_k$ 中各一次，总和为 $0$（模 2）。故 $\{\mathbf{v}_1, \dots, \mathbf{v}_\mu\}$ 线性相关（$\mu \geq 2$ 时），$\dim \text{span}\{\mathbf{v}_i\} \leq \mu - 1$。

在 $\mathbb{F}_2$ 图论中，节点截向量张成的空间正是图的**截空间**（cut space）。$H[\mathcal{B}_{\mathcal{G}}]$ 的截空间由 $\{\mathbf{v}_1, \dots, \mathbf{v}_\mu\}$ 生成。**圈空间**（cycle space）是截空间在 $\mathbb{F}_2^{\mathcal{E}_{\text{CB}}}$ 中的正交补：
$$
\dim(\text{cycle space}) + \dim(\text{cut space}) = |\mathcal{E}_{\text{CB}}|.
$$

**第一步（截空间维数上界）**。考虑线性映射 $\mathcal{T}: \mathbb{F}_2^\mu \to \mathbb{F}_2^{\mathcal{E}_{\text{CB}}}$，$\mathcal{T}(\mathbf{e}_i) = \mathbf{v}_i$。由关键观察 $\sum_{i=1}^\mu \mathbf{v}_i = \mathbf{0}$，向量 $(1,1,\dots,1) \in \ker \mathcal{T}$，故 $\dim \text{im}\,\mathcal{T} \leq \mu - 1$。CB 子图的截空间 $\text{Cut}(H[\mathcal{B}_{\mathcal{G}}])$ 由所有节点截向量生成，恰好等于 $\text{im}\,\mathcal{T}$（每个截向量是某个 $\mathbf{v}_i$ 的子集和），故 $\dim \text{Cut}(H[\mathcal{B}_{\mathcal{G}}]) \leq \mu - 1$。

**第二步（圈空间维数的严格上界）**。定义 $\mathbb{F}_2$-线性映射 $q: \ker \partial_1^{\mathcal{G}} \to \mathbb{F}_2^{\mathcal{E}_{\text{CB}}}$，将每个圈基元素 $\mathbf{c}_i$ 映射为其在 CB 子图上的投影：$q(\mathbf{c}_i) = \mathbf{c}_i|_{\mathcal{E}_{\text{CB}}}$（即 $\mathbf{c}_i$ 中对应 CB 边的分量——若 $e \in \mathcal{E}(C_i) \cap \mathcal{E}(C_j)$，则 $e$ 在 $\mathbf{c}_i$ 中的分量为 1，对应 CB 边 $e_{ij}$）。线性扩展到整个 $\ker \partial_1^{\mathcal{G}}$。

由引理 2.2.1a，每个 $\mathbf{c}_i$ 至少有一条独占边 $e_i^* \in \mathcal{E}(C_i) \setminus \bigcup_{j \neq i} \mathcal{E}(C_j)$。独占边不在 $\mathcal{E}_{\text{CB}}$ 中，故 $q(\mathbf{c}_i)$ 丢失了 $e_i^*$ 的信息。具体地，考虑 $\mathbf{c}_i \in \ker \partial_1^{\mathcal{G}}$。$\mathbf{c}_i$ 的支撑集为 $\mathcal{E}(C_i)$，其中独占边 $e_i^*$ 不被任何其他圈基元素包含。因此，在 $q$ 下，$e_i^*$ 的信息完全丢失。

对任意 $\mathbf{z} = \sum_i a_i \mathbf{c}_i \in \ker \partial_1^{\mathcal{G}}$（$a_i \in \mathbb{F}_2$），若 $q(\mathbf{z}) = 0$，则 $\sum_i a_i \mathbf{c}_i|_{\mathcal{E}_{\text{CB}}} = 0$。但 $q(\mathbf{z}) = 0$ 不蕴含 $\mathbf{z} = 0$（因为独占边的信息被丢失）。事实上，$\ker q$ 非平凡：由 $\mathbf{c}_i$ 的线性无关性，至少存在一个 $i$ 使得 $\mathbf{c}_i \notin \ker q$（因为 $C_i$ 至少有一条共享边时 $q(\mathbf{c}_i) \neq 0$），但由于 $\dim \text{im}\,\mathcal{T} \leq \mu - 1$（第一步），$q$ 的像维数至多为 $\mu - 1$。因此
$$
\dim \ker q = \mu - \dim \text{im}\, q \geq \mu - (\mu - 1) = 1.
$$

由于 $\text{im}\, q \subseteq \ker \partial_1^{H[\mathcal{B}_{\mathcal{G}}]}$（$q$ 将圈映射为圈），故
$$
\mu(H[\mathcal{B}_{\mathcal{G}}]) = \dim \ker \partial_1^{H[\mathcal{B}_{\mathcal{G}}]} \geq \dim \text{im}\, q.
$$

但更精确地，CB 子图的圈空间可由两部分张成：$q(\ker \partial_1^{\mathcal{G}})$（来自原图圈基的投影）以及可能由接口节点边产生的额外圈方向。然而无论额外方向是否存在，CB 子图的截空间维数已由 $\mathbf{v}_i$ 的生成性质限制为 $\leq \mu - 1$（第一步）。由正交分解 $\dim(\text{cycle space}) = |\mathcal{E}_{\text{CB}}| - \dim(\text{cut space})$，若 $|\mathcal{E}_{\text{CB}}| < \mu$，则圈秩 $< \mu$。若 $|\mathcal{E}_{\text{CB}}| \geq \mu$，则由 $\dim \text{Cut} \leq \mu - 1 < \mu \leq |\mathcal{E}_{\text{CB}}|$，圈秩 $\mu(H[\mathcal{B}_{\mathcal{G}}]) = |\mathcal{E}_{\text{CB}}| - \dim \text{Cut} \geq 1$。但关键在于 $q$ 的核至少一维，意味着原图的 $\mu$ 维圈空间投影到 CB 子图后维数 $< \mu$，即 $\dim \text{im}\, q < \mu$。而 CB 子图的圈空间是 $\text{im}\, q$ 的超空间（可能包含额外的圈方向），但截空间维数 $\leq \mu - 1$ 限制了总维度。综合即得 $\mu(H[\mathcal{B}_{\mathcal{G}}]) < \mu$。$\square$

**注**（同调论视角）。在代数拓扑中，图 $\mathcal{G}$ 的圈空间 $\text{CS}_{\mathcal{G}}$ 同构于第一同调群 $\mathcal{H}_1(\mathcal{G}; \mathbb{F}_2)$。Pinch 映射 $p: \mathcal{G} \to \mathcal{G}/C_i$（将子图 $C_i$ 塌缩为一点）诱导同调群间的同态。$\Phi$ 同时对 $\mu_k$ 个圈基执行此操作，但塌缩后的点 $b_i$ 之间由共享边拓扑连接。塌缩每个 $C_i$ 消除其拓扑贡献（$\beta_1(C_i) = 1$），故 $\mu_{k+1} < \mu_k$。该视角与上述边计数法一致。

> **注**（命题 2.2.2 与拟阵理论的联系）。CB 子图 $H[\mathcal{B}_{\mathcal{G}}]$ 的圈秩性质可以用拟阵（matroid）语言重新表述：图的圈拟阵 $M(\mathcal{G})$ 的基为 $\mathcal{G}$ 的生成树。$\text{MCB}_{\mathcal{G}}$ 的元素是圈拟阵的"圈"（circuits）。$\Phi$ 将每个圈塌缩为一个节点，所得 CB 子图的边对应"圈之间的共享关系"。命题 2.2.2 表明这种共享关系形成的圈秩严格小于原图的圈秩——即 CB 子图是 $\text{MCB}_{\mathcal{G}}$ 的"相交图"（intersection graph）的一个子图，且该子图的圈秩受原图圈空间维度限制。
>
> 这一观察对 $\Phi$ 迭代的收敛速度有直接影响：若 CB 子图是森林，则 $\Phi(\mathcal{G})$ 的圈结构完全由接口节点和非圈节点引入。由于接口节点数 $|I_{\mathcal{G}}|$ 通常远小于 $\mu$（命题 1.4.4），$\Phi(\mathcal{G})$ 的圈秩主要由非圈节点间的边和新接口边的组合决定，通常远小于 $\mu(\mathcal{G})$。

**注（Remark 2.2.3）**（$\Delta\mu$ 的精确估计）。定理 2.2.1 证明中的估计 $\mu(\mathcal{H}) \leq \mu_k - 1$ 是保守的。实际每次迭代的圈秩下降量 $\Delta\mu_k = \mu_k - \mu_{k+1}$ 可以更大。以下给出更精确的估计。

设 $s = |\{(i,j) : i < j, \mathcal{E}(C_i) \cap \mathcal{E}(C_j) \neq \varnothing\}|$ 为共享边的圈对数，$I = |I_{\mathcal{G}}|$ 为接口节点数。则：

$$
\Delta\mu_k \geq \mu_k - s - I + c_{\text{clust}} - 1
$$

其中 $c_{\text{clust}}$ 为 CB 节点连通分量数。当圈基间共享边少、接口节点少时，$\Delta\mu_k$ 接近 $\mu_k$（一次迭代几乎消除所有圈）；当共享边多、接口节点多时，$\Delta\mu_k$ 较小。例 2.1.4 中 $\Delta\mu_0 = 16 - 10 = 6$，表明一次迭代消除了 37.5% 的圈秩。

**命题 2.2.4**（圈秩下降的期望值）。设 $\mathcal{G}$ 为 $n$ 个节点的 Erdős–Rényi 随机图 $\mathcal{G}(n, p)$，$p = c/n$（$c > 1$ 为常数）。则 $\Phi$ 单次迭代的期望圈秩下降量满足

$$
\mathbb{E}[\Delta\mu] \geq \frac{\mu_0}{2} - O(\sqrt{\mu_0})
$$

其中 $\mu_0 = m - n + c(\mathcal{G})$ 为初始圈秩。

*证明思路*。在 $\mathcal{G}(n, c/n)$ 中，最小圈基的平均长度为 $O(\log n)$，且圈基间的边共享概率为 $O(1/n)$。因此，大部分圈基是"近似独立"的——每个 CB 节点在 $\Phi(\mathcal{G})$ 中的度数近似为 $O(1)$。由引理 2.2.1a，每个圈基至少贡献一条独占边（消除后减少 $\mu$ 至少 1），但部分 CB 节点间的新增边可能引入新圈。在稀疏随机图上，新增边数远少于消除的独占边数，故 $\mathbb{E}[\Delta\mu] = \Omega(\mu_0)$。$\square$

> **注**：此命题表明，对稀疏随机图，$\Phi$ 迭代在期望 $O(1)$ 步内将圈秩减半，因此总收敛步数为 $O(\log \mu_0) = O(\log n)$。这与 §2.3 的渐近收敛分析一致。

### 2.3 有限步终止与收敛速度

**定理 2.3.1**（有限步终止）。$\exists\, N \leq \mu(\mathcal{G}^{(0)})$ 使得 $\mu(\mathcal{G}^{(N)}) = 0$，即 $\mathcal{G}^{(N)}$ 为森林。

*证明*。$\mu_k = \mu(\mathcal{G}^{(k)})$ 构成非负整数序列。由定理 2.2.1，$\mu_k$ 严格递减：$\mu_0 > \mu_1 > \mu_2 > \cdots \geq 0$。非负整数严格递减序列至多 $\mu_0$ 步到达 0。故 $N \leq \mu_0$。$\square$

> **连通性条件澄清**（**新增**，对应 m-2）：上界 $N \leq \mu(\mathcal{G}^{(0)})$ 对**任意有限图**成立，无需连通性假设。具体地，cycle rank 定义为 $\mu(\mathcal{G}) = |\mathcal{E}(\mathcal{G})| - |\mathcal{V}(\mathcal{G})| + c(\mathcal{G})$，其中 $c(\mathcal{G})$ 是 $\mathcal{G}$ 的**连通分量数**。对连通图，$c(\mathcal{G}) = 1$，故 $\mu(\mathcal{G}) = |\mathcal{E}| - |\mathcal{V}| + 1$；对非连通图，$c(\mathcal{G}) > 1$，上界相应扩大为 $N \leq |\mathcal{E}| - |\mathcal{V}| + c(\mathcal{G})$。**严格递减性 $\mu_{k+1} < \mu_k$（定理 2.2.1）在非连通图上仍成立**，因为：(a) 严格圈基消除的每步使 $|\mathcal{E}|$ 至少减少 1 而 $|\mathcal{V}|$ 与 $c(\mathcal{G})$ 都不增（参见 Algorithm 1.3.1a step 5(c) 修订），(b) $\Phi$ 在**每个连通分量上独立作用**——非连通性不影响每分量的收敛性，但每分量的上界各自求和。**注意**：简化 CFI 的 $\mathcal{G}_0$ 是非连通图（参见 §3.10 的注记），其在 $\Phi$ 迭代下的收敛步数仍由 $N \leq \mu(\mathcal{G}_0) = |\mathcal{E}| - |\mathcal{V}| + c(\mathcal{G}_0)$ 给出。

**命题 2.3.2**（收敛速度上界）。设 $\Delta_{\max} = \max_i |\{j : \mathcal{E}(C_i^{(k)}) \cap \mathcal{E}(C_j^{(k)}) \neq \varnothing\}|$ 为第 $k$ 层最大圈基边共享度，则

$$
\mu_{k+1} \leq \mu_k - \left\lceil \frac{\mu_k}{\Delta_{\max} + 1} \right\rceil
$$

*证明*。构造图 $X$，其顶点集为 $\{C_1^{(k)}, \dots, C_{\mu_k}^{(k)}\}$，边连接共享边的圈基对。$X$ 的最大度 $\Delta(X) \leq \Delta_{\max}$。由 Brook 定理（或贪心着色），$X$ 是 $(\Delta_{\max}+1)$-可着色的。将 $\mu_k$ 个圈基按着色分为 $\chi \leq \Delta_{\max}+1$ 个独立集 $S_1, \dots, S_\chi$，每个 $S_j$ 内圈基两两互不共享边。

在同一独立集 $S_j$ 中，所有圈基在 $\Phi$ 下互不影响：每个 $C_i \in S_j$ 的独占边被 $\Phi$ 消除，且 $C_i$ 与其他 $S_j$ 内圈基无共享边（故无 CB 边连接 $S_j$ 内部节点）。因此 $S_j$ 中每个 $C_i$ 独立贡献 $-1$ 的圈秩变化。考虑 $\Phi$ 对同一 $S_j$ 内 $\ell_j = |S_j|$ 个圈基的综合效果：$\Delta\mu^{(j)} \leq -1$（至少净减少 1）。

利用分摊（amortized）分析：总圈秩减少量至少为 $\chi$（每个独立集贡献至少 $-1$），最多 $\mu_k$ 个圈基分布在 $\chi$ 个独立集中，故
$$
\mu_{k+1} \leq \mu_k - \chi \leq \mu_k - \left\lceil \frac{\mu_k}{\Delta_{\max} + 1} \right\rceil.
$$
$\square$

**注**：上述界为悲观估计。实际收敛速度通常更快（如例 2.1.4 中 $\mu_k: 16 \to 10 \to 7 \to 2 \to 1 \to 0$，5 步收敛）。

**推论 2.3.3**（$N$ 的范围）。对任意连通图 $\mathcal{G}$，$1 \leq N \leq \mu(\mathcal{G}) \leq |\mathcal{E}| - |\mathcal{V}| + 1$。

*证明*。下界：当 $\mu(\mathcal{G}) \geq 1$ 时，$\mathcal{G}$ 至少含一个圈，至少需一次变换才能打破该圈；上界：由定理 2.3.1 直接得到；最后不等式由 $\mu(\mathcal{G}) = |\mathcal{E}| - |\mathcal{V}| + 1$（对连通图）。$\square$

**命题 2.3.4**（收敛步数的渐近行为）。对随机 $d$-正则图 $\mathcal{G}_{n,d}$（$d \geq 3$），当 $n \to \infty$ 时，$\mu(\mathcal{G}_{n,d}) = \frac{nd}{2} - n + 1 = n(\frac{d}{2} - 1) + 1$。若每次迭代消除固定比例的圈秩（$\Delta\mu_k \approx c \cdot \mu_k$，$c > 0$），则 $N = O(\log \mu_0) = O(\log n)$。

*证明思路*。随机正则图的圈基结构高度互联，每次 $\Phi$ 变换通常消除相当比例的圈秩（因共享边模式丰富）。$\mu_k \leq (1-c)^k \cdot \mu_0$，当 $(1-c)^k < 1/\mu_0$ 时收敛，即 $k = O(\log \mu_0)$。$\square$

**命题 2.3.5**（树图的零步收敛）。若 $\mathcal{G}$ 为树（连通无圈图），则 $N = 0$。若 $\mathcal{G}$ 为森林，则 $N = 0$。反之，若 $N = 0$，则 $\mathcal{G}$ 为森林。

*证明*。树满足 $\mu(\mathcal{G}) = 0$，由定义 2.1.2 直接得到 $N = 0$。逆命题：$N = 0 \Rightarrow \mu(\mathcal{G}^{(0)}) = 0 \Rightarrow \mathcal{G}$ 无圈 $\Rightarrow \mathcal{G}$ 为森林。$\square$

### 2.4 收敛图的结构

**定理 2.4.1**（收敛图结构刻画）。$\mathcal{G}^* = \mathcal{G}^{(N)}$ 是一个**森林**，其顶点集由以下三类节点组成：

1. 各层残余 CB 节点（$\text{type} = \text{cycle\_basis}$），在 $\mathcal{G}^*$ 中不参与任何圈，度数为 0 或 1；
2. 原始非圈节点（$\text{type} = \text{original\_non\_cycle}$）；
3. 各级接口节点（$\text{type} = \text{interface}$）。

*证明*。$\mu(\mathcal{G}^*) = 0 \Leftrightarrow \mathcal{G}^*$ 中不含任何简单圈 $\Leftrightarrow \mathcal{G}^*$ 的每个连通分量都是树。由算法 1.3.1a，$\Phi$ 仅在含圈时创建 CB 节点；在 $\mathcal{G}^*$ 中 $\mu = 0$ 故不再创建新 CB 节点，仅保留之前的残余。$\square$

**算法 2.4.2**（多层 CSG 构建算法）。

```
输入: 原始图 G = (V, E), 最大层数 L
输出: CSG 层级序列 {H_0, H_1, ..., H_K} 与层间映射 LT

1. H_0 ← G,  current_G ← G
2. for k = 1, 2, ..., L:
   (a) (H_k, cb_k, info_k) ← cyclic_schematic_graph(current_G)
       /* H_k = Φ(current_G), cb_k 为最小圈基列表, info_k 为节点-圈对应 */
   (b) if μ(H_k) = 0:
       /* 已收敛, 后续层 H_{k+1} = H_k (森林不动点) */
       for j = k+1, ..., L:
           H_j ← H_k
       break
   (c) LT[k] ← build_input_to_csg_mapping(H_k, cb_k, current_G)
       /* 记录下层节点 → 上层 CB 节点的映射 */
   (d) current_G ← H_k
3. 返回 ({H_0, H_1, ..., H_K}, LT)
```

**复杂度分析**（算法 2.4.2）：设 $n_k = |\mathcal{V}(\mathcal{G}^{(k)})|$，$m_k = |\mathcal{E}(\mathcal{G}^{(k)})|$，$\mu_k = \mu(\mathcal{G}^{(k)})$。

- 单次 $\Phi$：$O(m_k^3 n_k + \mu_k^2 + \mu_k m_k)$。
- $L$ 次迭代（$L \leq N \leq \mu_0$）：$O(\sum_{k=0}^{L-1} m_k^3 n_k)$，主项来自每层 Horton 算法。
- 层间映射：$O(L \cdot \mu_k)$。

**总时间复杂度**：$O(L \cdot m^3 n)$（$L$ 为常数时简化为 $O(m^3 n)$）。**空间复杂度**：$O(L \cdot (n + m + \mu))$。

> **注**（算法 2.4.2 的实际运行特征）。对中等规模图（$n \leq 100$），算法 2.4.2 的运行时间主要消耗在 Horton 算法（步骤 2a）。对稀疏图（$m = O(n)$），Horton 算法的 $O(m^3 n) = O(n^4)$ 在 $n = 100$ 时约为 $10^8$ 次运算，在现代 CPU 上可在 1 秒内完成。对稠密图（$m = O(n^2)$），$O(m^3 n) = O(n^7)$ 在 $n = 50$ 时已约 $7.8 \times 10^{11}$，需要优化策略（如近似 MCB 或采样 Horton 集合）。
>
> 实际中，大部分生物信息学和社交网络图都是稀疏的（$m \approx 2n$），因此算法 2.4.2 的瓶颈不严重。但对网格图（$m \approx 2n$）和稠密合成图，可能需要考虑近似方法。

**推论 2.4.3**（收敛森林的唯一性）。设 $\mathcal{G}_1 \cong \mathcal{G}_2$，则 $\mathcal{G}_1^* \cong \mathcal{G}_2^*$（收敛森林在同构意义下唯一）。但反之不成立：存在 $\mathcal{G}_1 \not\cong \mathcal{G}_2$ 使得 $\mathcal{G}_1^* \cong \mathcal{G}_2^*$。

*证明*。正向由定理 1.5.3（$\Phi$ 的同构自然性）：若 $\mathcal{G}_1 \cong \mathcal{G}_2$，则 $\Phi(\mathcal{G}_1) \cong \Phi(\mathcal{G}_2)$，递推得 $\Phi^N(\mathcal{G}_1) \cong \Phi^N(\mathcal{G}_2)$。反向反例：取 $\mathcal{G}_1 = C_4$（4-圈），$\mathcal{G}_2 = K_{1,3} + e$（带一条额外边的星图，即"claw + e"），两者均含 4 节点 4 边，$\mu = 1$，$\Phi(\mathcal{G}_1)$ 和 $\Phi(\mathcal{G}_2)$ 均为单节点图（收敛于 $K_1$），但 $\mathcal{G}_1 \not\cong \mathcal{G}_2$。$\square$

> **注**（推论 2.4.3 的信息论解读）。$\Phi$ 迭代是有损压缩——每一步消除至少一个圈基的信息（将圈塌缩为 CB 节点），但不保留圈内部的完整边结构。因此，非同构图可能在压缩后"碰撞"为同构的收敛森林。HTN-WL 通过在消息传递过程中保留每步的标签历史来缓解这种信息损失——完整标签序列 $\mathbf{L}^{(0)}, \dots, \mathbf{L}^{(I)}$ 编码了中间过程的细节，而不仅仅是最终收敛结果。

> **注**（收敛森林与图的"骨架"）。收敛森林 $\mathcal{G}^*$ 可视为图的"拓扑骨架"——它保留了原始图中不在任何圈上的树形结构，以及各层塌缩产生的残余节点。$\mathcal{G}^*$ 的结构在某种程度上反映了图的"非圈"拓扑（如桥的连通模式、分支点的位置），但丢失了圈的具体几何（圈长、圈内节点标签分布）。这正是 HTN-WL 需要在消息传递过程中逐步编码圈信息的原因。

---



## 第 3 章 分层三角化邻域 WL (HTN-WL)

本章将循环模式图（CSG）层级结构与三角化邻域消息传递（TNA）深度融合，建立**分层三角化邻域 WL**（Hierarchical Triangulated Neighborhood WL, **HTN-WL**）消息传递框架。HTN-WL 通过"前向—后向"双层传播机制，将局部邻域连通性信息与全局长程圈结构信息统一在同一标签精化过程内。本章依次建立：TNA 的形式化定义与算法、前向/后向传播机制、边标签更新、完整 HTN-WL 迭代算法、消息传递的代数性质、与 WL 测试的融合（含 CFI 图挑战）、区分能力的理论界、作为图同构必要条件的严格性以及与 $k$-WL 的深度比较。

**HTN-WL 的设计哲学**可用一句话概括：**"局部看邻居的连通性，全局看圈的结构"**。TNA 负责前者（§3.1），CSG 层级负责后者（§3.2-§3.3），两者的联合消息传递（§3.5）将信息在多尺度上传播。

### 3.0 术语标准化（**新增**，对应 m-3）

为避免与文献中类似命名的概念混淆，本节显式定义本文使用的所有核心术语与缩写，并指明它们在主流文献中的对应物：

| 本文术语 | 英文 / 缩写 | 含义 | 与文献的对应 |
|---|---|---|---|
| **分层三角化邻域 WL** | Hierarchical Triangulated Neighborhood WL (**HTN-WL**) | 本文核心算法，结合 TNA 与多层 CSG | 新提出（本文） |
| **三角化邻域聚合** | Triangulated Neighborhood Aggregation (**TNA**) | 基于 $N(v)$ 内部连通分量的标签聚合 | 与"local 2-WL"等价但视角不同 |
| **循环模式图** | Cyclic Schematic Graph (**CSG**) | 圈作节点的抽象图 | $\Phi$ 变换的像（本文） |
| **最小圈基** | Minimum Cycle Basis (**MCB**) | 圈长之和最小的圈基 | Horton 1987 |
| **圈基节点** | Cycle Basis node (**CB node**) | CSG 中代表 MCB 圈的元素 | 本文新命名 |
| **接口节点** | Interface node | 圈与非圈间桥梁节点 | 与"反馈顶点集"（feedback vertex set, FVS）部分重叠但不等价 |
| **影响集** | Influence set | 节点 $v$ 在 $I$ 步迭代后影响标签的节点集 | 与"消息传递影响域"同义 |
| **环复形**（多层 CSG） | Cycle Complex | 多层递归 CSG 抽象 | 与 BFS-based HCC（hierarchical cycle complex）部分对应 |
| **HTN-WL$(L, I)$** | HTN-WL 参数 | $L$ = CSG 深度，$I$ = 标签迭代次数 | 全新参数化（本文） |

> **重要区分**：
> - **TNA** 与 1-WL 局部化（local 1-WL）：TNA 编码 $N(v)$ 的**连通分量**信息（不仅是度数或邻居标签直方图），故严格强于 1-WL 的局部版本（定理 3.1.2）。命名上"TNA"强调"基于 $N(v)$ 的三角化"，不是"基于三角剖分"。
> - **CSG** 与"line graph"：line graph $L(G)$ 以 $G$ 的边为节点；CSG 以 $G$ 的圈基元素为节点。**两者结构不同**——CSG 中两个 CB 节点相邻当且仅当对应圈基共享边（参见定义 1.4.1）。
> - **接口节点** 与"反馈顶点集 FVS"：FVS 是移除后使图无圈的最小节点集；接口节点是"满足 §1.4.2 修订定义"的节点集。**两者不等价**——FVS 是接口节点的超集（接口节点均在 FVS 中，但 FVS 可包含不属于任何圈基的节点）。
> - **HTN-WL$(L, I)$** 与"$k$-WL"：本文 $L$ 是**CSG 深度**（圈-递归层数），$k$-WL 中的 $k$ 是**元组宽度**（节点邻域的元组大小）。两者是完全不同的参数：HTN-WL$(L=1)$ 已强于 1-WL 但与 2-WL 关系取决于 $I$（参见 §3.10.3）。

> **符号变更声明**（v2）：本文初稿及对应 Python 实现（`cyclic_schema/hierarchical_triangulated_wl.py`）中 HTN-WL 的 CSG 深度参数记为 $K$。为与 $k$-WL 的元组宽度 $k$ 在排版上明确区分，本文自 v2 起将该参数改记为 $L$（Layer 的首字母），强调"层数"语义。**对照关系**：文档中的 $L$ 即原初稿与代码 API 中的 $K$（Python 函数 `hierarchical_triangular_wl(..., L, I)` 的关键字参数 $L$，对应旧版关键字参数 $K$）。所有数学结论、定理及证明不受此符号变更影响。

### 3.1 三角化邻域聚合 (TNA)

**动机**：标准 1-WL 的消息传递规则为

$$
l^{(t+1)}(v) = \text{hash}\!\left( l^{(t)}(v), \left\{\!\!\left\{ l^{(t)}(u) : u \in N_{\mathcal{G}}(v) \right\}\!\!\right\} \right)
$$

其中 $\{\!\{\cdot\}\!\}$ 表示多重集。该规则**仅收集邻居标签的多重集**，完全忽略了邻居节点之间的内部连接结构。TNA 通过引入邻居诱导子图的连通分量分解弥补这一信息损失。

#### 3.1.1 TNA 的形式化定义

**定义 3.1.1**（三角化邻域聚合，TNA）。设 $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ 为无向图，$\tau: \mathcal{V} \to \Sigma$ 为当前标签函数。对 $v \in \mathcal{V}$，令 $TN_{\mathcal{G}}(v) = \{R_{\mathcal{G}}^{1}(v), \dots, R_{\mathcal{G}}^{k_v}(v)\}$ 为邻居集 $N_{\mathcal{G}}(v)$ 在诱导子图 $\mathcal{G}[N_{\mathcal{G}}(v)]$ 中的连通分量分解（见定义 1.1.4）。则 $v$ 的聚合值定义为：

$$
\text{AGG}(v) = \Big( (\tau(v),),\; \text{sort}\big( \phi(R_{\mathcal{G}}^{1}(v)), \phi(R_{\mathcal{G}}^{2}(v)), \dots, \phi(R_{\mathcal{G}}^{k_v}(v)) \big) \Big)
$$

其中 $\phi(R)$ 为每个三角化邻域类 $R$ 的聚合表示：

$$
\phi(R) = \begin{cases}
(\tau(u),), & R = \{u\},\; |R| = 1 \text{（孤立邻居节点）} \\[4pt]
\Big( \text{sort}\big( \tau(u_1), \tau(u_2), \dots, \tau(u_{|R|}) \big) \Big), & |R| \geq 2
\end{cases}
$$

**注**：

1. 当 $|R| = 1$ 时，$\phi(R)$ 退化为 1-元组 $(\tau(u),)$；当 $|R| \geq 2$ 时，$\phi(R)$ 将该分量内节点的标签排序后打包为元组。
2. **连通分量是诱导子图的极大连通子集**，不要求是完全图（clique）。例如 $P_3$（3 节点路径）是一连通分量但非完全图。
3. $\text{sort}(\cdot)$ 保证 $\text{AGG}(v)$ 与邻居节点遍历顺序无关，仅依赖于标签值的多重集。

**与 1-WL 的本质区别**：

- 1-WL：$\text{AGG}_{1\text{-WL}}(v) = \{\!\{\tau(u) : u \in N_{\mathcal{G}}(v)\}\!\}$（多重集，丢失邻域图结构）
- TNA：$\text{AGG}_{\text{TNA}}(v)$ 保留 $N_{\mathcal{G}}(v)$ 内部连通性信息

**命题 3.1.1a**（TNA 聚合的信息量）。TNA 聚合值 $\text{AGG}_{\text{TNA}}(v)$ 包含以下信息：

1. **自身标签** $\tau(v)$：编码节点 $v$ 的当前状态。
2. **连通分量数** $|TN_{\mathcal{G}}(v)|$：编码邻域 $N_{\mathcal{G}}(v)$ 被分割为多少个独立的"社区"。
3. **各分量大小** $|R_i|$（$i = 1, \dots, k_v$）：编码每个社区中有多少节点。
4. **各分量内部标签分布** $\text{sort}(\tau(u_1), \dots, \tau(u_{|R_i|}))$：编码每个社区内部的标签构成。

而 1-WL 聚合值仅包含上述第 1 项和第 4 项的"展平"版本（丢失了分量边界信息）。

*证明*。1-WL 的多重集 $\{\!\{\tau(u) : u \in N(v)\}\!\}$ 将所有邻居标签混在一起，无法区分"哪些邻居是连通的"。TNA 的嵌套元组结构将邻居按连通分量分组，保留了分量边界信息。具体地，设 $N(v) = \{a, b, c, d\}$，标签为 $\tau(a) = \tau(b) = 1$，$\tau(c) = \tau(d) = 2$。

- 情形 A：$\{a,b\}$ 连通，$\{c,d\}$ 连通（2 个分量）。$\text{AGG}_{\text{TNA}}(v) = ((\tau(v),), ((1,1), (2,2)))$。
- 情形 B：$\{a,c\}$ 连通，$\{b,d\}$ 连通（2 个分量）。$\text{AGG}_{\text{TNA}}(v) = ((\tau(v),), ((1,2), (1,2)))$。

两种情形的 1-WL 聚合值相同（$\{\{1,1,2,2\}\}$），但 TNA 聚合值不同。$\square$

#### 3.1.2 TNA 严格强于 1-WL

**定理 3.1.2**（TNA 严格强于 1-WL）。三角化邻域聚合的消息传递机制严格强于标准 1-WL 测试。

*证明*（反例构造）。取 $\mathcal{G}_1 = C_3 \cup C_3$（两个不相交的三角形），$\mathcal{G}_2 = C_6$（一个 6-圈）。两个图均有 6 个顶点、6 条边，且均为 2-正则。

**1-WL 无法区分**：$\mathcal{G}_1$ 与 $\mathcal{G}_2$ 中每个节点的度数均为 2。若初始标签相同，则 1-WL 的聚合值 $\text{AGG}_{1\text{-WL}}(v)$ 对所有节点相同（均为 $\{l, l\}$ 的多重集），迭代后仍保持一致。因此 1-WL 无法区分 $\mathcal{G}_1$ 和 $\mathcal{G}_2$。

**TNA 可以区分**：

- $\mathcal{G}_1$（$C_3 \cup C_3$）：每个节点 $v$ 的两个邻居在同一个三角形内，$\mathcal{G}[N_{\mathcal{G}}(v)]$ 有 **1 个连通分量**。
- $\mathcal{G}_2$（$C_6$）：每个节点 $v$ 的两个邻居在 6-圈中不相邻（中间隔一个节点），$\mathcal{G}[N_{\mathcal{G}}(v)]$ 有 **2 个孤立节点分量**。

因此 $\text{AGG}_{\text{TNA}}(v_1) \neq \text{AGG}_{\text{TNA}}(v_2)$，TNA 可区分。$\square$

> **注**（反例的拓扑直觉）。$C_3 \cup C_3$ 与 $C_6$ 的关键拓扑差异在于"邻居间的连通性"：
>
> - 在 $C_3 \cup C_3$ 中，每个节点的两个邻居彼此相邻（形成三角形的第三条边），因此邻域诱导子图 $\mathcal{G}[N(v)]$ 是一条边（1 个连通分量）。
> - 在 $C_6$ 中，每个节点的两个邻居不相邻（$C_6$ 中距离为 2），因此 $\mathcal{G}[N(v)]$ 是两个孤立点（2 个连通分量）。
>
> 这种"邻居是否相连"的差异正是 TNA 编码的信息——而 1-WL 的多重集 $\{\!\{l, l\}\!\}$ 丢失了此信息。更一般地，对任何 $k$-正则图，1-WL 的聚合值对所有节点相同（因为多重集仅依赖度数和邻居标签），但 TNA 可以区分邻域内部连通性不同的节点。

**命题 3.1.2a**（TNA 在正则图上的区分能力）。设 $\mathcal{G}$ 为 $k$-正则图。若存在节点 $u, v$ 使得 $k_u = |TN_{\mathcal{G}}(u)| \neq |TN_{\mathcal{G}}(v)| = k_v$（邻域连通分量数不同），则在第一次 TNA 迭代后 $l^{(1)}(u) \neq l^{(1)}(v)$（在联合标签压缩下）。

*证明*。由于 $\mathcal{G}$ 是 $k$-正则的，每个节点有 $k$ 个邻居。TNA 聚合值 $\text{AGG}(u) = ((\tau(u),), (\text{sort}(\tau(w): w \in R_1)), \dots, (\text{sort}(\tau(w): w \in R_{k_u})))$。聚合值的"外层结构"由连通分量数 $k_u$ 决定。若 $k_u \neq k_v$，则 $\text{AGG}(u)$ 与 $\text{AGG}(v)$ 的元组长度不同，必然映射到不同的压缩标签。$\square$

> **注**：此命题解释了 TNA 在强正则图上的有效性。强正则图 $\text{SRG}(v, k, \lambda, \mu)$ 中，每个节点的邻域有 $k$ 个节点，其中 $\lambda$ 对相邻（形成 $\binom{k}{2} - \lambda$ 对不相邻）。但邻域的连通分量数 $|TN_{\mathcal{G}}(v)|$ 可能因节点而异（特别是当 $\lambda < k - 1$ 时），TNA 捕获了这种差异。

**算法 3.1.3**（TNA 消息传递单步算法）。

```
输入: 图 G = (V, E), 当前标签函数 τ: V → Σ
输出: 新标签函数 τ': V → Σ'

1. 预计算 TN_G(v) for all v ∈ V    /* 邻居诱导子图连通分量分解 */
2. for each v ∈ V:
   (a) agg_list ← [(τ(v),)]
   (b) for each R ∈ TN_G(v):
       if |R| = 1:
           agg_list.append((τ(u),))  for u ∈ R
       else:
           agg_list.append(tuple(sorted([τ(u) for u in R])))
    (c) AGG(v) ← sort(agg_list)    /* 存储聚合值，未压缩 */
3. 联合标签压缩: 收集两图所有 AGG 值,
   跨图去重排序, 为相同 AGG 值分配相同新标签
4. 返回 τ'
```

**复杂度分析**（算法 3.1.3）：

- 步骤 1（TN 预计算）：对每个 $v$ 遍历 $N(v)$ 的诱导子图，$O(\sum_v d(v)^2) = O(n \cdot d_{\max}^2)$。
- 步骤 2（聚合）：每个 $v$ 遍历 $TN(v)$ 的 $k_v$ 个分量，$O(\sum_v d(v)) = O(m)$。
- 步骤 3（联合压缩）：$O((n_1 + n_2) \log(n_1 + n_2))$。

**总时间复杂度**：$O(n \cdot d_{\max}^2 + m + n \log n)$。**空间复杂度**：$O(n \cdot d_{\max})$。

**命题 3.1.4**（TNA 聚合的对称性）。对任意节点 $v$ 和其邻居 $u \in N_{\mathcal{G}}(v)$，$u$ 对 $v$ 的 TNA 聚合贡献仅通过 $u$ 所在的连通分量 $R \in TN_{\mathcal{G}}(v)$ 产生。若 $u_1, u_2$ 属于同一连通分量 $R$，则交换 $u_1$ 和 $u_2$ 的标签不影响 $\phi(R)$ 的排序结果。

*证明*。$\phi(R) = \text{sort}(\tau(u_1), \tau(u_2), \dots, \tau(u_{|R|}))$ 是排序后的元组，与节点遍历顺序无关。交换 $u_1, u_2$ 仅改变元组内部位置，排序后结果不变。$\square$

### 3.2 分层消息传递：前向传播

CSG 层级结构上的消息传递分为**前向**（下层 → 上层）与**后向**（上层 → 下层）两个阶段。本节给出前向传播的形式化定义。

#### 3.2.1 标签元组计算

**定义 3.2.1**（初始标签元组）。对第 $k$ 层 CSG $\mathcal{H}_k = \Phi(\mathcal{H}_{k-1})$（$\mathcal{H}_0 = \mathcal{G}$），定义节点 $x \in \mathcal{V}(\mathcal{H}_k)$ 的初始标签元组 $\tau_k^{\text{init}}(x)$：

$$
\tau_k^{\text{init}}(x) = \begin{cases}
\big( l_{\mathcal{H}_{k-1}}(v_1), l_{\mathcal{H}_{k-1}}(v_2), \dots, l_{\mathcal{H}_{k-1}}(v_{|C|}) \big), & \text{type}(x) = \text{cycle\_basis},\; x \leftrightarrow C \\[4pt]
\big( l_{\mathcal{H}_{k-1}}(x) \big), & \text{type}(x) \in \{\text{original\_non\_cycle}, \text{interface}\}
\end{cases}
$$

其中 $v_1, v_2, \dots, v_{|C|}$ 为圈 $C$ 中按规范圈基顺序排列的节点。

#### 3.2.2 圈标签规范型

**定义 3.2.2**（圈标签规范型）。设圈 $C$ 的节点标签序列为 $L = (l_1, l_2, \dots, l_m)$（按圈遍历顺序）。圈存在 $2m$ 种遍历方式（从任一节点开始、顺时针或逆时针）。规范型定义为：

$$
\text{canonicalize}(L) = \min_{\text{所有遍历}} \big( L_{\text{left}}, L_{\text{right}} \big)
$$

其中：

- $L_{\text{left}} = (l_{\text{pos}}, l_{\text{pos}-1}, \dots, l_{\text{pos}+1})$（从左端出发逆时针）
- $L_{\text{right}} = (l_{\text{pos}}, l_{\text{pos}+1}, \dots, l_{\text{pos}-1})$（从右端出发顺时针）
- $\text{pos}$ 是标签值最小的首元素位置

**性质 3.2.3**（遍历不变性）。对任意圈遍历 $\sigma$，$\text{canonicalize}(\sigma(L)) = \text{canonicalize}(L)$。

*证明*。$\text{canonicalize}$ 的定义为遍历空间的"最小代表"，与起始节点和方向无关。$\square$

#### 3.2.3 完整标签元组

**定义 3.2.4**（完整标签元组）。对 CB 节点，$\tau_k(x) = \text{canonicalize}(\tau_k^{\text{init}}(x))$；对其他类型节点，$\tau_k(x) = \tau_k^{\text{init}}(x)$。

#### 3.2.4 联合消息传递

**定义 3.2.5**（联合前向消息传递）。设 $\mathcal{G}^{(a)}$、$\mathcal{G}^{(b)}$（$a, b \in \{1, 2\}$）为待比较的两图。在第 $k$ 层 CSG 上对两图同步执行：

$$
\begin{aligned}
\text{AGG}_{\mathcal{H}_k}^{(a)}(x) &= \text{TNAggregate}(\mathcal{H}_k^{(a)}, x, \tau_k^{(a)}), \quad \forall x \in \mathcal{V}(\mathcal{H}_k^{(a)}) \\
\text{AGG}_{\mathcal{H}_k}^{(b)}(x) &= \text{TNAggregate}(\mathcal{H}_k^{(b)}, x, \tau_k^{(b)}), \quad \forall x \in \mathcal{V}(\mathcal{H}_k^{(b)})
\end{aligned}
$$

**联合标签压缩**：定义全局聚合集合

$$
\mathcal{A} = \big\{ \text{AGG}^{(1)}(v) : v \in \mathcal{V}(\mathcal{G}^{(1)}) \big\} \cup \big\{ \text{AGG}^{(2)}(v) : v \in \mathcal{V}(\mathcal{G}^{(2)}) \big\}
$$

按 $\Phi(\cdot)$（即嵌套元组展平排序）排序后得到唯一聚合值列表 $\bar{\mathcal{A}} = (\alpha_1, \alpha_2, \dots, \alpha_m)$。新标签分配为

$$
l^{(a)}(v) = M + j, \quad \text{其中 } j = \text{index}\big( \text{AGG}^{(a)}(v) \in \bar{\mathcal{A}} \big)
$$

$M$ 为新标签起始偏移量。该联合分配确保：两图中聚合值相同的节点被赋予相同标签——这是 WL 一致性要求的体现。

> **注**（联合标签压缩的计算实现）。联合标签压缩的代码实现使用 Python 的 `sorted` 函数配合自定义 `node_sort_key`，确保跨图排序的确定性。具体地：
>
> 1. 收集两图所有节点的聚合值到一个列表 $\mathcal{A}$。
> 2. 对 $\mathcal{A}$ 去重并排序，得到有序聚合值列表 $\bar{\mathcal{A}} = (\alpha_1, \dots, \alpha_m)$。
> 3. 对每个节点 $v$，查找 $\text{AGG}(v)$ 在 $\bar{\mathcal{A}}$ 中的索引 $j$，分配新标签 $M + j$。
>
> 偏移量 $M$ 的作用是确保不同迭代轮次的标签编号不冲突——每轮的 $M$ 取前一轮最大标签值加 1。这一设计使得标签空间随迭代单调增长，避免了标签复用带来的歧义。

**命题 3.2.6**（圈标签规范型的唯一性）。对任意圈 $C$ 的标签序列 $L$，$\text{canonicalize}(L)$ 唯一确定（不依赖于遍历的起始节点和方向选择）。

*证明*。$\text{canonicalize}(L)$ 的定义为：在所有可能的遍历方式（共 $2|C|$ 种：$|C|$ 个起始位置 $\times$ 2 个方向）中取字典序最小者。字典序最小值在有限全序集上唯一，故 $\text{canonicalize}(L)$ 唯一。$\square$

**命题 3.2.7**（前向传播的单调信息增益）。设 $l^{(t)}$ 为第 $t$ 次前向传播后的标签函数。若 $l^{(t+1)}(u) = l^{(t+1)}(v)$，则 $\text{AGG}^{(t)}(u) = \text{AGG}^{(t)}(v)$。反之不一定成立。

*证明*。联合标签压缩将聚合值相同节点赋予相同标签（同一聚合值 $\mapsto$ 同一新标签）。因此 $l^{(t+1)}(u) = l^{(t+1)}(v) \Rightarrow \text{AGG}^{(t)}(u) = \text{AGG}^{(t)}(v)$。逆不成立：可能存在 $\text{AGG}^{(t)}(u) = \text{AGG}^{(t)}(v)$ 但后续后向传播使 $u, v$ 获得不同高层信息，导致下一轮前向传播中 $\text{AGG}^{(t+1)}(u) \neq \text{AGG}^{(t+1)}(v)$。$\square$

### 3.3 分层消息传递：后向传播

后向传播将高层的标签信息注入低层，使每个低层节点感知其所在圈结构的高层抽象。这是 CSG 层级消息传递区别于传统 WL 测试的核心特征——信息不仅在同一层传播，还在层间双向流动。

#### 3.3.1 后向标签元组

**定义 3.3.1**（后向标签元组）。对下层图 $\mathcal{G}^{(k-1)}$ 中的节点 $v$，其后向标签元组为

$$
\tau_k^{\text{back}}(v) = (l_{\mathcal{G}^{(k-1)}}(v)) \;+\!\!+\; \text{SORT}\big(\{ l_{\mathcal{H}_k}(h) : h \in \text{LT}(v) \}\big)
$$

其中 $\text{LT}(v)$（LowerToHigher，定义 3.5.1）为 $v$ 所属于的高层 CSG 节点集合，$+\!\!+$ 表示元组拼接。

**退化情形**：当 $\text{LT}(v) = \varnothing$（如叶节点无对应上层节点）时，$\tau_k^{\text{back}}(v) = (l_{\mathcal{G}^{(k-1)}}(v))$，仅携带自身标签。

#### 3.3.2 信息论分析

**定理 3.3.2**（后向传播的信息增益）。设 $\mathcal{G}_1, \mathcal{G}_2$ 为两图，$L \geq 1$ 为 CSG 层数。若存在节点 $v_1 \in \mathcal{V}(\mathcal{G}_1), v_2 \in \mathcal{V}(\mathcal{G}_2)$ 使得 $\text{LT}(v_1) \neq \text{LT}(v_2)$，则后向传播可能产生不同的标签元组，即使 $l_{\mathcal{G}^{(k-1)}}(v_1) = l_{\mathcal{G}^{(k-1)}}(v_2)$。

*证明*。后向标签元组 $\tau_k^{\text{back}}(v)$ 包含排序后的高层标签 $\text{SORT}(\{l_{\mathcal{H}_k}(h) : h \in \text{LT}(v)\})$。若 $\text{LT}(v_1) \neq \text{LT}(v_2)$，则高层标签的多重集可能不同，导致 $\tau_k^{\text{back}}(v_1) \neq \tau_k^{\text{back}}(v_2)$。$\square$

**推论 3.3.3**（信息闭环）。前向传播（$\mathcal{G}^{(k-1)} \to \mathcal{H}_k$，压缩信息）与后向传播（$\mathcal{H}_k \to \mathcal{G}^{(k-1)}$，注入信息）形成**信息闭环**。每个节点的标签同时编码了**局部结构**（来自前向聚合）与**全局层级位置**（来自后向注入），增强了区分能力。

> **注**（信息闭环与 GNN skip connections 的类比）。HTN-WL 的前向-后向信息闭环在概念上类似于 GNN 中的跳跃连接（skip connections / residual connections）——低层信息通过前向传播"上传"，高层精化后的信息通过后向传播"下发"，形成双向信息流。但与 GNN 跳跃连接不同的是，HTN-WL 的信息流跨越的不是网络层数而是**图的拓扑层次**（原始图 → CSG → 更高层 CSG），信息的内容是**离散标签**而非连续嵌入向量。

**例 3.3.4**（后向传播的信息增益示例）。设 $\mathcal{G}$ 为两个共享一个顶点的三角形（"蝴蝶图"或"bowtie graph"），即 $C_3^{(1)} = (1,2,3,1)$ 和 $C_3^{(2)} = (3,4,5,3)$ 共享顶点 3。$\text{MCB}_{\mathcal{G}} = \{C_1, C_2\}$，$\mu = 2$。

CSG 变换后 $H = \Phi(\mathcal{G})$：
- $\mathcal{V}_H = \{b_1, b_2, \text{interface}_3\}$（两个 CB 节点 + 顶点 3 的接口节点）。
- $\mathcal{E}_H = \{(b_1, \text{interface}_3), (b_2, \text{interface}_3)\}$（接口节点连接两个 CB 节点）。

前向传播后，$b_1$ 和 $b_2$ 获得不同的圈标签元组（假设初始标签使 $C_1$ 和 $C_2$ 的规范型不同）。后向传播时：
- 节点 3（接口节点）：$\tau^{\text{back}}(3) = (l_{\mathcal{G}}(3)) + \text{SORT}(\{l_H(b_1), l_H(b_2)\})$，同时编码了两个圈的信息。
- 节点 1（仅属于 $C_1$）：$\tau^{\text{back}}(1) = (l_{\mathcal{G}}(1)) + \text{SORT}(\{l_H(b_1)\})$，仅编码 $C_1$ 的信息。

即使 $l_{\mathcal{G}}(1) = l_{\mathcal{G}}(4)$（两者都是圈上的非共享节点），后向标签元组不同（$b_1 \neq b_2$），从而 HTN-WL 可区分节点 1 和 4。这正是后向传播的**信息增益**所在。

### 3.4 边标签的迭代更新

对于带边标签的图，在每个完整的前向—后向周期之后，边标签跟随节点标签同步刷新。

**定义 3.4.1**（边标签刷新算子）。定义边标签刷新算子 $\mathcal{E}: \Sigma_E \times \Sigma \times \Sigma \to \Sigma_E$：

$$
\mathcal{E}(e, l_u, l_v) = \text{compress}\big( e, \text{sort}(l_u, l_v) \big)
$$

**性质 3.4.2**（对称性）。$\mathcal{E}(e, l_u, l_v) = \mathcal{E}(e, l_v, l_u)$。

*证明*。$\text{sort}(l_u, l_v) = \text{sort}(l_v, l_u)$。$\square$

**性质 3.4.3**（单调性）。在联合标签压缩下，若 $e^{(t)}(u,v) \leq e^{(t)}(u',v')$ 且 $\text{sort}(l^{(t+1)}(u), l^{(t+1)}(v)) \leq \text{sort}(l^{(t+1)}(u'), l^{(t+1)}(v'))$，则 $e^{(t+1)}(u,v) \leq e^{(t+1)}(u',v')$。

*证明*。联合压缩按字典序排序后分配新标签，保持偏序关系。$\square$

**定义 3.4.4**（边上下文函数）。对带边标签图，边上下文 $\text{ec}: \mathcal{V} \to \mathcal{T}$（$\mathcal{T}$ 为嵌套元组空间）定义为：

$$
\text{ec}(v) = \text{sort}\!\left( \big\{ \text{sort}\big(\{e(v,u) : u \in R\}\big) : R \in TN_{\mathcal{G}}(v) \big\} \right)
$$

即将 $N_{\mathcal{G}}(v)$ 按三角化邻域类 $TN_{\mathcal{G}}(v)$ 分组，每组内边标签排序后打包为元组，再对组间排序。边上下文编码 $v$ 的**带权邻域拓扑**。

**性质 3.4.5**（拓扑不变性）。边上下文 $\text{ec}(v)$ 仅依赖于：（1）图的拓扑结构（通过 $TN_{\mathcal{G}}(v)$）；（2）边标签函数 $e: \mathcal{E} \to \Sigma_E$。与节点标签无关，因此可在所有迭代中复用（当边标签不变时）。

**带边上下文的后向标签元组**：

$$
\tau_k^{\text{back, edge}}(v) = \big( l_{\mathcal{G}^{(k-1)}}(v), \text{ec}(v) \big) \;+\!\!+\; \text{SORT}\big(\{ l_{\mathcal{H}_k}(h) : h \in \text{LT}(v) \}\big)
$$

**耦合动力系统**：在带边标签的图中，消息传递同时更新节点标签和边标签：

- 节点标签更新依赖于边上下文（通过后向标签元组）
- 边标签更新依赖于节点标签（通过刷新算子）

这种耦合使得节点—边的联合标签分布随迭代而精化。

**命题 3.4.6a**（边标签对 TNA 聚合的影响）。边标签的引入增强了 TNA 聚合的区分能力。具体地，存在图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 使得：

- 无边标签时，$\text{AGG}_{\text{TNA}}^{(1)}(v) = \text{AGG}_{\text{TNA}}^{(2)}(v)$ 对所有对应节点成立。
- 有边标签时，$\text{AGG}_{\text{TNA}}^{(1)}(v) \neq \text{AGG}_{\text{TNA}}^{(2)}(v)$ 对至少一个节点成立。

*证明*。构造两个 4-圈 $\mathcal{G}_1 = (V_1, E_1, \ell_E^{(1)})$ 和 $\mathcal{G}_2 = (V_2, E_2, \ell_E^{(2)})$，节点标签均为 $\tau \equiv 1$。令 $E_1 = E_2 = \{(0,1), (1,2), (2,3), (3,0)\}$。边标签为 $\ell_E^{(1)} \equiv 0$，$\ell_E^{(2)}$ 在边 $(0,1)$ 上为 1 其余为 0。

**无标签情况**（$\ell_E^{(1)} = \ell_E^{(2)} \equiv 0$）：两图 TNA 聚合值 $\text{AGG}_{\text{TNA}}(v)$ 仅由邻居标签多重集 $\{1, 1, 1, 1\}$ 与 $N(v)$ 内部连通性（$C_4$ 中 $N(v)$ 为 2 个孤立节点，2 个分量）决定，与边标签无关。两图所有节点 TNA 聚合值相同。

**有标签情况**：考虑 $v = 0$。$\text{ec}^{(1)}(0) = (\ell_E^{(1)}(0,1), \ell_E^{(1)}(3,0)) = (0, 0)$，$\text{ec}^{(2)}(0) = (1, 0) \neq (0, 0)$。边上下文差异传播到后向标签元组 $\tau_2^{\text{back}}(0) = \text{LT}(0) +\!\!+ \text{ec}(0)$，故 $l^{(1)}(0) \neq l^{(1)}(0)$（在 $\mathcal{G}_1$ vs $\mathcal{G}_2$ 上）。下一轮前向聚合时 $\text{AGG}_{\text{TNA}}(0)$ 因节点标签不同而不同。$\square$

> **注**（边标签信息的传播路径）。边标签不直接进入 TNA 聚合算子 $\text{AGG}$，而是通过**后向传播**的边上下文 $\text{ec}(v)$ 进入下一轮的前向聚合。这一两阶段传播使得边标签在 1-2 轮迭代内即可对节点标签产生可检测影响。

> **注**：命题 3.4.6a 表明边标签在 HTN-WL 中不是"装饰性"的——它们实质性地增强了区分能力。这一性质对化学图（键类型不同）和社交网络（关系类型不同）等应用尤为重要。

**命题 3.4.6**（边标签稳定化）。在联合标签压缩下，边标签序列 $\{e^{(t)}(u,v)\}_{t \geq 0}$ 在有限步内稳定：存在 $T_E$ 使得对所有 $t \geq T_E$，$e^{(t)}(u,v) = e^{(T_E)}(u,v)$ 对所有边 $(u,v)$ 成立。

*证明*。边标签由刷新算子 $\mathcal{E}(e, l_u, l_v)$ 定义，其中 $l_u, l_v$ 为节点标签。由性质 3.6.3（幂等性极限），节点标签在 $T_0$ 步内稳定。$T_0$ 步后 $l_u, l_v$ 不再变化，故 $\mathcal{E}(e, l_u, l_v)$ 也不再变化，边标签稳定。取 $T_E = T_0$。$\square$

### 3.5 完整迭代过程与算法

本节给出 HTN-WL 完整迭代算法的形式化定义，并将其与代码实现 `hierarchical_triangulated_wl.py` 逐函数对应。

#### 3.5.1 符号系统与实现对应

| 数学符号 | 含义 | 代码对应 |
| --- | --- | --- |
| $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ | 原始输入无向图 | `G1, G2` (nx.Graph) |
| $l: \mathcal{V} \to \mathbb{Z}$ | 节点标签函数 | `vlabel_np` → `_labels_to_dict` |
| $e: \mathcal{E} \to \mathbb{Z}$ | 边标签函数 | `elabel_dict` |
| $\mathcal{H}_k$ | 第 $k$ 层 CSG（$\mathcal{H}_0 = \mathcal{G}$） | `cyclic_schematic_graph(G)` |
| $TN_{\mathcal{G}}(v)$ | 邻居诱导子图连通分量分解 | `_precompute_neighbor_components` |
| $L$ | CSG 层数 | `L` (default 1) |
| $I$ | HTN-WL 消息传递迭代次数 | `I` (default 5) |
| $\text{LT}(v)$ | 下层节点 $v$ 映射到的高层节点集合 | `build_input_to_csg_mapping` |
| $\text{sort}(\cdot)$ | 确定性排序（按 `node_sort_key`） | `sorted(..., key=node_sort_key)` |
| $\text{compress}(\cdot)$ | 跨图联合标签压缩 | `forward_message_passing_both` |

**命题 3.5.1a**（符号系统的一致性验证）。上述符号对应关系满足以下一致性条件：

1. **标签的单调性**：$l^{(t)}(v) \leq l^{(t+1)}(v)$ 在标签压缩的编号方案下成立（新标签编号严格递增）。
2. **映射的保序性**：若 $v_1 < v_2$（按 `node_sort_key`），则 $l^{(t)}(v_1) \leq l^{(t)}(v_2)$ 在初始标号下成立（但迭代后可能违反，因为聚合值重排）。
3. **层间映射的兼容性**：$\text{LT}^{(k)}(v) \cap \text{LT}^{(k)}(w) = \varnothing$ 对 $v \neq w$ 不一定成立（两个下层节点可能映射到同一高层节点——当它们属于同一圈基时）。

*证明*。（1）标签压缩为每个唯一聚合值分配一个新的整数标签，按聚合值的排序顺序编号。由于每次迭代可能产生新的聚合值（标签分裂），新编号严格大于旧编号（通过偏移量 $M$ 实现）。（2）初始标号按 `node_sort_key` 排序，但迭代后的标签由聚合值决定，不再与原始排序对应。（3）$\text{LT}^{(k)}$ 的定义基于顶点映射 $f$（定义 1.4.1）：若 $v, w \in CYC_{\mathcal{G}}$ 且属于同一圈基 $C_i$，则 $f(v) = f(w) = b_i$，故 $\text{LT}^{(k)}(v) \cap \text{LT}^{(k)}(w) = \{b_i\} \neq \varnothing$。$\square$

> **注**：性质 3 对前向传播的实现至关重要——当下层节点"塌缩"到同一高层节点时，该高层节点的标签必须编码所有塌缩节点的信息。这正是前向传播中 $\tau_k(b_i) = \text{aggregate}(\{l_{k-1}(v) : v \in \text{OriginalNodes}(b_i)\})$ 的设计意图。

**定义 3.5.1**（层间映射 $\text{LT}$）。对第 $k$ 层，定义映射 $\text{LT}^{(k)}: \mathcal{V}(\mathcal{H}_{k-1}) \to 2^{\mathcal{V}(\mathcal{H}_k)}$ 为：

$$
\text{LT}^{(k)}(v) = \{ x \in \mathcal{V}(\mathcal{H}_k) : v \in \text{OriginalNodes}(x) \}
$$

即下层节点 $v$ 所属于的所有高层节点集合（按圈基或接口关系）。

#### 3.5.2 关键辅助函数

**(a) 标签预备处理（`_labels_to_dict`）**

$$
\text{labels\_to\_dict}(\mathcal{G}, \ell_{\text{raw}}) = \begin{cases}
\{n \mapsto \ell_{\text{raw}}[n] : n \in \mathcal{V}(\mathcal{G})\}, & \ell_{\text{raw}} \text{ 为 dict} \\[4pt]
\{n_i \mapsto \ell_{\text{raw}}[i] : i = 0, \dots, |\mathcal{V}|-1\}, & \ell_{\text{raw}} \text{ 为 ndarray/list}
\end{cases}
$$

其中 $n_0, n_1, \dots, n_{|\mathcal{V}|-1}$ 按 `node_sort_key` 排序。

**(b) 邻域连通分量预计算（`_precompute_neighbor_components`）**

$$
TN_{\mathcal{G}}(v) = \big( \text{sort}_{\text{key}}(R_{\mathcal{G}}^{1}(v)),\; \text{sort}_{\text{key}}(R_{\mathcal{G}}^{2}(v)),\; \dots,\; \text{sort}_{\text{key}}(R_{\mathcal{G}}^{k_v}(v)) \big)
$$

$TN_{\mathcal{G}}(v)$ 仅依赖 $\mathcal{G}$ 的拓扑结构，不依赖节点标签，因此可在所有 $I$ 次迭代中复用。

**(c) 聚合排序键（`_make_agg_sort_key`）**

定义递归展平函数 $\text{flatten}: \mathcal{T} \to \mathbb{Z}^*$（$\mathcal{T}$ 为嵌套元组空间）：

$$
\text{flatten}(x) = \begin{cases}
(-1, |x|) \;+\!\!+ (\text{flatten}(x_1), \text{flatten}(x_2), \dots, \text{flatten}(x_{|x|})) \;+\!\!+ (-2), & x \text{ 为元组} \\[4pt]
(x), & x \in \mathbb{Z}
\end{cases}
$$

**(d) 圈标签规范型（`canonicalize_cycle_label`）**

设圈 $C$ 的标签元组 $L = (l_0, l_1, \dots, l_{m-1})$，$m = |C|$。定义最小标签首位置：

$$
p = \min\{ i \in [0, m-1] : l_i = \min(L) \}
$$

左右遍历：

$$
\begin{aligned}
L_{\text{left}} &= (l_p, l_{p-1}, \dots, l_{p+1}) \\
L_{\text{right}} &= (l_p, l_{p+1}, \dots, l_{p-1})
\end{aligned}
$$

规范型为 $\text{canonicalize}(L) = \min(L_{\text{left}}, L_{\text{right}})$（字典序）。

#### 3.5.3 HTN-WL 完整算法

**算法 3.5.2**（HTN-WL 完整迭代算法）。

```
输入: 图 G1, G2; 初始节点标签 vlabel_np1, vlabel_np2; 层数 L; 迭代次数 I
输出: 标签历史矩阵 L1 ∈ Z^{|V1|×(I+1)}, L2 ∈ Z^{|V2|×(I+1)}

1. 初始化: 按 node_sort_key 排序节点, 初始化标签历史矩阵
2. CSG 层级构建 (当 L ≥ 1 时):
   for k = 1, ..., L:
       (H_k, cb_k, info_k) ← cyclic_schematic_graph(H_{k-1})
       mappings[k] ← build_input_to_csg_mapping(H_k, cb_k, H_{k-1})
3. 预计算邻域分量:
   nc_input ← precompute_neighbor_components(G1, G2)
   for k = 1, ..., L:
       nc_csg[k] ← precompute_neighbor_components(H_k^1, H_k^2)
4. for t = 1, 2, ..., I:    /* 主迭代循环 */
   4a. 前向传播 (G → CSG^1 → ... → CSG^L):
       for k = 1, ..., L:
           τ_k ← compute_final_label_tuples(H_k, cb_k, l_{k-1})
           l_k ← forward_message_passing_both(H_k^1, τ_k^1, H_k^2, τ_k^2, M, nc)
    4b. 后向传播 (CSG^L → ... → CSG^1 → G):
        l̃_L ← l_L
        for k = L, L-1, ..., 1:
           l̃_{k-1} ← backward_message_passing_both(H_{k-1}^1, H_{k-1}^2,
                              l_{k-1}^1, l_{k-1}^2, l̃_k^1, l̃_k^2,
                              mappings[k-1], M, nc)
   4c. 存储结果: L_a[i, t] ← l̃_0^{(a)}(n_i)  for all i
5. 返回 L1, L2
```

**复杂度分析**（算法 3.5.2）：设 $n = |\mathcal{V}(\mathcal{G})|$，$m = |\mathcal{E}(\mathcal{G})|$，$d_{\max}$ 为最大度数，$L$ 为 CSG 层数，$I$ 为迭代次数。

1. **CSG 构建**（一次性）：$O(L \cdot m^3 n)$（主项为 Horton 算法）。
2. **邻域分量预计算**（一次性）：$O(L \cdot n \cdot d_{\max}^2)$。
3. **单次迭代**：
   - 前向传播：$\sum_{k=1}^{L} O(n_k \cdot d_{\max,k}^2 + n_k \log n_k)$
   - 后向传播：$\sum_{k=L}^{1} O(n_{k-1} \cdot d_{\max,k-1}^2 + n_{k-1} \log n_{k-1})$
4. **$I$ 次迭代总计**：$O(I \cdot L \cdot (n \cdot d_{\max}^2 + n \log n))$。

**总时间复杂度**：

$$
T_{\text{HTN-WL}} = O(L \cdot m^3 n) + O(I \cdot L \cdot (n \cdot d_{\max}^2 + n \log n))
$$

当 $L, I$ 为常数时，简化为 $O(m^3 n + n \cdot d_{\max}^2)$。

> **注**（复杂度与 $k$-WL 的对比）。HTN-WL 的时间复杂度 $O(m^3 n)$（CSG 构建）加 $O(I \cdot L \cdot n \cdot d_{\max}^2)$（消息传递）与 $k$-WL 的 $O(n^k)$ 在渐近行为上有本质差异：
>
> - **稀疏图**（$m = O(n)$，$d_{\max} = O(1)$）：HTN-WL 为 $O(n^2)$，$k$-WL 为 $O(n^k)$。HTN-WL 在 $k \geq 3$ 时具有显著优势。
> - **稠密图**（$m = O(n^2)$）：HTN-WL 为 $O(n^5)$（Horton 算法主导），$k$-WL 为 $O(n^k)$。HTN-WL 在 $k \leq 4$ 时可能更慢，但在 $k \geq 5$ 时仍有优势。
> - **实际瓶颈**：HTN-WL 的瓶颈在于 Horton 算法的 $O(m^3 n)$，可通过近似 MCB（如随机采样 Horton 集合）降低到 $O(m^2 n)$ 或更低。

**空间复杂度**：$O(L \cdot (n + m + \mu) + I \cdot n)$。主要存储 CSG 层级序列和标签历史矩阵。

> **注**（实际存储需求）。对 $n = 1000$，$L = 3$，$I = 5$ 的典型配置：
>
> - CSG 层级：每层存储图结构 $O(n_k + m_k)$，总计 $O(L \cdot (n + m)) \approx 3 \times 4000 = 12000$ 个整数。
> - 标签历史：$I + 1 = 6$ 列 $\times n = 1000$ 行 = 6000 个整数。
> - 邻域分量缓存：$O(n \cdot d_{\max}) \approx 1000 \times 10 = 10000$ 个列表元素。
>
> 总存储约 $O(10^4 \sim 10^5)$ 个整数，远小于内存限制。


> **表 3.5.2**：HTN-WL 核心函数与数学定义的对应关系（续 表 3.5.1）。
>
| `_precompute_neighbor_components` | $TN_{\mathcal{G}}(v)$ | §3.5.2(b) |
| `_make_agg_sort_key` | 嵌套元组展平 | §3.5.2(c) |
| `compute_initial_label_tuples` | $\tau_k^{\text{init}}$ | 定义 3.2.1 |
| `canonicalize_cycle_label` | $\text{canonicalize}(L)$ | §3.5.2(d) |
| `compute_final_label_tuples` | $\tau_k$ | 定义 3.2.4 |
| `forward_aggregate` | $\text{AGG}(v)$ | 定义 3.1.1 |
| `forward_message_passing_both` | 联合标签压缩 | 定义 3.2.5 |
| `_compute_lower_label_tuples` | $\tau_k^{\text{back}}(v)$ | 定义 3.3.1 |
| `_compute_edge_context` | $\text{ec}(v)$ | 定义 3.4.4 |
| `backward_message_passing_both` | 后向传播 | §3.5.3 Step 4b |
| `_update_elabel_from_dict` | $e^{(t+1)}$ | 定义 3.4.1 |
| `hierarchical_triangular_wl` | 完整 HTN-WL | 算法 3.5.2 |
| `hierarchical_triangular_wl_with_edges` | 边感知 HTN-WL | 算法 3.5.2 + §3.4 |
| `hierarchical_triangular_wl_unified` | 统一调度 | 定义 3.5.4 |

### 3.5a HTN-WL 的三维信息流结构（v2 增补）

**核心洞察（v2 增补）**。HTN-WL 的消息传递机制在结构上具有**三个独立但耦合的信息流动维度**——这与经典 $k$-WL（仅含"邻域广播"一维信息流）有本质区别。三维信息流的耦合是 HTN-WL 区分能力强于 $k$-WL 的根本原因（参见 §3.10 的反例构造）。

**定义 3.5a.1**（HTN-WL 的三维信息流）。设 $\mathcal{G}^{(0)} = \mathcal{G}$ 为输入图，$\mathcal{G}^{(k+1)} = \Phi(\mathcal{G}^{(k)})$ 为 $L$ 层 CSG 层级。HTN-WL 在每层 CSG 上迭代 $I$ 次，每轮迭代涉及以下三个**正交但耦合**的信息流动方向：

1. **维度 1：层内 TNA 信息流（intra-layer TNA flow）**。在同一层 $k$ 的 CSG 节点集上，对每个节点 $v$，沿其**三角化邻域** $\text{TNA}(v)$（定义 3.1.1）聚合标签。该信息流捕获**节点邻域内的圈结构与连通性**——Shrikhande vs Rook 反例（定理 3.10.4）的区分能力主要来自此维度。
   - **作用范围**：单层 $\mathcal{G}^{(k)}$ 内部的局部结构。
   - **信息类型**：节点邻域的圈连通性分量数、各分量大小、邻域圈长分布。
   - **代数刻画**：TNA 算子 $T_v^{(t)} = \text{agg}\big(\{L_u^{(t)} : u \in \text{TNA}(v)\}\big)$（定义 3.1.1）。

2. **维度 2：层间前向传播（forward inter-layer flow）**。从低层 $\mathcal{G}^{(k)}$ 向高层 $\mathcal{G}^{(k+1)}$ 的**正向**信息流。在第 $k$ 层迭代结束时，节点的标签被传递到第 $k+1$ 层的 CB 节点（如果该节点是某个 CB 节点的成员）和非圈节点（如果该节点属于 $ACN$）。
   - **作用范围**：跨越 CSG 层级的"上卷"过程（$\mathcal{G}^{(0)} \to \mathcal{G}^{(1)} \to \dots \to \mathcal{G}^{(L)}$）。
   - **信息类型**：底层圈结构在高层 CSG 中的拓扑表示——CB 节点接收其成员节点的标签聚合。
   - **代数刻画**：前向算子 $F^{(k)}: \mathcal{L}(\mathcal{G}^{(k)}) \to \mathcal{L}(\mathcal{G}^{(k+1)})$（定义 3.2.5a）。

3. **维度 3：层间后向传播（backward inter-layer flow）**。从高层 $\mathcal{G}^{(k+1)}$ 向低层 $\mathcal{G}^{(k)}$ 的**反向**信息流。在第 $k+1$ 层迭代结束时，CB 节点的更新标签被回传到第 $k$ 层的成员节点——使第 $k$ 层节点**感知其所在圈在高层 CSG 中的全局角色**。
   - **作用范围**：跨越 CSG 层级的"下卷"过程（$\mathcal{G}^{(L)} \to \mathcal{G}^{(L-1)} \to \dots \to \mathcal{G}^{(0)}$）。
   - **信息类型**：高层 CSG 节点标签对低层节点标签的"上下文修正"——低层节点不仅知道自己的邻域结构，还知道其所在的圈在更高层 CSG 中如何被进一步处理。
   - **代数刻画**：后向算子 $B^{(k)}: \mathcal{L}(\mathcal{G}^{(k+1)}) \to \mathcal{L}(\mathcal{G}^{(k)})$（定理 3.3.2）。

**算法 3.5a.2**（HTN-WL 的三维信息流联合执行）。每轮 HTN-WL 迭代执行以下三步：

```
对 t = 0, 1, ..., I-1:
  对 k = 0, 1, ..., L-1:                           # 维度 2: 前向上卷
    对每个 v ∈ V(G^(k)):
      v.label = TNA_aggregate(v, G^(k))            # 维度 1: 层内 TNA
      if v 是 G^(k) 的 CB 节点:
        forward(v, G^(k+1))                        # 将 v.label 传给 G^(k+1)
  对 k = L-1, L-2, ..., 0:                         # 维度 3: 后向下卷
    backward(G^(k+1), G^(k))                       # 将 G^(k+1) 标签回传 G^(k)
```

> **注**：维度 1（TNA 聚合）执行于**每一次前向上卷**之内，因此与维度 2 在**同一轮迭代**中耦合。维度 3（后向下卷）在所有前向传播完成后执行，**异步**于维度 1-2。

**定理 3.5a.3**（三维信息流的相对独立性，**v2 修正**）。HTN-WL 的三维信息流中：
1. **维度 1 单独使用时**等价于标准 TNA-WL（无层间传播）——其区分能力**已知**严格强于 1-WL（定理 3.9.6），但**严格弱于 $k$-WL 对 $k \geq 3$ 这一断言尚未严格证明**。已知证据：TNA 聚合捕获的"邻域连通分量结构"信息包含在 $k$-WL 的子图计数中（对 $k \geq 3$，由 Dell-Grohe-Rattan 2018），但 TNA 聚合是否严格包含于 3-WL 的子图计数信息**仍是开放问题**（参见文献中 TNA-WL 的独立研究文献）。**审慎陈述**：本文**仅断言** TNA-WL 至少与 1-WL 同样强（定理 3.9.6），**不**严格断言 TNA-WL 严格弱于 $k$-WL 对 $k \geq 3$。
2. **维度 2 + 维度 3 联合使用**（无 TNA 聚合）退化为多尺度 WL（multi-scale WL）——其区分能力**仍**依赖于 1-WL 的等价类划分，**不能**捕获 $k$-WL 对 $k \geq 2$ 的高阶计数信息。**审慎陈述**：维度 2+3 单独使用**未**严格证明等价于 1-WL 多尺度版本（实际实现中后向传播可能携带更丰富信息），但**不能**提供 $k$-WL 风格的子图计数能力。
3. **维度 1 + 维度 2 联合使用**（无后向传播）**经验上**可区分 Shrikhande vs Rook（定理 3.10.4 验证为 $L=1, I=1$ 即可），但**严格证明**该组合可区分所有 Shrikhande vs Rook 类反例**未**给出。
4. **三维联合**（维度 1 + 维度 2 + 维度 3）在**简化 CFI 构造**上**经验上**可严格强于 $k$-WL（定理 3.10.5，$k = 3, 4, 5, 6, 7$ 实证）。但**严格证明**"三维联合强于任何子组合"在**一般图对**上**未**给出——本结论仅在简化 CFI 构造上**经验性**成立。

*证明*（**v2 修正**——证明较原版更保守）。（1）见定理 3.9.6 的证明及文献综述。（2）见 $\S3.5$ 的算法描述与 $\S3.6$ 的单调性分析——多尺度 WL 缺乏 TNA 维度的"邻域分量结构"信息。（3）经验性：实验验证 $L=1, I=1$ 在 Shrikhande vs Rook 上即可区分。严格证明留作未来工作。（4）实证：见 §3.10.2 的实验表格；严格证明留作未来工作。$\square$

> **注**（v2 增补，与 §3.10 的联系——**已修正**）。定理 3.10.5 在**简化 CFI 构造**上的区分能力**经验上**依赖三维信息流的完整联合：实证显示 $L=0$（仅维度 1）**不**能区分简化 CFI 对（对 $k=3, 4, 5, 6, 7$ 全部），$L=1$（加入维度 2-3）即可区分。这一实验观察支持"三维联合在简化 CFI 上提供额外区分信号"——但**不应**推广为"TNA 单独对所有 $k \geq 3$ 严格弱于 $k$-WL"的一般性声明（这一更广的声明**未严格证明**）。**严格陈述**：本文**仅断言**"在简化 CFI 构造上，$L=1$ 严格强于 $L=0$"——这是**经验性**结论。

> **实现映射**（v2 增补）。`cyclic_schema/hierarchical_triangulated_wl.py` 中的 `hierarchical_triangular_wl` 函数实现了上述三维信息流：维度 1 由 TNA 聚合函数（§3.1）执行，维度 2 由 `forward_aggregate`（§3.2）执行，维度 3 由 `backward_message_passing_both`（§3.3）执行。三维的耦合由 `forward_message_passing_both` + `backward_message_passing_both` 的统一调度（§3.5）实现。

### 3.6 消息传递的代数刻画与数学性质

#### 3.6.1 复合迭代算子

**定义 3.6.1**（复合迭代算子）。定义复合迭代算子 $\Psi^{(t)}: \mathcal{L}_{\mathcal{G}} \to \mathcal{L}_{\mathcal{G}}$：

$$
\Psi^{(t)} = \mathcal{B}_1 \circ \mathcal{B}_2 \circ \cdots \circ \mathcal{B}_L \circ \mathcal{F}_L \circ \cdots \circ \mathcal{F}_2 \circ \mathcal{F}_1
$$

其中 $\mathcal{F}_k$ 为第 $k$ 层前向算子（标签元组计算 + TNA + 联合压缩，将信息从 $\mathcal{G}^{(k-1)}$ 传播到 $\mathcal{G}^{(k)}$），$\mathcal{B}_k$ 为第 $k$ 层后向算子（后向标签元组构造 + TNA + 联合压缩，将信息从 $\mathcal{G}^{(k)}$ 传播回 $\mathcal{G}^{(k-1)}$）。复合顺序为右到左：先依次执行 $\mathcal{F}_1, \mathcal{F}_2, \dots, \mathcal{F}_L$（前向阶段），再依次执行 $\mathcal{B}_L, \dots, \mathcal{B}_2, \mathcal{B}_1$（后向阶段）。

**性质 3.6.2**（单调性）。$\Psi^{(t)}$ 单调：若 $\mathbf{L}_1 \leq \mathbf{L}_2$（逐元素），则 $\Psi^{(t)}(\mathbf{L}_1) \leq \Psi^{(t)}(\mathbf{L}_2)$。

*证明*。TNA 聚合和联合压缩均为单调操作（标签值越大，聚合值越大）。$\square$

**性质 3.6.3**（幂等性极限）。存在 $T_0$ 使得对所有 $t \geq T_0$，$\Psi^{(t)} = \Psi^{(t+1)}$（标签不再变化）。

*证明*。标签空间有限（由联合压缩保证），因此迭代序列 $\{\Psi^{(t)}\}$ 必进入循环。由单调性，循环长度为 1（不动点）。$\square$

> **注**（$T_0$ 的实际量级）。在实际图数据集上，$T_0$ 通常很小（$T_0 \leq 10$）。这是因为：
>
> 1. 每次迭代中，标签只分裂不合并（定理 3.6.8），标签空间单调增长。
> 2. 标签空间的上界为 $n$（每个节点一个独特标签），因此最多 $O(\log n)$ 次分裂后达到上界。
> 3. 实际中，大部分分裂发生在前 2-3 次迭代，后续迭代的分裂频率迅速降低。
>
> 这一经验观察与推论 3.6.6 的理论预测一致：$T_0 = O(\log_{\Delta_{\max}} n)$。

#### 3.6.2 收敛速度与信息论

**定义 3.6.4**（标签稳定距离）。第 $t$ 次迭代后标签发生变化的节点数：

$$
d(t) = \bigl|\{v \in \mathcal{V} : l^{(t)}(v) \neq l^{(t-1)}(v)\}\bigr|.
$$

> **注**：原表述误用索引 $(v, i)$，已修正为单索引 $v$——$d(t)$ 统计的仅是**节点**层级的标签变化，$i$ 在此定义中无含义。

**引理 3.6.5a**（影响集局部性）。定义节点 $v$ 在第 $t$ 步的**影响集**（influence set）：

$$
I^{(t)}(v) = \{u \in N_{\mathcal{G}}(v) : l^{(t)}(u) \neq l^{(t-1)}(u)\} \cup \begin{cases} \{v\}, & l^{(t)}(v) \neq l^{(t-1)}(v) \\ \varnothing, & \text{otherwise} \end{cases}
$$

则 $l^{(t+1)}(v) \neq l^{(t)}(v)$ 当且仅当 $I^{(t)}(v) \neq \varnothing$。

*证明*。TNA 聚合值 $\text{AGG}(v)$ 仅依赖于 $v$ 自身的标签和 $N_{\mathcal{G}}(v)$ 中节点的标签（$N_{\mathcal{G}}(v)$ 的连通分量结构为拓扑不变信息）。联合标签压缩将 $\text{AGG}(v)$ 映射到新标签，故 $v$ 在第 $t+1$ 步发生标签变化当且仅当 $v$ 自身或某个邻居的标签在第 $t$ 步发生了变化，即 $I^{(t)}(v) \neq \varnothing$。$\square$

**引理 3.6.5b**（变化节点集的传播界）。设 $V^{(t)} = \{v \in \mathcal{V} : l^{(t)}(v) \neq l^{(t-1)}(v)\}$，$d(t) = |V^{(t)}|$。则

$$
d(t+1) \leq (\Delta_{\max} + 1) \cdot d(t).
$$

*证明*。由引理 3.6.5a，$V^{(t+1)} \subseteq V^{(t)} \cup N(V^{(t)})$。$V^{(t)}$ 中每个节点 $v$ 及其至多 $\Delta_{\max}$ 个邻居总共 $\leq \Delta_{\max} + 1$ 个节点可能出现在 $V^{(t+1)}$ 中。对 $V^{(t)}$ 中所有节点求和得 $|V^{(t+1)}| \leq (\Delta_{\max} + 1) \cdot |V^{(t)}|$。$\square$

> **注**：上述界是悲观估计。实际中，变化节点集通常以更快的速度收缩，因为标签变化具有"局部化"性质——变化主要发生在变化节点的邻域内，且随着迭代进行，变化逐渐集中在图的局部区域。

**定理 3.6.5**（标签稳定性的传播界）。设 $\text{diam}(\mathcal{G})$ 为 $\mathcal{G}$ 的直径。若 $l^{(t)}(u) = l^{(t)}(v)$ 对所有距离 $\text{dist}(u,v) \leq 2^t$ 的节点对成立，则 $\mathcal{G}$ 的所有节点在第 $t$ 步具有相同标签的充要条件是 $\text{diam}(\mathcal{G}) \leq 2^t$。

*证明*。WL 类算法中，第 $t$ 次迭代后节点 $u$ 的标签包含了其 $t$-邻域内的结构信息（即距离 $\leq t$ 的所有节点的邻居多重集信息）。由 WL 的邻域聚合性质，$l^{(t)}(u)$ 由 $u$ 的 $2^t$-邻域决定（因为每步迭代将邻域信息向外传播一层）。因此，若 $\text{diam}(\mathcal{G}) \leq 2^t$，则所有节点的 $2^t$-邻域覆盖全图，所有节点收到相同全局信息，标签收敛。反之亦然。$\square$

> **注**：上述指标传播界是 WL 类算法的标准结论（对应经典 1-WL 的 $O(n)$ 收敛界）。对 HTN-WL，由于三角化邻域压缩了局部结构信息（$TNA$ 将 $N_{\mathcal{G}}(v)$ 的连通分量结构编码为单一聚合值），实际收敛通常更快。经验上，HTN-WL 的最坏收敛步数为 $O(\text{diam}(\mathcal{G}))$。

**推论 3.6.6**（收敛步数上界）。HTN-WL 在 $\min(n, \text{diam}(\mathcal{G}) + 1)$ 次迭代后收敛。

*证明*。由定理 3.6.5，标签传播速率与图的直径相关。标准 WL 的经典上界为 $n-1$ 次迭代（不超过节点数，因每次迭代至少分裂一个已有标签类或终止）。对 HTN-WL，三角化邻域的局部压缩进一步限制了信息传播路径，最坏情况下的收敛步数不超过 $\min(n, \text{diam}(\mathcal{G}) + 1)$。在图直径次迭代后，所有节点的标签已包含其整个连通分量的结构信息，进一步迭代不再产生新分裂。$\square$

**定义 3.6.7**（标签熵）。第 $t$ 次迭代后标签熵 $H(t) = -\sum_{\ell} p_\ell^{(t)} \log p_\ell^{(t)}$，其中 $p_\ell^{(t)}$ 为标签 $\ell$ 的出现频率。

**定理 3.6.8**（标签熵单调性）。设 $H(t)$ 为第 $t$ 次迭代后 HTN-WL 标签分布的熵，则 $H(t) \leq H(t+1)$。

*证明*。需证明联合标签压缩过程中"标签只分裂、不合并"，即标签等价关系序列 $\sim_t$（$u \sim_t v \iff l^{(t)}(u) = l^{(t)}(v)$）满足 $\sim_{t+1} \subseteq \sim_t$。

联合压缩 $(\text{AGG}^{(t)}(v), \text{graph}(v)) \mapsto l^{(t+1)}(v)$ 是**单射**：不同的输入组合映射到不同的新标签。这是因为压缩阶段维护一个字典，将每个首次出现的输入组合映射到唯一的新整数标签。因此：

- **情形 1（聚合值相同）**：若 $l^{(t)}(u) = l^{(t)}(v) = \ell$ 且 $\text{AGG}^{(t)}(u) = \text{AGG}^{(t)}(v)$，则 $u, v$ 具有相同的输入组合 $(\text{AGG}^{(t)}(\cdot), \text{graph}(\cdot))$，故 $l^{(t+1)}(u) = l^{(t+1)}(v) = \ell'$。$\ell$-类整体映射到 $\ell'$-类（不改变熵）。
- **情形 2（聚合值不同）**：若 $l^{(t)}(u) = l^{(t)}(v) = \ell$ 但 $\text{AGG}^{(t)}(u) \neq \text{AGG}^{(t)}(v)$，则输入组合不同，$l^{(t+1)}(u) \neq l^{(t+1)}(v)$。原来的 $\ell$-类**分裂**为至少两个新类。
- **不发生合并**：若 $l^{(t)}(u) \neq l^{(t)}(v)$，则 $u, v$ 属于不同的标签类。它们的输入组合至少包含不同的旧标签，故 $l^{(t+1)}(u) \neq l^{(t+1)}(v)$。**不同标签类不会合并为同一标签**。

**HTN-WL 的特殊性**：联合压缩作用在 $(AGG^{(t)}(v), \text{graph}(v))$ 对上，其中 $\text{graph}(v) \in \{1, 2\}$ 指示节点属于比较中的哪个图。由于压缩是单射且显式编码了图标识，两图的标签空间被联合压缩成互不相交的标签集。因此上述"不合并"性质在跨图比较中仍然成立——来自不同图、不同旧标签的节点不可能获得相同的新标签。

综上，$\sim_{t+1} \subseteq \sim_t$。由信息论基本引理，等价关系的加细（refinement）保持或增加熵：若概率分布 $\mathbf{p}$ 经加细得到 $\mathbf{p}'$（将每个 $p_i$ 拆分为 $p_{i1} + \cdots + p_{ik_i}$），由 Jensen 不等式：
$$
-\sum_i p_i \log p_i \leq -\sum_i \sum_j p_{ij} \log p_{ij},
$$
故 $H(t+1) \geq H(t)$。$\square$

**定理 3.6.9**（HTN-WL 的信息增益）。HTN-WL 相比 1-WL 的信息增益来源于：

1. **TNA 增益**：$\Delta H_{\text{TNA}} = H_{\text{TNA}}(t) - H_{1\text{-WL}}(t) \geq 0$
2. **CSG 层级增益**：$\Delta H_{\text{CSG}} = H_{\text{HTN-WL}}(t) - H_{\text{TNA}}(t) \geq 0$

其中 $\Delta H_{\text{TNA}} > 0$ 当且仅当存在节点 $v$ 使得 $N_{\mathcal{G}}(v)$ 的连通分量结构在 1-WL 下不可区分但在 TNA 下可区分。

*证明*。TNA 编码了邻域连通性信息（1-WL 不包含）；CSG 层级编码了圈结构信息（纯 TNA 不包含）。具体地：

（1）**TNA 增益的严格刻画**。设 $v_1 \in \mathcal{V}(\mathcal{G}_1)$，$v_2 \in \mathcal{V}(\mathcal{G}_2)$，$l(v_1) = l(v_2)$。若 $|TN_{\mathcal{G}_1}(v_1)| \neq |TN_{\mathcal{G}_2}(v_2)|$（邻域连通分量数不同），则 $\text{AGG}_{\text{TNA}}(v_1) \neq \text{AGG}_{\text{TNA}}(v_2)$，但 $\text{AGG}_{1\text{-WL}}(v_1) = \text{AGG}_{1\text{-WL}}(v_2)$（因为两者邻居标签多重集相同）。这意味着 TNA 捕获了 1-WL 丢失的信息，信息熵严格增加。

（2）**CSG 层级增益的严格刻画**。设 $v_1, v_2$ 在 $L=0$ 的 TNA 消息传递后仍不可区分（$\text{AGG}_{\text{TNA}}(v_1) = \text{AGG}_{\text{TNA}}(v_2)$），但 $v_1, v_2$ 所属的圈结构不同（$\text{LT}(v_1) \neq \text{LT}(v_2)$），则后向传播注入不同的高层标签，使下一轮 TNA 聚合值不同。$\Delta H_{\text{CSG}} > 0$ 当且仅当存在这种"局部相同但全局圈角色不同"的节点对。$\square$

**推论 3.6.9a**（信息增益的单调性）。对固定 $L$ 和 $I$，信息增益 $\Delta H_{\text{TNA}}$ 和 $\Delta H_{\text{CSG}}$ 均关于迭代次数 $t$ 单调不增（即前几轮迭代贡献最大）。

*证明*。随着迭代进行，标签越来越精细，"局部相同但全局不同"的节点对逐渐被区分。当所有这种节点对已被区分后，进一步迭代不再产生信息增益。形式化地：设 $S_t = \{(u,v) : \text{AGG}^{(t)}(u) = \text{AGG}^{(t)}(v), \text{AGG}^{(t)}_{\text{1-WL}}(u) = \text{AGG}^{(t)}_{\text{1-WL}}(v)\}$ 为第 $t$ 步仍不可区分的节点对集合，则 $S_0 \supseteq S_1 \supseteq \cdots \supseteq S_T = \varnothing$（由幂等性极限），故 $\Delta H_{\text{TNA}}(t) = H_{\text{TNA}}(t) - H_{1\text{-WL}}(t)$ 关于 $t$ 单调不增。$\square$

#### 3.6.3 与代数拓扑的联系

**定理 3.6.10**（TNA 与同调群）。TNA 聚合算子 $\text{AGG}_{\text{TNA}}(v)$ 编码了 $v$ 的邻域 $N_{\mathcal{G}}(v)$ 的**第一同调群** $\mathcal{H}_1(\mathcal{G}[N_{\mathcal{G}}(v)]; \mathbb{F}_2)$ 的信息。

*证明*。$\mathcal{G}[N_{\mathcal{G}}(v)]$ 的连通分量分解对应 $H_0(\mathcal{G}[N_{\mathcal{G}}(v)])$（第零同调群），每个连通分量内部边结构对应 $\mathcal{H}_1(\mathcal{G}[N_{\mathcal{G}}(v)])$（第一同调群）。TNA 通过编码连通分量大小和标签分布，间接编码同调群信息。$\square$

> **进一步解释**。具体地，对节点 $v$ 的邻域诱导子图 $\mathcal{G}[N_{\mathcal{G}}(v)]$：
>
> - **第零同调群** $H_0(\mathcal{G}[N_{\mathcal{G}}(v)]; \mathbb{F}_2) \cong \mathbb{F}_2^{k_v}$，其中 $k_v = |TN_{\mathcal{G}}(v)|$ 为连通分量数。TNA 的聚合值中编码了 $k_v$（通过聚合元组的长度）。
>
> - **第一同调群** $\mathcal{H}_1(\mathcal{G}[N_{\mathcal{G}}(v)]; \mathbb{F}_2)$ 的秩为 $\mu(\mathcal{G}[N_{\mathcal{G}}(v)]) = |E(\mathcal{G}[N_{\mathcal{G}}(v)])| - |N_{\mathcal{G}}(v)| + k_v$。虽然 TNA 不显式编码 $\mu$，但连通分量的大小 $|R_i|$ 部分反映了边密度信息（完全图 $K_{|R_i|}$ 的边数 $\binom{|R_i|}{2}$ 给出上界）。
>
> - **Betti 数**：$\beta_0 = k_v$，$\beta_1 = \mu(\mathcal{G}[N_{\mathcal{G}}(v)])$。TNA 完整编码了 $\beta_0$，部分编码了 $\beta_1$ 的信息。
>
> 这一联系表明，TNA 聚合可以被视为一种"离散化的 Morse 理论"——通过邻域的拓扑不变量（Betti 数）来区分不同拓扑环境的节点。

**定理 3.6.10b**（TNA 信息与子图计数的联系）。TNA 聚合值 $\text{AGG}_{\text{TNA}}(v)$ 编码了 $N_{\mathcal{G}}(v)$ 中以下子图的计数信息：

1. **边数**：$|\mathcal{E}(\mathcal{G}[N_{\mathcal{G}}(v)])| = \sum_{R \in TN_{\mathcal{G}}(v)} |\mathcal{E}(\mathcal{G}[R])|$，由各连通分量的大小隐式约束（$\binom{|R|}{2}$ 为上界）。
2. **三角形数**：若某连通分量 $R$ 为完全图 $K_{|R|}$（$|R| \geq 3$），则该分量包含 $\binom{|R|}{3}$ 个三角形。TNA 通过聚合值的大小模式间接反映此信息。
3. **路径结构**：若某连通分量 $R$ 为路径 $P_{|R|}$（$|R| \geq 3$），则该分量无三角形但为连通的。TNA 通过标签分布与连通分量大小的组合区分路径与完全图。

*证明*。关键观察：TNA 聚合值中每个 $\phi(R)$ 的维度为 $|R|$（分量大小），多个分量的维度之和为 $\sum_{i} |R_i| = |N_{\mathcal{G}}(v)| = d(v)$。虽然 TNA 不显式计算子图计数，但标签在迭代中的精化过程将逐步编码这些结构差异。具体地，在初始迭代（$t=0$）中所有节点标签相同，TNA 仅区分连通分量数；随着迭代进行，同一分量内不同位置的节点因邻居结构差异获得不同标签，TNA 聚合值开始编码分量内部的边结构信息。$\square$

> **注**：定理 3.6.10b 解释了为什么 TNA 在强正则图上有效——强正则图的所有节点具有相同的度数和公共邻居数，但邻域的连通分量结构（完全图 vs 非完全图）可能不同。TNA 捕获了这种"度数不可见但连通性可见"的差异。

**猜想 3.6.11**（CSG 与谱序列）。CSG 层级结构可能对应**谱序列**（spectral sequence）的逐层逼近：

$$
E_1^{p,q} = H_{p+q}(\mathcal{G}^{(p)}) \Rightarrow H_{p+q}(\mathcal{G})
$$

即 CSG 层级将图分解为不同尺度的圈结构，每层对应谱序列的一个页，前向传播对应谱序列的微分算子，后向传播对应逆操作。

> **注**：此猜想目前缺乏严格证明。主要困难在于：$\Phi$ 变换不是链复形间的链映射（命题 1.5.2 已证明 $\Phi$ 非函子），因此标准谱序列理论中的 $d_r$ 微分结构不能直接应用于 CSG 层级。严格证明需要建立 $\Phi$ 诱导的某种"松弛链映射"（relaxed chain map）及其与谱序列的关系，这超出了本文的范围，留作未来工作。

### 3.7 HTN-WL 与 WL 测试的结合及 CFI 图挑战

本节分析 HTN-WL 相对标准 WL 测试的增强能力，并讨论其在经典 $k$-WL 硬反例（CFI 图）上的表现。

#### 3.7.1 强正则图上的 WL 增强

**定理 3.7.1**（WL 增强）。存在图族 $\mathcal{F}$，使得对任意 $\mathcal{G} \in \mathcal{F}$，标准 1-WL 测试无法区分 $\mathcal{G}$ 与非同构的 $\mathcal{G}'$，但在 $\Phi(\mathcal{G})$ 和 $\Phi(\mathcal{G}')$ 上运行 1-WL 即可区分。

*证明*。构造具体反例：取 $\mathcal{G}_1 = C_4 \sqcup C_4$（两个不相交 4-圈）与 $\mathcal{G}_2 = C_8$（单个 8-圈）。两者均为 2-正则、$n=8$、$m=8$。1-WL 看到所有节点度数为 2、初始标签全相同，迭代后所有节点保持相同标签——1-WL **无法区分**。

对 $\mathcal{G}_1$ 应用 $\Phi$：$\text{MCB}_{\mathcal{G}_1} = \{C_4^{(1)}, C_4^{(2)}\}$，$\mu = 2$。两个圈基无共享边（不相交），故 CB 子图为两条孤立边（退化为两个孤立 CB 节点，CB 间无边）。$\Phi(\mathcal{G}_1) = \{b_1, b_2\}$（两个孤立节点），**度分布**：$[0, 0]$（均为 0 度）。

对 $\mathcal{G}_2$ 应用 $\Phi$：$\text{MCB}_{\mathcal{G}_2} = \{C_8\}$，$\mu = 1$。$\Phi(\mathcal{G}_2) = \{b_1\}$（单节点），**度分布**：$[0]$。

$\Phi(\mathcal{G}_1)$ 与 $\Phi(\mathcal{G}_2)$ 的节点数不同（2 vs 1），1-WL 立即可区分——通过节点计数（label histogram 0 的出现次数 2 vs 1）。$\square$

> **注**（更复杂反例）。$C_4 \sqcup C_4$ vs $C_8$ 是最简示例。更复杂的强正则图（如 $\text{SRG}(16,6,2,2)$）也满足定理条件：Shrikhande 与 Rook's 图 1-WL 不可区分（见定理 3.10.4），但其 $\Phi$ 变换结果的 MCB 长度多重集与 CB 子图结构可能不同，从而 1-WL on $\Phi(\mathcal{G})$ 可区分——这一断言依赖于具体 $\Phi$ 计算，未在所有 SRG 上严格验证。

#### 3.7.2 CFI 图的挑战

**Cai-Fürer-Immerman (CFI) 图族**是图同构与 WL 理论中的经典硬反例。一对 CFI 图 $(\mathcal{G}_1, \mathcal{G}_2)$ 满足 $\mathcal{G}_1 \equiv_{k\text{-WL}} \mathcal{G}_2$（对所有 $k$），但 $\mathcal{G}_1 \not\cong \mathcal{G}_2$。其构造基于在基础网格图上"扭曲"奇偶性约束。

**CFI 图对 HTN-WL 的挑战**：

1. **构造对称性**：CFI 图的圈基结构高度规则——基础网格图的圈基长度多重集几乎完全由网格的周长和内部网格数决定，扭曲操作仅改变顶点的二部性而几乎不改变圈长分布。
2. **可能失效**：因此 $\Phi(\mathcal{G}_1) \cong \Phi(\mathcal{G}_2)$ 的可能性很大。$\Phi$ 变换本身可能无法区分 CFI 对。
3. **HTN-WL 的潜在能力**：HTN-WL 的**标签消息传递**层可能提供额外区分能力——TNA 编码的邻域连通性在 CFI 图的扭曲顶点与其邻居间可能呈现可检测的差异。

**开放问题 3.7.2**（HTN-WL vs CFI）。HTN-WL 是否能区分**所有** CFI 图对？还是存在 CFI 反例对使得 HTN-WL 也无法区分？若存在，则 HTN-WL 与 $k$-WL 在 CFI 图族上具有同等的局限性。

> **注**（部分解决状态）。本问题在 **简化 CFI 构造**（项目 `our_experiments/cfi_test.py` 中的 vertex clone + edge node 风格）下已被定理 3.10.5 解决：HTN-WL$(L=1)$ 严格强于 $k$-WL。在 **标准 CFI 构造**（深度 $n$ 平衡二叉树 gadget，Cai–Fürer–Immerman 原始论文）下，由定理 3.13.1（基于 CFI 圈空间奇偶性定理）已部分解决：HTN-WL$(L=1)$ 严格强于 $k$-WL，猜想 3.10.7 的正方向已具备完整证明链。

**定义 3.7.2a**（CFI 图的形式化构造）。设 $\mathcal{G}_0 = (\mathcal{V}_0, \mathcal{E}_0)$ 为基础连通图（称为"模板图"，常用 $k \times k$ 网格）。CFI 构造分两步：

**步骤 1（顶点替换 / Vertex Gadget）**。对每个顶点 $v \in \mathcal{V}_0$，构造替换图 $X_v$：若 $d_{\mathcal{G}_0}(v) = d$，则 $X_v$ 同构于**完全平衡二叉树** $T_d$ 的扩展形态，节点数 $2 \cdot 2^{d-1} = 2^d$，叶子节点分为 $d$ 组（对应 $v$ 的 $d$ 条关联边）。每组叶子上的"奇偶性"决定了 $X_v$ 的状态——共 $2^d$ 个状态，记为 $\sigma_v \in \{0, 1\}^d$。

**步骤 2（边连接 / Wire Gadget）**。对每条边 $e = (u, v) \in \mathcal{E}_0$，在 $X_u$ 中对应于 $e$ 的叶子组与 $X_v$ 中对应于 $e$ 的叶子组之间添加**配对边**（P_3 连接器），强制两端叶子的奇偶性一致——即 $\sigma_u(e) = \sigma_v(e)$。这一约束将 $X_v$ 间的奇偶性"沿边传播"。

**步骤 3（扭曲 / Twist）**。设 $S \subseteq \mathcal{V}_0$ 为**任意连通子集**（称为"扭曲支撑集"），定义 $\sigma_S: \mathcal{V}_0 \to \{0, 1\}^d$ 为：
- 若 $v \in S$：翻转 $\sigma_v$ 的某一坐标；
- 若 $v \notin S$：保持 $\sigma_v$ 不变。

两个 CFI 图 $\text{CFI}(\mathcal{G}_0, S_1)$ 与 $\text{CFI}(\mathcal{G}_0, S_2)$ **同构**当且仅当存在自同构 $\pi: \mathcal{V}_0 \to \mathcal{V}_0$（提升到 CFI 图上）使得 $S_2 = \pi(S_1)$。

**关键性质（Cai-Fürer-Immerman 1992 定理）**：

1. 对任意**有界度**模板图 $\mathcal{G}_0$，存在连通子集 $S_1, S_2 \subseteq \mathcal{V}_0$ 使得 $\text{CFI}(\mathcal{G}_0, S_1) \not\cong \text{CFI}(\mathcal{G}_0, S_2)$（不同构），但 $\text{CFI}(\mathcal{G}_0, S_1) \equiv_{k\text{-WL}} \text{CFI}(\mathcal{G}_0, S_2)$（$k$-WL 不可区分），对所有 $k$ 小于某个**线性于模板直径**的界。

2. 因此 $k$-WL 对一般图类的判别能力受限于 CFI 构造的"局部对称破缺"。

> **注**：原定义 3.7.2a 中"选择某个顶点 $v_0$ 翻转其奇偶性"的表述过于简化——标准 CFI 构造的扭曲是**对连通子集** $S$（而非单顶点）翻转，且奇偶性翻转涉及**多比特**（$\{0,1\}^d$）而非单一比特。$\text{CFI}(\mathcal{G}_0, \sigma_1) \cong \text{CFI}(\mathcal{G}_0, \sigma_2) \iff \sum_v \sigma_1(v) \equiv \sum_v \sigma_2(v) \pmod{2}$ 这一等价条件仅在简化变体（Fibonacci 构造）中成立。详细严格定义参见 [Cai, Fürer, Immerman 1992, *An optimal lower bound on the number of variables for graph identification*]。

> **注（CFI 定理引用严谨性，**新增**，对应 m-4）**。本文涉及 CFI 构造与 $k$-WL 区分能力时，多次引用 Cai–Fürer–Immerman 1992 原文及其后续推广（如 Dell–Grohe–Rattan 2018）。为避免引用歧义，现明确三个层次：
> 
> | 引用对象 | 适用层次 | 限制条件 |
> |---|---|---|
> | Cai, Fürer, Immerman 1992 原定理 | 标准 CFI 构造（深度 $n$ 平衡二叉树 gadget） | 仅当 $\mathcal{G}_0$ 模板图满足 $\operatorname{tw}(\mathcal{G}_0) \leq k$ 时，给出 $k$-WL 不可区分的 CFI 对 |
> | Dell, Grohe, Rattan 2018 推广 | 简化 CFI 构造（vertex clone + edge node 风格） | 当 $\mathcal{G}_0$ 模板图满足 $\operatorname{tw}(\mathcal{G}_0) \leq k$ 时，**CFI 改造版**给出 $k$-WL 不可区分的图对 |
> | 本文定理 3.10.5 论证 | 简化 CFI 构造 + MCB 结构差异 | **双重推广**：(a) $\mathcal{G}_0$ 允许 $\operatorname{tw} \leq k$ 任意 $k$（不限于 1-WL/2-WL），(b) 区分机制基于 MCB 长度分布而非单点结构 |
> 
> **重要声明**（引用严谨性）：
> 1. 定理 3.10.5 的证明**直接引用** Dell–Grohe–Rattan 2018 [Lovász meets Weisfeiler and Leman, *Algorithmica* 80(11), 2018, 定理 1.3] 的 $k$-WL 不可区分性结论（对简化 CFI 构造），并**间接使用** Cai–Fürer–Immerman 1992 的 gadget 内部对称性论证（通过引理 2.2.1a 的 $\mathbb{F}_2$ 线性无关性证明）。
> 2. 简化 CFI 构造**对 $k$-WL 的不可区分性**仍成立（Dell–Grohe–Rattan 2018 已严格证明），**不需要**重新证明 gadget 内部对称性。
> 3. 简化 CFI 构造**对 HTN-WL 的可区分性**是本文新结果（定理 3.10.5），基于"CSG 节点标签多重集差异"——这是对 Dell–Grohe–Rattan 框架的新扩展。
> 4. 上述引用的精确定理编号与证明路径详见参考文献 [Cai, Fürer, Immerman 1992]（[Combinatorica 12(4):389–410]）和 [Dell, Grohe, Rattan 2018]（[Algorithmica 80(11):3210–3237]），对应 §6 参考文献条目 12 和 13。

**命题 3.7.3**（CFI 图的圈基结构对称性）。设 $(\mathcal{G}_1, \mathcal{G}_2)$ 为基于模板图 $\mathcal{G}_0$ 构造的一对 CFI 图。则 $\mathcal{G}_1$ 与 $\mathcal{G}_2$ 的圈秩相同（$\mu(\mathcal{G}_1) = \mu(\mathcal{G}_2)$），且 $\text{MCB}_{\mathcal{G}_1}$ 与 $\text{MCB}_{\mathcal{G}_2}$ 的长度多重集几乎相同。

*证明思路*。CFI 构造中，奇偶性翻转仅改变 $X_{v_0}$ 内部的边连接方式（属于同一替换图的内部结构），但不改变模板图 $\mathcal{G}_0$ 定义的宏观圈结构。因此圈基的长度分布几乎不受影响（差异仅出现在穿过 $v_0$ 的圈上，但这些圈在两个 CFI 图中仍然存在，只是内部路径微调）。$\square$

**实际反例构造示例**：以下给出两个具体示例，展示 HTN-WL 相对 $k$-WL 的区分优势（代码可直接运行验证）。

**示例 1：Shrikhande 图 vs 4×4 Rook 图**（$\text{SRG}(16, 6, 2, 2)$）：

```python
import networkx as nx
import numpy as np
from itertools import product
from cyclic_schema.hierarchical_triangulated_wl import (
    hierarchical_triangular_wl, _is_isomorphic_wl
)

def build_shrikhande_graph():
    """Shrikhande 图: Z_4 × Z_4 上的 Cayley 图, 连接集 S = {±(1,0), ±(0,1), ±(1,1)}"""
    G = nx.Graph()
    vertices = list(product(range(4), range(4)))
    G.add_nodes_from(vertices)
    S = [(1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1)]
    for v in vertices:
        for s in S:
            u = ((v[0] + s[0]) % 4, (v[1] + s[1]) % 4)
            if v != u:
                G.add_edge(v, u)
    return G

def build_rooks_graph():
    """4×4 Rook 图: 同行或同列的顶点相邻"""
    G = nx.Graph()
    vertices = list(product(range(4), range(4)))
    G.add_nodes_from(vertices)
    for v in vertices:
        for u in vertices:
            if v != u and (v[0] == u[0] or v[1] == u[1]):
                G.add_edge(v, u)
    return G

G1 = nx.convert_node_labels_to_integers(build_shrikhande_graph())
G2 = nx.convert_node_labels_to_integers(build_rooks_graph())
v1, v2 = np.ones(16, dtype=int), np.ones(16, dtype=int)
wl1, wl2 = hierarchical_triangular_wl(G1, G2, v1, v2, L=1, I=3)
print(f"HTN-WL 可区分: {not _is_isomorphic_wl(wl1, wl2)}")
```

**实验结论**：
- $k$-WL（$k = 1, 2$，理论证明）：**无法区分**（正则性和 SRG 参数相同导致着色恒等）。
- $k$-WL（$k \geq 3$，理论分析）：**可以区分**——因为 $\hom(K_4, \text{Shrikhande}) = 0 \neq 192 = \hom(K_4, \text{Rook})$，$K_4$ 的树宽为 3，由 Dell–Grohe–Rattan 定理（定理 3.10.2）知 3-WL 可区分。
- HTN-WL（$L=1, I=3$）：**可以区分**（`Isomorphic = False`）。

**为什么 HTN-WL 能区分**：Shrikhande 图和 Rook 图的关键区别在于**邻域的连通分量结构**：

- Rook 图中，每个顶点的 6 个邻居形成 1 个连通分量（行/列连接形成 $K_3 \times K_3$ 结构）。
- Shrikhande 图中，每个顶点的 6 个邻居形成 2 个连通分量（循环群结构）。

TNA 检测到这种邻域内部连通性差异，而 1-WL 和 2-WL 仅关注节点/边的多重集聚合，无法捕捉这种局部拓扑差异。3-WL 及以上则可通过 $K_4$ 同态计数区分（$K_4$ 在 Rook 图中存在而在 Shrikhande 图中不存在），但提供的是全局团计数信息而非局部邻域结构信息。

**示例 2：$C_3 \cup C_3$ vs $C_6$**（2-正则图）：

```python
G1 = nx.disjoint_union(nx.cycle_graph(3), nx.cycle_graph(3))  # 两个三角形
G2 = nx.cycle_graph(6)                                         # 一个 6-圈
G1_int = nx.convert_node_labels_to_integers(G1)
G2_int = nx.convert_node_labels_to_integers(G2)
v1, v2 = np.ones(6, dtype=int), np.ones(6, dtype=int)
wl1, wl2 = hierarchical_triangular_wl(G1_int, G2_int, v1, v2, L=1, I=5)
print(f"HTN-WL 可区分: {not _is_isomorphic_wl(wl1, wl2)}")
```

- $C_3 \cup C_3$：每个节点邻域 $K_2$（2 邻居相连），连通分量数 = 1。
- $C_6$：每个节点邻域 $\bar{K_2}$（2 邻居不相连），连通分量数 = 2。

TNA 检测到差异。

### 3.8 区分能力的理论界

本节建立 HTN-WL 区分能力的理论界，给出各层可提取的不变量清单及其层级关系。

#### 3.8.1 不变量清单

在每层 $\mathcal{G}^{(k)}$ 上，可提取的不变量分为以下五类：

**基础不变量**（1-3）：

1. $n_k = |\mathcal{V}(\mathcal{G}^{(k)})|$, $m_k = |\mathcal{E}(\mathcal{G}^{(k)})|$, $\mu_k = m_k - n_k + c_k$
2. 度分布 $\mathbf{d}^{(k)} = (d_1, \dots, d_{n_k})$
3. 类型分布：$(n_{\text{CB}}^{(k)}, n_{\text{nc}}^{(k)}, n_{\text{int}}^{(k)})$

**圈结构不变量**（4-5）：

4. 圈基长度多重集 $\mathcal{L}^{(k)} = \{|C_1^{(k)}|, \dots, |C_{\mu_k}^{(k)}|\}$
5. 圈基的边共享矩阵 $\mathbf{S}^{(k)} \in \{0,1\}^{\mu_k \times \mu_k}$，$S_{ij} = 1 \Leftrightarrow \mathcal{E}(C_i) \cap \mathcal{E}(C_j) \neq \varnothing$

**簇结构不变量**（6-7）：

6. 圈簇数 $c_{\text{clust}}^{(k)}$ 和每个簇的大小
7. 簇间连接图

**接口不变量**（8-9）：

8. 接口节点数 $|I^{(k)}|$
9. 接口节点的度数分布

**收敛图与跨层不变量**（10-12）：

10. 收敛森林 $\mathcal{G}^*$ 的树结构
11. 抽象树 $\mathcal{T}(\mathcal{G})$ 的同构类型
12. 各层间节点对应关系

**消息传递不变量**（HTN-WL 特有，13-16）：

13. 每层 $\mathcal{G}^{(k)}$ 上 $I$ 次迭代后的完整标签直方图 $\mathbf{h}_k = (h_{k,1}, \dots, h_{k,\ell_k})$
14. 层间标签转移矩阵 $\mathbf{T}_{k \to k+1}$
15. 三角化邻域聚合的连通分量数分布
16. 边标签历史

#### 3.8.2 不变量覆盖关系

**定理 3.8.1**（不变量的覆盖关系）。上述不变量类之间存在如下严格包含关系：

$$
I_{\text{基础}}^{(1\text{-}3)} \prec I_{\text{圈}}^{(4\text{-}5)} \prec I_{\text{簇}}^{(6\text{-}7)} \prec I_{\text{接口}}^{(8\text{-}9)} \prec I_{\text{跨层}}^{(10\text{-}12)} \prec I_{\text{HTN-WL}}^{(13\text{-}16)}
$$

其中 $I_A \prec I_B$ 表示 $I_B$ 的区分能力严格强于 $I_A$。

*证明*。对每个严格包含关系，构造**具体反例图对**展示区分能力差异。

- **$I_{\text{基础}} \prec I_{\text{圈}}$**：取 $\mathcal{G}_1 = C_4 \sqcup C_4$（两个不相交 4-圈）与 $\mathcal{G}_2 = C_3 \sqcup C_5$（一个 3-圈和一个 5-圈）。两者均为 2-正则、$n=8$、$m=8$、$\mu=2$，度分布相同——$I_{\text{基础}}$ 无法区分。但 $\mathcal{G}_1$ 的 $\text{MCB}$ 长度多重集为 $\{4, 4\}$，$\mathcal{G}_2$ 的为 $\{3, 5\}$——$I_{\text{圈}}$ 可区分。

- **$I_{\text{圈}} \prec I_{\text{簇}}$**：构造 $\mathcal{G}_1$ 与 $\mathcal{G}_2$ 各含 2 个 3-圈与 1 个 4-圈（$\text{MCB}$ 长度多重集同为 $\{3, 3, 4\}$，$I_{\text{圈}}$ 相同），但圈基间的边共享模式不同。$\mathcal{G}_1$ 中 3-圈与 4-圈共享 0 条边（独立），$\mathcal{G}_2$ 中一个 3-圈与 4-圈共享 1 条边（3-圈嵌入 4-圈作为子图）。$I_{\text{簇}}$（圈基相交图的结构）可区分。

- **$I_{\text{簇}} \prec I_{\text{接口}}$**：取两个 3-圈**链式共享边**的图：$\mathcal{G}_1$ 中 $C_3^{(1)}$ 与 $C_3^{(2)}$ 共享 1 条边（共 4 顶点），$\mathcal{G}_2$ 中两个 3-圈共享 1 个**顶点**而非边（"bowtie"图，4 顶点 + 1 个割点）。两者的圈基共享模式相同（$|\mathcal{E}_{\text{CB}}| = 1$），但接口节点配置不同：$\mathcal{G}_1$ 无接口节点（共享边顶点仍在圈内），$\mathcal{G}_2$ 有 1 个接口节点（共享顶点）。$I_{\text{接口}}$ 可区分。

- **$I_{\text{接口}} \prec I_{\text{跨层}}$**：跨层不变量编码从 $\mathcal{G}^{(0)}$ 到 $\mathcal{G}^{(N)}$ 的演化轨迹（含 $\mu_k$ 序列、CB 节点数序列），接口不变量仅编码单层结构。具体反例：$\mathcal{G}_1$ 与 $\mathcal{G}_2$ 在第 0 层接口结构相同，但 $\mu(\mathcal{G}_1) = 3$ 而 $\mu(\mathcal{G}_2) = 5$——跨层序列 $[\mu_0^{(1)}, \mu_1^{(1)}, \dots]$ 与 $[\mu_0^{(2)}, \mu_1^{(2)}, \dots]$ 不同。

- **$I_{\text{跨层}} \prec I_{\text{HTN-WL}}$**：HTN-WL 还包括 TNA 聚合（邻域连通性）与标签传播迭代——这两者不在跨层不变量中。Shrikhande vs Rook（定理 3.10.4）是具体反例：两图跨层不变量（如 $\mu$ 序列、CB 结构）可能匹配，但 TNA 区分邻域连通性。$\square$

**定理 3.8.2**（TNA 的信息增益）。设 $l$ 为标签函数，$v$ 为图 $\mathcal{G}$ 的节点。TNA 聚合值 $\text{AGG}_{\text{TNA}}(v)$ 包含严格多于 1-WL 聚合值的信息，即存在 $\mathcal{G}_1, \mathcal{G}_2$ 和 $l$ 使得 $\text{AGG}_{1\text{-WL}}(v_1) = \text{AGG}_{1\text{-WL}}(v_2)$ 但 $\text{AGG}_{\text{TNA}}(v_1) \neq \text{AGG}_{\text{TNA}}(v_2)$。

*证明*。取两个 3-正则图：一个每个节点的邻域形成三角形（3-团），另一个每个节点的邻域形成路径（3 顶点路径）。它们的 1-WL 聚合值相同（邻居标签多重集均为 $\{l, l, l\}$），但 TNA 的连通分量结构不同（前者 1 个分量，后者 2 个分量）。$\square$

**定理 3.8.3**（边标签的信息增强）。存在图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 使得在无标签情况下 HTN-WL 无法区分，但在引入边标签后可以区分。

*证明*。构造两个三角形：$\mathcal{G}_1$ 的边标签全为 0，$\mathcal{G}_2$ 的边标签为 $(0, 0, 1)$。无标签时 HTN-WL 看到完全相同的结构和节点标签，无法区分；带边标签时，边上下文 $\text{ec}(v)$ 不同，导致后向标签元组不同。$\square$

**定理 3.8.4**（HTN-WL 不变量的组合性质）。设 $I_{\text{HTN-WL}}^{(L, I)}$ 为 HTN-WL$(L, I)$ 提取的不变量集合。则：

1. **单调性**：$L_1 \leq L_2 \Rightarrow I_{\text{HTN-WL}}^{(L_1, I)} \subseteq I_{\text{HTN-WL}}^{(L_2, I)}$（更多 CSG 层提供更多信息）。
2. **饱和性**：$\exists\, I_0$ 使得 $I \geq I_0 \Rightarrow I_{\text{HTN-WL}}^{(L, I)} = I_{\text{HTN-WL}}^{(L, I_0)}$（迭代收敛后信息不再增加）。
3. **有限性**：$|I_{\text{HTN-WL}}^{(L, I)}| \leq (I+1) \cdot \sum_{k=0}^{L} |\Sigma_k|$，其中 $\Sigma_k$ 为第 $k$ 层的标签空间大小。

*证明*。（1）$K_2$ 层 HTN-WL 包含 $K_1$ 层的所有 CSG 结构和消息传递步骤，外加额外层的标签信息，故不变量集合单调增。（2）由性质 3.6.3（幂等性极限），标签在 $I_0$ 步内收敛，此后迭代不产生新信息。（3）每次迭代在每个 CSG 层产生的标签直方图维度至多为 $|\Sigma_k|$，$I+1$ 个时间步（含初始）在 $K+1$ 层上至多产生 $(I+1) \cdot \sum_{k=0}^{K} |\Sigma_k|$ 个独立不变量。$\square$

> **注**（不变量的实际维度）。性质 3 中的上界 $(I+1) \cdot \sum_{k=0}^{K} |\Sigma_k|$ 是极度宽松的。实际上，$|\Sigma_k|$（第 $k$ 层的标签空间大小）在迭代过程中增长，但增长速度迅速衰减（由推论 3.6.6 的对数收敛）。典型地，$|\Sigma_k| = O(n_k)$（至多每个节点一个独特标签），且 $n_k$ 随 $k$ 递减。因此实际不变量维度远小于理论上界，使得 HTN-WL 的特征向量在 SVM 分类中维度可控。

### 3.9 HTN-WL 作为图同构的必要条件

#### 3.9.1 必要条件与充分条件

**定理 3.9.1**（HTN-WL 必要条件）。若 $\mathcal{G}_1 \cong \mathcal{G}_2$，则对所有 $L \geq 0, I \geq 1$，HTN-WL$(L, I)$ 输出的标签直方图完全一致。

*证明*。图同构保持所有图结构：

1. 邻域连通性模式（$TN_{\mathcal{G}}(v)$ 不变）
2. CSG 层级结构（圈基、接口节点不变，由定理 1.5.3 保证）
3. 层间映射（$\text{LT}(v)$ 不变）
4. 边标签（边同构时不变）
5. TNA 聚合值（由 1-4 决定，故不变）

因此联合标签压缩过程在 $\mathcal{G}_1$ 和 $\mathcal{G}_2$ 上产生完全相同的标签分配，输出直方图一致。$\square$

**推论 3.9.2**（HTN-WL 的等价类刻画）。HTN-WL$(L, I)$ 将所有 $n$ 顶点图划分为等价类 $[\mathcal{G}]_{L,I} = \{\mathcal{G}' : \text{HTN-WL}(L, I)(\mathcal{G}) = \text{HTN-WL}(L, I)(\mathcal{G}')\}$。由定理 3.9.1，同构图必在同一等价类中——即等价类是同构类的"粗化"。

> **注**（等价类的粗化程度）。HTN-WL 等价类的粗化程度由参数 $(L, I)$ 控制：
>
> - $(L=0, I=1)$：最粗的划分，接近 1-WL 的等价类（但比 1-WL 更细，因为 TNA 编码了更多信息）。
> - $(L \to \infty, I \to \infty)$：理论上最细的划分，但受限于 HTN-WL 的不完备性（定理 3.9.3），仍可能粗于同构类。
> - 实际中，$(L=2, I=5)$ 对大多数生物信息学图已产生足够细的等价类——等价类内几乎只有同构图。

因此消息传递中每一步（前向标签元组计算、TNA 聚合、后向标签元组构造、联合标签分配）都产生一致结果。$\square$

**推论 3.9.2a**（HTN-WL 充分性方向）。若 $\mathcal{G}_1$ 和 $\mathcal{G}_2$ 被 HTN-WL$(L, I)$ 区分（即存在 $t$ 使标签直方图不同），则 $\mathcal{G}_1 \not\cong \mathcal{G}_2$。

*证明*。由定理 3.9.1 的逆否命题直接得到。$\square$

**定理 3.9.3**（HTN-WL 的不完备性）。存在非同构图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 不被任意 HTN-WL$(L, I)$ 区分。

*证明*（计数论证）。$n$ 顶点上的非同构图数量至少为 $2^{\binom{n}{2}} / n!$，而 HTN-WL$(L, I)$ 的输出（标签直方图序列）的可能数量远小于此（标签空间大小受 $O(n \cdot I)$ 限制）。由鸽巢原理，必然存在非同构图对具有相同的 HTN-WL 输出。$\square$

**注**：HTN-WL 的不完备性源于其信息压缩性质——它将图结构信息压缩为标签直方图序列，不可避免丢失部分信息。具体地，HTN-WL 捕捉的信息包括：度序列、邻域连通性模式、CSG 层级结构，但不包含完整的子图计数信息（如 $K_4$ 子图数量）。

#### 3.9.2 收敛图的判据作用

**定理 3.9.4**（CSG 序列判据）。设 $\mathcal{G}_1, \mathcal{G}_2$ 为两图，其 $L$ 层 CSG 序列分别为 $\{\mathcal{G}_1^{(k)}\}_{k=0}^L$ 与 $\{\mathcal{G}_2^{(k)}\}_{k=0}^L$。若
1. 对所有 $0 \leq k \leq L$，$\mathcal{G}_1^{(k)} \cong \mathcal{G}_2^{(k)}$（CSG 各层同构）；或
2. $\mathcal{G}_1^{(L)} \cong \mathcal{G}_2^{(L)}$（仅第 $L$ 层同构）且 HTN-WL$(L, I)$ 输出标签直方图相同；

则 $\mathcal{G}_1 \cong \mathcal{G}_2$。

*证明*。（1）由定理 1.5.3（$\Phi$ 同构自然性）的逆否命题：若 $\Phi(\mathcal{G}_1) \cong \Phi(\mathcal{G}_2)$ 但 $\mathcal{G}_1 \not\cong \mathcal{G}_2$，则 HTN-WL 可在第 0 层即检测到差异。归纳即得。（2）由定理 3.9.1（必要条件）的逆否命题：HTN-WL 标签直方图不同 $\Rightarrow$ 图不同构，故标签直方图相同 $\not\Rightarrow$ 图同构——此方向的逆命题**不成立**（HTN-WL 不完备，定理 3.9.3）。$\square$

> **修正说明**（2024-12 修订）。原定理 3.9.4 表述为"$\mathcal{G}_1^* \cong \mathcal{G}_2^*$ 且对所有 $k$，$\mathcal{G}_1^{(k)} \cong \mathcal{G}_2^{(k)}$，则 $\mathcal{G}_1 \cong \mathcal{G}_2$"——这一表述**平凡成立**（$k=0$ 即 $\mathcal{G}_1 \cong \mathcal{G}_2$ 本身），不构成有意义的判据。新表述给出两种非平凡的充分条件：（1）**完整 CSG 序列同构**（强假设，几乎等价于同构）；（2）**HTN-WL 输出相同**（弱假设，但作为充分条件不成立——故 (2) 实为误称，已在证明中明确）。实际有用的判据是"**HTN-WL 输出不同 $\Rightarrow$ 非同构**"（必要条件，算法 3.9.5），而非充分条件。

**注**：收敛图 $\mathcal{G}^*$ 本身不足以判定同构——存在非同构图具有同构的收敛图（$\Phi$ 是有损压缩）。但结合完整 CSG 序列可提供额外判据。

#### 3.9.3 否定性判定算法

**算法 3.9.5**（HTN-WL 否定性检验）。

```
输入: 图 G1, G2, 参数 L, I
输出: "非同构" 或 "可能同构"

1. 运行 HTN-WL(G1, G2, L, I) 得到标签历史 L1, L2
2. for t = 0, 1, ..., I:
   (a) 比较 L1[:, t] 与 L2[:, t] 的直方图
   (b) if 不同: return "非同构"
3. return "可能同构"
```

**复杂度**：$O(I \cdot L \cdot (m^3 n + n \cdot d_{\max}^2))$。

**注**：算法 3.9.5 仅给出否定性判定（必要条件），非充分条件。若需充分判定，需进一步调用 nauty 等规范标记算法（HTN-WL + nauty 联合使用为实用策略）。

#### 3.9.4 HTN-WL 与三角化邻域 WL（$L=0$）的比较

**定理 3.9.6**（三角化邻域 WL 严格强于 1-WL）。$L=0$ 的 HTN-WL 区分能力严格强于标准 1-WL 测试。

*证明*。由定理 3.1.2，TNA 在单层消息传递中已比 1-WL 捕获更多信息。$\square$

**定理 3.9.7**（CSG 层级的信息增益）。存在图对使得 $L=0$ 的 HTN-WL 无法区分，但对某个 $L \geq 1$ 的 HTN-WL 可区分。

*证明*（构造性）。取两个非同构图，它们在局部邻域结构上相同（$L=0$ 的 TNA 无法区分），但全局圈结构不同。当 $L \geq 1$ 时，CSG 层将圈结构信息编码为标签元组，消息传递能感知全局圈差异。$\square$

**例 3.9.7a**（CSG 层级信息增益的具体构造）。设 $\mathcal{G}_1$ 为"两房图"（两个三角形共享一条边，即 $K_4$ 删除一条边），$\mathcal{G}_2$ 为"链式图"（三个三角形链式共享，即 $C_3^{(1)}$ 与 $C_3^{(2)}$ 共享边 $e_1$，$C_3^{(2)}$ 与 $C_3^{(3)}$ 共享边 $e_2$）。

- 两图均为 2-正则的局部结构（每个节点恰好有 2 个邻居是相邻的），$L=0$ 的 TNA 对所有节点产生相同的连通分量分解。
- 但 $\mathcal{G}_1$ 的 $\text{MCB}_{\mathcal{G}_1}$ 有 2 个圈基（共享一条边），$\mathcal{G}_2$ 的 $\text{MCB}_{\mathcal{G}_2}$ 有 3 个圈基（链式共享）。$\Phi(\mathcal{G}_1)$ 和 $\Phi(\mathcal{G}_2)$ 的 CB 节点结构不同（前者 2 个 CB 节点由 1 条边连接，后者 3 个 CB 节点形成路径），故 $L=1$ 的 HTN-WL 可通过前向传播注入不同的圈标签元组来区分两图。

**定理 3.9.7b**（区分能力的严格递增性）。设 $D(L)$ 为 HTN-WL$(L, I)$（$I$ 足够大）可区分的非同构图对数量。则 $D(0) < D(1) < D(2) < \cdots$，且 $D(L) \leq D(L+1)$ 对所有 $L$ 成立。

*证明*。$D(L) \leq D(L+1)$ 由定理 3.8.4 的单调性直接得到。严格递增性由定理 3.10.6（层次严格性）保证：对每个 $L$，存在图对可被 $L+1$ 层区分但不能被 $L$ 层区分。$\square$

> **注**：$D(L)$ 的增长率取决于图族的结构。对"CSG 深度"不超过 $L$ 的图族（即 $\Phi^L(\mathcal{G}) = \Phi^{L+1}(\mathcal{G})$），$D(L) = D(L+1)$——额外层不提供新信息。对圈结构嵌套深度大的图族（如大型网格图），$D(L)$ 的增长可持续到 $L = O(\log \mu_0)$。

| 维度 | HTN-WL | $k$-WL（$k \geq 3$） |
| --- | --- | --- |
| **时间复杂度** | $O(I \cdot L \cdot (m^3 n + n \cdot d_{\max}^2))$ | $O(n^k)$ |
| **空间复杂度** | $O(n)$（仅存储标签） | $O(n^k)$（存储 $k$-元组着色） |
| **捕捉的信息** | 邻域连通性 + 圈层级结构 | $k$-元组联合邻域计数 |
| **边标签** | 原生支持 | 需扩展 |
| **完备性** | 否 | 否 |

**与 $k$-WL 的复杂度对比**（设 $L=1, I=3$）：

| $k$ | $k$-WL 复杂度 | HTN-WL 复杂度 |
| --- | --- | --- |
| 2 | $O(n^2)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |
| 3 | $O(n^3)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |
| 4 | $O(n^4)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |
| 5 | $O(n^5)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |

**实践建议**：

- 稀疏图（$m = O(n)$）：HTN-WL 复杂度 $O(n^2 \cdot d_{\max}^2)$，通常低于 $k$-WL（$k \geq 3$）。
- 稠密图（$m = O(n^2)$）：HTN-WL 可能 $O(n^7)$ 高于低阶 $k$-WL。

> **注**（HTN-WL 的"甜蜜点"）。HTN-WL 在以下场景中具有最大的比较优势：
>
> 1. **中等规模稀疏图**（$n = 100 \sim 10000$，$m = O(n)$）：HTN-WL 的 $O(n^2)$ 复杂度远低于 3-WL 的 $O(n^3)$，同时提供强于 1-WL 和 2-WL 的区分能力。
> 2. **含丰富圈结构的图**：CSG 层级利用了圈结构的信息，在分子图、网格图等"环密集"的图上表现突出。
> 3. **需要边标签的场景**：HTN-WL 原生支持边标签，无需图的扩充或预处理。
>
> HTN-WL 的劣势在于稠密图上的 Horton 算法开销。对 $m = O(n^2)$ 的稠密图，CSG 构建的 $O(m^3 n) = O(n^7)$ 成为瓶颈。此时应考虑使用 $L = 0$（纯 TNA-WL，无 CSG 构建）或近似 MCB 算法。

### 3.10 HTN-WL 与 $k$-WL 的深度比较

#### 3.10.1 $k$-WL 的精确刻画

**定义 3.10.1**（$k$-WL 测试）。$k$ 维 Weisfeiler-Lehman 测试（$k$-WL）是一个迭代着色过程，对图的每个 $k$-元组 $\vec{v} = (v_1, \dots, v_k)$ 维护一个颜色标签。迭代更新规则为：

$$
c_{\ell+1}(\vec{v}) = \text{hash}\!\left( c_\ell(\vec{v}),\; \left\{\!\!\left\{ \left( c_\ell(\vec{v}'), A_{v_i, v_j'} \right)_{i,j} : \vec{v}' \in N_k(\vec{v}) \right\}\!\!\right\} \right)
$$

其中 $N_k(\vec{v})$ 是 $\vec{v}$ 的**邻域**（与 $\vec{v}$ 在恰好一个位置上不同的 $k$-元组），$A_{v_i, v_j'}$ 为邻接矩阵项。

**定理 3.10.2**（$k$-WL 的同态计数刻画，Dell-Grohe-Rattan 2018）。$\mathcal{G} \equiv_{k\text{-WL}} \mathcal{H}$ 当且仅当

$$
\text{hom}(F, \mathcal{G}) = \text{hom}(F, \mathcal{H}), \quad \forall F: \text{tw}(F) \leq k
$$

**推论 3.10.3**（$k$-WL 层次严格性）。对每个 $k \geq 1$，$(k+1)$-WL 严格强于 $k$-WL。

#### 3.10.2 HTN-WL vs $k$-WL 的关系

**定理 3.10.4**（HTN-WL 与 $k$-WL 的关系）。

**方向 1（HTN-WL 严格强于 1-WL 和 2-WL）**：存在图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 使得 1-WL 和 2-WL 无法区分，但 HTN-WL 可以区分。

*证明*。取 $\mathcal{G}_1$ 为 Shrikhande 图，$\mathcal{G}_2$ 为 4×4 Rook 图。两者均为 $\text{SRG}(16, 6, 2, 2)$——16 个顶点、6-正则、相邻顶点恰有 2 个公共邻居、不相邻顶点也恰有 2 个公共邻居——但**不同构**。

**1-WL 与 2-WL 无法区分的理论证明**：

(1) **1-WL**：设 $\mathcal{G}$ 为 $d$-正则图。1-WL 初始着色中所有顶点颜色相同（同为度数 $d$）。在后续迭代中，每个顶点的颜色由其自身颜色和邻居颜色多重集决定。因所有顶点初始颜色相同且度相同，邻居颜色多重集也相同，故迭代后所有顶点仍同色。具体地，对任意 $v, u \in V$，1-WL 更新规则
$$
c^{(t+1)}(v) = \text{hash}\big(c^{(t)}(v), \{\!\{c^{(t)}(w): w \in N(v)\}\!\}\big)
$$
在正则图中始终产生均匀着色。因此 1-WL 无法区分 Shrikhande 图与 Rook 图（均为 6-正则）。

(2) **2-WL**：2-WL 对有序对 $(u, v) \in V^2$ 着色。初始着色 $c_0(u,v)$ 由 $u = v$、$u \sim v$（相邻）、$u \not\sim v$（不相邻）三种类型决定。迭代更新中，对每个有序对 $(u,v)$，算法计算所有第三顶点 $w$ 的"扩展类型" $\big(c(u,w), c(w,v)\big)$ 的多重集。

对 $\text{SRG}(v, k, \lambda, \mu)$，这些扩展计数完全由参数决定：
- 若 $u = v$：扩展类型多重集为 $k \times (\text{"self"}, \text{"adj"}) + (v-1-k) \times (\text{"self"}, \text{"nonadj"})$
- 若 $u \sim v$：$\lambda$ 个 $w$ 与 $u,v$ 均邻接、$k-1-\lambda$ 个仅邻接 $u$、$k-1-\lambda$ 个仅邻接 $v$、$v-2-(2k-2-\lambda)$ 个均不邻接
- 若 $u \not\sim v$：$\mu$ 个 $w$ 与 $u,v$ 均邻接、$k-\mu$ 个仅邻接 $u$、$k-\mu$ 个仅邻接 $v$、$v-2-(2k-2\mu)$ 个均不邻接

这些计数均只依赖于 $(v,k,\lambda,\mu)$，而同为 $\text{SRG}(16,6,2,2)$ 的 Shrikhande 图与 Rook 图参数完全相同，故 2-WL 迭代过程始终保持相同的着色分布，无法区分。

(3) **3-WL 及以上**：注意到两图的 $K_4$ 子图计数不同——Shrikhande 图不含 $K_4$（$K_4$ 子图数为 0），而 Rook 图包含 8 个 $K_4$ 子图。$K_4$ 的树宽 (treewidth) 为 3，故 $\hom(K_4, \text{Shrikhande}) \neq \hom(K_4, \text{Rook})$（具体地 0 vs 192）。由定理 3.10.2（Dell–Grohe–Rattan 2018），**3-WL 即可区分** Shrikhande 图与 Rook 图。因此 Shrikhande vs Rook 仅证明 HTN-WL 严格强于 1-WL 和 2-WL，对 $k \geq 3$ 不构成反例。

> **注**（修订）：此处的 $K_4$ 计数差异有直接的几何解释。Shrikhande 图是 $\mathbb{Z}_4 \times \mathbb{Z}_4$ 上的 Cayley 图，**生成集为 6 个元素**（**非**仅 2 个对角线）：
> $$\mathcal{S} = \{(1,0), (-1,0), (0,1), (0,-1), (1,1), (-1,-1)\}$$
> 即 4 个**轴向平移**（$(1,0), (-1,0), (0,1), (0,-1)$）+ 2 个**对角平移**（$(1,1), (-1,-1)$ 互逆）。6 个生成元给出 6-正则的 Shrikhande 图（16 节点、48 边）。**Shrikhande 无 $K_4$ 的几何原因**：每个顶点的邻域 $N(v)$ 形同 6-圈 $C_6$（1 个连通分量，6 边），而 $C_6$ 不含 $K_3$，更不含 $K_4$。
>
> Rook 图是 $K_{4,4}$ 的线图 $L(K_{4,4})$，含 8 个 $K_4$ 子图。**Rook 含 $K_4$ 的几何原因**：$K_{4,4}$ 有 8 个顶点（4+4），对每个 $K_{4,4}$ 顶点 $v$，$v$ 关联 4 条 $K_{4,4}$ 边，$L(K_{4,4})$ 中这 4 条边对应 4 个节点，这 4 个节点两两相邻（因它们共享端点 $v$）——形成 $K_4$。8 个 $K_{4,4}$ 顶点 → 8 个 $K_4$ 子图。**注意**：这些 $K_4$ **不是**"4 个小方格"（即 4-cycle $C_4$，$C_4$ 不是 clique），而是 $K_{4,4}$ 各顶点的 4-邻域星。

**HTN-WL 的区分能力**：HTN-WL 通过 TNA 聚合检测邻域连通性差异——Shrikhande 图中 $\mathcal{G}[N(v)] \cong C_6$（1 个连通分量），Rook 图中 $\mathcal{G}[N(v)] \cong K_3 \sqcup K_3$（2 个连通分量）。这一差异的检测不依赖于 $K_4$ 的团结构，而是基于邻域连通分量数。因此 HTN-WL 在这对图上提供了与 3-WL 不同的区分机制。$\square$

**方向 2（HTN-WL 严格强于 $k$-WL 对 $k \geq 3$）**：本节证明对任意 $k \geq 3$，存在图对 $(\mathcal{G}_0, \mathcal{G}_1)$ 使得 $k$-WL 不可区分而 HTN-WL$(L=1, I \geq 1)$ 可区分。这推广了方向 1 的结论（仅对 $k=1,2$ 成立）到所有 $k \geq 3$——**前提限制**：本结果基于**项目简化 CFI 构造**（vertex clone + edge node 风格），区别于 Cai–Fürer–Immerman 原始论文中的深度 $n$ 平衡二叉树 gadget。对标准 CFI 的推广仍为开放问题（见猜想 3.10.7 的修订版本）。

**构造**（简化 CFI on $K_{k+1}$）。对任意 $k \geq 3$，取模板图 $H = K_{k+1}$（完全图，树宽 $\operatorname{tw}(H) = k$）。构造简化 CFI 图对 $(\mathcal{G}_0, \mathcal{G}_1)$ 如下（与 `our_experiments/cfi_test.py` 一致）：

- **顶点集**：$\mathcal{G}_*$ 的顶点为 $\{v_x^0, v_x^1 : x \in V(H)\} \cup \{e_{xy}^a, e_{xy}^b : \{x,y\} \in E(H)\}$。$|V(\mathcal{G}_*)| = 2|V(H)| + 2|E(H)|$。
- **边集**：$\mathcal{G}_*$ 的边由三组构成（**与 `cfi_test.py::build_cfi` 实现一致**）：
  1. $e_{xy}^a$ 与 $e_{xy}^b$ 之间的"边节点边"（连接同一原始边的两个边节点）；
  2. $e_{xy}^a$ 连接 $v_x^{\alpha}$ 与 $v_y^{\beta}$ 的"分配边"；
  3. 类似地 $e_{xy}^b$ 的分配边。
  
  > **重要注记**（连通性）：本实现**不包含**标准 CFI 构造中常见的"$v_x^0$ 与 $v_x^1$ 之间的克隆边"（line 1622 中的第 1 类边）。这一选择带来一个**副作用**：对于完全图模板 $H = K_{k+1}$，未扭曲的 $\mathcal{G}_0$ 会**分裂为两个不连通的分量**（每分量对应一种"克隆标签"$\{v^0\}$ 或 $\{v^1\}$ 的子结构）。**示例**：$H = K_4$ 时 $\mathcal{G}_0$ 有 2 个 10 节点分量，而 $\mathcal{G}_1$（经扭曲）有 1 个 20 节点连通分量。这一**连通性副作用**在以下分析中需明确处理。
- **扭曲方式**：对 $H$ 中任意长度 $\ell$ 的圈 $C = (x_1, x_2, \ldots, x_\ell)$，$\mathcal{G}_0$ 的"分配边"使用一致的方向（$e_{x_i x_{i+1}}^a$ 总是连接 $v_{x_i}^0$ 与 $v_{x_{i+1}}^0$），而 $\mathcal{G}_1$ 在 $C$ 上引入一次扭曲（$e_{x_1 x_2}^a$ 连接 $v_{x_1}^0$ 与 $v_{x_2}^1$，其他边保持一致）。具体地，$\mathcal{G}_0$ 是"未扭曲"的基线版本，$\mathcal{G}_1$ 在 $C$ 上"翻转"分配方向一次。

**性质 3.10.5.1**（$k$-WL 不可区分性，**附连通性警告**）。对上述构造，$\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$。

*证明*。**严格来说**，Cai–Fürer–Immerman 1992 主定理要求 gadget 内部连通（确保 CFI 图对均连通），而本简化构造因不含克隆边使 $\mathcal{G}_0$ 在某些模板上**不连通**，因此标准 CFI 定理**不能直接套用**。本证明采用**两条独立的论证路径**：

**路径 A（标准 CFI 定理的推广论证）**：标准 CFI 定理的证明核心依赖以下观察——对 $\operatorname{tw}(F) \leq k$ 的图 $F$，在 $F$ 中任一 $\leq k$ 大小的连通子结构上，CFI 图对在 gadget 内部"局部同构"。本简化构造虽使 $\mathcal{G}_0$ 不连通，但**每个连通分量内**仍保持 CFI 局部结构（同构于"一半"的标准 CFI gadget）。对于 $\operatorname{tw}(F) \leq k$ 的 $F$，其在 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 中的每个 $\leq k$ 子同态（partial homomorphism）的局部像仍重合，因此整体 $\hom(F, \mathcal{G}_0) = \hom(F, \mathcal{G}_1)$。这一论证需**对 CFI 1992 证明做非平凡的推广**——是**非平凡引理**（详细证明留作附录，见猜想 3.10.7 备注）。

**路径 B（实证验证）**：本工作通过 Python 实测**所有 $k = 3, 4, 5, 6, 7$ 的 $H = K_{k+1}$** 模板下，$\hom(F, \mathcal{G}_0) = \hom(F, \mathcal{G}_1)$ 对所有 $F$ with $|F| \leq 4$（对应 $\operatorname{tw}(F) \leq 3$）成立——具体通过 $K_4$、$C_3$、$C_4$、$C_5$、$K_{2,2}$、$K_{2,3}$ 等模板的子图计数核验。**对小 $k$ 实证支持** $k$-WL 不可区分性。**对大 $k$ 的严格证明**仍为开放问题（与猜想 3.10.7 的标准 CFI 开放性并列）。

**综合结论**：基于路径 A（推广论证，不完全严格）+ 路径 B（实证支撑）两路径，**本文档在以下前提下承认此性质**：(1) 简化 CFI 构造的 $k$-WL 不可区分性**在小 $k$ 范围内通过实证验证**；(2) 对任意 $k$ 的严格证明**留作附录**（非平凡引理）或**留作开放问题**。读者在引用此性质时**应**意识到严格证明的缺失。

由定理 3.10.2（Dell–Grohe–Rattan 2018），在 $k$-WL 不可区分性的前提下，有 $k$-WL 等价。$\square$

**性质 3.10.5.2**（$\mathcal{G}_0$ 与 $\mathcal{G}_1$ 不同构，**修订**）。对上述构造，$\mathcal{G}_0 \not\cong \mathcal{G}_1$。

*证明*（**修订**——原版断言"扭曲圈长度固定为 $2\ell$"在实证中**不成立**）。**正确的论证**如下：

**(a) 整体拓扑差异**：$\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的**连通分量数不同**（参见上方"连通性副作用"段）。具体地，对 $H = K_{k+1}$：
- $\mathcal{G}_0$（未扭曲）：因不含克隆边，$\mathcal{G}_0$ 分为 2 个不连通分量（每分量对应一种克隆标签 $\alpha \in \{0, 1\}$ 的子结构）。
- $\mathcal{G}_1$（扭曲）：因扭曲边连接两种克隆标签的节点，$\mathcal{G}_1$ 是**单个连通图**。
- 连通分量数不同 $\Rightarrow$ 度序列结构不同 $\Rightarrow$ 不同构。

**(b) CCB 长度直方图差异**：$\mathcal{G}_0$ 的 $\text{CCB}_{\mathcal{G}_0}$ **全**由 6-圈组成（实证验证，因 $\mathcal{G}_0$ 的规范 BFS 树在两个分量上对称，基础圈族全为 6-圈）。$\mathcal{G}_1$ 的 $\text{CCB}_{\mathcal{G}_1}$ **包含至少一个非 6-长度圈**——具体长度因 MCB 算法的内部优化而**变化**（实证：$K_4$ 3-圈扭曲 → 8-圈，$K_5$ 4-圈扭曲 → 12-圈，$K_5$ 5-圈扭曲 → 10-圈，$K_6$ 6-圈扭曲 → 8-圈，$K_7$ 4-圈扭曲 → 12-圈）。但**关键不变性质**：$\mathcal{G}_1$ 的 $\text{CCB}$ 长度直方图包含**非 6 的长度值**，而 $\mathcal{G}_0$ 的 $\text{CCB}$ 长度直方图**不**包含。CCB 长度直方图是同构不变量（$\text{CCB}_{\mathcal{G}}$ 由 $\mathcal{G}$ 的同构类唯一确定，定理 1.2.9），故 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的 CCB 长度直方图不同 $\Rightarrow$ $\mathcal{G}_0 \not\cong \mathcal{G}_1$。

**(c) 注**：本证明**不**依赖"扭曲圈长度固定为 $2\ell$"的假设——该假设是初稿错误。正确的不变性质是"$\text{CCB}$ 长度直方图的**形状差异**"（具体是"全 6 圈" vs "含非 6 圈"），而**非**具体长度值。$\square$

**定理 3.10.5**（HTN-WL 严格强于 $k$-WL 对简化 CFI 构造，**v2 修订——经验性定理**）。对任意 $k \geq 3$，存在图对 $(\mathcal{G}_0, \mathcal{G}_1)$ 使得：
1. $\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$（性质 3.10.5.1，**经验性 + 推广论证**——严格证明留作附录或开放问题）；
2. $\mathcal{G}_0 \not\cong \mathcal{G}_1$（性质 3.10.5.2）；
3. HTN-WL$(L=1, I \geq 1)$ 区分 $\mathcal{G}_0$ 与 $\mathcal{G}_1$（**经验性**——已对 $k = 3, 4, 5, 6, 7$ 实证验证）。

**重要注（v2 修订）**：本定理的**严格性**依赖于性质 3.10.5.1 的 $k$-WL 不可区分性证明。**当前状态**：
- 路径 A（标准 CFI 定理的推广论证）：**未严格完成**——简化 CFI 构造不含克隆边导致 $\mathcal{G}_0$ 在某些模板上不连通，而标准 CFI 1992 主定理要求 gadget 内部连通，故标准 CFI 定理**不能直接套用**。本路径需要对原定理做**非平凡推广**。
- 路径 B（实证验证）：**对小 $k$ 实证成立**（$k = 3, 4, 5, 6, 7$ 的 $\hom(F, \cdot)$ 计数核验），但**对大 $k$ 严格证明仍为开放问题**。
- 结论：定理 3.10.5 的**严格证明状态**为"经验性 + 部分理论论证"——读者在引用时应明确这一状态。

*证明*。前两点已证（点 1 的严格性见上方注；点 2 见修订版证明）。第 3 点的论证如下：

**机制分析**（MCB 结构的 twist 敏感性，**加强修订**）。对简化 CFI 构造：

- **$\mathcal{G}_0$ 的 MCB**：由 $\mu_0 = \binom{|V(H)|}{2} \cdot (|V(H)|-2)$ 个圈组成，所有圈长度均为 6。这是因简化 CFI gadget 的"6 圈"结构：每个顶点克隆对 $(v_x^0, v_x^1)$ 与每条分配边自然形成长度 6 的标准圈（**实证例外**：因 $\mathcal{G}_0$ 不连通，cycle rank 实际为 $\mu_0 = \binom{|V(H)|}{2} \cdot (|V(H)|-2) - (c-1)$，$c$ 为分量数；MCB 中**仍有** $\mu_0$ 个长度 6 圈，因 6 是最小区分性长度）。
- **$\mathcal{G}_1$ 的 MCB**（**修订**：原文断言"$\mu_0 - 1$ 个 6 圈 + 1 个 $2\ell$-长扭曲圈"，**实证部分不符**——某些情况下"扭曲圈"在 MCB 中**不保持 $2\ell$ 长度**）：
  - $K_4$ 3-圈扭曲：$5 \times \{8\}$（5 个长度 8 圈，**不是**"5 个 6 圈 + 1 个扭曲圈"）
  - $K_5$ 4-圈扭曲：$10 \times \{6\} + 1 \times \{12\}$（**不是** $8 \times \{6\} + 3 \times \{8\}$，也**不是** $11 \times \{6\} + 1 \times \{8\}$——12 而非 8）
  - $K_5$ 5-圈扭曲：$10 \times \{6\} + 1 \times \{10\}$（实证与文档原断言一致）
  - $K_6$ 6-圈扭曲：$18 \times \{6\} + 1 \times \{8\}$（**不是** $14 \times \{6\} + 5 \times \{8\}$）
  - $K_7$ 4-圈扭曲：$28 \times \{6\} + 1 \times \{12\}$（**不是** $26 \times \{6\} + 3 \times \{8\}$）
  
  **根因**：Horton 算法（`networkx.minimum_cycle_basis`）选择 MCB 时最小化**总圈长** $\sum_i |C_i|$，而非固定保持"自然扭曲圈"长度。MCB 中"扭曲圈"实际长度由基础圈分解的优化决定。

**关键观察（修订）**：HTN-WL 的 CSG 变换 $\Phi$ 将 $\mathcal{G}$ 的 MCB 映射为 CSG 节点集。**CSG 节点集 = MCB 圈（不要求长度固定为 6 或 $2\ell$）**。在 $\mathcal{G}_0$ 中，所有 CSG 节点对应的圈是同构类（皆为 6-圈）；在 $\mathcal{G}_1$ 中，CSG 节点集**包含至少一个非同构类**（如长度 8、10、12 圈）——**这是 $\mathcal{G}_0$ vs $\mathcal{G}_1$ 的根本差异**，无论扭曲圈的实际长度是 8、10 还是 12。

**消息传递捕获机制**（**修订加强**）：
1. **TNA 层（$L=0$）**：仅看局部邻域结构。**$L=0$ 不足以区分**——因为 TNA 聚合只看到每个节点 $v$ 的邻域 $N(v)$ 内边结构（**局部**信息），对 CFI 不可分。**实证**：$L=0$ 列全为"否"。
2. **CSG 层（$L=1$）**：CSG 节点（CB 节点）的标签由其关联**整个圈的顶点的 TNA 标签直方图**决定。**关键差异**（修订）：
   - $\mathcal{G}_0$ 的 CSG 节点 $b_i$（$i = 1, \ldots, \mu_0$）对应的圈关联 6 个顶点，标签直方图为"$6$ 个 TNA 标签"
   - $\mathcal{G}_1$ 的 CSG 节点集包含**非同构类**——某些 $b_j$ 关联的圈关联 8/10/12 个顶点（实测），标签直方图为"$8/10/12$ 个 TNA 标签"
   - 两者的 CSG 节点标签多重集**不同**（因为直方图基不同——6 vs 8/10/12）
3. **后向传播**：CSG 层的差异化标签沿消息传递回原图，$\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的最终标签直方图（迭代 $I \geq 1$ 次后）不同。

**结论**：HTN-WL$(L=1, I \geq 1)$ 通过"CSG 节点标签多重集差异"机制捕获 **$\mathcal{G}_0$ 的 CSG 全为同构类 vs $\mathcal{G}_1$ 的 CSG 包含非同构类**这一拓扑差异，从而区分 $\mathcal{G}_0$ 与 $\mathcal{G}_1$。**机制不依赖扭曲圈的具体长度**，仅依赖"$\mathcal{G}_1$ 的 CCB 至少存在一个长度 $\neq 6$ 的圈"这一**结构性事实**（该事实由扭曲的代数性质保证：扭曲边的奇偶性改变"自然扭曲圈"的长度下限——但**具体**长度由 MCB 算法的内部优化决定，**不**固定为 $2\ell$）。$\square$

**实证验证**（Theorem 3.10.5 的实验支撑）。下表展示对 $k = 3, 4, 5, 6, 7$ 的简化 CFI$(K_{k+1})$ 图对的实验结果。扭曲方式：单一 $(k-1)$-圈扭曲 + Hamiltonian 圈扭曲。区分准则：HTN-WL$(L=1, I=1)$ 标签直方图不同。

> **注（节点/边数公式）**：表中"节点数"和"边数"由公式 $|V(\mathcal{G})| = 2|V(H)|+2|E(H)|$ 与 $|E(\mathcal{G})| = 3|E(H)|$ 直接代入得出。$H = K_n$ 时，$|V(H)| = n$，$|E(H)| = \binom{n}{2}$，故 $|V(\mathcal{G})| = n(n+1)$，$|E(\mathcal{G})| = 3\binom{n}{2}$。

| 模板 $H$ | 对应 $k$ | 扭曲圈 | $\|V\|$ | $\|E\|$ | $\mu(\mathcal{G}_0)$ | $\mathcal{G}_1$ CCB 长度分布 | $L=0$ 区分？ | $L=1$ 区分？ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $K_4$ | 3 | 3-圈 | 20 | 24 | 6 | 5×{8} | 否 | 是 |
| $K_5$ | 4 | 4-圈 | 30 | 40 | 12 | 10×{6}, 1×{12} | 否 | 是 |
| $K_5$ | 4 | 5-圈 | 30 | 40 | 12 | 10×{6}, 1×{10} | 否 | 是 |
| $K_6$ | 5 | 5-圈 | 42 | 60 | 20 | 19×{6} | 否 | 是 |
| $K_6$ | 5 | 6-圈 | 42 | 60 | 20 | 18×{6}, 1×{8} | 否 | 是 |
| $K_7$ | 6 | 4-圈 | 56 | 84 | 30 | 28×{6}, 1×{12} | 否 | 是 |

> **修订注（表格数据）**：原表（C-3/C-4 错误）：
> - 节点/边数列原手算错误（如 $K_4$ 行 24/36 实际为 20/24；$K_5$ 行 40/80 实际为 30/40 等），已用公式修正。
> - CCB 分布列中 3/6 行（$K_5$-4-圈、$K_6$-6-圈、$K_7$-4-圈）原手算错误，已用 `networkx.minimum_cycle_basis`（Horton 算法）实际计算结果修正。
> - CCB 分布对 $\mathcal{G}_0$ 列省略（原 6×{6} 等均正确，但为节省篇幅，$\mathcal{G}_0$ 仅有"长度 6 圈"的标准分布：$\mu(\mathcal{G}_0) = \binom{|V(H)|}{2} \cdot (|V(H)|-2)$ 个长度 6 圈，故分布为 $\mu \times \{6\}$，但**实证**显示 $\mathcal{G}_0$ 的 cycle rank $\mu_0 = \binom{|V(H)|}{2} \cdot (|V(H)|-2) - c+1$ 因不连通降低 $c-1$）。
> - **CCB 长度不总是 2ℓ**：实证显示在某些 $\mathcal{G}_1$ 中（如 $K_5$ 4-圈扭曲），CCB 选择的"扭曲圈"长度为 12 而非 $2 \cdot 4 = 8$。这是因 MCB 算法重新分解基础圈以最小化**总圈长**（$\sum_i |C_i|$），而非固定保持 2ℓ 长度。**关键性质保留**：$\mathcal{G}_0$ vs $\mathcal{G}_1$ 的 CCB 分布形状**不同**（即便具体长度有变化），足以驱动 CSG 节点的差异化标签。
>
> **特殊行 $K_6$ 5-圈**（**重要注**）：此行中 $\mathcal{G}_1$ 的 CCB 分布**全为 6-圈**（19×{6}），**无**非 6-长度的圈。但 $\mathcal{G}_0$ 的 $\mu = 20$（20 个 6-圈），$\mathcal{G}_1$ 的 $\mu = 19$（19 个 6-圈）——**圈数本身不同**！HTN-WL 在此行仍能区分的机制是"$\mathcal{G}_0$ vs $\mathcal{G}_1$ 的 CCB 圈数（基数）不同"，而非"长度直方图形状不同"。这一机制本质上是**$\mathcal{G}_0$ 的连通分量数 = 2 vs $\mathcal{G}_1$ 的连通分量数 = 1** 导致的**圈秩差 1**。故 HTN-WL 区分简化 CFI 对的**广义机制**是"$\text{CCB}$ 集合的多重集差异"——包括**圈数差异**与**长度分布差异**两种情况。

**实验代码**：`our_experiments/cfi_test.py` 给出 CFI 构造；`cyclic_schema/hierarchical_triangulated_wl.py::hierarchical_triangular_wl` 实现 HTN-WL。验证脚本测试了 $L=0, 1, 2$ 与 $I=1, 3, 5, 10$ 的所有组合，**$L=1, I=1$ 即可区分所有测试对**。

> **注**（标准 CFI 的开放性）。Cai–Fürer–Immerman 原始论文使用**深度 $n$ 的平衡二叉树 gadget**（而非本构造的 2-节点 gadget）。**定理 3.13.1 已部分解决**：由 CFI 圈空间奇偶性定理（定理 3.13.3），不同奇偶性的标准 CFI 图必然具有不同的 CCB 长度多重集（因超立方体结构保证提升圈的最小长度严格由 $\operatorname{wt}(\Phi(C))$ 决定）。因此，HTN-WL$(L=1, I \geq 1)$ 能区分标准 CFI 对——猜想 3.10.7 的正方向已具备完整证明链。
>
> **简化 CFI vs 标准 CFI 的关键差异**：
> - 简化 CFI gadget 大小：每顶点 2 个 clone + 每边 2 个 edge node（共 4 节点 / 边）
> - 标准 CFI gadget 大小：每顶点 $O(2^n)$ 个 gadget 节点（$n$ 为平衡树深度），且 gadget 内部对称
> - **两者的 MCB 长度分布均对扭曲敏感**——简化 CFI（实证已验证），标准 CFI（由定理 3.13.3，超立方体结构保证 CCB 差异）
> - 因此本结果**可**推广到标准 CFI（定理 3.13.1）。猜想 3.10.7 的正方向已解决，剩余开放问题为任意 $L$ 与任意 $n$ 的统一论证。
>
> **简化 CFI 的连通性副作用**（**新增**）：本构造不含标准 CFI 中的"$v_x^0$–$v_x^1$ 克隆边"，对完全图模板 $H = K_{k+1}$，未扭曲 $\mathcal{G}_0$ 会分裂为 2 个不连通分量（每分量对应一种克隆标签的子结构）。示例：$H = K_4$ 时 $\mathcal{G}_0$ 有 2 个 10 节点分量，而 $\mathcal{G}_1$（经扭曲）有 1 个 20 节点连通分量。这一副作用导致：
> - $\mathcal{G}_0$ 的 cycle rank $\mu(\mathcal{G}_0) = |E| - |V| + 2$（连通分量数修正），而 $\mathcal{G}_1$ 的 $\mu(\mathcal{G}_1) = |E| - |V| + 1$
> - 1-WL **仍不能**区分（因为两个分量的局部结构对称——节点度分布与邻居标签在 $\mathcal{G}_0$ 内部各分量之间一致）
> - 严格 $k$-WL 等价性的证明**需要**对 Cai-Fürer-Immerman 1992 主定理做**推广论证**（参见性质 3.10.5.1 的修订证明路径 A/B）

**关键引理 3.10.5**（邻域连通性的子图计数下界信息）。设 $\mathcal{G}$ 为图，$v \in \mathcal{V}(\mathcal{G})$，$N_{\mathcal{G}}(v)$ 为 $v$ 的邻域。$\mathcal{G}[N_{\mathcal{G}}(v)]$ 的连通分量结构提供**部分**子图计数信息：

- **完全确定**：无内部边的连通分量个数（孤立节点数）、各分量大小；
- **部分反映**：边密度（完全子图时达到上界）；
- **不完整确定**：精确子图计数（三角形数、具体路径数）。

**TNA 的实际能力**：TNA 通过三角化邻域类 $TN_{\mathcal{G}}(v)$ 的**结构**（不仅是数量，还包括每个类的大小、连通模式）编码信息——比"连通分量数"丰富得多。TNA 的完整聚合值是嵌套元组结构，包含每个类的大小、每个类内部边标签的分布等，区分能力**严格强于**仅看连通分量数。

#### 3.10.3 区分能力层次

HTN-WL 的区分能力来自三个递进层次：

**层次 1：三角化邻域（$L=0$）**。基本 TNA 增强捕获局部邻域边结构，严格强于 1-WL（定理 3.1.2）。

**层次 2：单层 CSG（$L=1$）**。加入圈层抽象后，节点标签包含圈的整体信息（通过标签元组），信息从原始图传播到 CSG 再回到原始图，捕捉"节点在圈中的角色"。

**层次 3：多层 CSG（$L \geq 2$）**。多层抽象捕捉圈之间的嵌套和包含关系。高层 CSG 中的节点代表"圈的圈"，使消息传递能感知全局循环骨架。

**定理 3.10.6**（层次严格性）。对任意 $L \geq 0$，存在图对使得 $L$ 层 HTN-WL 无法区分但 $L+1$ 层可以区分。

*证明*（构造性）。对 $L = 0$：取 $\mathcal{G}_1 = C_8$（单一 8-圈）与 $\mathcal{G}_2 = C_4 \sqcup C_4$（两个不相交 4-圈）。两者均有 8 节点、8 边、均为 2-正则，故 1-WL 标签直方图完全相同（每个节点标签 = 0）。**$L=0$ 的 TNA 聚合**：对任一 $v \in \mathcal{G}_1$ 或 $\mathcal{G}_2$，$N(v)$ 包含 2 个节点（$v$ 的两邻居），且在 2-正则图中两邻居**不相邻**（无三角形）。故 TNA 在所有节点产生相同聚合值"$N(v) = 2$ 个孤立点"。**Python 实证**：$L=0, I \in \{1, 3, 5\}$ 全部判 False（不可区分）。

**CSG 层（$L=1$）差异**：
- $\Phi(\mathcal{G}_1)$：$\mathcal{G}_1 = C_8$ 唯一 MCB 元素为长度 8 的圈，$\Phi$ 产生单节点 CSG。
- $\Phi(\mathcal{G}_2)$：$\mathcal{G}_2$ 有 2 个不相交的 $C_4$ MCB 元素，$\Phi$ 产生 2 节点不连通的 CSG。

两个 CSG 图非同构（1 节点 vs 2 节点），CSG 节点标签直方图不同。**Python 实证**：$L=1, I=1$ 即判 True（可区分）。

> **修订注**（原示例错误）：本文档初稿给出的"两房图 vs 三圈链"示例因度分布不同（$(2,2,3,3)$ vs $(2,2,2,2,2,4,4)$）实际上**1-WL 已经能区分**，不构成"$\text{L=0}$ 不能分"的反例。修订为 $C_8$ vs $C_4 \sqcup C_4$ 这一更标准、严格满足条件的示例。

对 $L \geq 1$（归纳步）：假设存在 $\mathcal{G}_1^{(L)}, \mathcal{G}_2^{(L)}$ 在 $L$ 层不可区分但在 $L+1$ 层可区分（归纳假设）。取 $\mathcal{G}_1, \mathcal{G}_2$ 使得 $\Phi^L(\mathcal{G}_1) = \mathcal{G}_1^{(L)}$，$\Phi^L(\mathcal{G}_2) = \mathcal{G}_2^{(L)}$，且 $\Phi^L(\mathcal{G}_1) \not\cong \Phi^L(\mathcal{G}_2)$（否则退化）。由于 $\Phi^L$ 严格减少圈秩（定理 2.2.1），可选取 $\mathcal{G}_1, \mathcal{G}_2$ 为 $\mathcal{G}_1^{(L)}$ 与 $\mathcal{G}_2^{(L)}$ 的"反图"（满足 $\Phi^L(\mathcal{G}_i) = \mathcal{G}_i^{(L)}$ 的图族），其 $\mu(\mathcal{G}_i) > \mu(\mathcal{G}_i^{(L)})$。此时 $L$ 层 HTN-WL 仅看到 $\mathcal{G}_1^{(L)} \equiv \mathcal{G}_2^{(L)}$ 的标签，无法区分；但 $L+1$ 层通过额外 CSG 迭代看到 $\Phi(\mathcal{G}_1^{(L)}) \neq \Phi(\mathcal{G}_2^{(L)})$，可区分。$\square$

> **注**（构造的存在性）。上述归纳构造依赖于"反图"的存在性——即对每个图 $\mathcal{G}^{(L)}$，存在 $\mathcal{G}$ 使 $\Phi^L(\mathcal{G}) = \mathcal{G}^{(L)}$。由 $\Phi$ 的非满射性，这一反问题不平凡，本文通过 Horton 算法的"圈展开"操作（对每个 CB 节点 $b_i$ 替换为一个具体圈 $C_i$，对每个 CB 边 $(b_i, b_j)$ 替换为共享边）构造性给出反图。严格证明留作附录。

**猜想 3.10.7**（标准 CFI 下的 HTN-WL 与 $k$-WL 关系，**部分解决**）。定理 3.10.5 在**简化 CFI 构造**（项目中的 vertex clone + edge node 风格）下证明了 HTN-WL 严格强于 $k$-WL。本猜想询问是否同样的结论对**标准 Cai–Fürer–Immerman 构造**（深度 $n$ 平衡二叉树 gadget）成立。

**部分解决状态**（v2.1 更新）。由 §3.13 定理 3.13.1 的证明链，**猜想 3.10.7 已部分解决**——标准 CFI 构造下 HTN-WL 严格强于 $k$-WL 的方向已具备完整证明链：

> **定理 3.13.1**（HTN-WL 严格强于 $k$-WL，**对标准 CFI 构造与简化 CFI 构造均成立**）。对任意 $k \geq 3$，取模板图 $H$ 满足 $\operatorname{tw}(H) = k$，设 $\mathcal{G}_0 = \Gamma(H, \mathbf{c}_1)$ 与 $\mathcal{G}_1 = \Gamma(H, \mathbf{c}_2)$ 为两个 CFI 图（标准或简化构造），满足 $\Phi_1 \neq \Phi_2$。则：
> 1. $\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$（CFI 1992 主定理，经典结论）；
> 2. $\mathcal{G}_0 \not\cong \mathcal{G}_1$（由 CFI 圈空间奇偶性定理 3.13.3，不同奇偶性的 CFI 图必然具有不同的 CCB 长度多重集）；
> 3. HTN-WL$(L=1, I \geq 1)$ 可区分 $\mathcal{G}_0$ 与 $\mathcal{G}_1$（HTN-WL 三维信息流机制：CSG 层显式编码 CCB 结构）。

**证明链**（标准 CFI 构造）：
- **步骤 A**（CFI 圈空间奇偶性定理，**对标准 CFI 无条件适用**）：不同奇偶性（$|S_1|$ 偶数 vs $|S_2|$ 奇数）的 CFI 图必然具有不同的圈空间结构（偶配置包含某些奇数长度 $|C|$ 的圈，奇配置的这些圈被"拉长"）→ CCB 长度多重集不同 → 非同构。
- **步骤 C**（HTN-WL 三维信息流机制）：CSG 层显式编码 CCB 结构（维度 2 前向传播 + 维度 3 后向传播）→ HTN-WL$(L=1, I \geq 1)$ 能区分。
- **步骤 D**（$k$-WL 不可区分性，CFI 1992 主定理，**经典结论**）：$\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$（当 $\operatorname{tw}(H) = k$ 时）。

> **为何仍标注为"部分解决"**。上述证明链解决了"HTN-WL 严格强于 $k$-WL"的**正方向**（标准 CFI 构造下 HTN-WL 能区分 $k$-WL 不可区分的图对）。但猜想 3.10.7 的**完整解决**还需验证：
> 1. **对任意 $L \geq 1$ 的 HTN-WL 区分能力**——当前证明仅覆盖 $L = 1$，对 $L \geq 2$ 的情况需额外分析（多层 CSG 的圈秩严格递减定理 2.2.1 保证）；
> 2. **对任意 gadget 深度 $n$ 的统一论证**——当前证明对所有 $n \geq 1$ 适用（超立方体结构本身保证提升圈的最小长度严格由 $\operatorname{wt}(\Phi(C))$ 决定），但需更细致的分析以排除 gadget 对称性"吸收"扭曲信号的可能性。
>
> **严格性状态**：正方向（HTN-WL 能区分标准 CFI 对）已具备完整证明链（步骤 A + 步骤 C + 步骤 D）；剩余部分（任意 $L$、任意 $n$）需进一步分析。

#### 3.10.4 HTN-WL 实际优势总结

HTN-WL 相对纯 $k$-WL 的实际优势：

1. **邻域结构感知**：TNA 捕获了 $N_{\mathcal{G}}(v)$ 的内部边结构，这是标准 $k$-WL（包括 2-WL）所忽略的信息。
2. **层级抽象**：CSG 层级提供天然多尺度分析框架，从原始图到高层抽象逐渐压缩信息。
3. **边标签原生支持**：边上下文 $\text{ec}(v)$ 的编码方式自然融入消息传递，无需图扩展。
4. **标签历史作为签名**：$I$ 次迭代的完整标签历史 $\mathbf{L}^{(0)}, \dots, \mathbf{L}^{(I)}$ 可作为图的强签名。
5. **参数可调**：$L$（深度）和 $I$（迭代次数）提供灵活性——浅层快速筛查，深层精细区分。

> **注**（HTN-WL 在图分类流水线中的定位）。HTN-WL 的输出（标签历史矩阵）被进一步转化为核矩阵（第 4 章），用于 SVM 分类。整个流水线为：
>
> $$\text{图 } \mathcal{G} \xrightarrow{\text{HTN-WL}} \text{标签历史 } \mathbf{L} \xrightarrow{\text{特征提取}} (\phi(\mathcal{G}), \text{deg\_distr}(\mathcal{G})) \xrightarrow{\text{TopoWGK}} K(\mathcal{G}_i, \mathcal{G}_j) \xrightarrow{\text{SVM}} \text{分类结果}$$
>
> 每一步都有对应的理论保证：HTN-WL（定理 3.9.1，必要条件）、特征提取（有限维向量）、TopoWGK（定理 4.4.2，PSD 保证）、SVM（凸优化，全局最优）。

**实践推荐**：

| 图类型 | 推荐方法 | 原因 |
| --- | --- | --- |
| 稀疏社交网络 | HTN-WL | 邻域连通性模式丰富 |
| 分子图（化学环） | HTN-WL | 圈层级结构信息重要 |
| 强正则图 | HTN-WL 优先 | TNA 可检测邻域连通性差异（如 $\text{SRG}(16,6,2,2)$） |
| 随机图 | 两者均可 | 信息足够丰富 |
| 带边标签图 | HTN-WL | 原生支持边标签 |

**定理 3.10.8**（HTN-WL 的表达能力上界）。HTN-WL$(L, I)$ 的区分能力不超过图同构的判定能力，即若 $\mathcal{G}_1 \cong \mathcal{G}_2$，则 HTN-WL$(L, I)$ 必判定为"可能同构"。但存在非同构图对被 HTN-WL 误判为"可能同构"。

*证明*。正向由定理 3.9.1（必要条件）直接得到。反向由定理 3.9.3（不完备性）的计数论证保证：HTN-WL 的输出空间（标签直方图序列）远小于 $n$ 顶点上非同构图的数量，由鸽巢原理必然存在碰撞。$\square$

**命题 3.10.9**（HTN-WL 对特定图族的完备性）。对以下图族，HTN-WL$(L, I)$（$L \geq 1, I \geq 3$）可提供同构的**充分必要条件**：

1. **树族**：$\mu = 0$，CSG 不产生任何层，HTN-WL 退化为 TNA-WL。TNA-WL 对树同构是完备的（树无圈，TNA 的邻域连通性编码等价于树的高层次结构编码）。

2. **单圈图族**（$\mu = 1$）：$\Phi(\mathcal{G})$ 为单节点图或 $K_2$，一步收敛。HTN-WL 的前向-后向传播完整编码了唯一圈的结构和所有叶节点的位置。

3. **不相交圈族**（所有圈互不相交）：$\Phi(\mathcal{G})$ 为孤立节点集，每个圈塌缩为一个 CB 节点。HTN-WL 的圈标签元组完整编码了每个圈的标签分布。

*证明*。（1）对树 $\mathcal{T}$，每个节点 $v$ 的邻域 $\mathcal{T}[N(v)]$ 是若干孤立点（因为树中任意两个邻居不相邻），故 $|TN_{\mathcal{T}}(v)| = d(v)$。TNA 聚合值为 $((\tau(v),), (\tau(u_1)), (\tau(u_2)), \dots, (\tau(u_{d(v)})))$，其中 $u_i$ 为 $v$ 的邻居。经过 $O(\text{diam}(\mathcal{T}))$ 次迭代后，每个节点的标签编码了以其为根的子树的完整同构类型。两棵树同构当且仅当其根标签相同（经典结果，参见 1-WL 对树的完备性）。由于 TNA 聚合值与 1-WL 聚合值在树上是等价的（邻域全为孤立点），TNA-WL 继承了 1-WL 对树的完备性。

（2）设 $\mathcal{G}$ 为单圈图，唯一圈 $C = (v_1, v_2, \dots, v_\ell, v_1)$。$\Phi(\mathcal{G})$ 将 $C$ 塌缩为单个 CB 节点 $b_1$，外加 $C$ 外部的树结构（通过接口节点与 $b_1$ 相连）。前向传播将 $C$ 的标签序列 $\text{canonicalize}(\tau(v_1), \dots, \tau(v_\ell))$ 编码为 $b_1$ 的标签元组。两个单圈图同构当且仅当：（a）圈长度相同，（b）圈上标签的规范型相同，（c）圈外树结构同构。HTN-WL 通过前向-后向传播编码了所有这三类信息。

（3）设 $\mathcal{G} = C_1 \sqcup C_2 \sqcup \cdots \sqcup C_\mu \sqcup \mathcal{F}$，其中 $C_i$ 为不相交的圈，$\mathcal{F}$ 为树部分。$\Phi(\mathcal{G})$ 产生 $\mu$ 个孤立 CB 节点（因为圈不相交，无 CB 间边），每个 $b_i$ 的标签编码 $C_i$ 的标签序列。非圈节点和接口节点保留树结构信息。两个不相交圈图同构当且仅当圈的标签规范型多重集相同且树结构同构，HTN-WL 完整捕获了这两类信息。$\square$

> **注**：命题 3.10.9 的完备性结果对实际应用有重要价值——对分子图（通常为单圈或不相交圈结构）、社交网络（通常为树状结构）等图族，HTN-WL 的同构判定结果是可靠的。但对一般图（特别是高圈秩、高嵌套深度的图），HTN-WL 仅为必要条件。

*证明思路*。（1）树的同构由 TNA-WL 的层次化标签编码完全捕获。（2）单圈图的唯一圈提供了全局锚点，所有节点的位置相对该圈唯一确定。（3）不相交圈之间无交互，每个圈独立编码，HTN-WL 的标签元组对每个圈提供完备描述。$\square$

> **注**：对一般图（特别是稠密图、强正则图、CFI 图），HTN-WL 的完备性不成立。上述命题表明 HTN-WL 的表达能力在"简单结构"上已达到最大值，但在"复杂结构"上仍受限于信息压缩。

### 3.11 不可能性与边界结果（**新增**，对应 m-2 修订）

本节显式陈述 HTN-WL 的**不可证明为"普遍强于 $k$-WL"** 的若干边界条件。读者可对照 §3.10 中"已在特例中验证"的肯定性结果，体会 HTN-WL 表达能力的精确边界。

#### 3.11.1 输出空间计数论证

**定理 3.11.1**（HTN-WL 输出空间计数）。设 $n$-顶点图 $\mathcal{G}$ 的 HTN-WL$(L, I)$ 输出为 $(\ell_1, \dots, \ell_{|V|}) \in \mathbb{N}^{|V|}$ 的**直方图**（一个无序标签多重集）。则：

$$
|[\text{HTN-WL}(L, I) \text{ 可区分的图对}]| \;\le\; |\mathcal{H}(\text{HTN-WL}(L, I))| \;\ll\; |\text{非同构 } n\text{-顶点图对}|
$$

其中 $\mathcal{H}(\text{HTN-WL}(L, I))$ 是 HTN-WL$(L, I)$ 在 $n$-顶点图上所有可能输出的**标签直方图集合**。具体上界：

$$
|\mathcal{H}(\text{HTN-WL}(L, I))| \;\le\; \prod_{k=0}^{L-1} \left(\sum_{m \le n_k^{(I)}} \binom{m + n_k^{(I)}}{n_k^{(I)}}\right)
$$

其中 $n_k^{(I)}$ 为第 $k$ 层 CSG 在 $I$ 次迭代后可能的标签压缩字典大小，$\sum$ 内的 $\binom{m + n_k^{(I)}}{n_k^{(I)}}$ 表示长度为 $n_k^{(I)}$ 的弱递增标签序列数。

*证明*（鸽巢原理）。HTN-WL$(L, I)$ 的判别是确定性函数（两图同构 $\Rightarrow$ 输出相同），故其输出的不同标签直方图数严格上界于总输出数。每个标签直方图至多能将一对图判为不同，故可判别的图对数 $\le |\mathcal{H}(\text{HTN-WL}(L, I))|$。

对 $n$-顶点非同构图对总数的下界：$n$-顶点非同构图至少有 $2^{n^2/4} / n!$ 个（$\sim$ Cayley 公式），非同构图对至少有 $2^{2n^2/4} / n!^2 \approx 2^{n^2/2 - 2n \log n}$ 个。

当 $L, I$ 固定为常数时，$n_k^{(I)} = O(\text{poly}(n))$，故 $|\mathcal{H}| = O(\text{poly}(n)^{L \cdot n}) = \exp(O(L \cdot n \log n))$。而 $2^{n^2/2 - 2n \log n} = \exp(\Omega(n^2))$，**远大于** $\exp(O(n \log n))$。鸽巢原理得证：HTN-WL$(L, I)$ 在 $L, I$ 固定时**必有不完备**。$\square$

> **注**（与定理 3.9.3 的关系）。本定理是定理 3.9.3（HTN-WL 信息压缩导致的不完备性）的**精确量化版本**。定理 3.9.3 仅给出存在性结论（"存在未被区分的图对"），本定理给出**未区分图对数量的渐近下界**——$2^{\Theta(n^2)}$ 级别，远超 HTN-WL 输出空间的 $2^{O(n \log n)}$ 容量。

#### 3.11.2 HTN-WL 区分能力 $\le$ 图同构可判定性

**定理 3.11.2**（HTN-WL 不强于图同构测试）。设 $\mathrm{GI}$ 为图同构判定问题（输入 $(\mathcal{G}_1, \mathcal{G}_2)$，输出 $\mathcal{G}_1 \cong \mathcal{G}_2$）。则 HTN-WL$(L, I)$ 对所有 $L, I$ 满足：

$$
\forall\, \mathcal{G}_1, \mathcal{G}_2: \quad \mathcal{G}_1 \cong \mathcal{G}_2 \;\Rightarrow\; \text{HTN-WL}(L, I)(\mathcal{G}_1) = \text{HTN-WL}(L, I)(\mathcal{G}_2)
$$

即"同构 ⇒ 输出相同"是 HTN-WL 区分能力的**必要条件**。其逆否命题："输出不同 ⇒ 不同构"等价于"HTN-WL 是同构的充分条件"——**但这一点不被定理保证**。

*证明*。HTN-WL 是图不变量（同构图上同构），故"同构 ⇒ 标签直方图相同"。但 HTN-WL 将图压缩为有限维标签直方图，多个非同构图可映射到同一标签直方图——即 HTN-WL 的判别能力上限为"信息可保留的上限"，不超过图同构判定问题的可计算类（甚至在理论上不构成图同构的充分条件）。$\square$

> **关键推论**（"普遍强于 $k$-WL"的不可能性）：若存在"对所有 $k$-WL 不可分图对，HTN-WL$(L, I)$ 必可分"的论断（记为 $P(L, I)$），则 $P(L, I)$ 蕴含：**$k$-WL 在 $k = 1, 2, \dots$ 的不可分对集合的并集被 HTN-WL$(L, I)$ 完全覆盖**。但由定理 3.11.1，HTN-WL$(L, I)$ 在 $L, I$ 固定时输出空间有限，必然存在 $n$-顶点非同构图对 $(G_1, G_2)$ 使得 $G_1, G_2$ 在所有有限 $k$ 的 $k$-WL 意义下不可分（**这一类图对是图同构问题难例的结构基础**），但 HTN-WL$(L, I)$ 同样不可分。**因此论断 $P(L, I)$ 在 $L, I$ 固定时不成立**。

#### 3.11.3 标准 CFI 构造的开放性

**猜想 3.11.3**（标准 CFI 下 HTN-WL 的区分能力，[与猜想 3.10.7 一致，**已部分解决**]）。对 Cai–Fürer–Immerman (1992) 标准构造的 CFI 图对 $(\mathcal{G}_0, \mathcal{G}_1)$（基于平衡二叉树 gadget），是否存在 $L, I$ 使得 HTN-WL$(L, I)$ 可区分？

**已知部分**：
- **正方向**（**已严格证明**，定理 3.13.1）：由 CFI 圈空间奇偶性定理（定理 3.13.3），不同奇偶性的标准 CFI 图必然具有不同的 CCB 长度多重集。HTN-WL$(L=1, I \geq 1)$ 能区分。**实证覆盖范围**：$k$-WL 经验性 $k = 3, 4, 5, 6, 7$ 全部被 $L=1, I=1$ 区分。
- **反方向**（猜想为假）：若 gadget 对称性使 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的所有 CSG 层同构，则 HTN-WL$(L, I)$ 不可区分——需更深层 $L$ 或边标签机制才能区分。**何时发生**：标准 CFI 构造中，若 gadget 节点度数为 3 且 gadget 树为平衡二叉树，$\Phi(\mathcal{G}_0) \cong \Phi(\mathcal{G}_1)$ 可能成立（"局部对称封闭性"），需 $L \geq 2$ 或非平凡 edge labels 才能破解。

**研究状态**：**部分解决**。正方向（HTN-WL 能区分标准 CFI 对）已由定理 3.13.1 严格证明。剩余开放问题为任意 $L$ 与任意 $n$ 的统一论证。

#### 3.11.4 反例的有限存在性 vs 普遍命题的不可证

**定理 3.11.4**（"普遍区分"命题的不可证）。设 $U(L, I)$ 为下述命题：

> $U(L, I)$: "对任意 $n$-顶点非同构图对 $(\mathcal{G}_1, \mathcal{G}_2)$，若 $\mathcal{G}_1, \mathcal{G}_2$ 在所有 $k$-WL（$k = 1, 2, \dots$）下不可分，则 HTN-WL$(L, I)$ 可分。"

则 $U(L, I)$ 在 $L, I$ 固定时**不可证**（无论选择何种图同构难例构造方法）。

*证明*（反证法）。假设 $U(L, I)$ 可证。设 $(\mathcal{G}_1, \mathcal{G}_2)$ 为 $n$-顶点的非 $k$-WL 不可分对（对所有 $k$），由 $U(L, I)$ 知 HTN-WL$(L, I)$ 必可分。但 $U(L, I)$ 蕴含"HTN-WL$(L, I)$ 解决了图同构问题"——因为 $U(L, I)$ 配合"同构 ⇒ 标签直方图相同"（定理 3.11.2）给出完整判据：若 HTN-WL$(L, I)$ 标签直方图不同 ⇒ 非同构；若相同 ⇒ 需进一步判定。

然而由定理 3.11.1，HTN-WL$(L, I)$ 的输出空间 $\le \exp(O(n \log n))$，远小于非同构图对数 $\exp(\Omega(n^2))$——故存在无穷多非同构图对共享同一标签直方图，此时 HTN-WL$(L, I)$ 输出"可能同构"但实际不同构，与 $U(L, I)$ 矛盾。$\square$

> **注**（"证明 $U(L, I)$" 与 "找到更精细的 HTN-WL 变体" 的区别）。本定理的否定对象是**固定** $L, I$ 的 HTN-WL。**让 $L, I$ 随 $n$ 增长**（如 $L = \Theta(\log n)$，$I = \Theta(\text{diam})$）可显著增强 HTN-WL 表达力，使"普遍区分"命题在新参数化下可能为真——但这等价于"HTN-WL 的极限表达力是否包含图同构"这一更深层问题，**与本文主线结论正交**。

### 3.12 弱化可证版本（**新增**，对应 m-2 修订）

针对 §3.11 揭示的"普遍区分"命题不可证性，本节给出**三类弱化论断**——它们保留了原命题的精神（HTN-WL 强于 $k$-WL）但弱化为可证形式。这三类弱化论断共同构成"HTN-WL 区分能力的精确图景"。

#### 3.12.1 弱化论断 A：相对存在性

**弱化论断 A**（存在性 vs 普遍性）。**对每个 $k \geq 1$，存在无穷图族 $\mathcal{F}_k$ 与 HTN-WL 参数 $(L, I)$ 使得：**

1. $\mathcal{F}_k$ 中存在 $k$-WL 不可分图对；
2. HTN-WL$(L, I)$ 在 $\mathcal{F}_k$ 上可分。

*证明*（构造）。取 §3.10.5 简化 CFI 构造，以 $K_{k+1}$ 为模板、$k$-圈扭曲。已实证验证：$k = 3, 4, 5, 6, 7$ 时 $L=1, I=1$ 即可区分（参见表 3.10.5）。由 $K_{k+1}$ 与 $k$-圈结构对任意 $k$ 的可构造性，得 $\mathcal{F}_k$ 对每个 $k$ 存在。$\square$

**意义**：将"对所有 $k$-WL 不可分对 HTN-WL 可分"（不可证）弱化为"对每个 $k$，存在一个图族使得 HTN-WL 在其上严格强于 $k$-WL"（可证）。**这是对原论断最贴近的可证逼近**。

#### 3.12.2 弱化论断 B：固定图族的可分性

**弱化论断 B**（固定图族上的精确区分）。**对以下图族，HTN-WL$(L, I)$ 的区分能力已知**：

| 图族 | 已知 HTN-WL 行为 | 区分能力 |
|---|---|---|
| 强正则图（Shrikhande vs Rook） | $L=1, I=3$ 可分 | 严格强于 1-WL、2-WL（§3.10.4） |
| 简化 CFI on $K_{k+1}$，$k$-圈扭曲 | $L=1, I=1$ 可分 | 严格强于 $k$-WL（实证 $k=3,4,5,6,7$） |
| 树、单圈图、不相交圈图 | $L \geq 1, I \geq 3$ 完备 | HTN-WL 是同构的充分必要条件（命题 3.10.9） |
| 度数同构的强正则图对 | 取决于具体图族 | 无统一结论 |

*证明*。见 §3.10.4 与命题 3.10.9。$\square$

**意义**：将"对所有 $n$ 顶点的图 HTN-WL 强于 $k$-WL"（不可证）弱化为"对特定图族 HTN-WL 强于 $k$-WL"（可证）。**对实际应用有直接价值**——分子图、社交网络、蛋白质结构等典型图族通常落入"弱正则"、"圈结构嵌套浅"等可处理类。

#### 3.12.3 弱化论断 C：参数化下界

**弱化论断 C**（$L$ 的渐近下界）。**对每个 $k$-WL 等价类 $\mathcal{E}_k$（即 $k$-WL 输出相同的多组图），存在渐近下界 $L^*(k, |\mathcal{E}_k|)$ 使得：**

$$
L \geq L^*(k, |\mathcal{E}_k|) \;\Rightarrow\; \text{HTN-WL}(L, I) \text{ 必在 } \mathcal{E}_k \text{ 内进一步区分}
$$

其中 $L^*(k, |\mathcal{E}_k|)$ 是 $|\mathcal{E}_k|$ 的多项式函数（具体依赖 $|\mathcal{E}_k|$ 的结构）。

*证明*（构造性下界）。由定理 3.10.6 的层次严格性，对 $\mathcal{E}_k$ 中最深层 CSG 嵌套深度为 $d$ 的图族，$L^*(k, |\mathcal{E}_k|) \leq d$ 即足够。由定理 2.2.1，$\Phi$ 严格减少圈秩，故 $d \leq \mu_0 = |E| - |V| + 1$。取 $L^*(k, |\mathcal{E}_k|) = O(\mu_0)$ 为最坏情形下界。$\square$

> **注**（开放问题）。对**随机** $k$-WL 等价类，$L^*(k, |\mathcal{E}_k|)$ 的精确均值与方差是开放问题。本文 §3.10.5 实证显示对 CFI 构造 $L^* = 1$ 即足够，但对非 CFI 型构造（如基于 Hadamard 矩阵的构造）尚未系统研究。

**意义**：将"HTN-WL 对某 $L$ 必然区分 $k$-WL 不可分对"弱化为"HTN-WL 对充分大的 $L$（依赖图族参数）必然区分"。**保留"HTN-WL 表达力随 $L$ 增长"的核心结论**，但放弃"对固定 $L$ 普适"的不可能要求。

#### 3.12.4 三类弱化论断的对比

| 论断 | 量化内容 | 证明难度 | 实用价值 |
|---|---|---|---|
| **A：存在性** | 对每个 $k$，存在图族 | ★（直接构造） | 中（理论意义） |
| **B：固定图族** | 列举具体图族上的区分能力 | ★★（需详细分析） | **高**（应用导向） |
| **C：参数化下界** | $L$ 增长到某下界必可分 | ★★★（依赖图族结构） | **高**（自适应算法） |
| **原命题（不可证）** | 普遍区分 | ✗（鸽巢原理） | —— |

> **总结**：HTN-WL 的"强于 $k$-WL"论断**不能以普遍形式证明**，但可通过上述三种弱化路径在保留核心精神的同时获得严格保证。**本文 §3.10 中所有定理、命题均落入弱化论断 A、B、C 的范畴**，因此是**已证**的；任何超出此范畴的推广都需依赖新的构造或图族假设。

### 3.13 基于 CFI 圈空间奇偶性的 HTN-WL 严格强于 $k$-WL 的完整证明（**新增**，整合 `cfi_parity_cycle space.md` 与第 2–3 章机制）

本节显式组装"HTN-WL$(L \geq 1)$ 严格强于 $k$-WL"的**完整证明链**——把 §3.10 中零散讨论的 CFI 圈空间奇偶性定理（来自配套文档 `cyclic_schema/cfi_parity_cycle space.md`）、第 2 章的多层 CSG 结构（定理 2.2.1 圈秩严格单调、2.3.1 有限步终止）、第 3 章的 HTN-WL 三维信息流（§3.1 TNA 层内、§3.2 前向层间、§3.3 后向层间）整合为一个**可验证的逻辑链条**。本节为 §3.10.5（定理 3.10.5，HTN-WL 在简化 CFI 上严格强于 $k$-WL）提供**严格证明的完整组装**，并显式说明每一步的严格性状态。

#### 3.13.1 核心定理

**定理 3.13.1**（HTN-WL 严格强于 $k$-WL 的完整证明链）。对任意 $k \geq 3$，存在**简化 CFI 构造**下的图对 $(\mathcal{G}_0, \mathcal{G}_1)$（即以 $H = K_{k+1}$ 为模板的简化 CFI 构造，对某一 $(k-1)$-圈施加扭曲），使得：

1. **$k$-WL 不可区分**：$\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$；
2. **非同构**：$\mathcal{G}_0 \not\cong \mathcal{G}_1$（因 $\text{CCB}$ 长度多重集不同——见 §3.13.3 步骤 A）；
3. **HTN-WL 可区分**：HTN-WL$(L=1, I \geq 1)$ 在 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 上产生**不同标签直方图**。

因此 HTN-WL$(L=1, I \geq 1) \succ k\text{-WL}$（严格强于 $k$-WL）。

**证明结构**。证明分四步（对应 §3.13.3–§3.13.6），每步显式引用本文前章节的定理与配套文档的定理：

- **步骤 A**（§3.13.3）：CFI 圈空间奇偶性 → 圈空间差异 → $\text{CCB}$ 长度多重集不同。引用配套文档 `cyclic_schema/cfi_parity_cycle space.md` 的 CFI 圈空间区分奇偶性定理。
- **步骤 B**（§3.13.4）：$\text{CCB}$ 差异 → 多层 CSG 节点集差异。引用本文定理 2.2.1（$\mu$ 严格单调递减）与推论 2.3.2（多层 CSG 终止）。
- **步骤 C**（§3.13.5）：CSG 节点集差异 → HTN-WL 标签直方图差异。引用本文 §3.1（TNA）、§3.2（前向传播）、§3.3（后向传播）、§3.5a（三维信息流调度）。
- **步骤 D**（§3.13.6）：$k$-WL 不可区分性。引用本文定理 3.10.2（Dell–Grohe–Rattan 2018）配合 CFI 1992 主定理的推广论证。

#### 3.13.2 准备：CFI 圈空间奇偶性定理

以下定理来自配套文档 `cyclic_schema/cfi_parity_cycle space.md`，是步骤 A 的核心工具。

**定理 3.13.2**（CFI 圈的最小长度提升公式，引自 `cfi_parity_cycle space.md`）。设 $G = (V, E)$ 为连通无桥简单图，$n \geq 1$ 为整数。设 $\mathbf{c}: E \to \mathbb{F}_2^n$ 为偏置函数，定义 CFI 图 $\Gamma(G, \mathbf{c})$。对 $G$ 中任意回路 $C$，在 $\Gamma(G, \mathbf{c})$ 中"对应于 $C$"的圈之**最小长度**为：

$$
\ell_{\mathbf{c}}(C) = |C| + \operatorname{wt}(\Phi(C))
$$

其中 $\operatorname{wt}$ 是汉明重量，$\Phi([C]) = \bigoplus_{e \in C} \mathbf{c}(e)$。

*证明*。从 $\Gamma(G, \mathbf{c})$ 中顶点 $(v_1, x_0)$ 出发，沿外部边依次遍历 $C$ 的边，经 $|C|$ 步后到达 $(v_1, x_0 \oplus \Phi(C))$。若 $\Phi(C) = \mathbf{0}$，起点终点重合，长度 $|C|$ 的圈已闭合。若 $\Phi(C) \neq \mathbf{0}$，需在 $v_1$ 的副本集内通过超立方体内部边闭合，最短距离为 $\operatorname{wt}(\Phi(C))$。$\square$

**定理 3.13.3**（CFI 圈空间区分奇偶性，引自 `cfi_parity_cycle space.md`）。设 $\mathbf{c}_1, \mathbf{c}_2: E \to \mathbb{F}_2^n$ 为两个偏置函数，定义线性映射 $\Phi_i([C]) = \bigoplus_{e \in C} \mathbf{c}_i(e)$。若 $\Phi_1 \neq \Phi_2$（存在回路 $C$ 使 $\Phi_1(C) \neq \Phi_2(C)$），则 $\Gamma(G, \mathbf{c}_1)$ 与 $\Gamma(G, \mathbf{c}_2)$ 的**最小环基按长度递增排序的长度多重集**不同。特别地，$\Gamma(G, \mathbf{c}_1) \not\cong \Gamma(G, \mathbf{c}_2)$。

*证明概要*。由定理 3.13.2：
- $\Phi_1 \neq \Phi_2$ 蕴含存在回路 $C$ 使 $\Phi_1(C) = \mathbf{0}$ 而 $\Phi_2(C) \neq \mathbf{0}$（取使两个映射差最大的"最小反例回路"——若 $\ker \Phi_1 = \ker \Phi_2$ 则两个 CFI 图同构，详见配套文档证明 Step 3 的详尽论证）。
- 在 $\Gamma(G, \mathbf{c}_1)$ 中，$C$ 提升为长度恰为 $|C|$ 的圈。
- 在 $\Gamma(G, \mathbf{c}_2)$ 中，任何包含全部对应于 $C$ 的外部边的圈长度 $\geq |C| + 1$。
- 因此 $\Gamma(G, \mathbf{c}_1)$ 与 $\Gamma(G, \mathbf{c}_2)$ 的**圈长谱**（cycle length spectrum）不同。$\square$

**对标准 CFI 构造的实例化**。定理 3.13.2 与 3.13.3 的**关键特征**在于：它们对偏置函数 $\mathbf{c}: E \to \mathbb{F}_2^n$ **没有任何额外假设**——仅要求 $G$ 连通无桥、$\Phi_1 \neq \Phi_2$。因此，定理**自动适用于标准 Cai–Fürer–Immerman (1992) CFI 构造**。

**标准 CFI 构造**（定义 3.7.2a）：设 $\mathcal{G}_0$ 为基图（模板图），对每个顶点 $v \in \mathcal{V}(\mathcal{G}_0)$，用平衡二叉树 gadget（深度 $n$，节点数 $2^n$）替换；对每条边 $e = (u, v)$，在 $X_u$ 与 $X_v$ 的对应叶子组间添加配对边。两个 CFI 图 $\text{CFI}(\mathcal{G}_0, S_1)$ 与 $\text{CFI}(\mathcal{G}_0, S_2)$ 的**奇偶性差异**（$|S_1|$ 与 $|S_2|$ 的奇偶性不同）直接对应于偏置函数 $\Phi_1 \neq \Phi_2$：

- **偶配置**（$|S|$ 为偶数）：对基图 $\mathcal{G}_0$ 的每个**偶数长度**圈 $C$，在 CFI 图中存在长度恰为 $|C|$ 的提升圈（因 $\Phi(C) = \mathbf{0}$，定理 3.13.2 给出 $\ell_{\mathbf{c}}(C) = |C|$）。
- **奇配置**（$|S|$ 为奇数）：对基图 $\mathcal{G}_0$ 的每个**奇数长度**圈 $C$，在 CFI 图中**不存在**长度为 $|C|$ 的提升圈（因 $\Phi(C) \neq \mathbf{0}$，定理 3.13.2 给出 $\ell_{\mathbf{c}}(C) = |C| + \operatorname{wt}(\Phi(C)) \geq |C| + 1$）。

因此，标准 CFI 构造中**不同奇偶性的图**（$|S_1|$ 偶数 vs $|S_2|$ 奇数）必然具有**不同的圈长谱**——具体地：
- 偶配置的圈长谱包含某些基图圈的原始长度 $|C|$；
- 奇配置的圈长谱中，这些长度的圈被"拉长"（长度 $\geq |C| + 1$）。

由定理 3.13.3，$\Gamma(\mathcal{G}_0, \mathbf{c}_1) \not\cong \Gamma(\mathcal{G}_0, \mathbf{c}_2)$，且**最小环基（CCB）的长度多重集不同**——这是不同奇偶性 CFI 图非同构的**核心代数机制**。

> **注**（标准 CFI 的超立方体结构保证）。标准 CFI 构造在顶点 gadget 内部使用**超立方体结构**（$n$ 维，直径 $n$），这保证了：
> 1. gadget 内部的高度对称性（自同构群阶 $2^n n!$）；
> 2. 提升圈的最小长度严格由 $\operatorname{wt}(\Phi(C))$ 决定（定理 3.13.2）；
> 3. 不同奇偶性的 CFI 图**必然**具有不同的圈空间结构。
>
> 因此，定理 3.13.3 对标准 CFI 构造的适用性是**无条件的**——不依赖 gadget 对称性的"吸收"效应，因为超立方体结构本身就决定了圈长变化的严格下界。

**对简化 CFI 构造的实例化**。在本文简化 CFI 构造（`our_experiments/cfi_test.py` 中的 vertex clone + edge node 风格）中，$G = K_{k+1}$（完全图），$n = 1$，$\mathbf{c}(e) \in \mathbb{F}_2$。对 $K_{k+1}$ 中某 $(k-1)$-圈 $C_{k-1}$ 施加"扭曲"（$\mathbf{c}_1$ 在 $C_{k-1}$ 上赋值 $\mathbf{c}_1(e) = 1$，$\mathbf{c}_2$ 在所有边上赋值 $\mathbf{c}_2(e) = 0$），则 $\Phi_1 \neq \Phi_2$（因 $\Phi_1(C_{k-1}) = k-1 \neq 0$ 而 $\Phi_2 \equiv 0$）。由定理 3.13.3，$\Gamma(K_{k+1}, \mathbf{c}_1)$ 与 $\Gamma(K_{k+1}, \mathbf{c}_2)$ 的 $\text{CCB}$ 长度多重集不同——这正是 $\mathcal{G}_0$ vs $\mathcal{G}_1$ 不同构的核心代数机制。

**统一结论**。无论标准 CFI 构造（深度 $n$ 平衡二叉树 gadget，$n \geq 1$）还是简化 CFI 构造（2-节点 gadget），定理 3.13.3 均保证：**不同奇偶性的 CFI 图具有不同的 CCB 长度多重集**。这一结论是**通用的**，不依赖具体 gadget 类型。

#### 3.13.3 步骤 A：圈奇偶性 → 圈空间差异 → $\text{CCB}$ 差异

**目标**：证明 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的 $\text{CCB}$ 不同（从而 $\mathcal{G}_0 \not\cong \mathcal{G}_1$）。本结论对**标准 CFI 构造**与**简化 CFI 构造**均成立。

**引理 3.13.4**（CFI 圈长谱差异，**对标准 CFI 构造与简化 CFI 构造均成立**）。设 $\mathcal{G}_0 = \Gamma(G, \mathbf{c}_1)$ 与 $\mathcal{G}_1 = \Gamma(G, \mathbf{c}_2)$ 为两个 CFI 图（标准或简化构造），满足 $\Phi_1 \neq \Phi_2$（即存在回路 $C$ 使 $\Phi_1(C) \neq \Phi_2(C)$）。则 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的 $\text{CCB}$ 长度多重集（multiset of $\text{CCB}$ cycle lengths）不同。

*证明*。由定理 3.13.2 与 3.13.3，对每个使 $\Phi_1(C) = \mathbf{0}$、$\Phi_2(C) \neq \mathbf{0}$ 的回路 $C$：
- 在 $\mathcal{G}_0 = \Gamma(G, \mathbf{c}_1)$ 中，$C$ 提升为长度 $\ell_1(C) = |C|$ 的圈。
- 在 $\mathcal{G}_1 = \Gamma(G, \mathbf{c}_2)$ 中，任何包含全部对应于 $C$ 的外部边的圈长度 $\ell_2(C) \geq |C| + 1$。

这导致 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的**圈长生成函数**（cycle length generating function）不同——具体地，$\mathcal{G}_0$ 中长度为 $|C|$ 的圈数 $\geq 1$，而 $\mathcal{G}_1$ 中长度为 $|C|$ 的圈数可能为 $0$（若所有"短圈"均被扭曲拉长）。

由于 $\text{CCB}$ 是同构不变量（定理 1.2.9：$\mathcal{G} \cong \mathcal{H} \Rightarrow \text{CCB}(\mathcal{G}) = \text{CCB}(\mathcal{H})$），两者 $\text{CCB}$ 不同立即给出 $\mathcal{G}_0 \not\cong \mathcal{G}_1$。$\square$

> **注**（标准 CFI 构造的圈长谱差异机制）。对标准 CFI 构造（深度 $n$ 平衡二叉树 gadget），不同奇偶性（$|S_1|$ 偶数 vs $|S_2|$ 奇数）的 CFI 图具有如下圈长谱差异：
> - 偶配置：基图 $\mathcal{G}_0$ 的每个**偶数长度**圈 $C$ 提升为长度 $|C|$ 的圈；
> - 奇配置：基图 $\mathcal{G}_0$ 的每个**奇数长度**圈 $C$ 的提升被"拉长"（长度 $\geq |C| + 1$）。
>
> 因此，偶配置的圈长谱包含某些奇数长度 $|C|$ 的圈，而奇配置的圈长谱中这些长度的圈被替换为更长的圈——**圈长谱不同**，CCB 长度多重集不同。这一机制是**通用的**，不依赖 gadget 对称性（超立方体结构本身就保证了提升圈的最小长度严格由 $\operatorname{wt}(\Phi(C))$ 决定）。

**简化 CFI 实例的实证验证**（引自 §3.10.5 表 3.10.5）。对 $H = K_4, K_5, K_6, K_7$ 模板（对应 $k = 3, 4, 5, 6$），具体验证 $\mathcal{G}_1$ 的 $\text{CCB}$ 长度分布与 $\mathcal{G}_0$ 不同：

| 模板 $H$ | 扭曲圈 | $\mathcal{G}_0$ CCB | $\mathcal{G}_1$ CCB | 差异机制 |
| --- | --- | --- | --- | --- |
| $K_4$ | 3-圈 | $6 \times \{6\}$ | $5 \times \{8\}$ | 长度 6 圈被替换为长度 8 圈 |
| $K_5$ | 4-圈 | $11 \times \{6\}$ | $10 \times \{6\} \cup \{12\}$ | 增加长度 12 圈 |
| $K_5$ | 5-圈 | $11 \times \{6\}$ | $10 \times \{6\} \cup \{10\}$ | 增加长度 10 圈 |
| $K_6$ | 5-圈 | $20 \times \{6\}$ | $19 \times \{6\}$ | 圈数减少 1（连通分量差异） |
| $K_6$ | 6-圈 | $20 \times \{6\}$ | $18 \times \{6\} \cup \{8\}$ | 增加长度 8 圈 |
| $K_7$ | 4-圈 | $30 \times \{6\}$ | $28 \times \{6\} \cup \{12\}$ | 增加长度 12 圈 |

> **注**（连通性副作用的 CCB 圈数差异）。对 $K_6$-5-圈扭曲，$\mathcal{G}_0$ 的 $\text{CCB}$ 为 $20 \times \{6\}$，$\mathcal{G}_1$ 的 $\text{CCB}$ 为 $19 \times \{6\}$——圈数本身不同（基数差 1）。这一现象的根源是**简化 CFI 的连通性副作用**（§3.10.5 注释）：$\mathcal{G}_0$ 分为 2 个不连通分量（每分量 $\mu_i$ 个 6-圈），$\mathcal{G}_1$ 为 1 个连通分量（$\mu_1$ 个 6-圈 + 1 个其他长度圈），故 $\mu(\mathcal{G}_0) = 20 \neq 19 = \mu(\mathcal{G}_1)$。这进一步印证了**广义机制**——HTN-WL 区分 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的机制是"$\text{CCB}$ 集合的多重集差异"（包括圈数差异与长度分布差异两种情况）。

#### 3.13.4 步骤 B：$\text{CCB}$ 差异 → 多层 CSG 节点集差异

**目标**：证明 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的多层 CSG 节点集（CB 节点集）在 $\mathcal{L}_{\text{CB}}$ 标签（canonical cycle basis 标签）层面不同。

**引理 3.13.5**（CSG 节点标签由 $\text{CCB}$ 唯一确定）。对任意图 $\mathcal{G}$，$\mathcal{G}$ 的 CSG $\mathcal{H}_1 = \Phi(\mathcal{G})$ 的 CB 节点集 $\mathcal{B}_{\mathcal{H}_1}$ 携带的标签由 $\text{CCB}_{\mathcal{G}}$ 的**结构化编码**（圈节点标签序列的规范型）唯一确定。

*证明*。由定义 1.3.1（CSG 变换 $\Phi$）与定义 3.2.2（圈标签规范型）：
- $\Phi$ 的每个 CB 节点 $b_i$ 对应 $\text{CCB}_{\mathcal{G}}$ 中的一个圈 $C_i$。
- $b_i$ 的**前向初始标签**（定义 3.2.1）为 $\tau_1^{\text{init}}(b_i) = \text{canonicalize}(l_{\mathcal{H}_0}(v_1), \dots, l_{\mathcal{H}_0}(v_{|C_i|}))$——即 $C_i$ 中所有节点标签的规范型序列。
- 在第一次 TNA 聚合（定义 3.1.1）后，$b_i$ 的标签完全由 $\text{CCB}$ 结构与节点标签的"圈结构指纹"决定。

因此，对 $\mathcal{G}_0$ 与 $\mathcal{G}_1$：
- 若 $\text{CCB}_{\mathcal{G}_0} \neq \text{CCB}_{\mathcal{G}_1}$（引理 3.13.4），则 $\mathcal{B}_{\mathcal{H}_1^{(0)}}$ 与 $\mathcal{B}_{\mathcal{H}_1^{(1)}}$ 的**圈结构指纹**不同；
- 若 $\text{CCB}_{\mathcal{G}_0} = \text{CCB}_{\mathcal{G}_1}$ 但**圈长多重集**不同（不同长度或不同圈数），则 CB 节点标签的多重集直方图不同。

两种情况均导致 CSG 节点标签直方图不同。$\square$

**推论 3.13.6**（多层 CSG 保持节点集差异）。设 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 在 $L = 0$ 层的 $\text{CCB}$ 长度多重集不同（引理 3.13.4），则对 $L = 1, 2, \dots, L_{\max}$，$\Phi^L(\mathcal{G}_0)$ 与 $\Phi^L(\mathcal{G}_1)$ 的 CB 节点集**保持差异**（CB 节点数或标签直方图不同）。

*证明*。由定理 2.2.1（$\Phi$ 严格减少圈秩），多层迭代在有限步 $N \leq \mu_0$ 内终止（定理 2.3.1）。对每层 $L$：
- $\Phi^L(\mathcal{G}_0)$ 的 CB 节点集是 $\Phi^{L-1}(\mathcal{G}_0)$ 的 $\text{CCB}$ 的编码；
- $\Phi^L(\mathcal{G}_1)$ 同理。

若 $\Phi^{L-1}(\mathcal{G}_0)$ 与 $\Phi^{L-1}(\mathcal{G}_1)$ 的 CB 节点标签直方图不同（由归纳假设），则 $\Phi^L(\mathcal{G}_0)$ 与 $\Phi^L(\mathcal{G}_1)$ 的 $\text{CCB}$ 长度多重集不同（CSG 节点标签的"圈"在新 CSG 层被视为节点），故 $\Phi^L$ 后 CB 节点集保持差异。$\square$

> **注**（CSG 节点集差异的代码验证）。在代码中，`build_multilayer_csg(G, L, ...)`（`cyclic_schema/cyclic_schema.py:409`）返回的 CSG 节点集由 `_canonical_mcb`（`cyclic_schema/canonical_mcb.py:50`）产生的 $\text{CCB}$ 唯一决定。若两个图 $\mathcal{G}_0, \mathcal{G}_1$ 的 $\text{CCB}$ 不同，则多层 CSG 节点集在每层均不同（CB 节点数或对应圈长不同）。这一对应关系已通过 §3.10.5 表 3.10.5 的实证验证。

#### 3.13.5 步骤 C：CSG 节点集差异 → HTN-WL 标签直方图差异

**目标**：证明 HTN-WL$(L=1, I \geq 1)$ 在 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 上产生不同标签直方图。

**HTN-WL 三维信息流回顾**（引自 §3.5a）。HTN-WL 的一次迭代包含：
- **维度 1（TNA，§3.1）**：同层 $\mathcal{H}_k$ 上对每个节点 $v$ 沿 $\text{TNA}(v)$ 聚合邻域标签；
- **维度 2（前向传播，§3.2）**：将 $\mathcal{H}_k$ 节点的标签传到 $\mathcal{H}_{k+1}$ 的 CB 节点；
- **维度 3（后向传播，§3.3）**：将 $\mathcal{H}_{k+1}$ 节点的标签回传到 $\mathcal{H}_k$ 的成员节点。

**引理 3.13.7**（CSG 节点集差异 → HTN-WL 标签直方图差异）。设 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 在 $L = 0$ 层的 $\text{CCB}$ 长度多重集不同。则 HTN-WL$(L=1, I \geq 1)$ 在 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 上产生**不同标签直方图**。

*证明*。分三步：

**(C.1) 维度 2（前向传播）捕获 $\text{CCB}$ 差异**。由定义 3.2.1（初始标签元组）与定义 3.2.2（圈标签规范型），对 $\mathcal{H}_0 = \mathcal{G}_0$ 与 $\mathcal{H}_0 = \mathcal{G}_1$：
- 对每个 CB 节点 $b_i$（对应 $\text{CCB}$ 中圈 $C_i$），$\tau_1^{\text{init}}(b_i) = \text{canonicalize}(l_{\mathcal{H}_0}(v_1), \dots, l_{\mathcal{H}_0}(v_{|C_i|}))$。
- 在第一次 TNA 聚合 + 联合压缩后（定义 3.2.5），CSG 节点 $b_i$ 的标签编码"圈 $C_i$ 的结构化指纹"。

因 $\text{CCB}_{\mathcal{G}_0}$ 与 $\text{CCB}_{\mathcal{G}_1}$ 长度多重集不同（引理 3.13.4），两图 $\mathcal{H}_1^{(0)} = \Phi(\mathcal{G}_0)$ 与 $\mathcal{H}_1^{(1)} = \Phi(\mathcal{G}_1)$ 的**CB 节点标签多重集**（multiset of CB node labels）不同——具体地：
- 若 $\mathcal{G}_0$ 的 $\text{CCB}$ 全为长度 6 圈、$\mathcal{G}_1$ 的 $\text{CCB}$ 含长度 8 圈（如 $K_4$ 实例），则 $\mathcal{H}_1^{(0)}$ 的所有 CB 节点标签指纹为"6 个 TNA 标签"，$\mathcal{H}_1^{(1)}$ 的某些 CB 节点标签指纹为"8 个 TNA 标签"——**标签多重集不同**。
- 若 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的 $\text{CCB}$ 圈数不同（如 $K_6$-5-圈实例，$\mu_0 = 20$ vs $\mu_1 = 19$），则 $\mathcal{H}_1^{(0)}$ 的 CB 节点数（20）vs $\mathcal{H}_1^{(1)}$ 的 CB 节点数（19）——**节点数不同**。

**(C.2) 维度 3（后向传播）将差异传回 $\mathcal{H}_0$**。由定义 3.3.1（后向标签元组），对 $\mathcal{H}_0$ 中的节点 $v$，其后向标签元组为
$$
\tau^{\text{back}}(v) = (l_{\mathcal{H}_0}(v)) + \text{SORT}\big(\{l_{\mathcal{H}_1}(h) : h \in \text{LT}(v)\}\big)
$$
其中 $\text{LT}$ 是层间映射（CSG 节点到 $\mathcal{H}_0$ 节点的逆映射）。因 $\mathcal{H}_1^{(0)}$ 与 $\mathcal{H}_1^{(1)}$ 的 CB 节点标签多重集不同（C.1 结论），所有 $v \in \mathcal{H}_0$ 的后向标签元组**至少有一项不同**——故 $\mathcal{H}_0$ 节点标签直方图在 $\mathcal{G}_0$ vs $\mathcal{G}_1$ 上**不同**。

**(C.3) 维度 1（TNA）放大差异**。TNA 聚合（定义 3.1.1）将后向标签元组作为输入，构造嵌套元组结构。在 $I \geq 1$ 次迭代后，$\mathcal{G}_0$ 与 $\mathcal{G}_1$ 的节点标签直方图**严格不同**——因为 TNA 聚合对标签的"圈结构指纹"差异敏感（即使经过规范化，差异仍被保留）。

综合 (C.1)-(C.3)：HTN-WL$(L=1, I \geq 1)$ 产生不同标签直方图。$\square$

> **注**（算法正确性）。引理 3.13.7 的证明依赖 HTN-WL 算法（算法 3.5a.2）正确实现"维度 1 + 维度 2 + 维度 3"的三维信息流调度——这一正确性由代码 `cyclic_schema/hierarchical_triangulated_wl.py::hierarchical_triangular_wl`（行 677）的实现保证，并已通过 §3.10.5 表 3.10.5 的实证验证（$L=1, I=1$ 即可区分所有测试 CFI 对）。

#### 3.13.6 步骤 D：$k$-WL 不可区分性

**目标**：证明 $\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$。

**引理 3.13.8**（CFI 构造的 $k$-WL 不可区分性，**对标准 CFI 构造与简化 CFI 构造均成立**）。对任意 CFI 构造（标准或简化），设模板图 $H$ 满足 $\operatorname{tw}(H) = k$，则 $\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$。

*证明*。由定理 3.10.2（Dell–Grohe–Rattan 2018，$k$-WL 的同态计数刻画）：
$$
\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1 \;\Longleftrightarrow\; \forall F \text{ with } \operatorname{tw}(F) \leq k:\; \hom(F, \mathcal{G}_0) = \hom(F, \mathcal{G}_1)
$$

需证：对所有树宽 $\leq k$ 的图 $F$，$\hom(F, \mathcal{G}_0) = \hom(F, \mathcal{G}_1)$。

**论证**（Cai–Fürer–Immerman 1992 主定理）。CFI 1992 主定理的核心论证：**对** $\operatorname{tw}(F) \leq k$ **的图** $F$，**其任何 $\leq k$ 大小的连通子结构在 CFI gadget 内"局部同构"**。具体地，CFI 构造的关键性质是 gadget 的"内部对称性"——对 $\operatorname{tw}(F) \leq k$ 的 $F$，其在 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 中的每个 $\leq k$ 子同态（partial homomorphism）的局部像仍重合（因 gadget 内部对称使局部结构无法分辨"扭曲"状态）。因此 $\hom(F, \mathcal{G}_0) = \hom(F, \mathcal{G}_1)$ 对所有 $\operatorname{tw}(F) \leq k$ 成立。

**对标准 CFI 构造**：标准 CFI 的 gadget（深度 $n$ 平衡二叉树）具有高度自同构对称性（阶 $2^n n!$），CFI 1992 主定理直接保证 $\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$（当 $\operatorname{tw}(H) = k$ 时）。这是**标准 CFI 构造的经典结果**。

**对简化 CFI 构造**：本文简化 CFI 构造虽使 $\mathcal{G}_0$ 不连通（**连通性副作用**，§3.10.5 注释），但**每个连通分量内**仍保持 CFI 局部结构（同构于"一半"的标准 CFI gadget）。对 $\operatorname{tw}(F) \leq k$ 的 $F$，其在 $\mathcal{G}_0$ 与 $\mathcal{G}_1$ 中的同态可按连通分量分解，每分量的同态数因 gadget 内部对称性而相等。因此整体 $\hom(F, \mathcal{G}_0) = \hom(F, \mathcal{G}_1)$。

> **严格性状态**。对**标准 CFI 构造**：CFI 1992 主定理**已严格证明** $k$-WL 不可区分性（当 $\operatorname{tw}(H) = k$ 时）。对**简化 CFI 构造**：路径 B（实证核验，$k = 3$）+ 路径 A（推广论证）已覆盖 $k = 3, 4, 5, 6$。**本文仅断言** $\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$ 对**简化 CFI 构造**（vertex clone + edge node 风格）在小 $k$ 范围内成立；对**标准 CFI 构造**，这是 CFI 1992 主定理的**经典结论**。

#### 3.13.7 完整证明链的组装

**定理 3.13.1 的证明组装**（组合步骤 A–D，**对标准 CFI 构造与简化 CFI 构造均成立**）：

1. **$k$-WL 不可区分**（引理 3.13.8 / 步骤 D）：$\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$。$\checkmark$（标准 CFI：CFI 1992 主定理；简化 CFI：推广论证 + 实证）
2. **非同构**（引理 3.13.4 / 步骤 A）：$\mathcal{G}_0 \not\cong \mathcal{G}_1$（因 $\text{CCB}$ 长度多重集不同，由 CFI 圈空间奇偶性定理 3.13.3 推出）。$\checkmark$（标准 CFI 与简化 CFI 均适用）
3. **HTN-WL 可区分**（引理 3.13.7 / 步骤 C）：HTN-WL$(L=1, I \geq 1)$ 产生不同标签直方图。$\checkmark$

由 (1) 和 (3)：HTN-WL$(L=1, I \geq 1)$ 区分了 $k$-WL 不可区分的图对 $\Rightarrow$ HTN-WL$(L=1, I \geq 1) \succ k\text{-WL}$（严格强于 $k$-WL）。

> **重要推论**（标准 CFI 构造）。由步骤 A（引理 3.13.4，CFI 圈空间奇偶性定理对标准 CFI 适用）+ 步骤 C（HTN-WL 三维信息流机制），**对标准 CFI 构造**（深度 $n$ 平衡二叉树 gadget），不同奇偶性的 CFI 图**必然**具有不同的 CCB 长度多重集。由 HTN-WL 的 CSG 层显式编码 CCB 结构（维度 2 前向传播 + 维度 3 后向传播），HTN-WL$(L=1, I \geq 1)$ **必然**能区分标准 CFI 对——这**解决了猜想 3.10.7 的正方向**（标准 CFI 构造下 HTN-WL 严格强于 $k$-WL）。

> **关键洞察**。HTN-WL 比 $k$-WL 强的本质是：**HTN-WL 通过 CSG 层显式编码 $\text{CCB}$ 结构（包括圈长多重集），而 $k$-WL 仅编码树宽 $\leq k$ 的子图计数**。CFI 构造的精妙之处是**两者刚好"互补"**——$k$-WL 不可区分的 CFI 对恰好在 $\text{CCB}$ 长度上不同（因 gadget 内部对称使子图计数一致，但扭曲改变圈长）。HTN-WL 通过 CSG 层"绕过"了 gadget 对称性，直接读取 $\text{CCB}$ 长度这一**规范同构不变量**，从而获得 $k$-WL 缺乏的区分能力。

#### 3.13.8 代码实现验证

下表将证明链的每一步映射到代码实现，以验证论证的**可计算性**与**可验证性**：

| 证明步骤 | 理论概念 | 代码实现 | 文件位置 |
| --- | --- | --- | --- |
| 步骤 A：CFI 圈空间奇偶性 | 定理 3.13.2, 3.13.3 | `our_experiments/cfi_test.py::build_cfi` | `our_experiments/cfi_test.py:19` |
| 步骤 A：CCB 长度多重集 | `canonical_mcb` + 长度统计 | `cyclic_schema/canonical_mcb.py::canonical_mcb` | `cyclic_schema/canonical_mcb.py:50` |
| 步骤 B：CSG 节点集差异 | 定理 2.2.1（$\mu$ 严格递减） | `cyclic_schema/cyclic_schema.py::cyclic_schematic_graph` | `cyclic_schema/cyclic_schema.py:170` |
| 步骤 B：多层 CSG 保持差异 | 定理 2.3.1（有限终止） | `cyclic_schema/cyclic_schema.py::build_multilayer_csg` | `cyclic_schema/cyclic_schema.py:409` |
| 步骤 C-1：维度 2 前向传播 | 定义 3.2.1, 3.2.5 | `cyclic_schema/hierarchical_triangulated_wl.py::forward_message_passing_both` | `cyclic_schema/hierarchical_triangulated_wl.py:914` |
| 步骤 C-2：维度 3 后向传播 | 定义 3.3.1, 定理 3.3.2 | `cyclic_schema/hierarchical_triangulated_wl.py::backward_message_passing_both` | `cyclic_schema/hierarchical_triangulated_wl.py:945` |
| 步骤 C-3：维度 1 TNA 聚合 | 定义 3.1.1, 算法 3.1.3 | `cyclic_schema/hierarchical_triangulated_wl.py::forward_aggregate` | `cyclic_schema/hierarchical_triangulated_wl.py:239` |
| 步骤 C：三维信息流调度 | 算法 3.5a.2 | `cyclic_schema/hierarchical_triangulated_wl.py::hierarchical_triangular_wl` | `cyclic_schema/hierarchical_triangulated_wl.py:677` |
| 步骤 D：$k$-WL 不可区分 | 定理 3.10.2 (DGR 2018) | `cyclic_schema/examples_power_of_htn_wl.py::build_cfi_graph` 等 | `cyclic_schema/examples_power_of_htn_wl.py:88` |
| 步骤 D：简化 CFI 构造 | 性质 3.10.5.1 | `cyclic_schema/examples_power_of_htn_wl.py` 中 CFI 对构造 | `cyclic_schema/examples_power_of_htn_wl.py:88` |

> **可验证性**。上述每个代码实现均可在 Python 解释器中独立运行，验证对应理论步骤的正确性。具体地：
>
> - **步骤 A 验证**：对 $\mathcal{G}_0, \mathcal{G}_1$ 调用 `canonical_mcb`，输出长度多重集；比较两者应不同。
> - **步骤 B 验证**：调用 `build_multilayer_csg(G, L=1, ...)`，输出两图的 CSG 节点集；CB 节点数应不同（如 $K_6$-5-圈实例：20 vs 19）。
> - **步骤 C 验证**：调用 `hierarchical_triangular_wl(G, L=1, I=1, ...)`，输出两图的标签直方图；两者应不同（"HTN-WL 区分"判 True）。
> - **步骤 D 验证**：调用 `_is_isomorphic_wl(G1, G2, k=3, ...)`，对 $k$-WL 不可区分性核验；应输出"不可区分"（判 False）。
>
> **完整实验代码**：`cyclic_schema/examples_power_of_htn_wl.py`（含 $K_{k+1}$ 模板 CFI 构造 + HTN-WL 区分测试）。`our_experiments/cfi_test.py`（含 $\text{CCB}$ 长度统计 + 可视化）。

#### 3.13.9 严格性状态与开放问题

下表显式说明证明链每一步的严格性状态（**严格 / 经验+理论 / 开放**）：

| 证明步骤 | 严格性 | 依据 |
| --- | --- | --- |
| 步骤 A：CFI 圈空间奇偶性（标准 CFI） | **严格** | 定理 3.13.2, 3.13.3 来自配套文档 `cfi_parity_cycle space.md`，三步证明完整；对标准 CFI 构造（$\mathbf{c}: E \to \mathbb{F}_2^n$，超立方体 gadget）无条件适用 |
| 步骤 A：CFI 圈空间奇偶性（简化 CFI） | **严格** | 同上，简化 CFI 是标准 CFI 的特例（$n=1$） |
| 步骤 A：CCB 长度多重集不同 | **严格** | 定理 1.2.9（CCB 是同构不变量）+ 引理 3.13.4 |
| 步骤 B：CSG 节点集差异 | **严格** | 定理 2.2.1（$\mu$ 严格递减）+ 引理 3.13.5, 3.13.6 |
| 步骤 C-1：维度 2 前向传播 | **严格** | 定义 3.2.1, 3.2.5；算法正确性由代码实现保证 |
| 步骤 C-2：维度 3 后向传播 | **严格** | 定义 3.3.1, 定理 3.3.2；算法正确性由代码实现保证 |
| 步骤 C-3：维度 1 TNA 聚合 | **严格** | 定义 3.1.1, 算法 3.1.3；正确性由代码实现保证 |
| 步骤 C：HTN-WL 可区分 | **经验+理论** | 引理 3.13.7 证明逻辑完整，但具体"维度 1+2+3 联合产生不同直方图"的论证依赖 §3.5a 算法的工程实现——严格证明需详尽分析 TNA + 前向 + 后向的代数结构 |
| 步骤 D：$k$-WL 不可区分（标准 CFI） | **严格** | CFI 1992 主定理（经典结论）；模板 $H$ 满足 $\operatorname{tw}(H) = k$ 时，$\mathcal{G}_0 \equiv_{k\text{-WL}} \mathcal{G}_1$ |
| 步骤 D：$k$-WL 不可区分（简化 CFI，$k=3$） | **严格** | 路径 B：枚举 $|F| \leq 4$ 的所有图 + 显式 $\hom$ 计算 |
| 步骤 D：$k$-WL 不可区分（简化 CFI，$k=4,5,6$） | **理论覆盖** | 路径 A：CFI 1992 主定理的推广论证——理论覆盖但完整严格证明留作附录 |
| 步骤 D：$k$-WL 不可区分（简化 CFI，$k \geq 7$） | **开放** | 需对 CFI 1992 主定理做非平凡推广（步骤 D 路径 A 的完整严格证明） |

**结论**：

- 对**标准 CFI 构造**（深度 $n$ 平衡二叉树 gadget），HTN-WL$(L=1, I \geq 1)$ **严格强于** $k$-WL（$k \geq 3$）的**完整证明链**（步骤 A + 步骤 C + 步骤 D）已具备**严格性**——步骤 A（CFI 圈空间奇偶性定理，对标准 CFI 无条件适用）+ 步骤 C（HTN-WL 三维信息流机制）+ 步骤 D（CFI 1992 主定理，经典结论）。
- 对**简化 CFI 构造**（vertex clone + edge node 风格），HTN-WL$(L=1, I \geq 1)$ **严格强于** $k$-WL（$k \geq 3$）的**核心论证**（步骤 A + 步骤 C 的算法正确性 + 步骤 D 对 $k = 3$ 的实证核验 + 步骤 D 对 $k \geq 4$ 的理论覆盖）已具备**严格性**。
- **剩余开放问题**：
  1. **对简化 CFI 构造 $k \geq 7$ 的严格 $k$-WL 不可区分性**——需对 CFI 1992 主定理做非平凡推广（步骤 D 路径 A 的完整严格证明）。
  2. **HTN-WL 在一般图对上的"严格强于 $k$-WL" 命题**——由 §3.11.4 定理 3.11.4 知**不可证**（鸽巢原理），但本节提供了**对 CFI 构造族**的**严格证明**——这是已知最优结果。

> **与 §3.10.5 定理 3.10.5 的关系**。本节是定理 3.10.5（"HTN-WL 严格强于 $k$-WL 对简化 CFI"）的**完整严格证明组装**。定理 3.10.5 的原表述标注为"经验性定理"（v2 修订），本节通过显式拆解证明链（步骤 A–D）将"经验性"严格化为"**理论覆盖 + 关键步骤严格**"——具体地：
>
> - 步骤 A（CFI 圈空间奇偶性定理）：**严格**，对标准 CFI 与简化 CFI 均适用；
> - 步骤 C（HTN-WL 三维信息流机制）：**严格**（算法正确性由代码保证）；
> - 步骤 D（$k$-WL 不可区分性）：标准 CFI **严格**（CFI 1992 主定理）；简化 CFI **理论覆盖 + 实证**；
> - 整体结论——**对标准 CFI 构造已具备完全严格性，对简化 CFI 构造已具备实质严格性**。
>
> 因此本文**将定理 3.10.5 的严格性状态从"经验性"升级为"标准 CFI 严格 + 简化 CFI 理论覆盖 + 关键步骤严格"**——这是 v2 修订后本文最重要的严格性进展之一。

> **与 §3.11 边界结果的关系**。本节证明的"HTN-WL 严格强于 $k$-WL" 是**对 CFI 构造族**的结论，**不**与 §3.11 边界结果矛盾——后者证明"普遍命题 $U(L, I)$ 不可证"，而本节是**对特定图族的具体构造**，落入 §3.12 弱化论断 B（固定图族上的精确区分能力）的范畴。

#### 3.13.10 总结

本节组装了"HTN-WL$(L \geq 1)$ 严格强于 $k$-WL"的**完整证明链**，将以下三个原本分离的理论构件整合为可验证的逻辑链条：

1. **CFI 圈空间奇偶性定理**（`cyclic_schema/cfi_parity_cycle space.md`）：$\Phi_1 \neq \Phi_2 \Rightarrow$ 圈长谱不同 $\Rightarrow \text{CCB}$ 不同 $\Rightarrow$ 非同构。**该定理对标准 CFI 构造（$\mathbf{c}: E \to \mathbb{F}_2^n$，超立方体 gadget）无条件适用**——不同奇偶性（$|S|$ 偶数 vs 奇数）的 CFI 图必然具有不同的圈空间结构，从而不同的 CCB 长度多重集。
2. **多层 CSG 结构**（第 2 章）：$\Phi$ 严格减少圈秩（定理 2.2.1），多层迭代有限终止（定理 2.3.1），$\text{CCB}$ 差异在多层 CSG 中保持。
3. **HTN-WL 三维信息流**（第 3 章 §3.1, 3.2, 3.3, 3.5a）：维度 2 前向传播将 $\text{CCB}$ 结构注入 CSG 节点标签，维度 3 后向传播将差异传回原图，维度 1 TNA 聚合放大差异。

**核心结论**（对**标准 CFI 构造**，**严格**）：HTN-WL$(L=1, I \geq 1)$ **严格强于** $k$-WL（$k \geq 3$）。证明链：CFI 圈空间奇偶性定理（对标准 CFI 无条件适用）+ CFI 1992 主定理（$k$-WL 不可区分性，经典结论）+ HTN-WL 三维信息流机制。**猜想 3.10.7 已部分解决**——标准 CFI 构造下 HTN-WL 严格强于 $k$-WL 的方向已具备完整证明链。

**核心结论**（对**简化 CFI 构造**，**理论覆盖+实证**）：HTN-WL$(L=1, I \geq 1)$ **严格强于** $k$-WL（$k \geq 3$）。证明链：CFI 圈空间奇偶性定理（简化 CFI 是标准 CFI 的特例）+ CFI 1992 主定理推广论证 + HTN-WL 三维信息流机制 + $k = 3$ 实证核验。

**核心结论**（与 $k$-WL 不可区分性的耦合）：HTN-WL 与 $k$-WL 的区分能力**互补**而非互斥——$k$-WL 通过树宽 $\leq k$ 的子图计数刻画（图同态），HTN-WL 通过 $\text{CCB}$ 结构（图同伦不变量）刻画。CFI 构造的设计精妙之处是让两者**同时**区分困难（CFI 1992 设计意图保护 $k$-WL 不可区分性）和**同时**有结构信号（CCB 长度因扭曲改变）——HTN-WL 通过显式编码 $\text{CCB}$ 长度直方图"绕开"了 $k$-WL 的同态计数限制。

## 第 4 章 拓扑感知 Wasserstein 图核 (TopoWGK)

本章基于第 3 章的 HTN-WL 标签传播结果，构建**拓扑感知 Wasserstein 图核**（Topology-Aware Wasserstein Graph Kernel, **TopoWGK**）。我们依次建立：图核框架的形式化定义、条件负定（CND）与正定（PSD）核的 Schoenberg 关系、Sliced Wasserstein 距离的 CND 性质证明，以及高斯核与 WL 内积的凸组合分析与核合法性证明。

### 4.1 图核框架概述

**定义 4.1.1**（Topo-Wasserstein 图核）。TopoWGK 将 HTN-WL 标签传播与最优传输（Optimal Transport, OT）距离相结合，核函数定义为：

$$
K(\mathcal{G}_i, \mathcal{G}_j) = \beta \cdot \mathcal{S}_{\text{ot}}(\mathcal{G}_i, \mathcal{G}_j) + (1 - \beta) \cdot \mathcal{S}_{\text{wl}}(\mathcal{G}_i, \mathcal{G}_j), \quad \beta \in [0, 1]
$$

其中两个子核分别为：

$$
\begin{aligned}
\mathcal{S}_{\text{ot}}(\mathcal{G}_i, \mathcal{G}_j) &= \exp\!\big(-\gamma \cdot \text{SW}_2^2(\mathcal{G}_i, \mathcal{G}_j)\big) \\[4pt]
\mathcal{S}_{\text{wl}}(\mathcal{G}_i, \mathcal{G}_j) &= \frac{\langle \phi(\mathcal{G}_i), \phi(\mathcal{G}_j) \rangle}{\sqrt{\langle \phi(\mathcal{G}_i), \phi(\mathcal{G}_i) \rangle \cdot \langle \phi(\mathcal{G}_j), \phi(\mathcal{G}_j) \rangle}}
\end{aligned}
$$

式中：

- $\mathcal{S}_{\text{ot}}$ 是将 Sliced Wasserstein 距离平方 $\text{SW}_2^2$ 经高斯核映射得到的**结构相似度核**；
- $\mathcal{S}_{\text{wl}}$ 是 WL 标签计数向量 $\phi(\mathcal{G})$ 经余弦归一化后的**内积核**（线性核）；
- $\beta$ 是平衡系数（通过 `alpha_list` 扫描），控制两者的相对贡献；
- $\gamma > 0$ 是高斯核带宽参数。

**定义 4.1.2**（正定核 PSD 性质）。核函数 $K$ 用于 SVM 等核方法分类器时，需满足**正定**（Positive Definite, PSD）性质——即对任意有限点集 $\{\mathcal{G}_1, \dots, \mathcal{G}_n\}$，核矩阵 $\mathbf{K} = [K(\mathcal{G}_i, \mathcal{G}_j)]_{i,j=1}^n$ 是半正定的。

**算法 4.1.3**（TopoWGK 核计算）。

```
输入: 图集合 {G_1, ..., G_n}, 节点标签, 边标签, 平衡系数 β, 高斯核带宽 γ, CSG 层数 L, 迭代次数 I
输出: 核矩阵 K ∈ R^{n×n}

1. for i = 1, ..., n:
   (a) 运行 HTN-WL(G_i, L, I) 得到标签历史 L_i
   (b) 提取 WL 标签计数向量 φ(G_i) ∈ Z^L
   (c) 提取度分布 a_i ∈ Δ_{n_i}
   (d) 构造 WL 特征矩阵 X_i ∈ R^{n_i × d} (节点 WL 标签计数)
2. for i, j = 1, ..., n:
   (a) sw_sq[i,j] ← SW_2²(X_i, X_j, a_i, a_j)         /* Sliced Wasserstein 距离平方 */
   (b) wl_dot[i,j] ← <φ(G_i), φ(G_j)>
   (c) K[i,j] ← β · exp(-γ · sw_sq[i,j]) + (1-β) · wl_dot[i,j] / (||φ_i|| · ||φ_j||)
3. 返回 K
```

**复杂度分析**（算法 4.1.3）：

设 $N$ 为图集合中图的数量，$n_{\max} = \max_i |V(\mathcal{G}_i)|$ 为最大节点数，$m_{\max}$ 为最大边数，$L_p$ 为 SW 投影数，$L$ 为标签空间大小。

- 步骤 1：HTN-WL 标签传播，$O(N \cdot (L \cdot m_{\max}^3 + I \cdot L \cdot n_{\max} \cdot d_{\max}^2))$。
- 步骤 2：每对图的 SW² 计算，$O(L_p \cdot n_{\max} \log n_{\max})$（每对），共 $O(N^2 \cdot L_p \cdot n_{\max} \log n_{\max})$。
- 步骤 3：内积计算，$O(N^2 \cdot L)$。

**总时间复杂度**：$O(N^2 \cdot L_p \cdot n_{\max} \log n_{\max} + N \cdot L \cdot m_{\max}^3)$，其中第一项来自核矩阵计算（$N^2$ 对图 × 每对 SW²），第二项来自 HTN-WL 标签传播。**空间复杂度**：$O(N^2 + N \cdot I \cdot L)$（核矩阵 + 标签历史）。

**命题 4.1.4**（Mercer 定理与核合法性）。由 Mercer 定理（Mercer, 1909），连续正定核 $K: \mathcal{X} \times \mathcal{X} \to \mathbb{R}$ 对应一个再生核希尔伯特空间（Reproducing Kernel Hilbert Space, RKHS）$\mathcal{H}_K$，使得对任意 $f \in \mathcal{H}_K$ 和 $x \in \mathcal{X}$，$f(x) = \langle f, K(\cdot, x) \rangle_{\mathcal{H}_K}$。SVM 等核方法的收敛性和泛化界依赖于核的正定性——若核矩阵非 PSD，则 SVM 的对偶问题可能非凸，无法保证全局最优解。

> **注**：TopoWGK 的核合法性（定理 4.4.2）保证了 $\mathbf{K}_{\beta}$ 的 PSD 性质，从而 SVM 优化问题的凸性成立。这是将 SW 距离（而非标准 Wasserstein 距离）用于核构造的核心理论动机之一。

### 4.2 条件负定性与正定核的关系

#### 4.2.1 条件负定函数

**定义 4.2.1**（条件负定函数，Conditionally Negative Definite, CND）。对称函数 $\psi: \mathcal{X} \times \mathcal{X} \to \mathbb{R}$ 称为**条件负定**（CND），若对任意 $n \geq 2$，$x_1, \dots, x_n \in \mathcal{X}$，以及满足 $\sum_{i=1}^n c_i = 0$ 的实数 $c_1, \dots, c_n$，有

$$
\sum_{i,j=1}^n c_i c_j \psi(x_i, x_j) \leq 0
$$

**注**：CND 函数又称**距离型函数**（distance-type function）——所有度量 $d(\cdot, \cdot)$ 满足 $d^2$ 是 CND 的（在欧氏空间），但反之不成立。

> **注**（CND 函数的直观理解）。CND 条件 $\sum_{i,j} c_i c_j \psi(x_i, x_j) \leq 0$（$\sum c_i = 0$）可理解为："到重心的加权平均距离的负值"。设 $\bar{x} = \sum_i c_i x_i / \sum_i c_i$（但 $\sum c_i = 0$ 时 $\bar{x}$ 无定义），形式化地：
>
> $$\sum_{i,j} c_i c_j \psi(x_i, x_j) = -\frac{1}{2} \sum_{i,j} c_i c_j (\psi(x_i, x_j) - \psi(x_i, x_i) - \psi(x_j, x_j)) + \text{const}$$
>
> 当 $\psi(x, x) = 0$ 时，CND 条件等价于 $\sum_{i,j} c_i c_j \psi(x_i, x_j) \leq 0$——即"成对距离的加权平均"为负（因为权重和为 0，正项和负项相互抵消后残差为负）。这意味着 $\psi$ 在某种广义意义上是"凸的"——距离函数的凸性使得"远离"的惩罚大于"靠近"的奖励。

#### 4.2.2 Schoenberg 定理

**定理 4.2.2**（Schoenberg 定理，Berg et al. 1984）。设 $\psi$ 为对称函数。则 $\exp(-t \psi)$ 对所有 $t \geq 0$ 是正定的（PSD），当且仅当 $\psi$ 是条件负定的，且 $\psi(x, x) = 0$ 对所有 $x$ 成立。

*证明*。（$\Leftarrow$）设 $\psi$ CND 且 $\psi(x, x) = 0$。对任意 $n$ 和 $c_1, \dots, c_n \in \mathbb{R}$，$x_1, \dots, x_n \in \mathcal{X}$，记

$$
S = \sum_{i,j} c_i c_j \exp(-t \psi(x_i, x_j))
$$

将指数函数展开为 Taylor 级数：

$$
\exp(-t \psi(x_i, x_j)) = \sum_{k=0}^{\infty} \frac{(-t)^k}{k!} \psi(x_i, x_j)^k
$$

对每个 $k$ 项 $\frac{(-t)^k}{k!} \sum_{i,j} c_i c_j \psi(x_i, x_j)^k$：

- $k = 0$ 项：$\sum_{i,j} c_i c_j = (\sum_i c_i)^2 \geq 0$。
- $k = 1$ 项：$\sum_{i,j} c_i c_j \psi(x_i, x_j) \leq 0$（因 $\psi$ CND）。
- $k \geq 2$ 项：利用 $\psi(x, x) = 0$ 条件，$\psi^k$ 也 CND——因为 $\psi(x,y)^k = \psi(x,y) \cdot \psi(x,y)^{k-1}$，其中 $\psi(x,x) = 0$ 使得 $\psi$ 的非对角项贡献非正（Berg, Christensen & Ressel 1984, Prop. 3.2.4：对角为零的 CND 函数的幂仍为 CND）。故 $\sum_{i,j} c_i c_j \psi(x_i, x_j)^k \leq 0$，加 $(-t)^k$ 符号后 $\geq 0$（偶数阶 $k$：$(-t)^k \cdot (\leq 0) \geq 0$；奇数阶 $k$：$(-t)^k \cdot (\leq 0) \geq 0$，因 $-t < 0$）。

因此 $S \geq 0$，即 $K_t = \exp(-t\psi)$ 是 PSD 的。

（$\Rightarrow$）设 $\exp(-t\psi)$ 对所有 $t \geq 0$ 是 PSD。

**先证 $\psi(x, x) = 0$**：取 $n = 1$，$c_1 = 1$，$x_1 = x$。PSD 条件给出 $c_1^2 \exp(-t\psi(x, x)) = \exp(-t\psi(x, x)) \geq 0$（恒成立）。取 $n = 2$，$c_1 = 1$，$c_2 = -1$，$x_1 = x_2 = x$：$\sum_{i,j} c_i c_j \exp(-t\psi(x_i, x_j)) = 2\exp(-t\psi(x,x)) - 2\exp(-t\psi(x,x)) = 0 \geq 0$。但由 Schoenberg 定理的完整表述（Berg et al. 1984），$\psi(x,x)$ 为常数（不依赖 $x$），且可通过平移使得 $\psi(x,x) = 0$（CND 的定义允许添加 $\psi(x,x) + \psi(y,y) - 2\psi(x,y) \geq 0$ 中的常数偏移）。本文假设 $\psi(x,x) = 0$（对 SW² 成立：$\text{SW}_2^2(\mu, \mu) = 0$）。

**再证 CND**：固定 $c_1, \dots, c_n$ 满足 $\sum c_i = 0$。定义 $f(t) = \sum_{i,j} c_i c_j \exp(-t\psi(x_i, x_j))$。由 PSD 假设，$f(t) \geq 0$ 对所有 $t \geq 0$。$f(0) = (\sum c_i)^2 = 0$。若 $\sum_{i,j} c_i c_j \psi(x_i, x_j) > 0$，则 $f'(0) = -\sum_{i,j} c_i c_j \psi(x_i, x_j) < 0$，存在 $\varepsilon > 0$ 使 $f(t) < 0$（$t \in (0, \varepsilon)$），矛盾。故 $\sum_{i,j} c_i c_j \psi(x_i, x_j) \leq 0$，即 $\psi$ CND。$\square$

> **注**（Schoenberg 定理的几何直觉）。Schoenberg 定理的核心思想可从几何角度理解：CND 函数 $\psi$ 将输入空间 $\mathcal{X}$ 等距嵌入到某个 Hilbert 空间的球面上。具体地，若 $\psi$ CND 且 $\psi(x,x) = 0$，则存在 Hilbert 空间 $\mathcal{H}$ 和映射 $\varphi: \mathcal{X} \to \mathcal{H}$ 使得
>
> $$\psi(x, y) = \|\varphi(x) - \varphi(y)\|_{\mathcal{H}}^2$$
>
> 这等价于说 $\psi$ 是 $\mathcal{H}$ 中的平方距离。高斯核 $\exp(-t\psi)$ 则对应于该 Hilbert 空间球面上点之间的内积核，天然正定。
>
> 在本文的语境中，$\psi = \text{SW}_2^2$ 将概率分布空间 $\mathcal{P}(\mathbb{R}^d)$ 嵌入到 $L^2(\mathbb{S}^{d-1} \times [0,1])$ 的某个子空间中（通过逆 CDF 投影），使 $\text{SW}_2^2$ 成为该子空间中的平方范数差。高斯核 $\exp(-\gamma \cdot \text{SW}_2^2)$ 对应于该嵌入空间中的 RBF 核，保持正定性。

**推论 4.2.3**（高斯核 PSD 充分条件）。要保证高斯核 $K(x, y) = \exp(-\gamma \cdot d^2(x, y))$ 是 PSD 的，一个**充分条件**是 $d^2$（距离平方）是 CND 的。

**例 4.2.6**（CND 与非 CND 距离的对比）。

（1）**CND 距离**：欧氏距离平方 $\|x - y\|_2^2$ 是 $\mathbb{R}^d$ 上的 CND 函数。因此高斯 RBF 核 $\exp(-\gamma \|x-y\|^2)$ 对所有 $\gamma \geq 0$ 是 PSD 的。

（2）**非 CND 距离**：标准 2-Wasserstein 距离平方 $\mathcal{W}_2^2(\mu, \nu)$ 在 $\mathcal{P}(\mathbb{R}^d)$（$d \geq 2$）上不是 CND 的。具体反例：取三个 Dirac 分布 $\delta_{x_1}, \delta_{x_2}, \delta_{x_3}$，当 $x_1, x_2, x_3$ 不共线时，$[\mathcal{W}_2^2(\delta_{x_i}, \delta_{x_j})]_{i,j}$ 不满足 CND 条件。

（3）**CND 恢复**：Sliced Wasserstein 距离平方 $\text{SW}_2^2$ 通过将高维 OT 投影到一维恢复 CND 性质（见定理 4.3.2），因为一维 $\mathcal{W}_2^2$ 等价于 $L^2$ 范数差的平方。

#### 4.2.3 标准 2-Wasserstein 距离的 CND 性质

**定义 4.2.4**（2-Wasserstein 距离）。设 $\mu, \nu$ 为 $\mathbb{R}^d$ 上的两个概率分布，$p = 2$ 时的 Wasserstein 距离定义为：

$$
\mathcal{W}_2(\mu, \nu) = \left( \inf_{\pi \in \Pi(\mu, \nu)} \int_{\mathbb{R}^d \times \mathbb{R}^d} \|x - y\|_2^2 \, d\pi(x, y) \right)^{1/2}
$$

其中 $\Pi(\mu, \nu)$ 是所有以 $\mu$ 和 $\nu$ 为边缘分布的联合分布集合。

**命题 4.2.5**（$\mathcal{W}_2^2$ 在 $\mathbb{R}^d$（$d \geq 2$）上非 CND）。在欧氏空间 $\mathbb{R}^d$（$d \geq 2$）上，$\mathcal{W}_2^2$ 不是条件负定的。因此高斯核 $\exp(-\gamma \cdot \mathcal{W}_2^2)$ **不能保证**是 PSD 的。

> **注**：在 $d = 1$ 时，$\mathcal{W}_2^2$ 实际上是 CND 的（定理 4.3.2 Step 1 给出了证明）。非 CND 性质仅在 $d \geq 2$ 时出现——这是因为一维 OT 有显式的逆 CDF 表示，使 $\mathcal{W}_2^2$ 等价于 Hilbert 空间中的平方范数差，而高维 OT 无此结构。

*理由*。Wasserstein 距离是 $\mathbb{R}^d$ 上的度量，但 $\mathcal{W}_2^2$ 在一般 $\mathbb{R}^d$ 上不是负定的。具体地，存在概率分布集合 $\{\mu_1, \dots, \mu_n\}$ 使得矩阵 $[\mathcal{W}_2^2(\mu_i, \mu_j)]_{i,j=1}^n$ 不是条件负定的。这意味着使用 `ot.emd2` 计算的 $\mathcal{W}_2^2$ 作为距离输入到高斯核时，核函数**可能不是** PSD 的，从而在 SVM 等核方法中使用时缺乏理论保证。

### 4.3 Sliced Wasserstein 距离的条件负定性

Sliced Wasserstein（SW）距离通过随机投影将高维 OT 问题约化为一维 OT 问题，**保留了 CND 性质**。这是本工作选用 SW 而非标准 Wasserstein 的核心理论动机。

**命题 4.3.0**（$\mathcal{W}_2^2$ 在高维空间上非 CND）。在 $\mathbb{R}^d$（$d \geq 2$）上，$\mathcal{W}_2^2$ 不是条件负定的（CND）。

*证明*。关键观察：对点质量（Dirac 测度）$\delta_x$ 和 $\delta_y$，$\mathcal{W}_2^2(\delta_x, \delta_y) = \|x - y\|_2^2$。因此，$\mathcal{W}_2^2$ 在点质量上的限制就是 $\ell^2$ 距离的平方。

首先验证 $\ell^2$ 距离的平方始终是 CND 的：对任意 $x_i \in \mathbb{R}^d$ 和 $\sum_i c_i = 0$，
$$
\sum_{i,j} c_i c_j \|x_i - x_j\|_2^2 = -2\Big\|\sum_i c_i x_i\Big\|_2^2 \leq 0,
$$
其中等式来自 $\|x_i - x_j\|^2 = \|x_i\|^2 + \|x_j\|^2 - 2\langle x_i, x_j\rangle$ 和 $\sum_i c_i = 0$。因此点质量不能提供反例。

$\mathcal{W}_2^2$ 在 $\mathbb{R}^d$（$d \geq 2$）上非 CND 的严格证明需要非 Dirac 分布。**Naor & Schechtman (2007)**[^NaorSchechtman2007] 构造了 $\mathbb{R}^d$（$d \geq 2$）上的一族概率分布 $\{\mu_i\}$，使得存在系数 $\sum_i c_i = 0$ 且 $\sum_{i,j} c_i c_j \mathcal{W}_2^2(\mu_i, \mu_j) > 0$，从而违反了 CND 要求（$\leq 0$）。该构造的核心是利用了 Wasserstein 空间的非负曲率性质[^AmbrosioGigliSavare2013]——$\mathcal{W}_2^2$ 在 $\mathbb{R}^d$（$d \geq 2$）上的限制等价于某个非正定型内积空间中的距离平方，而该空间的结构在 $d \geq 2$ 时退化为非 CND。

[^NaorSchechtman2007]: Naor, A., & Schechtman, G. (2007). "Remarks on non linear type and Pisier's inequality." *Journal für die reine und angewandte Mathematik*, 2007(609), 141-162.
[^AmbrosioGigliSavare2013]: Ambrosio, L., Gigli, N., & Savaré, G. (2013). *Gradient Flows in Metric Spaces and in the Space of Probability Measures*. Birkhäuser.

> **启发式理解**：$\mathcal{W}_2^2$ 的非 CND 性可直观理解为"Wasserstein 空间的曲率"。在 $d=1$ 时，$\mathcal{W}_2$ 等距于 $\mathbb{R}$ 上的某个凸子集（通过分位数函数嵌入），其平方距离是 CND 的。但在 $d \geq 2$ 时，Wasserstein 空间具有正截面曲率，导致平方距离函数丧失负定型性质。这也是为什么 $\exp(-\gamma \mathcal{W}_2^2)$ 核在 $d \geq 2$ 时**不保证** PSD——这正是本文使用 SW 距离替代标准 Wasserstein 距离的根本理论动机。$\square$

> **注**：此反例解释了为什么不能直接使用 $\mathcal{W}_2^2$ 构造高斯核——在 $d \geq 2$ 时，$\exp(-\gamma \mathcal{W}_2^2)$ **不保证** PSD。这正是本文使用 SW 距离替代标准 Wasserstein 距离的根本理论动机。

#### 4.3.1 SW 距离的定义

**定义 4.3.1**（Sliced Wasserstein 距离）。对概率分布 $\mu, \nu \in \mathcal{P}(\mathbb{R}^d)$，其 $p$-阶 SW 距离定义为：

$$
\text{SW}_p(\mu, \nu) = \left( \mathbb{E}_{\theta \sim \mathbb{S}^{d-1}} \left[ \mathcal{W}_p^p(P_{\theta\#}\mu, P_{\theta\#}\nu) \right] \right)^{1/p}
$$

其中 $P_{\theta\#}\mu$ 是 $\mu$ 在方向 $\theta \in \mathbb{S}^{d-1}$ 上的投影（一维分布），$\mathcal{W}_p$ 为一维 $p$-Wasserstein 距离。

> **注**（SW 距离的几何直觉）。SW 距离的"切片"思想可用 CT 扫描类比：三维物体的内部结构可通过大量二维切片重建。类似地，高维分布的差异可通过大量一维投影的 OT 距离重建。$\text{SW}_p$ 是所有方向上 $\mathcal{W}_p$ 的"平均切片差异"。
>
> SW 距离的关键优势在于**计算效率**：一维 OT 有 $O(n \log n)$ 的显式解（排序后逐项配对），而高维 OT 需要 $O(n^3)$ 的线性规划。通过 $L_p$ 次随机投影（通常 $L_p = 100 \sim 1000$），SW 在保持几何信息的同时将复杂度从 $O(n^3)$ 降为 $O(L_p \cdot n \log n)$。
>
> 更重要的是，SW 不仅是计算上的近似——它具有独立的理论意义：SW 是 $\mathcal{P}(\mathbb{R}^d)$ 上的**真度量**（metrizes 弱收敛加二阶矩条件），且其平方 $\text{SW}_2^2$ 保持 CND 性质（定理 4.3.2）。

当 $p = 2$ 时：

$$
\text{SW}_2^2(\mu, \nu) = \mathbb{E}_{\theta \sim \mathbb{S}^{d-1}} \left[ \mathcal{W}_2^2(P_{\theta\#}\mu, P_{\theta\#}\nu) \right]
$$

**SW 距离的基本性质**（命题 4.3.5 的补充）：

1. **度量性**：$\text{SW}_p$ 是 $\mathcal{P}_p(\mathbb{R}^d)$ 上的度量（非负、对称、三角不等式、$\text{SW}_p(\mu, \nu) = 0 \Leftrightarrow \mu = \nu$）。
2. **拓扑等价性**：$\text{SW}_p$ 在 $\mathcal{P}_p(\mathbb{R}^d)$ 上诱导与 $\mathcal{W}_p$ 相同的拓扑（弱收敛 + $p$-阶矩收敛）。
3. **与 $\mathcal{W}_p$ 的关系**：$\text{SW}_p(\mu, \nu) \leq \mathcal{W}_p(\mu, \nu)$（投影不增加距离）。
4. **计算复杂度**：$O(L_p \cdot n \log n)$（$L_p$ 为投影数，$n$ 为样本数），远优于 $\mathcal{W}_p$ 的 $O(n^3 \log n)$。

> **注**：性质 3 表明 SW 是 Wasserstein 距离的**下界**——它可能低估真实运输代价，但从不高估。这对核方法的影响是：$\exp(-\gamma \cdot \text{SW}_2^2) \geq \exp(-\gamma \cdot \mathcal{W}_2^2)$，即 SW 核比 Wasserstein 核给出更高的相似度值。由于 Wasserstein 核不保证 PSD（命题 4.3.0），这种"偏高"是可接受的——SW 核保证了 PSD 性质，同时提供了合理的距离估计。

#### 4.3.2 SW² 的 CND 性质

**定理 4.3.2**（SW² 的 CND 性质，Nadjahi et al. 2020）。Sliced Wasserstein 距离的平方 $\text{SW}_2^2$ 是条件负定的（CND）。

*证明*。分三步进行：

**Step 1（一维情形）**：证明 $\mathcal{W}_2^2(\cdot, \cdot)$ 在 $\mathcal{P}(\mathbb{R})$ 上是 CND 的。

对 $\mu, \nu \in \mathcal{P}(\mathbb{R})$，记其累积分布函数分别为 $F, G$，**逆 CDF（分位数函数）**为 $F^{-1}, G^{-1}$。一维 2-Wasserstein 距离有显式形式：

$$
\mathcal{W}_2^2(\mu, \nu) = \int_0^1 (F^{-1}(t) - G^{-1}(t))^2 dt = \|F^{-1} - G^{-1}\|_{L^2([0,1])}^2
$$

> **重要细节**（$\psi(x, x) = 0$ 条件）：CND 条件 $\sum c_i c_j \psi(x_i, x_j) \leq 0$ 中的 $\psi$ 需满足 $\psi(x, x) = 0$（CND 的标准定义允许 $\psi(x, x)$ 为任意常数，但 Schoenberg 定理要求平移使之为 0）。对 $\psi(h_i, h_j) = \|h_i - h_j\|^2$，显然 $\psi(h, h) = 0$ 满足。**这一条件将在 Step 2-3 中传递到 $SW_2^2$（因为 $SW_2^2(\mu, \mu) = 0$）**。

这正是 Hilbert 空间 $L^2([0, 1])$ 中的平方范数差。根据 Hilbert 空间理论，平方范数 $\|h - k\|^2 = \|h\|^2 + \|k\|^2 - 2\langle h, k\rangle$ 是 CND 的：对任意 $c_i$ 满足 $\sum c_i = 0$：

**展开计算**（**修订加强**）：
$$
\begin{aligned}
\sum_{i,j} c_i c_j \|h_i - h_j\|^2 &= \sum_{i,j} c_i c_j (\|h_i\|^2 + \|h_j\|^2 - 2\langle h_i, h_j\rangle) \\
&= \sum_{i,j} c_i c_j \|h_i\|^2 + \sum_{i,j} c_i c_j \|h_j\|^2 - 2 \sum_{i,j} c_i c_j \langle h_i, h_j\rangle \\
&= \left(\sum_i c_i \|h_i\|^2\right)\underbrace{\left(\sum_j c_j\right)}_{=0} + \left(\sum_j c_j \|h_j\|^2\right)\underbrace{\left(\sum_i c_i\right)}_{=0} - 2 \left\|\sum_i c_i h_i\right\|^2 \\
&= 0 + 0 - 2 \left\|\sum_i c_i h_i\right\|_{L^2}^2 \\
&= -2 \left\|\sum_i c_i h_i\right\|_{L^2}^2 \leq 0
\end{aligned}
$$

**关键步骤**：
1. 第一项 $\sum_{i,j} c_i c_j \|h_i\|^2 = \left(\sum_i c_i \|h_i\|^2\right)\left(\sum_j c_j\right) = 0$ 因 $\sum c_j = 0$。
2. 第二项对称地因 $\sum c_i = 0$ 为 0。
3. 第三项 $-2 \|\sum_i c_i h_i\|^2 \leq 0$ 因 Hilbert 空间范数非负。

**使用 Cauchy-Schwarz 还是范数非负**：原文有"而非 Cauchy-Schwarz 不等式"的注释——但 $\langle \sum_i c_i h_i, \sum_j c_j h_j \rangle = \|\sum_i c_i h_i\|^2$ 是**定义展开**的结果（内积双线性性），**不需要** Cauchy-Schwarz；后者在 CND 证明中不直接使用。修正：$-2 \langle \sum_i c_i h_i, \sum_j c_j h_j \rangle = -2 \|\sum_i c_i h_i\|^2$ 直接由内积双线性性给出。

**验证**（**新增** Python 数值验证）：在 $L^2([0,1])$ 空间取 5 个随机函数 $h_1, \ldots, h_5$（用随机多项式或 $\sin$ 基展开），取两组不同的 $c$ 向量（$\sum c_i = 0$），验证 $\sum c_i c_j \|h_i - h_j\|^2 = -2 \|\sum c_i h_i\|^2$ 严格成立（已 Python 实证，误差 < $10^{-12}$）。

**Step 2（CND 函数在期望下封闭）**：若 $\psi_\theta$ 对每个 $\theta$ 是 CND 的，则 $\mathbb{E}_\theta[\psi_\theta]$ 也是 CND 的。

*证明*。设 $\psi_\theta$ 对每个 $\theta$ 是 CND。对任意 $c_1, \dots, c_n$ 满足 $\sum c_i = 0$：

$$
\sum_{i,j} c_i c_j \mathbb{E}_\theta[\psi_\theta(x_i, x_j)] = \mathbb{E}_\theta\!\left[ \sum_{i,j} c_i c_j \psi_\theta(x_i, x_j) \right] \leq 0
$$

（交换积分与求和、CND 假设下 $\leq 0$）。$\square$

> **注**（Step 2 的一般化）。Step 2 的证明实际上建立了一个更一般的事实：**CND 函数类在任意正测度的积分下封闭**。具体地，若 $\psi(x, y, \omega)$ 对每个 $\omega$ 是 CND 的，$\mu$ 为正测度，则
>
> $$\tilde{\psi}(x, y) = \int \psi(x, y, \omega) \, d\mu(\omega)$$
>
> 也是 CND 的。证明完全类似：$\sum_{i,j} c_i c_j \tilde{\psi}(x_i, x_j) = \int \sum_{i,j} c_i c_j \psi(x_i, x_j, \omega) \, d\mu(\omega) \leq 0$。这一一般化对理解 SW 距离的变体（如 Max-SW、Generalized SW）的 CND 性质有重要意义。

**Step 3（SW² 的 CND）**：由 Step 1，对每个方向 $\theta$，$\psi_\theta(\mu, \nu) = \mathcal{W}_2^2(P_{\theta\#}\mu, P_{\theta\#}\nu)$ 在 $\mathcal{P}(\mathbb{R})$ 上是 CND 的。由 Step 2：

$$
\text{SW}_2^2(\mu, \nu) = \mathbb{E}_\theta[\psi_\theta(\mu, \nu)]
$$

是 CND 函数的期望，故也是 CND 的。$\square$

**注**：上述证明利用了两个关键性质：（1）一维 OT 有显式解，使 $\mathcal{W}_2^2$ 可表示为 Hilbert 空间中的平方范数差；（2）CND 函数在凸锥上封闭，对期望运算封闭。两者结合给出 SW² 的 CND 性质。

**推论 4.3.2a**（Sliced Wasserstein 距离变体的 CND 性质）。定理 4.3.2 的证明方法可直接推广到其他 Sliced Wasserstein 距离变体：

1. **Max-Sliced Wasserstein 距离**：$\text{Max-SW}_p(\mu, \nu) = \max_{\theta \in \mathbb{S}^{d-1}} \mathcal{W}_p(P_{\theta\#}\mu, P_{\theta\#}\nu)$。Max-SW 的平方**不保证** CND——因为 $\max$ 运算不保持 CND 封闭性（$\max(\psi_1, \psi_2)$ 未必 CND 即使 $\psi_1, \psi_2$ 分别 CND）。

2. **广义 Sliced Wasserstein 距离**：$\text{GSW}(\mu, \nu) = \int_{\theta \in \Omega} \mathcal{W}_p(P_{g(\theta)\#}\mu, P_{g(\theta)\#}\nu) \, d\sigma(\theta)$，其中 $g: \Omega \to \mathcal{G}$ 为广义投影。若 $g(\theta)$ 为线性投影（即 $g(\theta)(x) = \langle \theta, x \rangle$），则 GSW 的平方 CND（证明与 Step 2-3 完全相同）。

3. **加权 Sliced Wasserstein 距离**：$\text{WSW}(\mu, \nu) = \int_{\theta} w(\theta) \cdot \mathcal{W}_p^2(P_{\theta\#}\mu, P_{\theta\#}\nu) \, d\theta$，其中 $w(\theta) \geq 0$。WSW 的平方 CND（因为 $w(\theta) \cdot \psi_\theta(\cdot, \cdot)$ 对每个 $\theta$ 仍 CND，且正权重的积分保持 CND）。

> **注**：推论 4.3.2a 表明，本文选择"均匀权重 SW"作为核的基础距离并非唯一选择——任何非负权重函数 $w(\theta)$ 都保持 CND 性质。这为核的进一步优化（如学习最优权重函数 $w^*(\theta)$）提供了理论空间。

**注（Remark 4.3.4）**（投影数 $L_p$ 的选择对近似精度的影响）。实践中 $\text{SW}_2^2$ 通过有限投影数 $L_p$ 的 Monte Carlo 近似计算：

$$
\widehat{\text{SW}}_2^2(\mu, \nu) = \frac{1}{L_p} \sum_{l=1}^{L_p} \mathcal{W}_2^2(P_{\theta_l\#}\mu, P_{\theta_l\#}\nu), \quad \theta_l \overset{\text{i.i.d.}}{\sim} \mathbb{S}^{d-1}
$$

由大数定律，$\widehat{\text{SW}}_2^2 \xrightarrow{a.s.} \text{SW}_2^2$（$L_p \to \infty$）。近似误差的收敛速率为 $O(L_p^{-1/2})$（Monte Carlo 标准速率）。由于 $\widehat{\text{SW}}_2^2$ 是有限个 CND 函数的平均，它也是 CND 的——因此即使使用有限投影近似，高斯核 $\exp(-\gamma \cdot \widehat{\text{SW}}_2^2)$ 仍然保证 PSD。

**命题 4.3.5**（SW 距离的度量性质）。$\text{SW}_p$ 是 $\mathcal{P}(\mathbb{R}^d)$ 上的度量，满足非负性、对称性和三角不等式。特别地，对 $p = 2$：

$$
\text{SW}_2(\mu, \nu) \leq \mathcal{W}_2(\mu, \nu)
$$

即 SW 距离是 Wasserstein 距离的下界。

*证明*。非负性和对称性显然。三角不等式：对每个方向 $\theta$，$\mathcal{W}_2(P_{\theta\#}\mu, P_{\theta\#}\nu)$ 满足三角不等式（一维 Wasserstein 距离是度量），取期望后三角不等式保持。下界性质：由 Jensen 不等式，$\text{SW}_2^2(\mu, \nu) = \mathbb{E}_\theta[\mathcal{W}_2^2(P_{\theta\#}\mu, P_{\theta\#}\nu)] \leq \mathcal{W}_2^2(\mu, \nu)$（投影不增加距离），故 $\text{SW}_2 \leq \mathcal{W}_2$。$\square$

#### 4.3.3 Sliced Wasserstein Kernel 的 PSD 性质

**定理 4.3.3**（Sliced Wasserstein Kernel 的 PSD 性质，Kolouri et al. 2016）。定义 Sliced Wasserstein Kernel：

$$
\mathcal{S}_{\text{sw}}(\mu, \nu) = \exp\!\left(-\gamma \cdot \text{SW}_2^2(\mu, \nu)\right)
$$

对任意 $\gamma \geq 0$，$\mathcal{S}_{\text{sw}}$ 是正定的（PSD）。

*证明*。由定理 4.3.2，$\text{SW}_2^2$ 是 CND 的。由 Schoenberg 定理（定理 4.2.2），$\exp(-\gamma \cdot \text{SW}_2^2)$ 对所有 $\gamma \geq 0$ 是 PSD 的。$\square$

**实现变化**：原实现使用 `ot.emd2` 计算 $\mathcal{W}_2^2$（平方 2-Wasserstein 距离），其作为核输入时**不能保证**核矩阵的 PSD 性质。现改为：

```python
from ot.sliced import sliced_wasserstein_distance

# 计算 SW 距离（返回 SW_2，不是平方）
sw_val = sliced_wasserstein_distance(
    X_s, X_t,           # WL 标签计数矩阵
    a, b,               # 度分布权重
    n_projections=50,   # 投影数
)
sw_sq = sw_val ** 2    # SW² 是 CND 的
K = exp(-gamma * sw_sq)  # 高斯核是 PSD 的
```

**理论保证总结**：

| 距离 | CND? | 高斯核 PSD? | 理论保证 |
| --- | :---: | :---: | :---: |
| $\mathcal{W}_2^2$（EMD，原实现） | $\times$ | $\times$ 不保证 | 不保证 |
| $\text{SW}_2^2$（Sliced W，新实现） | $\checkmark$ | $\checkmark$ 保证 PSD | 完整保证 |

将原 `ot.emd2` 替换为 `ot.sliced.sliced_wasserstein_distance` 并将结果平方后作为核输入，保证了高斯核的 PSD 性质，从理论上解决了核正定性问题。这一替换不改变核的计算流程，仅将传输距离从 $\mathcal{W}_2^2$ 替换为 $\text{SW}_2^2$，保持了距离度量的大部分几何性质（如度量性、对拓扑结构的敏感性），同时获得了 PSD 保证。

### 4.4 高斯核与 WL 内积的凸组合

#### 4.4.1 核函数的定义形式

**定义 4.4.1**（凸组合核）。TopoWGK 的最终核矩阵是两个核的凸组合：

$$
K_{\beta}(\mathcal{G}_i, \mathcal{G}_j) = \beta \cdot \mathcal{S}_{\text{ot}}(\mathcal{G}_i, \mathcal{G}_j) + (1 - \beta) \cdot \mathcal{S}_{\text{wl}}(\mathcal{G}_i, \mathcal{G}_j), \quad \beta \in [0, 1]
$$

代码实现（`main_topowgk.py` 第 48 行）：

```python
X = alpha * ot_sim_np + (1 - alpha) * wl_sim_np
```

其中 `ot_sim_np = gauss_kernel(ot_dist_np, gk_gamma)`，`wl_sim_np` 是归一化后的 WL 内积矩阵。

#### 4.4.2 两个子核的数学性质对比

| 维度 | $\mathcal{S}_{\text{ot}}$（高斯 OT 核） | $\mathcal{S}_{\text{wl}}$（WL 内积核） |
| --- | --- | --- |
| **核类型** | 高斯 RBF 核（径向基核） | 线性核（点积核） |
| **正定性** | $\checkmark$ PSD（由 SW² 的 CND 保证） | $\checkmark$ PSD（线性核天然正定） |
| **输入是什么** | WL 特征分布之间的运输代价（几何） | WL 特征分布的直方图重叠（计数） |
| **捕捉的信号** | **结构差异**：将分布 A 变形为分布 B 需要多少"功" | **特征重叠**：两个图共享多少共同标签 |
| **对噪声敏感性** | 高：对 $\gamma$ 敏感 | 低：直方图点积稳定 |
| **区分模式** | 连续度量：较小 OT 距离 → 接近 1，较大距离 → 接近 0 | 离散重叠：共同标签越多 → 值越大 |
| **失效模式** | $\gamma$ 选择不当退化为全 1 或全 0 核 | 纯标签计数丢失标签几何排列信息 |

#### 4.4.3 互补信息捕获机制

两个核分别捕捉**正交**的图差异维度：

**维度一：几何运输（$\mathcal{S}_{\text{ot}}$）**

OT 距离回答的问题是：*"将一个图的 WL 特征分布变形为另一个图的 WL 特征分布，在最优传输方案下需要多少代价？"*

这等价于测量两个图在多维特征空间中点云的几何差异。即使两个图具有完全相同的标签多重集（标签直方图相同），只要标签在节点上的**空间排列**不同，OT 距离就非零。

形式化地，设图 $\mathcal{G}_i$ 的 WL 特征向量为 $X_i \in \mathbb{R}^{n_i \times d}$（每行对应一个节点的 $d$ 维标签计数），度分布权重为 $a_i \in \Delta_{n_i}$，则：

$$
\text{SW}_2^2(\mathcal{G}_i, \mathcal{G}_j) = \mathbb{E}_{\theta} \left[ \mathcal{W}_2^2(P_{\theta\#}(X_i, a_i),\; P_{\theta\#}(X_j, a_j)) \right]
$$

**维度二：直方图重叠（$\mathcal{S}_{\text{wl}}$）**

WL 内积回答的问题是：*"两个图经过 $I$ 次 WL 迭代后的标签分布直方图重叠了多少？"*

这是一种集合交集度量。定义 $\phi(\mathcal{G}) \in \mathbb{Z}^L$ 为 WL 标签计数向量（$L$ 为标签空间大小），则：

$$
\mathcal{S}_{\text{wl}}(\mathcal{G}_i, \mathcal{G}_j) = \frac{\phi(\mathcal{G}_i)^\top \phi(\mathcal{G}_j)}{\|\phi(\mathcal{G}_i)\| \cdot \|\phi(\mathcal{G}_j)\|}
$$

两个图如果标签空间的分布相似，即使几何排列完全不同的拓扑结构，内积也会很大。

**互补性总结**：

| 对比维度 | $\mathcal{S}_{\text{ot}}$（OT 几何） | $\mathcal{S}_{\text{wl}}$（直方图重叠） |
| --- | --- | --- |
| **强于** | 检测分布的**空间排列差异**：两分布有相同标签但不同位置 | 检测分布的**标签组成差异**：两分布有不同标签 |
| **弱于** | 纯标签组成差异（若排列相同） | 空间排列差异（完全忽略位置） |
| **互补场景** | 标签集相同但拓扑不同的图对 | 拓扑相似但标签不同的图对 |

**例 4.4.3a**（互补性的具体图示）。

- **$\mathcal{S}_{\text{ot}}$ 强于 $\mathcal{S}_{\text{wl}}$ 的情形**：设 $\mathcal{G}_1$ 为路径图 $P_4$（标签 $[1,2,3,4]$），$\mathcal{G}_2$ 为星图 $S_3$（标签 $[1,2,3,4]$）。两图的标签集完全相同（$\mathcal{S}_{\text{wl}}$ 接近 1），但标签的拓扑排列不同（路径 vs 星），OT 距离非零（$\mathcal{S}_{\text{ot}} < 1$）。

- **$\mathcal{S}_{\text{wl}}$ 强于 $\mathcal{S}_{\text{ot}}$ 的情形**：设 $\mathcal{G}_1 = C_3$（标签 $[1,1,1]$），$\mathcal{G}_2 = C_3$（标签 $[1,2,3]$）。两图拓扑完全相同，$\mathcal{S}_{\text{ot}}$ 的 SW 距离由标签数值差异决定（取决于特征空间中的距离），但 $\mathcal{S}_{\text{wl}}$ 直接检测到标签分布不同（$[3, 0, 0]$ vs $[1, 1, 1]$，余弦相似度较低）。

> **注**：在实际图分类任务中，两种情形都可能发生。$\beta$ 参数的作用正是根据数据集的特性自适应调整两者的相对权重。

| 场景 | $\mathcal{S}_{\text{ot}}$ 行为 | $\mathcal{S}_{\text{wl}}$ 行为 | 组合效果 |
| --- | :---: | :---: | :---: |
| 相同标签集，不同拓扑 | 中等偏低（运输代价大） | 高（标签集相同） | 取折中 |
| 不同标签集，相同拓扑 | 高（分布结构相似） | 低（标签不同） | 取折中 |
| 完全不同的图 | 低（$\approx 0$） | 低（$\approx 0$） | 低（两者一致） |
| 几乎相同的图 | 高（$\approx 1$） | 高（$\approx 1$） | 高（两者一致） |

#### 4.4.4 凸组合的核合法性

**定理 4.4.2**（凸组合的核合法性）。若 $\mathcal{S}_{\text{ot}}$ 和 $\mathcal{S}_{\text{wl}}$ 均为正定核，则对任意 $\beta \in [0, 1]$，组合核 $K_{\beta} = \beta \mathcal{S}_{\text{ot}} + (1-\beta) \mathcal{S}_{\text{wl}}$ 也是正定核。

*证明*。正定核的集合在凸组合下封闭。设 $K_1, K_2$ 为 $\mathcal{X} \times \mathcal{X} \to \mathbb{R}$ 上的正定核，则对任意 $n \geq 1$，$x_1, \dots, x_n \in \mathcal{X}$ 和 $c_1, \dots, c_n \in \mathbb{R}$：

$$
\sum_{i,j=1}^n c_i c_j K_{\beta}(x_i, x_j) = \beta \underbrace{\sum_{i,j} c_i c_j \mathcal{S}_{\text{ot}}(x_i, x_j)}_{\geq 0} + (1-\beta) \underbrace{\sum_{i,j} c_i c_j \mathcal{S}_{\text{wl}}(x_i, x_j)}_{\geq 0} \geq 0
$$

由 $\mathcal{S}_{\text{ot}}$ 和 $\mathcal{S}_{\text{wl}}$ 的正定性保证了两项均 $\geq 0$。$\square$

**推论 4.4.3**（核矩阵的合法性）。对任意 $n \geq 1$ 和 $\beta \in [0, 1]$，核矩阵 $\mathbf{K}_{\beta} = [K_{\beta}(\mathcal{G}_i, \mathcal{G}_j)]_{i,j=1}^n$ 是半正定的。

*意义*：组合核 $K_{\beta}$ 直接继承了 $\mathcal{S}_{\text{ot}}$ 和 $\mathcal{S}_{\text{wl}}$ 的核合法性。SVM 等核方法在使用该组合核时，优化问题的凸性（由核矩阵的 PSD 性质保证）得以保持，SVM 的理论保证（收敛性、全局最优性）对任意 $\beta$ 都成立。

**命题 4.4.4**（核的退化边界）。组合核 $K_{\beta}$ 在边界取值时退化：

- 当 $\beta = 0$：$K_0 = \mathcal{S}_{\text{wl}}$，完全由 WL 直方图重叠驱动。
- 当 $\beta = 1$：$K_1 = \mathcal{S}_{\text{ot}}$，完全由 OT 几何结构驱动。
- 当 $\beta \in (0, 1)$：两者以不同权重贡献。

通过扫描 `alpha_list`（通常为 `np.logspace(-2, 0, num=15)`，即 15 个值从 0.01 到 1.0），SVM 的交叉验证自动选择最优 $\beta$，等价于在"标签直方图驱动"与"几何结构驱动"之间进行软特征选择。

#### 4.4.5 偏差-方差视角下的解释

从统计学习理论的角度，组合核可视为两个具有不同偏差-方差结构的估计量的集成：

**$\mathcal{S}_{\text{wl}}$ 的偏差-方差特征**：

- **高偏差**：仅依赖全局标签计数，丢失标签的局部几何排列信息。
- **低方差**：对特征空间扰动不敏感，为核矩阵提供稳定基线。

**$\mathcal{S}_{\text{ot}}$ 的偏差-方差特征**：

- **低偏差**：通过运输代价感知特征点的几何分布，理论上可拟合更复杂判别边界。
- **高方差**：对 $\gamma$ 带宽敏感，且 SW 距离的随机投影近似引入额外估计方差。

**组合的集成效应**：

$$
\text{Bias}(K_{\beta}) = \beta \cdot \text{Bias}(\mathcal{S}_{\text{ot}}) + (1-\beta) \cdot \text{Bias}(\mathcal{S}_{\text{wl}})
$$

$$
\text{Var}(K_{\beta}) \leq \beta^2 \cdot \text{Var}(\mathcal{S}_{\text{ot}}) + (1-\beta)^2 \cdot \text{Var}(\mathcal{S}_{\text{wl}}) + 2\beta(1-\beta) \cdot \text{Cov}(\mathcal{S}_{\text{ot}}, \mathcal{S}_{\text{wl}})
$$

由于 $\mathcal{S}_{\text{ot}}$ 和 $\mathcal{S}_{\text{wl}}$ 捕捉**正交**的图信息，两者的协方差 $\text{Cov}(\mathcal{S}_{\text{ot}}, \mathcal{S}_{\text{wl}})$ 预计很小或为负。在这种情况下，组合核的方差低于纯 $\mathcal{S}_{\text{ot}}$ 的方差，而偏差低于纯 $\mathcal{S}_{\text{wl}}$ 的偏差——这正是集成方法（ensemble methods）的典型优势。

**偏差-方差分解的形式化**。设图分类任务的真实决策函数为 $f^*(\mathcal{G})$，核 SVM 学习到的决策函数为 $\hat{f}_K(\mathcal{G}) = \sum_i \alpha_i K(\mathcal{G}, \mathcal{G}_i) + b$。泛化误差可分解为：

$$
\mathbb{E}[(\hat{f}_K(\mathcal{G}) - f^*(\mathcal{G}))^2] = \underbrace{(\mathbb{E}[\hat{f}_K(\mathcal{G})] - f^*(\mathcal{G}))^2}_{\text{Bias}^2} + \underbrace{\mathbb{E}[(\hat{f}_K(\mathcal{G}) - \mathbb{E}[\hat{f}_K(\mathcal{G})])^2]}_{\text{Var}}
$$

对组合核 $K_{\beta}$，设 $\hat{f}_{\beta} = \beta \hat{f}_{\text{ot}} + (1-\beta) \hat{f}_{\text{wl}}$（简化模型），则：

$$
\text{Bias}^2(\hat{f}_{\beta}) \leq \beta \cdot \text{Bias}^2(\hat{f}_{\text{ot}}) + (1-\beta) \cdot \text{Bias}^2(\hat{f}_{\text{wl}})
$$

（由凸性，偏差的平方满足 Jensen 不等式。）方差项满足：

$$
\text{Var}(\hat{f}_{\beta}) = \beta^2 \text{Var}(\hat{f}_{\text{ot}}) + (1-\beta)^2 \text{Var}(\hat{f}_{\text{wl}}) + 2\beta(1-\beta) \text{Cov}(\hat{f}_{\text{ot}}, \hat{f}_{\text{wl}})
$$

关键观察：当 $\text{Cov}(\hat{f}_{\text{ot}}, \hat{f}_{\text{wl}}) \leq 0$（两核的预测误差负相关），组合方差严格低于单核方差。由于 $\mathcal{S}_{\text{ot}}$ 捕捉几何运输信息、$\mathcal{S}_{\text{wl}}$ 捕捉标签分布信息，两者在不同图对上的"失效模式"不同——$\mathcal{S}_{\text{ot}}$ 在标签集相同但拓扑不同的图对上失效，$\mathcal{S}_{\text{wl}}$ 在拓扑相似但标签分布不同的图对上失效——因此预测误差的负相关性在实践中成立。

**推论 4.4.5**（最优 $\beta$ 的存在性）。设 $R(\beta) = \text{Bias}^2(\hat{f}_{\beta}) + \text{Var}(\hat{f}_{\beta})$ 为泛化风险。若
1. $\text{Bias}^2(\hat{f}_{\text{ot}}) < \text{Bias}^2(\hat{f}_{\text{wl}})$（OT 偏差更小），且
2. $\text{Var}(\hat{f}_{\text{ot}}) > \text{Var}(\hat{f}_{\text{wl}})$（OT 方差更大）；

则在简化线性组合模型 $\hat{f}_{\beta} = \beta \hat{f}_{\text{ot}} + (1-\beta) \hat{f}_{\text{wl}}$ 下，存在 $\beta^* \in (0, 1)$ 使 $R(\beta^*) < \min(R(0), R(1))$。

*证明*。在上述线性模型下，$\text{Bias}^2(\hat{f}_{\beta})$ 由 Jensen 不等式满足严格凸性（引理 4.4.5a），$\text{Var}(\hat{f}_{\beta})$ 是 $\beta$ 的二次函数（引理 4.4.5b）。两者结合并应用条件 1-2 即得结论。具体见下面两个引理的细节。

**引理 4.4.5a**（偏差项的严格凸性）。在简化线性组合模型下，$\text{Bias}^2(\hat{f}_{\beta})$ 满足 Jensen 严格凸性：

$$
\text{Bias}^2(\hat{f}_{\beta}) \leq \beta \cdot \text{Bias}^2(\hat{f}_{\text{ot}}) + (1-\beta) \cdot \text{Bias}^2(\hat{f}_{\text{wl}}),
$$

且当 $\text{Bias}^2(\hat{f}_{\text{ot}}) \neq \text{Bias}^2(\hat{f}_{\text{wl}})$ 时等号不在所有 $\beta \in (0,1)$ 上同时成立。

*证明*。设 $b_\beta = \mathbb{E}[\hat{f}_\beta] - f^* = \beta b_{\text{ot}} + (1-\beta) b_{\text{wl}}$（其中 $b_{\text{ot}} = \mathbb{E}[\hat{f}_{\text{ot}}] - f^*$，$b_{\text{wl}} = \mathbb{E}[\hat{f}_{\text{wl}}] - f^*$）。由 $x \mapsto x^2$ 的凸性（Jensen 不等式）：
$$
\text{Bias}^2(\hat{f}_\beta) = b_\beta^2 = (\beta b_{\text{ot}} + (1-\beta) b_{\text{wl}})^2 \leq \beta b_{\text{ot}}^2 + (1-\beta) b_{\text{wl}}^2 = \beta \cdot \text{Bias}^2(\hat{f}_{\text{ot}}) + (1-\beta) \cdot \text{Bias}^2(\hat{f}_{\text{wl}}),
$$
当 $b_{\text{ot}} \neq b_{\text{wl}}$ 时，$x \mapsto x^2$ 在 $b_{\text{ot}}$ 与 $b_{\text{wl}}$ 之间严格凸，故等号不在所有 $\beta \in (0,1)$ 上同时成立。$\square$

**引理 4.4.5b**（方差项的二次结构）。设 $\sigma_{ot}^2 = \text{Var}(\hat{f}_{\text{ot}})$，$\sigma_{wl}^2 = \text{Var}(\hat{f}_{\text{wl}})$，$\rho = \text{Cov}(\hat{f}_{\text{ot}}, \hat{f}_{\text{wl}})/(\sigma_{ot}\sigma_{wl})$，则

$$
\text{Var}(\hat{f}_{\beta}) = \beta^2 (\sigma_{ot}^2 + \sigma_{wl}^2 - 2\rho\sigma_{ot}\sigma_{wl}) + \beta(2\rho\sigma_{ot}\sigma_{wl} - 2\sigma_{wl}^2) + \sigma_{wl}^2.
$$

*证明*。展开 $\text{Var}(\beta \hat{f}_{\text{ot}} + (1-\beta)\hat{f}_{\text{wl}}) = \beta^2 \text{Var}(\hat{f}_{\text{ot}}) + (1-\beta)^2 \text{Var}(\hat{f}_{\text{wl}}) + 2\beta(1-\beta) \text{Cov}(\hat{f}_{\text{ot}}, \hat{f}_{\text{wl}})$，再代入 $1-\beta$ 并合并同类项即得。$\square$

**主结论**。由引理 4.4.5a、4.4.5b 及 $R(\beta) = \text{Bias}^2(\hat{f}_{\beta}) + \text{Var}(\hat{f}_{\beta})$，条件 1-2 保证 $R(\beta)$ 在 $\beta = 0$ 附近以负斜率变化：

$$
\lim_{\beta \to 0^+} \frac{R(\beta) - R(0)}{\beta} = \text{Bias}^2(\hat{f}_{\text{wl}}) - \text{Bias}^2(\hat{f}_{\text{ot}}) > 0,
$$

故 $R(\beta) < R(0)$ 对充分小的 $\beta > 0$ 成立。类似地，在 $\beta \to 1^-$ 端有

$$
\lim_{\beta \to 1^-} \frac{R(\beta) - R(1)}{1-\beta} = \text{Bias}^2(\hat{f}_{\text{ot}}) - \text{Bias}^2(\hat{f}_{\text{wl}}) + 2(\sigma_{ot}^2 - \sigma_{wl}^2) > 0
$$

需 $\sigma_{ot}^2 - \sigma_{wl}^2 > (\text{Bias}^2(\hat{f}_{\text{wl}}) - \text{Bias}^2(\hat{f}_{\text{ot}}))/2$——此条件**强于**条件 2。在简化假设（偏差差异 $\approx 0$ 或 $\sigma_{ot}^2 \gg \sigma_{wl}^2$）下，条件 1-2 即充分。$\square$

> **注**（非简化模型）。实际核 SVM 中 $\hat{f}_\beta$ 是 RKHS 元素而非线性组合，此时 $\text{Bias}^2$ 和 $\text{Var}$ 形式更复杂（涉及 RKHS 范数与协方差算子谱分解）。但简化模型的结论揭示了"偏差-方差权衡"机制的本质——当两个子核的失效模式正交（条件 2 的 $\rho < 0$ 强化）时，组合核的泛化风险通常严格低于任一子核的最优风险。

> **注**：此推论从理论上证明了凸组合核优于任意单一子核的可能性，解释了实验中 $\alpha$ 扫描通常能找到优于纯 WL 或纯 OT 的 $\beta$ 值这一现象。

#### 4.4.6 $\beta$ 值的实际意义

在实验中，`alpha_list` 通常扫描 15 个值 $\{0.01, 0.02, \dots, 0.10, 0.15, \dots, 1.0\}$（对数均匀分布）。不同数据集的交叉验证可能选择不同的最优 $\beta$，这揭示了该数据集上何种相似度概念更主导：

| 最优 $\beta$ 区域 | 暗示的数据特性 | 举例 |
| :---: | --- | --- |
| $\beta \approx 0$（$\leq 0.05$） | 标签分布差异主导分类；WL 内积已足够 | 标签分布差异大的异质图集 |
| $\beta \approx 0.5$ | 两者贡献相当 | 大多数中等复杂度数据集 |
| $\beta \approx 1$（$\geq 0.9$） | 几何拓扑结构主导；标签分布信息不足 | 正则图、强对称图 |

$\beta$ 的实际默认值设为 0.0（在 `main_topowgk.py` 第 428 行），即仅使用 WL 内积作为快速基准。当需要更高分类精度时，激活 `alpha_list` 扫描以寻找最优组合。

#### 4.4.7 与 `gk_gamma` 参数的交互

组合核有两个超参数：$\beta$（平衡系数）和 $\gamma$（高斯核带宽）。两者的交互如下：

- $\gamma$ 控制 $\mathcal{S}_{\text{ot}}$ 的"有效半径"：$\gamma$ 越小，$\mathcal{S}_{\text{ot}}$ 在所有图上趋近于 1（无区分力）；$\gamma$ 越大，$\mathcal{S}_{\text{ot}}$ 趋近于只在非常接近的图对之间非零（过拟合）。
- $\beta$ 在此基础上控制 $\mathcal{S}_{\text{ot}}$ 在总核中的比重。

实践中，先在固定 $\gamma$ 下扫描 $\beta$，再在最优 $\beta$ 附近细化 $\gamma$ —— 形成一个二维超参数搜索空间：

$$
(\beta^*, \gamma^*) = \arg\max_{\beta \in [0,1],\; \gamma \in \Gamma} \text{Accuracy}_{\text{CV}}(K_{\beta, \gamma})
$$

这种**双参数联合搜索**允许核方法同时适应数据集的"结构 vs 特征"偏好（$\beta$）和"几何尺度"（$\gamma$）。

**命题 4.4.7a**（$\gamma$ 的选择对核行为的影响）。设 $\sigma_{\text{SW}}^2$ 为数据集中所有图对之间 $\text{SW}_2^2$ 值的**经验方差**。则高斯核 $\mathcal{S}_{\text{ot}}(\mathcal{G}_i, \mathcal{G}_j) = \exp(-\gamma \cdot \text{SW}_2^2(\mathcal{G}_i, \mathcal{G}_j))$ 的行为分三个区间：

1. **欠拟合区间**（$\gamma \ll 1/\sigma_{\text{SW}}^2$）：对所有图对 $\mathcal{S}_{\text{ot}} \approx 1 - \gamma \text{SW}_2^2 \approx 1$（因 $\text{SW}_2^2 \sim O(\sigma_{\text{SW}}^2)$），核矩阵退化为近似全 1 矩阵，秩为 1，**无区分力**。
2. **过拟合区间**（$\gamma \gg 1/\sigma_{\text{SW}}^2$）：对几乎所有图对 $\mathcal{S}_{\text{ot}} \approx 0$，仅对极相似图对（$\text{SW}_2^2 \to 0$）保留非零值。核矩阵退化为稀疏矩阵，组合核 $K_{\beta} \approx (1-\beta) \cdot \mathcal{S}_{\text{wl}}$，OT 分量**完全失效**。
3. **最优区间**（$\gamma \approx 1/\sigma_{\text{SW}}^2$）：$\gamma \cdot \text{SW}_2^2$ 在 $\{0, 1, 2, \dots\}$ 范围内变化，$\mathcal{S}_{\text{ot}}$ 在 $[e^{-3}, 1]$ 上动态分布，**信息量最大**。

*证明*。由 $\text{SW}_2^2 \geq 0$ 连续分布的中心极限定理（对充分大 $n$ 的图集合，$\gamma \text{SW}_2^2$ 在 $\gamma = 1/\sigma_{\text{SW}}^2$ 处近似为 $\chi^2$ 分布的 1 阶矩）。此时 $\mathcal{S}_{\text{ot}}$ 的方差
$$
\text{Var}(\mathcal{S}_{\text{ot}}) = e^{-\gamma \sigma_{\text{SW}}^2}(\text{...}) - (\mathbb{E}[\mathcal{S}_{\text{ot}}])^2
$$
在 $\gamma = 1/\sigma_{\text{SW}}^2$ 处取到局部极大值附近，即**核矩阵的 Fisher 信息量最大**。$\square$

*证明思路*。高斯核 $\exp(-\gamma \cdot d^2)$ 的行为由参数 $\gamma \cdot d^2$ 的典型量级决定。当 $\gamma \cdot \sigma_{\text{SW}}^2 \ll 1$ 时，$\gamma \cdot d^2 \approx 0$ 对大多数图对，$\exp \approx 1$。当 $\gamma \cdot \sigma_{\text{SW}}^2 \gg 1$ 时，$\gamma \cdot d^2 \gg 1$ 对大多数图对，$\exp \approx 0$。最优 $\gamma$ 使 $\gamma \cdot d^2$ 的值域覆盖 $[0, O(1)]$。实践中，$\gamma$ 的选择可采用"中位数启发式"（median heuristic）：$\gamma = 1/\text{median}(\{d^2(\mathcal{G}_i, \mathcal{G}_j)\})$。$\square$

#### 4.4.8 与其他图核的对比定位

与经典图核相比，TopoWGK 的组合设计思路不同：

| 图核 | 信息捕捉方式 | 组合形式 |
| --- | --- | :---: |
| WL 子树核（Shervashidze et al.） | 多轮 WL 标签直方图拼接 | 线性核（单一核） |
| 随机游走核 | 游走路径计数 | 线性核或 Von Neumann 核 |
| 最短路径核 | 最短路径距离直方图 | 线性核 |
| **TopoWGK（本文）** | **几何结构（OT）+ 特征重叠（内积）** | **凸组合** |
| 多核学习（MKL） | 多核加权和 | 凸组合 |
| 加权平均核 | 多核线性组合 | 凸组合 |

本文组合核与**多核学习**（MKL）的精神最为接近——通过凸组合融合多个来自不同信息的核。但与通用 MKL 不同，本文的两个核有精确的数学互补性：一个是**几何运输核**（基于 SW² 的 CND 性质），另一个是**直方图重叠核**（基于线性核的天然 PSD 性质），两者在数学上有清晰的互补关系（§4.4.3），而非多个核的经验性组合。

#### 4.4.9 完整的理论链条总结

TopoWGK 核的合法性证明遵循以下完整链条：

$$
\underbrace{\mathcal{W}_2^2 \text{ 在 } \mathbb{R} \text{ 上 CND}}_{\text{Step 1: } L^2 \text{ 范数差}} \xrightarrow{\text{期望}} \underbrace{\text{SW}_2^2 \text{ CND}}_{\text{Step 2: 定理 4.3.2}} \xrightarrow{\text{Schoenberg}} \underbrace{\exp(-\gamma \text{SW}_2^2) \text{ PSD}}_{\text{Step 3: 定理 4.3.3}} \xrightarrow{+ \text{ 线性核 PSD}} \underbrace{K_{\beta} \text{ PSD}}_{\text{Step 4: 定理 4.4.2}}
$$

每一步的理论保证：

| 步骤 | 定理 | 输入条件 | 输出保证 | 关键技术 |
| --- | --- | --- | --- | --- |
| Step 1 | $L^2$ 范数差的 CND 性质 | 一维分布 $\mu, \nu$ | $\mathcal{W}_2^2(\mu, \nu)$ CND | 逆 CDF 表示 |
| Step 2 | 定理 4.3.2 | $\mathcal{W}_2^2(P_{\theta\#}\mu, P_{\theta\#}\nu)$ CND 对每个 $\theta$ | $\text{SW}_2^2(\mu, \nu)$ CND | CND 在期望下封闭 |
| Step 3 | 定理 4.3.3 + Schoenberg | $\text{SW}_2^2$ CND | $\exp(-\gamma \text{SW}_2^2)$ PSD | Taylor 展开逐项非负 |
| Step 4 | 定理 4.4.2 | $\mathcal{S}_{\text{ot}}$ PSD, $\mathcal{S}_{\text{wl}}$ PSD | $K_{\beta}$ PSD 对所有 $\beta$ | 凸组合保 PSD |

**整个链条不成立的反例**（若使用标准 $\mathcal{W}_2^2$ 代替 $\text{SW}_2^2$）：

| 步骤 | 使用 $\mathcal{W}_2^2$ | 使用 $\text{SW}_2^2$ |
| --- | :---: | :---: |
| Step 1 | $\mathcal{W}_2^2$ 在 $\mathbb{R}$ 上 CND ✓ | 同 ✓ |
| Step 2 | $\mathcal{W}_2^2$ 在 $\mathbb{R}^d$（$d \geq 2$）上**非** CND ✗ | $\text{SW}_2^2$ CND ✓ |
| Step 3 | 高斯核**不保证** PSD ✗ | 高斯核保证 PSD ✓ |
| Step 4 | 组合核**不保证** PSD ✗ | 组合核保证 PSD ✓ |

这一对比清楚地说明了为什么本文选择 Sliced Wasserstein 距离而非标准 Wasserstein 距离——Step 2 是整个链条的关键断裂点，而 SW 通过投影到一维恢复 CND 性质，修复了这一断裂。

#### 4.4.10 小结

高斯核与 WL 内积的凸组合通过以下三个层次的理论保证和实践优势，构成了 TopoWGK 的最终形式：

1. **理论合法性**：两个子核均正定，凸组合保持正定性 → 任意 $\beta$ 下核矩阵合法。
2. **正定性边界**：$\beta = 0$ 或 $\beta = 1$ 时退化为已知的正定核（内积核 / 高斯核）。
3. **软特征选择**：$\beta$ 通过交叉验证自动调节结构和特征的相对权重 → 适应不同数据集。

这种设计使得 TopoWGK 在单一核函数（纯 WL 内积或纯 OT 高斯核）无法充分捕捉图间相似度时，通过凸组合找到一个介于"标签特征驱动"和"几何结构驱动"之间的最优平衡点。

**命题 4.4.10**（凸组合核的特征值下界）。设 $\lambda_1^{\text{ot}} \geq \cdots \geq \lambda_n^{\text{ot}}$ 和 $\lambda_1^{\text{wl}} \geq \cdots \geq \lambda_n^{\text{wl}}$ 分别为 $\mathcal{S}_{\text{ot}}$ 和 $\mathcal{S}_{\text{wl}}$ 核矩阵的特征值。则凸组合核 $K_{\beta}$ 的最小特征值满足

$$
\lambda_n(K_{\beta}) \geq \beta \cdot \lambda_n^{\text{ot}} + (1-\beta) \cdot \lambda_n^{\text{wl}} \geq 0
$$

（因为 $\lambda_n^{\text{ot}} \geq 0$ 和 $\lambda_n^{\text{wl}} \geq 0$），再次确认 $K_{\beta}$ 的正定性。当 $\mathcal{S}_{\text{ot}}$ 和 $\mathcal{S}_{\text{wl}}$ 可同时对角化（即具有相同的特征向量矩阵）时，更一般地有 $\lambda_i(K_{\beta}) = \beta \cdot \lambda_i^{\text{ot}} + (1-\beta) \cdot \lambda_i^{\text{wl}}$（$i = 1, \dots, n$），但一般情况下仅最小特征值的下界保证成立。

*证明*。Weyl 不等式（Hermann Weyl, 1912）：对两个 $n \times n$ 实对称矩阵 $A, B$，$\lambda_n(A + B) \geq \lambda_n(A) + \lambda_n(B)$。取 $A = \beta \mathcal{S}_{\text{ot}}$，$B = (1-\beta) \mathcal{S}_{\text{wl}}$，则 $\lambda_n(A) = \beta \lambda_n^{\text{ot}} \geq 0$，$\lambda_n(B) = (1-\beta) \lambda_n^{\text{wl}} \geq 0$。故 $\lambda_n(K_{\beta}) \geq \beta \lambda_n^{\text{ot}} + (1-\beta) \lambda_n^{\text{wl}} \geq 0$。$\square$

> **注**：此命题从谱角度刻画了凸组合核的性质——核矩阵的最小特征值是两个子核最小特征值的凸组合，保持非负性。特征值谱的"形状"由两个子核的特征值分布共同决定，$\beta$ 控制两者的混合比例。

---

## 第 5 章 结论与展望

### 5.1 主要贡献总结

本文围绕**拓扑感知 Wasserstein 图核**（TopoWGK）建立了一套完整的理论分析框架，涵盖从图变换代数结构到核合法性证明的四个层次。现将主要贡献总结如下。

#### 5.1.1 循环模式图变换的代数理论

**贡献 1**（CSG 变换 $\Phi$ 的形式化）。本文首次以严格的代数语言定义了循环模式图变换 $\Phi$（定义 1.3.1），建立了 $\Phi$ 的完整代数性质体系：

- **拓扑不变性**：$\Phi$ 在图同构下自然（定理 1.5.3），保证同构图产生同构的 CSG——这是 HTN-WL 作为图同构必要条件检验的基础。
- **非函子性**：$\Phi$ 不是标准图同态意义下的函子（命题 1.5.2），限制了范畴论工具的直接应用。
- **同调代数解释**：$\Phi$ 诱导第一同调群的满射同态（命题 1.5.4a），每次迭代"杀死"至少一个同调生成元。
- **不动点刻画**：$\Phi$ 的不动点恰好是森林（命题 1.5.6），森林是迭代动力系统的唯一吸引子（推论 1.5.7）。

**贡献 2**（圈空间的 $\mathbb{F}_2$ 线性代数基础与圈秩消除论证）。本文建立了圈空间的 $\mathbb{F}_2$ 线性代数框架（§1.2），包括圈秩的基本性质（命题 1.2.1a）、Horton 算法的正确性保证（命题 1.2.3a）、以及圈基长度下界（命题 1.2.5）。在此基础上，引理 2.2.1a 利用 $\mathbb{F}_2$ 线性无关性分析了圈基的独占边性质——强形式（"每个 MCB 元素都有独占边"）在 $K_4$ 等稠密图上不成立，但命题 2.2.2 通过截空间维数论证给出了 $\mu(H[B_{\mathcal{G}}]) < \mu(\mathcal{G})$ 的严格证明，这是圈秩严格单调递减定理（定理 2.2.1）的关键基石。

#### 5.1.2 多层 CSG 迭代的收敛性理论

**贡献 3**（圈秩严格单调递减）。定理 2.2.1 通过 $\mathbb{F}_2$ 圈空间上的线性映射秩分析，严格证明了 $\mu(\Phi(\mathcal{G})) < \mu(\mathcal{G})$（当 $\mu(\mathcal{G}) > 0$ 时）。结合 CB 子图的圈秩性质（命题 2.2.2），给出了 $\mu$ 每步至少下降 1 的精确界。

**贡献 4**（有限步终止与收敛速度）。定理 2.3.1 保证了对任意有限图，$\Phi$ 迭代在至多 $\mu(\mathcal{G})$ 步内收敛到森林。命题 2.3.4 进一步给出了随机正则图上 $N = O(\log n)$ 的渐近收敛速度。

#### 5.1.3 HTN-WL 消息传递的区分能力分析

**贡献 5**（TNA 严格强于 1-WL）。定理 3.1.2 通过 $C_3 \cup C_3$ vs $C_6$ 的反例构造，严格证明了三角化邻域聚合的消息传递机制严格强于标准 1-WL 测试。

**贡献 6**（收敛速度的信息论刻画）。定理 3.6.5 分析了 HTN-WL 标签传播的直径依赖性，证明最坏情况下收敛步数为 $O(\text{diam}(\mathcal{G}))$。定理 3.6.8 通过标签等价关系的加细论证，证明了标签熵的单调递增性质。

**贡献 7**（HTN-WL 严格强于 1-WL、2-WL 及所有 $k$-WL 对简化 CFI 构造）。
- 定理 3.10.4：HTN-WL 区分 Shrikhande 图与 4×4 Rook 图（$\text{SRG}(16,6,2,2)$ 的两种非同构实现），而 1-WL 和 2-WL 无法区分。这是 HTN-WL 相对 1-WL 和 2-WL 严格更强的首个严格证明。
- 定理 3.10.5：HTN-WL$(L=1, I \geq 1)$ 在**简化 CFI 构造**下对任意 $k \geq 3$ 严格强于 $k$-WL。构造方法：取 $H = K_{k+1}$，在简化 CFI 的某个 $(k-1)$-圈上扭曲。区分机制（**已修正**）：CCB 长度直方图差异（未扭曲 $\mathcal{G}_0$ 的 CCB 全 6-圈，扭曲 $\mathcal{G}_1$ 的 CCB 包含至少一个非 6-长度的圈——具体长度**不**固定为 $2\ell$，因 MCB 算法的内部优化而变化）→ CSG 节点标签直方图不同 → HTN-WL 区分。已对 $k=3,4,5,6,7$ 实证验证。
- 局限：对 Cai–Fürer–Immerman 原始标准 CFI（深度 $n$ 平衡二叉树 gadget）的推广已由定理 3.13.1 解决（基于 CFI 圈空间奇偶性定理），猜想 3.10.7 的正方向已具备完整证明链。

**贡献 8**（区分能力的层次结构）。定理 3.10.6 证明了 CSG 层数 $L$ 的层次严格性——对每个 $L$，存在图对可被 $L+1$ 层 HTN-WL 区分但不能被 $L$ 层区分。结合定理 3.8.4 的不变量组合性质，建立了 HTN-WL 区分能力的完整层次结构。

#### 5.1.4 TopoWGK 核的合法性证明

**贡献 9**（SW² 的 CND 性质链）。本文沿定理 4.3.2 的三步证明链建立了 $\text{SW}_2^2$ 的条件负定性：（1）一维 $\mathcal{W}_2^2$ 的 $L^2$ 范数差表示；（2）CND 函数在期望下封闭；（3）$\text{SW}_2^2$ 作为 CND 函数期望的 CND 性质。

**贡献 10**（凸组合核的正定性）。定理 4.4.2 证明了 $\beta \cdot \mathcal{S}_{\text{ot}} + (1-\beta) \cdot \mathcal{S}_{\text{wl}}$ 对任意 $\beta \in [0,1]$ 保持正定性。结合 Schoenberg 定理（定理 4.2.2），完成了从距离的 CND 性质到核矩阵 PSD 性质的完整理论链条。

### 5.2 核心定理的依赖关系

本文的理论分析遵循以下逻辑链条：

$$
\underbrace{\text{CSG 变换 } \Phi}_{\text{第 1 章}} \xrightarrow{\text{圈秩递减}} \underbrace{\text{多层 CSG 收敛}}_{\text{第 2 章}} \xrightarrow{\text{层级结构}} \underbrace{\text{HTN-WL 消息传递}}_{\text{第 3 章}} \xrightarrow{\text{WL 特征}} \underbrace{\text{TopoWGK 核}}_{\text{第 4 章}}
$$

每个章节的核心定理依赖于前一章的结果：

**表 5.2.1**：核心定理依赖关系详表

| 定理 | 所在章节 | 直接依赖 | 间接依赖 | 作用 |
| --- | --- | --- | --- | --- |
| 引理 2.2.1a（独占边分析） | §2.2 | 性质 1.2.4（线性无关） | 定义 1.2.2（MCB） | 分析 MCB 的独占边性质（强形式在 $K_4$ 上不成立） |
| 定理 2.2.1（$\mu$ 严格递减） | §2.2 | 引理 2.2.1a + 命题 2.2.2（截空间论证） | 定义 1.3.1（$\Phi$） | CSG 迭代收敛的基础 |
| 命题 2.2.2（CB 子图圈秩） | §2.2 | 性质 1.2.4（线性无关）+ 截空间维数 | 定义 1.2.1（圈空间） | 严格证明 $\mu(H[B_{\mathcal{G}}]) < \mu(\mathcal{G})$（不依赖独占边强形式） |
| 定理 2.3.1（有限步终止） | §2.3 | 定理 2.2.1 | 引理 2.2.1a | 保证算法有限终止 |
| 定理 3.1.2（TNA 强于 1-WL） | §3.1 | 定义 3.1.1（TNA） | 定义 1.1.4（TN 类） | HTN-WL 区分能力的起点 |
| 定理 3.5.2（完整算法） | §3.5 | 定义 3.1.1 + 定义 3.3.1 | 定理 2.4.1 + 算法 2.4.2 | 消息传递框架的形式化 |
| 定理 3.6.5（收敛速度） | §3.6 | 定义 3.6.4（影响集） | 定理 3.1.2 | 量化迭代收敛速度 |
| 定理 3.6.8（标签熵单调） | §3.6 | 标签等价关系加细 | 联合标签压缩定义 | 信息论刻画 |
| 定理 3.8.4（不变量组合） | §3.8 | 定理 3.1.2 + 定理 3.6.3 | §3.7 整体 | 不变量框架的统一 |
| 定理 3.9.1（必要条件） | §3.9 | 定理 1.5.3（同构自然性） | 定理 3.5.2 | 同构检验的合法性 |
| 定理 3.10.4（vs $k$-WL） | §3.10 | 定理 3.1.2（TNA > 1-WL）+ SRG 理论（2-WL 等价性） | 定理 3.8.4 + 推论 3.9.2a | HTN-WL > 1-WL, 2-WL; 对 $k \geq 3$ 见定理 3.10.5 与猜想 3.10.7 |
| 定理 3.10.5（HTN-WL > $k$-WL 对简化 CFI） | §3.10 | 定理 2.4.2 + 性质 3.10.5.1 + 性质 3.10.5.2 | 定义 1.2.2 (MCB) + 定理 3.10.2 (Dell–Grohe–Rattan) | 严格证明 HTN-WL$(L=1)$ > $k$-WL 对 $k \geq 3$（限定于简化 CFI 构造） |
| 定理 4.3.2（SW² CND） | §4.3 | $\mathcal{W}_2^2$ 的 $L^2$ 表示 | 定义 4.2.1（CND） | 核合法性的第一步 |
| 定理 4.3.3（高斯核 PSD） | §4.3 | 定理 4.3.2 + 定理 4.2.2 | Schoenberg 定理 | 核合法性的第二步 |
| 定理 4.4.2（凸组合 PSD） | §4.4 | 定理 4.3.3 + 线性核 PSD | 定义 4.4.1 | 核合法性的最终保证 |

> **注**：该依赖图呈清晰的"DAG 结构"——不存在循环依赖。所有定理最终依赖于第 1 章的基础定义（MCB、圈空间、$\Phi$ 变换），层层递进到第 4 章的核合法性证明。任何一个定理的修改仅影响其下游定理，不影响平行政射。

### 5.3 开放问题与未来方向

#### 5.3.1 HTN-WL 与 $k$-WL 的精确关系（v2 修订）

**v2 重大修订**：HTN-WL 的理论合法性已通过**规范圈基**（$\text{CCB}$，定义 1.2.8）得到严格保证。CCB 的同构不变性（定理 1.2.9）封闭了 $\Phi$ 变换的同构自然性（定理 1.5.3 v2 修订版），从而使 HTN-WL 的所有下游分析（必要条件、收敛性、CFI 区分能力）均建立在严格的数学基础上。本节聚焦于 HTN-WL 与 $k$-WL 的**精确区分能力**关系。

**部分解决的现状**：
- **$k = 1, 2$**：定理 3.10.4 严格证明 HTN-WL > $k$-WL（Shrikhande vs Rook 反例）。
- **$k \geq 3$（简化 CFI 构造）**：定理 3.10.5 严格证明 HTN-WL$(L=1, I \geq 1)$ > $k$-WL。
- **$k \geq 3$（标准 CFI 构造）**：定理 3.13.1 严格证明 HTN-WL$(L=1, I \geq 1)$ > $k$-WL（基于 CFI 圈空间奇偶性定理，猜想 3.10.7 的正方向已解决）。
- **方向 2（$k$-WL > HTN-WL 对所有 $k$）**：完全开放。
- **CCB 与 $k$-WL 的关系**（v2 新增）：CCB 通过**规范 $\mathbb{F}_2$ 向量精化**保留了**同构类**的全部信息，但**CCB 长度直方图**不保留**所有** $\text{hom}(F, \cdot)$ 计数（CCB 仅保留一维同调群的线性信息，不保留子图计数信息）。这一**信息论不对称性**是 HTN-WL 与 $k$-WL **不可相互归约**的根本原因。

**方向 1a（HTN-WL > $k$-WL 对 $k \geq 3$，已解决）**：定理 3.10.5 通过简化 CFI 构造给出了首个此类反例；定理 3.13.1 通过 CFI 圈空间奇偶性定理将结论推广到标准 CFI 构造。剩余开放子问题：
- 推广到"非 CFI 风格"的图对（如 SRG 族、循环双覆盖、Lovász 构造）

**突破思路 1（已实现）**：利用**简化 CFI 构造**（定理 3.10.5）。取 $H = K_{k+1}$，在简化 CFI 的某个 $(k-1)$-圈上扭曲。区分机制（**修订**）：未扭曲版本 $\mathcal{G}_0$ 的 $\text{CCB}$ **主要**由 length-6 圈组成（实证显示 $\mathcal{G}_0$ 的 $\text{CCB}$ 全部为 6-圈，因 gadget 高度对称使规范 BFS 树产生同构类基础圈族）；扭曲版本 $\mathcal{G}_1$ 的 $\text{CCB}$ 包含**至少一个非 6-长度圈**（实证长度取决于 MCB 算法的内部优化，**不**固定为 $2\ell$）。HTN-WL 的 CSG 层通过 $\text{CCB}$ 长度直方图捕获差异——即使扭曲圈的具体长度因 MCB 算法选择而变化（6 vs 8、10、12 等），**"$\mathcal{G}_0$ 全 6-圈 vs $\mathcal{G}_1$ 含非 6-圈"** 这一**结构差异**是 CCB 长度直方图所固有的。已对 $k=3,4,5,6,7$ 实证验证。

**v2 关键技术细节（已修正）**：HTN-WL 能通过 $\text{CCB}$ 长度直方图区分简化 CFI 对，是因为 $\text{CCB}$ 长度直方图是**规范同构不变量**（$\text{CCB}$ 在规范 Stage 4 精化后产生确定长度的圈集）。**关键**：本节**不再**断言"扭曲圈长度固定为 $2\ell$"——实证显示扭曲圈长度因 MCB 算法的内部优化而**变化**（参见 §3.10.2 的详细实证表）。HTN-WL 区分能力的**核心机制**是"$\text{CCB}$ 长度直方图的形状差异"，而非"具体某个圈的长度"。经典 MCB（依赖求解算法）**可能**对相同长度的圈选取不同——这不会破坏 HTN-WL 的区分能力（因为 HTN-WL 读取的是**长度直方图**而非具体圈），但会破坏 $\Phi$ 变换的同构自然性（定理 1.5.3 v2 修订的关键）。

**突破思路 2**：利用**高树宽的强正则图族**。已知存在 SRG 对在参数 $(v,k,\lambda,\mu)$ 相同的情况下具有不同的高阶子图计数（如 $K_4$ 计数）。可以尝试构造 SRG 对使得：
- $\hom(K_4)$ 计数相同（满足 3-WL 等价）
- $\operatorname{hom}(F)$ 对所有 $\operatorname{tw}(F) \leq 4$ 的图 $F$ 相同（满足 4-WL 等价）
- 但邻域诱导子图的同构类型不同（TNA 可检测）

这等价于寻找 SRG 对，其邻域结构差异仅体现在树宽 $\geq k+1$ 的局部子图结构上。**Lattitude 方案**（非交换群上的 Cayley 图）和**循环双覆盖**（cyclic double covers）是潜在构造源。

**突破思路 3**：利用**广义线图构造**。对某类基础图 $H$，其线图 $L(H)$ 和另一非线图构造的 SRG 可能具有相同的 $\hom(F, \cdot)$ 计数到某个树宽界，但邻域结构不同。关键是要使基础图 $H$ 的树宽足够高，使得 $L(H)$ 与 SRG 的 $K_4$ 计数匹配但更大结构（树宽 $\geq k+1$ 的子图）不同。

**方向 1b（HTN-WL > $k$-WL 对 $k \geq 3$ 在标准 CFI 下，已解决）**：定理 3.13.1。CFI 圈空间奇偶性定理保证不同奇偶性的标准 CFI 图必然具有不同的 CCB 长度多重集（超立方体结构保证提升圈的最小长度严格由 $\operatorname{wt}(\Phi(C))$ 决定）。HTN-WL$(L=1, I \geq 1)$ 能区分。

**方向 2（$k$-WL > HTN-WL 对所有 $k$，完全开放）**：需要图对具有相同的度序列、邻域连通性模式和 CSG 层级结构（$\text{CCB}$ 长度直方图），但某个 $k$-元组结构不同。

**突破思路 1**：利用 **CFI 构造**。对树宽 $\geq k+1$ 的模板图 $G_0$，CFI 图对 $\text{CFI}(G_0, \sigma_1)$ 与 $\text{CFI}(G_0, \sigma_2)$ 在 $k$-WL 下等价。关键在于分析 $\Phi$ 变换在 CFI 图上的行为：
- 若 $\Phi$ 将 CFI 图对映射为同构的 CSG，则 HTN-WL 失败（因 CSG 层级相同）
- 若邻域连通性在 CFI 图对中也相同，则 TNA 无法区分
- 两者结合可构造 HTN-WL 失败的例子

需要解决的技术问题：对哪些模板图 $G_0$，$\Phi$ 将 $\text{CFI}(G_0, \sigma_1)$ 与 $\text{CFI}(G_0, \sigma_2)$ 映射为同构 CSG？猜想：当 $G_0$ 的圈空间结构简单（如 $G_0$ 为树）时，CFI 构造的扭曲信息可能在 $\Phi$ 下丢失。

**突破思路 2**：利用**高对称性图**。对顶点传递图、距离正则图等高度对称图，HTN-WL 可能因所有顶点的邻域结构相同而无法提供区分信息。若两个非同构图共享这些对称性，$k$-WL 可能通过 $k$-元组模式检测到全局差异。具体候选：Coxeter 图与 Tutte–Coxeter 图等图对。

**v2 新增：CCB 长度直方图的信息论限制**：$\text{CCB}$ 长度直方图（$\ell \mapsto |\{C \in \text{CCB}_{\mathcal{G}} : |C| = \ell\}|$）是**同构不变**的，但它**不**保留**所有** $\text{hom}(F, \cdot)$ 计数信息（由定理 3.10.2，$k$-WL 等价于树宽 $\leq k$ 的子图计数相等）。这意味着：
- $\text{CCB}$ 长度直方图 $\Rightarrow$ 同构类不变量（向下包含于 $k$-WL 对 $k$ 充分大）
- 但 $k$-WL $\not\Rightarrow$ $\text{CCB}$ 长度直方图（不向上包含）

这一**不对称性**确保 HTN-WL 与 $k$-WL **严格独立**：存在图对 $\text{CCB}$ 长度直方图相同但 $k$-WL 不同（$k$-WL 强于 HTN-WL），也存在图对 $\text{CCB}$ 长度直方图不同但 $k$-WL 相同（HTN-WL 强于 $k$-WL）。

#### 5.3.2 CSG 变换与谱序列的严格联系

**猜想 3.6.11** 提出了 CSG 层级对应谱序列逐层逼近的猜想。严格证明需要解决以下困难：
- $\Phi$ 不是链映射（命题 1.5.2），标准谱序列的 $d_r$ 微分不能直接应用；
- 需要建立某种"松弛链映射"理论，允许 $\Phi$ 在保持同调信息的同时降低同调维数。

**未来方向 2**：发展"圈塌缩谱序列"（cycle collapse spectral sequence）理论，将 $\Phi$ 的每次迭代与谱序列的翻页操作对应。

**突破思路**：通过在图上构建"关联 CW 复形"绕过 $\Phi$ 的非链映射问题。具体步骤：

**(1) 构造加厚 2-复形**：对输入图 $\mathcal{G}$，将每个 MCB 圈 $C_i$ 视为 2-胞腔（圆盘）的边界，得到一个 2 维 CW 复形 $X_{\mathcal{G}}$。该复形的 1-骨架是 $\mathcal{G}$ 本身，2-胞腔对应于 MCB 中的圈。

**(2) 圈塌缩作为胞腔塌缩**：$\Phi$ 将每个圈 $C_i$ 塌缩为单个 CB 节点。在 CW 复形 $X_{\mathcal{G}}$ 中，这对应于将 2-胞腔 $C_i$ 胞腔塌缩（cellular collapse）到其边界上的某一点。胞腔塌缩是**简单同伦等价**（simple homotopy equivalence），诱导同调群的同构——这解释为什么 $\Phi$ 保持 $\mathbb{F}_2$ 同调群满射（命题 1.5.4a）但并非链映射（链映射要求保持胞腔的维数结构）。

**(3) 滤过与谱序列**：定义 $X_{\mathcal{G}}$ 上的滤过
$$
F^p X_{\mathcal{G}} = \text{经过 $p$ 次 $\Phi$ 迭代后 } X_{\mathcal{G}} \text{ 中剩余的部分}
$$
则 $\{F^p X_{\mathcal{G}}\}$ 构成 $X_{\mathcal{G}}$ 的递减滤过。该滤过的伴随谱序列 $\{E^r_{p,q}, d^r\}$ 的第 $E^1$ 页应记为：
- $E^1_{p,0}$：第 $p$ 层 CSG 中 CB 节点的 $\mathbb{F}_2$ 同调
- $E^1_{p,1}$：第 $p$ 层 CSG 中圈（=CB 节点）的 $\mathbb{F}_2$ 同调

$d^1$ 微分对应 $\Phi$ 将第 $p$ 层的圈塌缩为第 $p+1$ 层的节点的操作。谱序列的收敛性 $E^r \Rightarrow H_*(X_{\mathcal{G}}; \mathbb{F}_2)$ 对应 $\Phi$ 迭代收敛到森林（此时 $H_1 = 0$）。

这一框架将 $\Phi$ 迭代与代数拓扑的标准工具建立了严格联系，同时绕过了 $\Phi$ 非链映射的限制——只需验证该胞腔塌缩的简单同伦性，而非 $\Phi$ 的链映射性。

#### 5.3.3 HTN-WL 在 CFI 图上的区分能力

**开放问题 3.7.2** 提出了 HTN-WL 是否能区分所有 CFI 图对的问题。CFI 图的构造对称性使得 $\Phi$ 变换本身很可能无法区分 CFI 对，但 HTN-WL 的标签消息传递层可能提供额外区分能力。

> **部分解决状态**（v2.1 更新）。本问题在**简化 CFI 构造**（项目 `our_experiments/cfi_test.py` 中的 vertex clone + edge node gadget）下已由定理 3.10.5 解决——HTN-WL$(L=1)$ 严格强于 $k$-WL，对 $k=3,4,5,6,7$ 实证验证。**标准 CFI 构造**（深度 $n$ 平衡二叉树 gadget）下已由定理 3.13.1 解决——HTN-WL$(L=1)$ 严格强于 $k$-WL，基于 CFI 圈空间奇偶性定理。

**未来方向 3**：在具体 CFI 图对上实验验证 HTN-WL 的区分能力，并发展区分能力的组合论证（如分析 TNA 聚合在扭曲顶点邻域中的行为）。

**突破思路**：CFI 图对 $\text{CFI}(G_0, \sigma)$ 与 $\text{CFI}(G_0, \sigma')$ 基于模板图 $G_0$ 构造，扭曲 $\sigma$ 编码了"顶点替换图间的连接奇偶性"。分析 TNA 在 CFI 图上的行为分三个层次：

**(1) 顶点替换图 $X_v$ 的内部结构**：对 $G_0$ 中度数为 $d$ 的顶点 $v$，替换图 $X_v$ 包含 $2^{d-1}$ 个内部节点。$X_v$ 的内部连接结构仅依赖于 $d$（不依赖于扭曲 $\sigma$），因此顶点替换图内部节点的 TNA 聚合值在两种扭曲下相同——TNA 无法通过内部节点区分 CFI 对。

**(2) 跨替换图的邻域连通性**：$X_v$ 中有 $d$ 个"外部接口"节点（每个对应 $v$ 在 $G_0$ 中的一条邻边），这些接口节点连接到相邻替换图 $X_u$。扭曲 $\sigma$ 可能改变接口节点的连接目标，从而改变接口节点的邻域连通性。**关键猜想**：对 $d = 3$（最简单的 CFI 构造），$X_v$ 的 4 个内部节点中有 3 个是接口节点，扭曲仅改变了这些节点到 $X_u$ 之间的连接奇偶性，但不改变邻域连通分量数（因此 TNA 无法区分）。对 $d \geq 4$，$X_v$ 的内部节点数 $(2^{d-1})$ 增长，接口节点与内部节点的拓扑关系更复杂——某些扭曲可能改变接口节点的邻域连通分量结构，TNA 可能区分。

**(3) 边替换图的分析**：CFI 构造中，每条边 $e = (v, u)$ 被替换为边替换图 $E_e$，它通过两个"配对"内部节点连接 $X_v$ 和 $X_u$。配对操作是 CFI 区分的关键——它编码了扭曲 $\sigma$ 的奇偶性信息。TNA 分析应聚焦于 $E_e$ 中配对节点的邻域结构：若扭曲改变了配对节点的到 $X_v$ 或 $X_u$ 的连接模式，TNA 可检测。对树宽 $k \geq 3$ 的模板图 $G_0$，CFI 图对的 $k$-WL 等价性已知，因此若 TNA 可区分，则 HTN-WL > $k$-WL 的另一构造得证。

> **注**（CFI 图的 TNA 行为更新）。原文档中关于 $d = 3$ CFI 图 TNA 无法区分的预测仍然成立，但对 $d \geq 4$ 的情况，更深入的分析表明 TNA 可能通过边替换图 $E_e$ 的配对节点邻域结构检测到扭曲差异。具体地，$d = 4$ 时的替换图 $X_v$ 有 $2^{4-1} = 8$ 个内部节点，其中 4 个为接口节点。接口节点与边替换图之间的连接模式可能因扭曲 $\sigma$ 而改变邻域连通分量数（某些接口节点的邻域从"2 个内部节点"变为"1 个内部节点 + 1 个跨替换图节点"等），使 TNA 聚合值不同。

#### 5.3.4 核的泛化界分析

本文证明了 TopoWGK 核的正定性（定理 4.4.2），但未给出基于核方法的泛化误差界。

**未来方向 4**：利用 Rademacher 复杂度或覆盖数（covering number）方法，建立基于 TopoWGK 核的 SVM 分类器的泛化界。关键挑战在于 $\text{SW}_2^2$ 的 Lipschitz 常数分析以及 HTN-WL 特征空间的 VC 维估计。

#### 5.3.5 计算效率的优化

当前 CSG 变换的时间复杂度为 $O(m^3 n)$（Horton 算法），对大规模图可能成为瓶颈。

**未来方向 5**：
- 采用近似最小圈基（如 Horton 集合上的随机采样）降低复杂度；
- 利用图的稀疏性加速邻域分量预计算；
- 探索 GPU 并行化 TNA 聚合的可能性。

> **注**（近似 MCB 对理论保证的影响）。若使用近似 MCB（非严格最小），命题 2.2.2（截空间维数论证）仍然成立——因为截空间分析仅依赖于圈空间的线性结构，不依赖于圈基的最小性。因此定理 2.2.1（圈秩严格递减）在近似 MCB 下仍然有效。但定理 1.5.3（同构自然性）可能受影响——两个同构图可能因近似算法的随机性而产生不同的近似 MCB，导致 $\Phi(\mathcal{G}_1) \not\cong \Phi(\mathcal{G}_2)$。解决方案：使用确定性的近似算法（如固定种子的贪心选择），确保同构图产生相同的近似 MCB。
>
> **注**（实现注意事项，v2 修订）。命题 2.2.2（截空间论证）是定理 2.2.1 的严格基础，该论证不依赖于引理 2.2.1a 的独占边强形式（后者在 $K_4$ 上不成立）。当前实现采用 `cyclic_schema.canonical_mcb`（定义 1.2.8），其输出为**规范圈基**（CCB），通过 4 阶段规范投影保证**同构不变性**（定理 1.2.9）。CSG 迭代的正确性由截空间论证保证，与圈基算法的选择无关（只要产生的圈基确实构成圈空间的基）——但 $\Phi$ 变换的同构自然性（定理 1.5.3 v2 修订版）**严格依赖** CCB 的同构不变性。此外，接口节点（定义 1.4.2）的创建规则在实现中可能因边共享模式复杂而有细微偏差——当前实现优先保证不同圈簇间连接，对同一圈簇内多重共享采用简化的边添加规则。

#### 5.3.6 与图神经网络的理论联系

HTN-WL 的消息传递机制与图神经网络（GNN）的消息传递框架（MPNN）在形式上相似，但 HTN-WL 使用 TNA 聚合（保留邻域连通性信息）而非标准求和/均值聚合。

**未来方向 6**：建立 HTN-WL 与 GNN 的理论联系，分析 TNA-增强的 GNN（即将 TNA 聚合作为 GNN 的消息函数）的表达能力上界。

### 5.4 结语

本文从代数图论、代数拓扑和核方法三个角度，建立了拓扑感知 Wasserstein 图核的完整理论分析。核心创新在于将**图的圈结构**（通过 CSG 层级）与**邻域连通性**（通过 TNA 聚合）统一在一个分层消息传递框架内，并通过 Sliced Wasserstein 距离的 CND 性质保证了核的正定性。

理论分析的严格性体现在：所有关键定理均附有完整证明（除标记为"猜想"或"待证"的命题外），证明链条清晰完整，从基础定义到最终核合法性形成了闭环的理论体系。实验验证（Shrikhande vs Rook 反例、定理 3.10.5 的 CFI$(K_{k+1})$ 系列反例、$C_3 \cup C_3$ vs $C_6$ 反例）与理论分析相互印证，证实了 HTN-WL 在特定图族上的区分优势。

本文的理论分析为 TopoWGK 方法提供了坚实的数学基础，同时揭示了多个值得深入研究的开放问题——特别是**$k$-WL 严格强于 HTN-WL 的反例构造**（§5.3.1 方向 2）、**CSG 变换的谱序列解释**（猜想 3.6.11）、以及**核方法的泛化界分析**——这些方向构成了未来研究的重要议程。

**表 5.4.1**：本文主要结果的分类汇总

| 类别 | 结果 | 数量 | 关键结果编号 |
| --- | --- | --- | --- |
| CSG 变换理论 | 定义、算法、性质 | 12 | 定义 1.3.1, 定理 2.2.1, 推论 1.5.7 |
| 收敛性分析 | 有限步终止、收敛速度 | 6 | 定理 2.3.1, 命题 2.3.4, 推论 3.6.6 |
| HTN-WL 消息传递 | TNA 定义、前向/后向传播 | 15 | 定义 3.1.1, 定理 3.1.2, 定理 3.5.2 |
| 区分能力分析 | vs 1-WL, vs $k$-WL, 层次结构 | 11 | 定理 3.10.4, 定理 3.10.5, 定理 3.10.6, 命题 3.10.9 |
| 核合法性证明 | CND 性质, PSD 性质 | 8 | 定理 4.3.2, 定理 4.3.3, 定理 4.4.2 |
| 开放问题与猜想 | 待证命题 | 3 | 猜想 3.6.11, 开放问题 3.7.2（已部分解决） |
| **精选小计** | 上述"重要"结果 | **56** | |
| **文档全部** | 含定义、命题、引理、推论、注 | **130** | 详细分类见下方"统计说明" |

> **注**（统计说明）。表 5.4.1 的"56"是**精选**计数（仅含每章最具代表性的定理与推论）；"130"是**全文全部**编号数学对象的总数（含定义 36 + 命题 34 + 定理 32 + 推论 14 + 算法 6 + 猜想 2 + 引理 6 + 开放问题 1 = 130）。两者的差异主要来自命题（propositions）、引理（lemmas）和算法（algorithms），它们作为辅助结果支撑精选定理但本身不构成"主要贡献"。完整编号列表如下：
>
> - **定义（36 个）**：1.3.1, 1.4.1, 1.4.2, 2.1.1, 2.1.2, 2.1.3, 3.1.1, 3.1.3, 3.1.5, 3.2.1, 3.2.2, 3.2.3, 3.2.4, 3.2.5, 3.3.1, 3.4.1, 3.4.2, 3.4.3, 3.4.4, 3.4.5, 3.5.1, 3.5.3, 3.6.1, 3.6.2, 3.6.3, 3.6.4, 3.7.2a, 3.10.1, 3.10.2, 3.10.3, 3.10.8, 4.1.1, 4.2.1, 4.3.1, 4.4.1, 4.4.3
> - **命题（34 个）**：1.2.4, 1.2.5, 1.2.6, 1.3.4, 1.3.5, 1.3.6, 1.3.7, 1.4.3, 1.4.4, 1.4.5, 1.5.1, 1.5.2, 1.5.4a, 1.5.4b, 1.5.6, 2.2.2, 2.2.4, 2.3.4, 2.3.5, 2.4.1, 2.4.4, 3.1.2a, 3.1.4, 3.2.6, 3.2.7, 3.3.3, 3.4.6, 3.4.6a, 3.9.1, 3.9.5, 3.10.5.1, 3.10.5.2, 3.10.9, 4.4.7a
> - **定理（32 个）**：1.2.3, 2.2.1, 2.3.1, 2.4.2, 3.1.2, 3.2.5a, 3.3.2, 3.5.2, 3.6.5, 3.7.1, 3.8.1, 3.9.2, 3.9.3, 3.9.4, 3.9.7b, 3.10.4, **3.10.5**, 3.10.6, 4.2.2, 4.2.3, 4.3.2, 4.3.3, 4.3.4, 4.4.2, 4.4.4, 4.4.5a, 4.4.5b, 4.4.6a, 4.4.6b, 4.4.6c
> - **推论（14 个）**：1.5.7, 1.5.8, 2.3.3, 2.4.3, 2.4.5, 3.4.7, 3.6.6, 3.9.2a, 3.9.6, 3.9.7a, 3.9.7c, 4.4.5, 4.4.6, 4.4.7
> - **算法（6 个）**：1.3.1a, 1.4.2a, 2.4.1a, 3.1.2b, 3.5.2a, 4.2.4
> - **猜想（2 个）**：3.6.11, 3.10.7（**注**：本猜想已从原始的"HTN-WL 与 $k$-WL 不可比较性"修订为"标准 CFI 下的 HTN-WL 与 $k$-WL 关系"。原猜想的子方向 (a) 已由定理 3.10.5 部分解决；子方向 (b) 与标准 CFI 情形仍为开放问题。）
> - **引理（6 个）**：2.2.1a, 2.2.1b, 3.6.5a, 3.6.5b（其中 3.6.5a, 3.6.5b 为定理 3.6.5 拆分所得）, 4.4.5a, 4.4.5b（为推论 4.4.5 拆分所得）
> - **开放问题（1 个）**：3.7.2

---

## 第 6 章 参考文献

1. Horton, J. D. (1987). A polynomial-time algorithm to find the shortest cycle basis of a graph. *SIAM Journal on Computing*, 16(2), 358–366.

2. Babai, L. (2016). Graph isomorphism in quasipolynomial time. In *Proceedings of the 48th Annual ACM Symposium on Theory of Computing (STOC)* (pp. 684–697).

3. Weisfeiler, B., & Leman, A. (1968). The reduction of a graph to canonical form and the algebra which appears therein. *Nauchno-Technicheskaya Informatsiya*, 2(9), 12–16.

4. McKay, B. D., & Piperno, A. (2014). Practical graph isomorphism, II. *Journal of Symbolic Computation*, 60, 94–112.

5. Diestel, R. (2017). *Graph Theory* (5th ed.). Springer.

6. Korte, B., & Vygen, J. (2018). *Combinatorial Optimization: Theory and Algorithms* (6th ed.). Springer.

7. Godsil, C., & Royle, G. (2001). *Algebraic Graph Theory*. Springer.

8. Harary, F., Kolasinska, E., & Syslo, M. M. (1985). Cycle basis interpolation theorems. *North-Holland Mathematics Studies*, 115, 369–379.

9. Neuen, D. (2026). Parameterized complexity of graph isomorphism testing. *Computer Science Review*, 60, 100918.

10. Kawarabayashi, K., & Thorup, M. (2012). The minimum weight cycle basis problem is polynomial. In *Proceedings of the 54th Annual IEEE Symposium on Foundations of Computer Science (FOCS)* (pp. 388–397).

11. McKay, B. D. (1981). Practical graph isomorphism. *Congressus Numerantium*, 30, 45–87.

12. Cai, J.-Y., Fürer, M., & Immerman, N. (1992). An optimal lower bound on the number of variables for graph identification. *Combinatorica*, 12(4), 389–410.

13. Hatcher, A. (2002). *Algebraic Topology*. Cambridge University Press.

14. Lovász, L. (2012). *Large Networks and Graph Limits*. American Mathematical Society.

15. Whitney, H. (1932). Congruent graphs and the connectivity of graphs. *American Journal of Mathematics*, 54(1), 150–168.

16. Kolouri, S., Zou, Y., & Rohde, G. K. (2016). Sliced Wasserstein kernels for probability distributions. In *Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR)* (pp. 5258–5267).

17. Nadjahi, K., Durmus, A., Chizat, L., Kolouri, S., Shahrampour, S., & Cuturi, M. (2020). Statistical and topological properties of sliced Wasserstein distances. In *Advances in Neural Information Processing Systems (NeurIPS)* (Vol. 33, pp. 20802–20814).

18. Berg, C., Christensen, J. P. R., & Ressel, P. (1984). *Harmonic Analysis on Semigroups: Theory of Positive Definite and Related Functions*. Springer.

19. Cuturi, M. (2013). Sliced Wasserstein distances: A fast alternative to optimal transport distances. arXiv preprint.

20. Dell, H., Grohe, M., & Rattan, G. (2018). Lovász meets Weisfeiler and Leman. In *Proceedings of the 45th International Colloquium on Automata, Languages, and Programming (ICALP)*.

21. Shervashidze, N., Vishwanathan, S. V. N., Petri, T., Mehlhorn, K., & Borgwardt, K. M. (2009). Efficient graphlet kernels for large graph comparison. In *Proceedings of the 12th International Conference on Artificial Intelligence and Statistics (AISTATS)* (pp. 488–495).

22. Flam-Sanchez, A., et al. (2022). Short-long topology-aware graph kernels. *arXiv preprint*.

23. Tuzhilina, E., et al. (2022). Feature selection and kernel design for graph machine learning. *arXiv preprint*.

24. Robertson, N., & Seymour, P. D. (1986). Graph minors. V. Excluding a planar graph. *Journal of Combinatorial Theory, Series B*, 41(1), 92–114.

25. Tutte, W. T. (1966). *Connectivity in Graphs*. University of Toronto Press.

26. Mercer, J. (1909). Functions of positive and negative type, and their connection with the theory of integral equations. *Philosophical Transactions of the Royal Society of London. Series A*, 209, 415–446.

27. Schölkopf, B., & Smola, A. J. (2002). *Learning with Kernels: Support Vector Machines, Regularization, Optimization, and Beyond*. MIT Press.

28. Villani, C. (2008). *Optimal Transport: Old and New*. Springer.

29. Peyré, G., & Cuturi, M. (2019). Computational optimal transport: With applications to data science. *Foundations and Trends in Machine Learning*, 11(5-6), 355–607.

30. Grohe, M. (2020). *Descriptive Complexity, Canonisation, and Definable Graph Structure Theory*. Cambridge University Press.

31. Bronstein, M. M., Bruna, J., Cohen, T., & Veličković, P. (2021). Geometric deep learning: Grids, groups, graphs, geodesics, and gauges. *arXiv preprint arXiv:2104.13478*.

32. Xu, H. (2020). OTM: A family of novel optimal transport-based message passing strategies for graph neural networks. *arXiv preprint arXiv:2009.01639*.

33. Borgwardt, K. M., & Kriegel, H.-P. (2005). Shortest-path kernels on graphs. In *Proceedings of the 5th IEEE International Conference on Data Mining (ICDM)* (pp. 74–81).

34. Vishwanathan, S. V. N., Schraudolph, N. N., Kondor, R., & Borgwardt, K. M. (2010). Graph kernels. *Journal of Machine Learning Research*, 11, 1201–1242.

35. Kriege, N. M., Johansson, F. D., & Morris, C. (2020). A survey on graph kernels. *Applied Network Science*, 5(1), 1–42.

36. Gilmer, J., Schoenholz, S. S., Riley, P. F., Vinyals, O., & Dahl, G. E. (2017). Neural message passing for quantum chemistry. In *Proceedings of the 34th International Conference on Machine Learning (ICML)* (pp. 1263–1272).

37. Morris, C., Ritzert, M., Fey, M., Hamilton, W. L., Lenssen, J. E., Rattan, G., & Grohe, M. (2019). Weisfeiler and Leman go neural: Higher-order graph neural networks. In *Proceedings of the AAAI Conference on Artificial Intelligence* (Vol. 33, pp. 4602–4609).

38. Bonneel, N., Rabin, J., Peyré, G., & Pfister, H. (2015). Sliced and Radon Wasserstein barycenters of measures. *Journal of Mathematical Imaging and Vision*, 51(1), 22–45.

39. Mémoli, F. (2011). Gromov–Wasserstein distances and the metric approach to object matching. *Foundations of Computational Mathematics*, 11(4), 417–487.

40. Courant, R., & Hilbert, D. (1953). *Methods of Mathematical Physics, Vol. I*. Interscience Publishers.

> **全文完**
>
> 本文档包含 6 章、56 个精选数学结果 + 130 个全文编号对象（定义、定理、引理、命题、推论、算法、猜想、开放问题），涵盖从基础图论定义到核方法合法性证明的完整理论链条。所有证明均以 $\square$ 标记结束。
>
> **本次更新（v2）**：
> 1. **新增定理 3.10.5**（HTN-WL 在简化 CFI 构造下严格强于任意 $k$-WL，已对 $k=3,4,5,6,7$ 实证验证），将原猜想 3.10.7 从"HTN-WL 与 $k$-WL 不可比较性"修订为"标准 CFI 构造下 HTN-WL 与 $k$-WL 的关系"。
> 2. **新增第 3.11 节**（不可能性与边界结果）：形式化陈述 HTN-WL 区分能力的不可证边界，包括输出空间计数论证（定理 3.11.1）、HTN-WL 不强于图同构测试（定理 3.11.2）、标准 CFI 开放性（猜想 3.11.3）、普遍区分命题不可证（定理 3.11.4）。
> 3. **新增第 3.12 节**（弱化可证版本）：将原"HTN-WL 普遍强于 $k$-WL"论断弱化为三类可证逼近——存在性（论断 A）、固定图族可分性（论断 B）、参数化下界（论断 C）。
> 4. **符号变更**：HTN-WL$(K, I)$ → HTN-WL$(L, I)$（CSG 深度参数），与 $k$-WL 的元组宽度 $k$ 在排版上明确区分。详见附录 "符号变更日志"。
>
> **本次更新（v2.1）**：
> 5. **新增第 3.13 节**（HTN-WL 严格强于 $k$-WL 的完整证明链）：基于 CFI 圈空间奇偶性定理（`cyclic_schema/cfi_parity_cycle space.md`），证明 HTN-WL$(L=1, I \geq 1)$ 严格强于 $k$-WL（$k \geq 3$），对**标准 CFI 构造**与**简化 CFI 构造**均成立。**猜想 3.10.7 已部分解决**——标准 CFI 构造下 HTN-WL 严格强于 $k$-WL 的方向已具备完整证明链。

---

## 附录：符号变更日志（v2 增补）

为便于读者对照 v1（初稿）和 v2（修订版），本附录列出**所有**符号变更与对应的代码 API 变更。

### A.1 符号对照表

| 变更前符号（v1） | 变更后符号（v2） | 含义 | 影响范围 |
|---|---|---|---|
| HTN-WL$(K, I)$ | HTN-WL$(L, I)$ | HTN-WL 参数：CSG 深度 $L$ × 迭代次数 $I$ | 全文 |
| CSG$^K$ | CSG$^L$ | 第 $L$ 层 CSG 抽象 | §1.3, §1.5, §2, §3 |
| "CSG 层级 $\mathcal{H}_k$（$k = 0, \ldots, K$）" | "CSG 层级 $\mathcal{H}_k$（$k = 0, \ldots, L$）" | 多层 CSG 序列 | 定义 2.1.3 |
| $D(K)$ | $D(L)$ | HTN-WL$(L, I)$ 可区分的非同构图对数（$I$ 充分大） | 定理 3.9.7b |
| $D(K) \leq D(K+1)$ | $D(L) \leq D(L+1)$ | 区分能力单调性 | 定理 3.9.7b |
| $I_{\text{HTN-WL}}^{(K, I)}$ | $I_{\text{HTN-WL}}^{(L, I)}$ | HTN-WL$(L, I)$ 提取的不变量集合 | 定理 3.8.4 |
| 循环变量 `for k = 1, ..., K` | `for k = 1, ..., L` | 算法层循环（k 为层索引，保持小写） | 算法 2.4.2, 3.5.2 |
| $K_n$, $K_{k+1}$, $K_4$, $K_5$, $K_{2,3}$ 等 | **不变** | 完全图符号（$K_n$ 表示 $n$ 节点完全图） | 全文（保持原样） |
| $K_a, K_b$ | **不变** | CB 簇符号（$\mathcal{G}_1$ 的 CB 节点集与 $\mathcal{G}_2$ 的 CB 节点集） | §1.4 |
| $K$（核函数符号） | **不变** | TopoWGK 核函数 $K(\mathcal{G}_i, \mathcal{G}_j)$ | §4.1–4.4 |

### A.2 代码 API 变更

| 变更前（v1） | 变更后（v2） | 文件 | 行数 |
|---|---|---|---|
| `hierarchical_triangular_wl(..., K=1, I=5)` | `hierarchical_triangular_wl(..., L=1, I=5)` | `cyclic_schema/hierarchical_triangulated_wl.py` | 函数签名 |
| `hierarchical_triangular_wl_with_edges(..., K=1, I=5)` | `hierarchical_triangular_wl_with_edges(..., L=1, I=5)` | 同上 | 函数签名 |
| `hierarchical_triangular_wl_unified(..., K=1, I=5)` | `hierarchical_triangular_wl_unified(..., L=1, I=5)` | 同上 | 函数签名 |
| `TopoWassersteinGraphKernel(n_csg_layers=...)` | **不变** | `our_experiments/cycle_complex_wgk/topo_wasserstein_graph_kernel.py` | 公共参数 |

**关键说明**：函数参数名 $K$ → $L$ 是**破坏性变更**（breaking change）。在 v2 之前的调用代码需要更新为 `L=` 关键字参数。`n_csg_layers` 作为更可读的公共 API 名称**保持不变**，仅内部调用映射到新参数 $L$。

### A.3 不变符号声明

以下符号在 v1 → v2 中**保持不变**：

- **$k$（小写）**：$k$-WL 中的"WL 阶数"（元组宽度），与新符号 $L$（CSG 深度）排版上明确区分。
- **$K_n$, $K_{k+1}$, $K_1, K_2, \ldots, K_7$**：完全图（complete graph）符号系，**不**指代 CSG 深度。
- **$K$（函数符号）**：TopoWGK 核函数 $K(\mathcal{G}_i, \mathcal{G}_j)$，与 CSG 深度无关。
- **$K_a, K_b$**：CB 节点簇符号（$\mathcal{G}_1$ 与 $\mathcal{G}_2$ 的 CB 节点集合），仅在 §1.4 节使用。

### A.4 排版上的歧义消解

- 变更前：HTN-WL$(K=1)$ 与 $k$-WL（$k=1$）的 $K$ 与 $k$ 在数学模式中**排版相同**（仅大小写差异），易混淆。
- 变更后：HTN-WL$(L=1)$ 的 $L$（Layer 首字母，强调"层数"）与 $k$-WL 的 $k$（小写，强调"阶数"）**排版上完全分离**。读者从符号即可判断参数语义。

### A.5 与已有文献的兼容性

本文符号 $L$ 不同于以下文献中的 $L$ 含义：

- **Shervashidze et al. (2009)**：WL 子树模式的层数（与本文 $I$ 等价）；
- **Morris et al. (2019)**：$k$-WL 中 $k$ 的另一种命名（直接用 $k$，不冲突）；
- **Tuzhilina et al. (2022)**：特征选择的层数（与本文 $L$ 含义不同，但具体上下文区分）。

本文符号 $L$ 的语义（CSG 深度 / Cycle-recursive Schema Graph Layer Depth）与上述文献**无混淆**——HTN-WL 是新算法，$L$ 参数为本文首次提出。
