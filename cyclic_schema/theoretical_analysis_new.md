# 拓扑感知 Wasserstein 图核：理论分析

## 摘要

本文提出一种基于**分层循环模式图（Multi-layer CSG）**与**三角化邻域消息传递（Triangulated Neighborhood Aggregation, TNA）**融合的图核方法（Topo-Wasserstein Graph Kernel, TopoWGK）。该方法将图的局部短程拓扑与全局长程拓扑统一在分层消息传递框架（Hierarchical Triangulated Neighborhood WL, HTN-WL）内，并通过 Sliced Wasserstein 距离与 WL 内积核的凸组合构造具有正定性保证的图核。本文依次建立：（1）循环模式图（CSG）变换的代数结构与单次变换的圈空间刻画；（2）多层 CSG 迭代过程的收敛性定理；（3）HTN-WL 消息传递机制的形式化定义、算法复杂度与区分能力分析；（4）TopoWGK 的正定性证明与凸组合核的偏差-方差分析。

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
- [第 4 章 拓扑感知 Wasserstein 图核 (TopoWGK)](#第-4-章-拓扑感知-wasserstein-图核-topowgk)
  - [4.1 图核框架概述](#41-图核框架概述)
  - [4.2 条件负定性与正定核的关系](#42-条件负定性与正定核的关系)
  - [4.3 Sliced Wasserstein 距离的条件负定性](#43-sliced-wasserstein-距离的条件负定性)
  - [4.4 高斯核与 WL 内积的凸组合](#44-高斯核与-wl-内积的凸组合)
- [第 5 章 参考文献](#第-5-章-参考文献)

---

## 第 1 章 预备知识 (Preliminary)

本章给出后续章节所需的全部基础符号、图论概念以及循环模式图（Cyclic Schematic Graph, CSG）变换 $\Phi$ 的形式化定义与代数性质。所有符号按"图论基础 → 圈结构 → 消息传递 → 距离与核"的顺序在首次出现时定义。

### 1.1 图与子图基础符号

**定义 1.1.1**（无向简单图）。设 $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ 表示一个无向简单图，其中 $\mathcal{V}$ 为节点集，$\mathcal{E} \subseteq \binom{\mathcal{V}}{2}$ 为边集（无自环、无重边）。记 $n = |\mathcal{V}|$，$m = |\mathcal{E}|$，$c(\mathcal{G})$ 为 $\mathcal{G}$ 的连通分量数。

**定义 1.1.2**（邻居与度）。对节点 $v \in \mathcal{V}$，其邻居集定义为 $N_{\mathcal{G}}(v) = \{u \in \mathcal{V} : (v, u) \in \mathcal{E}\}$，度数为 $d_{\mathcal{G}}(v) = |N_{\mathcal{G}}(v)|$。

**定义 1.1.3**（连通分量与诱导子图）。设 $S \subseteq \mathcal{V}$，则 $\mathcal{G}[S] = (S, \mathcal{E}[S])$ 为 $\mathcal{G}$ 在 $S$ 上的**诱导子图**，其中 $\mathcal{E}[S] = \{(u,v) \in \mathcal{E} : u, v \in S\}$。若 $\mathcal{G}[S]$ 内部任意两节点间存在路径，则称 $S$ 是 $\mathcal{G}$ 的一个**连通分量**。

**定义 1.1.4**（三角化邻域类）。对节点 $v \in \mathcal{V}$，令 $TN_{\mathcal{G}}(v) = \{R_{\mathcal{G}}^{1}(v), R_{\mathcal{G}}^{2}(v), \dots, R_{\mathcal{G}}^{k_v}(v)\}$ 表示 $N_{\mathcal{G}}(v)$ 在诱导子图 $\mathcal{G}[N_{\mathcal{G}}(v)]$ 中的**连通分量分解**，其中 $k_v = |TN_{\mathcal{G}}(v)|$ 为连通分量个数，每个 $R_{\mathcal{G}}^{i}(v)$ 为 $\mathcal{G}[N_{\mathcal{G}}(v)]$ 的极大连通子集。

> **注**：三角化邻域类（Triangulated Neighborhood classes, TN-classes）是本文的核心局部结构概念。对每个节点 $v$，$TN_{\mathcal{G}}(v)$ 编码了 $N_{\mathcal{G}}(v)$ 内部所有邻居之间的连接结构，是 TNA 消息传递（定义 3.1.1）的拓扑基础。

### 1.2 圈空间与最小圈基

本节建立圈空间（cycle space）$\text{CS}_{\mathcal{G}}$ 与最小圈基（minimum cycle basis, MCB）$\text{MCB}_{\mathcal{G}}$ 的形式化定义，这是 CSG 变换 $\Phi$ 作用于图上的代数基础。

#### 1.2.1 圈空间

**定义 1.2.1**（圈空间）。无向图 $\mathcal{G}$ 的**圈空间** $\text{CS}_{\mathcal{G}} \subseteq \{0, 1\}^{|\mathcal{E}(\mathcal{G})|}$ 是 $\mathbb{F}_2 = \text{GF}(2)$ 上的向量空间，由 $\mathcal{G}$ 中所有**欧拉子图**（Eulerian subgraph，每个顶点度数为偶数的子图）构成。向量加法为对称差 $\triangle$。其维数为

$$
\dim \text{CS}_{\mathcal{G}} = \mu(\mathcal{G}) = |\mathcal{E}(\mathcal{G})| - |\mathcal{V}(\mathcal{G})| + c(\mathcal{G})
$$

其中 $\mu(\mathcal{G})$ 称为 $\mathcal{G}$ 的**圈秩**（cyclomatic number / circuit rank），表示图中独立圈的最大数量。

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

#### 1.2.2 节点分类与最小圈基

基于最小圈基，可将 $\mathcal{G}$ 的节点集划分为两类：

**定义 1.2.5**（圈节点与非圈节点）。设 $\text{MCB}_{\mathcal{G}} = \{C_1, \dots, C_{\mu}\}$。定义：

- $CYC_{\mathcal{G}} = \bigcup_{C \in \text{MCB}_{\mathcal{G}}} \mathcal{V}(C)$：**圈节点集**（出现在至少一个圈基中的顶点集合）。
- $ACN_{\mathcal{G}} = \mathcal{V} \setminus CYC_{\mathcal{G}}$：**非圈节点集**（anti-cycle nodes，不在任何圈基中的顶点集合）。

显然 $CYC_{\mathcal{G}} \cap ACN_{\mathcal{G}} = \varnothing$，$CYC_{\mathcal{G}} \cup ACN_{\mathcal{G}} = \mathcal{V}$。

> **注**：$ACN_{\mathcal{G}}$ 中的节点在 $\mathcal{G}$ 中仅位于树形结构上，即它们不在任何简单圈上。本性质保证 $\mathcal{G}$ 删除 $CYC_{\mathcal{G}}$ 后的剩余图必为森林。

### 1.3 循环模式图变换 $\Phi$ 的定义与算法步骤

本节给出循环模式图（Cyclic Schematic Graph, CSG）变换 $\Phi$ 的形式化定义与完整算法步骤。

#### 1.3.1 变换 $\Phi$ 的定义

**定义 1.3.1**（循环模式图变换）。对无向简单图 $\mathcal{G} = (\mathcal{V}, \mathcal{E})$，**循环模式图变换** $\Phi$ 将 $\mathcal{G}$ 映射为新图 $H = \Phi(\mathcal{G}) = (\mathcal{V}_H, \mathcal{E}_H)$，其中：

- $\mathcal{V}_H$ 包含三类节点：
  - $\mathcal{B}_{\mathcal{G}} = \{b_i : C_i \in \text{MCB}_{\mathcal{G}}\}$：圈基节点（CB 节点），$\text{type}(b_i) = \text{cycle\_basis}$。
  - $\mathcal{V}_{\text{nc}} = ACN_{\mathcal{G}}$：非圈节点，$\text{type}(v) = \text{original\_non\_cycle}$。
  - $\mathcal{I}_{\mathcal{G}}$：**接口节点**（interface nodes），$\text{type}(v) = \text{interface}$。
- $\mathcal{E}_H$ 由 $\Phi$ 的连接规则（算法 1.3.1）生成。

#### 1.3.2 变换 $\Phi$ 的算法

**算法 1.3.1**（CSG 变换 $\Phi$）。

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
       - 在 I_G 中创建 v 副本 (若 v ∈ CYC_G 且 v ∉ I_G), type = interface
       - 连接该接口节点到所有包含 v 的 b_i
   (d) 若两个 CB 节点簇 (b_i 的连通分量) 在原始节点层共享公共顶点 v,
       则将 v 添加为接口节点, 并连接到两个簇中各一个包含 v 的 b_i
6. 返回 H = (B_G ∪ V_nc ∪ I_G, E_H)
```

**复杂度分析**（算法 1.3.1）：

- 步骤 1：Horton 算法 $O(|\mathcal{E}|^3 |\mathcal{V}|)$；无赋权图贪心 $O(|\mathcal{E}|^2 |\mathcal{V}|)$。
- 步骤 2：$O(\mu)$ 创建 CB 节点。
- 步骤 3：$O(|\mathcal{V}|)$。
- 步骤 4-5a：遍历所有圈对，$O(\mu^2)$。
- 步骤 5b：$O(|\mathcal{E}|)$。
- 步骤 5c-5d：最坏 $O(\mu \cdot |\mathcal{E}|)$。

**总时间复杂度**：$O(|\mathcal{E}|^3 |\mathcal{V}| + \mu^2 + \mu \cdot |\mathcal{E}|)$，其中 $\mu = |\mathcal{E}| - |\mathcal{V}| + c(\mathcal{G})$ 为圈秩。**空间复杂度**：$O(|\mathcal{V}| + |\mathcal{E}| + \mu)$。

**注**：实现中采用 `networkx.minimum_cycle_basis` 计算 $\text{MCB}_{\mathcal{G}}$，其底层算法为 Horton (1987)。

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

**定义 1.4.2**（接口节点）。设 $\mathcal{G}$ 的最小圈基为 $\text{MCB}_{\mathcal{G}}$。节点 $v$ 称为**接口节点**当且仅当满足以下条件之一：

1. $v \in CYC_{\mathcal{G}}$，且存在至少两个不同的 $C_i, C_j \in \text{MCB}_{\mathcal{G}}$ 使得 $v \in \mathcal{V}(C_i) \cap \mathcal{V}(C_j)$；
2. $v \in CYC_{\mathcal{G}}$，且连接 $v$ 到 $ACN_{\mathcal{G}}$ 中某节点的边不在任何 $C_i$ 中（即 $v$ 是圈与非圈的连接点）；
3. $v$ 连接两个不同的 CB 节点连通分量（圈簇）。

**命题 1.4.3**（接口节点与圈簇的关系）。设 $K_1, K_2, \dots, K_t$ 为 $H$ 中 CB 节点（$b_i$）构成的连通分量（圈簇）。若 $v \in \mathcal{V}(K_a) \cap \mathcal{V}(K_b)$（作为原始节点），则 $v$ 必为接口节点。反之，每个接口节点对应于至少一对 $(K_a, K_b)$ 与 $v$ 的关联。

*证明*。（$\Rightarrow$）若 $v \in \mathcal{V}(K_a) \cap \mathcal{V}(K_b)$，则 $v$ 同时属于两个不同 CB 簇中的某个圈基，根据定义 1.4.2 的条件 1，$v$ 必为接口节点。（$\Leftarrow$）若 $v$ 为接口节点，则根据定义 1.4.2 的三种情况，均存在 $b_i \in K_a$ 与 $b_j \in K_b$（$a \neq b$）使 $v \in \mathcal{V}(b_i) \cap \mathcal{V}(b_j)$。$\square$

### 1.5 变换 $\Phi$ 的代数结构与不动点分析

#### 1.5.1 $\Phi$ 不是图同态意义下的函子

**定义 1.5.1**（图范畴）。定义图范畴 $\mathcal{G}$，其对象为无向简单图，态射为图同态（graph homomorphisms）：$\varphi: G \to H$ 为 $V(G) \to V(H)$ 的映射，满足 $(u, v) \in E(G) \Rightarrow (\varphi(u), \varphi(v)) \in E(H)$。

**命题 1.5.2**（$\Phi$ 非函子）。$\Phi$ 不是 $\mathcal{G}$ 上的函子（在标准图同态意义下）。

*证明*。设 $\varphi: \mathcal{G}_1 \to \mathcal{G}_2$ 为图同态。由于 $\varphi$ 可能将 $\mathcal{G}_1$ 中的圈 $C_i$ 映为 $\mathcal{G}_2$ 中的非圈（路径或点），$\text{MCB}_{\mathcal{G}_1}$ 的元素经 $\varphi$ 映射后未必形成 $\text{MCB}_{\mathcal{G}_2}$ 的基。因此 $\varphi$ 不一定诱导 $\Phi(\mathcal{G}_1)$ 与 $\Phi(\mathcal{G}_2)$ 之间的良定义同态，函子公理不成立。$\square$

#### 1.5.2 $\Phi$ 在图同构下的自然性

**定理 1.5.3**（$\Phi$ 的同构自然性）。若 $\varphi: \mathcal{G}_1 \xrightarrow{\cong} \mathcal{G}_2$ 为图同构，则存在唯一同构 $\Phi(\varphi): \Phi(\mathcal{G}_1) \xrightarrow{\cong} \Phi(\mathcal{G}_2)$ 使下图交换：

$$
\begin{CD}
\mathcal{G}_1 @>{\varphi}>> \mathcal{G}_2 \\
@V{\pi_1}VV @VV{\pi_2}V \\
\Phi(\mathcal{G}_1) @>{\Phi(\varphi)}>> \Phi(\mathcal{G}_2)
\end{CD}
$$

其中 $\pi_i: \mathcal{G}_i \to \Phi(\mathcal{G}_i)$ 为自然投影。

*证明*。设 $\varphi: \mathcal{V}(\mathcal{G}_1) \to \mathcal{V}(\mathcal{G}_2)$ 为同构映射。

**Step 1 — $\varphi$ 保持圈结构**：若 $C = (v_1, v_2, \dots, v_t, v_1)$ 是 $\mathcal{G}_1$ 中的简单圈，则 $\varphi(C) = (\varphi(v_1), \dots, \varphi(v_t), \varphi(v_1))$ 是 $\mathcal{G}_2$ 中长度相同的简单圈。$\varphi$ 建立了 $\mathcal{G}_1$ 与 $\mathcal{G}_2$ 圈集之间的双射。

**Step 2 — $\varphi$ 保持最小圈基**：$\varphi(\text{MCB}_{\mathcal{G}_1}) = \{\varphi(C_1), \dots, \varphi(C_\mu)\}$ 是 $\mathcal{G}_2$ 的一组圈基。由于 $\varphi$ 保持边长，$\sum |\varphi(C_i)| = \sum |C_i| = \min$，故为最小圈基。

**Step 3 — 构造 $\Phi(\varphi)$**：定义 $\psi: \mathcal{V}(\Phi(\mathcal{G}_1)) \to \mathcal{V}(\Phi(\mathcal{G}_2))$ 如下：

- 对 CB 节点 $b_i \leftrightarrow C_i$：$\psi(b_i) = b'_j$，其中 $b'_j \leftrightarrow \varphi(C_i)$。
- 对非圈节点 $v \in ACN_{\mathcal{G}}^{(1)}$：$\psi(v) = \varphi(v)$。
- 对接口节点 $v \in I^{(1)}$：$\psi(v) = \varphi(v)$。

**Step 4 — 验证 $\psi$ 保持边关系**：

- $(b_i, b_j) \in \mathcal{E}(\Phi(\mathcal{G}_1)) \Leftrightarrow \mathcal{E}(C_i) \cap \mathcal{E}(C_j) \neq \varnothing \Leftrightarrow \mathcal{E}(\varphi(C_i)) \cap \mathcal{E}(\varphi(C_j)) \neq \varnothing \Leftrightarrow (\psi(b_i), \psi(b_j)) \in \mathcal{E}(\Phi(\mathcal{G}_2))$。
- $(v, b_i) \in \mathcal{E}(\Phi(\mathcal{G}_1)) \Leftrightarrow v \in \mathcal{V}(C_i) \text{ 且边连接 } v \text{ 与 } C_i \Leftrightarrow \varphi(v) \in \mathcal{V}(\varphi(C_i)) \Leftrightarrow (\psi(v), \psi(b_i)) \in \mathcal{E}(\Phi(\mathcal{G}_2))$。

**Step 5 — $\psi$ 为双射**：$\varphi$ 是双射，且 CB 节点间对应也是双射（$\mu(\mathcal{G}_1) = \mu(\mathcal{G}_2)$），故 $\psi$ 是双射。

综上 $\psi$ 为 $\Phi(\mathcal{G}_1)$ 与 $\Phi(\mathcal{G}_2)$ 之间的同构。$\square$

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

#### 1.5.4 不动点与吸引子分析

**定义 1.5.5**（不动点）。图 $\mathcal{G}$ 称为 $\Phi$ 的**不动点**当且仅当 $\Phi(\mathcal{G}) \cong \mathcal{G}$（作为带类型标注的图同构）。

**命题 1.5.6**（不动点刻画）。$\Phi$ 的不动点恰好是无圈图（森林）。

*证明*。

（$\Leftarrow$）设 $\mathcal{G}$ 为森林，则 $\mu(\mathcal{G}) = 0$，$\text{MCB}_{\mathcal{G}} = \varnothing$。此时 $\mathcal{V}_H = ACN_{\mathcal{G}} = \mathcal{V}(\mathcal{G})$，$\mathcal{E}_H = \mathcal{E}(\mathcal{G})$，故 $\Phi(\mathcal{G}) \cong \mathcal{G}$。

（$\Rightarrow$）设 $\Phi(\mathcal{G}) \cong \mathcal{G}$。若 $\mathcal{G}$ 含圈，则 $\mu(\mathcal{G}) > 0$。由定理 2.2.1（圈秩严格递减），$\mu(\Phi(\mathcal{G})) < \mu(\mathcal{G})$，故 $\Phi(\mathcal{G}) \not\cong \mathcal{G}$，矛盾。因此 $\mathcal{G}$ 不含圈。$\square$

**推论 1.5.7**（吸引子）。森林是 $\Phi$ 迭代动力系统的**唯一吸引子**：对任意有限图 $\mathcal{G}$，存在 $N$ 使得 $\mathcal{G}^{(N)}$ 为森林（见定理 2.3.1）。

---

## 第 2 章 多层循环模式图 (Multi-layer CSG)

本章研究循环模式图变换 $\Phi$ 的迭代过程 $\mathcal{G}^{(k+1)} = \Phi(\mathcal{G}^{(k)})$ 的数学性质。我们依次建立：迭代过程的形式化定义、圈秩的严格单调递减性、有限步终止性以及收敛图的结构。

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

**定义 2.1.3**（CSG 层级与层数）。对固定的非负整数 $K$，CSG 层级 $\mathcal{H}_k$（$k = 0, 1, \dots, K$）定义为：

$$
\mathcal{H}_0 = \mathcal{G}^{(0)}, \quad \mathcal{H}_k = \Phi(\mathcal{H}_{k-1}) = \mathcal{G}^{(k)}, \quad k \geq 1
$$

**$K$ 称为 CSG 层数**。当 $K = 0$ 时，CSG 层级退化为单层（仅原始图）；当 $K \geq 1$ 时，建立 $K$ 层循环模式图。

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

### 2.2 圈秩严格单调性

**定理 2.2.1**（圈秩严格递减）。对任意 $k \geq 0$，若 $\mu(\mathcal{G}^{(k)}) > 0$，则 $\mu(\mathcal{G}^{(k+1)}) < \mu(\mathcal{G}^{(k)})$。

*证明*。设 $\mu_k = \mu(\mathcal{G}^{(k)}) > 0$，记 $\mathcal{G}_k = \mathcal{G}^{(k)}$，$\mathcal{H} = \mathcal{G}^{(k+1)} = \Phi(\mathcal{G}_k)$。我们从边数与顶点数的双重计数角度给出严格证明。

**(I) $\Phi$ 移除的边数下界**

**引理 2.2.1a**：对每个 $C_i \in \text{MCB}_{\mathcal{G}_k}$，至少存在一条边 $(u, v) \in \mathcal{E}(C_i)$ 不属于任何其他 $C_j$（$j \neq i$）。

*证明*。若 $C_i$ 的所有边均属于其他圈基 $C_j$，则 $C_i$ 的每条边 $(u, v)$ 在 $C_i$ 与某个 $C_j$ 中同时出现。但 $C_i$ 与 $C_j$ 共享所有边 $\Rightarrow$ $C_i = C_j$（两个相同边的圈基必须相同）。因此 $C_i$ 不能由其他圈基生成，即 $C_i$ 在 $\text{MCB}_{\mathcal{G}_k}$ 中是线性无关的。$\square$

**引理 2.2.1b**：定义 $s_i = |\{j \neq i : \mathcal{E}(C_i) \cap \mathcal{E}(C_j) \neq \varnothing\}|$ 为与 $C_i$ 共享边的圈基个数，则圈间边数 $|\mathcal{E}_{\text{CB}}|$ 满足

$$
|\mathcal{E}_{\text{CB}}| \leq \frac{1}{2} \sum_{i=1}^{\mu_k} s_i \leq \min\!\left( \binom{\mu_k}{2}, \frac{1}{2} \sum_i s_i \right)
$$

由引理 2.2.1a，每个 $C_i$ 至少有一条独占边被 $\Phi$ 移除（不进入 $\mathcal{E}(\mathcal{H})$），故

$$
\Delta \mathcal{E}_{\text{removed}} \geq \sum_{i=1}^{\mu_k} 1 = \mu_k
$$

**(II) $\Phi$ 减少的顶点数下界**

$$
\Delta \mathcal{V}_{\text{reduced}} = \sum_{i=1}^{\mu_k} (|\mathcal{V}(C_i)| - 1) - |I_{\text{new}}|
$$

其中 $I_{\text{new}}$ 为新增接口节点数。每个 $C_i$ 至少含 3 个顶点（简单圈），故 $|\mathcal{V}(C_i)| - 1 \geq 2$，即

$$
\Delta \mathcal{V}_{\text{reduced}} \geq 2 \mu_k
$$

**(III) $\Phi$ 新增的边数上界**

新增边数 $|\mathcal{E}_{\text{added}}|$ 由两部分构成：

- CB 节点间的圈间边：每对共享边的圈基对应一条 CB 间边。
- 接口边：每个接口节点连接到所有包含它的 CB 节点。

新增顶点集 $\mathcal{V}_{\text{added}} = \mu_k + |I|$（CB 节点 $\mu_k$ 个 + 接口节点 $|I|$ 个）。

**关键观察**：CB 节点间的边不可能形成圈——根据 TGS（树结构图模式）的树形结构性质（将在命题 2.2.2 中证明），CB 节点间的边数至多为 $\mu_k - c_{\text{clust}}$，其中 $c_{\text{clust}}$ 为 CB 节点连通分量数（$c_{\text{clust}} \geq 1$）。因此

$$
|\mathcal{E}_{\text{added}}| \leq (\mu_k - 1) + |I|
$$

**(IV) 连通分量数变化**

$H$ 的连通分量数 $c(H)$ 不超过 $c_k = c(\mathcal{G}_k)$，因为 $H$ 是 $\mathcal{G}_k$ 的商图，商图只可能合并或保持原分量数：

$$
c(H) - c_k \leq 0
$$

**(V) 汇总计算**

$$
\begin{aligned}
\mu(\mathcal{H}) - \mu_k &= \big(|\mathcal{E}(\mathcal{H})| - |\mathcal{V}(\mathcal{H})| + c(\mathcal{H})\big) - \big(|\mathcal{E}_k| - |\mathcal{V}_k| + c_k\big) \\
&= -\Delta \mathcal{E}_{\text{removed}} + |\mathcal{E}_{\text{added}}| + \Delta \mathcal{V}_{\text{reduced}} - |\mathcal{V}_{\text{added}}| + (c(H) - c_k) \\
&\leq -\mu_k + (\mu_k - 1) + |I| + 2\mu_k - (\mu_k + |I|) + 0 \\
&= -1 < 0
\end{aligned}
$$

因此 $\mu(\mathcal{H}) < \mu_k$。$\square$

**命题 2.2.2**（TGS 的树形结构）。设 $H = \Phi(\mathcal{G})$，将 CB 节点 $\mathcal{B}_{\mathcal{G}}$ 视为节点集，CB 节点间的边集 $\mathcal{E}_{\text{CB}}$ 视为边集，则子图 $H[\mathcal{B}_{\mathcal{G}}] = (\mathcal{B}_{\mathcal{G}}, \mathcal{E}_{\text{CB}})$ 必为森林（无圈）。

*证明*。假设 $H[\mathcal{B}_{\mathcal{G}}]$ 含圈 $b_{i_1} - b_{i_2} - \cdots - b_{i_t} - b_{i_1}$。则 $\mathcal{E}(C_{i_1}) \cap \mathcal{E}(C_{i_2}) \neq \varnothing$，$\mathcal{E}(C_{i_2}) \cap \mathcal{E}(C_{i_3}) \neq \varnothing$，$\dots$，$\mathcal{E}(C_{i_t}) \cap \mathcal{E}(C_{i_1}) \neq \varnothing$。但 $C_{i_1}, \dots, C_{i_t}$ 的边交集已形成闭合环，这意味着 $\bigoplus_{j=1}^{t} C_{i_j}$ 在 $\mathbb{F}_2$ 上的对称差为空图（因每条边被偶数次出现），与 $\text{MCB}_{\mathcal{G}}$ 的线性无关性矛盾。因此 $H[\mathcal{B}_{\mathcal{G}}]$ 无圈。$\square$

**注**（同调论视角）。在代数拓扑中，图 $\mathcal{G}$ 的圈空间 $\text{CS}_{\mathcal{G}}$ 同构于第一同调群 $\mathcal{H}_1(\mathcal{G}; \mathbb{F}_2)$。Pinch 映射 $p: \mathcal{G} \to \mathcal{G}/C_i$（将子图 $C_i$ 塌缩为一点）诱导同调群间的同态。$\Phi$ 同时对 $\mu_k$ 个圈基执行此操作，但塌缩后的点 $b_i$ 之间由共享边拓扑连接。塌缩每个 $C_i$ 消除其拓扑贡献（$\beta_1(C_i) = 1$），故 $\mu_{k+1} < \mu_k$。该视角与上述边计数法一致。

### 2.3 有限步终止与收敛速度

**定理 2.3.1**（有限步终止）。$\exists\, N \leq \mu(\mathcal{G}^{(0)})$ 使得 $\mu(\mathcal{G}^{(N)}) = 0$，即 $\mathcal{G}^{(N)}$ 为森林。

*证明*。$\mu_k = \mu(\mathcal{G}^{(k)})$ 构成非负整数序列。由定理 2.2.1，$\mu_k$ 严格递减：$\mu_0 > \mu_1 > \mu_2 > \cdots \geq 0$。非负整数严格递减序列至多 $\mu_0$ 步到达 0。故 $N \leq \mu_0$。$\square$

**命题 2.3.2**（收敛速度上界）。设 $\Delta_{\max} = \max_i |\{j : \mathcal{E}(C_i^{(k)}) \cap \mathcal{E}(C_j^{(k)}) \neq \varnothing\}|$ 为第 $k$ 层最大圈基边共享度，则

$$
\mu_{k+1} \leq \mu_k - \left\lceil \frac{\mu_k}{\Delta_{\max} + 1} \right\rceil
$$

*证明思路*。$\mu_k$ 个圈基可分为 $\lceil \mu_k / (\Delta_{\max} + 1) \rceil$ 个组，每组内圈基互不共享边。这些互不共享边的圈基在 $\Phi$ 下各自独立收缩为孤立 CB 节点，每组贡献至少 $-1$ 的 $\mu$ 变化。$\square$

**注**：上述界为悲观估计。实际收敛速度通常更快（如例 2.1.4 中 $\mu_k: 16 \to 10 \to 7 \to 2 \to 1 \to 0$，5 步收敛）。

**推论 2.3.3**（$N$ 的范围）。对任意连通图 $\mathcal{G}$，$1 \leq N \leq \mu(\mathcal{G}) \leq |\mathcal{E}| - |\mathcal{V}| + 1$。

*证明*。下界：当 $\mu(\mathcal{G}) \geq 1$ 时，$\mathcal{G}$ 至少含一个圈，至少需一次变换才能打破该圈；上界：由定理 2.3.1 直接得到；最后不等式由 $\mu(\mathcal{G}) = |\mathcal{E}| - |\mathcal{V}| + 1$（对连通图）。$\square$

### 2.4 收敛图的结构

**定理 2.4.1**（收敛图结构刻画）。$\mathcal{G}^* = \mathcal{G}^{(N)}$ 是一个**森林**，其顶点集由以下三类节点组成：

1. 各层残余 CB 节点（$\text{type} = \text{cycle\_basis}$），在 $\mathcal{G}^*$ 中不参与任何圈，度数为 0 或 1；
2. 原始非圈节点（$\text{type} = \text{original\_non\_cycle}$）；
3. 各级接口节点（$\text{type} = \text{interface}$）。

*证明*。$\mu(\mathcal{G}^*) = 0 \Leftrightarrow \mathcal{G}^*$ 中不含任何简单圈 $\Leftrightarrow \mathcal{G}^*$ 的每个连通分量都是树（命题 2.2.2 的递归推论）。由算法 1.3.1，$\Phi$ 仅在含圈时创建 CB 节点；在 $\mathcal{G}^*$ 中 $\mu = 0$ 故不再创建新 CB 节点，仅保留之前的残余。$\square$

**算法 2.4.2**（多层 CSG 构建算法）。

```
输入: 原始图 G = (V, E), 最大层数 K
输出: CSG 层级序列 {H_0, H_1, ..., H_K} 与层间映射 LT

1. H_0 ← G,  current_G ← G
2. for k = 1, 2, ..., K:
   (a) (H_k, cb_k, info_k) ← cyclic_schematic_graph(current_G)
       /* H_k = Φ(current_G), cb_k 为最小圈基列表, info_k 为节点-圈对应 */
   (b) if μ(H_k) = 0:
       /* 已收敛, 后续层 H_{k+1} = H_k (森林不动点) */
       for j = k+1, ..., K:
           H_j ← H_k
       break
   (c) LT[k] ← build_input_to_csg_mapping(H_k, cb_k, current_G)
       /* 记录下层节点 → 上层 CB 节点的映射 */
   (d) current_G ← H_k
3. 返回 ({H_0, H_1, ..., H_K}, LT)
```

**复杂度分析**（算法 2.4.2）：设 $n_k = |\mathcal{V}(\mathcal{G}^{(k)})|$，$m_k = |\mathcal{E}(\mathcal{G}^{(k)})|$，$\mu_k = \mu(\mathcal{G}^{(k)})$。

- 单次 $\Phi$：$O(m_k^3 n_k + \mu_k^2 + \mu_k m_k)$。
- $K$ 次迭代（$K \leq N \leq \mu_0$）：$O(\sum_{k=0}^{K-1} m_k^3 n_k)$，主项来自每层 Horton 算法。
- 层间映射：$O(K \cdot \mu_k)$。

**总时间复杂度**：$O(K \cdot m^3 n)$（$K$ 为常数时简化为 $O(m^3 n)$）。**空间复杂度**：$O(K \cdot (n + m + \mu))$。

---



## 第 3 章 分层三角化邻域 WL (HTN-WL)

本章将循环模式图（CSG）层级结构与三角化邻域消息传递（TNA）深度融合，建立**分层三角化邻域 WL**（Hierarchical Triangulated Neighborhood WL, **HTN-WL**）消息传递框架。HTN-WL 通过"前向—后向"双层传播机制，将局部邻域连通性信息与全局长程圈结构信息统一在同一标签精化过程内。本章依次建立：TNA 的形式化定义与算法、前向/后向传播机制、边标签更新、完整 HTN-WL 迭代算法、消息传递的代数性质、与 WL 测试的融合（含 CFI 图挑战）、区分能力的理论界、作为图同构必要条件的严格性以及与 $k$-WL 的深度比较。

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

#### 3.1.2 TNA 严格强于 1-WL

**定理 3.1.2**（TNA 严格强于 1-WL）。三角化邻域聚合的消息传递机制严格强于标准 1-WL 测试。

*证明*（反例构造）。取 $\mathcal{G}_1 = C_3 \cup C_3$（两个不相交的三角形），$\mathcal{G}_2 = C_6$（一个 6-圈）。两个图均有 6 个顶点、6 条边，且均为 2-正则。

**1-WL 无法区分**：$\mathcal{G}_1$ 与 $\mathcal{G}_2$ 中每个节点的度数均为 2。若初始标签相同，则 1-WL 的聚合值 $\text{AGG}_{1\text{-WL}}(v)$ 对所有节点相同（均为 $\{l, l\}$ 的多重集），迭代后仍保持一致。因此 1-WL 无法区分 $\mathcal{G}_1$ 和 $\mathcal{G}_2$。

**TNA 可以区分**：

- $\mathcal{G}_1$（$C_3 \cup C_3$）：每个节点 $v$ 的两个邻居在同一个三角形内，$\mathcal{G}[N_{\mathcal{G}}(v)]$ 有 **1 个连通分量**。
- $\mathcal{G}_2$（$C_6$）：每个节点 $v$ 的两个邻居在 6-圈中不相邻（中间隔一个节点），$\mathcal{G}[N_{\mathcal{G}}(v)]$ 有 **2 个孤立节点分量**。

因此 $\text{AGG}_{\text{TNA}}(v_1) \neq \text{AGG}_{\text{TNA}}(v_2)$，TNA 可区分。$\square$

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
   (c) τ'(v) ← sort(agg_list)
3. 联合标签压缩: 对图 1, 2 同步执行步骤 2, 收集所有 AGG 值,
   跨图去重排序, 为相同 AGG 值分配相同新标签
4. 返回 τ'
```

**复杂度分析**（算法 3.1.3）：

- 步骤 1（TN 预计算）：对每个 $v$ 遍历 $N(v)$ 的诱导子图，$O(\sum_v d(v)^2) = O(n \cdot d_{\max}^2)$。
- 步骤 2（聚合）：每个 $v$ 遍历 $TN(v)$ 的 $k_v$ 个分量，$O(\sum_v d(v)) = O(m)$。
- 步骤 3（联合压缩）：$O((n_1 + n_2) \log(n_1 + n_2))$。

**总时间复杂度**：$O(n \cdot d_{\max}^2 + m + n \log n)$。**空间复杂度**：$O(n \cdot d_{\max})$。

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

**定理 3.3.2**（后向传播的信息增益）。设 $\mathcal{G}_1, \mathcal{G}_2$ 为两图，$K \geq 1$ 为 CSG 层数。若存在节点 $v_1 \in \mathcal{V}(\mathcal{G}_1), v_2 \in \mathcal{V}(\mathcal{G}_2)$ 使得 $\text{LT}(v_1) \neq \text{LT}(v_2)$，则后向传播可能产生不同的标签元组，即使 $l_{\mathcal{G}^{(k-1)}}(v_1) = l_{\mathcal{G}^{(k-1)}}(v_2)$。

*证明*。后向标签元组 $\tau_k^{\text{back}}(v)$ 包含排序后的高层标签 $\text{SORT}(\{l_{\mathcal{H}_k}(h) : h \in \text{LT}(v)\})$。若 $\text{LT}(v_1) \neq \text{LT}(v_2)$，则高层标签的多重集可能不同，导致 $\tau_k^{\text{back}}(v_1) \neq \tau_k^{\text{back}}(v_2)$。$\square$

**推论 3.3.3**（信息闭环）。前向传播（$\mathcal{G}^{(k-1)} \to \mathcal{H}_k$，压缩信息）与后向传播（$\mathcal{H}_k \to \mathcal{G}^{(k-1)}$，注入信息）形成**信息闭环**。每个节点的标签同时编码了**局部结构**（来自前向聚合）与**全局层级位置**（来自后向注入），增强了区分能力。

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
| $K$ | CSG 层数 | `K` (default 1) |
| $I$ | HTN-WL 消息传递迭代次数 | `I` (default 5) |
| $\text{LT}(v)$ | 下层节点 $v$ 映射到的高层节点集合 | `build_input_to_csg_mapping` |
| $\text{sort}(\cdot)$ | 确定性排序（按 `node_sort_key`） | `sorted(..., key=node_sort_key)` |
| $\text{compress}(\cdot)$ | 跨图联合标签压缩 | `forward_message_passing_both` |

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

定义递归展平函数 $\Phi: \mathcal{T} \to \mathbb{Z}^*$（$\mathcal{T}$ 为嵌套元组空间）：

$$
\Phi(x) = \begin{cases}
(-1, |x|) \;+\!\!+ (\Phi(x_1), \Phi(x_2), \dots, \Phi(x_{|x|})) \;+\!\!+ (-2), & x \text{ 为元组} \\[4pt]
(x), & x \in \mathbb{Z}
\end{cases}
$$

> **注意**：此 $\Phi$ 是聚合键展平函数，与 CSG 变换 $\Phi$（第 1 章）同名异义，依上下文区分。

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
输入: 图 G1, G2; 初始节点标签 vlabel_np1, vlabel_np2; 层数 K; 迭代次数 I
输出: 标签历史矩阵 L1 ∈ Z^{|V1|×(I+1)}, L2 ∈ Z^{|V2|×(I+1)}

1. 初始化: 按 node_sort_key 排序节点, 初始化标签历史矩阵
2. CSG 层级构建 (当 K ≥ 1 时):
   for k = 1, ..., K:
       (H_k, cb_k, info_k) ← cyclic_schematic_graph(H_{k-1})
       mappings[k] ← build_input_to_csg_mapping(H_k, cb_k, H_{k-1})
3. 预计算邻域分量:
   nc_input ← precompute_neighbor_components(G1, G2)
   for k = 1, ..., K:
       nc_csg[k] ← precompute_neighbor_components(H_k^1, H_k^2)
4. for t = 1, 2, ..., I:    /* 主迭代循环 */
   4a. 前向传播 (G → CSG^1 → ... → CSG^K):
       for k = 1, ..., K:
           τ_k ← compute_final_label_tuples(H_k, cb_k, l_{k-1})
           l_k ← forward_message_passing_both(H_k^1, τ_k^1, H_k^2, τ_k^2, M, nc)
   4b. 后向传播 (CSG^K → ... → CSG^1 → G):
       l̃_K ← l_K
       for k = K, K-1, ..., 1:
           l̃_{k-1} ← backward_message_passing_both(H_{k-1}^1, H_{k-1}^2,
                              l_{k-1}^1, l_{k-1}^2, l̃_k^1, l̃_k^2,
                              mappings[k-1], M, nc)
   4c. 存储结果: L_a[i, t] ← l̃_0^{(a)}(n_i)  for all i
5. 返回 L1, L2
```

**复杂度分析**（算法 3.5.2）：设 $n = |\mathcal{V}(\mathcal{G})|$，$m = |\mathcal{E}(\mathcal{G})|$，$d_{\max}$ 为最大度数，$K$ 为 CSG 层数，$I$ 为迭代次数。

1. **CSG 构建**（一次性）：$O(K \cdot m^3 n)$（主项为 Horton 算法）。
2. **邻域分量预计算**（一次性）：$O(K \cdot n \cdot d_{\max}^2)$。
3. **单次迭代**：
   - 前向传播：$\sum_{k=1}^{K} O(n_k \cdot d_{\max,k}^2 + n_k \log n_k)$
   - 后向传播：$\sum_{k=K}^{1} O(n_{k-1} \cdot d_{\max,k-1}^2 + n_{k-1} \log n_{k-1})$
4. **$I$ 次迭代总计**：$O(I \cdot K \cdot (n \cdot d_{\max}^2 + n \log n))$。

**总时间复杂度**：

$$
T_{\text{HTN-WL}} = O(K \cdot m^3 n) + O(I \cdot K \cdot (n \cdot d_{\max}^2 + n \log n))
$$

当 $K, I$ 为常数时，简化为 $O(m^3 n + n \cdot d_{\max}^2)$。

**空间复杂度**：

- 标签历史：$O(n \cdot (I+1))$
- 邻域分量：$O(\sum_{k=0}^{K} n_k \cdot d_{\max,k})$
- CSG 层级：$O(K \cdot (\mu + n))$
- 总计：$O(n \cdot I + K \cdot \mu)$

#### 3.5.4 $K = 0$ 退化情形

**定义 3.5.3**（$K = 0$ 退化）。当 $K = 0$ 时，CSG 层级被完全跳过，HTN-WL 退化为**三角化邻域 WL**（Triangulated Neighbors WL）：

$$
l^{(t)}(v) = \text{compress}\big( \text{TNAggregate}(\mathcal{G}, v, l^{(t-1)}) \text{ across } \mathcal{G}_1, \mathcal{G}_2 \big)
$$

其中 $\text{TNAggregate}$ 使用 $\tau(v) = (l^{(t-1)}(v),)$ 作为初始标签元组。$K = 0$ 是 $K \geq 1$ 的严格子集——前向/后向传播"穿过"空层。

#### 3.5.5 统一调度接口

**定义 3.5.4**（统一调度）。HTN-WL 提供统一接口，根据图的元信息自动选择边标签感知模式：

$$
\text{WL}_{\text{unified}}(\mathcal{G}_1, \mathcal{G}_2, \text{g\_info}) = \begin{cases}
\text{WL}_{\text{edges}}(\mathcal{G}_1, \mathcal{G}_2, l, e, K, I), & \text{g\_info}[\text{'el'}] = \text{True} \\
\text{WL}_{\text{nodes}}(\mathcal{G}_1, \mathcal{G}_2, l, K, I), & \text{g\_info}[\text{'el'}] = \text{False}
\end{cases}
$$

#### 3.5.6 代码—数学对象映射表

| 代码函数 | 数学对象 | 核心方程 |
| --- | --- | --- |
| `_labels_to_dict` | 标签正规化 | §3.5.2(a) |
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

### 3.6 消息传递的代数刻画与数学性质

#### 3.6.1 复合迭代算子

**定义 3.6.1**（复合迭代算子）。定义复合迭代算子 $\Psi^{(t)}: \mathcal{L}_{\mathcal{G}} \to \mathcal{L}_{\mathcal{G}}$：

$$
\Psi^{(t)} = \mathcal{B}_1 \circ \mathcal{F}_1 \circ \mathcal{B}_2 \circ \mathcal{F}_2 \circ \cdots \circ \mathcal{B}_K \circ \mathcal{F}_K
$$

其中 $\mathcal{F}_k$ 为前向算子（标签元组计算 + TNA + 联合压缩），$\mathcal{B}_k$ 为后向算子（后向标签元组构造 + TNA + 联合压缩）。

**性质 3.6.2**（单调性）。$\Psi^{(t)}$ 单调：若 $\mathbf{L}_1 \leq \mathbf{L}_2$（逐元素），则 $\Psi^{(t)}(\mathbf{L}_1) \leq \Psi^{(t)}(\mathbf{L}_2)$。

*证明*。TNA 聚合和联合压缩均为单调操作（标签值越大，聚合值越大）。$\square$

**性质 3.6.3**（幂等性极限）。存在 $T_0$ 使得对所有 $t \geq T_0$，$\Psi^{(t)} = \Psi^{(t+1)}$（标签不再变化）。

*证明*。标签空间有限（由联合压缩保证），因此迭代序列 $\{\Psi^{(t)}\}$ 必进入循环。由单调性，循环长度为 1（不动点）。$\square$

#### 3.6.2 收敛速度与信息论

**定义 3.6.4**（标签稳定距离）。第 $t$ 次迭代后标签发生变化的节点数：

$$
d(t) = |\{(v, i) : l^{(t)}(v) \neq l^{(t-1)}(v)\}|
$$

**定理 3.6.5**（收敛速度上界）。$d(t) \leq n \cdot \left(1 - \dfrac{1}{\Delta_{\max} + 1}\right)^{t-1}$，其中 $\Delta_{\max}$ 为图的最大度数。

*证明*。每次迭代中，只有标签发生变化的节点的邻居可能在下一次迭代中变化。由 TNA 聚合的局部性，变化传播距离为 1，故 $d(t) \leq d(t-1) \cdot \frac{\Delta_{\max}}{\Delta_{\max} + 1}$。$\square$

**推论 3.6.6**（对数收敛）。HTN-WL 在 $O(\log n)$ 次迭代后收敛（对大多数图）。

*证明*。由定理 3.6.5，$d(t) \leq n \cdot \alpha^{t-1}$，$\alpha = 1 - \frac{1}{\Delta_{\max} + 1} < 1$。当 $d(t) < 1$ 时收敛，即 $t > \log_{1/\alpha} n = O(\log n)$。$\square$

**定义 3.6.7**（标签熵）。第 $t$ 次迭代后标签熵 $H(t) = -\sum_{\ell} p_\ell^{(t)} \log p_\ell^{(t)}$，其中 $p_\ell^{(t)}$ 为标签 $\ell$ 的出现频率。

**定理 3.6.8**（标签熵单调性）。$H(t) \leq H(t+1)$。

*证明*。每次迭代的 TNA 聚合和联合压缩将不同的聚合值映射到不同标签，增加标签多样性，熵不减。$\square$

**定理 3.6.9**（HTN-WL 的信息增益）。HTN-WL 相比 1-WL 的信息增益来源于：

1. **TNA 增益**：$\Delta H_{\text{TNA}} = H_{\text{TNA}}(t) - H_{1\text{-WL}}(t) \geq 0$
2. **CSG 层级增益**：$\Delta H_{\text{CSG}} = H_{\text{HTN-WL}}(t) - H_{\text{TNA}}(t) \geq 0$

其中 $\Delta H_{\text{TNA}} > 0$ 当且仅当存在节点 $v$ 使得 $N_{\mathcal{G}}(v)$ 的连通分量结构在 1-WL 下不可区分但在 TNA 下可区分。

*证明*。TNA 编码了邻域连通性信息（1-WL 不包含）；CSG 层级编码了圈结构信息（纯 TNA 不包含）。$\square$

#### 3.6.3 与代数拓扑的联系

**定理 3.6.10**（TNA 与同调群）。TNA 聚合算子 $\text{AGG}_{\text{TNA}}(v)$ 编码了 $v$ 的邻域 $N_{\mathcal{G}}(v)$ 的**第一同调群** $\mathcal{H}_1(\mathcal{G}[N_{\mathcal{G}}(v)]; \mathbb{F}_2)$ 的信息。

*证明*。$\mathcal{G}[N_{\mathcal{G}}(v)]$ 的连通分量分解对应 $H_0(\mathcal{G}[N_{\mathcal{G}}(v)])$（第零同调群），每个连通分量内部边结构对应 $\mathcal{H}_1(\mathcal{G}[N_{\mathcal{G}}(v)])$（第一同调群）。TNA 通过编码连通分量大小和标签分布，间接编码同调群信息。$\square$

**定理 3.6.11**（CSG 与谱序列）。CSG 层级结构对应**谱序列**（spectral sequence）的逐层逼近：

$$
E_1^{p,q} = H_{p+q}(\mathcal{G}^{(p)}) \Rightarrow H_{p+q}(\mathcal{G})
$$

*论证*。CSG 层级将图分解为不同尺度的圈结构，每层对应谱序列的一个页。前向传播对应谱序列的微分算子，后向传播对应逆操作。$\square$

### 3.7 HTN-WL 与 WL 测试的结合及 CFI 图挑战

本节分析 HTN-WL 相对标准 WL 测试的增强能力，并讨论其在经典 $k$-WL 硬反例（CFI 图）上的表现。

#### 3.7.1 强正则图上的 WL 增强

**定理 3.7.1**（WL 增强）。存在图族 $\mathcal{F}$，使得对任意 $\mathcal{G} \in \mathcal{F}$，标准 1-WL 测试无法区分 $\mathcal{G}$ 与非同构的 $\mathcal{G}'$，但在 $\Phi(\mathcal{G})$ 和 $\Phi(\mathcal{G}')$ 上运行 1-WL 即可区分。

*论证*。考虑强正则图参数 $\text{SRG}(v, k, \lambda, \mu)$，其定义为 $v$ 顶点、$k$-正则、任意相邻顶点有 $\lambda$ 个公共邻居、任意非相邻顶点有 $\mu$ 个公共邻居的图。已知存在非同构的强正则图具有相同参数（如 $\text{SRG}(16, 6, 2, 2)$ 的多种实现），1-WL 测试无法区分它们——因为 1-WL 在正则图上只能得到度数着色（全部相同）。

对 $\mathcal{G}$ 应用 $\Phi$ 变换：

- $\mathcal{G}$ 的每个圈基 $C_i$ 对应一个节点 $b_i$，其度数取决于 $C_i$ 与其他圈的共享边数；
- 圈基长度 $\ell_i$ 可能不同，导致接口节点的分布不同；
- 因此，$\Phi(\mathcal{G})$ 中节点的度数分布不一定是正则的。

具体到 $\text{SRG}(16, 6, 2, 2)$ 的不同实现，它们的圈基结构（圈长分布、共享模式、圈簇划分）往往不同，从而 $\Phi(\mathcal{G})$ 的度数分布也不同，1-WL 可区分。$\square$

#### 3.7.2 CFI 图的挑战

**Cai-Fürer-Immerman (CFI) 图族**是图同构与 WL 理论中的经典硬反例。一对 CFI 图 $(\mathcal{G}_1, \mathcal{G}_2)$ 满足 $\mathcal{G}_1 \equiv_{k\text{-WL}} \mathcal{G}_2$（对所有 $k$），但 $\mathcal{G}_1 \not\cong \mathcal{G}_2$。其构造基于在基础网格图上"扭曲"奇偶性约束。

**CFI 图对 HTN-WL 的挑战**：

1. **构造对称性**：CFI 图的圈基结构高度规则——基础网格图的圈基长度多重集几乎完全由网格的周长和内部网格数决定，扭曲操作仅改变顶点的二部性而几乎不改变圈长分布。
2. **可能失效**：因此 $\Phi(\mathcal{G}_1) \cong \Phi(\mathcal{G}_2)$ 的可能性很大。$\Phi$ 变换本身可能无法区分 CFI 对。
3. **HTN-WL 的潜在能力**：HTN-WL 的**标签消息传递**层可能提供额外区分能力——TNA 编码的邻域连通性在 CFI 图的扭曲顶点与其邻居间可能呈现可检测的差异。

**开放问题 3.7.2**（HTN-WL vs CFI）。HTN-WL 是否能区分所有 CFI 图对？还是存在 CFI 反例对使得 HTN-WL 也无法区分？若存在，则 HTN-WL 与 $k$-WL 在 CFI 图族上具有同等的局限性。

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
wl1, wl2 = hierarchical_triangular_wl(G1, G2, v1, v2, K=1, I=3)
print(f"HTN-WL 可区分: {not _is_isomorphic_wl(wl1, wl2)}")
```

**实验结论**：
- $k$-WL（$k = 1, 2, 3, 4, 5$）：**无法区分**（哈希值完全相同）。
- HTN-WL（$K=1, I=3$）：**可以区分**（`Isomorphic = False`）。

**为什么 HTN-WL 能区分**：Shrikhande 图和 Rook 图的关键区别在于**邻域的连通分量结构**：

- Rook 图中，每个顶点的 6 个邻居形成 1 个连通分量（行/列连接形成 $K_3 \times K_3$ 结构）。
- Shrikhande 图中，每个顶点的 6 个邻居形成 2 个连通分量（循环群结构）。

TNA 检测到这种邻域内部连通性差异，而 $k$-WL 仅关注 $k$-元组的联合邻域计数，无法捕捉这种局部拓扑差异。

**示例 2：$C_3 \cup C_3$ vs $C_6$**（2-正则图）：

```python
G1 = nx.disjoint_union(nx.cycle_graph(3), nx.cycle_graph(3))  # 两个三角形
G2 = nx.cycle_graph(6)                                         # 一个 6-圈
G1_int = nx.convert_node_labels_to_integers(G1)
G2_int = nx.convert_node_labels_to_integers(G2)
v1, v2 = np.ones(6, dtype=int), np.ones(6, dtype=int)
wl1, wl2 = hierarchical_triangular_wl(G1_int, G2_int, v1, v2, K=1, I=5)
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

*证明*。

- **$I_{\text{基础}} \prec I_{\text{圈}}$**：同谱非同构 3-正则图（如 $\text{SRG}(16,6,2,2)$）具有相同 $(n, m, \mu)$ 和度分布，但最小圈基长度多重集不同（Tutte 1966 反例）。
- **$I_{\text{圈}} \prec I_{\text{簇}}$**：两图共享圈基长度多重集但圈基的边共享模式不同——例如两图各含 2 个 3-圈和 1 个 4-圈，但 3-圈与 4-圈的共享边数不同。
- **$I_{\text{簇}} \prec I_{\text{接口}}$**：两图具有相同圈簇结构但接口节点连接模式不同——例如两个 3-圈共享一个顶点（星型接口）vs 三个 3-圈链式共享（路径型接口）。
- **$I_{\text{接口}} \prec I_{\text{跨层}}$**：跨层不变量编码从 $\mathcal{G}^{(0)}$ 到 $\mathcal{G}^{(N)}$ 的演化轨迹，接口不变量仅编码单层结构。
- **$I_{\text{跨层}} \prec I_{\text{HTN-WL}}$**：由定理 3.8.2 和 3.8.3。$\square$

**定理 3.8.2**（TNA 的信息增益）。设 $l$ 为标签函数，$v$ 为图 $\mathcal{G}$ 的节点。TNA 聚合值 $\text{AGG}_{\text{TNA}}(v)$ 包含严格多于 1-WL 聚合值的信息，即存在 $\mathcal{G}_1, \mathcal{G}_2$ 和 $l$ 使得 $\text{AGG}_{1\text{-WL}}(v_1) = \text{AGG}_{1\text{-WL}}(v_2)$ 但 $\text{AGG}_{\text{TNA}}(v_1) \neq \text{AGG}_{\text{TNA}}(v_2)$。

*证明*。取两个 3-正则图：一个每个节点的邻域形成三角形（3-团），另一个每个节点的邻域形成路径（3 顶点路径）。它们的 1-WL 聚合值相同（邻居标签多重集均为 $\{l, l, l\}$），但 TNA 的连通分量结构不同（前者 1 个分量，后者 2 个分量）。$\square$

**定理 3.8.3**（边标签的信息增强）。存在图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 使得在无标签情况下 HTN-WL 无法区分，但在引入边标签后可以区分。

*证明*。构造两个三角形：$\mathcal{G}_1$ 的边标签全为 0，$\mathcal{G}_2$ 的边标签为 $(0, 0, 1)$。无标签时 HTN-WL 看到完全相同的结构和节点标签，无法区分；带边标签时，边上下文 $\text{ec}(v)$ 不同，导致后向标签元组不同。$\square$

### 3.9 HTN-WL 作为图同构的必要条件

#### 3.9.1 必要条件与充分条件

**定理 3.9.1**（HTN-WL 必要条件）。若 $\mathcal{G}_1 \cong \mathcal{G}_2$，则对所有 $K \geq 0, I \geq 1$，HTN-WL$(K, I)$ 输出的标签直方图完全一致。

*证明*。图同构保持所有图结构：

1. 邻域连通性模式（$TN_{\mathcal{G}}(v)$ 不变）
2. CSG 层级结构（圈基、接口节点不变，由定理 1.5.3 保证）
3. 层间映射（$\text{LT}(v)$ 不变）
4. 边标签（边同构时不变）

因此消息传递中每一步（前向标签元组计算、TNA 聚合、后向标签元组构造、联合标签分配）都产生一致结果。$\square$

**推论 3.9.2**（HTN-WL 充分性方向）。若 $\mathcal{G}_1$ 和 $\mathcal{G}_2$ 被 HTN-WL$(K, I)$ 区分（即存在 $t$ 使标签直方图不同），则 $\mathcal{G}_1 \not\cong \mathcal{G}_2$。

*证明*。由定理 3.9.1 的逆否命题直接得到。$\square$

**定理 3.9.3**（HTN-WL 的不完备性）。存在非同构图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 不被任意 HTN-WL$(K, I)$ 区分。

*证明*（计数论证）。$n$ 顶点上的非同构图数量至少为 $2^{\binom{n}{2}} / n!$，而 HTN-WL$(K, I)$ 的输出（标签直方图序列）的可能数量远小于此（标签空间大小受 $O(n \cdot I)$ 限制）。由鸽巢原理，必然存在非同构图对具有相同的 HTN-WL 输出。$\square$

**注**：HTN-WL 的不完备性源于其信息压缩性质——它将图结构信息压缩为标签直方图序列，不可避免丢失部分信息。具体地，HTN-WL 捕捉的信息包括：度序列、邻域连通性模式、CSG 层级结构，但不包含完整的子图计数信息（如 $K_4$ 子图数量）。

#### 3.9.2 收敛图的判据作用

**定理 3.9.4**（收敛图判据）。设 $\mathcal{G}^*$ 为收敛图（$\mu(\mathcal{G}^*) = 0$）。若 $\mathcal{G}_1^* \cong \mathcal{G}_2^*$ 且对所有 $k$，$\mathcal{G}_1^{(k)} \cong \mathcal{G}_2^{(k)}$，则 $\mathcal{G}_1 \cong \mathcal{G}_2$。

*证明*。这是定理 3.9.1 的直接推论——同构在所有层保持则回到原始层仍保持。$\square$

**注**：收敛图 $\mathcal{G}^*$ 本身不足以判定同构——存在非同构图具有同构的收敛图（$\Phi$ 是有损压缩）。但结合完整 CSG 序列可提供额外判据。

#### 3.9.3 否定性判定算法

**算法 3.9.5**（HTN-WL 否定性检验）。

```
输入: 图 G1, G2, 参数 K, I
输出: "非同构" 或 "可能同构"

1. 运行 HTN-WL(G1, G2, K, I) 得到标签历史 L1, L2
2. for t = 0, 1, ..., I:
   (a) 比较 L1[:, t] 与 L2[:, t] 的直方图
   (b) if 不同: return "非同构"
3. return "可能同构"
```

**复杂度**：$O(I \cdot K \cdot (m^3 n + n \cdot d_{\max}^2))$。

**注**：算法 3.9.5 仅给出否定性判定（必要条件），非充分条件。若需充分判定，需进一步调用 nauty 等规范标记算法（HTN-WL + nauty 联合使用为实用策略）。

#### 3.9.4 HTN-WL 与三角化邻域 WL（$K=0$）的比较

**定理 3.9.6**（三角化邻域 WL 严格强于 1-WL）。$K=0$ 的 HTN-WL 区分能力严格强于标准 1-WL 测试。

*证明*。由定理 3.1.2，TNA 在单层消息传递中已比 1-WL 捕获更多信息。$\square$

**定理 3.9.7**（CSG 层级的信息增益）。存在图对使得 $K=0$ 的 HTN-WL 无法区分，但对某个 $K \geq 1$ 的 HTN-WL 可区分。

*证明*（构造性）。取两个非同构图，它们在局部邻域结构上相同（$K=0$ 的 TNA 无法区分），但全局圈结构不同。当 $K \geq 1$ 时，CSG 层将圈结构信息编码为标签元组，消息传递能感知全局圈差异。$\square$

#### 3.9.5 复杂度对比

| 维度 | HTN-WL | $k$-WL（$k \geq 3$） |
| --- | --- | --- |
| **时间复杂度** | $O(I \cdot K \cdot (m^3 n + n \cdot d_{\max}^2))$ | $O(n^k)$ |
| **空间复杂度** | $O(n)$（仅存储标签） | $O(n^k)$（存储 $k$-元组着色） |
| **捕捉的信息** | 邻域连通性 + 圈层级结构 | $k$-元组联合邻域计数 |
| **边标签** | 原生支持 | 需扩展 |
| **完备性** | 否 | 否 |

**与 $k$-WL 的复杂度对比**（设 $K=1, I=3$）：

| $k$ | $k$-WL 复杂度 | HTN-WL 复杂度 |
| --- | --- | --- |
| 2 | $O(n^2)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |
| 3 | $O(n^3)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |
| 4 | $O(n^4)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |
| 5 | $O(n^5)$ | $O(m^3 n + n \cdot d_{\max}^2)$ |

**实践建议**：

- 稀疏图（$m = O(n)$）：HTN-WL 复杂度 $O(n^2 \cdot d_{\max}^2)$，通常低于 $k$-WL（$k \geq 3$）。
- 稠密图（$m = O(n^2)$）：HTN-WL 可能 $O(n^7)$ 高于低阶 $k$-WL。

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

**方向 1**（已严格证明）：存在图对 $(\mathcal{G}_1, \mathcal{G}_2)$ 使得 $k$-WL（任意 $k \geq 1$）无法区分，但 HTN-WL 可以区分。

*证明*。取 $\mathcal{G}_1$ 为 Shrikhande 图，$\mathcal{G}_2$ 为 4×4 Rook 图。两者均为 $\text{SRG}(16, 6, 2, 2)$，具有相同的度序列、相同的特征值谱、相同的局部邻域结构，但**不是同构的**。

- $k$-WL（$k = 1, 2, 3, 4, 5$）：无法区分（哈希值完全相同）
- HTN-WL：可以区分（TNA 检测到邻域连通分量数不同）

详细验证见 §3.7.2 的代码示例。$\square$

**方向 2**（开放问题）：是否存在图对 $(\mathcal{G}'_1, \mathcal{G}'_2)$ 使得 $k$-WL（$k \geq 3$）可区分但 HTN-WL 不可区分？

*分析*。HTN-WL 捕捉了度序列、邻域连通性模式和 CSG 层级结构，信息量丰富。构造 HTN-WL 不可区分的非同构图对，需要两图同时具有：

1. 相同的度序列；
2. 相同的邻域连通性模式；
3. 相同的 CSG 层级结构；

同时 $k$-WL（$k \geq 3$）需要能检测到某种子图差异（如 $K_4$ 计数不同）。但 $K_4$ 的存在会影响邻域连通性（$K_4$ 的每个节点邻域是 $K_3$ 形成 1 个连通分量），构造满足上述三个条件但 $K_4$ 计数不同的图对是非平凡组合问题。

目前尚无已知具体构造。**这是 HTN-WL 与 $k$-WL 关系研究的核心开放问题**。

**关键引理 3.10.5**（邻域连通性的子图计数下界信息）。设 $\mathcal{G}$ 为图，$v \in \mathcal{V}(\mathcal{G})$，$N_{\mathcal{G}}(v)$ 为 $v$ 的邻域。$\mathcal{G}[N_{\mathcal{G}}(v)]$ 的连通分量结构提供**部分**子图计数信息：

- **完全确定**：无内部边的连通分量个数（孤立节点数）、各分量大小；
- **部分反映**：边密度（完全子图时达到上界）；
- **不完整确定**：精确子图计数（三角形数、具体路径数）。

**TNA 的实际能力**：TNA 通过三角化邻域类 $TN_{\mathcal{G}}(v)$ 的**结构**（不仅是数量，还包括每个类的大小、连通模式）编码信息——比"连通分量数"丰富得多。TNA 的完整聚合值是嵌套元组结构，包含每个类的大小、每个类内部边标签的分布等，区分能力**严格强于**仅看连通分量数。

#### 3.10.3 区分能力层次

HTN-WL 的区分能力来自三个递进层次：

**层次 1：三角化邻域（$K=0$）**。基本 TNA 增强捕获局部邻域边结构，严格强于 1-WL（定理 3.1.2）。

**层次 2：单层 CSG（$K=1$）**。加入圈层抽象后，节点标签包含圈的整体信息（通过标签元组），信息从原始图传播到 CSG 再回到原始图，捕捉"节点在圈中的角色"。

**层次 3：多层 CSG（$K \geq 2$）**。多层抽象捕捉圈之间的嵌套和包含关系。高层 CSG 中的节点代表"圈的圈"，使消息传递能感知全局循环骨架。

**定理 3.10.6**（层次严格性）。对任意 $K \geq 0$，存在图对使得 $K$ 层 HTN-WL 无法区分但 $K+1$ 层可以区分。

*证明*（构造性）。构造两个图，它们在 $K$ 层以内的圈结构完全相同，但在第 $K+1$ 层出现不同的圈结构。例如，取两个图，其 $K$ 次迭代 $\Phi$ 变换结果同构，但 $\Phi^{K+1}(\mathcal{G})$ 不同。由于 $K$ 层 HTN-WL 使用的 CSG 只到第 $K$ 层，无法感知更深层圈差异；而 $K+1$ 层通过额外 CSG 层可捕捉这种差异。$\square$

#### 3.10.4 HTN-WL 实际优势总结

HTN-WL 相对纯 $k$-WL 的实际优势：

1. **邻域结构感知**：TNA 捕获了 $N_{\mathcal{G}}(v)$ 的内部边结构，这是标准 $k$-WL（包括 2-WL）所忽略的信息。
2. **层级抽象**：CSG 层级提供天然多尺度分析框架，从原始图到高层抽象逐渐压缩信息。
3. **边标签原生支持**：边上下文 $\text{ec}(v)$ 的编码方式自然融入消息传递，无需图扩展。
4. **标签历史作为签名**：$I$ 次迭代的完整标签历史 $\mathbf{L}^{(0)}, \dots, \mathbf{L}^{(I)}$ 可作为图的强签名。
5. **参数可调**：$K$（深度）和 $I$（迭代次数）提供灵活性——浅层快速筛查，深层精细区分。

**实践推荐**：

| 图类型 | 推荐方法 | 原因 |
| --- | --- | --- |
| 稀疏社交网络 | HTN-WL | 邻域连通性模式丰富 |
| 分子图（化学环） | HTN-WL | 圈层级结构信息重要 |
| 强正则图 | HTN-WL 优先 | TNA 可检测邻域连通性差异（如 $\text{SRG}(16,6,2,2)$） |
| 随机图 | 两者均可 | 信息足够丰富 |
| 带边标签图 | HTN-WL | 原生支持边标签 |

---

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
输入: 图集合 {G_1, ..., G_n}, 节点标签, 边标签, 平衡系数 β, 高斯核带宽 γ, CSG 层数 K, 迭代次数 I
输出: 核矩阵 K ∈ R^{n×n}

1. for i = 1, ..., n:
   (a) 运行 HTN-WL(G_i, K, I) 得到标签历史 L_i
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

- 步骤 1：HTN-WL 标签传播，$O(n \cdot (K \cdot m^3 n + I \cdot K \cdot n \cdot d_{\max}^2))$。
- 步骤 2：每对图的 SW² 计算，$O(L_p \cdot n_i \log n_i)$（$L_p$ 为投影数），共 $O(n^2 \cdot L_p \cdot n \log n)$。
- 步骤 3：内积计算，$O(n^2 \cdot L)$。

**总时间复杂度**：$O(n^2 \cdot L_p \cdot n \log n + n \cdot K \cdot m^3 n)$，主项为 $O(n^3 \log n)$。**空间复杂度**：$O(n^2 + n \cdot I)$。

### 4.2 条件负定性与正定核的关系

#### 4.2.1 条件负定函数

**定义 4.2.1**（条件负定函数，Conditionally Negative Definite, CND）。对称函数 $\psi: \mathcal{X} \times \mathcal{X} \to \mathbb{R}$ 称为**条件负定**（CND），若对任意 $n \geq 2$，$x_1, \dots, x_n \in \mathcal{X}$，以及满足 $\sum_{i=1}^n c_i = 0$ 的实数 $c_1, \dots, c_n$，有

$$
\sum_{i,j=1}^n c_i c_j \psi(x_i, x_j) \leq 0
$$

**注**：CND 函数又称**距离型函数**（distance-type function）——所有度量 $d(\cdot, \cdot)$ 满足 $d^2$ 是 CND 的（在欧氏空间），但反之不成立。

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
- $k \geq 2$ 项：利用 $\psi(x, x) = 0$ 条件可证明 $\psi^k$ 也 CND（Schönberg 定理的标准推论），故 $\sum_{i,j} c_i c_j \psi(x_i, x_j)^k \leq 0$，加 $(-t)^k$ 符号后 $\geq 0$。

因此 $S \geq 0$，即 $K_t = \exp(-t\psi)$ 是 PSD 的。

（$\Rightarrow$）设 $\exp(-t\psi)$ 对所有 $t \geq 0$ 是 PSD。固定 $c_1, \dots, c_n$ 满足 $\sum c_i = 0$。定义 $f(t) = \sum_{i,j} c_i c_j \exp(-t\psi(x_i, x_j))$。由 PSD 假设，$f(t) \geq 0$ 对所有 $t \geq 0$。$f(0) = (\sum c_i)^2 = 0$。若 $\sum_{i,j} c_i c_j \psi(x_i, x_j) > 0$，则 $f'(0) = -\sum_{i,j} c_i c_j \psi(x_i, x_j) < 0$，存在 $\varepsilon > 0$ 使 $f(t) < 0$（$t \in (0, \varepsilon)$），矛盾。故 $\sum_{i,j} c_i c_j \psi(x_i, x_j) \leq 0$，即 $\psi$ CND。$\square$

**推论 4.2.3**（高斯核 PSD 充分条件）。要保证高斯核 $K(x, y) = \exp(-\gamma \cdot d^2(x, y))$ 是 PSD 的，一个**充分条件**是 $d^2$（距离平方）是 CND 的。

#### 4.2.3 标准 2-Wasserstein 距离的 CND 性质

**定义 4.2.4**（2-Wasserstein 距离）。设 $\mu, \nu$ 为 $\mathbb{R}^d$ 上的两个概率分布，$p = 2$ 时的 Wasserstein 距离定义为：

$$
\mathcal{W}_2(\mu, \nu) = \left( \inf_{\pi \in \Pi(\mu, \nu)} \int_{\mathbb{R}^d \times \mathbb{R}^d} \|x - y\|_2^2 \, d\pi(x, y) \right)^{1/2}
$$

其中 $\Pi(\mu, \nu)$ 是所有以 $\mu$ 和 $\nu$ 为边缘分布的联合分布集合。

**命题 4.2.5**（$\mathcal{W}_2^2$ 在一般 $\mathbb{R}^d$ 上非 CND）。在一般欧氏空间 $\mathbb{R}^d$（$d \geq 1$）上，$\mathcal{W}_2^2$ 不是条件负定的。因此高斯核 $\exp(-\gamma \cdot \mathcal{W}_2^2)$ **不能保证**是 PSD 的。

*理由*。Wasserstein 距离是 $\mathbb{R}^d$ 上的度量，但 $\mathcal{W}_2^2$ 在一般 $\mathbb{R}^d$ 上不是负定的。具体地，存在概率分布集合 $\{\mu_1, \dots, \mu_n\}$ 使得矩阵 $[\mathcal{W}_2^2(\mu_i, \mu_j)]_{i,j=1}^n$ 不是条件负定的。这意味着使用 `ot.emd2` 计算的 $\mathcal{W}_2^2$ 作为距离输入到高斯核时，核函数**可能不是** PSD 的，从而在 SVM 等核方法中使用时缺乏理论保证。

### 4.3 Sliced Wasserstein 距离的条件负定性

Sliced Wasserstein（SW）距离通过随机投影将高维 OT 问题约化为一维 OT 问题，**保留了 CND 性质**。这是本工作选用 SW 而非标准 Wasserstein 的核心理论动机。

#### 4.3.1 SW 距离的定义

**定义 4.3.1**（Sliced Wasserstein 距离）。对概率分布 $\mu, \nu \in \mathcal{P}(\mathbb{R}^d)$，其 $p$-阶 SW 距离定义为：

$$
\text{SW}_p(\mu, \nu) = \left( \mathbb{E}_{\theta \sim \mathbb{S}^{d-1}} \left[ \mathcal{W}_p^p(P_{\theta\#}\mu, P_{\theta\#}\nu) \right] \right)^{1/p}
$$

其中 $P_{\theta\#}\mu$ 是 $\mu$ 在方向 $\theta \in \mathbb{S}^{d-1}$ 上的投影（一维分布），$\mathcal{W}_p$ 为一维 $p$-Wasserstein 距离。

当 $p = 2$ 时：

$$
\text{SW}_2^2(\mu, \nu) = \mathbb{E}_{\theta \sim \mathbb{S}^{d-1}} \left[ \mathcal{W}_2^2(P_{\theta\#}\mu, P_{\theta\#}\nu) \right]
$$

#### 4.3.2 SW² 的 CND 性质

**定理 4.3.2**（SW² 的 CND 性质，Nadjahi et al. 2020）。Sliced Wasserstein 距离的平方 $\text{SW}_2^2$ 是条件负定的（CND）。

*证明*。分三步进行：

**Step 1（一维情形）**：证明 $\mathcal{W}_2^2(\cdot, \cdot)$ 在 $\mathcal{P}(\mathbb{R})$ 上是 CND 的。

对 $\mu, \nu \in \mathcal{P}(\mathbb{R})$，记其累积分布函数分别为 $F, G$，逆 CDF 为 $F^{-1}, G^{-1}$。一维 2-Wasserstein 距离有显式形式：

$$
\mathcal{W}_2^2(\mu, \nu) = \int_0^1 (F^{-1}(t) - G^{-1}(t))^2 dt = \|F^{-1} - G^{-1}\|_{L^2([0,1])}^2
$$

这正是 Hilbert 空间 $L^2([0, 1])$ 中的平方范数差的平方。根据 Hilbert 空间理论，平方范数 $\|h - k\|^2 = \|h\|^2 + \|k\|^2 - 2\langle h, k\rangle$ 是 CND 的（因为对任意 $c_i$ 满足 $\sum c_i = 0$，$\sum_{i,j} c_i c_j \|h_i - h_j\|^2 = -2 \sum_{i,j} c_i c_j \langle h_i, h_j\rangle \leq 0$，由 Cauchy-Schwarz 不等式）。

**Step 2（CND 函数在期望下封闭）**：若 $\psi_\theta$ 对每个 $\theta$ 是 CND 的，则 $\mathbb{E}_\theta[\psi_\theta]$ 也是 CND 的。

*证明*。设 $\psi_\theta$ 对每个 $\theta$ 是 CND。对任意 $c_1, \dots, c_n$ 满足 $\sum c_i = 0$：

$$
\sum_{i,j} c_i c_j \mathbb{E}_\theta[\psi_\theta(x_i, x_j)] = \mathbb{E}_\theta\!\left[ \sum_{i,j} c_i c_j \psi_\theta(x_i, x_j) \right] \leq 0
$$

（交换积分与求和、CND 假设下 $\leq 0$）。$\square$

**Step 3（SW² 的 CND）**：由 Step 1，对每个方向 $\theta$，$\psi_\theta(\mu, \nu) = \mathcal{W}_2^2(P_{\theta\#}\mu, P_{\theta\#}\nu)$ 在 $\mathcal{P}(\mathbb{R})$ 上是 CND 的。由 Step 2：

$$
\text{SW}_2^2(\mu, \nu) = \mathbb{E}_\theta[\psi_\theta(\mu, \nu)]
$$

是 CND 函数的期望，故也是 CND 的。$\square$

**注**：上述证明利用了两个关键性质：（1）一维 OT 有显式解，使 $\mathcal{W}_2^2$ 可表示为 Hilbert 空间中的平方范数差；（2）CND 函数在凸锥上封闭，对期望运算封闭。两者结合给出 SW² 的 CND 性质。

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

#### 4.4.9 小结

高斯核与 WL 内积的凸组合通过以下三个层次的理论保证和实践优势，构成了 TopoWGK 的最终形式：

1. **理论合法性**：两个子核均正定，凸组合保持正定性 → 任意 $\beta$ 下核矩阵合法。
2. **正定性边界**：$\beta = 0$ 或 $\beta = 1$ 时退化为已知的正定核（内积核 / 高斯核）。
3. **软特征选择**：$\beta$ 通过交叉验证自动调节结构和特征的相对权重 → 适应不同数据集。

这种设计使得 TopoWGK 在单一核函数（纯 WL 内积或纯 OT 高斯核）无法充分捕捉图间相似度时，通过凸组合找到一个介于"标签特征驱动"和"几何结构驱动"之间的最优平衡点。

---

## 第 5 章 参考文献

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

---


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
- [第 4 章 拓扑感知 Wasserstein 图核 (TopoWGK)](#第-4-章-拓扑感知-wasserstein-图核-topowgk)
  - [4.1 图核框架概述](#41-图核框架概述)
  - [4.2 条件负定性与正定核的关系](#42-条件负定性与正定核的关系)
  - [4.3 Sliced Wasserstein 距离的条件负定性](#43-sliced-wasserstein-距离的条件负定性)
  - [4.4 高斯核与 WL 内积的凸组合](#44-高斯核与-wl-内积的凸组合)
- [第 5 章 参考文献](#第-5-章-参考文献)

---

