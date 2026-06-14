# Hierarchical Triangulated WL Graph Kernel (Hierarchical_Tri_WGK)

[English](#english) | [中文](#中文)

---

## English

### Overview

This project implements a **Topology-Aware Graph Kernel** that combines short-range and long-range topological information from graphs using Weisfeiler-Lehman (WL) propagation and Optimal Transport (OT) distance.

### Core Innovation

Traditional WL kernels only capture local neighborhood structure. This project extends WL to incorporate:

1. **Short-Range Topology**: Triangulated neighbors (TN) capturing local chordal cycle structures
2. **Long-Range Topology**: Hierarchical Cycle Complexes (HCC) capturing global cycle context via BFS on the cycle graph

The combination is achieved via **Optimal Transport** to compute graph similarity.

### Project Structure

```
Hierarchical_Tri_WGK/
├── cyclic_schema/                       # Cyclic schema & hierarchical WL engine
│   ├── cyclic_schema.py                 # CSG construction, cycle basis, multi-layer CSG, mappings
│   ├── hierarchical_triangular_wl.py  # Hierarchical WL (TN + HCC unified dispatch)
│   ├── test_hierarchical_wl.py          # Comprehensive tests (100+ test cases)
│   ├── pyproject.toml                   # Package configuration
│   ├── cfi_3wl_analysis.md              # CFI 3-WL analysis documentation
│   ├── cfi_comparison.png               # CFI comparison figure
│   └── theoretic_analysis.md            # Theoretical background
├── hierarchical_tri_wl_tools/           # Core WL implementations
│   ├── classic_wl_tools.py                       # Classic Weisfeiler-Lehman label propagation
│   ├── classic_edge_wl_tools.py                  # Edge-aware classic WL
│   ├── local_structure_wl_tools.py               # Local structure propagation helpers
│   ├── local_structure_edge_wl_tools.py          # Edge-aware local structure helpers
│   ├── node_wl_via_triangulated_neighbors.py     # TN-based node WL
│   ├── edge_wl_via_triangulated_neighbors.py     # TN-based edge WL
│   ├── node_wl_via_cycle_complexes.py            # Cycle-complex-based node WL
│   ├── edge_wl_via_cycle_complexes.py            # Cycle-complex-based edge WL
│   ├── hierarchical_cycle_complex_bfs_neighbors.py  # BFS-based HCC construction
│   └── topo_aware_wl_test.py                     # Unified topology-aware WL test
├── our_experiments/
│   ├── cycle_complex_wgk/                # Main kernel implementation
│   │   ├── topo_wasserstein_graph_kernel.py  # Topo-Wasserstein graph kernel
│   │   ├── main_topowgk.py              # SVM training & evaluation pipeline
│   │   ├── process_dataset.py           # Dataset download & construction
│   │   ├── generate_strongly_regular_graphs.py  # SRG experiments
│   │   ├── graph_isomorphism_test.py     # Graph isomorphism testing
│   │   ├── hyperparameter_sensitivity.py  # Hyperparameter analysis
│   │   ├── postprocess_svm_results.py    # Results aggregation
│   │   ├── test_edge_label_integration.py  # Edge label integration tests
│   │   ├── test_sage.py                 # SageMath integration test
│   │   └── topowgk_outputs/             # Evaluation outputs
│   ├── cfi_test.py                      # CFI graph experiments
│   └── example.py                       # Quick-start example
├── eval_grakel_baselines/                # Baseline graph kernel comparisons
│   ├── test_grakel.py
│   ├── test_grakel_baselines.py
│   ├── test_grakel_baselines_scripts.sh
│   ├── test_wgk.py
│   ├── test_wgk_script.sh
│   ├── wgk_utils.py
│   └── outputs/                         # Baseline results
├── tests/                               # Unit & integration tests
│   ├── conftest.py                      # Shared pytest fixtures
│   ├── test_cyclic_schema.py            # CSG construction tests
│   ├── test_integration.py              # End-to-end pipeline tests
│   ├── test_kernel.py                   # Kernel fit/transform tests
│   ├── test_utils.py                    # Utility function tests
│   ├── test_wl_tools.py                 # WL tooling tests
│   └── pytest.ini                       # Test configuration
├── utils.py                             # Shared utilities
├── README.md                            # This file (bilingual)
└── README.en.md                         # English-only README
```

### Key Components

#### 1. Cyclic Schema (CSG)

The `cyclic_schema` package converts an input graph into a **Cyclic Schematic Graph (CSG)** — a higher-level abstraction where each cycle in the input graph becomes a node. This enables multi-layer hierarchical analysis:

- **Single-layer CSG** (L=1): One level of cycle abstraction
- **Multi-layer CSG** (L>1): Recursively build CSG on CSGs for deeper abstraction

The CSG construction identifies cycle bases, classifies nodes (cycle / non-cycle / interface), and builds input-to-CSG and CSG-to-input mappings for label propagation.

#### 2. Hierarchical WL Engine

`hierarchical_triangular_wl.py` is the central orchestrator combining two propagation strategies:

- **Triangulated Neighbors (TN)**: Short-range — looks at cycles formed within a node's neighbor set (chordal cycles), iteration count `I`
- **Hierarchical Cycle Complexes (HCC)**: Long-range — propagates labels through the multi-layer CSG, depth `K`

Supports both **node-only** and **edge-label-aware** propagation via a unified dispatch interface.

#### 3. Core WL Tools

Low-level WL implementations under `hierarchical_tri_wl_tools/`:

- **Classic WL**: Standard color refinement (vertex label → sorted neighbor labels → compression)
- **TN WL**: Propagate labels through triangulated (chordal cycle) neighborhoods
- **Cycle Complex WL**: Propagate labels through hierarchical cycle contexts
- **Edge-aware variants**: All WL variants support edge label propagation

#### 4. Topo-Wasserstein Graph Kernel

Main class: `TopoWassersteinGraphKernel`

```
K(G₁, G₂) = α × OT_distance + (1-α) × WL_similarity
```

- **WL similarity**: Inner product of WL label distribution matrices
- **OT distance**: Wasserstein distance between WL feature distributions (POT library)
- **Normalization**: Optional cosine normalization for WL similarity matrix

**Parameters:**
- `niter_tn`: Number of TN (short-range) WL iterations
- `niter_hcc`: Number of HCC (long-range) CSG layers

### Installation

```bash
# Core dependencies
pip install numpy networkx scipy pot  # Core
pip install scikit-learn pandas       # ML
pip install grakel                    # Baseline kernels (optional)
```

### Quick Usage

```python
from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import (
    TopoWassersteinGraphKernel
)
import networkx as nx
import numpy as np

# Build two small graphs
G1 = nx.Graph()
G1.add_edges_from([(0, 1), (1, 2), (0, 2)])
G2 = nx.Graph()
G2.add_edges_from([(0, 1), (1, 2)])

vlabel1 = np.array([1, 2, 3], dtype=np.int32)
vlabel2 = np.array([1, 2, 3], dtype=np.int32)

def deg_distr(g):
    d = np.array([g.degree(v) for v in g.nodes()], dtype=np.float32)
    return d / d.sum()

kernel = TopoWassersteinGraphKernel(niter_tn=3, niter_hcc=3, wl_normalized=True)
ot_dist, wl_sim, runtime = kernel.fit_transform(
    {'nl': True, 'el': False},
    [G1, G2],
    [vlabel1, vlabel2],
    [[(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]],
    [np.array([]), np.array([])],
    [deg_distr(G1), deg_distr(G2)],
)
```

### Running Tests

```bash
# Install test dependencies
pip install pytest networkx numpy

# Run all tests
pytest tests/ cyclic_schema/test_hierarchical_wl.py -v

# 154 test cases covering:
# - CSG construction & multi-layer abstraction
# - Classic WL propagation
# - Edge-aware WL propagation
# - Kernel fit/transform
# - End-to-end integration
```

### Datasets Evaluated

| Dataset | Type | Description |
|---------|------|-------------|
| IMDB-MULTI | Social Networks | Multi-class movie collaboration |
| IMDB-BINARY | Social Networks | Binary movie collaboration |
| PROTEINS | Bioinformatics | Protein structure classification |
| ENZYMES | Bioinformatics | Enzyme classification |
| NCI1 | Bioinformatics | Anti-cancer activity prediction |
| NCI109 | Bioinformatics | Anti-cancer activity prediction |
| Mutagenicity | Chemistry | Molecular mutagenicity |
| DHFR | Chemistry | Dihydrofolate reductase inhibition |
| BZR, COX2 | Chemistry | Molecular property prediction |
| MSRC_21, MSRC_9 | Image | Graph-based image segmentation |
| FIRSTMM_DB | Social Networks | Software bug prediction |
| SYNTHETIC | Synthetic | Strongly regular graphs, CFI graphs |

### References

- Weisfeiler, J. L., & Leman, A. (1968). The reduction of a graph to canonical form.
- Tuzhilina, E., et al. (2022). "Feature Selection and Kernel Design for Graph Machine Learning"
- Flam-Sanchez, et al. (2022). "Short-Long Topology-Aware Graph Kernels"
- GraKeL: A graph kernel library for machine learning

---

## 中文

### 项目概述

本项目实现了一个**拓扑感知图核 (Topology-Aware Graph Kernel)**，通过结合**短程**和**长程**拓扑信息，使用 Weisfeiler-Lehman (WL) 传播和最优传输 (Optimal Transport) 距离来计算图相似度。

### 核心创新

传统 WL 核仅捕获局部邻域结构。本项目扩展 WL 算法以整合：

1. **短程拓扑**: 三角化邻域 (Triangulated Neighbors, TN) — 捕获节点邻域内的弦环结构
2. **长程拓扑**: 层级环复形 (Hierarchical Cycle Complexes, HCC) — 通过环图上的 BFS 捕获全局环上下文

通过**最优传输**方法组合这些特征来计算图相似度。

### 项目结构

```
Hierarchical_Tri_WGK/
├── cyclic_schema/                       # 环式图构建和层级WL引擎
│   ├── cyclic_schema.py                 # CSG 构建、环基、多层CSG、映射
│   ├── hierarchical_triangular_wl.py  # 层级 WL（TN + HCC 统一调度的核心）
│   ├── test_hierarchical_wl.py          # 详尽测试（100+ 测试用例）
│   ├── cfi_3wl_analysis.md              # CFI 3-WL 分析文档
│   ├── cfi_comparison.png               # CFI 对比图
│   └── theoretic_analysis.md            # 理论分析
├── hierarchical_tri_wl_tools/           # 核心 WL 实现
│   ├── classic_wl_tools.py                       # 经典 Weisfeiler-Lehman
│   ├── classic_edge_wl_tools.py                  # 边感知经典 WL
│   ├── local_structure_wl_tools.py               # 局部结构传播辅助
│   ├── local_structure_edge_wl_tools.py          # 边感知局部结构传播辅助
│   ├── node_wl_via_triangulated_neighbors.py     # TN 节点 WL
│   ├── edge_wl_via_triangulated_neighbors.py     # TN 边 WL
│   ├── node_wl_via_cycle_complexes.py            # 环复形节点 WL
│   ├── edge_wl_via_cycle_complexes.py            # 环复形边 WL
│   ├── hierarchical_cycle_complex_bfs_neighbors.py  # BFS 构建 HCC
│   └── topo_aware_wl_test.py                     # 统一拓扑感知 WL 测试
├── our_experiments/
│   ├── cycle_complex_wgk/                # 主核实现
│   │   ├── topo_wasserstein_graph_kernel.py  # 拓扑 Wasserstein 图核
│   │   ├── main_topowgk.py              # SVM 训练与评估流程
│   │   ├── process_dataset.py           # 数据集下载与构建
│   │   ├── generate_strongly_regular_graphs.py  # 强正则图实验
│   │   ├── graph_isomorphism_test.py     # 图同构测试
│   │   ├── hyperparameter_sensitivity.py  # 超参数分析
│   │   ├── postprocess_svm_results.py    # 结果汇总
│   │   ├── test_edge_label_integration.py  # 边标签集成测试
│   │   ├── test_sage.py                 # SageMath 集成测试
│   │   └── topowgk_outputs/             # 评估结果
│   ├── cfi_test.py                      # CFI 图实验
│   └── example.py                       # 快速入门示例
├── eval_grakel_baselines/                # 基线图核对比
│   ├── test_grakel.py
│   ├── test_grakel_baselines.py
│   ├── test_grakel_baselines_scripts.sh
│   ├── test_wgk.py
│   └── outputs/
├── tests/                               # 单元测试和集成测试
│   ├── conftest.py                      # 共享 pytest fixtures
│   ├── test_cyclic_schema.py
│   ├── test_integration.py
│   ├── test_kernel.py
│   ├── test_utils.py
│   ├── test_wl_tools.py
│   └── pytest.ini
├── utils.py                             # 共享工具函数
├── README.md                            # 本文件（双语）
└── README.en.md                         # 英文版 README
```

### 核心组件

#### 1. 环式图 (Cyclic Schema, CSG)

`cyclic_schema` 包将输入图转换为**环式图 (CSG)** — 这是一种更高层次的抽象，其中输入图中的每个环都成为一个节点。这使得多层层级分析成为可能：

- **单层 CSG** (L=1)：一层环抽象
- **多层 CSG** (L>1)：递归地在 CSG 上构建 CSG，获得更深层次的抽象

CSG 构建过程包括识别环基、分类节点（环上节点/非环节点/接口节点），以及构建输入图到 CSG 的双向映射，用于标签传播。

#### 2. 层级 WL 引擎

`hierarchical_triangular_wl.py` 是核心编排器，组合两种传播策略：

- **三角化邻域 (TN)**：短程 — 查看节点邻域内形成的环（弦环），迭代次数 `I`
- **层级环复形 (HCC)**：长程 — 通过多层 CSG 传播标签，深度 `K`

通过统一调度接口支持**纯节点**和**边标签感知**两种传播模式。

#### 3. 核心 WL 工具

`hierarchical_tri_wl_tools/` 下的底层 WL 实现：

- **经典 WL**：标准颜色精化（顶点标签 → 排序的邻居标签 → 压缩）
- **TN WL**：通过三角化邻域传播标签
- **环复形 WL**：通过层级环上下文传播标签
- **边感知变体**：所有 WL 变体均支持边标签传播

#### 4. 拓扑 Wasserstein 图核

主类: `TopoWassersteinGraphKernel`

```
K(G₁, G₂) = α × OT 距离 + (1-α) × WL 相似度
```

- **WL 相似度**：WL 标签分布矩阵的内积
- **OT 距离**：WL 特征分布间的 Wasserstein 距离（POT 库）
- **归一化**：可选余弦归一化 WL 相似度矩阵

**参数：**
- `niter_tn`: TN（短程）WL 迭代次数
- `niter_hcc`: HCC（长程）CSG 层数

### 安装依赖

```bash
# 核心依赖
pip install numpy networkx scipy pot  # 核心库
pip install scikit-learn pandas       # 机器学习
pip install grakel                    # 基线图核（可选）
```

### 快速使用

```python
from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import (
    TopoWassersteinGraphKernel
)
import networkx as nx
import numpy as np

# 构建两个小图
G1 = nx.Graph()
G1.add_edges_from([(0, 1), (1, 2), (0, 2)])
G2 = nx.Graph()
G2.add_edges_from([(0, 1), (1, 2)])

vlabel1 = np.array([1, 2, 3], dtype=np.int32)
vlabel2 = np.array([1, 2, 3], dtype=np.int32)

def deg_distr(g):
    d = np.array([g.degree(v) for v in g.nodes()], dtype=np.float32)
    return d / d.sum()

kernel = TopoWassersteinGraphKernel(niter_tn=3, niter_hcc=3, wl_normalized=True)
ot_dist, wl_sim, runtime = kernel.fit_transform(
    {'nl': True, 'el': False},
    [G1, G2],
    [vlabel1, vlabel2],
    [[(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2)]],
    [np.array([]), np.array([])],
    [deg_distr(G1), deg_distr(G2)],
)
```

### 运行测试

```bash
# 安装测试依赖
pip install pytest networkx numpy

# 运行所有测试
pytest tests/ cyclic_schema/test_hierarchical_wl.py -v

# 共 154 个测试用例，覆盖：
# - CSG 构建与多层抽象
# - 经典 WL 传播
# - 边感知 WL 传播
# - 核的 fit/transform
# - 端到端集成测试
```

### 评估数据集

| 数据集 | 类型 | 描述 |
|--------|------|------|
| IMDB-MULTI | 社交网络 | 多类电影合作 |
| IMDB-BINARY | 社交网络 | 二类电影合作 |
| PROTEINS | 生物信息学 | 蛋白质结构分类 |
| ENZYMES | 生物信息学 | 酶分类 |
| NCI1 | 生物信息学 | 抗癌活性预测 |
| NCI109 | 生物信息学 | 抗癌活性预测 |
| Mutagenicity | 化学 | 分子致突变性 |
| DHFR | 化学 | 二氢叶酸还原酶抑制 |
| BZR, COX2 | 化学 | 分子性质预测 |
| MSRC_21, MSRC_9 | 图像 | 基于图的图像分割 |
| FIRSTMM_DB | 社交网络 | 软件缺陷预测 |
| SYNTHETIC | 合成 | 强正则图、CFI 图 |

### 参考文献

- Weisfeiler, J. L., & Leman, A. (1968). The reduction of a graph to canonical form.
- Tuzhilina, E., et al. (2022). "Feature Selection and Kernel Design for Graph Machine Learning"
- Flam-Sanchez, et al. (2022). "Short-Long Topology-Aware Graph Kernels"
- GraKeL: A graph kernel library for machine learning

---

### License

MIT License

### Contact

For questions and contributions, please open an issue.
