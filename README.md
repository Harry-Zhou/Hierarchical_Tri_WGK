# Short-Long Topology-Aware Graph Kernel (Short_Long_Topology_Aware_WGK)

[English](#english) | [中文](#中文)

---

## English

### Overview

This project implements a **Topology-Aware Graph Kernel** that combines short-range and long-range topological information from graphs using Weisfeiler-Lehman (WL) propagation and Optimal Transport (OT) distance.

### Core Innovation

Traditional WL kernels only capture local neighborhood structure. This project extends WL to incorporate:

1. **Short-Range Topology**: Triangulated neighbors capturing local cycle structures
2. **Long-Range Topology**: Hierarchical cycle complexes capturing global cycle context

The combination is achieved via **Optimal Transport** to compute graph similarity.

### Project Structure

```
Short_Long_Topology_Aware_WGK/
├── hierarchical_cycle_complex_wl_tools/   # Core WL implementations
│   ├── classic_wl_tools.py               # Classic Weisfeiler-Lehman
│   ├── node_wl_via_cycle_complexes.py    # Node-level WL via cycle complexes
│   ├── edge_wl_via_cycle_complexes.py    # Edge-level WL via cycle complexes  
│   ├── node_wl_via_triangulated_neighbors.py  # TN-based node WL
│   ├── edge_wl_via_triangulated_neighbors.py  # TN-based edge WL
│   ├── hierarchical_cycle_complex_bfs_neighbors.py  # Cycle complex construction
│   └── topo_aware_wl_test.py            # Topology-aware WL test
├── our_experiments/
│   └── cycle_complex_wgk/                # Main kernel implementation
│       ├── topo_wasserstein_graph_kernel.py  # Topo-Wasserstein kernel
│       ├── main_topowgk.py              # Training & evaluation script
│       └── process_dataset.py           # Dataset processing
├── eval_grakel_baselines/                # Baseline comparisons (GraKeL)
└── utils.py                              # Utility functions
```

### Key Components

#### 1. Hierarchical Cycle Complex WL Tools

- **Classic WL**: Basic color refinement algorithm
- **Cycle Complex WL**: Extends WL with cycle basis information
- **Triangulated Neighbor WL**: Captures short cycles (chordal cycles)
- **Topological-Aware WL**: Combines TN (short) + HCC (long) iterations

#### 2. Topo-Wasserstein Graph Kernel

Main class: `TopoWassersteinGraphKernel`

- Computes WL label distributions at each iteration
- Uses degree distribution as node features
- Computes transport cost via Optimal Transport (POT library)
- Combines OT distance with WL inner product similarity

### Algorithms

#### Triangulated Neighbors (TN)
- Extracts chordal cycles from node neighborhoods
- Captures **short-range topological information**

#### Hierarchical Cycle Contexts (HCC)  
- Builds maximum cycle complexes via BFS on cycle graph
- Captures **long-range topological information**
- Layer-by-layer propagation of cycle context

#### Topo-Wasserstein Kernel
```
K(G₁, G₂) = α × OT_distance + (1-α) × WL_similarity
```
- `niter_tn`: Number of TN iterations
- `niter_hcc`: Number of HCC iterations  
- `α`: Interpolation coefficient

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

### Installation

```bash
# Required dependencies
pip install numpy networkx scipy ot  # Core
pip install scikit-learn pandas      # ML
pip install grakel                   # Baseline kernels
```

### Usage Example

```python
from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import TopoWassersteinGraphKernel
from our_experiments.cycle_complex_wgk.process_dataset import download_dataset, construct_dataset

# Load dataset
dataset_dict = download_dataset('PROTEINS')
graph_list, vlabel_list, edges_list, elabel_list, \
    vtx_hcc_list, vtx_tri_list, deg_distr_list, y = construct_dataset(dataset_dict)

# Compute kernel
kernel = TopoWassersteinGraphKernel(niter_tn=3, niter_hcc=3, wl_normalized=True)
ot_dist, wl_sim, runtime = kernel.fit_transform(
    dataset_dict['dataset_info'], graph_list, vlabel_list, 
    edges_list, elabel_list, vtx_hcc_list, vtx_tri_list, deg_distr_list
)
```

### References

- Weisfeiler, J. L., & Leman, A. (1968). The reduction of a graph to canonical form.
- Tuzhilina, E., et al. (2022). "Feature Selection and Kernel Design for Graph Machine Learning"
- GraKeL: A graph kernel library for machine learning

---

## 中文

### 项目概述

本项目实现了一个**拓扑感知图核 (Topology-Aware Graph Kernel)**，通过结合**短程**和**长程**拓扑信息，使用Weisfeiler-Lehman (WL) 传播和最优传输 (Optimal Transport) 距离来计算图相似度。

### 核心创新

传统WL核仅捕获局部邻域结构。本项目扩展WL算法以整合：

1. **短程拓扑**: 三角化邻域 (Triangulated Neighbors) - 捕获局部环结构
2. **长程拓扑**: 层级环复形 (Hierarchical Cycle Complexes) - 捕获全局环上下文

通过**最优传输**方法组合这些特征来计算图相似度。

### 项目结构

```
Short_Long_Topology_Aware_WGK/
├── hierarchical_cycle_complex_wl_tools/   # 核心WL实现
│   ├── classic_wl_tools.py               # 经典Weisfeiler-Lehman
│   ├── node_wl_via_cycle_complexes.py    # 基于环复形的节点WL
│   ├── edge_wl_via_cycle_complexes.py    # 基于环复形的边WL
│   ├── node_wl_via_triangulated_neighbors.py  # TN节点WL
│   ├── edge_wl_via_triangulated_neighbors.py  # TN边WL
│   ├── hierarchical_cycle_complex_bfs_neighbors.py  # 环复形构建
│   └── topo_aware_wl_test.py            # 拓扑感知WL测试
├── our_experiments/
│   └── cycle_complex_wgk/                # 主核实现
│       ├── topo_wasserstein_graph_kernel.py  # 拓扑Wasserstein图核
│       ├── main_topowgk.py              # 训练与评估脚本
│       └── process_dataset.py           # 数据集处理
├── eval_grakel_baselines/                # 基线对比 (GraKeL)
└── utils.py                              # 工具函数
```

### 核心组件

#### 1. 层级环复形WL工具

- **经典WL**: 基本颜色精化算法
- **环复形WL**: 基于环基扩展WL
- **三角化邻域WL**: 捕获短环 (弦环)
- **拓扑感知WL**: 组合TN (短程) + HCC (长程) 迭代

#### 2. 拓扑Wasserstein图核

主类: `TopoWassersteinGraphKernel`

- 计算每次迭代的WL标签分布
- 使用度分布作为节点特征
- 通过最优传输 (POT库) 计算传输成本
- 组合OT距离与WL内积相似度

### 算法原理

#### 三角化邻域 (TN)
- 从节点邻域中提取弦环
- 捕获**短程拓扑信息**

#### 层级环复形 (HCC)
- 通过环图上的BFS构建最大环复形
- 捕获**长程拓扑信息**
- 逐层传播环上下文

#### 拓扑Wasserstein核
```
K(G₁, G₂) = α × OT距离 + (1-α) × WL相似度
```
- `niter_tn`: TN迭代次数
- `niter_hcc`: HCC迭代次数  
- `α`: 插值系数

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

### 安装依赖

```bash
# 核心依赖
pip install numpy networkx scipy ot  # 核心库
pip install scikit-learn pandas      # 机器学习
pip install grakel                   # 基线图核
```

### 使用示例

```python
from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import TopoWassersteinGraphKernel
from our_experiments.cycle_complex_wgk.process_dataset import download_dataset, construct_dataset

# 加载数据集
dataset_dict = download_dataset('PROTEINS')
graph_list, vlabel_list, edges_list, elabel_list, \
    vtx_hcc_list, vtx_tri_list, deg_distr_list, y = construct_dataset(dataset_dict)

# 计算核矩阵
kernel = TopoWassersteinGraphKernel(niter_tn=3, niter_hcc=3, wl_normalized=True)
ot_dist, wl_sim, runtime = kernel.fit_transform(
    dataset_dict['dataset_info'], graph_list, vlabel_list, 
    edges_list, elabel_list, vtx_hcc_list, vtx_tri_list, deg_distr_list
)
```

### 参考文献

- Weisfeiler, J. L., & Leman, A. (1968). The reduction of a graph to canonical form.
- Tuzhilina, E., et al. (2022). "Feature Selection and Kernel Design for Graph Machine Learning"
- GraKeL: A graph kernel library for machine learning

---

### License

MIT License

### Contact

For questions and contributions, please open an issue.