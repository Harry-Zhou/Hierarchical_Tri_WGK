# Short-Long Topology-Aware Graph Kernel (Short_Long_Topology_Aware_WGK)

## Overview

This project implements a **Topology-Aware Graph Kernel** that combines short-range and long-range topological information from graphs using Weisfeiler-Lehman (WL) propagation and Optimal Transport (OT) distance.

## Core Innovation

Traditional WL kernels only capture local neighborhood structure. This project extends WL to incorporate:

1. **Short-Range Topology**: Triangulated neighbors capturing local cycle structures
2. **Long-Range Topology**: Hierarchical cycle complexes capturing global cycle context

The combination is achieved via **Optimal Transport** to compute graph similarity.

## Project Structure

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

## Key Components

### 1. Hierarchical Cycle Complex WL Tools

- **Classic WL**: Basic color refinement algorithm
- **Cycle Complex WL**: Extends WL with cycle basis information
- **Triangulated Neighbor WL**: Captures short cycles (chordal cycles)
- **Topological-Aware WL**: Combines TN (short) + HCC (long) iterations

### 2. Topo-Wasserstein Graph Kernel

Main class: `TopoWassersteinGraphKernel`

- Computes WL label distributions at each iteration
- Uses degree distribution as node features
- Computes transport cost via Optimal Transport (POT library)
- Combines OT distance with WL inner product similarity

## Algorithms

### Triangulated Neighbors (TN)
- Extracts chordal cycles from node neighborhoods
- Captures **short-range topological information**

### Hierarchical Cycle Contexts (HCC)  
- Builds maximum cycle complexes via BFS on cycle graph
- Captures **long-range topological information**
- Layer-by-layer propagation of cycle context

### Topo-Wasserstein Kernel
```
K(G₁, G₂) = α × OT_distance + (1-α) × WL_similarity
```
- `niter_tn`: Number of TN iterations
- `niter_hcc`: Number of HCC iterations  
- `α`: Interpolation coefficient

## Datasets Evaluated

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

## Installation

```bash
# Required dependencies
pip install numpy networkx scipy ot  # Core
pip install scikit-learn pandas      # ML
pip install grakel                   # Baseline kernels
```

## Usage Example

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

## Key Features

- **Hybrid Topology**: Combines short-range (triangulated neighbors) and long-range (cycle complexes) topological features
- **Hierarchical Propagation**: Multi-layer cycle context propagation via BFS
- **Optimal Transport**: Uses Wasserstein distance for graph-level similarity
- **Flexible Configuration**: Adjustable iteration counts for TN and HCC (e.g., niter_tn=3, niter_hcc varies 0-7)

## References

- Weisfeiler, J. L., & Leman, A. (1968). The reduction of a graph to canonical form.
- Tuzhilina, E., et al. (2022). "Feature Selection and Kernel Design for Graph Machine Learning"
- GraKeL: A graph kernel library for machine learning

## License

MIT License

## Contact

For questions and contributions, please open an issue.