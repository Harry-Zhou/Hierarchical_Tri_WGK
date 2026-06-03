# Short-Long Topology-Aware Graph Kernel (Short_Long_Topology_Aware_WGK)

## Overview

This project implements a **Topology-Aware Graph Kernel** that combines short-range and long-range topological information from graphs using Weisfeiler-Lehman (WL) propagation and Optimal Transport (OT) distance.

## Core Innovation

Traditional WL kernels only capture local neighborhood structure. This project extends WL to incorporate:

1. **Short-Range Topology**: Triangulated neighbors (TN) capturing local chordal cycle structures
2. **Long-Range Topology**: Hierarchical Cycle Complexes (HCC) capturing global cycle context via BFS on the cycle graph

The combination is achieved via **Optimal Transport** to compute graph similarity.

## Project Structure

```
Short_Long_Topology_Aware_WGK/
├── cyclic_schema/                       # Cyclic schema & hierarchical WL engine
│   ├── cyclic_schema.py                 # CSG construction, cycle basis, multi-layer CSG, mappings
│   ├── hierarchical_triangulated_wl.py  # Hierarchical WL (TN + HCC unified dispatch)
│   ├── test_hierarchical_wl.py          # Comprehensive tests (100+ test cases)
│   ├── pyproject.toml                   # Package configuration
│   ├── cfi_3wl_analysis.md              # CFI 3-WL analysis
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
│   └── example.py                       # Usage example
├── eval_grakel_baselines/                # Baseline graph kernel comparisons
│   ├── test_grakel.py
│   ├── test_grakel_baselines.py
│   ├── test_grakel_baselines_scripts.sh
│   ├── test_wgk.py
│   ├── test_wgk_script.sh
│   ├── wgk_utils.py
│   └── outputs/
├── tests/                               # Unit & integration tests
│   ├── conftest.py                      # Shared pytest fixtures
│   ├── test_cyclic_schema.py            # CSG construction tests
│   ├── test_integration.py              # End-to-end pipeline tests
│   ├── test_kernel.py                   # Kernel fit/transform tests
│   ├── test_utils.py                    # Utility function tests
│   ├── test_wl_tools.py                 # WL tooling tests
│   └── pytest.ini                       # Test configuration
├── utils.py                             # Shared utilities
├── README.md                            # Bilingual README
└── README.en.md                         # This file
```

## Key Components

### 1. Cyclic Schema (CSG)

The `cyclic_schema` package converts an input graph into a **Cyclic Schematic Graph (CSG)** — a higher-level abstraction where each cycle in the input graph becomes a node. This enables multi-layer hierarchical analysis:

- **Single-layer CSG** (K=1): One level of cycle abstraction
- **Multi-layer CSG** (K>1): Recursively build CSG on CSGs for deeper abstraction

The CSG construction identifies cycle bases, classifies nodes (cycle / non-cycle / interface), and builds input-to-CSG and CSG-to-input mappings for label propagation.

### 2. Hierarchical WL Engine

`hierarchical_triangulated_wl.py` is the central orchestrator combining two propagation strategies:

- **Triangulated Neighbors (TN)**: Short-range — looks at cycles formed within a node's neighbor set (chordal cycles), iteration count `I`
- **Hierarchical Cycle Complexes (HCC)**: Long-range — propagates labels through the multi-layer CSG, depth `K`

Supports both **node-only** and **edge-label-aware** propagation via a unified dispatch interface.

### 3. Core WL Tools

Low-level WL implementations under `hierarchical_tri_wl_tools/`:

- **Classic WL**: Standard color refinement (vertex label → sorted neighbor labels → compression)
- **TN WL**: Propagate labels through triangulated (chordal cycle) neighborhoods
- **Cycle Complex WL**: Propagate labels through hierarchical cycle contexts
- **Edge-aware variants**: All WL variants support edge label propagation

### 4. Topo-Wasserstein Graph Kernel

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

## Installation

```bash
# Core dependencies
pip install numpy networkx scipy pot  # Core
pip install scikit-learn pandas       # ML
pip install grakel                    # Baseline kernels (optional)
```

## Quick Usage

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

## Running Tests

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
| SYNTHETIC | Synthetic | Strongly regular graphs, CFI graphs |

## Key Features

- **Hybrid Topology**: Combines short-range (triangulated neighbors) and long-range (cycle complexes) topological features
- **Hierarchical Propagation**: Multi-layer cycle context propagation via BFS
- **Optimal Transport**: Uses Wasserstein distance for graph-level similarity
- **Flexible Configuration**: Adjustable iteration counts for TN and HCC
- **Edge Label Support**: All WL variants support edge label propagation

## References

- Weisfeiler, J. L., & Leman, A. (1968). The reduction of a graph to canonical form.
- Tuzhilina, E., et al. (2022). "Feature Selection and Kernel Design for Graph Machine Learning"
- Flam-Sanchez, et al. (2022). "Short-Long Topology-Aware Graph Kernels"
- GraKeL: A graph kernel library for machine learning

## License

MIT License

## Contact

For questions and contributions, please open an issue.
