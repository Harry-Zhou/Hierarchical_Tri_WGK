"""
Multi-layer Cyclic Schematic Graph (CSG) construction.

This module provides a unified interface for building multi-layer CSG
hierarchies, designed for high cohesion, low coupling, and flat architecture.

The module supports:
1. Single-graph multi-layer CSG construction
2. Paired-graph multi-layer CSG construction (for WL/CSG-Transformer comparison)
3. Precomputation of neighbor components (topology-dependent, label-independent)
4. Consistent interface for both discrete WL and continuous Transformer message passing

Architecture:
- This module is the shared foundation for:
  * cyclic_schema/htn_wl.py (discrete WL message passing)
  * cyclic_schema/htn_wl.py (legacy alias, kept for backward compatibility)
  * cyclic_schema/csg_transformer.py (continuous Transformer message passing)
- Both modules import this module to get multi-layer CSG structures,
  avoiding code duplication and ensuring consistency.
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Hashable, List, Optional, Tuple

import networkx as nx

from .cyclic_schema import (
    build_csg_to_input_mapping,
    build_input_to_csg_mapping,
    build_multilayer_csg,
    cyclic_schematic_graph,
    node_sort_key,
)


# ============================================================================
# Data structures
# ============================================================================

@dataclass
class SingleGraphCSG:
    """
    Multi-layer CSG structure for a single graph.
    
    This is the unified data structure returned by build_multilayer_csg_single().
    It contains all information needed for message passing (WL or Transformer).
    
    Attributes
    ----------
    layers : list
        List of (H, cycle_basis, info) tuples, one per CSG layer.
        layers[0] is CSG^1, layers[1] is CSG^2, etc.
    mappings : list
        List of lower_to_csg mappings (dict: lower_node -> tuple of CSG nodes).
        mappings[k] maps from layer k to layer k+1 (where layer 0 = input graph).
    neighbor_components : list
        List of neighbor component dicts (node -> tuple of component tuples).
        neighbor_components[k] is for layer k (layer 0 = input graph).
    input_graph : nx.Graph
        The original input graph.
    L : int
        Number of CSG layers.
    """
    layers: List[Tuple[nx.Graph, List, Dict]]
    mappings: List[Dict[Hashable, Any]]
    neighbor_components: List[Dict[Hashable, Tuple]]
    input_graph: nx.Graph
    L: int


@dataclass
class PairedGraphCSG:
    """
    Multi-layer CSG structure for a pair of graphs (for comparison).
    
    This is the unified data structure returned by build_multilayer_csg_pair().
    It contains all information needed for joint message passing on two graphs.
    
    Attributes
    ----------
    csg1 : SingleGraphCSG
        Multi-layer CSG structure for graph 1.
    csg2 : SingleGraphCSG
        Multi-layer CSG structure for graph 2.
    """
    csg1: SingleGraphCSG
    csg2: SingleGraphCSG


# ============================================================================
# Core construction functions
# ============================================================================

def _precompute_neighbor_components(G: nx.Graph) -> Dict[Hashable, Tuple]:
    """
    Precompute neighbor components for all nodes in graph G.
    
    The induced-subgraph components of N(v) depend only on graph topology,
    not on node labels, so this is computed once and reused across all
    WL iterations.
    
    Parameters
    ----------
    G : nx.Graph
        Input graph.
    
    Returns
    -------
    dict
        Mapping from node to tuple of component tuples.
        Each component tuple contains neighbor nodes in deterministic order.
    """
    adj = {n: list(G.neighbors(n)) for n in G.nodes()}
    result: Dict[Hashable, Tuple] = {}
    
    for v, neighbors in adj.items():
        if not neighbors:
            result[v] = ()
            continue
        
        nbr_set = frozenset(neighbors)
        seen: set = set()
        components: List[Tuple] = []
        
        for u in neighbors:
            if u in seen:
                continue
            comp = [u]
            seen.add(u)
            stack = [u]
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if y in nbr_set and y not in seen:
                        seen.add(y)
                        comp.append(y)
                        stack.append(y)
            # Sort component nodes deterministically
            comp.sort(key=node_sort_key)
            components.append(tuple(comp))
        
        result[v] = tuple(components)
    
    return result


def build_multilayer_csg_single(
    G: nx.Graph,
    L: int = 1,
    cb_prefix_func: Optional[Callable[[int], str]] = None,
) -> SingleGraphCSG:
    """
    Build multi-layer CSG structure for a single graph.
    
    This is the primary interface for constructing multi-layer CSG
    hierarchies. It builds the CSG layers, inter-layer mappings, and
    precomputes neighbor components for all layers.
    
    Parameters
    ----------
    G : nx.Graph
        Input graph.
    L : int
        Number of CSG layers (default: 1).
    cb_prefix_func : callable, optional
        Function that takes layer index (1-based) and returns prefix.
        Default: lambda i: f"cb{i}"
    
    Returns
    -------
    SingleGraphCSG
        Complete multi-layer CSG structure with all necessary information.
    
    Examples
    --------
    >>> G = nx.cycle_graph(5)
    >>> csg = build_multilayer_csg_single(G, L=2)
    >>> len(csg.layers)  # 2 layers
    2
    >>> len(csg.mappings)  # 2 mappings (input->CSG^1, CSG^1->CSG^2)
    2
    >>> len(csg.neighbor_components)  # 3 graphs (input + 2 CSG layers)
    3
    """
    if cb_prefix_func is None:
        cb_prefix_func = lambda i: f"cb{i}"
    
    # Build CSG layers
    layers = build_multilayer_csg(G, L, cb_prefix_func)
    
    # Build inter-layer mappings
    mappings = []
    current_lower = G
    
    for k in range(L):
        H_k, cb_k, _ = layers[k]
        lower_to_csg = build_input_to_csg_mapping(H_k, cb_k, current_lower)
        mappings.append(lower_to_csg)
        current_lower = H_k
    
    # Precompute neighbor components for all layers
    neighbor_components = [_precompute_neighbor_components(G)]
    for k in range(L):
        H_k, _, _ = layers[k]
        neighbor_components.append(_precompute_neighbor_components(H_k))
    
    return SingleGraphCSG(
        layers=layers,
        mappings=mappings,
        neighbor_components=neighbor_components,
        input_graph=G,
        L=L,
    )


def build_multilayer_csg_pair(
    G1: nx.Graph,
    G2: nx.Graph,
    L: int = 1,
    cb_prefix_func: Optional[Callable[[int], str]] = None,
) -> PairedGraphCSG:
    """
    Build multi-layer CSG structure for a pair of graphs.
    
    This is the primary interface for constructing multi-layer CSG
    hierarchies for graph comparison tasks (WL test, CSG-Transformer, etc.).
    
    Parameters
    ----------
    G1, G2 : nx.Graph
        Input graphs.
    L : int
        Number of CSG layers (default: 1).
    cb_prefix_func : callable, optional
        Function that takes layer index (1-based) and returns prefix.
        Default: lambda i: f"cb{i}"
    
    Returns
    -------
    PairedGraphCSG
        Complete multi-layer CSG structure for both graphs.
    
    Examples
    --------
    >>> G1 = nx.cycle_graph(5)
    >>> G2 = nx.path_graph(5)
    >>> paired = build_multilayer_csg_pair(G1, G2, L=1)
    >>> paired.csg1.L  # Number of layers
    1
    """
    csg1 = build_multilayer_csg_single(G1, L, cb_prefix_func)
    csg2 = build_multilayer_csg_single(G2, L, cb_prefix_func)
    
    return PairedGraphCSG(csg1=csg1, csg2=csg2)


# ============================================================================
# Utility functions for accessing CSG structure
# ============================================================================

def get_layer_graph(csg: SingleGraphCSG, layer_idx: int) -> nx.Graph:
    """
    Get the graph at a specific layer.
    
    Parameters
    ----------
    csg : SingleGraphCSG
        Multi-layer CSG structure.
    layer_idx : int
        Layer index (0 = input graph, 1 = CSG^1, etc.).
    
    Returns
    -------
    nx.Graph
        The graph at the specified layer.
    """
    if layer_idx == 0:
        return csg.input_graph
    else:
        return csg.layers[layer_idx - 1][0]


def get_cycle_basis(csg: SingleGraphCSG, layer_idx: int) -> List:
    """
    Get the cycle basis at a specific CSG layer.
    
    Parameters
    ----------
    csg : SingleGraphCSG
        Multi-layer CSG structure.
    layer_idx : int
        CSG layer index (1-based, 1 = CSG^1, etc.).
    
    Returns
    -------
    list
        Cycle basis at the specified layer.
    
    Raises
    ------
    IndexError
        If layer_idx is 0 (input graph has no cycle basis in CSG structure).
    """
    if layer_idx == 0:
        raise IndexError("Input graph (layer 0) has no cycle basis in CSG structure")
    return csg.layers[layer_idx - 1][1]


def get_mapping(csg: SingleGraphCSG, layer_idx: int) -> Dict[Hashable, Any]:
    """
    Get the lower-to-CSG mapping at a specific layer.
    
    Parameters
    ----------
    csg : SingleGraphCSG
        Multi-layer CSG structure.
    layer_idx : int
        Layer index (1 = input->CSG^1, 2 = CSG^1->CSG^2, etc.).
    
    Returns
    -------
    dict
        Mapping from lower-layer nodes to tuple of CSG nodes.
    """
    return csg.mappings[layer_idx - 1]


def get_neighbor_components(csg: SingleGraphCSG, layer_idx: int) -> Dict[Hashable, Tuple]:
    """
    Get precomputed neighbor components at a specific layer.
    
    Parameters
    ----------
    csg : SingleGraphCSG
        Multi-layer CSG structure.
    layer_idx : int
        Layer index (0 = input graph, 1 = CSG^1, etc.).
    
    Returns
    -------
    dict
        Neighbor components for the specified layer.
    """
    return csg.neighbor_components[layer_idx]


# ============================================================================
# Backward compatibility helpers
# ============================================================================

def extract_htn_wl_args(
    paired_csg: PairedGraphCSG,
) -> Tuple[
    List[Tuple[nx.Graph, List, Dict]],
    List[Tuple[nx.Graph, List, Dict]],
    List[Dict[Hashable, Any]],
    List[Dict[Hashable, Any]],
    Dict[Hashable, Tuple],
    Dict[Hashable, Tuple],
    List[Tuple[Dict[Hashable, Tuple], Dict[Hashable, Tuple]]],
]:
    """
    Extract arguments in the format expected by htn_wl functions.
    
    This function provides backward compatibility for existing htn_wl.py code.
    It converts the new PairedGraphCSG structure to the old tuple format.
    
    Parameters
    ----------
    paired_csg : PairedGraphCSG
        Paired multi-layer CSG structure.
    
    Returns
    -------
    tuple
        (layers1, layers2, mappings1, mappings2, nc_input1, nc_input2, nc_csg)
        in the format expected by htn_wl functions.
    """
    csg1 = paired_csg.csg1
    csg2 = paired_csg.csg2
    
    # Extract layers
    layers1 = csg1.layers
    layers2 = csg2.layers
    
    # Extract mappings
    mappings1 = csg1.mappings
    mappings2 = csg2.mappings
    
    # Extract neighbor components
    nc_input1 = csg1.neighbor_components[0]
    nc_input2 = csg2.neighbor_components[0]
    
    # Build nc_csg list (pair of neighbor components for each CSG layer)
    nc_csg = []
    for k in range(csg1.L):
        nc_csg.append((
            csg1.neighbor_components[k + 1],
            csg2.neighbor_components[k + 1],
        ))
    
    return (
        layers1, layers2,
        mappings1, mappings2,
        nc_input1, nc_input2,
        nc_csg,
    )


# ============================================================================
# CSG-Transformer integration example
# ============================================================================

def csg_transformer_example():
    """
    Example showing how CSG-Transformer can use this module.
    
    This demonstrates the interface that csg_transformer.py should use.
    """
    import networkx as nx
    
    # Build multi-layer CSG for a pair of graphs
    G1 = nx.cycle_graph(5)
    G2 = nx.path_graph(5)
    
    paired_csg = build_multilayer_csg_pair(G1, G2, L=2)
    
    # Access CSG structure for graph 1
    csg1 = paired_csg.csg1
    
    # Get layers
    for layer_idx in range(csg1.L):
        layer_graph = get_layer_graph(csg1, layer_idx + 1)
        cycle_basis = get_cycle_basis(csg1, layer_idx + 1)
        mapping = get_mapping(csg1, layer_idx + 1)
        neighbor_comp = get_neighbor_components(csg1, layer_idx + 1)
        
        print(f"Layer {layer_idx + 1}: {layer_graph.number_of_nodes()} nodes")
    
    # For CSG-Transformer, you would:
    # 1. Use get_layer_graph() to get each layer's graph structure
    # 2. Use get_neighbor_components() for TNA-Attention computation
    # 3. Use get_mapping() for forward/backward cross-attention
    # 4. Initialize node embeddings using layer_graph.nodes()
    
    return paired_csg
