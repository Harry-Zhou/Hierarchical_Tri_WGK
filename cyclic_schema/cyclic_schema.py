"""
Cyclic Schematic Graph construction and utilities (Step-1).

This module provides functions for:
1. Detecting minimum cycle basis in graphs
2. Constructing cyclic schematic graphs (CSG)
3. Building multi-layer CSG hierarchies with inter-layer mappings
4. Building node correspondences between input graph and CSG (per Step-1)
5. Visualizing the input graph and the resulting CSG

Step-1 Correspondence:
  - For a CSG node of type 'cycle_basis': its corresponding input nodes
    are all input graph nodes that belong to that cycle basis (stored as a tuple).
  - For a CSG node of type 'original_non_cycle': its corresponding input node
    is the node itself (stored as a 1-tuple).
  - For a CSG node of type 'interface': its corresponding input node is the
    node itself (stored as a 1-tuple).

  - For an input graph node v:
    (1) If v appears in multiple cycle bases → maps to a tuple of all
        cycle-basis CSG nodes that contain v.
    (2) If v appears as original_non_cycle in the CSG → maps to a 1-tuple
        containing that original_non_cycle CSG node.
    (3) If v appears as interface in the CSG → maps to a 1-tuple containing
        that interface CSG node.
"""

from collections import defaultdict
from typing import Any, Callable, Dict, Hashable, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.figure import Figure
from matplotlib.patches import Patch


def node_sort_key(node: Any) -> Tuple:
    """
    Sort key for graph nodes to ensure deterministic ordering.

    Multi-layer CSGs mix int nodes (input graph ids) with str nodes
    (cycle_basis names like 'cb0_3'). Group ints by value, then strings
    lexically, so the resulting order is stable across graph constructions.
    """
    if isinstance(node, (int,)):  # Use (int,) for exact type match
        return (0, int(node), "")
    if isinstance(node, str):
        return (1, 0, node)
    # NumPy integers
    if hasattr(node, 'dtype') and 'int' in str(node.dtype):
        return (0, int(node), "")
    return (2, 0, str(node))


# Internal alias for backwards compatibility
_node_sort_key = node_sort_key


def _canonicalize_cycle(cycle: Sequence) -> Tuple:
    """
    Canonicalize a cycle's node ordering for deterministic representation.

    Rotates the cycle to start at the minimum element (by node_sort_key),
    then chooses the direction (forward or reversed) that gives the smaller
    tuple. This ensures isomorphic graphs produce the same cycle signature.
    """
    nodes = list(cycle)
    if not nodes:
        return tuple(nodes)

    # Rotate so that the minimum element comes first
    min_node = min(nodes, key=node_sort_key)
    min_idx = nodes.index(min_node)

    forward = tuple(nodes[min_idx:] + nodes[:min_idx])

    # Reverse the cycle: [v0, v1, ..., v_{n-1}] → [v0, v_{n-1}, ..., v1]
    rev_nodes = [nodes[0]] + nodes[:0:-1]
    rev_min_node = min(rev_nodes, key=node_sort_key)
    rev_min_idx = rev_nodes.index(rev_min_node)
    backward = tuple(rev_nodes[rev_min_idx:] + rev_nodes[:rev_min_idx])

    # Choose the direction with the smaller tuple
    return forward if forward <= backward else backward


def get_cycle_edges(cycle):
    """
    Convert cycle (node list) to edge set (each edge as frozenset, undirected).

    Parameters
    ----------
    cycle : list
        List of nodes in the cycle.

    Returns
    -------
    set
        Set of frozensets, each representing an edge.
    """
    edges = set()
    n = len(cycle)
    for i in range(n):
        u = cycle[i]
        v = cycle[(i + 1) % n]
        edges.add(frozenset([u, v]))
    return edges


def get_edge_in_cycle(cycle):
    """
    Convert cycle to list of (u, v) tuple edges (for edge membership checking).

    Parameters
    ----------
    cycle : list
        List of nodes in the cycle.

    Returns
    -------
    list
        List of (u, v) tuples with u <= v.
    """
    edges = []
    n = len(cycle)
    for i in range(n):
        u = cycle[i]
        v = cycle[(i + 1) % n]
        edges.append((u, v) if str(u) <= str(v) else (v, u))
    return edges


def cyclic_schematic_graph(G, cb_prefix="CB"):
    """
    Construct cyclic schematic graph from input graph G.

    Parameters
    ----------
    G : networkx.Graph
        Input undirected graph.
    cb_prefix : str, optional
        Prefix for cycle basis node names (default: "CB").

    Returns
    -------
    H : networkx.Graph
        Cyclic schematic graph. Node attributes:
          - type: 'cycle_basis' | 'original_non_cycle' | 'interface'
          - If type='cycle_basis': also has cycle_index, original_nodes attributes
    cycle_basis : list
        Minimum cycle basis list (each cycle is a canonicalized node tuple).
    info : dict
        Additional information including node/edge classification details.
    """
    # ============================================================
    # Step 1: Detect minimum cycle basis
    # ============================================================
    raw_basis = nx.minimum_cycle_basis(G)
    # Canonicalize each cycle and sort for deterministic ordering
    cycle_basis = []
    for c in raw_basis:
        canonical = _canonicalize_cycle(c)
        cycle_basis.append(canonical)
    cycle_basis.sort(key=lambda c: (len(c), tuple(node_sort_key(n) for n in c)))

    # Find nodes and edges in any cycle basis
    nodes_in_cycles = set()
    edges_in_cycles = set()  # stored as (min, max) tuples
    for cycle in cycle_basis:
        nodes_in_cycles.update(cycle)
        for u, v in get_edge_in_cycle(cycle):
            edges_in_cycles.add((u, v) if str(u) <= str(v) else (v, u))

    # Nodes not in any cycle basis
    nodes_not_in_cycles = set(G.nodes()) - nodes_in_cycles

    # Edges not in any cycle basis
    edges_not_in_cycles = []
    for u, v in G.edges():
        key = (u, v) if str(u) <= str(v) else (v, u)
        if key not in edges_in_cycles:
            edges_not_in_cycles.append((u, v))

    # ============================================================
    # Step 2: Build initial cyclic schematic graph
    # ============================================================
    H = nx.Graph()

    # Add cycle basis nodes
    num_bases = len(cycle_basis)
    for i in range(num_bases):
        H.add_node(
            f"{cb_prefix}_{i}",
            type="cycle_basis",
            cycle_index=i,
            original_nodes=cycle_basis[i],
        )

    # If two cycle bases share at least one edge, add edge between them
    basis_edge_sets = [get_cycle_edges(c) for c in cycle_basis]
    for i in range(num_bases):
        for j in range(i + 1, num_bases):
            if basis_edge_sets[i] & basis_edge_sets[j]:
                H.add_edge(f"{cb_prefix}_{i}", f"{cb_prefix}_{j}")

    # Add nodes not in any cycle basis
    for node in nodes_not_in_cycles:
        H.add_node(node, type="original_non_cycle")

    # Add edges not in any cycle basis (Type 1: both endpoints not in cycles)
    for u, v in edges_not_in_cycles:
        if u in nodes_not_in_cycles and v in nodes_not_in_cycles:
            H.add_edge(u, v)

    # ============================================================
    # Step 3: Identify interface nodes and build final CSG
    # ============================================================
    # Find connected components of cycle basis nodes (cycle clusters)
    cb_nodes = [f"{cb_prefix}_{i}" for i in range(num_bases)]
    if cb_nodes:
        cb_subgraph = H.subgraph(cb_nodes)
        clusters = list(nx.connected_components(cb_subgraph))
    else:
        clusters = []

    # Track added interface nodes
    interface_nodes_added = set()

    # --- Case 1: Common node between two clusters ---
    for a in range(len(clusters)):
        for b in range(a + 1, len(clusters)):
            # Get all original nodes in cluster a
            nodes_in_cluster_a = set()
            for cb_node in clusters[a]:
                idx = int(cb_node.split("_")[1])
                nodes_in_cluster_a.update(cycle_basis[idx])

            # Get all original nodes in cluster b
            nodes_in_cluster_b = set()
            for cb_node in clusters[b]:
                idx = int(cb_node.split("_")[1])
                nodes_in_cluster_b.update(cycle_basis[idx])

            common_nodes = nodes_in_cluster_a & nodes_in_cluster_b
            if len(common_nodes) == 1:
                v = list(common_nodes)[0]
                # Add interface node to CSG
                if v not in H:
                    H.add_node(v, type="interface")
                    interface_nodes_added.add(v)

                # Connect to cluster a (one edge to first CB containing v)
                for cb_node in clusters[a]:
                    idx = int(cb_node.split("_")[1])
                    if v in cycle_basis[idx]:
                        H.add_edge(v, cb_node)
                        break

                # Connect to cluster b (one edge to first CB containing v)
                for cb_node in clusters[b]:
                    idx = int(cb_node.split("_")[1])
                    if v in cycle_basis[idx]:
                        H.add_edge(v, cb_node)
                        break

    # --- Case 2 & 3: Process edges connecting cycles to non-cycles ---
    for u, v in edges_not_in_cycles:
        u_in_cycle = u in nodes_in_cycles
        v_in_cycle = v in nodes_in_cycles

        if u_in_cycle and not v_in_cycle:
            # u is interface node (in cycle), v is not
            if u not in H:
                H.add_node(u, type="interface")
                interface_nodes_added.add(u)
                # Connect to all cycle bases containing u
                for i, cycle in enumerate(cycle_basis):
                    if u in cycle:
                        H.add_edge(u, f"{cb_prefix}_{i}")
            H.add_edge(u, v)

        elif v_in_cycle and not u_in_cycle:
            # v is interface node (in cycle), u is not
            if v not in H:
                H.add_node(v, type="interface")
                interface_nodes_added.add(v)
                # Connect to all cycle bases containing v
                for i, cycle in enumerate(cycle_basis):
                    if v in cycle:
                        H.add_edge(v, f"{cb_prefix}_{i}")
            H.add_edge(u, v)

        elif u_in_cycle and v_in_cycle:
            # --- Case 3: Both endpoints in different cycle bases ---
            if u not in H:
                H.add_node(u, type="interface")
                interface_nodes_added.add(u)
                # Connect u to all cycle bases containing u
                for i, cycle in enumerate(cycle_basis):
                    if u in cycle:
                        H.add_edge(u, f"{cb_prefix}_{i}")
            if v not in H:
                H.add_node(v, type="interface")
                interface_nodes_added.add(v)
                # Connect v to all cycle bases containing v
                for i, cycle in enumerate(cycle_basis):
                    if v in cycle:
                        H.add_edge(v, f"{cb_prefix}_{i}")
            H.add_edge(u, v)

    # If a node in `nodes_in_cycles` never got added as an interface but IS
    # in some cycle basis yet NOT internal to exactly one cycle cluster,
    # it should still be reachable.  The edge connections above already
    # cover the standard cases — no extra action needed for pure internal
    # nodes (they only participate via their cycle_basis node).

    # ============================================================
    # Compile results
    # ============================================================
    info = {
        "cycle_basis": cycle_basis,
        "nodes_in_cycles": nodes_in_cycles,
        "nodes_not_in_cycles": nodes_not_in_cycles,
        "edges_in_cycles": edges_in_cycles,
        "edges_not_in_cycles": edges_not_in_cycles,
        "clusters": clusters,
        "interface_nodes_added": interface_nodes_added,
    }
    return H, cycle_basis, info


def build_example1_graph():
    """
    Build the test graph G1 from notebook example 1.

    Returns
    -------
    G1 : networkx.Graph
        Test graph with 26 nodes and 41 edges.
    """
    G1 = nx.Graph()
    G1.add_edges_from([
        (0, 2), (0, 3), (0, 4), (0, 18), (0, 19), (0, 25),
        (1, 2), (1, 5), (1, 7), (1, 4),
        (2, 3), (2, 5),
        (3, 4), (3, 7),
        (4, 5), (4, 6),
        (5, 14),
        (6, 7), (6, 14),
        (7, 14),
        (8, 11), (8, 15),
        (9, 11), (9, 12),
        (10, 12), (10, 13),
        (11, 12), (11, 16),
        (12, 13), (12, 16), (12, 17),
        (13, 17),
        (18, 19),
        (20, 21), (20, 25),
        (21, 25),
        (22, 23), (22, 25),
        (23, 24), (23, 25),
        (24, 25)
    ])
    return G1


def build_multilayer_csg(
    G, K=1, cb_prefix_func=None
):
    """
    Build multi-layer cyclic schematic graph hierarchy.

    For K >= 1, produces layers [CSG^1, CSG^2, ..., CSG^K].
    Each layer is a (H, cycle_basis, info) tuple.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    K : int
        Number of CSG layers (default: 1).
    cb_prefix_func : callable, optional
        Function that takes layer index (1-based) and returns prefix.
        Default: lambda i: f"cb{i}"

    Returns
    -------
    layers : list
        List of (H, cycle_basis, info) tuples, one per layer.
    """
    if cb_prefix_func is None:
        cb_prefix_func = lambda i: f"cb{i}"

    layers = []
    current_G = G

    for i in range(1, K + 1):
        H, cycle_basis, info = cyclic_schematic_graph(current_G, cb_prefix=cb_prefix_func(i))
        layers.append((H, cycle_basis, info))
        current_G = H

    return layers


def build_multilayer_csg_with_mappings(G, K=1, cb_prefix_func=None):
    """
    Build multi-layer CSG hierarchy AND the inter-layer mappings.

    For each adjacent pair (layer_{k-1}, layer_k), returns two mappings:
      - csg_to_lower: CSG node → tuple of lower-graph nodes
      - lower_to_csg: lower-graph node → tuple of CSG nodes

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    K : int
        Number of CSG layers (default: 1).
    cb_prefix_func : callable, optional
        Function that takes layer index (1-based) and returns prefix.

    Returns
    -------
    layers : list
        List of (H, cycle_basis, info) tuples, one per layer.
    mappings : list
        List of (csg_to_input, input_to_csg) tuples for each adjacent pair.
        Element k maps between layer k-1 and layer k
        (where layer 0 = input graph G).
    """
    layers = build_multilayer_csg(G, K, cb_prefix_func)

    mappings = []
    current_lower = G

    for k in range(K):
        H_k, cb_k, _ = layers[k]

        csg_to_lower = build_csg_to_input_mapping(H_k, cb_k)
        lower_to_csg = build_input_to_csg_mapping(H_k, cb_k, current_lower)

        mappings.append((csg_to_lower, lower_to_csg))
        current_lower = H_k

    return layers, mappings


def get_node_type(H, node):
    """
    Get the type of a node in the CSG.

    Parameters
    ----------
    H : networkx.Graph
        Cyclic schematic graph.
    node : hashable
        Node identifier.

    Returns
    -------
    str
        Node type: 'cycle_basis', 'original_non_cycle', or 'interface'.
    """
    return H.nodes[node].get('type', 'unknown')


def get_original_nodes_for_csg_node(H, node, cycle_basis):
    """
    Get the original graph nodes corresponding to a CSG node.

    Parameters
    ----------
    H : networkx.Graph
        Cyclic schematic graph.
    node : hashable
        CSG node identifier.
    cycle_basis : list
        List of cycles from the minimum cycle basis.

    Returns
    -------
    tuple
        Tuple of original graph nodes.
    """
    node_type = get_node_type(H, node)

    if node_type == 'cycle_basis':
        idx = H.nodes[node].get('cycle_index', int(str(node).split('_')[1]))
        return tuple(cycle_basis[idx])
    elif node_type in ('original_non_cycle', 'interface'):
        return (node,)
    else:
        return (node,)


def build_csg_to_input_mapping(H, cycle_basis):
    """
    Build mapping from CSG nodes to input graph nodes (Step-1 direction).

    For each CSG node, returns a tuple of the corresponding input graph node(s):
      - cycle_basis: tuple of all input nodes in that cycle
      - original_non_cycle: 1-tuple containing the node itself
      - interface: 1-tuple containing the node itself

    Parameters
    ----------
    H : networkx.Graph
        Cyclic schematic graph.
    cycle_basis : list
        List of cycles from the minimum cycle basis.

    Returns
    -------
    dict
        Mapping from CSG node to tuple of input graph nodes.
    """
    mapping = {}
    for node in H.nodes():
        mapping[node] = get_original_nodes_for_csg_node(H, node, cycle_basis)
    return mapping


def build_input_to_csg_mapping(H, cycle_basis, lower_G):
    """
    Build mapping from lower-graph nodes to CSG nodes (Step-1 reverse direction).

    For each node v in the 'lower' graph (input graph or lower CSG layer),
    returns a tuple of CSG nodes that 'contain' v:
      - If v appears in multiple cycle bases → tuple of all cycle-basis CSG nodes
        whose cycles contain v (sorted deterministically).
      - If v appears as an original_non_cycle node in the CSG → 1-tuple
        containing that original_non_cycle CSG node.
      - If v appears as an interface node in the CSG → 1-tuple containing
        that interface CSG node.
      - Nodes in exactly one cycle basis (internal to a cycle) → 1-tuple
        containing that single cycle-basis CSG node.

    The tuple ordering is deterministic (sorted by node_sort_key) to ensure
    reproducible label propagation.

    Parameters
    ----------
    H : networkx.Graph
        Cyclic schematic graph (the 'higher' layer).
    cycle_basis : list
        List of cycles from the minimum cycle basis used to build H.
    lower_G : networkx.Graph
        The 'lower' graph (input graph or previous CSG layer) whose nodes
        we want to map to CSG nodes.

    Returns
    -------
    dict
        Mapping from lower-graph node to tuple of CSG nodes (deterministic order).
    """
    # Reverse mapping: collect CSG nodes per lower-graph node
    input_to_csg: Dict[Hashable, Any] = {}

    for csg_node in H.nodes():
        # Get the lower-graph nodes that this CSG node corresponds to
        lower_nodes = get_original_nodes_for_csg_node(H, csg_node, cycle_basis)
        for lower_node in lower_nodes:
            if lower_node not in input_to_csg:
                input_to_csg[lower_node] = []
            input_to_csg[lower_node].append(csg_node)

    # Convert lists to sorted tuples for deterministic ordering
    for node in list(input_to_csg.keys()):
        input_to_csg[node] = tuple(sorted(input_to_csg[node], key=node_sort_key))

    # Ensure ALL lower-graph nodes are in the mapping
    for node in lower_G.nodes():
        if node not in input_to_csg:
            input_to_csg[node] = ()

    return input_to_csg


def draw_input_graph(G, ax=None, title="Input Graph", layout_seed=42):
    """
    Draw the input graph.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created.
    title : str
        Plot title.
    layout_seed : int
        Seed for spring_layout (ensures reproducible layouts).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    pos = nx.spring_layout(G, seed=layout_seed)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color="lightblue",
                           node_size=500, edgecolors="black")
    nx.draw_networkx_edges(G, pos, ax=ax, width=1.5, alpha=0.7)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_weight="bold")

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    return fig


def draw_cyclic_schematic(H, cycle_basis, info, ax=None,
                          title="Cyclic Schematic Graph", layout_seed=42):
    """
    Draw the cyclic schematic graph with node-type-aware styling.

    Different node types use distinct colors and shapes:
      - cycle_basis: blue circle
      - original_non_cycle: green square
      - interface: red diamond

    Parameters
    ----------
    H : networkx.Graph
        Cyclic schematic graph.
    cycle_basis : list
        List of cycles from the minimum cycle basis.
    info : dict
        Information dictionary returned by cyclic_schematic_graph.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created.
    title : str
        Plot title.
    layout_seed : int
        Seed for spring_layout.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    else:
        fig = ax.figure

    cb_nodes = [n for n, d in H.nodes(data=True) if d.get("type") == "cycle_basis"]
    orig_non_cycle = [n for n, d in H.nodes(data=True)
                      if d.get("type") == "original_non_cycle"]
    interface_nodes = [n for n, d in H.nodes(data=True)
                       if d.get("type") == "interface"]

    pos = nx.spring_layout(H, seed=layout_seed, k=2, iterations=50)

    nx.draw_networkx_edges(H, pos, ax=ax, width=1.2, alpha=0.5, edge_color="gray")

    if cb_nodes:
        nx.draw_networkx_nodes(
            H, pos, nodelist=cb_nodes, ax=ax,
            node_color="#4A90D9", node_size=800,
            edgecolors="#2C5F8A", linewidths=2,
        )

    if orig_non_cycle:
        nx.draw_networkx_nodes(
            H, pos, nodelist=orig_non_cycle, ax=ax,
            node_color="#7DD67D", node_size=600,
            edgecolors="#2E8B2E", linewidths=2,
            node_shape="s",
        )

    if interface_nodes:
        nx.draw_networkx_nodes(
            H, pos, nodelist=interface_nodes, ax=ax,
            node_color="#FF6B6B", node_size=700,
            edgecolors="#CC3333", linewidths=2,
            node_shape="d",
        )

    labels = {}
    for n in H.nodes():
        d = H.nodes[n]
        labels[n] = str(n)

    nx.draw_networkx_labels(H, pos, labels=labels, ax=ax,
                            font_size=9, font_weight="bold")

    legend_elements = [
        Patch(facecolor="#4A90D9", edgecolor="#2C5F8A", label="Cycle Basis (CB)"),
        Patch(facecolor="#7DD67D", edgecolor="#2E8B2E",
              label="Original Node (not in cycles)"),
        Patch(facecolor="#FF6B6B", edgecolor="#CC3333", label="Interface Node"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=10)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    return fig


def draw_input_graph_with_cycles(G, cycle_basis, ax=None,
                                 title="Input Graph (cycle-highlighted)",
                                 layout_seed=42):
    """
    Draw the input graph with cycle-basis highlighting.

    Each cycle basis uses an independent color; nodes belonging to multiple
    cycle bases are highlighted in gold with a red border.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    cycle_basis : list
        List of cycles from the minimum cycle basis.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on. If None, a new figure is created.
    title : str
        Plot title.
    layout_seed : int
        Seed for spring_layout.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing the plot.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    pos = nx.spring_layout(G, seed=layout_seed)

    node_cycles = defaultdict(list)
    for i, cycle in enumerate(cycle_basis):
        for node in cycle:
            node_cycles[node].append(i)

    cycle_cmap = plt.get_cmap("tab10")
    cycle_colors = [cycle_cmap(i) for i in range(cycle_cmap.N)]

    nx.draw_networkx_edges(G, pos, ax=ax, width=1.5, alpha=0.7, edge_color="gray")

    non_cycle_nodes = [n for n in G.nodes() if n not in node_cycles]
    if non_cycle_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=non_cycle_nodes, ax=ax,
                               node_color="lightblue", node_size=500,
                               edgecolors="black")

    single_cycle_nodes_by_color = defaultdict(list)
    multi_cycle_nodes = []
    for node, cycles in node_cycles.items():
        if len(cycles) == 1:
            single_cycle_nodes_by_color[cycles[0]].append(node)
        else:
            multi_cycle_nodes.append(node)

    for cycle_idx, nodes in single_cycle_nodes_by_color.items():
        color = cycle_colors[cycle_idx % len(cycle_colors)]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax,
                               node_color=[color], node_size=500,
                               edgecolors="black")

    if multi_cycle_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=multi_cycle_nodes, ax=ax,
                               node_color="gold", node_size=550,
                               edgecolors="red", linewidths=2, node_shape="d")

    nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_weight="bold")

    legend_elements = [
        Patch(facecolor="#4A90D9", edgecolor="#2C5F8A", label="Cycle Basis (CB)"),
        Patch(facecolor="#7DD67D", edgecolor="#2E8B2E",
              label="Original Node (not in cycles)"),
        Patch(facecolor="#FF6B6B", edgecolor="#CC3333", label="Interface Node"),
    ]

    ax.legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.8)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    return fig


def draw_side_by_side(G, H, cycle_basis, info, layout_seed=42):
    """
    Draw the input graph and the cyclic schematic graph side by side.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    H : networkx.Graph
        Cyclic schematic graph.
    cycle_basis : list
        List of cycles from the minimum cycle basis.
    info : dict
        Information dictionary returned by cyclic_schematic_graph.
    layout_seed : int
        Seed for spring_layout.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing both subplots.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))

    draw_input_graph(G, ax=ax1, title="(a) Input Graph", layout_seed=layout_seed)
    draw_cyclic_schematic(H, cycle_basis, info, ax=ax2,
                          title="(b) Cyclic Schematic Graph", layout_seed=layout_seed)

    plt.tight_layout()
    return fig


def draw_side_by_side_with_cycles(G, H, cycle_basis, info, layout_seed=42):
    """
    Draw the input graph (with cycle highlighting) and the CSG side by side.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    H : networkx.Graph
        Cyclic schematic graph.
    cycle_basis : list
        List of cycles from the minimum cycle basis.
    info : dict
        Information dictionary returned by cyclic_schematic_graph.
    layout_seed : int
        Seed for spring_layout.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure containing both subplots.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, 9))

    draw_input_graph_with_cycles(G, cycle_basis, ax=ax1,
                                 title="(a) Input Graph (cycle-highlighted)",
                                 layout_seed=layout_seed)
    draw_cyclic_schematic(H, cycle_basis, info, ax=ax2,
                          title="(b) Cyclic Schematic Graph",
                          layout_seed=layout_seed)

    plt.tight_layout()
    return fig


def print_analysis_info(G, H, cycle_basis, info):
    """
    Print a human-readable analysis report.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    H : networkx.Graph
        Cyclic schematic graph.
    cycle_basis : list
        List of cycles from the minimum cycle basis.
    info : dict
        Information dictionary returned by cyclic_schematic_graph.
    """
    print("=" * 60)
    print("Graph Analysis Report")
    print("=" * 60)

    print(f"\nInput graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    print(f"\n--- Step 1: Minimum Cycle Basis ---")
    print(f"Detected {len(cycle_basis)} minimum cycle basis cycles:")
    cb_nodes_h = [n for n, d in H.nodes(data=True) if d.get("type") == "cycle_basis"]
    cb_prefix = cb_nodes_h[0].rsplit("_", 1)[0] if cb_nodes_h else "CB"
    for i, cycle in enumerate(cycle_basis):
        print(f"  {cb_prefix}_{i}: {cycle}")
        ce = get_edge_in_cycle(cycle)
        print(f"         Edges: {ce}")

    print(f"\nNodes in at least one cycle basis: "
          f"{sorted(info['nodes_in_cycles'], key=str)}")
    print(f"Nodes not in any cycle basis: "
          f"{sorted(info['nodes_not_in_cycles'], key=str)}")
    print(f"Edges not in any cycle basis:")
    for e in info['edges_not_in_cycles']:
        print(f"  {e}")

    print(f"\n--- Step 2: Initial Cyclic Schematic Graph ---")
    print(f"H nodes: {H.number_of_nodes()}, edges: {H.number_of_edges()}")

    print(f"\n--- Step 3: Cycle Clusters and Interface Nodes ---")
    clusters = info['clusters']
    print(f"There are {len(clusters)} cycle cluster(s) (CB-only connected components):")
    for i, cluster in enumerate(clusters):
        nodes_in_cluster = set()
        for cb_node in cluster:
            idx = int(cb_node.split("_")[1])
            nodes_in_cluster.update(cycle_basis[idx])
        print(f"  Cluster {i}: CB nodes={sorted(cluster, key=str)}, "
              f"original nodes={sorted(nodes_in_cluster, key=str)}")

    if info['interface_nodes_added']:
        print(f"Interface nodes: "
              f"{sorted(info['interface_nodes_added'], key=str)}")
    else:
        print(f"Interface nodes: none")

    print(f"\nFinal Cyclic Schematic Graph:")
    print(f"  H nodes: {H.number_of_nodes()}, edges: {H.number_of_edges()}")
    print(f"  Nodes: {list(H.nodes())}")
    print(f"  Edges: {list(H.edges())}")
    print("=" * 60)


if __name__ == "__main__":
    # Example usage
    G1 = build_example1_graph()
    print(f"Example 1: {G1.number_of_nodes()} nodes, {G1.number_of_edges()} edges")

    # Build CSG
    H, cycle_basis, info = cyclic_schematic_graph(G1, cb_prefix="cb0")
    print(f"CSG: {H.number_of_nodes()} nodes, {H.number_of_edges()} edges")
    print(f"Cycle basis: {len(cycle_basis)} cycles")

    # Build mappings
    csg_to_input = build_csg_to_input_mapping(H, cycle_basis)
    input_to_csg = build_input_to_csg_mapping(H, cycle_basis, G1)

    print(f"\nCSG to input mapping:")
    for csg_node, input_nodes in sorted(csg_to_input.items(), key=lambda x: str(x[0])):
        print(f"  {csg_node}: {input_nodes}")

    print(f"\nInput to CSG mapping (first 5):")
    for input_node, csg_nodes in sorted(input_to_csg.items(), key=lambda x: str(x[0]))[:5]:
        print(f"  {input_node}: {csg_nodes}")

    # Demonstration of visualization (only when not in a headless context)
    try:
        import matplotlib
        if matplotlib.get_backend().lower() != "agg":
            print_analysis_info(G1, H, cycle_basis, info)
            draw_side_by_side_with_cycles(G1, H, cycle_basis, info)
            plt.show()
    except (ImportError, AttributeError):
        pass
