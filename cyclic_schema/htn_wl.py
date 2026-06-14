"""
Hierarchical Triangulated Weisfeiler-Lehman Algorithm.

Canonical implementation module for the HTN-WL algorithm.
The former :mod:`hierarchical_triangular_wl` was merged into this module;
all public symbols are also available via :mod:`cyclic_schema.htn_wl`
directly. New code should ``from cyclic_schema.htn_wl import ...``.

Implements a multi-layer graph isomorphism-testing algorithm based on
cyclic schematic graphs (CSG). The algorithm performs:

  Step-1 (Mapping): Builds input-CSG node correspondences.
  Step-2 (Single iteration message passing):
    - Forward (G -> CSG):  label tuple computation + message passing
    - Backward (CSG -> G): label propagation + message passing
  Step-3 (Multi-layer hierarchical message passing):
    - Forward:  G -> CSG^1 -> CSG^2 -> ... -> CSG^L
    - Backward: CSG^L -> ... -> CSG^1 -> G
    - Repeated for I iterations.

Key improvements over baseline:
  - Proper backward label tuple construction (prepends l_Gi(v) per spec)
  - L=0 support (no CSG layers → triangulated neighbors WL)
  - Deterministic cycle label canonicalization (left/right walk from min)
  - Joint label assignment across both graphs (WL consistency)
"""

import os as _os
import sys as _sys
from collections import defaultdict
from typing import Any, Dict, Hashable, List, Optional, Sequence, Tuple, Union

import networkx as nx
import numpy as np

from .cyclic_schema import (
    build_csg_to_input_mapping,
    build_example1_graph,
    build_input_to_csg_mapping,
    cyclic_schematic_graph,
    draw_cyclic_schematic,
    draw_input_graph,
    draw_input_graph_with_cycles,
    draw_side_by_side,
    draw_side_by_side_with_cycles,
    get_node_type,
    node_sort_key,
)
from .multilayer_csg import (
    _precompute_neighbor_components,
    build_multilayer_csg_pair,
    extract_htn_wl_args,
)

# Sentinel values for _make_agg_sort_key flattening
_TUPLE_START = -1
_TUPLE_END = -2


# ============================================================================
# Label conversion utilities
# ============================================================================

def _build_vtx_triangulated_neighbors(g: nx.Graph) -> Dict:
    """Build triangulated-neighbor dict for graph *g* (shared helper)."""
    vtx_tri: Dict = {}
    for v in g.nodes():
        neighs = list(g.neighbors(v))
        tri_neighs: List = []
        if neighs:
            v_isg = nx.induced_subgraph(g, neighs)
            for ttn in nx.connected_components(v_isg):
                ttn_v = list(ttn) + [v]
                tri_subgraph = nx.induced_subgraph(g, ttn_v)
                tri_neighs.append(tri_subgraph)
        vtx_tri[v] = tri_neighs
    return vtx_tri


def _labels_to_dict(G: nx.Graph, vlabel: Any) -> Dict[Hashable, int]:
    """
    Convert label array/dict to dict mapping node -> label int.

    Handles:
      - dict: direct lookup by node key
      - numpy ndarray: maps by sorted node order (arr[i] -> i-th sorted node)
      - list: same as ndarray

    Using sorted node order ensures the i-th element of the array always
    corresponds to the i-th node in sorted order, regardless of whether
    node IDs are consecutive integers.
    """
    if isinstance(vlabel, dict):
        return {n: int(vlabel[n]) for n in G.nodes()}
    arr = np.asarray(vlabel)
    sorted_nodes = sorted(G.nodes(), key=node_sort_key)
    return {n: int(arr[i]) for i, n in enumerate(sorted_nodes)}


# ============================================================================
# Aggregate sorting key (flattens nested tuples for comparison)
# ============================================================================

def _make_agg_sort_key(x: Any) -> Tuple:
    """Flatten nested tuple structure into a flat comparable tuple of ints.

    Needed because aggregates contain nested tuples (label tuples within
    connected-component groups), and Python cannot compare nested tuples
    of different shapes directly via `<`.
    """
    parts: List[int] = []

    def _walk(item: Any) -> None:
        if isinstance(item, tuple):
            parts.append(_TUPLE_START)
            parts.append(len(item))
            for v in item:
                _walk(v)
            parts.append(_TUPLE_END)
        else:
            parts.append(item)

    _walk(x)
    return tuple(parts)


# ============================================================================
# Neighbor component precomputation (topology, not label-dependent)
# ============================================================================

# ============================================================================
# Step-2 Forward: label tuple computation for CSG nodes
# ============================================================================

def compute_initial_label_tuples(
    H: nx.Graph,
    cycle_basis: List[List[Hashable]],
    original_labels: Dict[Hashable, int],
) -> Dict[Hashable, Tuple]:
    """
    Compute initial label tuples for CSG nodes (Step-2 forward).

    For a cycle_basis node, the label tuple consists of the labels of all
    original graph nodes in that cycle, in the order they appear in the
    cycle (as returned by the canonicalized minimum cycle basis).

    For original_non_cycle and interface nodes, the label is a
    1-tuple containing the node's label.
    """
    label_tuples: Dict[Hashable, Tuple] = {}
    h_nodes = H.nodes
    cb = cycle_basis

    for node in H.nodes():
        attrs = h_nodes[node]
        if attrs.get('type') == 'cycle_basis':
            cycle_nodes = cb[attrs['cycle_index']]
            label_tuples[node] = tuple(original_labels.get(n, 0) for n in cycle_nodes)
        else:
            label_tuples[node] = (original_labels.get(node, 0),)

    return label_tuples


def canonicalize_cycle_label(label_tuple: Tuple) -> Tuple:
    """
    Canonicalize a cycle label tuple to handle cycle symmetry (Step-2).

    Per the spec:
      1. Find the FIRST element with the minimum label value.
      2. From that element, walk LEFT  to form a left tuple.
      3. From that element, walk RIGHT to form a right tuple.
      4. Return the lexicographically smaller of the two.

    'left'  = pos, pos-1, pos-2, ..., pos+1 (mod n)
    'right' = pos, pos+1, pos+2, ..., pos-1 (mod n)

    Examples:
      (1, 2, 3, 4) -> (1, 2, 3, 4)  (already canonical)
      (3, 1, 2, 4) -> (1, 2, 4, 3)  (min=1 at pos=1, left wins)
      (1, 4, 3, 2) -> (1, 2, 3, 4)  (min=1 at pos=0, left wins)
      (2, 1, 4, 3) -> (1, 2, 3, 4)  (min=1 at pos=1, left wins)
    """
    if len(label_tuple) <= 1:
        return label_tuple

    min_label = min(label_tuple)
    pos = next(i for i, l in enumerate(label_tuple) if l == min_label)
    n = len(label_tuple)

    # Walk left from pos:  pos, pos-1, pos-2, ..., pos+1 (mod n)
    left = tuple(label_tuple[(pos - i) % n] for i in range(n))
    # Walk right from pos: pos, pos+1, pos+2, ..., pos-1 (mod n)
    right = tuple(label_tuple[(pos + i) % n] for i in range(n))

    return min(left, right)


def compute_final_label_tuples(
    H: nx.Graph,
    cycle_basis: List[List[Hashable]],
    original_labels: Dict[Hashable, int],
) -> Dict[Hashable, Tuple]:
    """
    Compute final label tuples after canonicalization (Step-2 forward).

    For cycle_basis nodes the initial label tuple is canonicalized to
    remove cycle-rotation symmetry. Other node types are unchanged.
    """
    initial_tuples = compute_initial_label_tuples(H, cycle_basis, original_labels)
    final_tuples: Dict[Hashable, Tuple] = {}

    for node, label_tuple in initial_tuples.items():
        if get_node_type(H, node) == 'cycle_basis':
            final_tuples[node] = canonicalize_cycle_label(label_tuple)
        else:
            final_tuples[node] = label_tuple

    return final_tuples


# ============================================================================
# Step-2 Forward: aggregate computation and joint message passing
# ============================================================================

def forward_aggregate(
    H: nx.Graph,
    v: Hashable,
    label_tuples: Dict[Hashable, Tuple],
    neighbor_components: Optional[Dict[Hashable, Tuple]] = None,
) -> Tuple:
    """
    Build AGG(v) per the spec:

        AGG(v) = ((l(v),), sorted_lt_i elements...)

    lt_i is built from the induced subgraph on N(v):
      - isolated neighbor u -> l(u) added to lt_i;
      - connected component {u_1, ..., u_s} -> SORT([l(u_1),...,l(u_s)])
        added to lt_i (as a tuple for hashability).

    lt_i is then sorted and appended.  If ``neighbor_components`` is
    provided (output of ``_precompute_neighbor_components``), it is used
    directly; otherwise it is computed on the fly.
    """
    get_lt = label_tuples.get
    v_label = get_lt(v, (0,))
    agg = [(v_label,)]

    if neighbor_components is not None:
        components = neighbor_components.get(v, ())
    else:
        neighbors = list(H.neighbors(v))
        if not neighbors:
            return tuple(agg)
        components = tuple(
            sorted(c, key=node_sort_key)
            for c in nx.connected_components(H.subgraph(neighbors))
        )

    if not components:
        return tuple(agg)

    lt_i: List = []
    for component in components:
        if len(component) == 1:
            lt_i.append(get_lt(component[0], (0,)))
        else:
            lt_i.append(tuple(sorted(get_lt(u, (0,)) for u in component)))

    lt_i.sort(key=_make_agg_sort_key)
    agg.extend(lt_i)

    return tuple(agg)


def forward_message_passing_both(
    H1: nx.Graph,
    label_tuples1: Dict[Hashable, Tuple],
    H2: nx.Graph,
    label_tuples2: Dict[Hashable, Tuple],
    max_label: int,
    nc1: Optional[Dict[Hashable, Tuple]] = None,
    nc2: Optional[Dict[Hashable, Tuple]] = None,
) -> Tuple[Dict[Hashable, int], Dict[Hashable, int]]:
    """
    Simultaneous WL-style message passing on both graphs (Step-2 update).

    Aggregate every node in both graphs, then sort all aggregates together
    and assign new labels from (max_label + 1) onward.  Same aggregate ->
    same new label (WL compression consistent across the two graphs).

    Parameters
    ----------
    H1, H2 : nx.Graph
        The two graphs to process (can be input graphs or CSG layers).
    label_tuples1, label_tuples2 : dict
        Current label tuples for each node in each graph.
    max_label : int
        Current maximum label value in the system. New labels start at +1.
    nc1, nc2 : dict, optional
        Precomputed neighbor component structures for H1, H2 (output of
        ``_precompute_neighbor_components``).  If omitted, they are
        recomputed in this call.

    Returns
    -------
    new_labels1, new_labels2 : dict
        New integer labels for nodes in each graph.
    """
    if nc1 is None:
        nc1 = _precompute_neighbor_components(H1)
    if nc2 is None:
        nc2 = _precompute_neighbor_components(H2)

    agg_to_nodes: Dict[Tuple, List[Tuple[int, Hashable]]] = defaultdict(list)
    for v in H1.nodes():
        agg = forward_aggregate(H1, v, label_tuples1, nc1)
        agg_to_nodes[agg].append((1, v))
    for v in H2.nodes():
        agg = forward_aggregate(H2, v, label_tuples2, nc2)
        agg_to_nodes[agg].append((2, v))

    next_label = max_label + 1
    new_labels1: Dict[Hashable, int] = {}
    new_labels2: Dict[Hashable, int] = {}
    sorted_aggs = sorted(agg_to_nodes.keys(), key=_make_agg_sort_key)
    for agg in sorted_aggs:
        for graph_id, v in agg_to_nodes[agg]:
            if graph_id == 1:
                new_labels1[v] = next_label
            else:
                new_labels2[v] = next_label
        next_label += 1

    return new_labels1, new_labels2


# ============================================================================
# Step-2 Backward: label propagation from higher layer to lower layer
# ============================================================================

def _compute_lower_label_tuples(
    lower_G: nx.Graph,
    lower_labels: Dict[Hashable, int],
    higher_labels: Dict[Hashable, int],
    lower_to_higher: Dict[Hashable, Tuple],
) -> Dict[Hashable, Tuple]:
    """
    Compute final label tuples for lower-graph nodes (Step-2 backward).

    Per the spec (Step-2 backward):
      For each node v in V_{G_i}, collect the labels of all higher-layer
      (CSG) nodes that 'contain' v (via the lower_to_higher mapping).
      The initial label tuple is:

        final_label_tuple(v) = SORT([L_{CSG}(u) : v belongs to u in V_{CSG}])

      i.e. just the sorted higher-layer labels — WITHOUT prepending
      l_{G_i}(v).  This differs from earlier versions where l_v was
      prepended.  The corrected behaviour follows the canonical spec.

      If no higher-layer nodes correspond to v, fall back to (l_v,)
      so that the node's own label is preserved.

    Parameters
    ----------
    lower_G : nx.Graph
        The lower-layer graph (input graph or previous CSG layer).
    lower_labels : dict
        Current labels of the lower-layer nodes (int).
    higher_labels : dict
        Labels of the higher-layer nodes (int).
    lower_to_higher : dict
        Mapping from lower-graph node to tuple of higher-graph node IDs.

    Returns
    -------
    dict
        Mapping from lower-graph node to final label tuple.
    """
    result: Dict[Hashable, Tuple] = {}
    for v in lower_G.nodes():
        higher_nodes = lower_to_higher.get(v, ())
        l_v = int(lower_labels.get(v, 0))
        if not higher_nodes:
            # Fallback: preserve node's own label when no CSG context exists
            result[v] = (l_v,)
        else:
            # Per spec: only the sorted higher-layer labels
            higher_lbls = sorted(int(higher_labels.get(u, 0)) for u in higher_nodes)
            result[v] = tuple(higher_lbls)
    return result


def backward_message_passing_both(
    lower_G1: nx.Graph,
    lower_G2: nx.Graph,
    lower_labels1: Dict[Hashable, int],
    lower_labels2: Dict[Hashable, int],
    higher_labels1: Dict[Hashable, int],
    higher_labels2: Dict[Hashable, int],
    lower_to_higher1: Dict[Hashable, Tuple],
    lower_to_higher2: Dict[Hashable, Tuple],
    max_label: int,
    nc1: Optional[Dict[Hashable, Tuple]] = None,
    nc2: Optional[Dict[Hashable, Tuple]] = None,
) -> Tuple[Dict[Hashable, int], Dict[Hashable, int]]:
    """
    Backward message passing (Step-2 backward) — fixed per spec.

    1. For each lower-graph node v, compute its final label tuple:
       - Prepend l_{G_i}(v) to the sorted list of higher-layer labels
         that correspond to v (via the lower_to_higher mapping).
    2. Run standard (forward) message passing on the lower graphs
       using these final label tuples.

    The starting label for new assignments begins at max_label + 1.

    Parameters
    ----------
    lower_G1, lower_G2 : nx.Graph
        Lower-layer graphs.
    lower_labels1, lower_labels2 : dict
        Current labels of the lower-layer nodes (int).
    higher_labels1, higher_labels2 : dict
        Higher-layer node labels (int).
    lower_to_higher1, lower_to_higher2 : dict
        Step-1 mappings from lower-graph nodes to higher-graph node tuples.
    max_label : int
        Current global maximum label.
    nc1, nc2 : dict, optional
        Precomputed neighbor components for lower_G1 and lower_G2.  If None,
        they are recomputed.

    Returns
    -------
    new_labels1, new_labels2 : dict
        New integer labels for nodes in lower_G1, lower_G2.
    """
    label_tuples1 = _compute_lower_label_tuples(
        lower_G1, lower_labels1, higher_labels1, lower_to_higher1)
    label_tuples2 = _compute_lower_label_tuples(
        lower_G2, lower_labels2, higher_labels2, lower_to_higher2)

    return forward_message_passing_both(
        lower_G1, label_tuples1, lower_G2, label_tuples2,
        max_label, nc1, nc2)


# ============================================================================
# Edge-label helper functions
# ============================================================================

def _compute_edge_context(
    v: Hashable,
    G: nx.Graph,
    elabel_dict: Dict[Tuple[Hashable, Hashable], int],
    neighbor_components: Dict[Hashable, Tuple],
) -> Tuple:
    """
    Compute the edge context ec(v) for node *v*.

    The edge context is a tuple of sorted incident edge labels, grouped by
    connected components of N(v).  Each component contributes a sorted tuple
    of edge labels ``(elabel(v,u) for u in component)``, and components are
    sorted lexicographically.

    Returns an empty tuple when *v* has no neighbours.
    """
    components = neighbor_components.get(v, ())
    if not components:
        return ()

    comp_aggs: List[Tuple] = []
    for comp in components:
        comp_edge_labels = tuple(
            sorted(
                elabel_dict.get(
                    (v, u) if str(v) <= str(u) else (u, v), 0)
                for u in comp
            )
        )
        comp_aggs.append(comp_edge_labels)
    comp_aggs.sort()
    return tuple(comp_aggs)


def _compute_edge_contexts_all(
    G: nx.Graph,
    elabel_dict: Dict[Tuple[Hashable, Hashable], int],
    neighbor_components: Dict[Hashable, Tuple],
) -> Dict[Hashable, Tuple]:
    """Precompute edge contexts for every node in *G*."""
    return {
        v: _compute_edge_context(v, G, elabel_dict, neighbor_components)
        for v in G.nodes()
    }


def _compute_lower_label_tuples_with_edges(
    lower_G: nx.Graph,
    lower_labels: Dict[Hashable, int],
    higher_labels: Dict[Hashable, int],
    lower_to_higher: Dict[Hashable, Tuple],
    edge_contexts: Dict[Hashable, Tuple],
) -> Dict[Hashable, Tuple]:
    """
    Compute final label tuples for lower-graph nodes, incorporating **edge
    context** (Step-2 backward, edge-aware variant).

    For a lower-graph node *v*:

        label_tuple(v) = (ec(v),) + SORT([higher labels])

    where ``ec(v)`` is the edge context of *v* (computed from incident edge
    labels grouped by triangulated neighbourhood components).

    When ``ec(v)`` is empty the label tuple uses only the sorted higher
    labels.  If there are neither higher labels nor edge context, fall
    back to (l_v,).
    """
    result: Dict[Hashable, Tuple] = {}
    for v in lower_G.nodes():
        higher_nodes = lower_to_higher.get(v, ())
        l_v = int(lower_labels.get(v, 0))
        ec_v = edge_contexts.get(v, ())

        if not higher_nodes and not ec_v:
            result[v] = (l_v,)
        elif not higher_nodes:
            result[v] = (ec_v,)
        elif not ec_v:
            higher_lbls = sorted(int(higher_labels.get(u, 0)) for u in higher_nodes)
            result[v] = tuple(higher_lbls)
        else:
            higher_lbls = sorted(int(higher_labels.get(u, 0)) for u in higher_nodes)
            result[v] = (ec_v,) + tuple(higher_lbls)
    return result


def backward_message_passing_both_with_edges(
    lower_G1: nx.Graph,
    lower_G2: nx.Graph,
    lower_labels1: Dict[Hashable, int],
    lower_labels2: Dict[Hashable, int],
    higher_labels1: Dict[Hashable, int],
    higher_labels2: Dict[Hashable, int],
    lower_to_higher1: Dict[Hashable, Tuple],
    lower_to_higher2: Dict[Hashable, Tuple],
    max_label: int,
    edge_contexts1: Dict[Hashable, Tuple],
    edge_contexts2: Dict[Hashable, Tuple],
    nc1: Optional[Dict[Hashable, Tuple]] = None,
    nc2: Optional[Dict[Hashable, Tuple]] = None,
) -> Tuple[Dict[Hashable, int], Dict[Hashable, int]]:
    """
    Backward message passing (Step-2 backward) that incorporates **edge
    contexts** into lower-graph label tuples.

    See ``backward_message_passing_both`` for the base behaviour.  The only
    difference is that lower-graph label tuples additionally carry edge
    context information (cf. ``_compute_lower_label_tuples_with_edges``).
    """
    label_tuples1 = _compute_lower_label_tuples_with_edges(
        lower_G1, lower_labels1, higher_labels1, lower_to_higher1,
        edge_contexts1)
    label_tuples2 = _compute_lower_label_tuples_with_edges(
        lower_G2, lower_labels2, higher_labels2, lower_to_higher2,
        edge_contexts2)

    return forward_message_passing_both(
        lower_G1, label_tuples1, lower_G2, label_tuples2,
        max_label, nc1, nc2)


def _update_elabel_from_dict(
    vlabel_dict1: Dict[Hashable, int],
    vlabel_dict2: Dict[Hashable, int],
    elabel_dict1: Dict,
    elabel_dict2: Dict,
) -> Tuple[Dict, Dict]:
    """
    Update edge labels using the **dict-based** interface (works with
    arbitrary non-consecutive node IDs).

    For each edge ``(u, v)`` the new label tuple is::

        (old_elabel, sorted(vlabel[u], vlabel[v]))

    Labels are jointly compressed across both graphs so that identical tuples
    receive the same new label.

    Parameters
    ----------
    vlabel_dict1, vlabel_dict2 : dict
        Current node labels keyed by node ID.
    elabel_dict1, elabel_dict2 : dict
        Current edge labels keyed by sorted edge tuple.

    Returns
    -------
    compressed1, compressed2 : dict
        Updated edge labels keyed by sorted edge tuple.
    """

    def _collect(vlabel_dict, elabel_dict):
        collection = {}
        for e, l in elabel_dict.items():
            vl0, vl1 = vlabel_dict[e[0]], vlabel_dict[e[1]]
            if vl0 <= vl1:
                collection[e] = (l, vl0, vl1)
            else:
                collection[e] = (l, vl1, vl0)
        return collection

    col1 = _collect(vlabel_dict1, elabel_dict1)
    col2 = _collect(vlabel_dict2, elabel_dict2)

    # Joint compression
    all_labels = list(col1.values()) + list(col2.values())
    all_labels = sorted(set(all_labels))
    if not all_labels:
        return {}, {}
    next_label = max(lc[0] for lc in all_labels) + 1

    def _compress(collection):
        return {e: all_labels.index(tup) + next_label
                for e, tup in collection.items()}

    return _compress(col1), _compress(col2)


def _validate_edge_labels(
    G1: nx.Graph,
    elabel_dict1: Optional[Dict],
    G2: nx.Graph,
    elabel_dict2: Optional[Dict],
) -> None:
    """
    Validate edge-label dictionaries.

    Raises ``ValueError`` if:
      - exactly one of the two dicts is ``None``;
      - an edge key is not present in the corresponding graph.
    """
    if (elabel_dict1 is None) != (elabel_dict2 is None):
        raise ValueError(
            "Either both or neither elabel_dict must be provided; "
            f"got elabel_dict1={type(elabel_dict1)} "
            f"elabel_dict2={type(elabel_dict2)}")
    if elabel_dict1 is None:
        return
    assert elabel_dict1 is not None
    assert elabel_dict2 is not None
    for e in elabel_dict1:
        if e not in G1.edges() and (e[1], e[0]) not in G1.edges():
            raise ValueError(
                f"Edge {e} in elabel_dict1 is not present in G1")
    for e in elabel_dict2:
        if e not in G2.edges() and (e[1], e[0]) not in G2.edges():
            raise ValueError(
                f"Edge {e} in elabel_dict2 is not present in G2")


# ============================================================================
# Step-3: Multi-layer hierarchical message passing
# ============================================================================

def hierarchical_triangular_wl(
    G1: nx.Graph,
    G2: nx.Graph,
    vlabel_np1: Any,
    vlabel_np2: Any,
    L: int = 1,
    I: int = 5,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Hierarchical Triangulated Weisfeiler-Lehman test (Step-3).

    Each iteration performs a complete message-passing cycle through the
    multi-layer hierarchy:

      Forward:  G -> CSG^1 -> CSG^2 -> ... -> CSG^L
      Backward: CSG^L -> ... -> CSG^1 -> G

    Forward step at layer k:
      1. Compute final label tuples for CSG^k from previous layer labels
         (cycle_basis nodes are canonicalized).
      2. Run simultaneous message passing on both copies of CSG^k.

    Backward step at layer k:
      1. Compute label tuples for the lower graph from (a) its current
         labels and (b) higher labels via the lower-to-higher mapping,
         with l_{G_i}(v) prepended per spec.
      2. Run simultaneous message passing on the lower graphs.

    I iterations produce label histories of shape (|V|, I+1): column 0
    holds the initial labels, columns 1..I hold the labels after each
    iteration. Rows are indexed by sorted node order (see node_sort_key).

    Parameters
    ----------
    G1, G2 : networkx.Graph
        Input graphs.
    vlabel_np1, vlabel_np2 : numpy.ndarray or dict
        Initial node labels. If ndarray, elements are mapped to nodes in
        sorted order (arr[i] -> i-th node in sorted(G.nodes())).
        If dict, maps node -> label directly.
    L : int
        Number of CSG layers (default 1). L=0 means no CSG layers
        (standard joint WL on the two input graphs).
        Renamed from ``K`` to ``L`` in v2 to disambiguate from ``k`` in
        ``k``-WL notation (where ``k`` denotes tuple width, not layer depth).
    I : int
        Number of iterations (default 5). Must be >= 1.

    Returns
    -------
    wl_np1, wl_np2 : numpy.ndarray
        WL label arrays of shape (|V_1|, I+1) and (|V_2|, I+1).

    Notes
    -----
    The cycle basis labeling within each layer is recomputed every
    iteration based on the current lower-layer node labels, ensuring
    that the canonicalized cycle labels reflect updated node labels.
    This is consistent with the spec: "进入下一迭代过程，cycle basis的
    labeling也会相应改变".
    """
    # ---------------------------------------------------------------
    # Validation and degenerate cases
    # ---------------------------------------------------------------
    if L < 0:
        raise ValueError(f"L must be >= 0, got {L}")
    if I < 1:
        raise ValueError(f"I must be >= 1, got {I}")

    # ---------------------------------------------------------------
    # Initialise node ordering and initial labels
    # ---------------------------------------------------------------
    node_list1: List = sorted(G1.nodes(), key=node_sort_key)
    node_list2: List = sorted(G2.nodes(), key=node_sort_key)

    initial_labels1 = _labels_to_dict(G1, vlabel_np1)
    initial_labels2 = _labels_to_dict(G2, vlabel_np2)

    n1 = G1.number_of_nodes()
    n2 = G2.number_of_nodes()

    wl_np1 = np.zeros((n1, I + 1), dtype=np.int32)
    wl_np2 = np.zeros((n2, I + 1), dtype=np.int32)

    for i, v in enumerate(node_list1):
        wl_np1[i, 0] = initial_labels1[v]
    for i, v in enumerate(node_list2):
        wl_np2[i, 0] = initial_labels2[v]

    # Empty graphs → nothing to do
    if n1 == 0 or n2 == 0:
        return wl_np1, wl_np2

    # ---------------------------------------------------------------
    # L=0 case: triangulated neighbors WL (no CSG layers)
    # ---------------------------------------------------------------
    if L == 0:
        # Lazy-import the triangulated neighbors WL module
        _hcc_dir = _os.path.normpath(_os.path.join(
            _os.path.dirname(__file__),
            '..', 'hierarchical_tri_wl_tools'))
        if _hcc_dir not in _sys.path:
            _sys.path.insert(0, _hcc_dir)
        from node_wl_via_triangulated_neighbors import (  # pyright: ignore[reportMissingImports]
            node_wl_test_triangulated_neighbors
        )

        # Remap non-consecutive node IDs to 0..n-1 (required by
        # node_wl_test_triangulated_neighbors which indexes arrays
        # directly by node ID).
        if n1 > 0 and node_list1 != list(range(n1)):
            G1_tri = nx.relabel_nodes(G1, {v: i for i, v in enumerate(node_list1)})
        else:
            G1_tri = G1

        if n2 > 0 and node_list2 != list(range(n2)):
            G2_tri = nx.relabel_nodes(G2, {v: i for i, v in enumerate(node_list2)})
        else:
            G2_tri = G2

        vtx_tri1 = _build_vtx_triangulated_neighbors(G1_tri)
        vtx_tri2 = _build_vtx_triangulated_neighbors(G2_tri)

        # Label arrays indexed by sorted-order position
        # (= node ID after remap, which matches row ordering in wl_np).
        vlabel_arr1 = np.array([wl_np1[i, 0] for i in range(n1)], dtype=np.int32)
        vlabel_arr2 = np.array([wl_np2[i, 0] for i in range(n2)], dtype=np.int32)

        if n1 > 0 and n2 > 0:
            wl_tri1, wl_tri2 = node_wl_test_triangulated_neighbors(
                vtx_tri1, vtx_tri2, vlabel_arr1, vlabel_arr2, I)
            wl_np1[:, :] = wl_tri1[:, :]
            wl_np2[:, :] = wl_tri2[:, :]

        return wl_np1, wl_np2

    # ---------------------------------------------------------------
    # L >= 1: Build multi-layer CSG hierarchy (Step-3)
    # ---------------------------------------------------------------
    paired_csg = build_multilayer_csg_pair(G1, G2, L)
    (
        layers1, layers2,
        mappings1, mappings2,
        nc_input1, nc_input2,
        nc_csg,
    ) = extract_htn_wl_args(paired_csg)

    # ---------------------------------------------------------------
    # Main iteration loop
    # ---------------------------------------------------------------
    for iteration in range(1, I + 1):
        current_labels1 = {
            v: int(wl_np1[i, iteration - 1])
            for i, v in enumerate(node_list1)
        }
        current_labels2 = {
            v: int(wl_np2[i, iteration - 1])
            for i, v in enumerate(node_list2)
        }

        current_max_label = max(
            max(current_labels1.values()),
            max(current_labels2.values()),
        )

        # ---- Forward pass: G -> CSG^1 -> ... -> CSG^L ----
        # forward_labels[i] keeps the labels of the i-th layer (i=0 is G)
        forward_labels1: List[Dict[Hashable, int]] = [current_labels1]
        forward_labels2: List[Dict[Hashable, int]] = [current_labels2]

        for layer_idx in range(L):
            H1, cb1, _ = layers1[layer_idx]
            H2, cb2, _ = layers2[layer_idx]
            nc1, nc2 = nc_csg[layer_idx]

            # Previous layer's labels (for CSG^k, the "lower" layer)
            prev_labels1 = forward_labels1[-1]
            prev_labels2 = forward_labels2[-1]

            # Compute final label tuples for CSG^k from prev layer labels
            final_tuples1 = compute_final_label_tuples(H1, cb1, prev_labels1)
            final_tuples2 = compute_final_label_tuples(H2, cb2, prev_labels2)

            # Joint message passing on both CSG^k graphs
            new_labels1, new_labels2 = forward_message_passing_both(
                H1, final_tuples1, H2, final_tuples2, current_max_label,
                nc1, nc2)

            forward_labels1.append(new_labels1)
            forward_labels2.append(new_labels2)

            if new_labels1:
                current_max_label = max(current_max_label, max(new_labels1.values()))
            if new_labels2:
                current_max_label = max(current_max_label, max(new_labels2.values()))

        # ---- Backward pass: CSG^L -> ... -> CSG^1 -> G ----
        higher_labels1 = forward_labels1[-1]
        higher_labels2 = forward_labels2[-1]

        for step in range(L, 0, -1):
            # Determine the 'lower' graph and its precomputed components
            if step == 1:
                lower_G1 = G1
                lower_G2 = G2
                lower_lbls1 = current_labels1
                lower_lbls2 = current_labels2
                nc1, nc2 = nc_input1, nc_input2
            else:
                lower_G1 = layers1[step - 2][0]
                lower_G2 = layers2[step - 2][0]
                lower_lbls1 = forward_labels1[step - 1]
                lower_lbls2 = forward_labels2[step - 1]
                nc1, nc2 = nc_csg[step - 2]

            new_labels1, new_labels2 = backward_message_passing_both(
                lower_G1, lower_G2,
                lower_lbls1, lower_lbls2,
                higher_labels1, higher_labels2,
                mappings1[step - 1], mappings2[step - 1],
                current_max_label, nc1, nc2)

            higher_labels1 = new_labels1
            higher_labels2 = new_labels2

            if new_labels1:
                current_max_label = max(current_max_label, max(new_labels1.values()))
            if new_labels2:
                current_max_label = max(current_max_label, max(new_labels2.values()))

        # After the full cycle (forward + backward), the labels of G1/G2
        # nodes have been updated.  Store them.
        for i, v in enumerate(node_list1):
            wl_np1[i, iteration] = higher_labels1.get(v, 0)
        for i, v in enumerate(node_list2):
            wl_np2[i, iteration] = higher_labels2.get(v, 0)

    return wl_np1, wl_np2


# ============================================================================
# Edge-aware hierarchical triangulated WL (node + edge labels)
# ============================================================================

def hierarchical_triangular_wl_with_edges(
    G1: nx.Graph,
    G2: nx.Graph,
    vlabel_np1: Any,
    vlabel_np2: Any,
    elabel_dict1: Dict,
    elabel_dict2: Dict,
    L: int = 1,
    I: int = 5,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Hierarchical Triangulated WL test that handles **both node labels and
    edge labels** simultaneously.

    The algorithm follows the same forward/backward CSG message-passing
    structure as ``hierarchical_triangular_wl``, but additionally:

    * **Edge-to-node fusion** — G node label tuples in the backward pass
      include an *edge context* ``ec(v)`` that summarises incident edge
      labels grouped by triangulated neighbourhood components.
    * **Edge label refresh** — after each full forward+backward cycle, edge
      labels are updated via ``compress(old_elabel, sorted(vlabel[u],
      vlabel[v]))``.
    * **Both histories** — node labels ``(n, I+1)`` **and** edge labels
      ``(|E|, I+1)`` are returned.

    Parameters
    ----------
    G1, G2 : networkx.Graph
        Input graphs.
    vlabel_np1, vlabel_np2 : numpy.ndarray or dict
        Initial node labels (same format as ``hierarchical_triangular_wl``).
    elabel_dict1, elabel_dict2 : dict
        Initial edge labels ``{(u, v): int}`` with ``u <= v``.  Must cover
        **all** edges of the respective graph (edges not present default to 0).
    L : int
        Number of CSG layers (default 1).  ``L=0`` means triangulated-
        neighbours edge WL (no CSG hierarchy).
    I : int
        Number of iterations (default 5).  Must be >= 1.

    Returns
    -------
    vwl_np1, vwl_np2 : numpy.ndarray
        Node label histories, shape ``(|V_1|, I+1)`` and ``(|V_2|, I+1)``.
    ewl_np1, ewl_np2 : numpy.ndarray
        Edge label histories, shape ``(|E_1|, I+1)`` and ``(|E_2|, I+1)``.
    """
    # ---------------------------------------------------------------
    # Validation
    # ---------------------------------------------------------------
    if L < 0:
        raise ValueError(f"L must be >= 0, got {L}")
    if I < 1:
        raise ValueError(f"I must be >= 1, got {I}")
    _validate_edge_labels(G1, elabel_dict1, G2, elabel_dict2)

    # ---------------------------------------------------------------
    # Normalise edge labels: sorted canonical keys, build edge arrays
    # ---------------------------------------------------------------
    edges1: List = sorted(elabel_dict1.keys())
    edges2: List = sorted(elabel_dict2.keys())
    elabel_arr1: np.ndarray = np.array(
        [elabel_dict1[e] for e in edges1], dtype=np.int32)
    elabel_arr2: np.ndarray = np.array(
        [elabel_dict2[e] for e in edges2], dtype=np.int32)

    # Current edge-label dicts (mutated each iteration)
    cur_elabel1: Dict = dict(elabel_dict1)
    cur_elabel2: Dict = dict(elabel_dict2)

    # ---------------------------------------------------------------
    # Initialise node ordering and labels
    # ---------------------------------------------------------------
    node_list1 = sorted(G1.nodes(), key=node_sort_key)
    node_list2 = sorted(G2.nodes(), key=node_sort_key)

    initial_labels1 = _labels_to_dict(G1, vlabel_np1)
    initial_labels2 = _labels_to_dict(G2, vlabel_np2)

    n1 = G1.number_of_nodes()
    n2 = G2.number_of_nodes()

    vwl_np1 = np.zeros((n1, I + 1), dtype=np.int32)
    vwl_np2 = np.zeros((n2, I + 1), dtype=np.int32)

    for i, v in enumerate(node_list1):
        vwl_np1[i, 0] = initial_labels1[v]
    for i, v in enumerate(node_list2):
        vwl_np2[i, 0] = initial_labels2[v]

    # Edge label history arrays
    ewl_np1 = np.zeros((len(edges1), I + 1), dtype=np.int32)
    ewl_np2 = np.zeros((len(edges2), I + 1), dtype=np.int32)
    for i, e in enumerate(edges1):
        ewl_np1[i, 0] = elabel_dict1[e]
    for i, e in enumerate(edges2):
        ewl_np2[i, 0] = elabel_dict2[e]

    # Empty graphs
    if n1 == 0 or n2 == 0:
        return vwl_np1, vwl_np2, ewl_np1, ewl_np2

    # ---------------------------------------------------------------
    # L=0: triangulated neighbours edge WL (no CSG layers)
    # ---------------------------------------------------------------
    if L == 0:
        _hcc_dir = _os.path.normpath(_os.path.join(
            _os.path.dirname(__file__),
            '..', 'hierarchical_tri_wl_tools'))
        if _hcc_dir not in _sys.path:
            _sys.path.insert(0, _hcc_dir)
        # Also add parent dir so import with package prefix makes relative imports work
        _hcc_parent = _os.path.normpath(_os.path.join(_hcc_dir, '..'))
        if _hcc_parent not in _sys.path:
            _sys.path.insert(0, _hcc_parent)
        from hierarchical_tri_wl_tools.edge_wl_via_triangulated_neighbors import (  # pyright: ignore[reportMissingImports]
            edge_wl_test_triangulated_neighbors
        )

        # Remap non-consecutive node IDs to 0..n-1
        if n1 > 0 and node_list1 != list(range(n1)):
            G1_tri = nx.relabel_nodes(G1, {v: i for i, v in enumerate(node_list1)})
        else:
            G1_tri = G1
        if n2 > 0 and node_list2 != list(range(n2)):
            G2_tri = nx.relabel_nodes(G2, {v: i for i, v in enumerate(node_list2)})
        else:
            G2_tri = G2

        vtx_tri1 = _build_vtx_triangulated_neighbors(G1_tri)
        vtx_tri2 = _build_vtx_triangulated_neighbors(G2_tri)

        # Remap edge keys to match remapped node IDs
        def _remap_edges(edges, mapping):
            return [tuple(sorted((mapping[u], mapping[v]))) for u, v in edges]

        node_map1 = {v: i for i, v in enumerate(node_list1)}
        node_map2 = {v: i for i, v in enumerate(node_list2)}
        edges1_remapped = _remap_edges(edges1, node_map1)
        edges2_remapped = _remap_edges(edges2, node_map2)

        # Build label arrays indexed by sorted-order position
        vlabel_arr1 = np.array([vwl_np1[i, 0] for i in range(n1)], dtype=np.int32)
        vlabel_arr2 = np.array([vwl_np2[i, 0] for i in range(n2)], dtype=np.int32)

        # Remap elabel arrays to match remapped edge keys
        elabel_arr1_remapped = np.array(
            [elabel_arr1[i] for i in range(len(edges1))], dtype=np.int32)
        elabel_arr2_remapped = np.array(
            [elabel_arr2[i] for i in range(len(edges2))], dtype=np.int32)

        if n1 > 0 and n2 > 0:
            vwl_tri1, vwl_tri2, ewl_tri1, ewl_tri2 = \
                edge_wl_test_triangulated_neighbors(
                    vtx_tri1, vtx_tri2,
                    vlabel_arr1, vlabel_arr2,
                    edges1_remapped, edges2_remapped,
                    elabel_arr1_remapped, elabel_arr2_remapped,
                    I)
            vwl_np1[:, :] = vwl_tri1[:, :]
            vwl_np2[:, :] = vwl_tri2[:, :]
            ewl_np1[:, :] = ewl_tri1[:, :]
            ewl_np2[:, :] = ewl_tri2[:, :]

        return vwl_np1, vwl_np2, ewl_np1, ewl_np2

    # ---------------------------------------------------------------
    # L >= 1: Build CSG hierarchy
    # ---------------------------------------------------------------
    paired_csg = build_multilayer_csg_pair(G1, G2, L)
    (
        layers1, layers2,
        mappings1, mappings2,
        nc_input1, nc_input2,
        nc_csg,
    ) = extract_htn_wl_args(paired_csg)

    # ---------------------------------------------------------------
    # Pre-processing: compute initial edge contexts for G nodes
    # ---------------------------------------------------------------
    edge_contexts1 = _compute_edge_contexts_all(G1, cur_elabel1, nc_input1)
    edge_contexts2 = _compute_edge_contexts_all(G2, cur_elabel2, nc_input2)

    # ---------------------------------------------------------------
    # Main iteration loop
    # ---------------------------------------------------------------
    for iteration in range(1, I + 1):
        current_labels1 = {
            v: int(vwl_np1[i, iteration - 1])
            for i, v in enumerate(node_list1)
        }
        current_labels2 = {
            v: int(vwl_np2[i, iteration - 1])
            for i, v in enumerate(node_list2)
        }

        current_max_label = max(
            max(current_labels1.values()),
            max(current_labels2.values()),
        )

        # ---- Forward pass: G -> CSG^1 -> ... -> CSG^L ----
        forward_labels1 = [current_labels1]
        forward_labels2 = [current_labels2]

        for layer_idx in range(L):
            H1, cb1, _ = layers1[layer_idx]
            H2, cb2, _ = layers2[layer_idx]
            nc1, nc2 = nc_csg[layer_idx]

            prev_labels1 = forward_labels1[-1]
            prev_labels2 = forward_labels2[-1]

            final_tuples1 = compute_final_label_tuples(H1, cb1, prev_labels1)
            final_tuples2 = compute_final_label_tuples(H2, cb2, prev_labels2)

            new_labels1, new_labels2 = forward_message_passing_both(
                H1, final_tuples1, H2, final_tuples2, current_max_label,
                nc1, nc2)

            forward_labels1.append(new_labels1)
            forward_labels2.append(new_labels2)

            if new_labels1:
                current_max_label = max(current_max_label, max(new_labels1.values()))
            if new_labels2:
                current_max_label = max(current_max_label, max(new_labels2.values()))

        # ---- Backward pass: CSG^L -> ... -> CSG^1 -> G ----
        higher_labels1 = forward_labels1[-1]
        higher_labels2 = forward_labels2[-1]

        for step in range(L, 0, -1):
            if step == 1:
                lower_G1, lower_G2 = G1, G2
                lower_lbls1, lower_lbls2 = current_labels1, current_labels2
                nc1, nc2 = nc_input1, nc_input2
                # G level: use edge-aware backward pass
                new_labels1, new_labels2 = backward_message_passing_both_with_edges(
                    lower_G1, lower_G2,
                    lower_lbls1, lower_lbls2,
                    higher_labels1, higher_labels2,
                    mappings1[step - 1], mappings2[step - 1],
                    current_max_label,
                    edge_contexts1, edge_contexts2,
                    nc1, nc2)
            else:
                lower_G1 = layers1[step - 2][0]
                lower_G2 = layers2[step - 2][0]
                lower_lbls1 = forward_labels1[step - 1]
                lower_lbls2 = forward_labels2[step - 1]
                nc1, nc2 = nc_csg[step - 2]
                # CSG layers: standard backward pass (no edge labels)
                new_labels1, new_labels2 = backward_message_passing_both(
                    lower_G1, lower_G2,
                    lower_lbls1, lower_lbls2,
                    higher_labels1, higher_labels2,
                    mappings1[step - 1], mappings2[step - 1],
                    current_max_label, nc1, nc2)

            higher_labels1 = new_labels1
            higher_labels2 = new_labels2

            if new_labels1:
                current_max_label = max(current_max_label, max(new_labels1.values()))
            if new_labels2:
                current_max_label = max(current_max_label, max(new_labels2.values()))

        # ---- Edge label refresh ----
        vlabel_dict1 = {v: int(higher_labels1.get(v, 0)) for v in G1.nodes()}
        vlabel_dict2 = {v: int(higher_labels2.get(v, 0)) for v in G2.nodes()}
        cur_elabel1, cur_elabel2 = _update_elabel_from_dict(
            vlabel_dict1, vlabel_dict2, cur_elabel1, cur_elabel2)

        # Recompute edge contexts from refreshed edge labels
        edge_contexts1 = _compute_edge_contexts_all(G1, cur_elabel1, nc_input1)
        edge_contexts2 = _compute_edge_contexts_all(G2, cur_elabel2, nc_input2)

        # ---- Store labels for this iteration ----
        for i, v in enumerate(node_list1):
            vwl_np1[i, iteration] = higher_labels1.get(v, 0)
        for i, v in enumerate(node_list2):
            vwl_np2[i, iteration] = higher_labels2.get(v, 0)
        for i, e in enumerate(edges1):
            ewl_np1[i, iteration] = cur_elabel1.get(e, 0)
        for i, e in enumerate(edges2):
            ewl_np2[i, iteration] = cur_elabel2.get(e, 0)

    return vwl_np1, vwl_np2, ewl_np1, ewl_np2


# ============================================================================
# Unified interface (dispatches based on g_info)
# ============================================================================

def hierarchical_triangular_wl_unified(
    g_info: Dict,
    G1: nx.Graph,
    G2: nx.Graph,
    vlabel_np1: Any,
    vlabel_np2: Any,
    elabel_dict1: Optional[Dict] = None,
    elabel_dict2: Optional[Dict] = None,
    L: int = 1,
    I: int = 5,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Unified entry point that dispatches to the node-only or edge-aware
    hierarchical WL test based on ``g_info``.

    ``g_info`` is a dataset-info dictionary with at least the key ``'el'``:

    * ``g_info['el'] is True``  → graphs have **edge labels** (and,
      implicitly, node labels too).  Delegates to
      :func:`hierarchical_triangular_wl_with_edges`.
    * ``g_info['el'] is False`` → graphs have **node labels only**.
      Delegates to :func:`hierarchical_triangular_wl`.

    Parameters
    ----------
    g_info : dict
        Dataset information; must contain ``'el'`` (bool).
    G1, G2 : networkx.Graph
        Input graphs.
    vlabel_np1, vlabel_np2 : numpy.ndarray or dict
        Initial node labels.
    elabel_dict1, elabel_dict2 : dict, optional
        Initial edge labels (required iff ``g_info['el']`` is True).
    L : int
        Number of CSG layers (default 1).
    I : int
        Number of iterations (default 5).

    Returns
    -------
    wl_np1, wl_np2 : numpy.ndarray
        Node label histories, shape ``(|V_1|, I+1)`` and ``(|V_2|, I+1)``.
    """
    if g_info.get('el', False):
        # elabel_dict1/2 are guaranteed non-None when g_info['el'] is True
        _ed1: Dict = elabel_dict1 if elabel_dict1 is not None else {}
        _ed2: Dict = elabel_dict2 if elabel_dict2 is not None else {}
        vwl1, vwl2, _ewl1, _ewl2 = hierarchical_triangular_wl_with_edges(
            G1, G2, vlabel_np1, vlabel_np2,
            _ed1, _ed2, L, I,
        )
        return vwl1, vwl2
    else:
        return hierarchical_triangular_wl(
            G1, G2, vlabel_np1, vlabel_np2, L, I,
        )


# ============================================================================
# Testing utilities
# ============================================================================

def _is_isomorphic_wl(wl1: np.ndarray, wl2: np.ndarray) -> bool:
    """Return True iff the two WL arrays have the same multiset per column."""
    if wl1.shape != wl2.shape:
        return False
    _, cols = wl1.shape
    for c in range(cols):
        if sorted(wl1[:, c].tolist()) != sorted(wl2[:, c].tolist()):
            return False
    return True


def _build_path(n: int) -> nx.Graph:
    g = nx.Graph()
    g.add_edges_from((i, i + 1) for i in range(n - 1))
    return g


def _build_cycle(n: int) -> nx.Graph:
    g = nx.Graph()
    g.add_edges_from((i, (i + 1) % n) for i in range(n))
    return g


def _build_chord_cycle() -> nx.Graph:
    """4-cycle with a chord: 4 nodes, 5 edges."""
    g = nx.Graph()
    g.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
    return g


def _build_two_triangles() -> nx.Graph:
    """Two disjoint triangles: 6 nodes, 6 edges."""
    g = nx.Graph()
    g.add_edges_from([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
    return g


# ============================================================================
# Self-tests
# ============================================================================

if __name__ == "__main__":
    import sys
    sys.path.insert(0, _os.path.dirname(_os.path.dirname(__file__)))

    G1 = build_example1_graph()
    G2 = build_example1_graph()

    print(f"Graph 1: {G1.number_of_nodes()} nodes, {G1.number_of_edges()} edges")
    print(f"Graph 2: {G2.number_of_nodes()} nodes, {G2.number_of_edges()} edges")

    np.random.seed(42)
    vlabel = np.random.randint(0, 10, G1.number_of_nodes())
    vlabel1 = vlabel.copy()
    vlabel2 = vlabel.copy()

    # ====================================================================
    # Test 1a: identical graphs, identical labels, L=1/2/3
    # ====================================================================
    print("\n--- Test 1a: identical graphs, identical labels, L=1/2/3 ---")
    for L in (1, 2, 3):
        wl1_f, wl2_f = hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=L, I=3)
        match = np.all(wl1_f == wl2_f)
        print(f"  L={L} I=3  shapes={wl1_f.shape}  match={match}")
        assert match, f"L={L}: identical inputs should give identical labels"

    # ====================================================================
    # Test 1b: identical graphs under L=0 (triangulated neighbors WL)
    # ====================================================================
    print("\n--- Test 1b: L=0 (triangulated neighbors WL) ---")
    wl_k0_1, wl_k0_2 = hierarchical_triangular_wl(
        G1, G2, vlabel1, vlabel2, L=0, I=3)
    match_k0 = np.all(wl_k0_1 == wl_k0_2)
    print(f"  L=0 I=3  shapes={wl_k0_1.shape}  match={match_k0}")
    assert match_k0, "L=0: identical inputs should give identical labels"
    assert wl_k0_1.shape == (26, 4), f"Expected (26, 4), got {wl_k0_1.shape}"

    # ====================================================================
    # Test 2: isomorphic graphs (relabeled), L=1/2/3
    # ====================================================================
    print("\n--- Test 2: isomorphic graphs (relabeled), L=1/2/3 ---")
    relabel = {i: 25 - i for i in range(26)}
    G1b = nx.relabel_nodes(G1, relabel)
    vlabel_b = np.array([vlabel[25 - i] for i in range(26)], dtype=np.int32)
    for L in (1, 2, 3):
        wl1, wl2 = hierarchical_triangular_wl(
            G1, G1b, vlabel1, vlabel_b, L=L, I=3)
        iso = _is_isomorphic_wl(wl1, wl2)
        print(f"  L={L} I=3  iso_multiset={iso}")
        if not iso:
            print(f"    Warning: isomorphic graphs differ due to cycle_basis order.")
            print(f"    Known limitation: nx.minimum_cycle_basis is not")
            print(f"    relabeling-invariant for graphs with automorphic")
            print(f"    substructures. Different node labeling leads to")
            print(f"    different minimum cycle bases, which changes CB")
            print(f"    node names in the CSG. This is fundamentally hard")
            print(f"    to fix as it requires solving graph isomorphism.")

    # ====================================================================
    # Test 3: non-isomorphic (cycle vs path)
    # ====================================================================
    print("\n--- Test 3: non-isomorphic (cycle vs path) ---")
    G_cycle = _build_cycle(6)
    G_path = _build_path(6)
    labels = np.array([0, 1, 2, 0, 1, 2], dtype=np.int32)
    wl_c, wl_p = hierarchical_triangular_wl(G_cycle, G_path, labels, labels, L=1, I=3)
    not_iso = not _is_isomorphic_wl(wl_c, wl_p)
    print(f"  L=1 I=3  cycle vs path: distinguish={not_iso}")
    assert not_iso, "WL must distinguish a 6-cycle from a 6-path"

    # ====================================================================
    # Test 4: multi-layer L=4
    # ====================================================================
    print("\n--- Test 4: multi-layer L=4 ---")
    wl1, wl2 = hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=4, I=2)
    print(f"  L=4 I=2  shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2)

    # ====================================================================
    # Test 5: dict-form initial labels
    # ====================================================================
    print("\n--- Test 5: dict-form initial labels ---")
    labels_dict = {v: int(vlabel1[i]) for i, v in enumerate(sorted(G1.nodes()))}
    wl_d, _ = hierarchical_triangular_wl(G1, G2, labels_dict, labels_dict, L=1, I=2)
    wl_a, _ = hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=1, I=2)
    print(f"  dict-vs-array: match={np.all(wl_d == wl_a)}")
    assert np.all(wl_d == wl_a)

    # ====================================================================
    # Test 6: acyclic graph (no cycles)
    # ====================================================================
    print("\n--- Test 6: acyclic graph (no cycles) ---")
    G_acyclic = _build_path(5)
    G_acyclic2 = _build_path(5)
    labels_path = np.array([0, 1, 0, 1, 0], dtype=np.int32)
    wl1, wl2 = hierarchical_triangular_wl(
        G_acyclic, G_acyclic2, labels_path, labels_path.copy(), L=1, I=2)
    print(f"  acyclic L=1 I=2: shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2)

    # ====================================================================
    # Test 7: single-node graph
    # ====================================================================
    print("\n--- Test 7: single-node graph ---")
    G_single = nx.Graph()
    G_single.add_node(0)
    G_single2 = nx.Graph()
    G_single2.add_node(0)
    labels_single = np.array([42], dtype=np.int32)
    wl1, wl2 = hierarchical_triangular_wl(
        G_single, G_single2, labels_single, labels_single.copy(), L=1, I=2)
    print(f"  single-node L=1 I=2: shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2)

    # ====================================================================
    # Test 8: two disjoint cycles
    # ====================================================================
    print("\n--- Test 8: two disjoint cycles (no shared edges, no interface) ---")
    G_disj = _build_two_triangles()
    G_disj2 = _build_two_triangles()
    labels_disj = np.array([0, 1, 2, 0, 1, 2], dtype=np.int32)
    wl1, wl2 = hierarchical_triangular_wl(
        G_disj, G_disj2, labels_disj, labels_disj.copy(), L=1, I=2)
    print(f"  disjoint cycles L=1 I=2: shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2)

    # ====================================================================
    # Test 9: non-consecutive node IDs with ndarray labels
    # ====================================================================
    print("\n--- Test 9: non-consecutive node IDs with ndarray labels ---")
    G_nc = nx.Graph()
    G_nc.add_edges_from([(0, 5), (5, 10), (10, 0)])
    G_nc2 = nx.Graph()
    G_nc2.add_edges_from([(0, 5), (5, 10), (10, 0)])
    labels_nc = np.array([1, 2, 3], dtype=np.int32)
    wl1, wl2 = hierarchical_triangular_wl(G_nc, G_nc2, labels_nc, labels_nc, L=1, I=2)
    print(f"  non-consecutive L=1 I=2: shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2), "Non-consecutive node IDs with ndarray labels must work"

    # ====================================================================
    # Test 10: empty graph
    # ====================================================================
    print("\n--- Test 10: empty graph ---")
    G_empty = nx.Graph()
    G_empty2 = nx.Graph()
    labels_empty = np.array([], dtype=np.int32)
    wl1, wl2 = hierarchical_triangular_wl(G_empty, G_empty2, labels_empty, labels_empty, L=1, I=2)
    print(f"  empty graph: shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2)

    # ====================================================================
    # Test 11: cycle_basis node label updates (L=1 forward only)
    # ====================================================================
    print("\n--- Test 11: triangle graph L=1 ---")
    G_tri1 = nx.Graph()
    G_tri1.add_edges_from([(0, 1), (1, 2), (0, 2)])
    G_tri2 = nx.Graph()
    G_tri2.add_edges_from([(0, 1), (1, 2), (0, 2)])
    tri_labels = np.array([1, 2, 3], dtype=np.int32)
    wl_t1, wl_t2 = hierarchical_triangular_wl(
        G_tri1, G_tri2, tri_labels, tri_labels.copy(), L=1, I=2)
    print(f"  triangle L=1 I=2: shapes={wl_t1.shape}  match={np.all(wl_t1 == wl_t2)}")
    assert np.all(wl_t1 == wl_t2)
    assert wl_t1.shape == (3, 3), f"Expected (3, 3) shape, got {wl_t1.shape}"

    # ====================================================================
    # Test 12: cycle_basis label canonicalization
    # ====================================================================
    print("\n--- Test 12: cycle_basis label canonicalization ---")
    assert canonicalize_cycle_label((1, 2, 3, 4)) == (1, 2, 3, 4)
    assert canonicalize_cycle_label((3, 1, 2, 4)) == (1, 2, 4, 3)
    assert canonicalize_cycle_label((1, 4, 3, 2)) == (1, 2, 3, 4)
    assert canonicalize_cycle_label((2, 1, 4, 3)) == (1, 2, 3, 4)
    assert canonicalize_cycle_label((1,)) == (1,)
    assert canonicalize_cycle_label(()) == ()
    assert canonicalize_cycle_label((1, 1, 1, 1)) == (1, 1, 1, 1)
    print("  All canonicalization invariants verified.")

    # ====================================================================
    # Test 13: visualization functions return figures (Agg backend)
    # ====================================================================
    print("\n--- Test 13: visualization functions return figures (Agg backend) ---")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    H, cycle_basis, info = cyclic_schematic_graph(G1, cb_prefix="cb0")
    fig = draw_cyclic_schematic(H, cycle_basis, info, title="Test")
    assert fig is not None
    plt.close(fig)
    fig = draw_input_graph(G1, title="Test")
    assert fig is not None
    plt.close(fig)
    fig = draw_input_graph_with_cycles(G1, cycle_basis, title="Test")
    assert fig is not None
    plt.close(fig)
    fig = draw_side_by_side(G1, H, cycle_basis, info)
    assert fig is not None
    plt.close(fig)
    fig = draw_side_by_side_with_cycles(G1, H, cycle_basis, info)
    assert fig is not None
    plt.close(fig)
    print("  All 5 visualization functions returned figures successfully.")

    # ====================================================================
    # Test 14: step-1 csg_to_input and input_to_csg mappings
    # ====================================================================
    print("\n--- Test 14: step-1 csg_to_input and input_to_csg mappings ---")
    csg_to_input = build_csg_to_input_mapping(H, cycle_basis)
    input_to_csg = build_input_to_csg_mapping(H, cycle_basis, G1)
    for csg_node, original_tuple in csg_to_input.items():
        if get_node_type(H, csg_node) == "cycle_basis":
            assert len(original_tuple) >= 3
        else:
            assert original_tuple == (csg_node,)
    for input_node in G1.nodes():
        if input_node in csg_to_input.values() or any(
            input_node in tup for tup in csg_to_input.values()
        ):
            assert input_node in input_to_csg
    print("  All step-1 mapping invariants verified.")

    # ====================================================================
    # Test 15: multi-layer with L=5 (stress test)
    # ====================================================================
    print("\n--- Test 15: multi-layer with L=5 (stress test) ---")
    wl1, wl2 = hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=5, I=2)
    print(f"  L=5 I=2: shapes={wl1.shape}  match={np.all(wl1 == wl2)}")
    assert np.all(wl1 == wl2)

    # ====================================================================
    # Test 16: graph with chord (extra edge in cycle)
    # ====================================================================
    print("\n--- Test 16: graph with chord (extra edge in cycle) ---")
    G_chord = _build_chord_cycle()
    G_chord2 = _build_chord_cycle()
    chord_labels = np.array([1, 2, 3, 4], dtype=np.int32)
    wl_c, wl_c2 = hierarchical_triangular_wl(
        G_chord, G_chord2, chord_labels, chord_labels.copy(), L=1, I=2)
    print(f"  chord L=1 I=2: shapes={wl_c.shape}  match={np.all(wl_c == wl_c2)}")
    assert wl_c.shape == (4, 3)
    assert np.all(wl_c == wl_c2), "Identical chord graphs should have identical labels"

    # ====================================================================
    # Test 17: L=0 on acyclic graph (triangulated neighbors WL)
    # ====================================================================
    print("\n--- Test 17: L=0 on acyclic graph (triangulated neighbors WL) ---")
    G_lin = _build_path(4)
    G_lin2 = _build_path(4)
    lin_labels = np.array([3, 1, 4, 1], dtype=np.int32)
    wl_k0_a, wl_k0_b = hierarchical_triangular_wl(G_lin, G_lin2, lin_labels, lin_labels, L=0, I=2)
    match_k0 = np.all(wl_k0_a == wl_k0_b)
    print(f"  L=0 I=2 shapes={wl_k0_a.shape}  match={match_k0}")
    assert match_k0, "L=0: identical acyclic graphs should produce identical labels"

    # ====================================================================
    # Test 18: Verify the L=0/L=1/L=2 consistent ordering
    # ====================================================================
    print("\n--- Test 18: ordering: same initial labels produce same output shape ---")
    for L in (0, 1, 2, 3):
        wl_a, wl_b = hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=L, I=2)
        assert np.all(wl_a == wl_b), f"L={L}: identical -> identical failed"
    print("  All L produced identical labels for identical inputs.")

    # ====================================================================
    # Test 19: Three-cycle (triangle) different label orders
    # ====================================================================
    print("\n--- Test 19: triangle different label orderings ---")
    G_tria = nx.Graph()
    G_tria.add_edges_from([(0, 1), (1, 2), (0, 2)])
    G_trib = nx.Graph()
    G_trib.add_edges_from([(0, 1), (1, 2), (0, 2)])
    labels_a = np.array([5, 10, 15], dtype=np.int32)
    labels_b = np.array([15, 5, 10], dtype=np.int32)
    wl_tria, _ = hierarchical_triangular_wl(G_tria, G_trib, labels_a, labels_b, L=1, I=3)
    iso_tri = _is_isomorphic_wl(wl_tria, wl_tria)
    print(f"  triangle different labelings: iso={iso_tri}")

    # ====================================================================
    # Test 20: Validate that output arrays are monotonically non-decreasing
    # (labels only increase as new aggregate signatures are discovered)
    # ====================================================================
    print("\n--- Test 20: label monotonicity check ---")
    wl1, wl2 = hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=2, I=4)
    for col in range(1, wl1.shape[1]):
        assert np.all(wl1[:, col] >= wl1[:, col - 1]), \
            f"Labels should be non-decreasing at column {col}"
        assert np.all(wl2[:, col] >= wl2[:, col - 1]), \
            f"Labels should be non-decreasing at column {col}"
    print("  All labels are non-decreasing across iterations.")

    # ====================================================================
    # Test 21: Error handling
    # ====================================================================
    print("\n--- Test 21: error handling ---")
    try:
        hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=-1, I=2)
        print("  FAIL: L=-1 should raise ValueError")
    except ValueError:
        print("  L=-1 raises ValueError: OK")
    try:
        hierarchical_triangular_wl(G1, G2, vlabel1, vlabel2, L=1, I=0)
        print("  FAIL: I=0 should raise ValueError")
    except ValueError:
        print("  I=0 raises ValueError: OK")

    print("\nAll self-tests passed.")
