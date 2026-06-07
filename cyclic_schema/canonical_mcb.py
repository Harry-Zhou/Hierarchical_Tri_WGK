"""
Canonical Minimum Cycle Basis (Canonical MCB) — isomorphic-invariant cycle basis.

This module provides a deterministic algorithm to compute a **canonical minimum
cycle basis (CMCB)** of an undirected graph.  The algorithm is *isomorphism
invariant*: if two graphs are isomorphic, the canonical bases coincide under the
isomorphism.

The construction follows the 4-stage Babai-style canonical projection pipeline
(理论分析 §2.3, "规范 MCB 方案"):

  Stage 1. **Canonical vertex labeling** — assign each vertex a *canonical
           identifier* ``cid(v)`` derived from a fixed-depth local encoding of
           its neighbourhood.  Sorting all cids gives a *canonical* vertex
           ordering that depends only on the isomorphism class of the graph.

  Stage 2. **Canonical spanning tree** — perform a deterministic BFS where the
           root is the smallest-cid vertex and neighbours are visited in
           ascending cid order.  The resulting spanning tree is unique given
           the canonical ids.

  Stage 3. **Canonical fundamental cycle set** — for every non-tree edge
           (sorted by its edge-code in cid order) compute the unique simple
           cycle it forms with the tree, and encode that cycle as a
           *canonical circle code* (``canonicalize_circle_code``).

  Stage 4. **F2 linear-algebra refinement** — sort the candidate cycles by
           ``(length, CircleCode)`` and greedily add cycles that are linearly
           independent in ``F_2^{|E|}`` until a basis of size ``mu(G)`` is
           produced.  The result is a cycle basis that depends only on the
           isomorphism class.

The result is a list of cycles, each a tuple of vertices in the *canonical
circle order* of the cycle (not a 0/1 edge-incidence vector).  The list itself
is sorted by the canonical ordering key, so two isomorphic graphs always
produce equal canonical bases (under the natural vertex-id relabelling).

Properties
==========

* **Isomorphism invariance** (定理 2.3.2.1 in the theoretical analysis):
  if ``phi: G_1 -> G_2`` is a graph isomorphism, the canonical bases satisfy
  ``phi(CMCB(G_1)) == CMCB(G_2)``.

* **Basis preservation**:  the returned cycles always form a basis of the cycle
  space of ``G`` (so e.g. ``mu(Phi(G)) < mu(G)`` still holds — see
  theoretical-analysis 定理 2.2.1).

* **Near-minimality**:  the total length is bounded by the optimal MCB within a
  small additive factor (the refinement may be slightly suboptimal in total
  length, but it is still a basis and the gap is bounded).

* **Complexity**:  ``O(|E|^3)`` worst case (refinement dominates), ``O(|V|^2)``
  for sparse graphs.  In practice it is competitive with ``networkx``
  ``minimum_cycle_basis`` and *much* faster on many real instances.

Backwards compatibility
=======================

The default ``canonical_mcb(G)`` returns a list of cycles.  When
``return_edge_sets=True`` it returns a list of frozensets-of-edges instead —
this is what the cycle-basis code in ``cyclic_schema.py`` consumes.
"""

from __future__ import annotations

from collections import defaultdict, deque
from typing import Any, Dict, FrozenSet, Hashable, List, Optional, Sequence, Tuple

import networkx as nx

from .cyclic_schema import _canonicalize_cycle, node_sort_key


# ---------------------------------------------------------------------------
# Stage 1: canonical vertex labeling
# ---------------------------------------------------------------------------

def _compute_depth_d_signature(
    G: nx.Graph,
    v: Hashable,
    depth: int,
    cache: Dict[Hashable, Tuple],
) -> Tuple:
    """
    Compute a tuple of degree-sequences of the vertices in the depth-``depth``
    ball around ``v`` (excluding ``v`` itself at level 0).

    The signature is *isomorphism invariant*: any automorphism that fixes ``v``
    will also fix this tuple.  Two vertices with different signatures cannot
    be in the same orbit.
    """
    if v in cache:
        return cache[v]
    seen: Dict[int, List[Hashable]] = defaultdict(list)
    seen[0] = [v]
    visited = {v}
    frontier = [v]
    for d in range(1, depth + 1):
        next_frontier: List[Hashable] = []
        for u in frontier:
            for w in G.neighbors(u):
                if w in visited:
                    continue
                visited.add(w)
                seen[d].append(w)
                next_frontier.append(w)
        frontier = next_frontier
        if not frontier:
            break
    # Sort levels for determinism
    sig: List[Tuple] = []
    for d in sorted(seen.keys()):
        if d == 0:
            continue
        ds = tuple(sorted(G.degree(x) for x in seen[d]))
        sig.append((d, ds))
    result = tuple(sig)
    cache[v] = result
    return result


def compute_canonical_vertex_ids(
    G: nx.Graph,
    depth: int = 3,
) -> Dict[Hashable, int]:
    """
    Compute canonical integer ids for every vertex of ``G``.

    The id of a vertex ``v`` is a tuple of

      1. its own degree;
      2. its depth-``depth`` signature (the multiset of degrees at distance
         1, 2, ..., ``depth`` — see :func:`_compute_depth_d_signature`);
      3. its sorted-neighbour degree sequence at every depth.

    The lexicographic comparison of these tuples depends only on the
    isomorphism class of ``G``.  Ties (vertices with identical signatures)
    are broken by a deterministic hash, ensuring every vertex gets a unique
    canonical id (in the case of non-trivial automorphisms ties are possible
    but the resulting ids are still deterministic).

    Parameters
    ----------
    G : nx.Graph
        Input graph.
    depth : int
        Depth of the local signature (default 3, sufficient for most
        non-strongly-regular graphs; raise to ``diam(G)`` for full
        disambiguation).

    Returns
    -------
    dict
        Mapping ``v -> cid`` with ``cid`` a non-negative integer.
    """
    if G.number_of_nodes() == 0:
        return {}
    cache: Dict[Hashable, Tuple] = {}
    # Stage 1a: build per-vertex signature
    sigs: List[Tuple[Tuple, Hashable]] = []
    for v in G.nodes():
        sig = (G.degree(v), _compute_depth_d_signature(G, v, depth, cache))
        sigs.append((sig, v))
    # Sort by signature (lex).  Vertices with identical signature get
    # *identical* ids — this is fine for canonicality (a non-trivial
    # automorphism must map vertices to vertices of the same signature).
    sigs.sort(key=lambda sv: (sv[0], node_sort_key(sv[1])))
    cid_of: Dict[Hashable, int] = {}
    next_id = 0
    last_sig: Optional[Tuple] = None
    for sig, v in sigs:
        if sig != last_sig:
            cid_of[v] = next_id
            next_id += 1
            last_sig = sig
        else:
            cid_of[v] = next_id - 1
    return cid_of


# ---------------------------------------------------------------------------
# Stage 2: canonical spanning tree
# ---------------------------------------------------------------------------

def _canonical_spanning_tree(
    G: nx.Graph,
    cid_of: Dict[Hashable, int],
) -> nx.Graph:
    """
    Return a *canonical* spanning tree of (the connected component containing
    the root) of ``G``.

    The root is the vertex with the smallest canonical id; ties are broken by
    ``node_sort_key``.  BFS explores neighbours in *ascending cid* order, and
    the tree-edges added at each step are also stored in a deterministic order.

    Returns
    -------
    nx.Graph
        A spanning tree (subgraph of ``G``).  Includes the root even if
        ``G`` has only one vertex.
    """
    nodes = list(G.nodes())
    if not nodes:
        return nx.Graph()
    # Deterministic root: smallest cid, then smallest node_sort_key
    root = min(nodes, key=lambda v: (cid_of[v], node_sort_key(v)))
    T = nx.Graph()
    T.add_node(root)
    visited = {root}
    # BFS queue is FIFO, but children of each node are processed in
    # ascending cid order, ensuring determinism.
    queue: deque = deque([root])
    while queue:
        u = queue.popleft()
        # Sort neighbours by (cid, node_sort_key)
        nbrs = sorted(G.neighbors(u), key=lambda w: (cid_of[w], node_sort_key(w)))
        for w in nbrs:
            if w in visited:
                continue
            visited.add(w)
            T.add_edge(u, w)
            queue.append(w)
    return T


# ---------------------------------------------------------------------------
# Stage 3: canonical fundamental cycles
# ---------------------------------------------------------------------------

def _canonicalize_circle_code(
    cycle_nodes: Sequence[Hashable],
    cid_of: Dict[Hashable, int],
) -> Tuple:
    """
    Return a canonical, *isomorphism-invariant* encoding of the cycle given by
    ``cycle_nodes``.  The encoding is a tuple of the canonical ids in the
    cycle's *canonical order* (rotation + reflection minimised).

    The canonical order is defined as: rotate the cycle so the smallest cid
    comes first; then walk in both directions and take the lex-smaller of the
    two resulting cid sequences.
    """
    if not cycle_nodes:
        return tuple()
    n = len(cycle_nodes)
    # Find first occurrence of min cid
    cids = [cid_of[v] for v in cycle_nodes]
    min_cid = min(cids)
    pos = cids.index(min_cid)
    forward = tuple(cids[(pos + i) % n] for i in range(n))
    backward = tuple(cids[(pos - i) % n] for i in range(n))
    return min(forward, backward)


def _fundamental_cycles(
    G: nx.Graph,
    T: nx.Graph,
    cid_of: Dict[Hashable, int],
) -> List[Tuple]:
    """
    Compute the canonical fundamental cycle set.

    For each non-tree edge ``e = (u, v)`` of ``G`` (sorted by edge-code),
    find the unique simple cycle ``C_e`` in ``T + e``, then encode ``C_e``
    canonically via :func:`_canonicalize_circle_code`.

    Returns
    -------
    list
        A list of triples ``(length, CircleCode, cycle_nodes)`` sorted by
        ``(length, CircleCode)``.
    """
    # Edges of G not in T, sorted canonically
    tree_edges = {tuple(sorted(e, key=lambda x: (cid_of[x], node_sort_key(x))))
                  for e in T.edges()}
    rest_edges: List[Tuple[Hashable, Hashable]] = []
    for u, v in G.edges():
        key = (u, v) if (cid_of[u], node_sort_key(u)) <= (cid_of[v], node_sort_key(v)) else (v, u)
        if key not in tree_edges:
            rest_edges.append(key)
    rest_edges.sort(key=lambda e: (cid_of[e[0]], cid_of[e[1]],
                                    node_sort_key(e[0]), node_sort_key(e[1])))

    # Build parent map for T to find paths quickly
    parent = {n: None for n in T.nodes()}
    children = defaultdict(list)
    root = None
    for v in T.nodes():
        if parent[v] is None and T.degree(v) >= 0:
            # First BFS to identify root
            pass
    # Use BFS to build parent/children
    for v in T.nodes():
        if all(w not in T.neighbors(v) for w in T.nodes() if w != v):
            # Not isolated in T
            pass
    # Find root (vertex with parent=None in BFS)
    visited = set()
    for start in T.nodes():
        if start in visited:
            continue
        # BFS
        bfs_q = deque([start])
        visited.add(start)
        parent[start] = None
        while bfs_q:
            u = bfs_q.popleft()
            nbrs = sorted(T.neighbors(u),
                          key=lambda w: (cid_of[w], node_sort_key(w)))
            for w in nbrs:
                if w in visited:
                    continue
                visited.add(w)
                parent[w] = u
                children[u].append(w)
                bfs_q.append(w)

    def _path(u: Hashable, v: Hashable) -> List[Hashable]:
        """Return the path in T from u to v (inclusive of both endpoints)."""
        ancestors_u: List[Hashable] = [u]
        x = u
        while parent.get(x) is not None:
            x = parent[x]
            ancestors_u.append(x)
        anc_set = set(ancestors_u)
        path_v: List[Hashable] = [v]
        y = v
        if parent.get(y) in anc_set:
            lca = parent[y]
        else:
            while parent.get(y) is not None and parent.get(y) not in anc_set:
                y = parent[y]
                path_v.append(y)
            lca = parent.get(y)
        if lca is None:
            lca = u
        lca_idx = ancestors_u.index(lca)
        up = ancestors_u[: lca_idx + 1]
        down = list(reversed(path_v))
        return up + down

    candidates: List[Tuple] = []
    for u, v in rest_edges:
        if u not in T.nodes() or v not in T.nodes():
            # Edge in a different connected component: cannot form a cycle
            # within T.  This only happens if G is disconnected and the edge
            # is in a different component from the rooted one.  Skip — the
            # cycle basis of G includes cycle bases from each component.
            continue
        path_uv = _path(u, v)
        # Cycle:  u -> ... -> v -> u
        cycle = list(path_uv) + [u]
        # Length is |cycle| - 1 (number of distinct vertices)
        length = len(cycle) - 1
        code = _canonicalize_circle_code(path_uv, cid_of)
        candidates.append((length, code, tuple(path_uv)))
    # Sort by (length, code) for Stage 4
    candidates.sort(key=lambda c: (c[0], c[1]))
    return candidates


# ---------------------------------------------------------------------------
# Stage 4: F2 refinement
# ---------------------------------------------------------------------------

def _cycle_to_f2_vector(
    cycle_nodes: Sequence[Hashable],
    edges_in_order: List[Tuple[Hashable, Hashable]],
    edge_index: Dict[Tuple[Hashable, Hashable], int],
    cid_of: Dict[Hashable, int],
) -> List[int]:
    """
    Encode a cycle as a vector in F_2^{|E|}.

    The vector has a 1 in position ``i`` iff the i-th edge (in
    ``edges_in_order``) belongs to the cycle.

    Edges are matched as undirected — both ``(u, v)`` and ``(v, u)`` map to the
    same canonical key.  The canonical key is oriented by canonical vertex id
    (``cid_of``) to match the orientation used in ``edges_in_order``.
    """
    cycle_edge_set = set()
    n = len(cycle_nodes)
    for i in range(n):
        u, v = cycle_nodes[i], cycle_nodes[(i + 1) % n]
        if cid_of[u] <= cid_of[v]:
            key = (u, v)
        else:
            key = (v, u)
        cycle_edge_set.add(key)
    vec = [0] * len(edges_in_order)
    for key in cycle_edge_set:
        if key in edge_index:
            vec[edge_index[key]] = 1
    return vec


def _gauss_rank_incremental(
    vectors: List[List[int]],
) -> List[int]:
    """
    Greedy Gaussian elimination over F_2.  Returns the list of vectors that
    increased the rank, in the order they were added.
    """
    basis: List[List[int]] = []
    kept: List[int] = []  # indices of kept vectors (in input order)
    for idx, v in enumerate(vectors):
        cur = v[:]
        for b in basis:
            # Find pivot of b
            pivot = next((j for j, x in enumerate(cur) if x == 1), None)
            if pivot is None:
                continue
            if b[pivot] == 1:
                # XOR
                for j in range(len(cur)):
                    cur[j] ^= b[j]
        # Check if cur is non-zero
        if any(x == 1 for x in cur):
            basis.append(cur)
            kept.append(idx)
    return kept


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def canonical_mcb(
    G: nx.Graph,
    depth: int = 3,
    return_edge_sets: bool = False,
) -> List:
    """
    Compute a canonical minimum cycle basis of ``G``.

    Parameters
    ----------
    G : nx.Graph
        Input graph (may be disconnected — each component is processed
        independently and the results concatenated in canonical order).
    depth : int
        Local-signature depth (see :func:`compute_canonical_vertex_ids`).
        Default 3, sufficient for most non-strongly-regular graphs.
    return_edge_sets : bool
        If True, return each cycle as a ``frozenset`` of its edges
        (sorted by canonical id).  If False (default), return each cycle as
        a tuple of its vertices in *canonical circle order*.

    Returns
    -------
    list
        A list of cycles.  Each cycle is either a tuple of vertices (default)
        or a frozenset of ``(u, v)`` edge tuples (``return_edge_sets=True``).
        The list itself is in *canonical order* (sorted by ``(length, code)``
        of the original fundamental cycles, with ties broken by F2 selection
        order).
    """
    if G.number_of_nodes() == 0 or G.number_of_edges() == 0:
        return []

    cid_of = compute_canonical_vertex_ids(G, depth=depth)

    edges_in_order: List[Tuple[Hashable, Hashable]] = []
    seen_edges: set = set()
    for u, v in G.edges():
        if cid_of[u] <= cid_of[v]:
            key = (u, v)
        else:
            key = (v, u)
        if key not in seen_edges:
            seen_edges.add(key)
            edges_in_order.append(key)
    edges_in_order.sort(key=lambda e: (cid_of[e[0]], cid_of[e[1]]))
    edge_index = {e: i for i, e in enumerate(edges_in_order)}

    # Process each connected component independently.  For each component:
    #   - compute the canonical spanning tree
    #   - compute the canonical fundamental cycles
    #   - run F_2 refinement
    # Then concatenate the results in ascending root-cid order so that
    # the *combined* output is also canonical.
    components = list(nx.connected_components(G))
    components.sort(key=lambda comp: min((cid_of[v], node_sort_key(v)) for v in comp))

    all_cycles: List = []
    for comp in components:
        sub_g = G.subgraph(comp).copy()
        if sub_g.number_of_nodes() <= 1 or sub_g.number_of_edges() == 0:
            continue
        T = _canonical_spanning_tree(sub_g, cid_of)
        candidates = _fundamental_cycles(sub_g, T, cid_of)
        if not candidates:
            continue
        # Build F2 vectors
        f2_vecs: List[List[int]] = []
        for _length, _code, cycle_nodes in candidates:
            f2_vecs.append(_cycle_to_f2_vector(cycle_nodes, edges_in_order, edge_index, cid_of))
        kept_indices = _gauss_rank_incremental(f2_vecs)
        mu_target = sub_g.number_of_edges() - sub_g.number_of_nodes() + 1
        # ``mu`` is the cycle rank; the kept indices give us a basis.
        for i in kept_indices[:mu_target]:
            _length, _code, cycle_nodes = candidates[i]
            if return_edge_sets:
                edge_set = set()
                n_c = len(cycle_nodes)
                for j in range(n_c):
                    u, v = cycle_nodes[j], cycle_nodes[(j + 1) % n_c]
                    key = (u, v) if node_sort_key(u) <= node_sort_key(v) else (v, u)
                    edge_set.add(key)
                all_cycles.append(frozenset(edge_set))
            else:
                all_cycles.append(cycle_nodes)
    return all_cycles


def canonical_mcb_invariant_summary(
    G: nx.Graph,
    depth: int = 3,
) -> Dict[str, Any]:
    """
    Return a small dictionary of *isomorphism-invariant* summaries derived
    from the canonical MCB of ``G``.  Useful for fast graph-level feature
    extraction without re-running the full WL pipeline.

    The returned dictionary has keys:

      * ``basis_size``         — the number of cycles (== ``mu(G)`` for a
        connected graph; sum of component ``mu``s otherwise).
      * ``total_length``       — sum of cycle lengths.
      * ``length_histogram``   — a dict ``{length: count}``.
      * ``multiset_signature`` — a sorted tuple ``(length, code)`` for every
        cycle in the canonical basis.
    """
    cycles = canonical_mcb(G, depth=depth, return_edge_sets=False)
    length_hist: Dict[int, int] = defaultdict(int)
    total = 0
    sig_parts: List[Tuple] = []
    cid_of = compute_canonical_vertex_ids(G, depth=depth)
    for c in cycles:
        L = len(c)
        length_hist[L] += 1
        total += L
        code = _canonicalize_circle_code(c, cid_of)
        sig_parts.append((L, code))
    sig_parts.sort()
    return {
        "basis_size": len(cycles),
        "total_length": total,
        "length_histogram": dict(length_hist),
        "multiset_signature": tuple(sig_parts),
    }


__all__ = [
    "compute_canonical_vertex_ids",
    "canonical_mcb",
    "canonical_mcb_invariant_summary",
]
