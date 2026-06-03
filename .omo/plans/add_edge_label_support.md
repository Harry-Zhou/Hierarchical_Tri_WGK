# Edge-Label Support for Hierarchical Triangulated WL — Design Document

> **Objective**: Extend `cyclic_schema/hierarchical_triangulated_wl.py` to handle graphs with **both node labels and edge labels simultaneously**, producing both node and edge label histories across iterations.

---

## 1. Motivation

The current `hierarchical_triangular_wl()` implements a powerful node-label-only WL test through a CSG (Cyclic Schematic Graph) hierarchy. However, many real-world graphs (e.g., molecular graphs with bond types, social networks with relationship types) carry significant information on their **edges**. Edge labels provide crucial distinguishing power that node labels alone cannot capture:

- **Molecular graphs**: bond type (single/double/triple/aromatic) determines chemical properties
- **Knowledge graphs**: relation type matters as much as entity type
- **Social networks**: "friend" vs "colleague" vs "family" edges carry different semantics

The existing flat (non-hierarchical) edge WL implementations in `hierarchical_cycle_complex_wl_tools/` already demonstrate how to propagate edge labels through triangulated neighbors and cycle complexes. This design extends that approach to the **hierarchical triangulated** framework, leveraging the CSG hierarchy to capture both edge-labeled local structure and long-range cycle information.

---

## 2. Current State

### 2.1 Node-Only Hierarchical WL (`cyclic_schema/hierarchical_triangulated_wl.py`)

```
Function: hierarchical_triangular_wl(G1, G2, vlabel_np1, vlabel_np2, K, I)
Returns : (wl_np1, wl_np2)            # shape: (|V|, I+1) — node labels only
```

**Iteration loop** (K≥1):
```
For each iteration:
  1. Forward pass:  G → CSG¹ → CSG² → ... → CSGᵏ
     — Compute CSG label tuples from previous-layer node labels
     — Joint message passing on CSG layers
  2. Backward pass: CSGᵏ → ... → CSG¹ → G
     — Compute lower-graph label tuples prepending l_{G_i}(v) + sorted([higher_labels])
     — Joint message passing on lower graphs
  3. Store G node labels for this iteration
```

**K=0 path**: Delegates to `node_wl_test_triangulated_neighbors` (node-only TN-WL).

**Key characteristics**:
- No edge label parameters or outputs
- CSG layers carry no edge information — they encode cycle membership only
- Message passing aggregates neighbor node labels only
- Labels are stored in arrays indexed by sorted node order

### 2.2 Existing Flat Edge WL Implementations

| File | Function | Edge Role |
|------|----------|-----------|
| `edge_wl_via_triangulated_neighbors.py` | `edge_wl_test_triangulated_neighbors()` | Edge labels used in TN node propagation + updated each iteration |
| `edge_wl_via_cycle_complexes.py` | `edge_wl_test_cycle_complexes()` | Edge labels used in HCC node propagation + updated each iteration |
| `local_structure_edge_wl_tools.py` | `update_elabel()` | Combines (old_elabel, sorted_endpoint_vlabels) → compress |

**Common pattern**: Each iteration:
1. Propagate node labels **using edge labels** → new node labels
2. Update edge labels `(old_elabel, vlabel_lhs, vlabel_rhs)` → compress → new edge labels

**Edge label composition** (`collect_vlabels_edge`):
```
edge_label_at_iter_t = compress(elabel_at_t-1, vlabel_lhs_at_t-1, vlabel_rhs_at_t-1)
```
This means each edge label at iteration t encodes:
- The edge's own label history
- The labels of its two endpoint nodes at iteration t-1

### 2.3 Top-Level Consumer (`topo_aware_wl_test.py`)

Currently, `topo_aware_edge_wl_test()` calls the flat edge WL functions but **discards the edge label history** — only returning concatenated node label arrays. This means even in the flat case, edge labels are only used internally to improve node discrimination; their history is not consumed by upper layers.

---

## 3. Design Goals

1. **Handle both node + edge labels simultaneously**: The algorithm processes labeled edges alongside labeled nodes.
2. **Edge labels participate in node label propagation**: Edge labels provide distinguishing power that affects WL color refinement.
3. **Edge labels are updated each iteration**: Edge labels evolve as node labels change, maintaining the coupled dynamics from the flat WL.
4. **Both node and edge label histories are output**: Downstream consumers can use both for kernel construction.
5. **Backward compatible**: Existing node-only callers work unchanged — edge labels are optional.
6. **Consistent joint label assignment**: Labels are compressed jointly across both graphs (same aggregate → same label), preserving WL consistency.
7. **Handles non-consecutive/non-integer node IDs**: Unlike the flat WL which assumes 0..n-1 integer nodes, the hierarchical WL must work with arbitrary hashable node IDs.

---

## 4. Edge Label Data Model

### 4.1 Representation

```
elabel_dict: Dict[Tuple[Hashable, Hashable], int]
```

- Keys are **undirected edge tuples** `(u, v)` with `u <= v` (sorted, canonical form)
- Values are integer labels
- Edge labels must be provided for all edges of both graphs (or not at all)
- The dictionary also serves as the edge list (keys = edges)

This dict-based representation works with arbitrary node IDs (strings, integers, tuples, etc.) — unlike array-based representations that require 0..n-1 consecutive integers.

### 4.2 Initialization

Edge labels can come from:
1. **Dataset attributes**: Molecular bond types, relation types in knowledge graphs
2. **Topological features**: Common-neighbor count, Jaccard coefficient, etc. (as used in the existing self-tests)
3. **Uniform labels**: All 1s if all edges are equivalent

### 4.3 Internal Evolution

Each iteration, edge labels are updated via:
```
new_elabel = compress(old_elabel, sorted(vlabel[u], vlabel[v]))
```

Where:
- `old_elabel` is the edge's label from the previous iteration
- `vlabel[u]`, `vlabel[v]` are the current node labels of the endpoints
- `compress` assigns the same new label to all edges with the same `(old_elabel, sorted_endpoint_labels)` tuple across both graphs

This is implemented by `update_elabel()` from `local_structure_edge_wl_tools.py`, adapted for dict-based node labels (see §5.3).

---

## 5. Algorithm Design

### 5.1 High-Level Flow (K≥1)

```
Input: G1, G2, vlabel_np1/2, elabel_dict1/2 (optional), K, I
Output: (vwl_np1, vwl_np2, ewl_np1, ewl_np2)   [if edge labels provided]
        (wl_np1, wl_np2)                        [if no edge labels]

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
For each iteration t = 1..I:
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 1: Forward pass (G → CSG¹ → ... → CSGᵏ)                 │
  │                                                               │
  │  Same as current node-only implemention:                      │
  │   - Compute CSG label tuples from current G node labels       │
  │   - Joint message passing on each CSG layer                   │
  │   - NO direct edge label involvement                          │
  │   (Edge labels affect this indirectly through G's enhanced    │
  │    node labels from the previous iteration's backward pass)   │
  └─────────────────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 2: Backward pass (CSGᵏ → ... → CSG¹ → G)                │
  │                                                               │
  │  For each layer k = K down to 1:                              │
  │    - Compute lower-graph label tuples:                        │
  │        if k == 1 (G level) and edge labels exist:             │
  │          label_tuple(v) = (l_G(v), edge_agg(v),               │
  │                           sorted([higher_labels]))            │
  │        else (CSG level or no edge labels):                    │
  │          label_tuple(v) = (l_G(v),                            │
  │                           sorted([higher_labels]))            │
  │    - Joint message passing on lower graph                     │
  └─────────────────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 3: Edge label update                                     │
  │                                                               │
  │  If edge labels exist:                                        │
  │    - Build vlabel dict from updated G node labels             │
  │    - Call update_elabel(vlabel_dict, elabel_dict)             │
  │    - Store compressed edge labels for this iteration          │
  └─────────────────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 4: Store G node labels for iteration t                   │
  │                                                               │
  │  vwl_np1[:, t] = labels from backward pass's G level          │
  │  vwl_np2[:, t] = labels from backward pass's G level          │
  └─────────────────────────────────────────────────────────────────┘
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

### 5.2 Edge Aggregate Computation (`edge_agg(v)`)

The core innovation is how edge labels are incorporated into G node label tuples. Following the design of `propagate_via_triangulated_neighboring_edges`, the edge label aggregate for node `v` is computed from its **triangulated neighborhood components**:

```
edge_agg(v) = (
    (elabel(v, u₁₁), elabel(v, u₁₂), ...),   # sorted edge labels for component 1
    (elabel(v, u₂₁), elabel(v, u₂₂), ...),   # sorted edge labels for component 2
    ...
)
```

Where:
- N(v) is split into connected components of the induced subgraph `G[N(v)]`
- For each component, we collect edge labels of edges `(v, u)` for each `u` in the component
- Each component's edge labels are sorted internally
- Components are sorted lexicographically by their edge label tuple

**Rationale**: This matches the structure of the existing node-level aggregation where each component contributes a sorted tuple of neighbor node labels. Here, each component contributes a sorted tuple of incident **edge** labels. This ensures that the edge information is structured in the same way as the neighborhood decomposition — preserving the topology-aware grouping.

**Implementation sketch**:
```python
def _edge_aggregate_for_node(v, G, elabel_dict, neighbor_components):
    """
    Compute edge label aggregate for node v.
    
    Parameters
    ----------
    v : Hashable
        Node ID
    G : nx.Graph
        The graph
    elabel_dict : dict
        {edge_tuple: int} — current edge labels
    neighbor_components : tuple of tuples
        Precomputed connected components of N(v)
    
    Returns
    -------
    tuple
        Edge label aggregate, empty tuple if no neighbors
    """
    components = neighbor_components.get(v, ())
    if not components:
        return ()
    
    comp_aggs = []
    for comp in components:
        # Collect edge labels of edges from v to each node in this component
        comp_edge_labels = tuple(sorted(
            elabel_dict.get(tuple(sorted((v, u))), 0)
            for u in comp
        ))
        comp_aggs.append(comp_edge_labels)
    
    comp_aggs.sort()
    return tuple(comp_aggs)
```

### 5.3 Dict-Based Edge Label Update

Since the hierarchical WL supports arbitrary node IDs, the standard `update_elabel` (which uses `vlabel_np[e[0]]` array indexing) must be adapted:

```python
def _update_elabel_from_dict(vlabel_dict1, vlabel_dict2, elabel_dict1, elabel_dict2):
    """
    Dict-based edge label update for non-consecutive node IDs.
    
    For each edge (u,v), the new label is:
        label = compress(old_elabel, sorted(vlabel_dict[u], vlabel_dict[v]))
    
    Returns new elabel dicts with updated labels.
    """
    # Collect edge labels with their current (old_elabel, sorted_endpoint_labels)
    elabel_collection1 = {}
    for e, l in elabel_dict1.items():
        vl0, vl1 = vlabel_dict1[e[0]], vlabel_dict1[e[1]]
        if vl0 <= vl1:
            elabel_collection1[e] = (l, vl0, vl1)
        else:
            elabel_collection1[e] = (l, vl1, vl0)
    
    elabel_collection2 = {}
    for e, l in elabel_dict2.items():
        vl0, vl1 = vlabel_dict2[e[0]], vlabel_dict2[e[1]]
        if vl0 <= vl1:
            elabel_collection2[e] = (l, vl0, vl1)
        else:
            elabel_collection2[e] = (l, vl1, vl0)
    
    # Joint compression across both graphs
    all_labels = list(elabel_collection1.values()) + list(elabel_collection2.values())
    all_labels = sorted(set(all_labels))
    next_label = max(lc[0] for lc in all_labels) + 1
    
    compressed1 = {}
    for e, tup in elabel_collection1.items():
        compressed1[e] = all_labels.index(tup) + next_label
    
    compressed2 = {}
    for e, tup in elabel_collection2.items():
        compressed2[e] = all_labels.index(tup) + next_label
    
    return compressed1, compressed2
```

### 5.4 K=0 Case: Triangulated Neighbors Edge WL

For K=0 with edge labels, we use `edge_wl_test_triangulated_neighbors` from the existing flat module, which requires:
- Remapped consecutive node IDs (0..n-1)
- Edge label arrays indexed by edge order

**Conversion**: `elabel_dict` → `edges_list` + `elabel_np`

```python
edges1 = list(elabel_dict1.keys())
elabel_np1 = np.array([elabel_dict1[e] for e in edges1], dtype=np.int32)
```

The function returns `(vwl_np1, vwl_np2, ewl_np1, ewl_np2)` — both node and edge label histories.

### 5.5 Cycle Complex Edge WL (HCC Layer)

For K≥1, the HCC-based propagation in `node_wl_via_cycle_complexes.py` uses node labels only. The cycle complexes are constructed from graph topology (BFS on cycle graph), and the node label propagation through these complexes uses node labels.

With edge labels, the HCC propagation should also incorporate edge labels, following the same pattern as `edge_wl_via_cycle_complexes.py`. This means:
- When building the node label collection from cycle complex neighbors, instead of collecting node labels of neighbor nodes, collect **edge labels** of edges within each cycle complex.
- Update edge labels after each HCC iteration.

**However**, since `hierarchical_triangular_wl` uses `compute_final_label_tuples` + `forward_message_passing_both` for CSG layers (not the `propagate_via_cycle_complexes` from the flat module), the edge label integration for cycle complexes follows the same enhanced-backward-pass approach: edge-enhanced node label tuples propagate edge information through the standard message passing mechanism.

For this initial design, edge labels in the CSG hierarchy are handled exclusively through:
- Enhanced G-level backward label tuples (incorporating edge aggregates)
- Edge label refresh after each complete iteration

CSG layer message passing remains unchanged — it operates on node label tuples that already encode edge information (from the G-level enhancement).

---

## 6. Function Interface

### 6.1 New Signature

```python
def hierarchical_triangular_wl(
    G1: nx.Graph,
    G2: nx.Graph,
    vlabel_np1: Any,
    vlabel_np2: Any,
    elabel_dict1: Optional[Dict] = None,
    elabel_dict2: Optional[Dict] = None,
    K: int = 1,
    I: int = 5,
) -> Union[
    Tuple[np.ndarray, np.ndarray],                                # node-only
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],        # with edge labels
]:
```

### 6.2 Parameter Semantics

| Parameter | Type | Description |
|-----------|------|-------------|
| `G1, G2` | `nx.Graph` | Input graphs (unchanged) |
| `vlabel_np1/2` | `ndarray` or `dict` | Initial node labels (unchanged) |
| `elabel_dict1/2` | `Optional[Dict]` | **NEW**: `{edge_tuple: int}` — initial edge labels. `None` = no edge labels |
| `K` | `int` | CSG layers (unchanged) |
| `I` | `int` | Iterations (unchanged) |

### 6.3 Return Values

| Return | Shape | Description |
|--------|-------|-------------|
| `vwl_np1, vwl_np2` (with or w/o edges) | `(|V|, I+1)` | Node label history (same as current `wl_np1/2`) |
| `ewl_np1, ewl_np2` (when edge labels provided) | `(|E|, I+1)` | **NEW**: Edge label history |

**Backward compatibility**: When `elabel_dict1` is `None`, returns `(wl_np1, wl_np2)` — identical to current behavior.

### 6.4 Error Conditions

- `elabel_dict1 is None` xor `elabel_dict2 is None` → `ValueError` (both or neither)
- `elabel_dict1` has edges not in `G1` → `ValueError`
- `elabel_dict1` missing edges from `G1` → accepted (unlabeled edges get label 0)
- Negative values for `K` or `I` → same as current (`ValueError`)

---

## 7. Precomputation Changes

### 7.1 New Precomputation: Neighbor Components for G

The `_precompute_neighbor_components` already computes the connected components of N(v) for each node. This is used in `forward_aggregate` for message passing.

For edge label aggregation, the **same** neighbor component structure can be reused — edge labels within each component are collected.

**No additional precomputation needed** for the edge aggregate — we reuse the existing `nc_input1`, `nc_input2`.

### 7.2 New Precomputation: Triangulated Neighbors for K=0 Edge WL

For K=0 with edge labels, the triangulated neighbor structure is built using `_build_vtx_triangulated_neighbors` (same as current, but now used with `edge_wl_test_triangulated_neighbors` instead of the node-only variant).

---

## 8. Detailed Modification Plan

### 8.1 Files to Modify

| File | Changes |
|------|---------|
| `cyclic_schema/hierarchical_triangulated_wl.py` | Main implementation changes (detailed below) |

### 8.2 Changes to `hierarchical_triangular_wl()`

#### 8.2.1 Signature (lines 463-470)
- Add `elabel_dict1=None, elabel_dict2=None` parameters
- Update docstring to document edge label behavior
- Add parameter validation for edge label consistency

#### 8.2.2 Validation (after existing K/I validation)
```python
if (elabel_dict1 is None) != (elabel_dict2 is None):
    raise ValueError("Either both or neither elabel_dict must be provided")
if elabel_dict1 is not None and elabel_dict2 is not None:
    for e in elabel_dict1:
        if e not in G1.edges() and (e[1], e[0]) not in G1.edges():
            raise ValueError(f"Edge {e} in elabel_dict1 not in G1")
    # ... same for G2 ...
```

#### 8.2.3 Edge Label Storage

Add edge label history arrays when edge labels are present:
```python
has_edge_labels = elabel_dict1 is not None
if has_edge_labels:
    edges1 = sorted(elabel_dict1.keys())
    edges2 = sorted(elabel_dict2.keys())
    ewl_np1 = np.zeros((len(edges1), I + 1), dtype=np.int32)
    ewl_np2 = np.zeros((len(edges2), I + 1), dtype=np.int32)
    # Store initial edge labels
    for i, e in enumerate(edges1):
        ewl_np1[i, 0] = elabel_dict1[e]
    for i, e in enumerate(edges2):
        ewl_np2[i, 0] = elabel_dict2[e]
```

#### 8.2.4 Modified K=0 Path

```python
if K == 0:
    if has_edge_labels:
        # Use edge-aware triangulated neighbors WL
        from edge_wl_via_triangulated_neighbors import edge_wl_test_triangulated_neighbors
        # ... build vtx_tri (same as current) ...
        # Convert elabel_dict to arrays
        edges1 = sorted(elabel_dict1.keys())
        edges2 = sorted(elabel_dict2.keys())
        elabel_arr1 = np.array([elabel_dict1[e] for e in edges1], dtype=np.int32)
        elabel_arr2 = np.array([elabel_dict2[e] for e in edges2], dtype=np.int32)
        # ... remap nodes (same as current) ...
        vwl_tri1, vwl_tri2, ewl_tri1, ewl_tri2 = edge_wl_test_triangulated_neighbors(
            vtx_tri1, vtx_tri2, vlabel_arr1, vlabel_arr2,
            edges1, edges2, elabel_arr1, elabel_arr2, I)
        vwl_np1[:, :] = vwl_tri1[:, :]
        vwl_np2[:, :] = vwl_tri2[:, :]
        ewl_np1[:, :] = ewl_tri1[:, :]
        ewl_np2[:, :] = ewl_tri2[:, :]
        return vwl_np1, vwl_np2, ewl_np1, ewl_np2
    else:
        # Same as current: use node_wl_test_triangulated_neighbors
        ...
```

#### 8.2.5 Modified K≥1 Backward Pass

When computing lower label tuples for G (step=1 in backward pass) and edge labels exist:

```python
def _compute_lower_label_tuples_with_edges(
    lower_G, lower_labels, higher_labels, lower_to_higher,
    elabel_dict, neighbor_components,
):
    result = {}
    for v in lower_G.nodes():
        higher_nodes = lower_to_higher.get(v, ())
        l_v = int(lower_labels.get(v, 0))
        
        # Edge aggregate
        edge_agg = _edge_aggregate_for_node(v, lower_G, elabel_dict, neighbor_components)
        
        higher_lbls = sorted(int(higher_labels.get(u, 0)) for u in higher_nodes) if higher_nodes else ()
        
        if not higher_nodes and not edge_agg:
            result[v] = (l_v,)
        elif not higher_nodes:
            result[v] = (l_v, edge_agg)
        elif not edge_agg:
            result[v] = (l_v,) + tuple(higher_lbls)
        else:
            result[v] = (l_v, edge_agg) + tuple(higher_lbls)
    return result
```

Then in the backward pass loop (step == 1 and has_edge_labels):
```python
if step == 1 and has_edge_labels:
    label_tuples1 = _compute_lower_label_tuples_with_edges(
        G1, current_labels1, higher_labels1,
        mappings1[0], current_elabel1, nc_input1)
    label_tuples2 = _compute_lower_label_tuples_with_edges(
        G2, current_labels2, higher_labels2,
        mappings2[0], current_elabel2, nc_input2)
else:
    label_tuples1 = _compute_lower_label_tuples(...)  # original
    label_tuples2 = _compute_lower_label_tuples(...)
```

#### 8.2.6 Edge Label Update After Each Iteration

After the backward pass completes:
```python
if has_edge_labels:
    # Build vlabel dicts from updated G node labels
    vlabel_dict1 = {v: int(new_labels_G1[v]) for v in G1.nodes()}
    vlabel_dict2 = {v: int(new_labels_G2[v]) for v in G2.nodes()}
    
    # Update edge labels
    current_elabel1, current_elabel2 = _update_elabel_from_dict(
        vlabel_dict1, vlabel_dict2, current_elabel1, current_elabel2)
    
    # Store edge labels for this iteration
    for i, e in enumerate(edges1_sorted):
        ewl_np1[i, iteration] = current_elabel1[e]
    for i, e in enumerate(edges2_sorted):
        ewl_np2[i, iteration] = current_elabel2[e]
```

### 8.3 Modularization

The following new helper functions are added to `hierarchical_triangulated_wl.py`:

| Function | Purpose |
|----------|---------|
| `_edge_aggregate_for_node(v, G, elabel_dict, nc)` | Compute edge label aggregate for node v |
| `_compute_lower_label_tuples_with_edges(...)` | Enhanced backward label tuples with edge info |
| `_update_elabel_from_dict(vlabel_dict1/2, elabel_dict1/2)` | Dict-based edge label update |
| `_validate_edge_labels(G1, elabel_dict1, G2, elabel_dict2)` | Validate edge label inputs |

---

## 9. Worked Example

### 9.1 Simple Test: Triangle with Colored Edges

```
Graph: Triangle (3 nodes, 3 edges)

Nodes:  all have vlabel = 1
Edges:  (0,1) elabel = 10, (1,2) elabel = 20, (0,2) elabel = 10

K=1, I=2 (1 CSG layer, 2 iterations)
```

**Iteration 1**:

1. **Forward**: G labels all 1 → CSG cycle basis node gets label tuple `(1, 1, 1)` (three node labels, canonicalized). CSG message passing doesn't differentiate.

2. **Backward** (G level, with edge aggregates):
   - Node 0: N(0) = {1, 2}, connected component = (1, 2), edges (0,1)=10, (0,2)=10
     - `edge_agg = ((10, 10),)` — one component with two edge labels
     - label_tuple = `(1, (10, 10), [CSG_labels])`
   - Node 1: N(1) = {0, 2}, edges (1,0)=10, (1,2)=20
     - `edge_agg = ((10, 20),)`
     - label_tuple = `(1, (10, 20), [CSG_labels])`
   - Node 2: N(2) = {0, 1}, edges (2,0)=10, (2,1)=20
     - `edge_agg = ((10, 20),)`
     - label_tuple = `(1, (10, 20), [CSG_labels])`

3. **Message passing on G**: Nodes 0 and 1/2 have different label tuples → they get different new labels.

4. **Edge label update**: 
   - (0,1): old=10, vlabels=(l(0), l(1))
   - (1,2): old=20, vlabels=(l(1), l(2))
   - (0,2): old=10, vlabels=(l(0), l(2))
   - After compression: edges with same (old_elabel, sorted_endpoint_labels) get same new label

**Without edge labels**: All 3 nodes would have identical label tuples → no differentiation → indistinguishable triangle.

**With edge labels**: Node 0 is distinguished from nodes 1 and 2 (different incident edge patterns).

### 9.2 Expected Test Outcomes

| Test Case | Expected Result |
|-----------|----------------|
| Identical graphs, identical edge labels | wl arrays match exactly (WL consistency) |
| Graphs differ in edge labels only | WL distinguishes them |
| K=0 edge WL matches existing `edge_wl_test_triangulated_neighbors` | Numerically identical |
| K≥1 with uniform edge labels | Same as node-only (no added distinguishing power) |

---

## 10. Verification Strategy

### 10.1 WL Consistency Checks

1. **Joint compression correctness**: Identical label tuples across both graphs get identical compressed labels.
2. **Iteration monotonicity**: Labels are non-decreasing within each iteration.
3. **Determinism**: Multiple runs with same inputs produce identical outputs.

### 10.2 Isomorphism Tests

| Test | Description |
|------|-------------|
| Identical graphs + labels | `vwl_np1 == vwl_np2` and `ewl_np1 == ewl_np2` |
| Isomorphic relabeling | Multiset of labels per column is the same (WL invariant property) |
| Edge-only differing graphs | WL distinguishes graphs that differ only in edge labels |
| K=0 vs flat edge WL | Outputs match `edge_wl_test_triangulated_neighbors` |
| K≥1 with uniform edges | Outputs match node-only `hierarchical_triangular_wl` |
| Empty edges | Degenerate case with no edges → identical to node-only |

### 10.3 Edge Label Validity

- All edge labels in output are non-negative integers
- Edge label arrays have correct shape `(|E|, I+1)`
- Initial edge labels match input
- Edge label progression is non-decreasing

---

## 11. Implementation Order

| Step | Description | Dependencies |
|------|-------------|-------------|
| 1 | Add `_edge_aggregate_for_node()` helper | None |
| 2 | Add `_update_elabel_from_dict()` helper | None |
| 3 | Add `_compute_lower_label_tuples_with_edges()` helper | Step 1 |
| 4 | Modify `hierarchical_triangular_wl()` signature + validation | None |
| 5 | Implement K=0 edge WL path | None (uses existing `edge_wl_test_triangulated_neighbors`) |
| 6 | Implement K≥1 backward pass edge integration | Steps 1-3 |
| 7 | Implement edge label update + storage each iteration | Step 2 |
| 8 | Add edge label return + backward-compatible wrapper | Step 4 |
| 9 | Add self-tests | Steps 1-8 |
| 10 | Integration test with `topo_wasserstein_graph_kernel.py` | Step 9 |

---

## 12. Open Questions

1. **Should edge label history be returned separately or merged into a single multi-dimensional array?** The current flat WL returns separate `vwl_np` and `ewl_np`. Following this convention for consistency.

2. **Should edge labels participate in CSG-layer message passing, or only at the G level?** For this design, only at G level. CSG layers encode cycle membership structure, and edges in CSGs represent cycle-sharing relationships, not the original graph edges. Edge labels of the original graph don't map naturally to CSG edges.

3. **How should the edge label update handle edges that appear in multiple cycles?** Edge labels are properties of the original graph edges, and each edge is updated once per iteration regardless of how many cycles it belongs to. This is consistent with the flat edge WL.

4. **Should we add a pre-processing step to incorporate edge labels into G's node labels before the first forward pass?** Not necessary — the first backward pass will incorporate edge aggregates. However, adding one such step before iteration 1 would make edge labels effective from the very first forward pass. Design decision: skip pre-processing for simplicity; edge labels take effect from iteration 1's backward pass and affect the forward pass from iteration 2 onward.
