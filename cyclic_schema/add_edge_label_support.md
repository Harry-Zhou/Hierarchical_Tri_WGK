# Edge-Label Support for Hierarchical Triangulated WL

## Summary

This design extends `hierarchical_triangular_wl()` to process graphs with **both node labels and edge labels simultaneously** through the CSG hierarchy. The core mechanism is **edge-to-node fusion**: each G node's label tuple is augmented with an *edge context* that summarizes labels of incident edges (structured by triangulated neighborhood components). Once fused into G's node label tuples, edge information propagates through the standard hierarchical message passing:

- **Forward** (G → CSG¹ → ... → CSGᵏ): CSG cycle labels are computed from edge-enhanced G node labels → cycle labels naturally encode edge info.
- **Backward** (CSGᵏ → ... → CSG¹ → G): G's label tuples incorporate edge context + higher-layer labels → message passing spreads edge-aware signatures.
- **Edge refresh**: After each full forward+backward cycle, edge labels are refreshed via `compress(old_elabel, sorted(vlabel[u], vlabel[v]))`.

The result: **two label histories** — node labels of shape `(|V|, I+1)` and edge labels of shape `(|E|, I+1)` — with edge information flowing through all K layers of the hierarchy.

---

## 1. Core Idea: Edge-to-Node Fusion

Edge labels are not a separate stream that travels up the CSG hierarchy. Instead, they are **fused into G node label tuples** before node labels enter hierarchical message passing.

### How fusion works

For each G node `v`, its *edge context* `ec(v)` is computed from incident edge labels, grouped by triangulated neighborhood components (the same components used in `forward_aggregate`):

```
ec(v) = (
    (elabel(v, u₁₁), elabel(v, u₁₂), ...),   # sorted edge labels for component C₁
    (elabel(v, u₂₁), elabel(v, u₂₂), ...),   # sorted edge labels for component C₂
    ...
)
```

Where:
- N(v) is split into connected components of the induced subgraph G[N(v)] (precomputed by `_precompute_neighbor_components`).
- For component C = {u₁, u₂, ...}, collect edge labels of edges `(v, uᵢ)`.
- Each component's edge labels are sorted internally; components are sorted lexicographically.

**Then**, when G's label tuples are computed in the backward pass, each node's tuple becomes:

```
label_tuple(v) = (l_G(v), ec(v)) + sorted([higher_layer_labels])
# Node-only (baseline):  (l_G(v),) + sorted([higher_layer_labels])
```

This `ec(v)` component becomes part of v's identity during message passing — when a neighbor `w` aggregates v's label tuple, `ec(v)` is included in the aggregate. This is how edge information spreads through the graph via standard WL message passing.

### Why this works in the hierarchy

Once edge context is part of G's node label tuples:

1. **Forward pass** (G → CSG¹): `compute_final_label_tuples()` composes CSG cycle labels from G node labels. Since G node labels now carry edge context `ec(v)`, the cycle labels `(l_G(n₁), ec(n₁), l_G(n₂), ec(n₂), ...)` encode edge information.

2. **CSG message passing**: Standard `forward_aggregate` collects label tuples from neighbors. These tuples already contain edge-derived information, which propagates normally.

3. **Back to G**: After backward pass, G nodes receive new scalar labels reflecting the edge-enhanced signatures.

4. **Edge refresh**: `update_elabel()` recomputes each edge label as `compress(old_elabel, sorted(vlabel[u], vlabel[v]))`, using the freshly updated G node labels.

---

## 2. Algorithm Flow

### 2.1 K=0: Triangulated Neighbors (no CSG)

Uses the existing `edge_wl_via_triangulated_neighbors.py` — no changes needed:

```
For each iteration:
  1. vlabel = propagate_via_triangulated_neighboring_edges(vtx_tri, vlabel, elabel_dict)
     (edge labels used in triangulated component aggregation)
  2. vlabel = compress_and_relabel(vlabel) across both graphs
  3. elabel_dict = update_elabel(elabel_dict, vlabel)
```

Input: `elabel_dict1/2` → convert to `edges1/2` + `elabel_arr1/2` for the existing API.
Output: `(vwl_np1, vwl_np2, ewl_np1, ewl_np2)`.

### 2.2 K≥1: Hierarchical CSG Message Passing

```
Input:  G1, G2, vlabel_np1/2, elabel_dict1/2 (optional), K, I
Output: (vwl_np1, vwl_np2, ewl_np1, ewl_np2)   [if edge labels]
        (wl_np1, wl_np2)                         [if no edge labels]

┌─────────────────────────────────────────────────────────────────────┐
│ Pre-processing (only when edge labels provided):                    │
│   - Compute initial edge context ec(v) for each G node:            │
│     ec(v) = component-grouped incident edge labels                 │
│   - Store ec(v) for later use in backward label tuples             │
│                                                                     │
│   Effect: edge labels are "pre-fused" into G's label representation │
│   so that the FIRST forward pass already sees edge-aware node info │
└─────────────────────────────────────────────────────────────────────┘

For each iteration t = 1..I:
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 1: FORWARD pass (G → CSG¹ → CSG² → ... → CSGᵏ)           │
  │                                                                 │
  │  for each layer k = 0..K-1:                                    │
  │    a) Compute CSG label tuples from previous-layer labels       │
  │       (if edge labels exist, previous-layer labels are already  │
  │        edge-aware from pre-processing or prior backward pass)   │
  │    b) Joint message passing on CSGᵏ (same as current)          │
  │                                                                 │
  │  Edge labels do NOT directly participate in CSG computation.     │
  │  They affect it INDIRECTLY through G's edge-enhanced node labels │
  └─────────────────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 2: BACKWARD pass (CSGᵏ → ... → CSG¹ → G)                 │
  │                                                                 │
  │  for each layer k = K down to 1:                               │
  │    a) Compute lower-graph label tuples:                        │
  │                                                                 │
  │       if k == 1 (G level) AND edge labels provided:            │
  │         label_tuple(v) = (l_G(v), ec(v),                       │
  │                           sorted([higher_labels]))             │
  │       else (CSG level or no edge labels):                      │
  │         label_tuple(v) = (l_G(v),                              │
  │                           sorted([higher_labels]))             │
  │                                                                 │
  │    b) Joint message passing on lower graph                     │
  │       (edge context ec(v) is aggregated from neighbors         │
  │        → edge info spreads through standard WL propagation)    │
  └─────────────────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 3: EDGE LABEL REFRESH                                      │
  │                                                                 │
  │  After backward pass gives G new scalar labels:                 │
  │    vlabel_dict[v] = new_G_label[v]                              │
  │    for each edge (u,v):                                         │
  │      new_elabel = compress(old_elabel,                          │
  │                            sorted(vlabel_dict[u],               │
  │                                   vlabel_dict[v]))              │
  │    Store compressed edge labels for iteration t                 │
  │                                                                 │
  │  Also RECOMPUTE ec(v) using refreshed edge labels:              │
  │    ec(v) = component-grouped incident edge labels               │
  │    (ready for next iteration's backward pass)                   │
  └─────────────────────────────────────────────────────────────────┘
       │
       ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ Step 4: STORE G node labels for iteration t                    │
  │  vwl_np1[:, t] = labels from backward pass's G level           │
  └─────────────────────────────────────────────────────────────────┘
```

### Key Design Points

1. **Edge-to-node fusion happens at two points**:
   - **Pre-processing** (before iteration 1): Initial `ec(v)` computed from input edge labels.
   - **After edge refresh** (each iteration): `ec(v)` recomputed from updated edge labels.

2. **Forward pass is unchanged** — but it processes edge-enhanced node labels.

3. **Backward pass G-level is modified** — label tuples include `ec(v)`.

4. **CSG layers are completely unchanged** — they operate on node label tuples that already carry edge information from lower layers.

5. **Edge labels are refreshed after each full cycle** — capturing the coupled node-edge dynamics from the flat edge WL.

---

## 3. Concrete Example: Triangle with Colored Edges

```
Graph: triangle (3 nodes, 3 edges)
Nodes:  all vlabel = 1
Edges:  (0,1)=10, (1,2)=20, (0,2)=10
K=1, I=2
```

**Pre-processing**: Compute ec(v) for each node.

| Node | Neighbor Components | ec(v) |
|------|-------------------|-------|
| 0 | {1, 2} (one component) | `((10, 10),)` |
| 1 | {0, 2} (one component) | `((10, 20),)` |
| 2 | {0, 1} (one component) | `((10, 20),)` |

**Iteration 1**:

1. **Forward** (G → CSG¹):
   - G labels are `{(0: 1), (1: 1), (2: 1)}` — still scalars from init, edge-unaware
   - CSG cycle label: `(1, 1, 1)` → canonicalized → `(1, 1, 1)`
   - CSG message passing → all CSG nodes get same label → no discrimination

2. **Backward** (CSG¹ → G):
   - G label tuples with edge context:
     - Node 0: `(1, (10, 10), [CSG_labels])`  ← distinguished from node 1,2!
     - Node 1: `(1, (10, 20), [CSG_labels])`
     - Node 2: `(1, (10, 20), [CSG_labels])`
   - Message passing on G: nodes 0 vs 1/2 have different label tuples → different new labels
   - → e.g., node 0 gets label 100, nodes 1/2 get label 101

3. **Edge refresh**:
   - `(0,1)`: old=10, vlabels=(100, 101) → new compressed label
   - `(1,2)`: old=20, vlabels=(101, 101) → different from (0,1)
   - `(0,2)`: old=10, vlabels=(100, 101) → same as (0,1)
   - Recompute ec(0)=((10,10),), ec(1)=((10,20),), ec(2)=((10,20),) using refreshed edge labels

4. **Store iteration 1 G labels**: `[100, 101, 101]`

**Iteration 2**:

1. **Forward** (G → CSG¹):
   - G labels = `{0: 100, 1: 101, 2: 101}` — now edge-aware!
   - CSG cycle label: `(100, 101, 101)` → different from `(1, 1, 1)` in iter 1
   - CSG message passing produces distinct labels

2. **Backward** with refreshed ec(v) — further refinement...

3. **Result**: triangle with 10/20/10 edges is distinguished from triangle with all-10 edges (which would have all nodes get same labels throughout).

**Without edge labels**: All 3 nodes identical → no discrimination → triangle indistinguishable from any other triangle.

---

## 4. Function Interface

```python
def hierarchical_triangular_wl(
    G1: nx.Graph,
    G2: nx.Graph,
    vlabel_np1: Any,
    vlabel_np2: Any,
    elabel_dict1: Optional[Dict] = None,   # NEW
    elabel_dict2: Optional[Dict] = None,   # NEW
    K: int = 1,
    I: int = 5,
) -> Union[
    Tuple[np.ndarray, np.ndarray],                                    # node-only
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],            # with edge labels
]:
```

### Parameters

| Param | Type | Description |
|-------|------|-------------|
| `elabel_dict1/2` | `Optional[Dict[Tuple, int]]` | NEW. `{(u,v): label}` with `u<=v`. Both or neither. |

### Returns

| Return | Shape | Description |
|--------|-------|-------------|
| `vwl_np1/2` | `(|V|, I+1)` | Node label history (same as current `wl_np1/2`) |
| `ewl_np1/2` | `(|E|, I+1)` | **NEW**. Edge label history. Only when `elabel_dict` provided. |

### Edge Label Dict Format

```python
# Canonical undirected edge key: smaller node first
elabel_dict1 = {
    (0, 1): 10,   # edge (0,1) has label 10
    (1, 2): 20,   # edge (1,2) has label 20
    (0, 2): 10,   # edge (0,2) has label 10
}
```

Must use sorted node order for keys: `(u, v)` where `u <= v`.

---

## 5. New Helper Functions

All added to `hierarchical_triangular_wl.py`:

### 5.1 `_compute_edge_context(v, G, elabel_dict, nc)`

```
Compute ec(v) = component-grouped, sorted incident edge labels.

Input:
  v: node ID
  G: nx.Graph
  elabel_dict: {edge: int}
  nc: neighbor components from _precompute_neighbor_components

Output: tuple of tuples, e.g., ((10, 10), (20,)) or () if no neighbors
```

### 5.2 `_compute_lower_label_tuples_with_edges(...)`

```
Enhanced version of _compute_lower_label_tuples.
For G nodes (when edge labels exist):
  label_tuple(v) = (l_G(v), ec(v), sorted([higher_labels]))
For CSG nodes or when no edge labels: same as original.
```

### 5.3 `_update_elabel_from_dict(vlabel_dict1, vlabel_dict2, elabel_dict1, elabel_dict2)`

```
Dict-based edge label update (works with arbitrary node IDs).

For each edge (u,v):
  new_label = compress(old_elabel, sorted(vlabel[u], vlabel[v]))
Joint compression across both graphs.
Returns (new_elabel_dict1, new_elabel_dict2).
```

### 5.4 `_validate_edge_labels(G1, elabel_dict1, G2, elabel_dict2)`

```
Validate edge label dicts:
  - Both or neither provided
  - All edge keys exist in the corresponding graph
  - All keys are sorted (u <= v)
```

---

## 6. Implementation Plan

### Step 1: Add helper functions (no behavioral changes)

- `_compute_edge_context`
- `_update_elabel_from_dict`
- `_validate_edge_labels`

### Step 2: Modify `hierarchical_triangular_wl` signature

- Add `elabel_dict1=None, elabel_dict2=None`
- Add validation call
- Add `has_edge_labels` flag
- Initialize `ewl_np1/2` arrays when edge labels present

### Step 3: Modify K=0 path

- When `has_edge_labels`: use `edge_wl_test_triangulated_neighbors` (existing module)
- Convert `elabel_dict` to `edges_list` + `elabel_np`
- Return both node and edge label arrays

### Step 4: Modify K≥1 backward pass

- When `step == 1 and has_edge_labels`: use `_compute_lower_label_tuples_with_edges` instead of `_compute_lower_label_tuples`
- Precompute initial `ec(v)` before iteration loop (pre-processing)

### Step 5: Add edge label refresh after each iteration

- Call `_update_elabel_from_dict` after backward pass
- Recompute `ec(v)` from refreshed edge labels
- Store edge labels into `ewl_np1/2`

### Step 6: Add return logic

- If `has_edge_labels`: return `(vwl_np1, vwl_np2, ewl_np1, ewl_np2)`
- Else: return `(wl_np1, wl_np2)` (backward compatible)

### Step 7: Add tests

- Identical graphs with identical edge labels → full match
- Graphs differing only in edge labels → WL distinguishes
- K=0 with edge labels matches flat `edge_wl_test_triangulated_neighbors`
- K≥1 with uniform edge labels matches node-only result
- Edge label monotonicity (non-decreasing)
- Non-consecutive node IDs with edge labels
- Both node-labeled + edge-labeled vs edge-only vs node-only

---

## 7. Code Change Summary

| File | Change | Lines |
|------|--------|-------|
| `cyclic_schema/hierarchical_triangular_wl.py` | Add 4 helper functions | ~80 new |
| `cyclic_schema/hierarchical_triangular_wl.py` | Modify signature (lines 463-470) | ~5 changed |
| `cyclic_schema/hierarchical_triangular_wl.py` | Add validation after line 527 | ~15 new |
| `cyclic_schema/hierarchical_triangular_wl.py` | Init edge arrays (~line 542) | ~15 new |
| `cyclic_schema/hierarchical_triangular_wl.py` | Modify K=0 path (~lines 557-610) | ~20 changed |
| `cyclic_schema/hierarchical_triangular_wl.py` | Add pre-processing before iteration loop (~line 661) | ~15 new |
| `cyclic_schema/hierarchical_triangular_wl.py` | Modify backward pass at G level (~line 714-734) | ~10 changed |
| `cyclic_schema/hierarchical_triangular_wl.py` | Add edge refresh after iteration (~line 744) | ~20 new |
| `cyclic_schema/hierarchical_triangular_wl.py` | Modify return (~line 751) | ~5 changed |
| `cyclic_schema/hierarchical_triangular_wl.py` | Add edge label tests (~line 1103) | ~50 new |

Total: ~240 new/changed lines.

---

## 8. Verification

| Check | Method |
|-------|--------|
| WL consistency | Identical inputs → identical outputs |
| Edge-aware discrimination | Graphs differing only in edge labels produce different WL histories |
| K=0 correctness | Output matches existing `edge_wl_test_triangulated_neighbors` |
| K≥1 edge uniformity | All edges same label → output matches node-only |
| Monotonicity | All labels non-decreasing across iterations |
| Arbitrary node IDs | Test with string labels, non-consecutive ints |
| Both labels | Test with node label + edge label combinations |
