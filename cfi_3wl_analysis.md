# Why 3-WL Cannot Distinguish CFI(K4) and CFI'(K4)

## 1. How 3-WL Works

3-WL (3-dimensional Weisfeiler-Lehman) distinguishes nodes through iterative refinement:
- Iteration 1: Node's degree/neighbors
- Iteration 2: Neighbor's neighbors (2-hop)
- Iteration 3: Neighbor's neighbor's neighbors (3-hop)

If node colors no longer change after iterations, 3-WL considers them equivalent.

## 2. Local Structures Are Identical

| Feature | CFI(K4) | CFI'(K4) |
|---------|---------|----------|
| Nodes | 20 | 20 |
| Edges | 24 | 24 |
| Degree distribution | Identical | Identical |
| 1-hop neighbors | Identical | Identical |
| 2-hop neighbors | Identical | Identical |
| 3-hop neighbors | Identical | Identical |

**Reason**: CFI construction preserves symmetry—all vertex-clones (v*_0) have the same degree, and all edge-nodes (e*_a or e*_b) have the same degree.

## 3. What 3-WL Cannot Capture: Cycle Parity

The fundamental difference between the two graphs:

```
CFI(K4) triangle cycle: v1_0 → e12_a → v2_0 → e23_a → v3_0 → e13_a → v1_0
                         ↓        ↓        ↓
                        parallel parallel parallel → 0 twisted edges (even)

CFI'(K4) triangle cycle: v1_0 → e12_a → v2_0 → e23_a → v3_0 → e13_a → v1_0
                         ↓        ↓        ↓
                        twisted  twisted  twisted  → 3 twisted edges (odd)
```

**Problem**: Detecting this difference requires traversing the entire cycle (6-hop), but 3-WL only looks at 3-hop neighborhoods.

## 4. Root Cause: Treewidth

| Base Graph | Treewidth | k-WL Distinguishability Condition |
|------------|-----------|----------------------------------|
| K4 | 3 | k > 3 (requires 4-WL) |

- Treewidth of K4 = 3
- Distinguishability of CFI(K4) depends on detecting "circuits around K4"
- 3-WL's capability is equivalent to detecting structures with treewidth < 3
- K4 is exactly at the boundary, so 3-WL **cannot distinguish** it

## 5. Conclusion

**Root Cause**: 3-WL's expressive power is limited to 3-hop neighborhoods and cannot capture global topological information like "parity of twisted edges in a cycle" which requires a longer path to observe.

The difference in minimum cycle basis (6 cycles of length 6 vs 5 cycles of length 8) is exactly this global topological difference—it requires 4-WL to detect.

## 6. Supporting Evidence

### Minimum Cycle Basis Comparison

| Metric | CFI(K4) | CFI'(K4) |
|--------|---------|----------|
| Number of cycles | 6 | 5 |
| Cycle length | 6 | 8 |
| Vertex clones | Only v*_0 or v*_1 | Mix of v*_0 and v*_1 |

### CFI(K4) - 6 cycles (all length 6)

| # | Cycle Path |
|---|------------|
| 1 | e13_a → v2_0 → e23_a → e12_a → v1_0 → v3_0 |
| 2 | e14_a → v2_0 → v4_0 → e12_a → v1_0 → e24_a |
| 3 | e13_a → e14_a → v4_0 → v1_0 → v3_0 → e34_a |
| 4 | e23_b → e12_b → v3_1 → v2_1 → v1_1 → e13_b |
| 5 | e24_b → e12_b → e14_b → v2_1 → v1_1 → v4_1 |
| 6 | e14_b → v3_1 → v1_1 → e34_b → v4_1 → e13_b |

### CFI'(K4) - 5 cycles (all length 8)

| # | Cycle Path |
|---|------------|
| 1 | e24_b → v3_1 → e12_a → e13_a → v2_1 → v1_0 → e34_b → v4_1 |
| 2 | v3_1 → v4_0 → e13_a → e14_a → v2_0 → e23_a → v1_0 → e24_a |
| 3 | e23_b → v4_0 → e12_a → e34_a → e14_a → v2_1 → v1_0 → v3_0 |
| 4 | e12_b → v4_0 → e34_a → v2_0 → v1_1 → e24_a → v3_0 → e13_b |
| 5 | e24_b → e23_b → e14_b → v2_1 → v1_1 → v3_0 → v4_1 → e13_b |

The difference in cycle lengths (6 vs 8) and the mixing of vertex clone types in CFI'(K4) cycles directly reflects the twisted edge structure that 3-WL cannot detect.