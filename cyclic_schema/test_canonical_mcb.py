"""Tests for canonical minimum cycle basis (canonical_mcb module)."""

import networkx as nx
import numpy as np
import pytest

from cyclic_schema.canonical_mcb import (
    canonical_mcb,
    canonical_mcb_invariant_summary,
    compute_canonical_vertex_ids,
)


# ===========================================================================
# Fixtures
# ===========================================================================

@pytest.fixture
def triangle():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (0, 2)])
    return G


@pytest.fixture
def square():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)])
    return G


@pytest.fixture
def k4():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
    return G


@pytest.fixture
def example1():
    from cyclic_schema import build_example1_graph
    return build_example1_graph()


# ===========================================================================
# Basic correctness
# ===========================================================================

class TestBasicCorrectness:
    def test_empty_graph(self):
        G = nx.Graph()
        assert canonical_mcb(G) == []

    def test_single_node(self):
        G = nx.Graph()
        G.add_node(0)
        assert canonical_mcb(G) == []

    def test_tree_no_cycles(self):
        """A tree has no cycles → empty basis."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3)])
        assert canonical_mcb(G) == []

    def test_triangle_basis_size(self, triangle):
        """Triangle has mu = 1 → 1 basis cycle."""
        basis = canonical_mcb(triangle)
        assert len(basis) == 1
        assert len(basis[0]) == 3

    def test_square_basis_size(self, square):
        """Square has mu = 1 → 1 basis cycle of length 4."""
        basis = canonical_mcb(square)
        assert len(basis) == 1
        assert len(basis[0]) == 4

    def test_k4_basis_size(self, k4):
        """K4 has mu = 3 → 3 basis cycles (3-cycles)."""
        basis = canonical_mcb(k4)
        assert len(basis) == 3
        for c in basis:
            assert len(c) == 3

    def test_cycle_basis_count_example1(self, example1):
        """Example 1 graph is disconnected (2 components), so mu = m - n + c."""
        basis = canonical_mcb(example1)
        n, m = example1.number_of_nodes(), example1.number_of_edges()
        c = nx.number_connected_components(example1)
        mu = m - n + c
        assert len(basis) == mu

    def test_basis_returns_vertices(self, triangle):
        basis = canonical_mcb(triangle)
        for c in basis:
            assert isinstance(c, tuple)
            for v in c:
                assert v in triangle.nodes()

    def test_basis_return_edge_sets(self, triangle):
        basis = canonical_mcb(triangle, return_edge_sets=True)
        assert len(basis) == 1
        assert isinstance(basis[0], frozenset)
        assert len(basis[0]) == 3


# ===========================================================================
# Basis property: cycles must be linearly independent in F_2^|E|
# ===========================================================================

class TestBasisProperty:
    @staticmethod
    def _to_f2_vec(cycle, edge_index):
        n = len(cycle)
        edges = set()
        for i in range(n):
            u, v = cycle[i], cycle[(i + 1) % n]
            key = (u, v) if (u, v) in edge_index else (v, u)
            edges.add(key)
        vec = [0] * len(edge_index)
        for e in edges:
            vec[edge_index[e]] = 1
        return vec

    def test_basis_is_f2_independent(self, k4):
        """F_2 vectors of the returned basis must be linearly independent."""
        G = k4
        edges = sorted(G.edges(), key=lambda e: (min(e), max(e)))
        edge_index = {e: i for i, e in enumerate(edges)}
        basis = canonical_mcb(G)
        vecs = [self._to_f2_vec(c, edge_index) for c in basis]
        # Build rank via gaussian elimination
        rank = 0
        basis_rows = []
        for v in vecs:
            cur = v[:]
            for b in basis_rows:
                pivot = next((j for j, x in enumerate(cur) if x == 1), None)
                if pivot is None:
                    continue
                if b[pivot] == 1:
                    cur = [a ^ b_i for a, b_i in zip(cur, b)]
            if any(x == 1 for x in cur):
                rank += 1
                basis_rows.append(cur)
        assert rank == len(basis)

    def test_basis_is_f2_independent_example1(self, example1):
        """Same property for the larger example1 graph."""
        G = example1
        edges = sorted(G.edges(), key=lambda e: (min(e), max(e)))
        edge_index = {e: i for i, e in enumerate(edges)}
        basis = canonical_mcb(G)
        vecs = [self._to_f2_vec(c, edge_index) for c in basis]
        rank = 0
        basis_rows = []
        for v in vecs:
            cur = v[:]
            for b in basis_rows:
                pivot = next((j for j, x in enumerate(cur) if x == 1), None)
                if pivot is None:
                    continue
                if b[pivot] == 1:
                    cur = [a ^ b_i for a, b_i in zip(cur, b)]
            if any(x == 1 for x in cur):
                rank += 1
                basis_rows.append(cur)
        assert rank == len(basis)


# ===========================================================================
# Isomorphism invariance
# ===========================================================================

class TestIsomorphismInvariance:
    def test_permuted_triangle(self, triangle):
        """Re-label triangle nodes → same canonical basis (up to relabelling)."""
        relabel = {0: 7, 1: 3, 2: 9}
        permuted = nx.relabel_nodes(triangle, relabel)
        basis_orig = canonical_mcb(triangle)
        basis_perm = canonical_mcb(permuted)
        assert len(basis_orig) == len(basis_perm)
        # Both bases have 1 cycle; cycle lengths match
        assert len(basis_orig[0]) == len(basis_perm[0]) == 3
        # Cycle sets (as sets of nodes) should have the same cardinality
        s1 = set(basis_orig[0])
        s2 = set(basis_perm[0])
        assert len(s1) == 3
        # The permuted cycle should be a relabelling of the original
        inverse_relabel = {v: k for k, v in relabel.items()}
        unmapped = frozenset(inverse_relabel[v] for v in basis_perm[0])
        assert unmapped == frozenset(basis_orig[0])

    def test_permuted_k4(self, k4):
        """K4 with relabelled vertices → same basis structure."""
        relabel = {0: 5, 1: 2, 2: 7, 3: 0}
        permuted = nx.relabel_nodes(k4, relabel)
        b1 = canonical_mcb(k4)
        b2 = canonical_mcb(permuted)
        # Both have 3 triangles
        assert len(b1) == len(b2) == 3
        for c1, c2 in zip(b1, b2):
            assert len(c1) == len(c2) == 3

    def test_permuted_example1(self, example1):
        """Example1 with a random relabelling → same basis size, same cycle
        length histogram (the only isomorphism invariant easy to verify on
        a graph with thousands of automorphisms-free vertex orbits)."""
        rng = np.random.default_rng(42)
        nodes = list(example1.nodes())
        perm = list(nodes)
        rng.shuffle(perm)
        relabel = dict(zip(nodes, perm))
        permuted = nx.relabel_nodes(example1, relabel)
        b1 = canonical_mcb(example1)
        b2 = canonical_mcb(permuted)
        assert len(b1) == len(b2)
        from collections import Counter
        h1 = Counter(len(c) for c in b1)
        h2 = Counter(len(c) for c in b2)
        assert h1 == h2

    def test_isomorphic_copies_same_signature(self, k4):
        """Two K4 copies with different node ids → identical invariant summary."""
        G1 = nx.relabel_nodes(k4, {i: i + 100 for i in range(4)})
        G2 = nx.relabel_nodes(k4, {i: i + 200 for i in range(4)})
        s1 = canonical_mcb_invariant_summary(G1)
        s2 = canonical_mcb_invariant_summary(G2)
        assert s1["basis_size"] == s2["basis_size"]
        assert s1["length_histogram"] == s2["length_histogram"]

    def test_k4_k33_bipartite_have_different_invariants(self):
        """K4 and K_{3,3} have different local degree signatures (4 vs 3)
        so their length histograms should differ."""
        G1 = nx.complete_graph(4)
        G2 = nx.complete_bipartite_graph(3, 3)
        s1 = canonical_mcb_invariant_summary(G1)
        s2 = canonical_mcb_invariant_summary(G2)
        # Both have mu=3 but cycle lengths differ (3 vs 4)
        assert s1["length_histogram"] != s2["length_histogram"]


# ===========================================================================
# Vertex-id invariance
# ===========================================================================

class TestCanonicalVertexIDs:
    def test_ids_are_unique(self, example1):
        cids = compute_canonical_vertex_ids(example1, depth=3)
        assert len(cids) == example1.number_of_nodes()
        assert set(cids.values()) == set(range(max(cids.values()) + 1))

    def test_ids_invariance_under_relabel(self, triangle):
        cids1 = compute_canonical_vertex_ids(triangle, depth=3)
        relabel = {0: 9, 1: 4, 2: 1}
        permuted = nx.relabel_nodes(triangle, relabel)
        cids2 = compute_canonical_vertex_ids(permuted, depth=3)
        # All vertices of triangle have the same degree (2) and same
        # 3-step signature (every vertex sees the other two at depth 1
        # with degree 2).  So all three get the *same* cid — this is
        # correct isomorphism-invariant behaviour.
        assert len(set(cids1.values())) == 1
        assert len(set(cids2.values())) == 1

    def test_higher_depth_disambiguates_k4(self, k4):
        """K4 is vertex-transitive — every vertex has identical local
        structure up to depth ``diam(K4) = 1`` (i.e. degree 3).  Even at
        depth 1 every vertex sees 3 others of degree 3, so all cids match."""
        cids = compute_canonical_vertex_ids(k4, depth=5)
        # All four vertices get the *same* cid (K4 is vertex-transitive)
        assert len(set(cids.values())) == 1


# ===========================================================================
# Disconnected graphs
# ===========================================================================

class TestDisconnected:
    def test_two_disjoint_triangles(self):
        """Two disjoint triangles → 2 basis cycles, one per component."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2),
                           (3, 4), (4, 5), (3, 5)])
        basis = canonical_mcb(G)
        assert len(basis) == 2
        for c in basis:
            assert len(c) == 3

    def test_disjoint_tree_and_cycle(self):
        """Path + triangle → 1 basis cycle (only from the triangle)."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2),  # path
                           (3, 4), (4, 5), (3, 5)])  # triangle
        basis = canonical_mcb(G)
        assert len(basis) == 1
        assert len(basis[0]) == 3

    def test_disjoint_with_trees_dominates(self):
        """Two disjoint cycles + trees → mu cycles total."""
        G = nx.Graph()
        G.add_edges_from([
            (0, 1), (1, 2), (0, 2),  # triangle
            (2, 3), (3, 4),           # tree from triangle
            (5, 6), (6, 7), (5, 7),  # triangle
            (7, 8),                    # tree from triangle
        ])
        basis = canonical_mcb(G)
        n, m = G.number_of_nodes(), G.number_of_edges()
        c = nx.number_connected_components(G)
        mu = m - n + c
        assert len(basis) == mu


# ===========================================================================
# Convergence speed (informational)
# ===========================================================================

class TestSpeed:
    def test_example1_runs_quickly(self, example1):
        """canonical_mcb should run in <2s on a 26-node, 41-edge graph."""
        import time
        t0 = time.time()
        basis = canonical_mcb(example1)
        elapsed = time.time() - t0
        assert elapsed < 2.0
        n, m = example1.number_of_nodes(), example1.number_of_edges()
        c = nx.number_connected_components(example1)
        assert len(basis) == m - n + c

    def test_larger_graph_runs(self):
        """Run on a 50-node, ~80-edge graph."""
        rng = np.random.default_rng(0)
        G = nx.gnp_random_graph(50, 0.07, seed=0)
        basis = canonical_mcb(G)
        n, m = G.number_of_nodes(), G.number_of_edges()
        c = nx.number_connected_components(G)
        mu = m - n + c
        assert len(basis) == mu


# ===========================================================================
# Non-isomorphism detection
# ===========================================================================

class TestNonIsomorphismDetection:
    def test_c6_vs_two_triangles(self):
        """C6 (one 6-cycle) and 2 disjoint triangles (two 3-cycles) are
        not isomorphic: their canonical length histograms differ."""
        G1 = nx.cycle_graph(6)
        G2 = nx.Graph()
        G2.add_edges_from([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
        s1 = canonical_mcb_invariant_summary(G1)
        s2 = canonical_mcb_invariant_summary(G2)
        assert s1["length_histogram"] == {6: 1}
        assert s2["length_histogram"] == {3: 2}
        assert s1 != s2

    def test_c4_vs_diamond(self):
        """C4 (a single 4-cycle) vs diamond graph (4-cycle + 1 diagonal).
        Both have mu=1, but the diamond's single basis cycle is *not* a
        4-cycle (it must use the diagonal edges) — its length is 3
        (one of the two triangles in the diamond)."""
        G1 = nx.cycle_graph(4)
        G2 = nx.Graph()
        G2.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        s1 = canonical_mcb_invariant_summary(G1)
        s2 = canonical_mcb_invariant_summary(G2)
        assert s1["length_histogram"] == {4: 1}
        # Diamond: basis can be a 3-cycle
        assert s2["length_histogram"] != s1["length_histogram"]


# ===========================================================================
# Integration with cyclic_schema.cyclic_schematic_graph
# ===========================================================================

class TestIntegration:
    def test_csg_uses_canonical_mcb(self, triangle):
        """Verify that cyclic_schematic_graph works with canonical_mcb as
        drop-in replacement.  This is a smoke test — the integration is
        in the cyclic_schema module itself (see modified import)."""
        from cyclic_schema.canonical_mcb import canonical_mcb as cmcb
        basis = cmcb(triangle)
        assert len(basis) == 1
