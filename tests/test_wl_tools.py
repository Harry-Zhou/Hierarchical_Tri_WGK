"""Unit tests for hierarchical_tri_wl_tools sub-modules.

Tests classic_wl_tools and classic_edge_wl_tools.
"""

import networkx as nx
import numpy as np
import pytest

from hierarchical_tri_wl_tools.classic_wl_tools import (
    classic_wl_test,
    propagate_vlabel,
    compress_and_relabel_vlabel,
    gdv_WL_test,
)
from hierarchical_tri_wl_tools.classic_edge_wl_tools import (
    collect_vlabels,
    compress_and_relabel_elabel,
    classic_edge_WL_test,
)


# ============================================================================
# classic_wl_tools
# ============================================================================

class TestClassicWL:
    """Classic WL label propagation."""

    def test_identical_graphs(self, triangle):
        g1 = g2 = triangle
        v1 = np.array([1, 2, 3], dtype=np.int32)
        v2 = np.array([1, 2, 3], dtype=np.int32)
        wl1, wl2 = classic_wl_test(g1, g2, v1, v2, niter=3)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (3, 4)

    def test_non_isomorphic_different_sizes(self):
        g1 = nx.Graph()
        g1.add_edges_from([(0, 1), (1, 2)])
        g2 = nx.Graph()
        g2.add_edges_from([(0, 1)])
        v1 = np.array([1, 2, 3], dtype=np.int32)
        v2 = np.array([1, 2], dtype=np.int32)
        wl1, wl2 = classic_wl_test(g1, g2, v1, v2, niter=2)
        assert wl1.shape[0] != wl2.shape[0]

    def test_monotonic_labels(self, triangle):
        v = np.array([5, 10, 15], dtype=np.int32)
        wl1, wl2 = classic_wl_test(triangle, triangle, v, v.copy(), niter=4)
        for col in range(1, wl1.shape[1]):
            assert np.all(wl1[:, col] >= wl1[:, col - 1])

    def test_label_permutation_invariant(self, triangle):
        """Permuting labels should preserve WL histogram equivalence."""
        v1 = np.array([10, 20, 30], dtype=np.int32)
        v2 = np.array([30, 10, 20], dtype=np.int32)
        wl1, wl2 = classic_wl_test(triangle, triangle, v1, v2, niter=3)

        # Check histograms after propagation
        from collections import Counter
        for col in range(wl1.shape[1]):
            h1 = Counter(wl1[:, col].tolist())
            h2 = Counter(wl2[:, col].tolist())
            assert h1 == h2, f"Column {col} histograms differ"

    def test_zero_iterations(self, triangle):
        v = np.array([1, 2, 3], dtype=np.int32)
        wl1, wl2 = classic_wl_test(triangle, triangle, v, v.copy(), niter=0)
        assert np.all(wl1[:, 0] == v)
        assert np.all(wl2[:, 0] == v)
        assert wl1.shape == (3, 1)

    def test_single_node(self, single_node):
        v = np.array([42], dtype=np.int32)
        wl1, wl2 = classic_wl_test(single_node, single_node, v, v.copy(), niter=2)
        assert wl1.shape == (1, 3)
        assert np.all(wl1 == wl2)


class TestPropagateVLabel:
    """Node label propagation helper."""

    def test_triangle(self, triangle):
        node_neighs = {v: list(triangle.neighbors(v)) for v in triangle.nodes()}
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        result = propagate_vlabel(node_neighs, vlabel)
        assert len(result) == 3
        for v in triangle.nodes():
            assert isinstance(result[v], tuple)

    def test_isolated_node(self):
        node_neighs = {0: []}
        vlabel = np.array([5], dtype=np.int32)
        result = propagate_vlabel(node_neighs, vlabel)
        assert result[0] == (5,)


class TestCompressRelabel:
    """Label compression across two graphs."""

    def test_identical_labels(self, triangle):
        v1 = {0: (1, 2), 1: (2, 3), 2: (3, 1)}
        v2 = {0: (1, 2), 1: (2, 3), 2: (3, 1)}
        c1, c2 = compress_and_relabel_vlabel(v1, v2)
        assert len(c1) == 3
        assert len(c2) == 3

    def test_distinct_labels(self, triangle):
        v1 = {0: (1, 2), 1: (3, 4), 2: (5, 6)}
        v2 = {0: (7, 8), 1: (9, 10), 2: (11, 12)}
        c1, c2 = compress_and_relabel_vlabel(v1, v2)
        assert not np.any(c1 == c2)


# ============================================================================
# classic_edge_wl_tools
# ============================================================================

class TestEdgeWL:
    """Edge-aware WL propagation."""

    def test_identical_graphs(self, triangle):
        v = np.array([1, 2, 3], dtype=np.int32)
        edges = [(0, 1), (1, 2), (0, 2)]
        elabels = np.array([10, 20, 10], dtype=np.int32)
        vwl1, vwl2, ewl1, ewl2 = classic_edge_WL_test(
            triangle, triangle, v, v.copy(),
            edges, edges, elabels, elabels, niter=2)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)

    def test_different_edge_labels_distinguish(self, triangle):
        v = np.array([1, 1, 1], dtype=np.int32)
        edges = [(0, 1), (1, 2), (0, 2)]
        elabels_a = np.array([10, 10, 10], dtype=np.int32)
        elabels_b = np.array([10, 99, 10], dtype=np.int32)
        vwl_a, vwl_b, _, _ = classic_edge_WL_test(
            triangle, triangle, v, v.copy(),
            edges, edges, elabels_a, elabels_b, niter=2)
        assert not np.all(vwl_a == vwl_b)

    def test_uniform_edge_labels_match(self, triangle):
        v = np.array([1, 2, 3], dtype=np.int32)
        edges = [(0, 1), (1, 2), (0, 2)]
        elabels = np.array([0, 0, 0], dtype=np.int32)
        vwl1, vwl2, ewl1, ewl2 = classic_edge_WL_test(
            triangle, triangle, v, v.copy(),
            edges, edges, elabels, elabels, niter=2)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)

    def test_three_node_line(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2)])
        v = np.array([1, 2, 3], dtype=np.int32)
        edges = [(0, 1), (1, 2)]
        elabels = np.array([10, 20], dtype=np.int32)
        vwl1, vwl2, ewl1, ewl2 = classic_edge_WL_test(
            G, G,
            v, v.copy(), edges, edges, elabels, elabels, niter=2)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)


class TestCollectVLabels:
    """Edge label collection from node labels."""

    def test_triangle(self, triangle):
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        elabel_dict = {(0, 1): 10, (1, 2): 20, (0, 2): 30}
        result = collect_vlabels(vlabel, elabel_dict)
        assert len(result) == 3
        for e, t in result.items():
            assert len(t) == 3  # (edge_label, vlabel_lhs, vlabel_rhs)
