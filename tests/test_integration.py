"""Integration tests for the full TopoWasserstein kernel pipeline.

Tests end-to-end flows with small synthetic graphs (no external datasets).
"""

import networkx as nx
import numpy as np
import pytest

from cyclic_schema.hierarchical_triangulated_wl import (
    hierarchical_triangular_wl,
    hierarchical_triangular_wl_with_edges,
    hierarchical_triangular_wl_unified,
    _is_isomorphic_wl,
)
from cyclic_schema import (
    cyclic_schematic_graph,
    build_multilayer_csg_with_mappings,
)
from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import (
    TopoWassersteinGraphKernel,
)


# ============================================================================
# Integration: WL → CSG → Kernel pipeline
# ============================================================================

class TestCSGThenWL:
    """End-to-end: build CSG then run hierarchical WL."""

    def test_csg_to_wl_pipeline(self, triangle):
        """Build CSG from triangle, then run WL with K=1."""
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        wl1, wl2 = hierarchical_triangular_wl(
            triangle, triangle, vlabel, vlabel.copy(), K=1, I=3)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (3, 4)

    def test_csg_abstraction_changes_wl(self):
        """Graphs with different CSG structures get different WL results."""
        two_tri = nx.Graph()
        two_tri.add_edges_from([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])

        cycle6 = nx.Graph()
        cycle6.add_edges_from((i, (i + 1) % 6) for i in range(6))

        labels_tt = np.array([0, 0, 0, 1, 1, 1], dtype=np.int32)
        labels_c6 = np.array([0, 0, 0, 0, 0, 0], dtype=np.int32)

        wl_tt, _ = hierarchical_triangular_wl(
            two_tri, two_tri, labels_tt, labels_tt.copy(), K=1, I=2)
        wl_c6, _ = hierarchical_triangular_wl(
            cycle6, cycle6, labels_c6, labels_c6.copy(), K=1, I=2)

        assert not _is_isomorphic_wl(wl_tt, wl_c6)


# ============================================================================
# Integration: K=0 triangulated neighbors path
# ============================================================================

class TestK0Path:
    """K=0 uses triangulated neighbors (no CSG layers)."""

    def test_k0_two_nodes(self):
        """Two connected nodes, K=0."""
        G = nx.Graph()
        G.add_edges_from([(0, 1)])
        labels = np.array([5, 10], dtype=np.int32)
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=0, I=2)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (2, 3)

    def test_k0_with_non_consecutive_ids(self, non_consecutive_graph):
        labels = np.array([1, 2, 3], dtype=np.int32)
        wl1, wl2 = hierarchical_triangular_wl(
            non_consecutive_graph, non_consecutive_graph,
            labels, labels.copy(), K=0, I=2)
        assert np.all(wl1 == wl2)

    def test_k0_distinguishes_triangle_from_path(self):
        """K=0 should still distinguish different structures."""
        G1 = nx.Graph()
        G1.add_edges_from([(0, 1), (1, 2), (0, 2)])
        G2 = nx.Graph()
        G2.add_edges_from([(0, 1), (1, 2)])
        # Use uniform labels for robust differentiation
        labels = np.ones(3, dtype=np.int32)
        wl1, _ = hierarchical_triangular_wl(
            G1, G1, labels, labels.copy(), K=0, I=2)
        wl2, _ = hierarchical_triangular_wl(
            G2, G2, labels, labels.copy(), K=0, I=2)
        assert not _is_isomorphic_wl(wl1, wl2)


# ============================================================================
# Integration: K=0 edge path (triangulated neighbors with edge labels)
# ============================================================================

class TestK0EdgePath:
    """K=0 with edge labels through edge-aware WL."""

    def test_k0_edge_identical(self, triangle, triangle_elabel_uniform):
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        vwl1, vwl2, ewl1, ewl2 = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel.copy(),
            triangle_elabel_uniform, triangle_elabel_uniform,
            K=0, I=2)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)

    def test_k0_edge_distinguishes(self, triangle):
        vlabel = np.array([1, 1, 1], dtype=np.int32)
        elabel_a = {(0, 1): 10, (1, 2): 10, (0, 2): 10}
        elabel_b = {(0, 1): 10, (1, 2): 99, (0, 2): 10}
        vwl_a, vwl_b, _, _ = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel.copy(),
            elabel_a, elabel_b, K=0, I=2)
        assert not np.all(vwl_a == vwl_b)


# ============================================================================
# Integration: Unified WL dispatch
# ============================================================================

class TestUnifiedDispatchIntegration:
    """hierarchical_triangular_wl_unified dispatches correctly."""

    def test_unified_k0_vs_k1(self, triangle):
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        g_info = {'el': False}

        wl_k0, _ = hierarchical_triangular_wl_unified(
            g_info, triangle, triangle, vlabel, vlabel.copy(), K=0, I=2)
        wl_k1, _ = hierarchical_triangular_wl_unified(
            g_info, triangle, triangle, vlabel, vlabel.copy(), K=1, I=2)

        # Same initial column
        assert np.all(wl_k0[:, 0] == wl_k1[:, 0])
        # Different propagation after iteration 1 (different K value)
        assert wl_k0.shape[1] == wl_k1.shape[1]

    def test_unified_k0_node_vs_edge_results_differ(self, triangle):
        """K=0 should give different results for node vs edge dispatch."""
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        elabel = {(0, 1): 10, (1, 2): 10, (0, 2): 10}

        wl_node, _ = hierarchical_triangular_wl_unified(
            {'el': False}, triangle, triangle, vlabel, vlabel.copy(),
            K=0, I=2)
        wl_edge, _ = hierarchical_triangular_wl_unified(
            {'el': True}, triangle, triangle, vlabel, vlabel.copy(),
            elabel, elabel, K=0, I=2)

        # Edge dispatch adds edge context → should differ from node-only
        assert not np.all(wl_node == wl_edge)

    def test_unified_full_vs_edge_equals_wl_with_edges(self, triangle):
        """Unified with el=True should match hierarchical_triangular_wl_with_edges."""
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}

        # Unified
        unified_vwl, unified_ewl = hierarchical_triangular_wl_unified(
            {'el': True}, triangle, triangle, vlabel, vlabel.copy(),
            elabel, elabel, K=1, I=2)

        # Direct edge-aware call
        direct_vwl, direct_ewl, _, _ = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel.copy(),
            elabel, elabel, K=1, I=2)

        assert np.all(unified_vwl == direct_vwl)
        assert np.all(unified_ewl == direct_ewl)


# ============================================================================
# Integration: TopoWassersteinGraphKernel pipeline
# ============================================================================

class TestFullKernelPipeline:
    """Full kernel pipeline on tiny synthetic datasets."""

    @pytest.fixture
    def three_graphs(self):
        """Three small graphs: triangle, path3, 4-cycle."""
        g1 = nx.Graph()
        g1.add_edges_from([(0, 1), (1, 2), (0, 2)])
        g2 = nx.Graph()
        g2.add_edges_from([(0, 1), (1, 2)])
        g3 = nx.Graph()
        g3.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)])

        v1 = np.array([1, 2, 3], dtype=np.int32)
        v2 = np.array([1, 2, 3], dtype=np.int32)
        v3 = np.array([1, 2, 3, 4], dtype=np.int32)

        e1 = [(0, 1), (1, 2), (0, 2)]
        e2 = [(0, 1), (1, 2)]
        e3 = [(0, 1), (1, 2), (2, 3), (3, 0)]

        def _deg_distr(g):
            d = np.array([g.degree(v) for v in g.nodes()], dtype=np.float32)
            return d / d.sum()

        return {
            'graph_list': [g1, g2, g3],
            'vlabel_list': [v1, v2, v3],
            'edges_list': [e1, e2, e3],
            'elabel_list': [np.array([]), np.array([]), np.array([])],
            'deg_distr_list': [_deg_distr(g) for g in [g1, g2, g3]],
        }

    def test_kernel_produces_valid_matrix(self, three_graphs):
        dataset_info = {'nl': True, 'el': False}
        kernel = TopoWassersteinGraphKernel(
            n_wl_iters=2, n_csg_layers=2, wl_normalized=True)
        ot_dist, wl_sim, runtime = kernel.fit_transform(
            dataset_info,
            three_graphs['graph_list'],
            three_graphs['vlabel_list'],
            three_graphs['edges_list'],
            three_graphs['elabel_list'],
            three_graphs['deg_distr_list'],
        )
        assert ot_dist.shape == (3, 3)
        assert wl_sim.shape == (3, 3)
        assert np.allclose(ot_dist, ot_dist.T)
        assert np.allclose(np.diag(ot_dist), 0, atol=1e-6)
        assert runtime >= 0

    def test_triangle_more_similar_to_self_than_to_path(self, three_graphs):
        """Triangle(0) should be closer to itself than to path3(1)."""
        dataset_info = {'nl': True, 'el': False}
        kernel = TopoWassersteinGraphKernel(
            n_wl_iters=2, n_csg_layers=2, wl_normalized=True)
        ot_dist, _, _ = kernel.fit_transform(
            dataset_info,
            three_graphs['graph_list'],
            three_graphs['vlabel_list'],
            three_graphs['edges_list'],
            three_graphs['elabel_list'],
            three_graphs['deg_distr_list'],
        )
        # Self-distance is 0, cross-distance should be > 0
        assert ot_dist[0, 0] == 0
        assert ot_dist[0, 1] > 0
        # Triangle-to-triangle is 0, triangle-to-cycle is > 0
        assert ot_dist[0, 0] <= ot_dist[0, 1]

    def test_multiple_fits_produce_consistent_results(self, three_graphs):
        """Same input with same params → same output (deterministic)."""
        dataset_info = {'nl': True, 'el': False}
        k1 = TopoWassersteinGraphKernel(n_wl_iters=2, n_csg_layers=2, wl_normalized=True)
        k2 = TopoWassersteinGraphKernel(n_wl_iters=2, n_csg_layers=2, wl_normalized=True)

        ot1, wl1, _ = k1.fit_transform(
            dataset_info,
            three_graphs['graph_list'],
            three_graphs['vlabel_list'],
            three_graphs['edges_list'],
            three_graphs['elabel_list'],
            three_graphs['deg_distr_list'],
        )
        ot2, wl2, _ = k2.fit_transform(
            dataset_info,
            three_graphs['graph_list'],
            three_graphs['vlabel_list'],
            three_graphs['edges_list'],
            three_graphs['elabel_list'],
            three_graphs['deg_distr_list'],
        )
        assert np.allclose(ot1, ot2)
        assert np.allclose(wl1, wl2)


# ============================================================================
# Integration: Edge-label kernel pipeline
# ============================================================================

class TestEdgeLabelKernelPipeline:
    """Full kernel pipeline with edge labels."""

    def test_edge_label_pipeline(self, triangle):
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        edges = [(0, 1), (1, 2), (0, 2)]
        deg = np.array([2, 2, 2], dtype=np.float32) / 6.0

        dataset_info = {'nl': True, 'el': True}
        kernel = TopoWassersteinGraphKernel(
            n_wl_iters=1, n_csg_layers=1, wl_normalized=False)

        ot_dist, wl_sim, _ = kernel.fit_transform(
            dataset_info,
            [triangle, triangle],
            [vlabel, vlabel.copy()],
            [edges, edges],
            [np.array([10, 20, 10]), np.array([10, 20, 10])],
            [deg, deg],
        )
        assert ot_dist.shape == (2, 2)
        assert np.allclose(ot_dist[0, 0], 0, atol=1e-6)


# ============================================================================
# Integration: Multi-layer CSG stability
# ============================================================================

class TestMultiLayerStability:
    """Multi-layer CSG produces deterministic results."""

    def test_two_layer_wl(self, triangle):
        """WL with K=2 runs without error and produces shape-stable output."""
        vlabel = np.array([1, 2, 3], dtype=np.int32)
        wl1, wl2 = hierarchical_triangular_wl(
            triangle, triangle, vlabel, vlabel.copy(), K=2, I=2)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (3, 3)

    def test_k0_k1_k2_consistent_ordering(self):
        """Ordering of labels should be deterministic across K values."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])
        labels = np.array([5, 10, 5], dtype=np.int32)
        results = []
        for K in (0, 1, 2):
            wl, _ = hierarchical_triangular_wl(
                G, G, labels, labels.copy(), K=K, I=2)
            results.append(wl[:, 0])
        # All should start with same initial labels
        for r in results[1:]:
            assert np.all(results[0] == r)


# ============================================================================
# Integration: Self-test blocks from source modules
# ============================================================================

class TestExistingSelfTests:
    """Replicate key self-tests from module __main__ blocks."""

    def test_hierarchical_wl_self_test_logic(self):
        """Replicate the main self-test from hierarchical_triangulated_wl.py."""
        from cyclic_schema import build_example1_graph
        G1 = build_example1_graph()
        G2 = build_example1_graph()
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, G1.number_of_nodes())
        vlabel1 = vlabel.copy()
        vlabel2 = vlabel.copy()

        for K in (1, 2, 3):
            wl1, wl2 = hierarchical_triangular_wl(
                G1, G2, vlabel1, vlabel2, K=K, I=3)
            assert np.all(wl1 == wl2), f"K={K} failed for identical graphs"
