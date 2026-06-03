"""Unit tests for TopoWassersteinGraphKernel.

Tests static methods, initialization, and basic fit/transform
with small synthetic graphs (no external datasets).
"""

import networkx as nx
import numpy as np
import pytest

from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import (
    TopoWassersteinGraphKernel,
    _build_elabel_dict,
)


class TestBuildElabelDict:
    """Edge label dictionary construction."""

    def test_basic(self):
        edges = [(0, 1), (1, 2)]
        elabels = [10, 20]
        result = _build_elabel_dict(edges, elabels)
        assert result == {(0, 1): 10, (1, 2): 20}

    def test_empty(self):
        assert _build_elabel_dict([], []) == {}


class TestWLInnerProduct:
    """Static method for WL similarity computation."""

    def test_identical(self):
        a = np.array([[1, 2], [1, 3]], dtype=np.int32)
        b = a.copy()
        sim = TopoWassersteinGraphKernel.wl_inner_product(a, b)
        assert sim > 0

    def test_zero_for_disjoint_labels(self):
        a = np.array([[0, 0], [0, 0]], dtype=np.int32)
        b = np.array([[1, 1], [1, 1]], dtype=np.int32)
        sim = TopoWassersteinGraphKernel.wl_inner_product(a, b)
        assert sim == 0

    def test_symmetric(self):
        a = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int32)
        b = np.array([[3, 2, 1], [6, 5, 4]], dtype=np.int32)
        sim_ab = TopoWassersteinGraphKernel.wl_inner_product(a, b)
        sim_ba = TopoWassersteinGraphKernel.wl_inner_product(b, a)
        assert sim_ab == sim_ba

    def test_single_node(self):
        a = np.array([[5]], dtype=np.int32)
        b = np.array([[5]], dtype=np.int32)
        sim = TopoWassersteinGraphKernel.wl_inner_product(a, b)
        assert sim > 0


class TestKernelInitialization:
    """Kernel construction and parameter handling."""

    def test_default_params(self):
        kernel = TopoWassersteinGraphKernel(niter_tn=3, niter_hcc=3)
        assert kernel._niter_tn == 3
        assert kernel._niter_hcc == 3
        assert kernel._wl_normalized is True

    def test_custom_params(self):
        kernel = TopoWassersteinGraphKernel(
            niter_tn=5, niter_hcc=7, wl_normalized=False)
        assert kernel._niter_tn == 5
        assert kernel._niter_hcc == 7
        assert kernel._wl_normalized is False


class TestKernelFitTransform:
    """Minimal fit/transform with tiny synthetic graphs."""

    def _make_two_graphs(self):
        g1 = nx.Graph()
        g1.add_edges_from([(0, 1), (1, 2), (0, 2)])
        g2 = nx.Graph()
        g2.add_edges_from([(0, 1), (1, 2)])
        v1 = np.array([1, 2, 3], dtype=np.int32)
        v2 = np.array([1, 2, 3], dtype=np.int32)
        edges1 = [(0, 1), (1, 2), (0, 2)]
        edges2 = [(0, 1), (1, 2)]
        deg1 = np.array([2, 2, 2], dtype=np.float32)
        deg1 = deg1 / deg1.sum()
        deg2 = np.array([1, 2, 1], dtype=np.float32)
        deg2 = deg2 / deg2.sum()
        dataset_info = {'nl': True, 'el': False}
        return {
            'dataset_info': dataset_info,
            'graph_list': [g1, g2],
            'vlabel_list': [v1, v2],
            'edges_list': [edges1, edges2],
            'elabel_list': [np.array([]), np.array([])],
            'deg_distr_list': [deg1, deg2],
        }

    def test_fit_transform_returns_correct_shape(self):
        data = self._make_two_graphs()
        kernel = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=True)
        ot_dist, wl_sim, runtime = kernel.fit_transform(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        assert ot_dist.shape == (2, 2)
        assert wl_sim.shape == (2, 2)
        assert runtime >= 0

    def test_symmetric_ot_distance(self):
        data = self._make_two_graphs()
        kernel = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=True)
        ot_dist, _, _ = kernel.fit_transform(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        assert np.allclose(ot_dist, ot_dist.T)

    def test_zero_diagonal_ot(self):
        """OT distance on diagonal (graph vs itself) should be ~0."""
        data = self._make_two_graphs()
        kernel = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=True)
        ot_dist, _, _ = kernel.fit_transform(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        assert np.allclose(np.diag(ot_dist), 0, atol=1e-6)

    def test_wl_similarity_diagonal_positive(self):
        """WL similarity on diagonal should be > 0."""
        data = self._make_two_graphs()
        kernel = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=True)
        _, wl_sim, _ = kernel.fit_transform(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        assert np.all(np.diag(wl_sim) > 0)

    def test_fit_then_transform(self):
        """Explicit fit + transform should match fit_transform."""
        data = self._make_two_graphs()
        kernel = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=True)
        ot_ft, wl_ft, _ = kernel.fit_transform(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        kernel2 = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=True)
        kernel2.fit(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        ot_t, wl_t, _ = kernel2.transform()
        assert np.allclose(ot_ft, ot_t)
        assert np.allclose(wl_ft, wl_t)

    def test_node_only_vs_edge_label_dispatch(self):
        """Kernel should handle has_el=False and has_el=True gracefully."""
        data = self._make_two_graphs()

        # Node-only mode
        kernel_node = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=False)
        ot_node, wl_node, _ = kernel_node.fit_transform(
            data['dataset_info'], data['graph_list'],
            data['vlabel_list'], data['edges_list'],
            data['elabel_list'], data['deg_distr_list'],
        )
        assert ot_node.shape == (2, 2)

        # Edge-label mode (same data, but dispatch is different)
        data_el = data.copy()
        data_el['dataset_info'] = {'nl': True, 'el': True}
        data_el['elabel_list'] = [
            np.array([10, 10, 10]),
            np.array([10, 10]),
        ]
        kernel_el = TopoWassersteinGraphKernel(
            niter_tn=1, niter_hcc=1, wl_normalized=False)
        ot_el, wl_el, _ = kernel_el.fit_transform(
            data_el['dataset_info'], data_el['graph_list'],
            data_el['vlabel_list'], data_el['edges_list'],
            data_el['elabel_list'], data_el['deg_distr_list'],
        )
        assert ot_el.shape == (2, 2)
