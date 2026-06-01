"""Tests for hierarchical triangulated WL (Step-1/Step-2/Step-3 + edge labels)."""

import networkx as nx
import numpy as np
import pytest

from cyclic_schema.hierarchical_triangulated_wl import (
    canonicalize_cycle_label,
    compute_initial_label_tuples,
    compute_final_label_tuples,
    forward_aggregate,
    forward_message_passing_both,
    backward_message_passing_both,
    backward_message_passing_both_with_edges,
    hierarchical_triangular_wl,
    hierarchical_triangular_wl_with_edges,
    hierarchical_triangular_wl_unified,
    _compute_edge_context,
    _compute_lower_label_tuples,
    _compute_lower_label_tuples_with_edges,
    _update_elabel_from_dict,
    _labels_to_dict,
    _precompute_neighbor_components,
    _is_isomorphic_wl,
)
from cyclic_schema import (
    cyclic_schematic_graph,
    build_example1_graph,
    build_csg_to_input_mapping,
    build_input_to_csg_mapping,
    get_node_type,
)

# ===========================================================================
# Helper fixtures
# ===========================================================================

@pytest.fixture
def triangle():
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (0, 2)])
    return G

@pytest.fixture
def example1():
    return build_example1_graph()

@pytest.fixture
def path4():
    G = nx.Graph()
    G.add_edges_from((i, i + 1) for i in range(3))
    return G

# ===========================================================================
# Cycle canonicalization unit tests
# ===========================================================================

class TestCycleCanonicalization:
    def test_already_canonical(self):
        assert canonicalize_cycle_label((1, 2, 3, 4)) == (1, 2, 3, 4)

    def test_min_at_nonzero_position(self):
        assert canonicalize_cycle_label((3, 1, 2, 4)) == (1, 2, 4, 3)

    def test_reversed_cycle(self):
        assert canonicalize_cycle_label((1, 4, 3, 2)) == (1, 2, 3, 4)

    def test_min_middle(self):
        assert canonicalize_cycle_label((2, 1, 4, 3)) == (1, 2, 3, 4)

    def test_single_element(self):
        assert canonicalize_cycle_label((1,)) == (1,)

    def test_empty(self):
        assert canonicalize_cycle_label(()) == ()

    def test_all_equal(self):
        assert canonicalize_cycle_label((1, 1, 1, 1)) == (1, 1, 1, 1)

    def test_two_elements(self):
        result = canonicalize_cycle_label((2, 1))
        assert result == (1, 2) or result == (2, 1)

    def test_uniqueness_invariant(self):
        """Same multiset of labels on same cycle topology → same canonical form."""
        r1 = canonicalize_cycle_label((5, 3, 7, 1))
        r2 = canonicalize_cycle_label((5, 3, 7, 1))
        assert r1 == r2


# ===========================================================================
# Identical graph tests (K=0 through K=5)
# ===========================================================================

class TestIdenticalGraphs:
    def test_identical_k0(self, example1):
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        wl1, wl2 = hierarchical_triangular_wl(
            example1, example1, vlabel, vlabel.copy(), K=0, I=3)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (26, 4)

    @pytest.mark.parametrize("K", [1, 2, 3])
    def test_identical_k1_to_k3(self, example1, K):
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        wl1, wl2 = hierarchical_triangular_wl(
            example1, example1, vlabel, vlabel.copy(), K=K, I=3)
        assert np.all(wl1 == wl2)

    @pytest.mark.parametrize("K", [4, 5])
    def test_identical_k4_k5_stress(self, example1, K):
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        wl1, wl2 = hierarchical_triangular_wl(
            example1, example1, vlabel, vlabel.copy(), K=K, I=2)
        assert np.all(wl1 == wl2)

    def test_identical_triangle_label_permutation(self, triangle):
        """Triangle with different label orderings — same multiset per column."""
        labels_a = np.array([5, 10, 15])
        labels_b = np.array([15, 5, 10])
        wl_a, wl_b = hierarchical_triangular_wl(
            triangle, triangle, labels_a, labels_b, K=1, I=3)
        assert _is_isomorphic_wl(wl_a, wl_b)

    def test_dict_form_equals_array_form(self, example1):
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        labels_dict = {v: int(vlabel[i])
                       for i, v in enumerate(sorted(example1.nodes()))}
        wl_d, _ = hierarchical_triangular_wl(
            example1, example1, labels_dict, labels_dict, K=1, I=2)
        wl_a, _ = hierarchical_triangular_wl(
            example1, example1, vlabel, vlabel.copy(), K=1, I=2)
        assert np.all(wl_d == wl_a)

    def test_all_k_produce_identical_for_identical(self, example1):
        """All K=0..3 produce identical labels for identical inputs."""
        np.random.seed(1)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        for K in (0, 1, 2, 3):
            wl_a, wl_b = hierarchical_triangular_wl(
                example1, example1, vlabel, vlabel.copy(), K=K, I=2)
            assert np.all(wl_a == wl_b), f"K={K} failed"


# ===========================================================================
# Non-isomorphic graph discrimination
# ===========================================================================

class TestDiscrimination:
    def test_cycle_vs_path(self):
        """WL must distinguish 6-cycle from 5-path."""
        cycle6 = nx.Graph()
        cycle6.add_edges_from((i, (i + 1) % 6) for i in range(6))
        path5 = nx.Graph()
        path5.add_edges_from((i, i + 1) for i in range(4))
        labels6 = np.array([0, 1, 2, 0, 1, 2])
        labels5 = np.array([0, 1, 2, 0, 1])
        wl_c, wl_p = hierarchical_triangular_wl(
            cycle6, path5, labels6, labels5, K=1, I=3)
        assert not _is_isomorphic_wl(wl_c, wl_p)

    def test_triangle_vs_path3(self, triangle):
        """Triangle (3 nodes) vs path (3 nodes) — different cycle structure."""
        path3 = nx.Graph()
        path3.add_edges_from([(0, 1), (1, 2)])
        labels = np.array([1, 2, 3])
        wl_tri, wl_path = hierarchical_triangular_wl(
            triangle, path3, labels, labels, K=1, I=2)
        assert not _is_isomorphic_wl(wl_tri, wl_path)

    def test_two_disjoint_triangles_vs_6cycle(self):
        """Two disjoint triangles vs 6-cycle — both 6 nodes, different cycle structure."""
        two_tri = nx.Graph()
        two_tri.add_edges_from([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
        cycle6 = nx.Graph()
        cycle6.add_edges_from((i, (i + 1) % 6) for i in range(6))
        labels = np.array([0, 0, 0, 1, 1, 1])
        wl_tri, wl_cycle = hierarchical_triangular_wl(
            two_tri, cycle6, labels, labels, K=1, I=2)
        assert not _is_isomorphic_wl(wl_tri, wl_cycle)


# ===========================================================================
# Edge cases
# ===========================================================================

class TestEdgeCases:
    def test_acyclic_graph(self):
        """Path with no cycles should still work with K>=1."""
        path = nx.Graph()
        path.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4)])
        labels = np.array([0, 1, 0, 1, 0])
        wl1, wl2 = hierarchical_triangular_wl(
            path, path, labels, labels.copy(), K=1, I=2)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (5, 3)

    def test_single_node(self):
        G = nx.Graph()
        G.add_node(0)
        labels = np.array([42])
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=1, I=2)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (1, 3)

    def test_empty_graph(self):
        G = nx.Graph()
        labels = np.array([], dtype=np.int32)
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=1, I=2)
        assert wl1.shape == (0, 3)

    def test_non_consecutive_node_ids(self):
        G = nx.Graph()
        G.add_edges_from([(0, 5), (5, 10), (10, 0)])
        labels = np.array([1, 2, 3])
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=1, I=2)
        assert np.all(wl1 == wl2)

    def test_k0_on_acyclic(self):
        """K=0 with acyclic graph — triangulated neighbors WL."""
        path = nx.Graph()
        path.add_edges_from([(0, 1), (1, 2), (2, 3)])
        labels = np.array([3, 1, 4, 1])
        wl1, wl2 = hierarchical_triangular_wl(
            path, path, labels, labels.copy(), K=0, I=2)
        assert np.all(wl1 == wl2)

    def test_chord_cycle(self):
        """4-cycle with a chord — tests cycle bases that share edges."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        labels = np.array([1, 2, 3, 4])
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=1, I=2)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (4, 3)

    @pytest.mark.parametrize("K", [0, 1, 2])
    def test_k0_and_k1_consistent_ordering(self, K):
        """Identical graphs with same labels → identical output regardless of K."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])
        labels = np.array([5, 10, 5])
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=K, I=2)
        assert np.all(wl1 == wl2)


# ===========================================================================
# Label monotonicity (WL labels only increase)
# ===========================================================================

class TestMonotonicity:
    def test_labels_non_decreasing(self, example1):
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        wl1, wl2 = hierarchical_triangular_wl(
            example1, example1, vlabel, vlabel.copy(), K=2, I=4)
        for col in range(1, wl1.shape[1]):
            assert np.all(wl1[:, col] >= wl1[:, col - 1])
            assert np.all(wl2[:, col] >= wl2[:, col - 1])


# ===========================================================================
# Error handling
# ===========================================================================

class TestErrorHandling:
    def test_negative_k(self, triangle):
        vlabel = np.array([1, 2, 3])
        with pytest.raises(ValueError, match="K must be >= 0"):
            hierarchical_triangular_wl(triangle, triangle, vlabel, vlabel, K=-1, I=2)

    def test_zero_iterations(self, triangle):
        vlabel = np.array([1, 2, 3])
        with pytest.raises(ValueError, match="I must be >= 1"):
            hierarchical_triangular_wl(triangle, triangle, vlabel, vlabel, K=1, I=0)


# ===========================================================================
# Edge-label tests
# ===========================================================================

class TestEdgeLabels:
    def test_identical_edge_labels_match(self, triangle):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        vwl1, vwl2, ewl1, ewl2 = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel, elabel, elabel, K=1, I=3)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)
        assert vwl1.shape == (3, 4)
        assert ewl1.shape == (3, 4)

    def test_different_edge_labels_distinguish(self, triangle):
        vlabel = np.array([1, 1, 1])
        elabel_a = {(0, 1): 10, (1, 2): 10, (0, 2): 10}
        elabel_b = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        vwl_a, vwl_b, _, _ = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel, elabel_a, elabel_b, K=1, I=3)
        assert not np.all(vwl_a == vwl_b)

    @pytest.mark.parametrize("K", [0, 1, 2])
    def test_edge_k0_k1_k2(self, triangle, K):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        vwl1, vwl2, ewl1, ewl2 = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel, elabel, elabel, K=K, I=2)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)

    def test_uniform_edge_labels_match_node_only(self, triangle):
        """Uniform edge labels (all same) → should match node-only (no edge info)."""
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 0, (1, 2): 0, (0, 2): 0}
        vwl_edge, _, _, _ = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel, elabel, elabel, K=1, I=2)
        vwl_node, _ = hierarchical_triangular_wl(
            triangle, triangle, vlabel, vlabel, K=1, I=2)
        # Edge histories differ (present vs absent) but node histories should be identical
        assert np.all(vwl_edge == vwl_node)

    def test_edge_label_monotonicity(self, triangle):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        _, _, ewl1, ewl2 = hierarchical_triangular_wl_with_edges(
            triangle, triangle, vlabel, vlabel, elabel, elabel, K=1, I=3)
        for col in range(1, ewl1.shape[1]):
            assert np.all(ewl1[:, col] >= ewl1[:, col - 1])

    def test_non_consecutive_ids_with_edge_labels(self):
        G = nx.Graph()
        G.add_edges_from([(0, 5), (5, 10), (10, 0)])
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 5): 10, (5, 10): 20, (0, 10): 10}
        vwl1, vwl2, ewl1, ewl2 = hierarchical_triangular_wl_with_edges(
            G, G, vlabel, vlabel, elabel, elabel, K=1, I=2)
        assert np.all(vwl1 == vwl2)
        assert np.all(ewl1 == ewl2)

    def test_one_sided_edge_label_error(self, triangle):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        with pytest.raises(ValueError, match="Either both or neither"):
            hierarchical_triangular_wl_with_edges(
                triangle, triangle, vlabel, vlabel, elabel, None, K=1, I=2)

    def test_invalid_edge_key(self, triangle):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 99): 10}  # 99 not in triangle
        with pytest.raises(ValueError, match="not present"):
            hierarchical_triangular_wl_with_edges(
                triangle, triangle, vlabel, vlabel, elabel, elabel, K=1, I=2)


# ===========================================================================
# Unified interface tests
# ===========================================================================

class TestUnifiedInterface:
    def test_unified_no_el_key(self, triangle):
        vlabel = np.array([1, 2, 3])
        g_info = {}
        wl1, wl2 = hierarchical_triangular_wl_unified(
            g_info, triangle, triangle, vlabel, vlabel, K=1, I=2)
        assert np.all(wl1 == wl2)

    def test_unified_el_false(self, triangle):
        vlabel = np.array([1, 2, 3])
        g_info = {"el": False}
        wl1, wl2 = hierarchical_triangular_wl_unified(
            g_info, triangle, triangle, vlabel, vlabel, K=1, I=2)
        assert np.all(wl1 == wl2)

    def test_unified_el_true(self, triangle):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        g_info = {"el": True}
        wl1, wl2 = hierarchical_triangular_wl_unified(
            g_info, triangle, triangle, vlabel, vlabel, elabel, elabel, K=1, I=2)
        assert np.all(wl1 == wl2)

    @pytest.mark.parametrize("K", [0, 1])
    def test_unified_k0_k1(self, triangle, K):
        vlabel = np.array([1, 2, 3])
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 10}
        g_info = {"el": True}
        wl1, wl2 = hierarchical_triangular_wl_unified(
            g_info, triangle, triangle, vlabel, vlabel, elabel, elabel, K=K, I=2)
        assert np.all(wl1 == wl2)

    def test_unified_integration_with_example1(self, example1):
        """Real-sized graph through unified interface."""
        np.random.seed(42)
        vlabel = np.random.randint(0, 10, example1.number_of_nodes())
        g_info = {"el": False}
        wl1, wl2 = hierarchical_triangular_wl_unified(
            g_info, example1, example1, vlabel, vlabel.copy(), K=2, I=3)
        assert np.all(wl1 == wl2)
        assert wl1.shape == (26, 4)


# ===========================================================================
# Step-1 mapping tests
# ===========================================================================

class TestStep1Mappings:
    def test_csg_to_input_mapping(self, example1):
        H, cb, info = cyclic_schematic_graph(example1)
        csg_to_input = build_csg_to_input_mapping(H, cb)
        for csg_node, original_tuple in csg_to_input.items():
            if get_node_type(H, csg_node) == "cycle_basis":
                assert len(original_tuple) >= 3, \
                    f"Cycle basis {csg_node} should have >= 3 nodes"
            else:
                assert original_tuple == (csg_node,), \
                    f"Non-cycle node {csg_node} should map to itself"

    def test_input_to_csg_mapping_coverage(self, example1):
        H, cb, info = cyclic_schematic_graph(example1)
        input_to_csg = build_input_to_csg_mapping(H, cb, example1)
        for v in example1.nodes():
            assert v in input_to_csg, f"Node {v} missing from input_to_csg mapping"

    def test_multi_layer_mappings_consistent(self, example1):
        """Each CSG layer should have valid mappings from its lower graph."""
        from cyclic_schema import build_multilayer_csg_with_mappings
        layers, mappings = build_multilayer_csg_with_mappings(example1, K=2)
        assert len(layers) == 2
        assert len(mappings) == 2
        for k, (csg_to_lower, lower_to_csg) in enumerate(mappings):
            assert len(csg_to_lower) == layers[k][0].number_of_nodes()
            assert len(lower_to_csg) > 0


# ===========================================================================
# Forward aggregate unit tests
# ===========================================================================

class TestForwardAggregate:
    def test_isolated_node(self):
        """Node with no neighbors → AGG = (((label,),),)."""
        G = nx.Graph()
        G.add_node(0)
        nc = _precompute_neighbor_components(G)
        labels = {0: (5,)}
        agg = forward_aggregate(G, 0, labels, nc)
        # AGG(0) = ((l(v),),) = (((5,),),) — l(v) is itself a tuple
        assert agg == (((5,),),)

    def test_two_neighbors_connected(self, triangle):
        """Triangle: node with neighbors that are connected."""
        nc = _precompute_neighbor_components(triangle)
        labels = {0: (1,), 1: (2,), 2: (3,)}
        agg = forward_aggregate(triangle, 0, labels, nc)
        # AGG(0) = (((1,),), ((2,), (3,))) — two neighbors 1,2 in same component
        assert agg == (((1,),), ((2,), (3,)))

    def test_two_neighbors_disconnected(self):
        """Star: node with neighbors that are NOT connected (each isolated)."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2)])
        nc = _precompute_neighbor_components(G)
        labels = {0: (10,), 1: (20,), 2: (30,)}
        agg = forward_aggregate(G, 0, labels, nc)
        # AGG(0) = (((10,),), (20,), (30,)) — each neighbor isolated, labels are tuples
        assert agg == (((10,),), (20,), (30,))


# ===========================================================================
# Edge context unit tests
# ===========================================================================

class TestEdgeContext:
    def test_no_neighbors(self):
        G = nx.Graph()
        G.add_node(0)
        nc = _precompute_neighbor_components(G)
        ec = _compute_edge_context(0, G, {}, nc)
        assert ec == ()

    def test_triangle_uniform_edges(self, triangle):
        nc = _precompute_neighbor_components(triangle)
        elabel = {(0, 1): 10, (1, 2): 10, (0, 2): 10}
        ec0 = _compute_edge_context(0, triangle, elabel, nc)
        # Node 0 has neighbors {1, 2} in one component → ((10, 10),)
        assert ec0 == ((10, 10),)

    def test_triangle_distinct_edges(self, triangle):
        nc = _precompute_neighbor_components(triangle)
        elabel = {(0, 1): 10, (1, 2): 20, (0, 2): 30}
        ec0 = _compute_edge_context(0, triangle, elabel, nc)
        # Node 0 with neighbors {1,2} → ((10, 30),) — sorted edge labels
        assert ec0 == ((10, 30),)

    def test_star_with_edge_labels(self):
        """Star: center connected to two leaves (leaves are isolated)."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (0, 2)])
        nc = _precompute_neighbor_components(G)
        elabel = {(0, 1): 5, (0, 2): 15}
        ec0 = _compute_edge_context(0, G, elabel, nc)
        # Node 0: neighbors 1 and 2 are isolated → ((5,), (15,))
        # Components are sorted lexicographically
        assert ec0 == ((5,), (15,))


# ===========================================================================
# Edge label refresh unit tests
# ===========================================================================

class TestEdgeLabelRefresh:
    def test_basic_update(self, triangle):
        vlabel = {0: 10, 1: 20, 2: 30}
        elabel = {(0, 1): 1, (1, 2): 2, (0, 2): 3}
        new_e1, new_e2 = _update_elabel_from_dict(
            vlabel, vlabel, elabel, elabel)
        assert len(new_e1) == 3
        assert len(new_e2) == 3
        # All edges should have updated labels
        for e in elabel:
            assert new_e1[e] >= 0

    def test_update_identical_vlabels(self, triangle):
        """Same node labels → same edge label updates for all edges."""
        vlabel = {0: 5, 1: 5, 2: 5}
        elabel = {(0, 1): 1, (1, 2): 1, (0, 2): 1}
        new_e1, _ = _update_elabel_from_dict(
            vlabel, vlabel, elabel, elabel)
        # All edges should get the same new label
        vals = list(new_e1.values())
        assert len(set(vals)) == 1


# ===========================================================================
# Backward label tuple unit tests
# ===========================================================================

class TestBackwardLabelTuples:
    def test_no_higher_nodes(self):
        """No corresponding CSG nodes → just (l_G(v),)."""
        G = nx.Graph()
        G.add_node(0)
        lower_labels = {0: 42}
        higher_labels = {}
        lower_to_higher = {0: ()}
        result = _compute_lower_label_tuples(
            G, lower_labels, higher_labels, lower_to_higher)
        assert result[0] == (42,)

    def test_with_higher_labels(self, triangle):
        lower_labels = {0: 5, 1: 10, 2: 15}
        higher_labels = {"cb0": 100, "cb1": 200}
        lower_to_higher = {0: ("cb0", "cb1"), 1: ("cb0",), 2: ("cb1",)}
        result = _compute_lower_label_tuples(
            triangle, lower_labels, higher_labels, lower_to_higher)
        # Node 0: (5,) + sorted([100, 200]) = (5, 100, 200)
        assert result[0] == (5, 100, 200)
        # Node 1: (10,) + sorted([100]) = (10, 100)
        assert result[1] == (10, 100)

    def test_with_edge_context(self, triangle):
        lower_labels = {0: 5, 1: 10, 2: 15}
        higher_labels = {"cb0": 100}
        lower_to_higher = {0: ("cb0",), 1: (), 2: ()}
        edge_contexts = {0: ((10, 20),), 1: (), 2: ()}
        result = _compute_lower_label_tuples_with_edges(
            triangle, lower_labels, higher_labels,
            lower_to_higher, edge_contexts)
        # Node 0: (5, ((10, 20),), 100) — ec(v) is a tuple of tuples
        assert result[0] == (5, ((10, 20),), 100)


# ===========================================================================
# Multi-layer structural tests
# ===========================================================================

class TestMultiLayer:
    def test_two_disjoint_cycles(self):
        """Two disjoint triangles — each cycle is its own cluster."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
        labels = np.array([0, 1, 2, 0, 1, 2])
        wl1, wl2 = hierarchical_triangular_wl(
            G, G, labels, labels.copy(), K=1, I=2)
        assert np.all(wl1 == wl2)

    def test_k0_k1_same_initial(self, triangle):
        """K=0 and K=1 should both start with same initial labels."""
        vlabel = np.array([5, 10, 15])
        wl_k0, _ = hierarchical_triangular_wl(
            triangle, triangle, vlabel, vlabel, K=0, I=1)
        wl_k1, _ = hierarchical_triangular_wl(
            triangle, triangle, vlabel, vlabel, K=1, I=1)
        # Both should have same initial labels (column 0)
        assert np.all(wl_k0[:, 0] == wl_k1[:, 0])


# ===========================================================================
# Import verification
# ===========================================================================

class TestPackageExports:
    def test_cyclic_schema_imports(self):
        """Verify package-level exports."""
        from cyclic_schema import (
            hierarchical_triangular_wl,
            hierarchical_triangular_wl_with_edges,
            hierarchical_triangular_wl_unified,
            cyclic_schematic_graph,
            build_csg_to_input_mapping,
            build_input_to_csg_mapping,
        )
        assert hierarchical_triangular_wl is not None
        assert hierarchical_triangular_wl_with_edges is not None
        assert hierarchical_triangular_wl_unified is not None
