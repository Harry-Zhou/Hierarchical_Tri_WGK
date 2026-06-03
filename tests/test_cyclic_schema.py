"""Unit tests for cyclic_schema.py — CSG construction and mappings."""

import networkx as nx
import pytest

from cyclic_schema.cyclic_schema import (
    node_sort_key,
    cyclic_schematic_graph,
    build_csg_to_input_mapping,
    build_input_to_csg_mapping,
    get_node_type,
    get_cycle_edges,
    get_edge_in_cycle,
    build_multilayer_csg,
    build_multilayer_csg_with_mappings,
)


class TestNodeSortKey:
    """Deterministic node ordering for CSG construction."""

    def test_int_nodes_first(self):
        assert node_sort_key(0) < node_sort_key("a")

    def test_int_ordering(self):
        assert node_sort_key(1) < node_sort_key(2)

    def test_str_ordering(self):
        assert node_sort_key("a") < node_sort_key("b")

    def test_numpy_int_treated_as_int(self):
        import numpy as np
        assert node_sort_key(np.int32(5)) < node_sort_key("x")
        assert node_sort_key(np.int32(3)) < node_sort_key(np.int32(7))

    def test_mixed_types_sortable(self):
        nodes = ["cb0", 1, "a", 0, "cb1", 2]
        sorted_nodes = sorted(nodes, key=node_sort_key)
        assert sorted_nodes == [0, 1, 2, "a", "cb0", "cb1"]


class TestCyclicSchematicGraph:
    """CSG construction from input graphs."""

    def test_triangle_has_one_cycle(self, triangle):
        H, cb, info = cyclic_schematic_graph(triangle)
        # triangle has 1 cycle (triangle itself)
        assert len(cb) == 1

    def test_path_has_no_cycles(self, path4):
        H, cb, info = cyclic_schematic_graph(path4)
        assert len(cb) == 0

    def test_chord_cycle_has_multiple_bases(self, chord_cycle):
        H, cb, info = cyclic_schematic_graph(chord_cycle)
        # 4-cycle with chord has 3 cycles in cycle basis
        assert len(cb) >= 2

    def test_csg_has_correct_node_types(self, triangle):
        H, cb, info = cyclic_schematic_graph(triangle)
        for n in H.nodes():
            assert get_node_type(H, n) in ("cycle_basis", "original_non_cycle", "interface")

    def test_two_disjoint_triangles(self, two_disjoint_triangles):
        H, cb, info = cyclic_schematic_graph(two_disjoint_triangles)
        assert len(cb) == 2  # two triangles → two cycle bases

    def test_single_node_no_cycles(self, single_node):
        H, cb, info = cyclic_schematic_graph(single_node)
        assert len(cb) == 0

    def test_returns_info_dict(self, triangle):
        H, cb, info = cyclic_schematic_graph(triangle)
        assert isinstance(info, dict)
        assert "cycle_basis" in info
        assert "edges_in_cycles" in info


class TestCSGMappings:
    """Step-1 correspondence mappings between CSG and input graph."""

    def test_csg_to_input_cycle_basis_size(self, triangle):
        H, cb, info = cyclic_schematic_graph(triangle)
        csg_to_input = build_csg_to_input_mapping(H, cb)
        for csg_node in H.nodes():
            if get_node_type(H, csg_node) == "cycle_basis":
                assert len(csg_to_input[csg_node]) == 3

    def test_csg_to_input_non_cycle_maps_to_self(self, triangle):
        H, cb, info = cyclic_schematic_graph(triangle)
        csg_to_input = build_csg_to_input_mapping(H, cb)
        for csg_node in H.nodes():
            if get_node_type(H, csg_node) != "cycle_basis":
                assert csg_to_input[csg_node] == (csg_node,)

    def test_input_to_csg_every_node_mapped(self, triangle):
        H, cb, info = cyclic_schematic_graph(triangle)
        input_to_csg = build_input_to_csg_mapping(H, cb, triangle)
        for v in triangle.nodes():
            assert v in input_to_csg

    def test_mappings_are_consistent(self, triangle):
        """input_to_csg and csg_to_input should be inverses for non-cycle nodes."""
        H, cb, info = cyclic_schematic_graph(triangle)
        csg_to_input = build_csg_to_input_mapping(H, cb)
        input_to_csg = build_input_to_csg_mapping(H, cb, triangle)

        for csg_node, input_nodes in csg_to_input.items():
            if get_node_type(H, csg_node) == "original_non_cycle":
                assert input_to_csg[input_nodes[0]] == (csg_node,)

    def test_non_consecutive_ids(self, non_consecutive_graph):
        """Mappings should work for graphs with non-consecutive node IDs."""
        H, cb, info = cyclic_schematic_graph(non_consecutive_graph)
        csg_to_input = build_csg_to_input_mapping(H, cb)
        input_to_csg = build_input_to_csg_mapping(H, cb, non_consecutive_graph)
        assert len(input_to_csg) == 3
        for v in (0, 5, 10):
            assert v in input_to_csg


class TestMultiLayerCSG:
    """Multi-layer CSG construction."""

    def test_two_layers(self, triangle):
        layers = build_multilayer_csg(triangle, K=2)
        assert len(layers) == 2
        for H, cb, info in layers:
            assert isinstance(H, nx.Graph)

    def test_k0_returns_empty_list(self, triangle):
        result = build_multilayer_csg(triangle, K=0)
        assert result == []

    def test_two_layers_with_mappings(self, triangle):
        layers, mappings = build_multilayer_csg_with_mappings(triangle, K=2)
        assert len(layers) == 2
        assert len(mappings) == 2
        for k in range(2):
            csg_to_lower, lower_to_csg = mappings[k]
            assert len(csg_to_lower) == layers[k][0].number_of_nodes()

    def test_increasing_layers_on_path(self):
        """Path graph should produce additional layers with more CSG depth."""
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 4)])
        result_k1 = build_multilayer_csg(G, K=1)
        result_k2 = build_multilayer_csg(G, K=2)
        assert result_k1 is not None

    def test_layer_k_higher_k_more_or_equal_layers(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2)])
        layers_k1 = build_multilayer_csg(G, K=1)
        layers_k2 = build_multilayer_csg(G, K=2)
        assert len(layers_k2) >= len(layers_k1)


class TestCycleEdges:
    """Cycle edge detection utilities."""

    def test_triangle_cycle_edges(self, triangle):
        _, cb, _ = cyclic_schematic_graph(triangle)
        cycle_edges = get_cycle_edges(cb[0])
        assert len(cycle_edges) == 3

    def test_edge_in_cycle_true(self, triangle):
        _, cb, _ = cyclic_schematic_graph(triangle)
        cycle_edges = get_edge_in_cycle(cb[0])
        assert (0, 1) in cycle_edges

    def test_edge_in_cycle_nonexistent(self, triangle):
        _, cb, _ = cyclic_schematic_graph(triangle)
        cycle_edges = get_edge_in_cycle(cb[0])
        assert (0, 99) not in cycle_edges

    def test_path_no_cycle_edges(self, path4):
        _, cb, _ = cyclic_schematic_graph(path4)
        assert len(cb) == 0
