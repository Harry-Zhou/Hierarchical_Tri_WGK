"""
Full matrix test: TopoWassersteinGraphKernel with edge-label dispatch.

Validates:
1. Edge-label-aware kernel distinguishes graphs differing only in edge labels
2. Node-only kernel treats those graphs as identical
3. Edge-label kernel produces different results from node-only kernel
4. Unified dispatch works for both has_el=True and has_el=False cases
5. L=0 properly dispatches to triangulated neighbors WL for both node/edge cases
6. Full matrix comparison: edge-labeled vs non-edge-labeled outputs
7. g_info passage through from dataset processing to WL test
"""
import sys
import numpy as np
import networkx as nx

sys.path.insert(0, '.')
from utils import precomp_node_neighs
from our_experiments.cycle_complex_wgk.topo_wasserstein_graph_kernel import (
    TopoWassersteinGraphKernel,
)
from cyclic_schema.htn_wl import (
    hierarchical_triangular_wl_unified,
    hierarchical_triangular_wl,
    hierarchical_triangular_wl_with_edges,
)


def _make_path3_graphs():
    """Return 3 path-graphs with controlled edge labels."""
    graphs = []
    vlabels = []
    edges_list = []
    elabels_list = []
    deg_distr_list = []

    for _ in range(3):
        G = nx.path_graph(3)
        graphs.append(G)
        vl = np.array([G.degree(v) for v in G.nodes()])
        vlabels.append(vl)
        e = [(u, v) if u <= v else (v, u) for u, v in G.edges()]
        edges_list.append(e)
        deg = vl / vl.sum().astype(np.float32)
        deg_distr_list.append(deg)

    # G0 and G2: all edge labels = 1; G1: edge labels differ (1, 2)
    elabels_list.append(np.array([1, 1]))
    elabels_list.append(np.array([1, 2]))
    elabels_list.append(np.array([1, 1]))

    return (graphs, vlabels, edges_list, elabels_list, deg_distr_list)


def test_unified_wl_dispatch():
    """Test that hierarchical_triangular_wl_unified dispatches correctly."""
    G1 = nx.path_graph(3)
    G2 = nx.path_graph(3)
    vl1 = np.array([1, 2, 1], dtype=np.int32)
    vl2 = vl1.copy()

    # === Node-only (L=1) ===
    g_info_noel = {'el': False, 'nl': True}
    wl_n, _ = hierarchical_triangular_wl_unified(
        g_info_noel, G1, G2, vl1, vl2, None, None, L=1, I=2)
    assert wl_n.shape == (3, 3), f"Expected (3,3), got {wl_n.shape}"

    # === Node-only (L=0) ===
    wl_n0, _ = hierarchical_triangular_wl_unified(
        g_info_noel, G1, G2, vl1, vl2, None, None, L=0, I=2)
    assert wl_n0.shape == (3, 3), f"Expected (3,3), got {wl_n0.shape}"

    # === Edge-label (L=1) ===
    g_info_el = {'el': True, 'nl': True}
    ed1 = {(0, 1): 5, (1, 2): 10}
    ed2 = {(0, 1): 5, (1, 2): 10}
    wl_e, _ = hierarchical_triangular_wl_unified(
        g_info_el, G1, G2, vl1, vl2, ed1, ed2, L=1, I=2)
    assert wl_e.shape == (3, 3), f"Expected (3,3), got {wl_e.shape}"

    # === Edge-label (L=0) ===
    wl_e0, _ = hierarchical_triangular_wl_unified(
        g_info_el, G1, G2, vl1, vl2, ed1, ed2, L=0, I=2)
    assert wl_e0.shape == (3, 3), f"Expected (3,3), got {wl_e0.shape}"

    # Edge-label output should differ from node-only
    assert not np.all(wl_e == wl_n), \
        "Edge-label WL should differ from node-only WL"

    print("[PASS] Unified WL dispatch correct for all (el, L) combinations")


def test_kernel_node_only():
    """Test kernel with node-only (has_el=False) datasets."""
    graphs, vlabels, edges_list, elabels_list, deg_distr_list = _make_path3_graphs()

    # Test with n_csg_layers=0 (L=0, TN path)
    g_info = {'el': False, 'nl': True}
    kernel = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=0, wl_normalized=False,
    )
    kernel.fit(g_info, graphs, vlabels, edges_list, elabels_list, deg_distr_list)
    ot, wl, _ = kernel.transform()
    # All 3 graphs are identical (same path structure, same node labels)
    assert ot.max() < 1e-8, \
        "Node-only L=0: identical graphs should have zero OT distance"
    print("[PASS] Kernel node-only L=0 produces correct full matrix")

    # Test with n_csg_layers=1 (L=1, CSG path)
    kernel2 = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=1, wl_normalized=False,
    )
    kernel2.fit(g_info, graphs, vlabels, edges_list, elabels_list, deg_distr_list)
    ot2, wl2, _ = kernel2.transform()
    assert ot2.max() < 1e-8, \
        "Node-only L=1: identical graphs should have zero OT distance"
    print("[PASS] Kernel node-only L=1 produces correct full matrix")


def test_kernel_edge_label():
    """Test kernel with edge-label (has_el=True) datasets."""
    graphs, vlabels, edges_list, elabels_list, deg_distr_list = _make_path3_graphs()

    g_info = {'el': True, 'nl': True}

    # --- Test with n_csg_layers=0 (L=0, TN edge WL) ---
    kernel = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=0, wl_normalized=False,
    )
    kernel.fit(g_info, graphs, vlabels, edges_list, elabels_list, deg_distr_list)
    ot_el, wl_el, _ = kernel.transform()

    # G0 vs G1: different edge labels -> different
    assert ot_el[0, 0] != ot_el[0, 1] or wl_el[0, 0] != wl_el[0, 1], \
        "Edge-label L=0 kernel must distinguish G0 from G1"
    # G0 vs G2: same edge labels -> same
    assert ot_el[0, 2] < 1e-8, \
        "G0 and G2 (same edge labels) should have zero OT distance (L=0)"
    print("[PASS] Kernel edge-label L=0 distinguishes by edge labels")

    # --- Test with n_csg_layers=1 (L=1, CSG edge WL) ---
    kernel2 = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=1, wl_normalized=False,
    )
    kernel2.fit(g_info, graphs, vlabels, edges_list, elabels_list, deg_distr_list)
    ot2, wl2, _ = kernel2.transform()

    assert ot2[0, 0] != ot2[0, 1] or wl2[0, 0] != wl2[0, 1], \
        "Edge-label L=1 kernel must distinguish G0 from G1"
    assert ot2[0, 2] < 1e-8, \
        "G0 and G2 (same edge labels) should have zero OT distance (L=1)"
    print("[PASS] Kernel edge-label L=1 distinguishes by edge labels")


def test_full_matrix_comparison():
    """
    Full matrix comparison: edge-labeled vs non-edge-labeled outputs.
    Validates that the same dataset produces different kernel matrices
    when edge labels are present vs absent.
    """
    graphs, vlabels, edges_list, elabels_list, deg_distr_list = _make_path3_graphs()

    # Edge-label kernel
    g_info_el = {'el': True, 'nl': True}
    kernel_el = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=1, wl_normalized=False,
    )
    kernel_el.fit(g_info_el, graphs, vlabels, edges_list, elabels_list, deg_distr_list)
    ot_el, wl_el, _ = kernel_el.transform()

    # Same graphs, node-only
    g_info_noel = {'el': False, 'nl': True}
    kernel_noel = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=1, wl_normalized=False,
    )
    kernel_noel.fit(g_info_noel, graphs, vlabels, edges_list, elabels_list, deg_distr_list)
    ot_noel, wl_noel, _ = kernel_noel.transform()

    # The full matrices should differ
    diff_ot = np.abs(ot_el - ot_noel).max()
    assert diff_ot > 1e-6, \
        "Full OT matrix: edge-label should differ from node-only"
    print(f"[PASS] Full matrix OT difference: {diff_ot:.6f} > 1e-6")

    diff_wl = np.abs(wl_el - wl_noel).max()
    assert diff_wl > 1e-6, \
        "Full WL matrix: edge-label should differ from node-only"
    print(f"[PASS] Full matrix WL difference: {diff_wl:.6f} > 1e-6")


def test_g_info_pass_through():
    """
    Verify that g_info is properly passed through and affects the WL dispatch.
    The kernel should use different WL algorithms based on g_info['el'].
    """
    G1 = nx.path_graph(3)
    G2 = nx.path_graph(3)
    G3 = nx.path_graph(3)
    graphs = [G1, G2, G3]

    vl = np.array([1, 2, 1], dtype=np.int32)
    vlabels = [vl.copy() for _ in range(3)]

    edges = [(0, 1), (1, 2)]
    edges_list = [edges, edges, edges]

    elabels = [np.array([1, 1]), np.array([1, 1]), np.array([1, 1])]

    deg = vl / vl.sum().astype(np.float32)
    deg_distr_list = [deg, deg, deg]

    # Test with g_info['el']=True (all 3 graphs have same edge labels)
    g_info_el = {'el': True, 'nl': True}
    kernel_el = TopoWassersteinGraphKernel(
        n_wl_iters=2, n_csg_layers=1, wl_normalized=False,
    )
    kernel_el.fit(g_info_el, graphs, vlabels, edges_list, elabels, deg_distr_list)
    ot_el, wl_el, _ = kernel_el.transform()

    assert ot_el.max() < 1e-8, \
        "All graphs same (including edge labels): OT should be zero"
    print("[PASS] g_info pass-through: edge-label mode works for identical graphs")


def test_wl_label_decomposition():
    """
    Test that the unified WL decomposes correctly by comparing:
    - hierarchical_triangular_wl with el=None => same as direct call
    - hierarchical_triangular_wl_with_edges with el provided => edge-aware
    """
    G1 = nx.path_graph(3)
    G2 = nx.path_graph(3)
    vl = np.array([1, 2, 1], dtype=np.int32)

    g_info = {'el': False, 'nl': True}
    wl_unified, _ = hierarchical_triangular_wl_unified(
        g_info, G1, G2, vl, vl.copy(), None, None, L=1, I=2)
    wl_direct, _ = hierarchical_triangular_wl(
        G1, G2, vl, vl.copy(), L=1, I=2)

    assert np.all(wl_unified == wl_direct), \
        "Unified (el=False) should match direct call"
    print("[PASS] Unified WL decomposition matches direct WL call")


def run_all_tests():
    test_unified_wl_dispatch()
    test_kernel_node_only()
    test_kernel_edge_label()
    test_full_matrix_comparison()
    test_g_info_pass_through()
    test_wl_label_decomposition()
    print("\n=== All edge-label integration tests passed! ===")


if __name__ == '__main__':
    run_all_tests()
