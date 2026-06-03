import numpy as np
from .edge_wl_via_cycle_complexes import edge_wl_test_cycle_complexes
from .edge_wl_via_triangulated_neighbors import edge_wl_test_triangulated_neighbors
from .node_wl_via_triangulated_neighbors import node_wl_test_triangulated_neighbors
from .node_wl_via_cycle_complexes import node_wl_test_cycle_complexes

def topo_aware_node_wl_test(
    vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
    vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
    vlabel1, vlabel2, 
    niter_tn, niter_hcc
):
    tn_vwl1, tn_vwl2 = node_wl_test_triangulated_neighbors(
        vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
        vlabel1, vlabel2, niter_tn
    )
    if niter_hcc > 0:
        hcc_vwl1, hcc_vwl2 = node_wl_test_cycle_complexes(
            vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
            tn_vwl1[:, -1], tn_vwl2[:, -1], niter_hcc
        )
        return np.concatenate([tn_vwl1, hcc_vwl1], axis = 1), np.concatenate([tn_vwl2, hcc_vwl2], axis = 1)
    else:
        return tn_vwl1, tn_vwl2

def topo_aware_edge_wl_test(
    vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
    vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
    vlabel1, vlabel2, 
    edges1, edges2, elabel_np1, elabel_np2, 
    niter_tn, niter_hcc
):
    tn_vwl1, tn_vwl2, _, _ = edge_wl_test_triangulated_neighbors(
        vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
        vlabel1, vlabel2, edges1, edges2, elabel_np1, elabel_np2, 
        niter_tn
    )
    if niter_hcc > 0:
        hcc_vwl1, hcc_vwl2, _, _ = edge_wl_test_cycle_complexes(
            vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
            tn_vwl1[:, -1], tn_vwl2[:, -1], edges1, edges2, elabel_np1, elabel_np2, 
            niter_hcc
        )
        return np.concatenate([tn_vwl1, hcc_vwl1], axis = 1), np.concatenate([tn_vwl2, hcc_vwl2], axis = 1)
    else:
        return tn_vwl1, tn_vwl2

def topo_aware_wl_test(
    g_info, niter_tn, niter_hcc, 
    vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
    vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
    vlabel1, vlabel2, 
    edges1, edges2, elabel1, elabel2, 
):
    if g_info['el']:
        assert type(elabel1) == np.ndarray and type(elabel2) == np.ndarray
        assert elabel1.shape[0] == len(edges1) and elabel2.shape[0] == len(edges2)
        vwl_np1, vwl_np2 = topo_aware_edge_wl_test(
            vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
            vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
            vlabel1, vlabel2, 
            edges1, edges2, elabel1, elabel2, 
            niter_tn, niter_hcc
        )
    else:
        vwl_np1, vwl_np2 = topo_aware_node_wl_test(
            vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
            vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
            vlabel1, vlabel2, 
            niter_tn, niter_hcc
        )
    return vwl_np1, vwl_np2