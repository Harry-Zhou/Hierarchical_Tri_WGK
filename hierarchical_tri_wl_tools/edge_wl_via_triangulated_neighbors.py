import numpy as np
from functools import reduce
import os
import sys
cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)
from .local_structure_edge_wl_tools import update_elabel
# from .local_structure_wl_tools import compress_and_relabel_vlabel
from .node_wl_via_triangulated_neighbors import node_wl_test_triangulated_neighbors

def propagate_via_triangulated_neighboring_edges(vtx_triangulated_neighbors, vlabel_np, elabel_dict):
    vlabel_collection = {}
    for vtx, triangulated_neighbors in vtx_triangulated_neighbors.items():
        vlabels = []
        for tri_nei in triangulated_neighbors: # tri_nei是由中心节点，及其一组三角邻居诱导出的子图
            vlabels_tri = [elabel_dict[(v0, v1) if v0 <= v1 else (v1, v0)] for v0, v1 in tri_nei.edges()]
            # [elabel_dict[(vtx, vn) if vtx <= vn else (vn, vtx)] for vn in tri_nei]
            vlabels_tri.sort()
            vlabels.append(tuple(vlabels_tri))
        vlabels.sort()
        vlabel_collection[vtx] = tuple([(vlabel_np[vtx], ), ] + vlabels)
    return vlabel_collection

def gen_compressed_vlabel(vlabel_collection, total_label_collection, next_label):
    compressed_vlabel_np = np.zeros(
        (len(vlabel_collection), ), 
        dtype = np.int32
    )
    for vtx, vlabel_tuple in vlabel_collection.items():
        lidx = total_label_collection.index(vlabel_tuple)
        compressed_vlabel_np[vtx] = lidx + next_label
    return compressed_vlabel_np
 
def compress_and_relabel_vlabel(vlabel_collection1, vlabel_collection2):
    total_label_collection = list(vlabel_collection1.values())
    total_label_collection.extend(vlabel_collection2.values())
    total_label_collection = list(set(total_label_collection))
    total_label_collection.sort()
    next_label = max([lc[0][0] for lc in total_label_collection]) + 1
    node_compressed_label_np1 = gen_compressed_vlabel(
        vlabel_collection1, total_label_collection, next_label
    )
    node_compressed_label_np2 = gen_compressed_vlabel(
        vlabel_collection2, total_label_collection, next_label
    )
    return node_compressed_label_np1, node_compressed_label_np2

def edge_wl_test_triangulated_neighbors(
    vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
    vlabel_np1, vlabel_np2, edges1, edges2, 
    elabel_np1, elabel_np2, niter
):
    vwl_np1 = np.zeros((len(vtx_triangulated_neighbors1), niter + 1), dtype = np.int32)
    vwl_np2 = np.zeros((len(vtx_triangulated_neighbors2), niter + 1), dtype = np.int32)
    vwl_np1[:, 0] = vlabel_np1 # copy
    vwl_np2[:, 0] = vlabel_np2 # copy
    ewl_np1 = np.zeros((len(edges1), niter + 1), dtype = np.int32)
    ewl_np2 = np.zeros((len(edges2), niter + 1), dtype = np.int32)
    ewl_np1[:, 0], ewl_np2[:, 0] = elabel_np1, elabel_np2
    compressed_elabel1, compressed_elabel2 = dict(zip(edges1, elabel_np1)), dict(zip(edges2, elabel_np2))
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_via_triangulated_neighboring_edges(
            vtx_triangulated_neighbors1, vwl_np1[:, i - 1], compressed_elabel1
        )
        vlabel_collection2 = propagate_via_triangulated_neighboring_edges(
            vtx_triangulated_neighbors2, vwl_np2[:, i - 1], compressed_elabel2
        )
        temp_vlabel_np1, temp_vlabel_np2 = compress_and_relabel_vlabel(
            vlabel_collection1, vlabel_collection2
        )
        vwl_np1[:, i] = temp_vlabel_np1
        vwl_np2[:, i] = temp_vlabel_np2
        compressed_elabel1, compressed_elabel2 = update_elabel(
            temp_vlabel_np1, temp_vlabel_np2, 
            compressed_elabel1, compressed_elabel2
        )
        for eidx1, e1 in enumerate(edges1):
            ewl_np1[eidx1, i] = compressed_elabel1[e1]
        for eidx2, e2 in enumerate(edges2):
            ewl_np2[eidx2, i] = compressed_elabel2[e2]
        
    return vwl_np1, vwl_np2, ewl_np1, ewl_np2
