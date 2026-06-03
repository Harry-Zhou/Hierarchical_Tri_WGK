from functools import reduce
import numpy as np
import networkx as nx
import os
import sys
cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)
# from local_structure_wl_tools import compress_and_relabel_vlabel

def propagate_via_triangulated_neighbors(vtx_triangulated_neighbors, vlabel_np):
    vlabel_collection = {}
    for v, tri_neighbors in vtx_triangulated_neighbors.items():
        vlabel_collection_value = [(vlabel_np[v], ), ]
        vlabels = []
        for tri_nei in tri_neighbors: # tri_nei是由中心节点，及其一组三角邻居诱导出的子图
            vlabel_trin = [vlabel_np[trin] for trin in tri_nei.nodes() if trin != v]
            vlabel_trin.sort()
            vlabels.append(tuple(vlabel_trin))
        vlabels.sort()
        vlabel_collection_value.extend(vlabels)
        vlabel_collection[v] = tuple(vlabel_collection_value)
    return vlabel_collection # vtx: ((l, ...), ...)

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
    # total_label_collection = list(set(total_label_collection))
    total_label_collection.sort()
    next_label = max([lc[0][0] for lc in total_label_collection]) + 1
    node_compressed_label_np1 = gen_compressed_vlabel(
        vlabel_collection1, total_label_collection, next_label
    )
    node_compressed_label_np2 = gen_compressed_vlabel(
        vlabel_collection2, total_label_collection, next_label
    )
    return node_compressed_label_np1, node_compressed_label_np2

def node_wl_test_triangulated_neighbors(
    vtx_triangulated_neighbors1, vtx_triangulated_neighbors2, 
    vlabel_np1, vlabel_np2, niter
):
    wl_np1 = np.zeros((len(vtx_triangulated_neighbors1), niter + 1), dtype = np.int32)
    wl_np2 = np.zeros((len(vtx_triangulated_neighbors2), niter + 1), dtype = np.int32)
    wl_np1[:, 0] = vlabel_np1
    wl_np2[:, 0] = vlabel_np2
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_via_triangulated_neighbors(
            vtx_triangulated_neighbors1, wl_np1[:, i - 1]
        )
        vlabel_collection2 = propagate_via_triangulated_neighbors(
            vtx_triangulated_neighbors2, wl_np2[:, i - 1]
        )
        temp_vlabel_np1, temp_vlabel_np2 = compress_and_relabel_vlabel(
            vlabel_collection1, vlabel_collection2
        )
        wl_np1[:, i] = temp_vlabel_np1
        wl_np2[:, i] = temp_vlabel_np2
    return wl_np1, wl_np2
