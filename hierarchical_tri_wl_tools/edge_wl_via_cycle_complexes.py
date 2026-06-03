import numpy as np
from functools import reduce
import os
import sys
cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)
from .local_structure_edge_wl_tools import update_elabel
# from .local_structure_wl_tools import compress_and_relabel_vlabel
from .node_wl_via_cycle_complexes import node_wl_test_cycle_complexes

def propagate_via_cycle_complexes(vtx_hierarchical_cycle_contexts, vlabel_np, elabel_dict):
    """
    vtx_hierarchical_cycle_contexts: {vtx: [[cycle_complex1, cycle_complex2, ...], ...]}
    """
    # # 原始实现（恢复为主逻辑）：三层嵌套循环，行为与历史版本一致
    # vlabel_collection = {}
    # for vtx, hierarchical_cycle_contexts in vtx_hierarchical_cycle_contexts.items():
    #     vlabel_collection_value = [[[vlabel_np[vtx],], ], ] # [layer[cycle[...], ], ]
    #     vlabel_hccs = []
    #     for cycle_contexts in hierarchical_cycle_contexts:
    #         vlabel_ccs = []
    #         for cycle_context in cycle_contexts:
    #             vlabel_cc = [elabel_dict[(v0, v1) if v0 <= v1 else (v1, v0)] for v0, v1 in cycle_context.edges()]
    #             vlabel_cc.sort()
    #             vlabel_ccs.append(tuple(vlabel_cc))
    #         vlabel_ccs.sort()
    #         vlabel_hccs.append(tuple(vlabel_ccs))
    #     vlabel_collection_value.extend(vlabel_hccs)
    #     vlabel_collection[vtx] = tuple(vlabel_collection_value)
    # return vlabel_collection

    # Micro-optimizations similar to node version:
    # - bind local names to reduce attribute/name lookup cost
    # - use explicit list comprehensions and built-in sorted for clarity and speed
    vlabel_collection = {}
    elabel_loc = elabel_dict  # now expected to be mapping edge-index -> label
    vlabel_loc = vlabel_np
    for vtx, hierarchical_cycle_contexts in vtx_hierarchical_cycle_contexts.items():
        vlabel_collection_value = [[[vlabel_loc[vtx]], ], ]  # keep same nesting
    
        vlabel_hccs = []
        for cycle_contexts in hierarchical_cycle_contexts:
            # build list of sorted edge-label tuples for each cycle
            # cycle_context is now a tuple/list of edge indices (into the graph's edges list)
            vlabel_ccs_list = [
                tuple(sorted([int(elabel_loc[eidx]) for eidx in cycle_context]))
                for cycle_context in cycle_contexts
            ]
            vlabel_ccs = tuple(sorted(vlabel_ccs_list))
            vlabel_hccs.append(vlabel_ccs)
    
        vlabel_collection_value.extend(vlabel_hccs)
        vlabel_collection[vtx] = tuple(vlabel_collection_value)
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
    # total_label_collection = list(set(total_label_collection))
    total_label_collection.sort()
    next_label = max([lc[0][0][0] for lc in total_label_collection]) + 1
    node_compressed_label_np1 = gen_compressed_vlabel(
        vlabel_collection1, total_label_collection, next_label
    )
    node_compressed_label_np2 = gen_compressed_vlabel(
        vlabel_collection2, total_label_collection, next_label
    )
    return node_compressed_label_np1, node_compressed_label_np2

def collect_vlabels_edge(vlabel_np, elabel_dict):
    elabel_collection = {}
    for e, l in elabel_dict.items():
        ev_label_lhs, ev_label_rhs = vlabel_np[e[0]], vlabel_np[e[1]]
        if ev_label_lhs <= ev_label_rhs: # 边的两个端点的标签排序
            elabel_collection[e] = (l, ev_label_lhs, ev_label_rhs)
        else:
            elabel_collection[e] = (l, ev_label_rhs, ev_label_lhs)
    return elabel_collection

def gen_compressed_elabel(elabel_collection, total_elabel_list, next_label):
    compressed_elabel = {}
    for e, elabel_tuple in elabel_collection.items():
        lidx = total_elabel_list.index(elabel_tuple)
        compressed_elabel[e] = lidx + next_label
    return compressed_elabel

def compress_and_relabel_elabel(elabel_collection1, elabel_collection2):
    total_elabel_list = list(elabel_collection1.values())
    total_elabel_list.extend(elabel_collection2.values())
    total_elabel_list = list(set(total_elabel_list))
    total_elabel_list.sort()
    next_label = max([lc[0] for lc in total_elabel_list]) + 1
    compressed_elabel1 = gen_compressed_elabel(elabel_collection1, total_elabel_list, next_label)
    compressed_elabel2 = gen_compressed_elabel(elabel_collection2, total_elabel_list, next_label)
    return compressed_elabel1, compressed_elabel2

def update_elabel(
    vlabel_np1, vlabel_np2, 
    elabel_dict1, elabel_dict2
):
    # edge label collection
    elabel_collection1 = collect_vlabels_edge(vlabel_np1, elabel_dict1)
    elabel_collection2 = collect_vlabels_edge(vlabel_np2, elabel_dict2)
    # compression and relabeling
    return compress_and_relabel_elabel(elabel_collection1, elabel_collection2)

def edge_wl_test_cycle_complexes(
    vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, 
    vlabel_np1, vlabel_np2, edges1, edges2, 
    elabel_np1, elabel_np2, niter
):
    """
    vtx_hierarchical_cycle_contexts: {vtx: [[cycle_complex1, cycle_complex2, ...], ...]}
    """
    vwl_np1 = np.zeros((len(vtx_hierarchical_cycle_contexts1), niter + 1), dtype = np.int32)
    vwl_np2 = np.zeros((len(vtx_hierarchical_cycle_contexts2), niter + 1), dtype = np.int32)
    vwl_np1[:, 0] = vlabel_np1 # copy
    vwl_np2[:, 0] = vlabel_np2 # copy
    ewl_np1 = np.zeros((len(edges1), niter + 1), dtype = np.int32)
    ewl_np2 = np.zeros((len(edges2), niter + 1), dtype = np.int32)
    ewl_np1[:, 0], ewl_np2[:, 0] = elabel_np1, elabel_np2
    # compressed_elabelX maps edge-tuple -> label (used by update_elabel)
    compressed_elabel1, compressed_elabel2 = dict(zip(edges1, elabel_np1)), dict(zip(edges2, elabel_np2))
    # elabel_idx_map maps edge-index -> label (used by propagate_via_cycle_complexes)
    elabel_idx_map1 = dict(enumerate(elabel_np1))
    elabel_idx_map2 = dict(enumerate(elabel_np2))
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_via_cycle_complexes(
            vtx_hierarchical_cycle_contexts1, vwl_np1[:, i - 1], elabel_idx_map1
        )
        vlabel_collection2 = propagate_via_cycle_complexes(
            vtx_hierarchical_cycle_contexts2, vwl_np2[:, i - 1], elabel_idx_map2
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
        # rebuild index->label maps from updated compressed_elabel (which maps edge-tuple -> label)
        elabel_idx_map1 = {idx: compressed_elabel1[edges1[idx]] for idx in range(len(edges1))}
        elabel_idx_map2 = {idx: compressed_elabel2[edges2[idx]] for idx in range(len(edges2))}
        for eidx1, e1 in enumerate(edges1):
            ewl_np1[eidx1, i] = compressed_elabel1[e1]
        for eidx2, e2 in enumerate(edges2):
            ewl_np2[eidx2, i] = compressed_elabel2[e2]
        
    return vwl_np1, vwl_np2, ewl_np1, ewl_np2
