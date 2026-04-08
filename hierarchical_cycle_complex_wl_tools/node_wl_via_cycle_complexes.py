from functools import reduce
import numpy as np
import networkx as nx
import os
import sys
cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)
# from local_structure_wl_tools import compress_and_relabel_vlabel

def propagate_via_cycle_complexes(vtx_hierarchical_cycle_contexts, vlabel_np):
    """
    vtx_hierarchical_cycle_contexts: {vtx: [[cycle_complex1, cycle_complex2, ...], ...]}
    """
    # # 原始实现（保留为主逻辑）：三层嵌套循环，易于理解且与历史实现一致。
    # vlabel_collection = {}
    # for vtx, hierarchical_cycle_contexts in vtx_hierarchical_cycle_contexts.items():
    #     vlabel_collection_value = [[[vlabel_np[vtx],], ], ] # [layer[cycle[...], ], ]
    #     vlabel_hccs = []
    #     for cycle_contexts in hierarchical_cycle_contexts:
    #         vlabel_ccs = []
    #         for cycle_context in cycle_contexts:
    #             vlabel_cc = [vlabel_np[ccvtx] for ccvtx in cycle_context.nodes() if ccvtx != vtx]
    #             vlabel_cc.sort()
    #             vlabel_ccs.append(tuple(vlabel_cc))
    #         vlabel_ccs.sort()
    #         vlabel_hccs.append(tuple(vlabel_ccs))
    #     vlabel_collection_value.extend(vlabel_hccs)
    #     vlabel_collection[vtx] = tuple(vlabel_collection_value)
    # return vlabel_collection
    
    # Micro-optimizations:
    # - bind local variables to avoid global lookups
    # - use generator expressions + built-in sorted to reduce Python-level loops
    # - produce tuples directly to match original return structure
    vlabel_collection = {}
    vlabel_np_loc = vlabel_np
    for vtx, hierarchical_cycle_contexts in vtx_hierarchical_cycle_contexts.items():
        # keep the same nested structure for the base label
        vlabel_collection_value = [[[int(vlabel_np_loc[vtx])], ], ]
    
        vlabel_hccs = []
        for cycle_contexts in hierarchical_cycle_contexts:
            # For each cycle (graph) produce a sorted tuple of neighbor labels (excluding vtx).
            # Use explicit list comprehensions so empty cycles produce empty tuples and
            # sorting is unambiguous.
            vlabel_ccs_list = [
                tuple(sorted([int(vlabel_np_loc[n]) for n in cycle_context if n != vtx]))
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
    for vtx, vlabel_list in vlabel_collection.items():
        vlabel_tuple = tuple(vlabel_list)
        lidx = total_label_collection.index(vlabel_tuple)
        compressed_vlabel_np[vtx] = lidx + next_label
    return compressed_vlabel_np
 
def compress_and_relabel_vlabel(vlabel_collection1, vlabel_collection2):
    total_label_collection = [v for v in vlabel_collection1.values()]
    total_label_collection.extend([v for v in vlabel_collection2.values()])
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

def node_wl_test_cycle_complexes(
    vtx_hierarchical_cycle_contexts1, vtx_hierarchical_cycle_contexts2, vlabel_np1, vlabel_np2, niter
):
    """
    vtx_hierarchical_cycle_contexts: vtx_hierarchical_cycle_contexts {vtx: [[cycle_complex1, cycle_complex2, ...], ...]}
    """
    wl_np1 = np.zeros((len(vtx_hierarchical_cycle_contexts1), niter + 1), dtype = np.int32)
    wl_np2 = np.zeros((len(vtx_hierarchical_cycle_contexts2), niter + 1), dtype = np.int32)
    wl_np1[:, 0] = vlabel_np1
    wl_np2[:, 0] = vlabel_np2
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_via_cycle_complexes(
            vtx_hierarchical_cycle_contexts1, wl_np1[:, i - 1]
        )
        vlabel_collection2 = propagate_via_cycle_complexes(
            vtx_hierarchical_cycle_contexts2, wl_np2[:, i - 1]
        )
        temp_vlabel_np1, temp_vlabel_np2 = compress_and_relabel_vlabel(
            vlabel_collection1, vlabel_collection2
        )
        wl_np1[:, i] = temp_vlabel_np1
        wl_np2[:, i] = temp_vlabel_np2
    return wl_np1, wl_np2
