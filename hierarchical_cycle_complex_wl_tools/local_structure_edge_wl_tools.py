import numpy as np
from collections import Counter

import os
import sys
cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)
from local_structure_wl_tools import compress_and_relabel_vlabel, local_structure_wl_test

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
    
def propagate_elabels(g, vlabel_np, elabel_dict):
    node_neighs_dict = {v: g.neighbors(v) for v in g.nodes()}
    vlabel_collection = {}
    for v in node_neighs_dict:
        cur_vlabels = [(vlabel_np[v], ), ]
        v_neighs = node_neighs_dict[v]
        v_induced_sg = nx.induced_subgraph(g, v_neighs)
        e_labeled_cc = list(
            map(
                lambda x: tuple(sorted(x)), 
                [[elabel_dict[(v, v_nei) if v <= v_nei else (v_nei, v)] for v_nei in v_cc] for v_cc in nx.connected_components(v_induced_sg)]
            )
        )
        e_labeled_cc.sort()
        cur_vlabels.extend(e_labeled_cc)
        vlabel_collection[v] = tuple(cur_vlabels)
    return vlabel_collection

def local_structure_edge_WL_test(
    g1, g2, 
    vlabel_np1, vlabel_np2, 
    edges1, edges2, 
    elabel_np1, elabel_np2,  
    niter
):
    vwl_np1 = np.zeros((g1.number_of_nodes(), niter + 1), dtype = np.int32)
    vwl_np2 = np.zeros((g2.number_of_nodes(), niter + 1), dtype = np.int32)
    vwl_np1[:, 0] = vlabel_np1 # copy
    vwl_np2[:, 0] = vlabel_np2 # copy
    ewl_np1 = np.zeros((len(edges1), niter + 1), dtype = np.int32)
    ewl_np2 = np.zeros((len(edges2), niter + 1), dtype = np.int32)
    ewl_np1[:, 0], ewl_np2[:, 0] = elabel_np1, elabel_np2
    compressed_elabel1, compressed_elabel2 = dict(zip(edges1, elabel_np1)), dict(zip(edges2, elabel_np2))
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_elabels(
            g1, vwl_np1[:, i - 1], compressed_elabel1)
        vlabel_collection2 = propagate_elabels(
            g2, vwl_np2[:, i - 1], compressed_elabel2)
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

def gdv_local_structure_edge_WL_test(
    dataset_info, niter, 
    g1, g2, 
    gdv1, gdv2, 
    vlabel_np1 = np.array([]), vlabel_np2 = np.array([]), 
    edges1 = [], edges2 = [], 
    elabel1 = np.array([]), elabel2 = np.array([])
):
    gdv_dim = gdv1.shape[1]
    wl_np1 = np.zeros((gdv1.shape[0], niter + 1, gdv_dim + 1), dtype = np.int32)
    wl_np2 = np.zeros((gdv2.shape[0], niter + 1, gdv_dim + 1), dtype = np.int32)
    if dataset_info['nl'] and dataset_info['el']:
        assert vlabel_np1.size != 0 and vlabel_np2.size != 0
        assert len(edges1) != 0 and len(edges2) != 0 and elabel1.size != 0 and elabel2.size != 0
        # edge labeled
        wl_np1[:, :, 0], wl_np2[:, :, 0], _, _ = local_structure_edge_WL_test(
            g1, g2, 
            vlabel_np1, vlabel_np2, 
            edges1, edges2, 
            elabel1, elabel2, 
            niter
        )
        for gdv_idx in range(gdv_dim):
            wl_np1[:, :, gdv_idx + 1], wl_np2[:, :, gdv_idx + 1] = local_structure_wl_test(
                g1, g2, 
                gdv1[:, gdv_idx], gdv2[:, gdv_idx], 
                niter
            )
    elif ((not dataset_info['nl']) and dataset_info['el']):
        assert len(edges1) != 0 and len(edges2) != 0 and elabel1.size != 0 and elabel2.size != 0
        # edge labeled
        wl_np1[:, :, 0], wl_np2[:, :, 0], _, _ = local_structure_edge_WL_test(
            g1, g2, 
            gdv1[:, 0], gdv2[:, 0], 
            edges1, edges2, 
            elabel1, elabel2, 
            niter
        )
        for gdv_idx in range(gdv_dim):
            wl_np1[:, :, gdv_idx + 1], wl_np2[:, :, gdv_idx + 1] = local_structure_wl_test(
                g1, g2, 
                gdv1[:, gdv_idx], gdv2[:, gdv_idx], 
                niter
            )
    elif (dataset_info['nl'] and (not dataset_info['el'])):
        # node labeled and edge not labeled
        wl_np1[:, :, 0], wl_np2[:, :, 0] = local_structure_wl_test(
            g1, g2, 
            vlabel_np1, vlabel_np2, 
            niter
        )
        for gdv_idx in range(gdv_dim):
            wl_np1[:, :, gdv_idx + 1], wl_np2[:, :, gdv_idx + 1] = local_structure_wl_test(
                g1, g2, 
                gdv1[:, gdv_idx], gdv2[:, gdv_idx], 
                niter
            )
    else:
        # node not labeled and edge not labeled
        for gdv_idx in range(gdv_dim):
            wl_np1[:, :, gdv_idx + 1], wl_np2[:, :, gdv_idx + 1] = local_structure_wl_test(
                g1, g2, 
                gdv1[:, gdv_idx], gdv2[:, gdv_idx], 
                niter
            )
        wl_np1, wl_np2 = np.delete(wl_np1, 0, 2), np.delete(wl_np2, 0, 2) # [:, :, o] 都是0，删除
    return wl_np1, wl_np2

# def propagate_gdv_label(g, gdv):
#     gdv_label_collection = []
#     vtx_neighs = {v: g.neighbors(v) for v in g.nodes()}
#     for gdv_idx in range(gdv.shape[1]):
#         label_collection = dict()
#         for v, v_neighs in vtx_neighs.items():
#             label_collection[(gdv[v, gdv_idx], v)] = [gdv[v_neigh, gdv_idx] for v_neigh in v_neighs]
#         gdv_label_collection.append(label_collection)
#     return gdv_label_collection

# def gdv_1st_step_mp(g1, g2, gdv1, gdv2, gamma):
#     gdv_new1, gdv_new2 = np.zeros_like(gdv1), np.zeros_like(gdv2)
#     gdv_label_collection1 = propagate_gdv_label(g1, gdv1)
#     gdv_label_collection2 = propagate_gdv_label(g2, gdv2)
#     for gdv_idx in range(gdv1.shape[1]):
#         label_collection1 = gdv_label_collection1[gdv_idx]
#         label_collection2 = gdv_label_collection2[gdv_idx]
#         gdv_vtx_flag_list1 = [(gdv_vtx, vtx, 1) for (gdv_vtx, vtx) in label_collection1.keys()]
#         gdv_vtx_flag_list2 = [(gdv_vtx, vtx, 2) for (gdv_vtx, vtx) in label_collection2.keys()]
#         sorted_gdv_vtx_flag_list = gdv_vtx_flag_list1 + gdv_vtx_flag_list2
#         sorted_gdv_vtx_flag_list.sort()
#         cur_gdv_vtx = sorted_gdv_vtx_flag_list[0][0]
#         cur_gdv_vtx_flags = []
#         cur_labels = []
#         for gdv_vtx, vtx, flag in sorted_gdv_vtx_flag_list:
#             if gdv_vtx == cur_gdv_vtx:
#                 cur_gdv_vtx_flags.append((gdv_vtx, vtx, flag))
#                 if flag == 1:
#                     cur_labels.extend(label_collection1[(gdv_vtx, vtx)])
#                 else:
#                     cur_labels.extend(label_collection2[(gdv_vtx, vtx)])
#             else:
#                 weights = Counter(cur_labels)
#                 total_w = sum(weights.values())
#                 weights = {k: w / total_w for k, w in weights.items()}
#                 for cgv, cv, cf in cur_gdv_vtx_flags:
#                     if cf == 1:
#                         lc_weights1 = {}.fromkeys(label_collection1[(cgv, cv)])
#                         lc_weights1 = {k:weights[k] for k in lc_weights1.keys()}
#                         total_lcw1 = sum(lc_weights1.values())
#                         gdv_new1[cv, gdv_idx] = sum(
#                             [k * (v / total_lcw1) for k, v in lc_weights1.items()]
#                         )
#                     else:
#                         lc_weights2 = {}.fromkeys(label_collection2[(cgv, cv)])
#                         lc_weights2 = {k:weights[k] for k in lc_weights2.keys()}
#                         total_lcw2 = sum(lc_weights2.values())
#                         gdv_new2[cv, gdv_idx] = sum(
#                             [k * (v / total_lcw2) for k, v in lc_weights2.items()]
#                         )
#                 cur_gdv_vtx_flags.clear()
#                 cur_labels.clear()
#                 cur_gdv_vtx = gdv_vtx
#                 cur_gdv_vtx_flags.append((gdv_vtx, vtx, flag))
#                 if flag == 1:
#                     cur_labels.extend(label_collection1[(gdv_vtx, vtx)])
#                 else:
#                     cur_labels.extend(label_collection2[(gdv_vtx, vtx)])
#     return np.round(gamma * gdv1 + (1 - gamma) * gdv_new1), np.round(gamma * gdv2 + (1 - gamma) * gdv_new2)

# def gdv_mp(g1, g2, gdv1, gdv2, gamma, n):
#     gdv_list1, gdv_list2 = [gdv1, ], [gdv2, ]
#     gdv_new1, gdv_new2 = gdv1, gdv2
#     for _ in range(n):
#         gdv_new1, gdv_new2 = gdv_1st_step_mp(g1, g2, gdv_new1, gdv_new2, gamma)
#         gdv_list1.append(gdv_new1)
#         gdv_list2.append(gdv_new2)
#     return np.concatenate(gdv_list1, axis = 1), np.concatenate(gdv_list2, axis = 1)
                
if __name__ == '__main__':
    import networkx as nx
    g1 = nx.Graph()
    g1.add_edges_from(
        [
            (0, 1), (0, 2), (0, 6), (0, 7), 
            (1, 2), (1,3), (1, 7), 
            (2, 3), (2, 4), 
            (3, 4), (3, 5), 
            (4, 5), (4, 6), 
            (5, 6), (5, 7), 
            (6, 7)
        ]
    )
    vlabel_np1 = np.ones((8, ), dtype = np.int32)
    v0_neighs = set(g1.neighbors(0))
    v1_neighs = set(g1.neighbors(1))
    elabel_tuple10, elabel_tuple11 = {}, {}
    for vi, vj in g1.edges():
        l0 = len(list(nx.common_neighbors(g1, vi, vj)))
        l1 = len(set(g1.neighbors(vi)).union(set(g1.neighbors(vj)))) - l0
        elabel_tuple10[(vi, vj)] = l0
        elabel_tuple11[(vi, vj)] = l1
    g2 = nx.Graph()
    g2.add_edges_from(
        [
            (0, 1), (0, 3), (0, 5), (0, 7), 
            (1, 2), (1, 4), (1, 6), 
            (2, 3), (2, 5), (2, 7), 
            (3, 4), (3, 6), 
            (4, 5), (4, 7), 
            (5, 6), 
            (6, 7)
        ]
    )
    vlabel_np2 = np.ones((8, ), dtype = np.int32)
    elabel_tuple20, elabel_tuple21 = {}, {}
    for vi, vj in g2.edges():
        l0 = len(list(nx.common_neighbors(g2, vi, vj)))
        l1 = len(set(g2.neighbors(vi)).union(set(g2.neighbors(vj)))) - l0
        elabel_tuple20[(vi, vj)] = l0
        elabel_tuple21[(vi, vj)] = l1
    wl_np1, wl_np2 = local_structure_edge_WL_test(
        g1, 
        g2, 
        vlabel_np1, 
        vlabel_np2, 
        elabel_tuple11, 
        elabel_tuple21, 
        2
    )
    print(wl_np1)
    print(wl_np2)
