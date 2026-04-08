import numpy as np

def propagate_vlabel(node_neighs_dict, vlabel_np):
    vlabel_collection = {}
    for v in node_neighs_dict:
        v_neigh_label = [vlabel_np[nv] for nv in node_neighs_dict[v]]
        v_neigh_label.sort() # 邻居节点的标签排序
        cur_vlabels = [vlabel_np[v], ]
        cur_vlabels.extend(
            v_neigh_label
        )
        vlabel_collection[v] = tuple(cur_vlabels)
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
    next_label = max([lc[0] for lc in total_label_collection]) + 1
    node_compressed_label_np1 = gen_compressed_vlabel(
        vlabel_collection1, total_label_collection, next_label
    )
    node_compressed_label_np2 = gen_compressed_vlabel(
        vlabel_collection2, total_label_collection, next_label
    )
    return node_compressed_label_np1, node_compressed_label_np2

def classic_wl_test(
    g1, g2, 
    vlabel_np1, vlabel_np2, 
    niter
):
    node_neighs_dict1 = {v: g1.neighbors(v) for v in g1.nodes()}
    node_neighs_dict2 = {v: g2.neighbors(v) for v in g2.nodes()}
    wl_np1 = np.zeros((g1.number_of_nodes(), niter + 1), dtype = np.int32)
    wl_np2 = np.zeros((g2.number_of_nodes(), niter + 1), dtype = np.int32)
    wl_np1[:, 0] = vlabel_np1 # copy
    wl_np2[:, 0] = vlabel_np2 # copy
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_vlabel(
            node_neighs_dict1, wl_np1[:, i - 1]
        )
        vlabel_collection2 = propagate_vlabel(
            node_neighs_dict2, wl_np2[:, i-1]
        )
        temp_vlabel_np1, temp_vlabel_np2 = \
            compress_and_relabel_vlabel(
                vlabel_collection1, vlabel_collection2
            )
        wl_np1[:, i] = temp_vlabel_np1
        wl_np2[:, i] = temp_vlabel_np2
    return wl_np1, wl_np2

def gdv_WL_test(g1, g2, gdv1, gdv2, niter):
    node_neighs_dict1 = {v: g1.neighbors(v) for v in g1.nodes()}
    node_neighs_dict2 = {v: g2.neighbors(v) for v in g2.nodes()}
    gdv_wl1 = np.zeros((*gdv1.shape, niter + 1), dtype = np.int32)
    gdv_wl2 = np.zeros((*gdv2.shape, niter + 1), dtype = np.int32)
    ndim_gdv = gdv1.shape[1]
    for gdv_idx in range(ndim_gdv):
        gdv_wl1[:, gdv_idx, 0] = gdv1[:, gdv_idx] # copy
        gdv_wl2[:, gdv_idx, 0] = gdv2[:, gdv_idx] # copy
        for i in range(1, niter + 1):
            vlabel_collection1 = propagate_vlabel(
                node_neighs_dict1, gdv_wl1[:, gdv_idx, i - 1]
            )
            vlabel_collection2 = propagate_vlabel(
                node_neighs_dict2, gdv_wl2[:, gdv_idx, i - 1]
            )
            temp_node_compressed_label_np1, temp_node_compressed_label_np2 = \
                compress_and_relabel_vlabel(
                    vlabel_collection1, vlabel_collection2
                )
            gdv_wl1[:, gdv_idx, i] = temp_node_compressed_label_np1
            gdv_wl2[:, gdv_idx, i] = temp_node_compressed_label_np2
    return gdv_wl1, gdv_wl2

if __name__ == '__main__':
    import networkx as nx
    g1 = nx.Graph()
    g1.add_edges_from(
        [
            (0,1), (0,3), (0,2),
            (1,3),
            (2,3), (2,4), (2,5)
        ]
    )
    # g2 = nx.Graph()
    # g2.add_edges_from(
    #     [
    #         (0,1), (0,2),
    #         (1,2), (1,3),
    #         (2,3), (2,4), 
    #         (3,5)
    #     ]
    # )
    g2 = nx.relabel_nodes(g1, {0:1, 1:3, 2:4, 3:5, 4:0, 5:2})
    vlabel_np1 = np.array([5,2,4,3,1,1])
    vlabel_np2 = np.array([1,5,1,2,4,3])
    wl_np1, wl_np2 = classic_wl_test(
        g1, g2, 
        vlabel_np1, vlabel_np2, 
        5
    )
    print(wl_np1)
    print(wl_np2)
