import numpy as np
import networkx as nx

def propagate_connected_vlabel(g, vlabel_np):
    node_neighs_dict = {v: g.neighbors(v) for v in g.nodes()}
    vlabel_collection = {}
    i = 0
    for v in node_neighs_dict:
        cur_vlabels = [(vlabel_np[v], ), ]
        v_neighs = node_neighs_dict[v]
        v_induced_sg = nx.induced_subgraph(g, v_neighs)
        v_labeled_cc = list(
            map(
                lambda x: tuple(sorted(x)), 
                [[vlabel_np[v_nei] for v_nei in v_cc] for v_cc in nx.connected_components(v_induced_sg)]
            )
        )
        v_labeled_cc.sort()
        cur_vlabels.extend(v_labeled_cc)
        vlabel_collection[v] = tuple(cur_vlabels)
        i += 1
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

def local_structure_wl_test(
    g1, g2, 
    vlabel_np1, vlabel_np2, 
    niter
):
    wl_np1 = np.zeros((g1.number_of_nodes(), niter + 1), dtype = np.int32)
    wl_np2 = np.zeros((g2.number_of_nodes(), niter + 1), dtype = np.int32)
    wl_np1[:, 0] = vlabel_np1 # copy
    wl_np2[:, 0] = vlabel_np2 # copy
    for i in range(1, niter + 1):
        vlabel_collection1 = propagate_connected_vlabel(
            g1, wl_np1[:, i - 1]
        )
        vlabel_collection2 = propagate_connected_vlabel(
            g2, wl_np2[:, i-1]
        )
        temp_vlabel_np1, temp_vlabel_np2 = \
            compress_and_relabel_vlabel(
                vlabel_collection1, vlabel_collection2
            )
        wl_np1[:, i] = temp_vlabel_np1
        wl_np2[:, i] = temp_vlabel_np2
    return wl_np1, wl_np2

def gdv_local_structure_WL_test(g1, g2, gdv1, gdv2, niter):
    gdv_wl1 = np.zeros((*gdv1.shape, niter + 1), dtype = np.int32)
    gdv_wl2 = np.zeros((*gdv2.shape, niter + 1), dtype = np.int32)
    ndim_gdv = gdv1.shape[1]
    for gdv_idx in range(ndim_gdv):
        gdv_wl1[:, gdv_idx, 0] = gdv1[:, gdv_idx] # copy
        gdv_wl2[:, gdv_idx, 0] = gdv2[:, gdv_idx] # copy
        for i in range(1, niter + 1):
            vlabel_collection1 = propagate_connected_vlabel(
                g1, gdv_wl1[:, gdv_idx, i - 1]
            )
            vlabel_collection2 = propagate_connected_vlabel(
                g2, gdv_wl2[:, gdv_idx, i - 1]
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
    wl_np1, wl_np2 = local_structure_wl_test(
        g1, g2, 
        vlabel_np1, vlabel_np2, 
        5
    )
    print(wl_np1)
    print(wl_np2)
