import networkx as nx
import numpy as np

from hierarchical_cycle_complex_wl_tools.hierarchical_cycle_complex_bfs_neighbors import build_vtx_hierarchical_cycle_contexts

def get_vtx_cc_neighbors(g):
    node_neighs_dict = {v: g.neighbors(v) for v in g.nodes()}
    vtx_cc_neighbors = {}
    for v, neighs in node_neighs_dict.items():
        v_isg = nx.induced_subgraph(g, neighs)
        cc_neighs = []
        for cc in nx.connected_components(v_isg):
            cc_neighs.append(list(cc))
        vtx_cc_neighbors[v] = cc_neighs
    return vtx_cc_neighbors

if __name__ == '__main__':
    g1=nx.Graph()
    g1.add_edges_from(
        [
            (1,2), 
            (2,3), (2,5), (2,6), (2,8), 
            (3,5), (3,6), (3,9), 
            (4,5), 
            (5,8), (5,9), 
            (6,9), 
            (7,8), 
            (8,9)
        ]
    )
    g2=nx.Graph()
    g2.add_edges_from(
        [
            (1,2), 
            (2,3), (2,4), (2,7), 
            (3,4), (3,5), 
            (4,7), (4,8), 
            (5,8), 
            (6,7), 
            (7,8), 
        ]
    )
    dataset_info={'nl':True, 'el':False}
    graph_list=[g1,g2] 
    vlabel_list=[[1,5,4,1,5,3,1,4,4], [1,4,3,4,2,1,4,3]]
    edges_list=[list(g1.edges()), list(g2.edges())]
    elabel_list=[[0 for _ in g1.edges()], [0 for _ in g2.edges()]]
    vtx_cb_neighbors_list = []
    vtx_cc_neighbors_list = []
    deg_distr_list = []
    for g in graph_list:
        vtx_cb_neighbors = build_vtx_hierarchical_cycle_contexts(g)
        vtx_cc_neighbors = get_vtx_cc_neighbors(g)
        vtx_cb_neighbors_list.append(vtx_cb_neighbors)
        vtx_cc_neighbors_list.append(vtx_cc_neighbors)
        deg_distr = np.array([g.degree(v) for v in g.nodes()])
        deg_distr = deg_distr / deg_distr.sum().astype(np.float32)
        deg_distr_list.append(deg_distr)
    