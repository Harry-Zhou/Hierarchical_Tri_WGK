import networkx as nx
import numpy as np
import os

import sys
cur_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(cur_dir)
parent_dir = os.path.dirname(parent_dir)
sys.path.append(parent_dir)
from topo_wasserstein_graph_kernel import TopoWassersteinGraphKernel
from hierarchical_tri_wl_tools.hierarchical_cycle_complex_bfs_neighbors import build_vtx_hierarchical_cycle_contexts
from utils import precomp_node_neighs

def build_vtx_triangulated_neighbors(g):
    node_neighs_dict = {v: g.neighbors(v) for v in g.nodes()}
    vtx_triangulated_neighbors = {}
    for v, neighs in node_neighs_dict.items():
        v_isg = nx.induced_subgraph(g, neighs)
        tri_neighs = []
        for ttn in nx.connected_components(v_isg):
            ttn_v = list(ttn) + [v, ]
            tri_subgraph = nx.induced_subgraph(g, ttn_v)
            tri_neighs.append(tri_subgraph)
        vtx_triangulated_neighbors[v] = tri_neighs
    return vtx_triangulated_neighbors

def load_graph_from_edgelist(filepath):
    g = nx.read_edgelist(filepath, nodetype=int)
    g = nx.convert_node_labels_to_integers(g, first_label=0)
    return g

def prepare_graph_data(g):
    vlabel_np = np.array([g.degree(v) for v in g.nodes()], dtype=np.int32)
    edges = []
    for (vi, vj) in g.edges():
        if vi <= vj:
            edges.append((vi, vj))
        else:
            edges.append((vj, vi))
    edges_list = edges
    elabel_list = np.array([])
    
    node_neighs = precomp_node_neighs(g)
    vtx_hierarchical_cycle_contexts = build_vtx_hierarchical_cycle_contexts(g)
    vtx_triangulated_neighbors = build_vtx_triangulated_neighbors(g)
    deg_distr = np.array([g.degree(v) for v in g.nodes()])
    deg_distr = deg_distr / deg_distr.sum().astype(np.float32)
    
    return vlabel_np, edges_list, elabel_list, node_neighs, vtx_hierarchical_cycle_contexts, vtx_triangulated_neighbors, deg_distr

def test_graph_isomorphism(graph1_path, graph2_path):
    g1 = load_graph_from_edgelist(graph1_path)
    g2 = load_graph_from_edgelist(graph2_path)
    
    print(f"Graph 1: {g1.number_of_nodes()} nodes, {g1.number_of_edges()} edges")
    print(f"Graph 2: {g2.number_of_nodes()} nodes, {g2.number_of_edges()} edges")
    
    dataset_info = {'nl': True, 'el': False}
    
    vlabel1, edges1, elabel1, node_neighs1, vtx_hcc1, vtx_tri1, deg_distr1 = prepare_graph_data(g1)
    vlabel2, edges2, elabel2, node_neighs2, vtx_hcc2, vtx_tri2, deg_distr2 = prepare_graph_data(g2)
    
    graph_list = [g1, g2]
    vlabel_list = [vlabel1, vlabel2]
    edges_list = [edges1, edges2]
    elabel_list = [elabel1, elabel2]
    vtx_hierarchical_cycle_contexts_list = [vtx_hcc1, vtx_hcc2]
    vtx_triangulated_neighbors_list = [vtx_tri1, vtx_tri2]
    deg_distr_list = [deg_distr1, deg_distr2]
    
    topo_wgk = TopoWassersteinGraphKernel(n_wl_iters=3, n_csg_layers=3, wl_normalized=True)
    
    ot_dist_np, wl_sim_np, run_time = topo_wgk.fit_transform(
        dataset_info, graph_list, vlabel_list, edges_list, elabel_list,
        vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list,
        deg_distr_list
    )
    
    print(f"\nOptimal Transport distance matrix:")
    print(ot_dist_np)
    print(f"\nWL similarity matrix:")
    print(wl_sim_np)
    
    graph_dist = ot_dist_np[0, 1]
    
    if np.isclose(graph_dist, 0.0):
        print(f"\nGraph distance is {graph_dist:.6f} (approximately 0)")
        print("Conclusion: The two graphs are isomorphic (same structure)")
    else:
        print(f"\nGraph distance is {graph_dist:.6f} (not 0)")
        print("Conclusion: The two graphs are NOT isomorphic (different structures)")
    
    return graph_dist, ot_dist_np, wl_sim_np

if __name__ == '__main__':
    graph1_path = './our_experiments/cycle_complex_wgk/graph_isomorphism_testing/graph1.edgelist'
    graph2_path = './our_experiments/cycle_complex_wgk/graph_isomorphism_testing/graph2.edgelist'
    
    print("=" * 80)
    print("Testing Graph Isomorphism using TopoWassersteinGraphKernel")
    print("=" * 80)
    
    graph_dist, ot_dist_np, wl_sim_np = test_graph_isomorphism(graph1_path, graph2_path)
    
    print("=" * 80)
    print("Test completed")
    print("=" * 80)
