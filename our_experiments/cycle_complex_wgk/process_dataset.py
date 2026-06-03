from multiprocessing import Pool, cpu_count
import numpy as np
import os
from grakel.datasets import fetch_dataset, get_dataset_info
import networkx as nx
import sys
cur_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(cur_dir)
parent_dir = os.path.dirname(parent_dir)
sys.path.append(parent_dir)

from utils import get_start_end_indices
from hierarchical_tri_wl_tools.hierarchical_cycle_complex_bfs_neighbors import build_vtx_hierarchical_cycle_contexts

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

def dataset_worker(g_data_list, start_gidx, end_gidx):
    graph_list = []
    vlabel_list = []
    edges_list = []
    elabel_list = []
    vtx_hierarchical_cycle_contexts_list = []
    vtx_triangulated_neighbors_list = []
    deg_distr_list = []
    # 提取图的结构信息和 Tag 信息
    for g_info in g_data_list[start_gidx: end_gidx]:
        g = nx.Graph() # 注意不连通，甚至有孤立点的图
        g.add_nodes_from(g_info[1].keys())
        g.add_edges_from(g_info[0])
        # 重新设置节点索引
        old_new_label = {
            old_label: new_label for new_label, old_label \
                in enumerate(sorted(g.nodes()))
        }
        g = nx.relabel_nodes(g, old_new_label)
        graph_list.append(g)
        vlabel_np = np.array([g.degree(v) for v in g.nodes()])
        for v, l in g_info[1].items():
            vlabel_np[old_new_label[v]] = l
        vlabel_list.append(vlabel_np)
        edges = []
        elabel = []
        if g_info[2] == {}:
            for (vi, vj) in g.edges():
                if vi <= vj:
                    e = (vi, vj)
                else:
                    e = (vj, vi)
                edges.append(e)
        else:
            for (vi, vj), l in g_info[2].items():
                if vi <= vj:
                    e = (old_new_label[vi], old_new_label[vj])
                    elabel.append(l)
                else:
                    e = (old_new_label[vj], old_new_label[vi])
                    elabel.append(l)
                edges.append(e)
        edges_list.append(edges)
        elabel_list.append(np.array(elabel))
    for g in graph_list:
        vtx_hierarchical_cycle_contexts = build_vtx_hierarchical_cycle_contexts(g)
        vtx_hierarchical_cycle_contexts_list.append(vtx_hierarchical_cycle_contexts)
        vtx_triangulated_neighbors = build_vtx_triangulated_neighbors(g)
        vtx_triangulated_neighbors_list.append(vtx_triangulated_neighbors)
        deg_distr = np.array([g.degree(v) for v in g.nodes()])
        deg_distr = deg_distr / deg_distr.sum().astype(np.float32)
        deg_distr_list.append(deg_distr)
    return start_gidx, end_gidx, graph_list, vlabel_list, edges_list, \
        elabel_list, vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list, deg_distr_list

def download_dataset(dname):
    dataset = fetch_dataset(dname, verbose=False)
    dataset_info = get_dataset_info(dname) # {'nl': bool, 'el': bool, 'na': bool, 'ea': bool, 'link': str}
    dataset_info['nl'] = True
    print((dataset.keys()))
    g_data_list, y = dataset['data'], dataset['target']
    return {'g_data': g_data_list, 'target': y, 'dataset_info': dataset_info, 'dname': dname}

def construct_dataset(dataset_dict): 
    g_data_list, y, dataset_info = dataset_dict['g_data'], dataset_dict['target'], dataset_dict['dataset_info']
    # 从 0 开始重新设置标签 label
    label_list = list(np.unique(y))
    y = np.array([label_list.index(t) for t in y])
    num_graphs = len(g_data_list)
    graph_list = [None for _ in range(num_graphs)]
    vlabel_list = [None for _ in range(num_graphs)]
    edges_list = [None for _ in range(num_graphs)]
    elabel_list = [None for _ in range(num_graphs)]
    vtx_hierarchical_cycle_contexts_list = [None for _ in range(num_graphs)]
    vtx_triangulated_neighbors_list = [None for _ in range(num_graphs)]
    deg_distr_list = [None for _ in range(num_graphs)]
    start_end_indices = get_start_end_indices(num_graphs)
    with Pool(cpu_count()) as pool:
        future_results = []
        for start_gidx, end_gidx in start_end_indices:
            fr = pool.apply_async(
                dataset_worker, 
                (g_data_list, start_gidx, end_gidx)
            )
            future_results.append(fr)
        for fr in future_results:
            start_gidx, end_gidx, graph_seg, vlabel_seg, edges_seg, elabel_seg, \
            vtx_hierarchical_cycle_contexts_seg, vtx_triangulated_neighbors_seg, deg_distr_seg = fr.get()
            graph_list[start_gidx: end_gidx] = graph_seg
            vlabel_list[start_gidx: end_gidx] = vlabel_seg
            edges_list[start_gidx: end_gidx] = edges_seg
            elabel_list[start_gidx: end_gidx] = elabel_seg
            vtx_hierarchical_cycle_contexts_list[start_gidx: end_gidx] = vtx_hierarchical_cycle_contexts_seg
            vtx_triangulated_neighbors_list[start_gidx: end_gidx] = vtx_triangulated_neighbors_seg
            deg_distr_list[start_gidx: end_gidx] = deg_distr_seg
    
    # Convert vtx_hierarchical_cycle_contexts into lightweight representations
    # Node-WL: per cycle -> tuple(sorted(node_ids_excl_center))
    # Edge-WL: per cycle -> tuple(sorted(edge_idx_list)), where edge_idx is index in edges_list
    def make_light_node_hcc(vtx_hcc):
        out = {}
        for v, layers in vtx_hcc.items():
            layers_out = []
            for layer in layers:
                cycles_out = []
                for cycle_g in layer:
                    nodes = [int(n) for n in cycle_g.nodes() if n != v]
                    cycles_out.append(tuple(nodes))
                layers_out.append(tuple(cycles_out))
            out[int(v)] = tuple(layers_out)
        return out

    def make_light_edge_hcc(vtx_hcc, edge_to_idx):
        out = {}
        for v, layers in vtx_hcc.items():
            layers_out = []
            for layer in layers:
                cycles_out = []
                for cycle_g in layer:
                    eidxs = [edge_to_idx[(a, b) if a <= b else (b, a)] for a, b in cycle_g.edges()]
                    cycles_out.append(tuple(eidxs))
                layers_out.append(tuple(cycles_out))
            out[int(v)] = tuple(layers_out)
        return out

    # Build lightweight lists per graph, replacing networkx cycle objects with tuples
    light_vtx_hcc_list = []
    for gidx, vtx_hcc in enumerate(vtx_hierarchical_cycle_contexts_list):
        # edges_list[gidx] contains normalized edge tuples
        if dataset_info.get('el', False):
            # build edge_to_idx mapping
            edge_to_idx = {e: idx for idx, e in enumerate(edges_list[gidx])}
            light_hcc = make_light_edge_hcc(vtx_hcc, edge_to_idx)
        else:
            light_hcc = make_light_node_hcc(vtx_hcc)
        light_vtx_hcc_list.append(light_hcc)

    # replace original list with lightweight version (compatible with downstream WL functions after adaptation)
    vtx_hierarchical_cycle_contexts_list = light_vtx_hcc_list

    return graph_list, vlabel_list, edges_list, elabel_list, \
        vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list, deg_distr_list, y
