import networkx as nx
import numpy as np
import itertools
from collections import Counter
from functools import reduce
from multiprocessing import Pool, cpu_count

import os
import sys
cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)
parent_dir = os.path.dirname(cur_dir)
sys.path.append(parent_dir)

def extract_vtx_maximal_cycle_complexes(g, achordal_cycle_basis_list):
    """
    构造一个 Graph, 节点是 achordal cycle basis, 
    边是两个 achordal cycle basis 共享至少一条边, 
    通过寻找该图的连通分量来提取 maximal cycle complexes
    """
    achordal_cycle_graph = nx.Graph()
    achordal_cycle_graph.add_nodes_from(range(len(achordal_cycle_basis_list)))
    for i, achordal_cycle_basis in enumerate(achordal_cycle_basis_list):
        achordal_cycle_graph.nodes[i]['achordal_cycle_basis'] = achordal_cycle_basis
    for (i,j) in itertools.combinations(achordal_cycle_graph, 2):
        cbi, cbj = achordal_cycle_basis_list[i], achordal_cycle_basis_list[j]
        cbi = [tuple(sorted(e)) for e in g.subgraph(cbi).edges()]
        cbj = [tuple(sorted(e)) for e in g.subgraph(cbj).edges()]
        if set(cbi).intersection(set(cbj)) != set():
            achordal_cycle_graph.add_edge(i, j)
    maximal_cycle_complexes = []
    for cc in nx.connected_components(achordal_cycle_graph):
        maximal_cycle_complexes.append(achordal_cycle_graph.subgraph(cc))
    
    # vtx 与 maximal_cycle_complexes 的映射（子图）
    vtx_maximal_cycle_complexes = {vtx: [] for vtx in g.nodes()}
    for maximal_cycle_complex in maximal_cycle_complexes:
        for cbidx in maximal_cycle_complex.nodes():
            for vtx in achordal_cycle_graph.nodes[cbidx]['achordal_cycle_basis']:
                vtx_maximal_cycle_complexes[vtx].append(maximal_cycle_complex)
    
    return vtx_maximal_cycle_complexes

def extract_vtx_achordal_cycle_basis_idxlist(achordal_cycle_basis_list):
    vtx_achordal_cycle_basis_idxlist = {}
    for vtx_skeleton_cb, achordal_cycle_basis in enumerate(achordal_cycle_basis_list):
        for vtx in achordal_cycle_basis:
            if vtx not in vtx_achordal_cycle_basis_idxlist:
                vtx_achordal_cycle_basis_idxlist[vtx] = []
            vtx_achordal_cycle_basis_idxlist[vtx].append(vtx_skeleton_cb)
    return vtx_achordal_cycle_basis_idxlist

def init_vtx_hierarchical_cycle_contexts(g, achordal_cycle_basis_list):
    vtx_hierarchical_cycle_contexts = {vtx: {} for vtx in g.nodes()}
    vtx_achordal_cycle_basis_idxlist = extract_vtx_achordal_cycle_basis_idxlist(achordal_cycle_basis_list)
    for vtx, contextual_achordal_cycle_basis_idxlist in vtx_achordal_cycle_basis_idxlist.items():
        # 第 0 层 context 
        vtx_hierarchical_cycle_contexts[vtx][0] = [
            achordal_cycle_basis_list[vtx_skeleton_cb] for vtx_skeleton_cb in contextual_achordal_cycle_basis_idxlist
        ]
    for vtx in g.nodes():
        if vtx_hierarchical_cycle_contexts[vtx] == {}:
            vtx_hierarchical_cycle_contexts[vtx][0] = [g.subgraph([vtx, v_nei]) for v_nei in g.neighbors(vtx)]
        else:
            h0_cycle_contexts = vtx_hierarchical_cycle_contexts[vtx][0]
            acycle_neighbors = set(g.neighbors(vtx)).difference(reduce(lambda x,y: x.union(y), h0_cycle_contexts, set()))
            vtx_hierarchical_cycle_contexts[vtx][0] = [g.subgraph([x, vtx]) for x in acycle_neighbors] + h0_cycle_contexts
    
    return vtx_hierarchical_cycle_contexts, vtx_achordal_cycle_basis_idxlist

def build_vtx_hierarchical_cycle_contexts(g):
    achordal_cycle_basis_list = [g.subgraph(mcb) for mcb in nx.minimum_cycle_basis(g)]
    vtx_hierarchical_cycle_contexts, vtx_achordal_cycle_basis_idxlist = init_vtx_hierarchical_cycle_contexts(g, achordal_cycle_basis_list)
    vtx_maximal_cycle_complexes = extract_vtx_maximal_cycle_complexes(g, achordal_cycle_basis_list)
    for vtx in g.nodes():
        if vtx in vtx_achordal_cycle_basis_idxlist:
            achordal_cycle_basis_idxlist = vtx_achordal_cycle_basis_idxlist[vtx]
            for maximal_cycle_complex in vtx_maximal_cycle_complexes[vtx]:
                layer_achordal_cb_idxlist = nx.bfs_layers(
                    maximal_cycle_complex, 
                    [cbidx for cbidx in achordal_cycle_basis_idxlist if cbidx in maximal_cycle_complex]
                )
                for layer, achordal_cb_idxlist in enumerate(layer_achordal_cb_idxlist):
                    if layer > 0:
                        if layer not in vtx_hierarchical_cycle_contexts[vtx]:
                            vtx_hierarchical_cycle_contexts[vtx][layer] = [
                                maximal_cycle_complex.nodes[achordal_cb]['achordal_cycle_basis']
                                for achordal_cb in achordal_cb_idxlist
                            ]
                        else:
                            vtx_hierarchical_cycle_contexts[vtx][layer] += [
                                maximal_cycle_complex.nodes[achordal_cb]['achordal_cycle_basis']
                                for achordal_cb in achordal_cb_idxlist
                            ]
    for vtx, hierarchical_cycle_contexts in vtx_hierarchical_cycle_contexts.items():
        hierarchical_cycle_contexts_list = [None for _ in range(len(hierarchical_cycle_contexts))]
        for layer, cycle_contexts in hierarchical_cycle_contexts.items():
            hierarchical_cycle_contexts_list[layer] = cycle_contexts
        vtx_hierarchical_cycle_contexts[vtx] = hierarchical_cycle_contexts_list      
    return vtx_hierarchical_cycle_contexts

if __name__ == '__main__':
    import itertools
    import matplotlib.pyplot as plt

    g = nx.Graph()
    for cycle in[
        ('a', 'b', 'c', 'd', 'a'),
        ('b', 'c', 'e', 'b'),
        ('c', 'e', 'g', 'c'),
        ('c', 'd', 'g', 'c'),
        ('b', 'e', 'f', 'b'),
        ('f', 'h', 'i', 'f'),
        ('j', 'k', 'l', 'j'),
        ('m', 'n', 'o', 'm'),
    ]:
        g.add_edges_from([(u, v) for u, v in zip(cycle[0:-1], cycle[1:])])
    g.add_edges_from([('a', 'p'), ('f', 'q'), ('g', 'r'), ('g', 's'), ('s', 'm'), ('s', 't'), ('h', 'j'), ('i', 'k')])
    
    vtx_hierarchical_cycle_contexts = build_vtx_hierarchical_cycle_contexts(g)
    for vtx, hierarchical_cycle_contexts in vtx_hierarchical_cycle_contexts.items():
        print(f"Vertex {vtx}:")
        for layer, cycle_contexts in enumerate(hierarchical_cycle_contexts):
            print(f"  Layer {layer}: {[list(ccg.nodes()) for ccg in cycle_contexts]}, {[list(ccg.edges()) for ccg in cycle_contexts]}")
    
    # cbs = [cb for cb in nx.minimum_cycle_basis(g)]
    # for i, cb in enumerate(cbs):
    #     print(f"Cycle Basis {i}: {cb}")
    # cg = nx.Graph()
    # cg.add_nodes_from(range(len(cbs)))
    # for (i,j) in itertools.combinations(range(len(cbs)), 2):
    #     cbi, cbj = cbs[i], cbs[j]
    #     sgi, sgj = g.subgraph(cbi), g.subgraph(cbj)
    #     cbi = [tuple(sorted(e)) for e in sgi.edges()]
    #     cbj = [tuple(sorted(e)) for e in sgj.edges()]
    #     # print(f"Comparing Cycle Basis {i} and {j}: {cbi} vs {cbj}")
    #     if set(cbi).intersection(set(cbj)) != set():
    #         cg.add_edge(i, j)
    
    # # 设置图形布局
    # pos = nx.circular_layout(g)  # 可选择其他布局：circular_layout, random_layout...

    # # 绘制图形
    # plt.figure(figsize=(10, 8))
    # nx.draw_networkx_nodes(g, pos, node_size=100, node_color='lightblue')
    # nx.draw_networkx_edges(g, pos, width=0.5, alpha=0.6)
    # nx.draw_networkx_labels(g, pos, font_size=8, font_family='sans-serif')

    # # 显示图形
    # plt.axis('off')
    # plt.tight_layout()
    # plt.savefig("./cycle_basis_graph.jpg")
    