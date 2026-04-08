from collections import Counter
from itertools import product
from multiprocessing import Pool, cpu_count
import numpy as np
import ot
from scipy.spatial.distance import cdist
import os
import time

import sys
cur_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(cur_dir)
parent_dir = os.path.dirname(parent_dir)
sys.path.append(parent_dir)
from utils import get_start_end_indices, precomp_node_neighs, get_recommended_nproc
from hierarchical_cycle_complex_wl_tools.topo_aware_wl_test import topo_aware_wl_test

# Module-level globals populated in worker initializer to avoid passing large objects per-task
_GLOBAL_VTX_HCC = None
_GLOBAL_VTX_TRI = None
_GLOBAL_DEG_DISTR = None
_GLOBAL_NODE_NEIGHS = None
_GLOBAL_VLABEL_LIST = None
_GLOBAL_EDGES_LIST = None
_GLOBAL_ELABEL_LIST = None
_GLOBAL_DATASET_INFO = None
_GLOBAL_NITER_TN = None
_GLOBAL_NITER_HCC = None

def _worker_init(vtx_hcc_list, vtx_tri_list, deg_distr_list, node_neighs_list,
                 vlabel_list, edges_list, elabel_list, dataset_info, niter_tn, niter_hcc):
    global _GLOBAL_VTX_HCC, _GLOBAL_VTX_TRI, _GLOBAL_DEG_DISTR, _GLOBAL_NODE_NEIGHS
    global _GLOBAL_VLABEL_LIST, _GLOBAL_EDGES_LIST, _GLOBAL_ELABEL_LIST
    global _GLOBAL_DATASET_INFO, _GLOBAL_NITER_TN, _GLOBAL_NITER_HCC
    _GLOBAL_VTX_HCC = vtx_hcc_list
    _GLOBAL_VTX_TRI = vtx_tri_list
    _GLOBAL_DEG_DISTR = deg_distr_list
    _GLOBAL_NODE_NEIGHS = node_neighs_list
    _GLOBAL_VLABEL_LIST = vlabel_list
    _GLOBAL_EDGES_LIST = edges_list
    _GLOBAL_ELABEL_LIST = elabel_list
    _GLOBAL_DATASET_INFO = dataset_info
    _GLOBAL_NITER_TN = niter_tn
    _GLOBAL_NITER_HCC = niter_hcc

def _compute_block_worker(block):
    start_ridx, end_ridx, start_cidx, end_cidx = block
    ot_dist_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
    wl_sim_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
    temp_start_cidx = start_cidx
    is_diagonal = start_ridx == start_cidx and end_ridx == end_cidx
    print(f'\t The ({start_ridx}, {end_ridx}, {start_cidx}, {end_cidx}) graph similarity submatrix')
    for gridx in range(start_ridx, end_ridx):
        dridx = gridx - start_ridx
        r_vtx_hierarchical_cycle_contexts = _GLOBAL_VTX_HCC[gridx]
        r_vtx_triangulated_neighbors = _GLOBAL_VTX_TRI[gridx]
        r_deg_distr = _GLOBAL_DEG_DISTR[gridx]
        r_node_neighs = _GLOBAL_NODE_NEIGHS[gridx]
        rvlabel_np = _GLOBAL_VLABEL_LIST[gridx]
        redges = _GLOBAL_EDGES_LIST[gridx]
        relabel_np = _GLOBAL_ELABEL_LIST[gridx]
        if is_diagonal:
            temp_start_cidx = gridx + 1
        for gcidx in range(temp_start_cidx, end_cidx):
            c_vtx_hierarchical_cycle_contexts = _GLOBAL_VTX_HCC[gcidx]
            c_vtx_triangulated_neighbors = _GLOBAL_VTX_TRI[gcidx]
            c_deg_distr = _GLOBAL_DEG_DISTR[gcidx]
            c_node_neighs = _GLOBAL_NODE_NEIGHS[gcidx]
            cvlabel_np = _GLOBAL_VLABEL_LIST[gcidx]
            cedges = _GLOBAL_EDGES_LIST[gcidx]
            celabel_np = _GLOBAL_ELABEL_LIST[gcidx]
            rvwl_np, cvwl_np = topo_aware_wl_test(
                _GLOBAL_DATASET_INFO, _GLOBAL_NITER_TN, _GLOBAL_NITER_HCC,
                r_vtx_triangulated_neighbors, c_vtx_triangulated_neighbors,
                r_vtx_hierarchical_cycle_contexts, c_vtx_hierarchical_cycle_contexts,
                rvlabel_np, cvlabel_np,
                redges, cedges, relabel_np, celabel_np
            )
            # compute transport cost using the same helper logic as the class method
            from collections import Counter as _Counter
            def _wl_counter(node_neighs, vwl_np, vmax_label):
                vwl_counter = np.zeros((vwl_np.shape[0], vmax_label + 1))
                for ridx, neighs in node_neighs.items():
                    for l, n in _Counter(vwl_np[neighs, :].flatten('F')).items():
                        vwl_counter[ridx, l] = n
                return vwl_counter
            native_vmax_label = max(rvwl_np.max(), cvwl_np.max())
            native_vwl_counter1 = _wl_counter(_GLOBAL_NODE_NEIGHS[gridx], rvwl_np, native_vmax_label)
            native_vwl_counter2 = _wl_counter(_GLOBAL_NODE_NEIGHS[gcidx], cvwl_np, native_vmax_label)
            transport_cost = cdist(native_vwl_counter1, native_vwl_counter2, metric = 'sqeuclidean')

            dcidx = gcidx - start_cidx
            ot_dist_diag[dridx, dcidx] = ot.emd2(
                _GLOBAL_DEG_DISTR[gridx], _GLOBAL_DEG_DISTR[gcidx], transport_cost
            )
            wl_sim_diag[dridx, dcidx] = np.dot(
                (rvwl_np.flatten('F')==rvwl_np).sum(), (cvwl_np.flatten('F')==cvwl_np).sum()
            ) if False else TopoWassersteinGraphKernel.wl_inner_product(None, rvwl_np, cvwl_np)

        if is_diagonal:
            rvwl_np, cvwl_np = topo_aware_wl_test(
                _GLOBAL_DATASET_INFO, _GLOBAL_NITER_TN, _GLOBAL_NITER_HCC,
                r_vtx_triangulated_neighbors, r_vtx_triangulated_neighbors,
                r_vtx_hierarchical_cycle_contexts, r_vtx_hierarchical_cycle_contexts,
                rvlabel_np, rvlabel_np,
                redges, redges, relabel_np, relabel_np
            )
            wl_sim_diag[dridx, dridx] = TopoWassersteinGraphKernel.wl_inner_product(None, rvwl_np, cvwl_np)
    return (start_ridx, end_ridx, start_cidx, end_cidx, ot_dist_diag, wl_sim_diag)

class TopoWassersteinGraphKernel:
    def __init__(self, niter_tn, niter_hcc, wl_normalized = True):
        self._niter_tn = niter_tn
        self._niter_hcc = niter_hcc
        self._wl_normalized = wl_normalized
    
    def wl_inner_product(self, vwl_np1, vwl_np2):
        vmax_label = max(vwl_np1.max(), vwl_np2.max())
        label_distr1 = np.zeros((vmax_label + 1, ))
        label_distr2 = np.zeros_like(label_distr1)
        for l, n in Counter(vwl_np1.flatten('F')).items():
            label_distr1[l] = n
        for l, n in Counter(vwl_np2.flatten('F')).items():
            label_distr2[l ] = n
        return np.dot(label_distr1, label_distr2)
    
    def _comp_transport_cost(
        self, node_neighs1, node_neighs2, vwl_np1, vwl_np2
    ):
        """
        vwl_np1, vwl_np2: (num_nodes, niter + 1, gdv_dim + 1)
        """
        def _wl_counter(node_neighs, vwl_np, vmax_label):
            vwl_counter = np.zeros((vwl_np.shape[0], vmax_label + 1))
            for ridx, neighs in node_neighs.items():
                for l, n in Counter(vwl_np[neighs, :].flatten('F')).items():
                    vwl_counter[ridx, l] = n
            return vwl_counter
        
        native_vmax_label = max(vwl_np1.max(), vwl_np2.max())
        native_vwl_counter1 = _wl_counter(
            node_neighs1, vwl_np1, native_vmax_label
        )
        native_vwl_counter2 = _wl_counter(
            node_neighs2, vwl_np2, native_vmax_label
        )
        transport_cost = cdist(
            native_vwl_counter1, native_vwl_counter2, metric = 'sqeuclidean'
        )
        
        return transport_cost
    
    def _otdist_wlsim_worker(self, start_ridx, end_ridx, start_cidx, end_cidx):
        ot_dist_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
        wl_sim_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
        temp_start_cidx = start_cidx
        is_diagonal = start_ridx == start_cidx and end_ridx == end_cidx
        print(f'\t The ({start_ridx}, {end_ridx}, {start_cidx}, {end_cidx}) graph similarity submatrix')
        for gridx in range(start_ridx, end_ridx):
            dridx = gridx - start_ridx
            # gridx = indices_data[i]
            r_vtx_hierarchical_cycle_contexts = self.vtx_hierarchical_cycle_contexts_list[gridx]
            r_vtx_triangulated_neighbors = self.vtx_triangulated_neighbors_list[gridx]
            r_deg_distr = self.deg_distr_list[gridx]
            # r_graph = self.graph_list[gridx]
            r_node_neighs = self.node_neighs_list[gridx]
            rvlabel_np = self.vlabel_list[gridx]
            redges = self.edges_list[gridx]
            relabel_np = self.elabel_list[gridx]
            if is_diagonal:
                temp_start_cidx = gridx + 1
            for gcidx in range(temp_start_cidx, end_cidx):
                # gcidx = indices_data[j]
                c_vtx_hierarchical_cycle_contexts = self.vtx_hierarchical_cycle_contexts_list[gcidx]
                c_vtx_triangulated_neighbors = self.vtx_triangulated_neighbors_list[gcidx]
                c_deg_distr = self.deg_distr_list[gcidx]
                # c_graph = self.graph_list[gcidx]
                c_node_neighs = self.node_neighs_list[gcidx]
                cvlabel_np = self.vlabel_list[gcidx]
                cedges = self.edges_list[gcidx]
                celabel_np = self.elabel_list[gcidx]
                rvwl_np, cvwl_np = topo_aware_wl_test(
                    self._dataset_info, self._niter_tn, self._niter_hcc, 
                    r_vtx_triangulated_neighbors, c_vtx_triangulated_neighbors, 
                    r_vtx_hierarchical_cycle_contexts, c_vtx_hierarchical_cycle_contexts, 
                    rvlabel_np, cvlabel_np, 
                    redges, cedges, relabel_np, celabel_np
                )
                transport_cost = self._comp_transport_cost(
                    r_node_neighs, c_node_neighs, rvwl_np, cvwl_np
                )
                dcidx = gcidx - start_cidx
                ot_dist_diag[dridx, dcidx] = ot.emd2(
                    r_deg_distr, c_deg_distr, transport_cost
                )
                wl_sim_diag[dridx, dcidx] = self.wl_inner_product(rvwl_np, cvwl_np)
            
            if is_diagonal:
                rvwl_np, cvwl_np = topo_aware_wl_test(
                    self._dataset_info, self._niter_tn, self._niter_hcc, 
                    r_vtx_triangulated_neighbors, r_vtx_triangulated_neighbors, 
                    r_vtx_hierarchical_cycle_contexts, r_vtx_hierarchical_cycle_contexts, 
                    rvlabel_np, rvlabel_np, 
                    redges, redges, relabel_np, relabel_np
                ) # (num_nodes, niter)
                wl_sim_diag[dridx, dridx] = self.wl_inner_product(rvwl_np, cvwl_np)
        return (start_ridx, end_ridx, start_cidx, end_cidx, ot_dist_diag, wl_sim_diag)

    def fit(
        self, dataset_info, graph_list, vlabel_list, edges_list, elabel_list, 
        vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list, 
        deg_distr_list
    ):
        self._dataset_info = dataset_info
        self._nl, self._el = dataset_info['nl'], dataset_info['el']
        self.vtx_hierarchical_cycle_contexts_list = vtx_hierarchical_cycle_contexts_list
        self.vtx_triangulated_neighbors_list = vtx_triangulated_neighbors_list
        self.deg_distr_list = deg_distr_list
        self.node_neighs_list = [precomp_node_neighs(g) for g in graph_list]
        self.graph_list = graph_list
        self.num_graphs = len(graph_list)
        self.vlabel_list = vlabel_list
        self.edges_list = edges_list
        self.elabel_list = elabel_list
    
    def transform(self):
        ot_dist_diag = np.zeros((self.num_graphs, self.num_graphs))
        wl_sim_diag = np.zeros((self.num_graphs, self.num_graphs))
        start_end_indices = get_start_end_indices(self.num_graphs)
        # Build block tasks (each block is a submatrix defined by index ranges)
        blocks = []
        for (start_ridx, end_ridx), (start_cidx, end_cidx) in product(start_end_indices, start_end_indices):
            if start_ridx <= start_cidx and end_ridx <= end_cidx:
                blocks.append((start_ridx, end_ridx, start_cidx, end_cidx))

        start_time = time.time()
        # Initialize worker pool with precomputed, read-only data to benefit from fork+COW
        nproc = get_recommended_nproc()
        pool = Pool(
            processes = nproc,
            initializer = _worker_init,
            initargs = (
                self.vtx_hierarchical_cycle_contexts_list,
                self.vtx_triangulated_neighbors_list,
                self.deg_distr_list,
                self.node_neighs_list,
                self.vlabel_list,
                self.edges_list,
                self.elabel_list,
                self._dataset_info,
                self._niter_tn,
                self._niter_hcc,
            )
        )

        try:
            for result in pool.imap_unordered(_compute_block_worker, blocks, chunksize=1):
                start_ridx, end_ridx, start_cidx, end_cidx, ot_dist_slice, wl_sim_slice = result
                ot_dist_diag[start_ridx: end_ridx, start_cidx: end_cidx] = ot_dist_slice
                wl_sim_diag[start_ridx: end_ridx, start_cidx: end_cidx] = wl_sim_slice
        finally:
            pool.close()
            pool.join()
        ot_dist_np = ot_dist_diag + ot_dist_diag.transpose(1, 0)
        wl_sim_np = wl_sim_diag + wl_sim_diag.transpose(1, 0)
        if self._wl_normalized:
            wl_sim_np_diag = np.diag(wl_sim_np)
            wl_sim_np = wl_sim_np / np.sqrt(np.outer(wl_sim_np_diag, wl_sim_np_diag))
        end_time = time.time()
        print(f'TopoWassersteinGraphKernel transform finished, Duration: {end_time - start_time:.4f} seconds')
        return ot_dist_np, wl_sim_np, end_time - start_time
    
    def fit_transform(
        self, dataset_info, graph_list, vlabel_list, edges_list, elabel_list, 
        vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list, 
        deg_distr_list
    ):
        self.fit(
            dataset_info, graph_list, vlabel_list, edges_list, elabel_list, 
            vtx_hierarchical_cycle_contexts_list, vtx_triangulated_neighbors_list, 
            deg_distr_list
        )
        return self.transform()
