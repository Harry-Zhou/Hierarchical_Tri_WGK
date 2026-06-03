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
from cyclic_schema.hierarchical_triangulated_wl import hierarchical_triangular_wl_unified


def _build_elabel_dict(edges, elabels):
    return {tuple(e): int(elabels[i]) for i, e in enumerate(edges)}


# Module-level globals populated in worker initializer to avoid passing large objects per-task
_GLOBAL_DEG_DISTR = None
_GLOBAL_NODE_NEIGHS = None
_GLOBAL_VLABEL_LIST = None
_GLOBAL_EDGES_LIST = None
_GLOBAL_ELABEL_LIST = None
_GLOBAL_DATASET_INFO = None
_GLOBAL_NITER_TN = None
_GLOBAL_NITER_HCC = None
_GLOBAL_GRAPH_LIST = None

def _worker_init(deg_distr_list, node_neighs_list,
                 vlabel_list, edges_list, elabel_list, dataset_info, niter_tn, niter_hcc,
                 graph_list):
    global _GLOBAL_DEG_DISTR, _GLOBAL_NODE_NEIGHS
    global _GLOBAL_VLABEL_LIST, _GLOBAL_EDGES_LIST, _GLOBAL_ELABEL_LIST
    global _GLOBAL_DATASET_INFO, _GLOBAL_NITER_TN, _GLOBAL_NITER_HCC
    global _GLOBAL_GRAPH_LIST
    _GLOBAL_DEG_DISTR = deg_distr_list
    _GLOBAL_NODE_NEIGHS = node_neighs_list
    _GLOBAL_VLABEL_LIST = vlabel_list
    _GLOBAL_EDGES_LIST = edges_list
    _GLOBAL_ELABEL_LIST = elabel_list
    _GLOBAL_DATASET_INFO = dataset_info
    _GLOBAL_NITER_TN = niter_tn
    _GLOBAL_NITER_HCC = niter_hcc
    _GLOBAL_GRAPH_LIST = graph_list

def _compute_block_worker(block):
    # Globals are populated by _worker_init before any worker executes a block.
    # The asserts both document the contract and narrow the global types
    # for static analysis (Pylance/Pyright cannot trace Pool initializer side effects).
    assert (
        _GLOBAL_DEG_DISTR is not None and _GLOBAL_NODE_NEIGHS is not None
        and _GLOBAL_VLABEL_LIST is not None and _GLOBAL_EDGES_LIST is not None
        and _GLOBAL_ELABEL_LIST is not None
        and _GLOBAL_DATASET_INFO is not None
        and _GLOBAL_NITER_TN is not None and _GLOBAL_NITER_HCC is not None
        and _GLOBAL_GRAPH_LIST is not None
    ), "Worker globals not initialized; _worker_init must run before _compute_block_worker"
    _g_info = _GLOBAL_DATASET_INFO
    _graph_list = _GLOBAL_GRAPH_LIST
    _vlabel_list = _GLOBAL_VLABEL_LIST
    _edges_list = _GLOBAL_EDGES_LIST
    _elabel_list = _GLOBAL_ELABEL_LIST
    _node_neighs = _GLOBAL_NODE_NEIGHS
    _deg_distr = _GLOBAL_DEG_DISTR

    has_el = _g_info.get('el', False)

    def _run_wl_test(ridx, cidx):
        G1 = _graph_list[ridx]
        G2 = _graph_list[cidx]
        rv = _vlabel_list[ridx]
        cv = _vlabel_list[cidx]
        if has_el:
            ed1 = _build_elabel_dict(_edges_list[ridx], _elabel_list[ridx])
            ed2 = _build_elabel_dict(_edges_list[cidx], _elabel_list[cidx])
        else:
            ed1, ed2 = None, None
        return hierarchical_triangular_wl_unified(
            _g_info, G1, G2, rv, cv, ed1, ed2,
            K=_GLOBAL_NITER_HCC, I=_GLOBAL_NITER_TN,
        )

    def _wl_counter(node_neighs, vwl_np, vmax_label):
        from collections import Counter as _Counter
        vwl_counter = np.zeros((vwl_np.shape[0], vmax_label + 1))
        for ridx, neighs in node_neighs.items():
            for l, n in _Counter(vwl_np[neighs, :].flatten('F')).items():
                vwl_counter[ridx, l] = n
        return vwl_counter

    start_ridx, end_ridx, start_cidx, end_cidx = block
    ot_dist_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
    wl_sim_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
    temp_start_cidx = start_cidx
    is_diagonal = start_ridx == start_cidx and end_ridx == end_cidx
    print(f'\t The ({start_ridx}, {end_ridx}, {start_cidx}, {end_cidx}) graph similarity submatrix')
    for gridx in range(start_ridx, end_ridx):
        dridx = gridx - start_ridx
        r_deg_distr = _deg_distr[gridx]
        if is_diagonal:
            temp_start_cidx = gridx + 1
        for gcidx in range(temp_start_cidx, end_cidx):
            rvwl_np, cvwl_np = _run_wl_test(gridx, gcidx)
            native_vmax_label = max(rvwl_np.max(), cvwl_np.max())
            native_vwl_counter1 = _wl_counter(_node_neighs[gridx], rvwl_np, native_vmax_label)
            native_vwl_counter2 = _wl_counter(_node_neighs[gcidx], cvwl_np, native_vmax_label)
            transport_cost = cdist(native_vwl_counter1, native_vwl_counter2, metric='sqeuclidean')
            dcidx = gcidx - start_cidx
            sw_dist = float(ot.sliced.sliced_wasserstein_distance(
                native_vwl_counter1, native_vwl_counter2,
                r_deg_distr, _deg_distr[gcidx],
                n_projections=50, log=False,
            ))
            ot_dist_diag[dridx, dcidx] = sw_dist ** 2
            wl_sim_diag[dridx, dcidx] = TopoWassersteinGraphKernel.wl_inner_product(rvwl_np, cvwl_np)
        if is_diagonal:
            rvwl_np, _ = _run_wl_test(gridx, gridx)
            wl_sim_diag[dridx, dridx] = TopoWassersteinGraphKernel.wl_inner_product(rvwl_np, rvwl_np)
    return (start_ridx, end_ridx, start_cidx, end_cidx, ot_dist_diag, wl_sim_diag)

class TopoWassersteinGraphKernel:
    def __init__(self, niter_tn, niter_hcc, wl_normalized = True):
        self._niter_tn = niter_tn
        self._niter_hcc = niter_hcc
        self._wl_normalized = wl_normalized
    
    @staticmethod
    def wl_inner_product(vwl_np1, vwl_np2):
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
        
        return transport_cost, native_vwl_counter1, native_vwl_counter2
    
    def _run_wl_test(self, ridx, cidx):
        G1 = self.graph_list[ridx]
        G2 = self.graph_list[cidx]
        rv = self.vlabel_list[ridx]
        cv = self.vlabel_list[cidx]
        if self._dataset_info.get('el', False):
            ed1 = _build_elabel_dict(self.edges_list[ridx], self.elabel_list[ridx])
            ed2 = _build_elabel_dict(self.edges_list[cidx], self.elabel_list[cidx])
        else:
            ed1, ed2 = None, None
        return hierarchical_triangular_wl_unified(
            self._dataset_info, G1, G2, rv, cv, ed1, ed2,
            K=self._niter_hcc, I=self._niter_tn,
        )

    def _otdist_wlsim_worker(self, start_ridx, end_ridx, start_cidx, end_cidx):
        ot_dist_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
        wl_sim_diag = np.zeros((end_ridx - start_ridx, end_cidx - start_cidx))
        temp_start_cidx = start_cidx
        is_diagonal = start_ridx == start_cidx and end_ridx == end_cidx
        print(f'\t The ({start_ridx}, {end_ridx}, {start_cidx}, {end_cidx}) graph similarity submatrix')
        for gridx in range(start_ridx, end_ridx):
            dridx = gridx - start_ridx
            r_deg_distr = self.deg_distr_list[gridx]
            r_node_neighs = self.node_neighs_list[gridx]
            if is_diagonal:
                temp_start_cidx = gridx + 1
            for gcidx in range(temp_start_cidx, end_cidx):
                rvwl_np, cvwl_np = self._run_wl_test(gridx, gcidx)
                _, vwl_c1, vwl_c2 = self._comp_transport_cost(
                    r_node_neighs, self.node_neighs_list[gcidx], rvwl_np, cvwl_np
                )
                dcidx = gcidx - start_cidx
                sw_dist = float(ot.sliced.sliced_wasserstein_distance(
                    vwl_c1, vwl_c2,
                    r_deg_distr, self.deg_distr_list[gcidx],
                    n_projections=50, log=False,
                ))
                ot_dist_diag[dridx, dcidx] = sw_dist ** 2
                wl_sim_diag[dridx, dcidx] = self.wl_inner_product(rvwl_np, cvwl_np)

            if is_diagonal:
                rvwl_np, _ = self._run_wl_test(gridx, gridx)
                wl_sim_diag[dridx, dridx] = self.wl_inner_product(rvwl_np, rvwl_np)
        return (start_ridx, end_ridx, start_cidx, end_cidx, ot_dist_diag, wl_sim_diag)

    def fit(
        self, dataset_info, graph_list, vlabel_list, edges_list, elabel_list, 
        deg_distr_list
    ):
        self._dataset_info = dataset_info
        self._nl, self._el = dataset_info['nl'], dataset_info['el']
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
                self.deg_distr_list,
                self.node_neighs_list,
                self.vlabel_list,
                self.edges_list,
                self.elabel_list,
                self._dataset_info,
                self._niter_tn,
                self._niter_hcc,
                self.graph_list,
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
        deg_distr_list
    ):
        self.fit(
            dataset_info, graph_list, vlabel_list, edges_list, elabel_list, 
            deg_distr_list
        )
        return self.transform()
