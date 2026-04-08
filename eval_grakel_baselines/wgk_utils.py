from sklearn.model_selection import train_test_split
import numpy as np
import torch
import dgl
import collections
import ot

class WassersteinGraphKernel:
    def __init__(
        self, len_rw_trace, num_rw_trace, is_normalized, 
        p = 2, eigval_ndigits = 4, prob_ndigits = 6
    ):
        self._len_rw_trace = len_rw_trace
        self._num_rw_trace = num_rw_trace
        self._is_normalized = is_normalized
        self._p = p
        self._eigval_ndigits = eigval_ndigits
        self._prob_ndigits = prob_ndigits
        
    def comp_graph_matrix(self, dgl_graph):
        adj = dgl_graph.adj().to_dense()
        lap = torch.diag(adj.sum(dim = 1)) - adj
        
        return torch.block_diag(adj, lap)

    def comp_graph_eigval_distribution(self, dgl_graph, rw_trace):
        num_rws = rw_trace.shape[2]
        eigval_list = []
        for v in dgl_graph.nodes():
            v = v.item()
            for l in range(num_rws):
                ctx_nodes = torch.unique(rw_trace[v, :, l])
                ctx_subg = dgl.node_subgraph(dgl_graph, ctx_nodes)
                svd_vals = torch.linalg.svdvals(
                    self.comp_graph_matrix(ctx_subg)
                ).tolist()
                eigval_list.extend(svd_vals)
        eigval_distribution = collections.Counter(
            [round(ev, self._eigval_ndigits) for ev in eigval_list]
        )
        total_num = float(sum(eigval_distribution.values()))
        return {
            eigval: round(num/total_num, self._prob_ndigits) 
            for eigval, num in eigval_distribution.items()
        }

    def comp_graph_dist(
        self, src_eigval_list, src_prob_list, 
        tgt_eigval_list, tgt_prob_list
    ):
        return ot.wasserstein_1d(
            src_eigval_list, tgt_eigval_list, 
            src_prob_list, tgt_prob_list, 
            self._p
        ).item()
    
    def fit(self, tgt_graph_list):
        graph_eigval_distribution_list = []
        for g in tgt_graph_list:
            rw_trace = self.gen_ctx_trace_nodes(
                g, self._len_rw_trace, self._num_rw_trace
            )
            eigval_distribution = self.comp_graph_eigval_distribution(
                g, rw_trace, self._eigval_ndigits, self._prob_ndigits
            )
            self._graph_eigval_distribution_list.append(eigval_distribution)
        
        return graph_eigval_distribution_list
        
    def fit_transform(self, tgt_graph_list):
        self._tgt_graph_eigval_distribution_list = self.fit(tgt_graph_list)
        self._tgt_gnum = len(tgt_graph_list)
        graph_dist = torch.zeros(
            self._tgt_gnum, self._tgt_gnum
        )
        for ridx in range(self._tgt_gnum):
            row_eigval_distribution = \
                self._tgt_graph_eigval_distribution_list[ridx]
            row_eigval_list = np.array(
                list(row_eigval_distribution.keys())
            ).reshape(-1, 1)
            row_prob_list = np.array(
                list(row_eigval_distribution.values())
            ).reshape(-1, 1)
            for cidx in range(ridx + 1, self._tgt_gnum):
                col_eigval_distribution = \
                    self._tgt_graph_eigval_distribution_list[cidx]
                col_eigval_list = np.array(
                    list(col_eigval_distribution.keys())
                ).reshape(-1, 1)
                col_prob_list = np.array(
                    list(col_eigval_distribution.values())
                ).reshape(-1, 1)
                graph_dist[ridx, cidx] \
                    = graph_dist[cidx, ridx] \
                    = self.comp_graph_dist(
                        row_eigval_list, row_prob_list, 
                        col_eigval_list, col_prob_list
                    )

        gamma = 1.0 / self._tgt_gnum
        if self._is_normalized:
            row_sum = graph_dist.sum(axis = 1)
            graph_dist /= torch.sqrt(torch.outer(row_sum, row_sum))
        tgt_graph_sim_mat = torch.exp(
            -1.0 * gamma * torch.pow(graph_dist, self._p)
        )

        return tgt_graph_sim_mat

    def transform(self, src_graph_list):
        src_gnum = len(src_graph_list)
        src_graph_eigval_distribution_list = self.fit(src_graph_list)
        graph_dist = torch.zeros(
            src_gnum, self._tgt_gnum
        )
        for ridx in range(src_gnum):
            row_eigval_distribution = src_graph_eigval_distribution_list[ridx]
            row_eigval_list = np.array(
                list(row_eigval_distribution.keys())
            ).reshape(-1, 1)
            row_prob_list = np.array(
                list(row_eigval_distribution.values())
            ).reshape(-1, 1)
            for cidx in range(self._tgt_gnum):
                col_eigval_distribution = \
                    self._tgt_graph_eigval_distribution_list[cidx]
                col_eigval_list = np.array(
                    list(col_eigval_distribution.keys())
                ).reshape(-1, 1)
                col_prob_list = np.array(
                    list(col_eigval_distribution.values())
                ).reshape(-1, 1)
                graph_dist[ridx, cidx] \
                    = graph_dist[cidx, ridx] \
                    = self.comp_graph_dist(
                        row_eigval_list, row_prob_list, 
                        col_eigval_list, col_prob_list
                    )

        gamma = 1.0 / self._tgt_gnum
        if self._is_normalized:
            row_sum = graph_dist.sum(axis = 1)
            graph_dist /= torch.sqrt(torch.outer(row_sum, row_sum))
        graph_sim_mat = torch.exp(
            -1.0 * gamma * torch.pow(graph_dist, self._p)
        )
    
        return graph_sim_mat

class GromovWassersteinGraphDist:
    def __init__(self, normalized):
        self._normalized = normalized
    
    def fit(self, tgt_graph_list):
        self._node_deg_list = []
        for tgt_graph in tgt_graph_list:
            pass
    
    def fit_transform(self, tgt_graph_list):
        pass
    
    def transform(self, src_graph_list):
        pass
 
def gen_train_test_dataset(X, y):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size = 0.2, random_state = 42
    )
    return X_train, X_test, np.array(y_train), np.array(y_test)

def get_param_setting():
    param_setting = {
        'is_normalized': True, 
        'p': 2, 
        'eigval_ndigits': 6, 
        'prob_ndigits': 6
    }
    
    return param_setting
  