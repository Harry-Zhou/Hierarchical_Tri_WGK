import networkx as nx
from multiprocessing import cpu_count
try:
    import psutil
except Exception:
    psutil = None

def get_start_end_indices(num_samples, num_cells = cpu_count()):
    assert num_samples >= 1
    if num_samples <= num_cells:
        start_end_indices = [
            (start_idx, end_idx) for start_idx, end_idx in zip(range(num_samples), range(1, num_samples + 1))
        ]
    else:
        cell_size = num_samples // num_cells
        remainder = num_samples % num_cells
        if remainder == 0:
            start_end_indices = [
                (start_idx * cell_size, end_idx * cell_size) \
                    for start_idx, end_idx in zip(range(num_cells), range(1, num_cells + 1))
            ]
        else:
            header_start_end_indices = [
                (start_idx * (cell_size + 1), end_idx * (cell_size + 1)) \
                    for start_idx, end_idx in zip(range(remainder), range(1, remainder + 1))
            ]
            header_tail_idx = header_start_end_indices[-1][-1]
            tail_start_end_indices = []
            for i, j in zip(range(num_cells - remainder), range(1, num_cells - remainder + 1)):
                if header_tail_idx + j * cell_size <= num_samples:
                    tail_start_end_indices.append(
                        (header_tail_idx + i * cell_size, header_tail_idx + j * cell_size)
                    )
            tail_start_end_indices[-1][-1]
            if tail_start_end_indices[-1][-1] < num_samples:
                tail_start_end_indices.append((tail_start_end_indices[-1][-1], num_samples))
            start_end_indices = header_start_end_indices + tail_start_end_indices
        
    return start_end_indices

def precomp_ego_nets_list(graph_list, radius):
    ego_nets_list = []
    for g in graph_list:
        ego_nets = [None, ] * g.number_of_nodes()
        for v in g.nodes():
            ego_nets[v] = nx.ego_graph(g, v, radius = radius)
        ego_nets_list.append(ego_nets)
    return ego_nets_list

def precomp_node_neighs(g):
    node_neighs = {vtx: list(g.neighbors(vtx)) for vtx in g.nodes()}
    return node_neighs


def get_recommended_nproc(reserve_cores=1, est_mem_per_proc_gb=None):
    """Return a recommended number of worker processes.

    - reserve_cores: keep this many physical cores free for system/overhead
    - est_mem_per_proc_gb: if provided, bound processes by available RAM
    """
    try:
        physical = psutil.cpu_count(logical=False) or cpu_count()
    except Exception:
        physical = cpu_count()

    nproc = max(1, physical - reserve_cores)
    if est_mem_per_proc_gb is not None and psutil is not None:
        try:
            total_mem_gb = psutil.virtual_memory().total / (1024 ** 3)
            mem_limit = int(total_mem_gb // est_mem_per_proc_gb)
            if mem_limit >= 1:
                nproc = min(nproc, mem_limit)
        except Exception:
            pass
    return nproc

if __name__ == '__main__':
    start_end_indices = get_start_end_indices(10)
    for sei in start_end_indices:
        print(sei)
