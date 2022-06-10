from .graph_dataset import (
    GraphClassificationDataset,
    GraphClassificationDatasetLabeled,
    LoadBalanceGraphDataset,
    NodeClassificationDataset,
    NodeClassificationDatasetLabeled,
    worker_init_fn,
)
        
# GRAPH_CLASSIFICATION_DSETS = ["collab", "imdb-binary", "imdb-multi", "rdt-b", "rdt-5k"]
GRAPH_CLASSIFICATION_DSETS = ['reddit_multi_12k', 'reddit_multi_5k', 'reddit_binary', \
    'github_stargazers', 'collab', 'dd', 'proteins_full', 'uacc257h', 'mutagenicity', \
    'nci-h23h', 'p388h', 'pc-3h', 'sn12ch']

__all__ = [
    "GRAPH_CLASSIFICATION_DSETS",
    "LoadBalanceGraphDataset",
    "GraphClassificationDataset",
    "GraphClassificationDatasetLabeled",
    "NodeClassificationDataset",
    "NodeClassificationDatasetLabeled",
    "worker_init_fn",
]
