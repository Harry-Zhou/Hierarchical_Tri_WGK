"""
Dataset loading and preprocessing for CSG-Transformer evaluation.

Supports:
- Graph classification: MUTAG, PROTEINS, NCI1, NCI109, ENZYMES, D&D,
  IMDB-BINARY, IMDB-MULTI, COLLAB, REDDIT-BINARY
- Node classification: Cora, Citeseer, Pubmed, Computers, Photo, Squirrel, Chameleon
- Custom datasets through TUDataset
"""

import os
from typing import Any, Callable, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
import torch
from torch_geometric.datasets import TUDataset, Planetoid, Amazon, Coauthor, WebKB
from torch_geometric.utils import to_networkx


GRAPH_CLASSIFICATION_DATASETS = {
    'MUTAG': {'num_graphs': 188, 'num_classes': 2, 'avg_nodes': 17.9, 'in_dim': 7},
    'PROTEINS': {'num_graphs': 1113, 'num_classes': 2, 'avg_nodes': 39.1, 'in_dim': 3},
    'NCI1': {'num_graphs': 4110, 'num_classes': 2, 'avg_nodes': 29.8, 'in_dim': 37},
    'NCI109': {'num_graphs': 4110, 'num_classes': 2, 'avg_nodes': 29.6, 'in_dim': 38},
    'ENZYMES': {'num_graphs': 600, 'num_classes': 6, 'avg_nodes': 32.6, 'in_dim': 3},
    'D&D': {'num_graphs': 1178, 'num_classes': 2, 'avg_nodes': 284.3, 'in_dim': 89},
    'IMDB-BINARY': {'num_graphs': 1000, 'num_classes': 2, 'avg_nodes': 19.8, 'in_dim': 0},
    'IMDB-MULTI': {'num_graphs': 1500, 'num_classes': 3, 'avg_nodes': 13.0, 'in_dim': 0},
    'COLLAB': {'num_graphs': 5000, 'num_classes': 3, 'avg_nodes': 74.5, 'in_dim': 0},
    'REDDIT-BINARY': {'num_graphs': 2000, 'num_classes': 2, 'avg_nodes': 429.6, 'in_dim': 0},
}

NODE_CLASSIFICATION_DATASETS = {
    'Cora': {'num_nodes': 2708, 'num_edges': 5429, 'num_classes': 7, 'in_dim': 1433},
    'Citeseer': {'num_nodes': 3327, 'num_edges': 4732, 'num_classes': 6, 'in_dim': 3703},
    'Pubmed': {'num_nodes': 19717, 'num_edges': 44338, 'num_classes': 3, 'in_dim': 500},
    'Computers': {'num_nodes': 13752, 'num_edges': 245861, 'num_classes': 10, 'in_dim': 767},
    'Photo': {'num_nodes': 7650, 'num_edges': 119081, 'num_classes': 8, 'in_dim': 745},
    'Squirrel': {'num_nodes': 5201, 'num_edges': 217073, 'num_classes': 5, 'in_dim': 2089},
    'Chameleon': {'num_nodes': 2277, 'num_edges': 36052, 'num_classes': 5, 'in_dim': 3233},
}


def _pyg_to_networkx(data: Any) -> nx.Graph:
    """Convert PyG Data object to NetworkX graph."""
    if hasattr(data, 'edge_index') and data.edge_index.numel() > 0:
        G = nx.Graph()
        G.add_nodes_from(range(data.num_nodes))
        G.add_edges_from(data.edge_index.t().tolist())
        return G
    return nx.empty_graph(data.num_nodes)


def load_graph_classification_dataset(
    name: str, root: str = './data'
) -> Tuple[List[Tuple[nx.Graph, torch.Tensor, torch.Tensor]], Dict[str, Any]]:
    """
    Load graph classification dataset.
    
    Args:
        name: Dataset name (MUTAG, PROTEINS, etc.)
        root: Root directory for data storage
    
    Returns:
        data_list: List of (G, features, label) tuples
        stats: Dataset statistics dict
    """
    dataset = TUDataset(root=f'{root}/TUDataset', name=name)
    dataset = dataset.shuffle()
    
    data_list = []
    for data in dataset:
        G = _pyg_to_networkx(data)
        if hasattr(data, 'x') and data.x is not None and data.num_features > 0:
            features = data.x.float()
        else:
            # Featureless graph: use degree-based features
            deg = torch.tensor([G.degree(v) for v in G.nodes()], dtype=torch.float32).unsqueeze(1)
            # Also add normalized degree and clustering coefficient as features
            max_deg = deg.max().item() if deg.numel() > 0 and deg.max().item() > 0 else 1.0
            norm_deg = deg / max_deg
            features = torch.cat([deg, norm_deg], dim=1)
        label = data.y
        data_list.append((G, features, label))
    
    # Determine actual input dimension from features
    actual_in_dim = data_list[0][1].shape[1] if data_list else 1
    
    stats = {
        'num_graphs': len(dataset),
        'num_classes': dataset.num_classes,
        'in_dim': actual_in_dim,
        'avg_nodes': np.mean([d[0].number_of_nodes() for d in data_list]),
        'max_nodes': max(d[0].number_of_nodes() for d in data_list),
        'min_nodes': min(d[0].number_of_nodes() for d in data_list),
    }
    
    return data_list, stats


def load_node_classification_dataset(
    name: str, root: str = './data'
) -> Tuple[nx.Graph, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Load node classification dataset.
    
    Args:
        name: Dataset name (Cora, Citeseer, etc.)
        root: Root directory for data storage
    
    Returns:
        G: NetworkX graph
        features: Node features
        labels: Node labels
        train_mask: Training mask
        val_mask: Validation mask
        test_mask: Test mask
    """
    dataset_name_map = {
        'Cora': 'cora', 'Citeseer': 'citeseer', 'Pubmed': 'pubmed',
    }
    
    if name in ['Cora', 'Citeseer', 'Pubmed']:
        pyg_name = dataset_name_map[name]
        dataset = Planetoid(root=f'{root}/Planetoid', name=pyg_name)
    elif name in ['Computers', 'Photo']:
        dataset = Amazon(root=f'{root}/Amazon', name=name)
    elif name in ['Squirrel', 'Chameleon']:
        dataset = WebKB(root=f'{root}/WebKB', name=name)
    else:
        raise ValueError(f"Dataset {name} not supported")
    
    data = dataset[0]
    G = _pyg_to_networkx(data)
    
    features = data.x
    labels = data.y
    
    train_mask = data.train_mask if hasattr(data, 'train_mask') else None
    val_mask = data.val_mask if hasattr(data, 'val_mask') else None
    test_mask = data.test_mask if hasattr(data, 'test_mask') else None
    
    if train_mask is None:
        num_nodes = data.num_nodes
        perm = torch.randperm(num_nodes)
        train_mask = torch.zeros(num_nodes, dtype=torch.bool)
        val_mask = torch.zeros(num_nodes, dtype=torch.bool)
        test_mask = torch.zeros(num_nodes, dtype=torch.bool)
        train_mask[perm[:int(0.6 * num_nodes)]] = True
        val_mask[perm[int(0.6 * num_nodes):int(0.8 * num_nodes)]] = True
        test_mask[perm[int(0.8 * num_nodes):]] = True
    
    return G, features, labels, train_mask, val_mask, test_mask


def get_dataset_stats(
    data_list: List[Tuple],
    task: str = 'graph_classification'
) -> Dict[str, Any]:
    """Get dataset statistics from loaded data."""
    if task == 'graph_classification':
        num_graphs = len(data_list)
        num_classes = len(set(d[2].item() if torch.is_tensor(d[2]) else d[2] for d in data_list))
        node_counts = [d[0].number_of_nodes() for d in data_list]
        in_dim = data_list[0][1].size(-1) if len(data_list) > 0 else 1
        
        return {
            'num_graphs': num_graphs,
            'num_classes': num_classes,
            'in_dim': in_dim,
            'avg_nodes': float(np.mean(node_counts)),
            'max_nodes': max(node_counts),
            'min_nodes': min(node_counts),
        }
    else:
        G = data_list[0] if isinstance(data_list, list) and len(data_list) > 0 else data_list
        return {
            'num_nodes': G.number_of_nodes(),
            'in_dim': data_list[1].size(-1) if isinstance(data_list, list) and len(data_list) > 1 else 1,
        }
