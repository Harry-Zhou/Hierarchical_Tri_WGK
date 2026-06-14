"""
Dataset loading and preprocessing for CSG-Transformer evaluation.

Supports:
- Graph classification: MUTAG, PROTEINS, NCI1, etc.
- Node classification: Cora, Citeseer, Pubmed, etc.
"""

import torch
from torch_geometric.datasets import TUDataset, Planetoid
from torch_geometric.transforms import NormalizeFeatures
from torch_geometric.loader import DataLoader
import numpy as np


def load_graph_classification_dataset(name, root='./data'):
    """
    Load graph classification dataset.
    
    Args:
        name: Dataset name (MUTAG, PROTEINS, NCI1, etc.)
        root: Root directory for data storage
    
    Returns:
        train_dataset, val_dataset, test_dataset
    """
    dataset = TUDataset(root=f'{root}/TUDataset', name=name)
    
    dataset = dataset.shuffle()
    
    num_graphs = len(dataset)
    train_size = int(0.8 * num_graphs)
    val_size = int(0.1 * num_graphs)
    test_size = num_graphs - train_size - val_size
    
    train_dataset = dataset[:train_size]
    val_dataset = dataset[train_size:train_size + val_size]
    test_dataset = dataset[train_size + val_size:]
    
    return train_dataset, val_dataset, test_dataset


def load_node_classification_dataset(name, root='./data'):
    """
    Load node classification dataset.
    
    Args:
        name: Dataset name (Cora, Citeseer, Pubmed, etc.)
        root: Root directory for data storage
    
    Returns:
        dataset with train_mask, val_mask, test_mask
    """
    if name in ['Cora', 'Citeseer', 'Pubmed']:
        dataset = Planetoid(
            root=f'{root}/Planetoid',
            name=name,
            transform=NormalizeFeatures()
        )
    else:
        raise ValueError(f"Dataset {name} not supported")
    
    return dataset


def create_data_loaders(train_dataset, val_dataset, test_dataset, 
                        batch_size=64, task='graph_classification'):
    """
    Create data loaders for training, validation, and testing.
    
    Args:
        train_dataset: Training dataset
        val_dataset: Validation dataset
        test_dataset: Test dataset
        batch_size: Batch size
        task: 'graph_classification' or 'node_classification'
    
    Returns:
        train_loader, val_loader, test_loader
    """
    if task == 'graph_classification':
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    else:
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    return train_loader, val_loader, test_loader


def get_dataset_stats(dataset, task='graph_classification'):
    """
    Get dataset statistics.
    
    Args:
        dataset: Dataset
        task: 'graph_classification' or 'node_classification'
    
    Returns:
        Dictionary with dataset statistics
    """
    if task == 'graph_classification':
        num_graphs = len(dataset)
        num_classes = dataset.num_classes
        
        node_counts = [g.num_nodes for g in dataset]
        edge_counts = [g.num_edges for g in dataset]
        
        in_dim = dataset.num_features
        
        return {
            'num_graphs': num_graphs,
            'num_classes': num_classes,
            'in_dim': in_dim,
            'avg_nodes': np.mean(node_counts),
            'avg_edges': np.mean(edge_counts),
            'max_nodes': max(node_counts),
            'min_nodes': min(node_counts)
        }
    else:
        data = dataset[0]
        return {
            'num_nodes': data.num_nodes,
            'num_edges': data.num_edges,
            'num_classes': dataset.num_classes,
            'in_dim': data.x.shape[1],
            'train_mask': data.train_mask.sum().item(),
            'val_mask': data.val_mask.sum().item(),
            'test_mask': data.test_mask.sum().item()
        }


GRAPH_CLASSIFICATION_DATASETS = {
    'MUTAG': {'num_graphs': 188, 'num_classes': 2, 'avg_nodes': 17.9},
    'PROTEINS': {'num_graphs': 1113, 'num_classes': 2, 'avg_nodes': 39.1},
    'NCI1': {'num_graphs': 4110, 'num_classes': 2, 'avg_nodes': 29.8},
    'NCI109': {'num_graphs': 4110, 'num_classes': 2, 'avg_nodes': 29.6},
    'ENZYMES': {'num_graphs': 600, 'num_classes': 6, 'avg_nodes': 32.6},
    'D&D': {'num_graphs': 1178, 'num_classes': 2, 'avg_nodes': 284.3},
    'IMDB-B': {'num_graphs': 1000, 'num_classes': 2, 'avg_nodes': 19.8},
    'IMDB-M': {'num_graphs': 1500, 'num_classes': 3, 'avg_nodes': 13.0},
    'COLLAB': {'num_graphs': 5000, 'num_classes': 3, 'avg_nodes': 74.5},
    'REDDIT-B': {'num_graphs': 2000, 'num_classes': 2, 'avg_nodes': 429.6},
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
